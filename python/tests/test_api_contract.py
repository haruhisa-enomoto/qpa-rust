"""The shipped signatures must describe the extension that users import."""

import ast
import inspect
from pathlib import Path
import types

import qpa


STUB = Path(__file__).resolve().parents[1] / 'qpa' / '_qpa.pyi'


def test_stub_members_and_calling_conventions_match_runtime():
    problems = []
    for node in ast.parse(STUB.read_text()).body:
        if isinstance(node, ast.ClassDef):
            owner = getattr(qpa, node.name, None)
            if owner is None:
                problems.append(f'{node.name}: missing class')
                continue
            members = [(f'{node.name}.{m.name}', owner, m)
                       for m in node.body if isinstance(m, ast.FunctionDef)]
            declared = {m.name for m in node.body if isinstance(m, ast.FunctionDef)}
            declared |= {m.target.id for m in node.body
                         if isinstance(m, ast.AnnAssign) and isinstance(m.target, ast.Name)}
            undocumented = {name for name in owner.__dict__
                            if not name.startswith('_') and name not in declared}
            if undocumented:
                problems.append(f'{node.name}: runtime members missing from stub {sorted(undocumented)}')
            names = [m.name for m in node.body if isinstance(m, ast.FunctionDef)
                     and not any(ast.unparse(d).endswith('.setter') or ast.unparse(d) == 'overload'
                                 for d in m.decorator_list)]
            if len(names) != len(set(names)):
                problems.append(f'{node.name}: duplicate stub members')
        elif isinstance(node, ast.FunctionDef):
            members = [(node.name, qpa, node)]
        else:
            continue
        for name, owner, member in members:
            descriptor = inspect.getattr_static(owner, member.name, None)
            if descriptor is None:
                problems.append(f'{name}: missing member')
                continue
            decorators = {ast.unparse(d) for d in member.decorator_list}
            is_property = 'property' in decorators or any(d.endswith('.setter') for d in decorators)
            runtime_property = isinstance(descriptor, (property, types.GetSetDescriptorType))
            if is_property != runtime_property:
                problems.append(f'{name}: method/property mismatch')
            if is_property or member.name.startswith('__') and member.name != '__init__':
                # CPython slot wrappers use generic parameter names for operators.
                continue
            if inspect.isclass(owner) and ('staticmethod' in decorators) != isinstance(descriptor, staticmethod):
                problems.append(f'{name}: static/instance method mismatch')
            try:
                signature = inspect.signature(owner if member.name == '__init__' else getattr(owner, member.name))
            except (ValueError, TypeError):
                continue  # Some builtin methods do not expose text signatures.
            positional = member.args.posonlyargs + member.args.args
            required = len(positional) - len(member.args.defaults)
            expected = [(arg.arg, i < required) for i, arg in enumerate(positional)]
            expected += [(arg.arg, default is None) for arg, default in
                         zip(member.args.kwonlyargs, member.args.kw_defaults)]
            expected = [(n, r) for n, r in expected if n not in ('self', 'cls')]
            actual = [(p.name, p.default is inspect.Parameter.empty) for p in signature.parameters.values()
                      if p.name not in ('self', 'cls')]
            if expected != actual:
                problems.append(f'{name}: parameters/defaults {expected} != {actual}')
                continue
            kinds = {arg.arg: inspect.Parameter.POSITIONAL_ONLY for arg in member.args.posonlyargs}
            kinds.update({arg.arg: inspect.Parameter.POSITIONAL_OR_KEYWORD for arg in member.args.args})
            kinds.update({arg.arg: inspect.Parameter.KEYWORD_ONLY for arg in member.args.kwonlyargs})
            for parameter, _ in expected:
                if kinds[parameter] != signature.parameters[parameter].kind:
                    problems.append(f'{name}.{parameter}: positional/keyword convention mismatch')
            defaults = list(zip(positional[required:], member.args.defaults))
            defaults += list(zip(member.args.kwonlyargs, member.args.kw_defaults))
            for arg, default in defaults:
                if default is None or isinstance(default, ast.Constant) and default.value is Ellipsis:
                    continue
                if ast.literal_eval(default) != signature.parameters[arg.arg].default:
                    problems.append(f'{name}.{arg.arg}: default value mismatch')
    assert not problems, '\n'.join(problems)


def test_public_classes_and_functions_have_stub_entries():
    declared = {node.name for node in ast.parse(STUB.read_text()).body
                if isinstance(node, (ast.ClassDef, ast.FunctionDef))}
    actual = {name for name in qpa.__all__
              if inspect.isclass(getattr(qpa, name)) or inspect.isbuiltin(getattr(qpa, name))}
    assert not actual - declared
