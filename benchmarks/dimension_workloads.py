#!/usr/bin/env python3
"""Compare basis, scalar-dimension, and batch-dimension module APIs.

Run each mode in a separate process against the same release build. Consumer
fixtures are read-only: --mode dimensions/batch rewrites only len(hom_basis_to)
and ext1(...).dimension expressions in the two explicitly loaded source files.
No imported helpers or consumer files on disk are changed. Batch mode uses the
batch API in ext_grid; consumer expressions still use scalar dimension calls.
Timings include module construction. Imports and AST rewriting are excluded.
"""
import argparse
import ast
import hashlib
import json
import os
import platform
import statistics
import sys
import time
from pathlib import Path

import qpa
import qpa._qpa as extension
from module_workloads import kronecker


class Dimensions(ast.NodeTransformer):
    def visit_Call(self, node):
        self.generic_visit(node)
        if (isinstance(node.func, ast.Name) and node.func.id == 'len'
                and len(node.args) == 1 and not node.keywords):
            inner = node.args[0]
            if isinstance(inner, ast.Call) and isinstance(inner.func, ast.Attribute) and inner.func.attr == 'hom_basis_to':
                inner.func.attr = 'hom_dimension_to'
                return ast.copy_location(inner, node)
        return node

    def visit_Attribute(self, node):
        self.generic_visit(node)
        inner = node.value
        if (node.attr == 'dimension' and isinstance(inner, ast.Call)
                and isinstance(inner.func, ast.Attribute) and inner.func.attr == 'ext1'):
            inner.func.attr = 'ext1_dimension'
            return ast.copy_location(inner, node)
        return node


def cases(mode, root):
    def hom(field):
        module = kronecker(field, 24)
        result = len(module.hom_basis_to(module)) if mode == 'basis' else module.hom_dimension_to(module)
        assert result == 24
        return result

    def ext_grid():
        algebra = qpa.QuotientAlgebra.truncated(qpa.Field.prime(2), qpa.Quiver.indexed(1, [(1, 1, 'x'), (1, 1, 'y')]), 2)
        simple = algebra.simple_modules()[0]
        projective = algebra.indecomposable_projectives()[0]
        source = qpa.MatrixModule.direct_sum([simple] * 4 + [projective] * 2)
        targets = [qpa.MatrixModule.direct_sum([simple] * n + [projective] * p)
                   for n, p in [(4, 0), (2, 2), (0, 4), (3, 1), (1, 3), (2, 0)]]
        if mode == 'basis':
            return [source.ext1(target).dimension for target in targets]
        if mode == 'dimensions':
            return [source.ext1_dimension(target) for target in targets]
        return source.ext1_dimensions(targets)

    result = {'hom_rational': lambda: hom(qpa.Field.rational()),
              'hom_binary': lambda: hom(qpa.Field.prime(2)), 'ext_grid': ext_grid}
    if root:
        sys.path.insert(0, str(root))

        def load(relative):
            path = root / relative
            tree = ast.parse(path.read_text(), filename=str(path))
            if mode != 'basis':
                tree = ast.fix_missing_locations(Dimensions().visit(tree))
            scope = {'__name__': 'benchmark_fixture', '__file__': str(path)}
            exec(compile(tree, str(path), 'exec'), scope)
            return scope

        js_path = 'experiments/weak-gorenstein-symmetry-counterexample/verify_js29.py'
        js = load(js_path)
        certificate = json.loads((root / js_path).with_name('js29-input.json').read_text())
        graded = load('lib/graded_square_zero_qpa.py')

        def graded_profiles():
            specs = [graded['thin_complex_specs']()[-1]] + graded['kronecker_complex_specs'](1)
            data = graded['graded_profiles'](2, specs)
            assert data['dimension'] == 9
            assert data['simple_pd'] == [3, 2, 1, 0]
            assert [data['profiles'][i][i] for i in range(3)] == [[2, 1, 0], [0, 1, 0], [0, 1, 0]]
            return data

        result.update(consumer_js29_rational=lambda: js['verify'](certificate),
                      consumer_graded_ext_binary=graded_profiles)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['basis', 'dimensions', 'batch'], required=True)
    parser.add_argument('--repeat', type=int, default=3)
    parser.add_argument('--case', action='append', help='Run only these named cases')
    parser.add_argument('--consumer-root', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compare', type=Path)
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error('--repeat must be positive')
    workloads = cases(args.mode, args.consumer_root.resolve() if args.consumer_root else None)
    if args.case:
        if unknown := set(args.case) - workloads.keys():
            parser.error(f'unknown cases: {sorted(unknown)}')
        workloads = {name: workloads[name] for name in args.case}
    baseline = json.loads(args.compare.read_text()) if args.compare else None
    report = {'mode': args.mode, 'python': sys.version, 'platform': platform.platform(),
              'cpu_affinity': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None,
              'qpa_package': qpa.__file__, 'qpa_build_revision': qpa.PORT_GIT_COMMIT,
              'extension_sha256': hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest(),
              'cases': {}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for name, run in workloads.items():
        samples, values = [], []
        for _ in range(args.repeat):
            start = time.perf_counter()
            value = run()
            samples.append(time.perf_counter() - start)
            values.append(json.loads(json.dumps(value)))
        assert all(value == values[0] for value in values), name
        entry = {'seconds': samples, 'median_seconds': statistics.median(samples), 'result': values[0]}
        if baseline:
            assert entry['result'] == baseline['cases'][name]['result'], name
            entry['speedup'] = baseline['cases'][name]['median_seconds'] / entry['median_seconds']
        report['cases'][name] = entry
        args.output.write_text(json.dumps(report, indent=2) + '\n')
        print(f"{name}: {entry['median_seconds']:.6f}s", flush=True)


if __name__ == '__main__':
    main()
