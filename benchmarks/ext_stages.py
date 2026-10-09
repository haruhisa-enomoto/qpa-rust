#!/usr/bin/env python3
"""Profile the cover, Hom, and restriction stages of the graded Ext fixture.

Load the consumer's graded_profiles function with its final calculation replaced
in memory by `return placed`. This preserves its module construction without
changing consumer files. Timings use the existing Python APIs and include their
wrapper overhead. Quotient-complement construction is not measured here; use
module_workloads.py for complete Ext timings. Covers and maps are discarded
between pairs, as in the consumer's ext1(...).dimension calls.
"""
import argparse
import ast
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import qpa
import qpa._qpa as extension


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--consumer-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.consumer_root.resolve()
    sys.path.insert(0, str(root))
    path = root / 'lib/graded_square_zero_qpa.py'
    source = path.read_text()
    tree = ast.parse(source)
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'graded_profiles')
    index = next(i for i, n in enumerate(function.body)
                 if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'profiles' for t in n.targets))
    function.body = function.body[:index] + [ast.Return(ast.Name('placed', ast.Load()))]
    scope = {'__name__': 'benchmark_fixture', '__file__': str(path)}
    exec(compile(ast.fix_missing_locations(tree), str(path), 'exec'), scope)
    specs = [scope['thin_complex_specs']()[-1]] + scope['kronecker_complex_specs'](1)
    placed = scope['graded_profiles'](2, specs)
    timings = {key: 0.0 for key in ['cover', 'kernel', 'projective_hom', 'syzygy_hom', 'compose']}

    def timed(stage, operation):
        start = time.perf_counter()
        value = operation()
        timings[stage] += time.perf_counter() - start
        return value

    counts = []
    for sources in placed:
        for targets in placed:
            for target in reversed(targets):
                module = sources[0]
                cover = timed('cover', module.projective_cover)
                inclusion = timed('kernel', cover.kernel_inclusion)
                projective_hom = timed('projective_hom', lambda: cover.source.hom_basis_to(target))
                syzygy_hom = timed('syzygy_hom', lambda: inclusion.source.hom_basis_to(target))
                restrictions = timed('compose', lambda: [inclusion.compose(h) for h in projective_hom])
                counts.append([len(projective_hom), len(syzygy_hom), cover.source.dimension,
                               target.dimension, inclusion.source.dimension])
                del restrictions, projective_hom, syzygy_hom, inclusion, cover
    report = {'source_sha256': hashlib.sha256(source.encode()).hexdigest(),
              'qpa_build_revision': qpa.PORT_GIT_COMMIT,
              'extension_sha256': hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest(),
              'cpu_affinity': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None,
              'seconds': timings, 'total_seconds': sum(timings.values()), 'calls': len(counts),
              'hom_dimensions_and_module_sizes': counts}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(timings, indent=2))


if __name__ == '__main__':
    main()
