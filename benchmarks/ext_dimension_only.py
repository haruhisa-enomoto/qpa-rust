#!/usr/bin/env python3
"""Ordinary-API Ext dimension benchmark with retained covers, across fields.

Preparation is timed separately. No algebra-specific engine configuration or
supplied bases are used. Expected dimensions come from truncated-loop modules
and regular Kronecker modules; nonzero restriction images are included.
Run under an external wall-time guard.
"""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import resource
import statistics
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeat', type=int, default=3)
    parser.add_argument('--scale', type=int, default=24)
    parser.add_argument('--memory-mib', type=int, default=4096)
    args = parser.parse_args()
    assert args.repeat > 0 and args.scale >= 6
    resource.setrlimit(resource.RLIMIT_AS, (args.memory_mib * 1024**2,) * 2)
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    import qpa
    import qpa._qpa as extension

    def loop_case(field):
        algebra = qpa.QuotientAlgebra.truncated(field, qpa.Quiver.indexed(1, [(1, 1, 'x')]), 3)
        simple = algebra.simple_modules()[0]
        projective = algebra.indecomposable_projectives()[0]
        length_two = algebra.module_from_integer_matrices([2], [('x', [[0, 1], [0, 0]])])
        k = args.scale // 6
        direct_sum = qpa.MatrixModule.direct_sum
        source = direct_sum([simple] * (2*k) + [length_two] * k + [projective] * k)
        targets = [direct_sum([simple] * k), direct_sum([length_two] * (2*k)),
                   direct_sum([projective] * k + [simple] * k), direct_sum([projective] * k)]
        return source, targets, [3*k*k, 6*k*k, 3*k*k, 0]

    def kronecker_case(field):
        n = args.scale
        parent = qpa.PathAlgebra(field, qpa.Quiver.indexed(2, [(1, 2, 'a'), (1, 2, 'b')]))
        zero, one = field.element(0), field.element(1)
        eigenvalue = field.primitive_element() if (field.degree or 1) > 1 else zero
        identity = [[one if i == j else zero for j in range(n)] for i in range(n)]
        def regular(value):
            jordan = [[value if i == j else one if j == i+1 else zero
                       for j in range(n)] for i in range(n)]
            return parent.module_from_matrices([n, n], [('a', identity), ('b', jordan)])
        source = regular(eigenvalue)
        return source, [source, regular(one)], [n, 0]

    report = dict(command=sys.argv, scale=args.scale, repeat=args.repeat,
                  version=version('qpa-rust'), revision=qpa.PORT_GIT_COMMIT,
                  runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  extension_sha256=hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest(),
                  cases=[], status='running')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    def save():
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    save()
    for name, field in [('Q', qpa.Field.rational()), ('GF2', qpa.Field.prime(2)),
                        ('GF101', qpa.Field.prime(101)), ('GF4', qpa.Field.extension(2, 2))]:
        for family, build in [('truncated_loop_mixed', loop_case), ('regular_kronecker', kronecker_case)]:
            start = time.perf_counter()
            source, targets, expected = build(field)
            construction = time.perf_counter() - start
            start = time.perf_counter()
            cover = source.projective_cover()
            inclusion = cover.kernel_inclusion()
            preparation = time.perf_counter() - start
            # Keep both maps alive during every query, without retaining a
            # whole projective family or installing additional algebra caches.
            assert cover.is_surjective and inclusion.is_injective
            samples = []
            for _ in range(args.repeat):
                start = time.perf_counter()
                actual = source.ext1_dimensions(targets)
                samples.append(time.perf_counter() - start)
                assert actual == expected, (name, family, actual, expected)
            case = dict(field=name, family=family, dimensions=actual,
                        source_dimensions=source.dimension_vector,
                        target_dimensions=[target.dimension_vector for target in targets],
                        construction_seconds=construction, preparation_seconds=preparation,
                        query_seconds=samples, query_median_seconds=statistics.median(samples),
                        process_peak_rss_kib_so_far=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            report['cases'].append(case)
            save()
            print(json.dumps(case), flush=True)
            del source, targets, cover, inclusion
    report['status'] = 'complete'
    save()


if __name__ == '__main__':
    main()
