#!/usr/bin/env python3
"""Benchmark exact module calculations with fresh modules on every repetition.

Build with `.venv/bin/maturin develop --release`, then run:
    .venv/bin/python benchmarks/module_workloads.py --output target/baseline.json

Optional consumer fixtures are loaded from --consumer-root without running their
command-line entry points or writing research certificates. Use the same Python,
build profile, backend, cases and repetitions for before/after comparisons. The
JSON records every sample, result, extension hash and loaded package location.
Timings include construction and validation but exclude imports and JSON output.
Consumer fixtures retain their existing presentation caches between repetitions;
these are repeated in-process workloads, not fresh-process cold-start timings.
"""

import argparse
import hashlib
import json
import platform
import runpy
import statistics
import sys
import time
from pathlib import Path

import qpa
import qpa._qpa as extension


def kronecker(field, size):
    algebra = qpa.PathAlgebra(field, qpa.Quiver.indexed(2, [(1, 2, "a"), (1, 2, "b")]))
    identity = [[int(i == j) for j in range(size)] for i in range(size)]
    jordan = [[int(j == i + 1) for j in range(size)] for i in range(size)]
    module = algebra.module_from_integer_matrices([size, size], [("a", identity), ("b", jordan)])
    return module


def hom_rational():
    module = kronecker(qpa.Field.rational(), 24)
    dimension = len(module.hom_basis_to(module))
    assert dimension == 24
    return dimension


def kernel_rational():
    size, shift = 48, 16
    module = kronecker(qpa.Field.rational(), size)
    zero, one = module.field.element(0), module.field.element(1)
    matrix = [[one if j == i + shift else zero for j in range(size)] for i in range(size)]
    kernel = module.homomorphism_to(module, [matrix, matrix]).kernel()
    assert kernel.dimension_vector == [shift, shift]
    return kernel.dimension_vector


def ext_binary():
    algebra = qpa.QuotientAlgebra.truncated(
        qpa.Field.prime(2), qpa.Quiver.indexed(1, [(1, 1, "x"), (1, 1, "y")]), 2
    )
    module = algebra.module_from_integer_matrices([8], [])
    dimension = module.ext1(module).dimension
    assert dimension == 128
    return dimension


def additive_closure_prime():
    algebra = qpa.PathAlgebra(qpa.Field.prime(5), qpa.Quiver.indexed(1, []))
    candidate = algebra.module_from_integer_matrices([6], [])
    generator = algebra.module_from_integer_matrices([7], [])
    result = candidate.is_in_additive_closure_of(generator)
    assert result
    return result


def consumer_cases(root):
    sys.path.insert(0, str(root))

    def load(path):
        return runpy.run_path(str(root / path), run_name="benchmark_fixture")

    js_path = "experiments/weak-gorenstein-symmetry-counterexample/verify_js29.py"
    js = load(js_path)
    certificate = json.loads((root / js_path).with_name("js29-input.json").read_text())
    graded = load("tests/test_graded_square_zero_qpa.py")
    delooping = load("tests/test_delooping_qpa.py")
    tubes = load("experiments/reflexive-categories/verify_triangular_tube_lattice.py")
    e6 = load("experiments/brick-semibrick-bicompact/verify_e6_faithful_contraction.py")
    gsz = load("tests/test_bounded_global_dimension_qpa.py")

    def gsz_nonminimal():
        result = gsz["bounded_right_global_dimension_gsz"](
            gsz["dimension_sixteen_p9_nonminimal"](), bound=6
        )
        assert result.status == "exceeds_bound"
        return result.status

    return {
        "consumer_js29_rational": lambda: js["verify"](certificate),
        "consumer_graded_ext_binary": graded["test_middle_arrow_complexes_and_negative_self_hom"],
        "consumer_delooping": lambda: delooping[
            "test_second_syzygy_summand_need_not_be_a_third_syzygy_summand"
        ]("transpose"),
        "consumer_reflexive_tubes": lambda: tubes["calculate"](tubes["inputs"]()),
        "consumer_e6_rational": lambda: e6["verify"](qpa.Field.rational(), "Q"),
        "consumer_gsz_nonminimal": gsz_nonminimal,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--case", action="append", help="Run only these named cases")
    parser.add_argument("--consumer-root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compare", type=Path, help="Check results against an earlier JSON run")
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error("--repeat must be positive")
    cases = {
        "hom_rational": hom_rational,
        "kernel_rational": kernel_rational,
        "ext_binary": ext_binary,
        "additive_closure_prime": additive_closure_prime,
    }
    if args.consumer_root:
        cases.update(consumer_cases(args.consumer_root.resolve()))
    selected = args.case or list(cases)
    if unknown := set(selected) - cases.keys():
        parser.error(f"unknown cases: {sorted(unknown)}")
    baseline = json.loads(args.compare.read_text()) if args.compare else None
    if baseline and not set(selected) <= baseline["cases"].keys():
        parser.error("comparison run is missing a selected case")
    report = {
        "python": sys.version,
        "platform": platform.platform(),
        "qpa_package": qpa.__file__,
        "qpa_build_revision": qpa.PORT_GIT_COMMIT,
        "extension_sha256": hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest(),
        "repetitions": args.repeat,
        "cases": {},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for name in selected:
        samples = []
        results = []
        for _ in range(args.repeat):
            started = time.perf_counter()
            result = cases[name]()
            samples.append(time.perf_counter() - started)
            results.append(json.loads(json.dumps(result)))
        assert all(result == results[0] for result in results), f"unstable result: {name}"
        report["cases"][name] = {
            "seconds": samples, "median_seconds": statistics.median(samples), "result": results[0]
        }
        if baseline:
            assert results[0] == baseline["cases"][name]["result"], f"changed result: {name}"
            report["cases"][name]["speedup"] = (
                baseline["cases"][name]["median_seconds"] / statistics.median(samples)
            )
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(f"{name}: {statistics.median(samples):.6f}s", flush=True)


if __name__ == "__main__":
    main()
