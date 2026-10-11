#!/usr/bin/env python3
"""Bounded, synthetic presentation-to-basis benchmarks (never imports consumers).

Run with the release Python environment, e.g.:
  .venv/bin/python benchmarks/algebra_construction.py --output benchmarks/results/algebra-construction.json
Each sample is a fresh process, with a wall timeout and address-space limit.
Stage timings exclude interpreter/package startup; process wall time includes it.
Only standard-library modules are needed by the coordinator and GAP backend.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import resource
import signal
import statistics
import subprocess
import sys
import tempfile
import time


CASES = [f"monomial-{m}-{d}" for m, d in [(2, 3), (2, 4), (3, 3), (3, 4)]] + [
    f"{family}-{n}" for n in [2, 3, 4] for family in ["quantum", "dense"]
] + ["grid-2", "grid-3"]


def presentation(name):
    """Return integer coefficient/path lists shared verbatim by both engines."""
    parts = name.split("-")
    family, n = parts[0], int(parts[1])
    if family == "grid":
        vertices = n * n
        arrows, lookup = [], {}
        for i in range(n):
            for j in range(n):
                for di, dj in [(1, 0), (0, 1)]:
                    if i + di < n and j + dj < n:
                        source, target = i * n + j, (i + di) * n + j + dj
                        lookup[source, target] = len(arrows)
                        arrows.append((source + 1, target + 1, f"a{len(arrows)}"))
        relations = []
        for i in range(n - 1):
            for j in range(n - 1):
                a, b, c, d = i*n+j, (i+1)*n+j, i*n+j+1, (i+1)*n+j+1
                relations.append([(1, [lookup[a, b], lookup[b, d]]),
                                  (-1, [lookup[a, c], lookup[c, d]])])
        dimension = (n * (n + 1) // 2) ** 2
    else:
        vertices = 1
        arrows = [(1, 1, f"x{i}") for i in range(n)]
        if family == "monomial":
            degree = int(parts[2])
            relations = [[(1, list(word))] for word in itertools.product(range(n), repeat=degree)]
            dimension = sum(n ** i for i in range(degree))
        else:
            # y_i = x_i (quantum) or x_i + sum_j x_j (dense).
            # I + all-ones has determinant n+1, invertible over Q and GF(101).
            change = [[int(i == j) + int(family == "dense") for j in range(n)] for i in range(n)]
            def product(i, j, factor=1):
                return {(a, b): factor * change[i][a] * change[j][b]
                        for a in range(n) for b in range(n)}
            polynomials = [product(i, i) for i in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    p, q = product(j, i), product(i, j, -2)
                    polynomials.append({word: p[word] + q[word] for word in p})
            relations = [[(coefficient, list(word)) for word, coefficient in sorted(p.items()) if coefficient]
                         for p in polynomials]
            dimension = 2 ** n
    return {"vertices": vertices, "arrows": arrows, "relations": relations, "dimension": dimension}


def emit(event, **data):
    print("BENCH " + json.dumps({"event": event, **data}), flush=True)


def rust_worker(name, field, route):
    import qpa
    import qpa._qpa as extension
    emit("metadata", revision=qpa.PORT_GIT_COMMIT, package=qpa.__file__,
         extension_sha256=hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest())
    def timed(stage, operation):
        emit("start", stage=stage)
        start = time.perf_counter()
        value = operation()
        emit("timing", stage=stage, seconds=time.perf_counter() - start)
        return value
    spec = presentation(name)
    def setup():
        k = qpa.Field.rational() if field == "Q" else qpa.Field.prime(101)
        algebra = qpa.PathAlgebra(k, qpa.Quiver.indexed(spec["vertices"], spec["arrows"]))
        generators = [algebra.generator(a[2]) for a in spec["arrows"]]
        relations = []
        for relation in spec["relations"]:
            value = algebra.zero()
            for coefficient, word in relation:
                term = generators[word[0]]
                for index in word[1:]:
                    term = term * generators[index]
                value = value + term.scale_left(k.element(coefficient))
            relations.append(value)
        return algebra, relations
    algebra, relations = timed("presentation", setup)
    if route == "direct":
        quotient = timed("quotient_with_gb", lambda: algebra.quotient(relations))
        gb = quotient.groebner_basis
        count = len(gb.relations)
    else:
        # The ideal overload returns the completed relation list without wrapping
        # an extra GroebnerBasis, so attachment is measured separately.
        ideal = timed("ideal", lambda: algebra.ideal(relations))
        complete = timed("groebner", ideal.gbnp_groebner_basis)
        count = len(complete)
        # The direct relation-list constructor replaces the defining generators
        # by the completed list, then also checks admissibility. Mirror both.
        ideal = timed("completed_ideal", lambda: algebra.ideal(complete))
        timed("attach", lambda: ideal.with_complete_groebner_basis(complete))
        quotient = timed("quotient", ideal.quotient)
        assert timed("admissibility", ideal.is_admissible)
    dimension = int(timed("dimension", quotient.dimension))
    basis = timed("basis", quotient.canonical_basis)
    assert dimension == len(basis) == spec["dimension"], (name, dimension, len(basis))
    emit("result", dimension=dimension, basis_size=len(basis), gb_relations=count)


def gap_program(name, field, route):
    spec = presentation(name)
    arrows = ",".join(f'[{a},{b},"{label}"]' for a, b, label in spec["arrows"])
    relations = []
    for relation in spec["relations"]:
        relations.append("+".join(f"({c}*One(k))*" + "*".join(f"g[{i+1}]" for i in word)
                                  for c, word in relation))
    def timed(stage, code):
        return f'''Print("BENCH {{\\"event\\":\\"start\\",\\"stage\\":\\"{stage}\\"}}\\n");
t := NanosecondsSinceEpoch();;
{code}
Print("BENCH {{\\"event\\":\\"timing\\",\\"stage\\":\\"{stage}\\",\\"nanoseconds\\":", NanosecondsSinceEpoch()-t, "}}\\n");
'''
    program = '''if LoadPackage("qpa") <> true then FORCE_QUIT_GAP(1); fi;
SetPrintFormattingStatus("*stdout*", false);;
Print("BENCH {\\"event\\":\\"metadata\\",\\"gap\\":\\"", GAPInfo.Version,
      "\\",\\"qpa\\":\\"", PackageInfo("qpa")[1].Version,
      "\\",\\"gbnp\\":\\"", PackageInfo("gbnp")[1].Version, "\\"}\\n");
'''
    program += timed("presentation", f'''k := {"Rationals" if field == "Q" else "GF(101)"};;
Q := Quiver({spec['vertices']}, [{arrows}]);;
A := PathAlgebra(k,Q);;
g := List([{','.join(json.dumps(a[2]) for a in spec['arrows'])}], s -> A.(s));;
rels := [{','.join(relations)}];;''')
    if route == "direct":
        program += timed("quotient_with_gb", "B := A/rels;;")
        # Fetch below, after dimension has guaranteed the basis is available.
    else:
        program += timed("ideal", "I := Ideal(A,rels);;")
        program += timed("groebner", "complete := GBNPGroebnerBasis(I);;")
        program += timed("completed_ideal", "I := Ideal(A,complete);;")
        program += timed("attach", "gb := GroebnerBasis(I,complete);;")
        program += timed("quotient", "B := A/I;;")
        program += timed("admissibility", "admissible := IsAdmissibleIdeal(I);;")
        program += "if not admissible then FORCE_QUIT_GAP(3); fi;\n"
    program += timed("dimension", "dim := Dimension(B);;")
    program += timed("basis", "bas := BasisVectors(Basis(B));;")
    if route == "direct":
        program += "complete := AsList(GroebnerBasisOfIdeal(IdealOfQuotient(B)));;\n"
    program += f'''if dim <> {spec['dimension']} or Length(bas) <> dim then FORCE_QUIT_GAP(2); fi;
Print("BENCH {{\\"event\\":\\"result\\",\\"dimension\\":", dim,
      ",\\"basis_size\\":", Length(bas), ",\\"gb_relations\\":", Length(complete), "}}\\n");
QUIT;
'''
    return program


def run_sample(args, case, field, backend, route):
    with tempfile.TemporaryDirectory(prefix="qpa-construction-") as directory:
        directory = Path(directory)
        if backend == "rust":
            command = [sys.executable, str(Path(__file__).resolve()), "--worker", case, field, route]
        else:
            script = directory / "case.g"
            script.write_text(gap_program(case, field, route))
            command = [args.gap, "-q", "-b", "-T", "-m", "64m", "-o", "512m", "-K", "768m", str(script)]
        rss_file = directory / "rss"
        command = ["/usr/bin/time", "-f", "%M", "-o", str(rss_file), *command]
        def limit():
            cap = args.memory_mib * 1024 * 1024
            resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
        start = time.perf_counter()
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   text=True, start_new_session=True, preexec_fn=limit)
        status = "ok"
        try:
            stdout, stderr = process.communicate(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            status = "timeout"
            os.killpg(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
        wall = time.perf_counter() - start
        events = [json.loads(line[6:]) for line in stdout.splitlines() if line.startswith("BENCH ")]
        timings = {e["stage"]: e.get("seconds", e.get("nanoseconds", 0) / 1e9)
                   for e in events if e["event"] == "timing"}
        result = next((e for e in events if e["event"] == "result"), None)
        if status == "ok" and (process.returncode != 0 or result is None):
            status = "error"
        rss_lines = rss_file.read_text().splitlines() if rss_file.exists() else []
        sample = {"status": status, "returncode": process.returncode, "wall_seconds": wall,
                  "stage_seconds": timings, "result": result,
                  "metadata": next((e for e in events if e["event"] == "metadata"), None),
                  "peak_rss_kib": int(rss_lines[-1]) if rss_lines and rss_lines[-1].isdigit() else None,
                  "last_started_stage": next((e["stage"] for e in reversed(events) if e["event"] == "start"), None)}
        if status != "ok" or stderr:
            sample.update(stdout=stdout, stderr=stderr)
        return sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", nargs=3, metavar=("CASE", "FIELD", "ROUTE"), help=argparse.SUPPRESS)
    parser.add_argument("--case", action="append", choices=CASES)
    parser.add_argument("--field", action="append", choices=["Q", "GF101"])
    parser.add_argument("--backend", action="append", choices=["rust", "gap"])
    parser.add_argument("--route", action="append", choices=["staged", "direct"])
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=10)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--gap", default="gap")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.worker:
        rust_worker(*args.worker)
        return
    if not args.output or args.repeat < 1 or args.timeout <= 0 or args.memory_mib < 1:
        parser.error("supply --output and positive repeat/timeout/memory limits")
    cases = args.case or CASES
    specs = {case: presentation(case) for case in cases}
    report = {"platform": platform.platform(), "python": sys.version,
              "cpu_affinity": sorted(os.sched_getaffinity(0)),
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "timeout_seconds": args.timeout, "address_space_limit_mib": args.memory_mib,
              "repeat": args.repeat, "presentations": specs, "runs": {}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for case, field, backend, route in itertools.product(cases, args.field or ["Q", "GF101"],
                                                        args.backend or ["rust", "gap"],
                                                        args.route or ["staged"]):
        key = "/".join([case, field, backend, route])
        samples = []
        for _ in range(args.repeat):
            sample = run_sample(args, case, field, backend, route)
            samples.append(sample)
            if sample["status"] != "ok":
                break  # Do not repeat a timeout or failure.
        entry = {"samples": samples}
        if all(s["status"] == "ok" for s in samples):
            assert all(s["result"] == samples[0]["result"] for s in samples), key
            entry["median_stage_seconds"] = {stage: statistics.median(s["stage_seconds"][stage] for s in samples)
                                              for stage in samples[0]["stage_seconds"]}
            entry["median_pipeline_seconds"] = statistics.median(sum(s["stage_seconds"].values()) for s in samples)
        report["runs"][key] = entry
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(key, samples[-1]["status"], entry.get("median_pipeline_seconds", ""), flush=True)
    # Cross-engine/route verification of invariant dimension and basis count.
    # GB counts are retained for diagnosis, not assumed equal across engines.
    for case, field in itertools.product(cases, args.field or ["Q", "GF101"]):
        results = [s["result"] for key, entry in report["runs"].items()
                   if key.startswith(f"{case}/{field}/") for s in entry["samples"] if s["status"] == "ok"]
        assert all(r["dimension"] == r["basis_size"] == specs[case]["dimension"] for r in results)


if __name__ == "__main__":
    main()
