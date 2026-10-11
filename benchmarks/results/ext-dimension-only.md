# Priority 1: genuinely dimension-only Ext

Follow-up: [the completed four-priority pass](native-performance.md) implements
priorities 2–4 and records their additional gains and remaining tradeoffs. This
report retains the historical priority-1 measurements.

This implements priority 1 of the
[native performance plan](../../docs/development/performance-plan.md).
Priorities 2–4 (projective construction, relation checking, adaptive elimination)
are recorded but **not implemented in this pass**. No version number was changed.

## Implementation and exactness

`ExtDimensionSource::dimension` no longer constructs `Hom(P,N)` basis maps,
restriction morphisms, or a matrix of their images. For the native projective
cover sequence `0 -> Omega M -> P -> M -> 0`, exactness gives

```
image_rank = dim Hom(P,N) - dim Hom(M,N)
dim Ext^1(M,N) = dim Hom(Omega M,N) - image_rank
```

The native cover has one `e_i A` summand for each vertex-uniform minimal
generator. These generators are already cached by cover construction; their
vertex multiplicities are counted once per prepared source. Evaluation at the
idempotent identifies `Hom(e_i A,N)` with `N_i`, so

```
dim Hom(P,N) = sum_i multiplicity_i * dim(N_i).
```

Only `Hom(Omega M,N)` and `Hom(M,N)` need dimension/rank computations. No field
approximations or algebra-specific information enter the calculation. The
basis-returning `ext1` / `ext_over_algebra` routine is unchanged, including its
ordered representatives and coefficient closure.

The zero-syzygy and zero-syzygy-Hom cases retain their fast exits. The implementation
does not introduce a target-zero shortcut that would suppress unsupported
projective-cover errors. Products and sums are checked for address-size overflow;
subtractions check the exact-sequence inequalities. Subtracting the image rank
avoids an unnecessary intermediate sum that could overflow even when the answer
fits.

Files:

- `crates/qpa-gbnp/src/module/ext_and_iyama.rs`
- `crates/qpa-gbnp/src/module/ext_dimension_tests.rs`
- `crates/qpa-gbnp/src/module/errors.rs`
- `python/tests/test_module_dimensions.py`

Public function signatures are unchanged. The measured development builds added
`ExtDimensionOverflow` and `InconsistentExtDimensions`. During [0.2.1 release
preparation](../../docs/development/release-0.2.1.md), these unpublished variants
were replaced by existing `Module(DimensionOverflow)` and
`InternalExtHomRejected` errors, preserving downstream exhaustive Rust matches.
Python continues to use the existing exception mapping. The historical binary
hashes and measurements below still refer to those development builds.

## Synthetic before/after measurements

Runner: [`../ext_dimension_only.py`](../ext_dimension_only.py), scale 24,
three repeated queries per case in each fresh process. Cover/kernel preparation
is timed separately; their native maps are retained in ordinary local variables
so the query measurements exclude rebuilding them. Each query is a target batch.

The mixed-loop family uses modules over `k[x]/(x^3)` with simple, length-two,
and projective summands. Its expected output is `[48,96,48,0]`, including nonzero
restriction images. The Kronecker family uses size-24 regular modules in the
same and different eigenvalue families, with expected output `[24,0]`. The GF(4)
case includes a nonprime eigenvalue. Every output agreed before and after.

Median **query-only milliseconds** (three samples; not an end-to-end or universal
speedup claim):

| Field | Family | Before | After | Before / after |
|---|---|---:|---:|---:|
| Q | mixed truncated loop | 50.729 | 4.026 | 12.6x |
| Q | regular Kronecker | 292.108 | 1.601 | 182.5x |
| GF(2) | mixed truncated loop | 3.173 | 0.308 | 10.3x |
| GF(2) | regular Kronecker | 7.917 | 0.478 | 16.6x |
| GF(101) | mixed truncated loop | 3.848 | 0.386 | 10.0x |
| GF(101) | regular Kronecker | 12.241 | 0.555 | 22.1x |
| GF(4) | mixed truncated loop | 76.295 | 5.770 | 13.2x |
| GF(4) | regular Kronecker | 350.658 | 3.460 | 101.3x |

The complete eight-case worker, including preparation, startup and JSON output,
went from **2.54 s / 162,444 KiB peak RSS** to
**0.18 s / 53,124 KiB peak RSS**. These are single process-level observations.
The per-case RSS values in JSON are explicitly cumulative process high-water
marks, not independent case peaks. The workloads are structured examples;
dense fill-in performance remains part of priority 4, not a solved general claim.

Raw reports: `ext-dimension-only-before.json`, `ext-dimension-only-after.json`,
and their matching `.stdout.log` / `.stderr.log` files. Both runs used the same
runner SHA-256:

```
fe3b1baaff24479a6aa09b6f3a83ee05e394f127b4eeda61984fa5db63271e83
```

## Large ordinary-API regression

The before run is [`arc-ext-ordinary.json`](arc-ext-ordinary.json), documented in
[`arc-native-api.md`](arc-native-api.md). The after run starts from the same raw
presentation and witness, through the same ordinary quotient and module APIs,
without a supplied completed basis.

The after worker completed successfully. Both selected dimensions remain zero,
both native short exact sequences pass their checks, and the input hashes and
module/projective/syzygy dimension vectors agree with the before run.

| Measurement | Before | After |
|---|---:|---:|
| Degree-one Ext dimension call | 10.906 s | 4.906 s |
| Degree-two Ext dimension call | 142.659 s | 50.118 s |
| End-to-end ordinary-API worker | 540.65 s | 467.40 s |
| Peak RSS | 4,566,792 KiB (4.36 GiB) | 2,860,504 KiB (2.73 GiB) |

Thus the degree-two query is about **2.85x faster** and process peak RSS is
about **37% lower** in this comparison. Total wall time falls by about **14%**.
These are single large-process observations, not medians. The unchanged cover
calls took 132.131 and 185.108 seconds in the after run, versus 120.952 and
179.084 seconds before; they now account for about 68% of total runtime. No
improvement in cover construction is claimed. This confirms why priority 2
remains important rather than treating query-only speedups as whole-workflow
speedups.

After records: `arc-ext-dimension-only.json` and matching `.stdout.log` /
`.stderr.log`. The worker starts from raw input, uses no checkpoint or special
optimization flag, and stays within its 900-second / 8-GiB guards.

## Validation

Passed on the implementation:

- **224 qpa-core + 215 qpa-gbnp tests**;
- **603 Python tests**;
- strict Clippy for both Rust crates.

New Rust tests compare projective multiplicity counts with actual projective
Hom basis dimensions, and compare Ext dimensions with the unchanged basis API,
over Q, GF(2), GF(101), GF(4). The source has repeated summands at several
vertices: its multiplicities `[3,3,0]` differ from both its dimension vector
`[3,4,0]` and its cover's vector `[3,6,3]`. Arithmetic tests cover multiplication
and addition overflow, inconsistent ranks, and fitting answers whose naive
intermediate sum would overflow.

Existing grids and new Python boundaries cover zero modules, projectives,
repeated batch targets, nonzero restriction images, nonprime field coefficients,
different algebras/fields, infinite path algebras and nonadmissible quotients.
The source-zero cases still succeed while unsupported nonzero-source covers
still fail, even with a zero target.

## Reproduction and provenance

Both measurements use release builds of qpa-rust 0.2.0 at base revision
`6207cd941fefcadab5a695001be3aac573b1dacf-dirty`, on the same Linux/WSL2 machine.
They are different modified working-tree builds, **not the base commit alone**.
Extension SHA-256:

```
before 10c34ab8e5793782b184a6345d03e16e56b5367a1fe9b868327b48853a8e53e2
after  a84d9a1be30948cfd4b1fe7ec48141213ab38756aa7c1a44794cba98ea2430b9
```

From the repository root:

```bash
VIRTUAL_ENV="$PWD/.venv" .venv/bin/maturin develop --release

/usr/bin/time -v timeout -k 2 300 .venv/bin/python benchmarks/ext_dimension_only.py \
  --repeat 3 --scale 24 --output benchmarks/results/ext-dimension-only-reproduced.json

/usr/bin/time -v timeout -k 2 900 .venv/bin/python benchmarks/arc_ext.py \
  --degree 2 --output benchmarks/results/arc-ext-dimension-only-reproduced.json

cargo test -p qpa-core -p qpa-gbnp --lib
cargo clippy -p qpa-core -p qpa-gbnp --all-targets -- -D warnings
.venv/bin/python -m pytest python/tests -q
```

The synthetic runner is self-contained. The large runner uses the read-only
presentation and witness paths and hashes recorded in the previous ordinary-API
report. Reproducing **before** requires the pre-priority-1 implementation, not
just rerunning the command with the new library.
