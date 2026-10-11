# Native performance work before the next version bump

Scope agreed with the user: complete **all four priorities**, beginning with
priority 1. All four are now implemented and validated. Version/release changes
are left to the maintainer. All optimizations
must work through ordinary QPA APIs, without algebra-specific names, dimensions,
manual relation partitions, or supplied checkpoints.

## 1. Genuinely dimension-only Ext — completed

The previous `ExtDimensionSource::dimension` constructed `Hom(P,N)` basis maps,
composed their restrictions, and built a dense image-rank matrix even though
only a dimension was requested. This work has been replaced using the exact sequence:

```
0 -> Hom(M,N) -> Hom(P,N) -> Hom(Omega M,N) -> Ext^1(M,N) -> 0
```

Thus `dim Ext^1 = dim Hom(Omega M,N) - dim Hom(P,N) + dim Hom(M,N)`.
For the native cover `P = direct_sum_i (e_i A)^m_i`, compute
`dim Hom(P,N) = sum_i m_i dim(N_i)` from the cover's generator multiplicities.
Keep the basis-returning Ext API and its ordered representatives unchanged.

Completed checks:
- [x] Dimension-only evaluation does not build projective Hom bases, restriction
  maps, or a restriction matrix.
- [x] Preserve existing algebra/field/domain checks, zero cases, and batch order.
- [x] Check integer arithmetic; retain short-circuiting when `Hom(Omega M,N)=0`.
- [x] Compare against the unchanged basis-returning API over Q, prime fields and
  extension fields, including nonzero restriction images and repeated summands.
- [x] Run core/GBNP/Python regressions and record ordinary-API time and peak memory
  before/after on multiple workloads, including the existing large regression.

Evidence: [priority-1 implementation and benchmark report](../../benchmarks/results/ext-dimension-only.md).
439 Rust tests, 603 Python tests and strict Clippy pass. Eight synthetic query
cases agree before/after. The large ordinary-API run drops from 540.65 to
467.40 seconds and from 4.36 to 2.73 GiB peak RSS; the degree-two dimension call
drops from 142.659 to 50.118 seconds. The basis-returning API is unchanged.
Checked arithmetic guards overflow/inconsistent dimensions; no public function
signature or Python exception class changes. The development-only error variants
were replaced with existing 0.2.0 errors during [0.2.1 release
preparation](release-0.2.1.md), preserving downstream exhaustive Rust matches.

## 2. Projective-cover construction — completed

- [x] Profile the previous cover stages: full-family construction costs
  116.017 / 114.388 seconds; direct sums cost 7.074 / 61.448 seconds.
- [x] Construct only selected projectives and reuse individual live weak-cache
  members when siblings expire; retain identity and cycle-safe lifetimes.
- [x] Extract product coordinates directly into local vertex blocks using lazy
  per-projective word indexes, rather than full-algebra vectors.
- [x] Assemble cover maps as ordered row blocks instead of projection composites.
- [x] Make direct-sum component maps lazy without losing live-map identity.
- [x] Compare matrices with the previous constructions across fields, including
  nonprime coefficients, repeated summands, and weak-cache expiry.

Ordinary cover calls fall from 132.131 / 185.108 to 4.631 / 41.288 seconds in the
large regression. See the [complete pass report](../../benchmarks/results/native-performance.md).

## 3. Module relation checking — completed

- [x] Reuse common path prefixes within bounded, local validation caches.
- [x] Borrow input actions, avoid zero-product allocations, and accumulate scaled
  matrices in place instead of making scaled and added copies.
- [x] Reuse validation only for structurally certified uniform direct sums;
  retain full checks for entered actions, projectives, nonuniform relations,
  and foreign-runtime-field NC cases.
- [x] Compare with the old evaluator across fields, zero spaces, and failures.

No public validation-skip flag or uncertified-presentation shortcut was added.
The diagnostic reference direct sums fall to 0.114 / 0.898 seconds; standalone
checked loop modules also improve in all four measured fields.

## 4. Elimination after sparse fill-in — completed

- [x] Measure fill-in and promote individual rows to dense storage under a
  shared slot budget, without allocating a large Cartesian-product matrix.
- [x] Pack sufficiently dense GF(2) inputs directly from sparse storage, with
  memory and density guards; do not confuse GF(4) with GF(2).
- [x] Avoid exporting unused basis copies for dense rank-only queries; retain
  transformation-free sparse rank evaluation.
- [x] Compare exact ordered nullspaces and ranks across Q, GF(2), GF(101), GF(4),
  including late/early pivots, dependent rows, and guard boundaries.
- [x] Measure dense and sparse workloads; disclose remaining sparse-query
  overhead instead of claiming every operation became faster.

Dense Hom cases improve across all four fields. A sparse-GF(2) packing regression
was mitigated by a density guard; residual small-query regressions, a cancelled
rational coefficient-growth experiment, and policy limits are documented in
[the report](../../benchmarks/results/native-performance.md).

## Baseline and release gates

The original ordinary-API baseline is
[`arc-native-api.md`](../../benchmarks/results/arc-native-api.md):
540.65 seconds and 4,566,792 KiB peak RSS. After priority 1: 467.40 seconds and
2,860,504 KiB. After all four: **180.29 seconds and 2,861,104 KiB**.
The global memory reduction is from priority 1; priorities 2–4 substantially
reduce runtime and earlier-stage memory, but do not lower the final global peak.

Completed release gates:
- [x] 446 Rust unit tests, 4 integration tests, 2 doctests, 611 Python tests;
  strict Clippy, formatting and syntax checks.
- [x] Frozen QPA fixtures and exact/ordered reference comparisons preserved.
- [x] 24 multi-field general benchmark cases and 8 dimension-only Ext cases
  agree before/after; raw samples and time/memory measurements retained.
- [x] Large ordinary-API results and both short exact sequences agree; original
  inputs remain unchanged and no checkpoint or special optimization flag is used.
- [x] Final binary/input/runner hashes checked, separate call and end-to-end
  timings reported, and regressions/limits disclosed.
- [x] No version bump or release change performed.

Authoritative evidence: [native-performance.md](../../benchmarks/results/native-performance.md).
