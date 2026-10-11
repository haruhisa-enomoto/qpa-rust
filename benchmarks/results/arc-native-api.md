# Ordinary QPA API: generic optimization, no example-specific preconditioner

Follow-ups: [dimension-only Ext](ext-dimension-only.md) and the
[complete four-priority performance pass](native-performance.md) reduce the
same ordinary-API workload further. The timings below are the baseline **before**
those follow-ups; the generic-constructor workflow remains valid.

This supersedes the caller-side construction workflow in
[the earlier optimization report](arc-optimization.md). The goal here is QPA
engine performance, using the ARC input only as a regression workload.

## What changed in the native library

The normal checked GBNP adapter now performs **automatic structural
preconditioning**:

1. Uniformize the input relations as before.
2. Compute strongly connected components of the quiver with an iterative,
   linear-time graph traversal.
3. If both kinds are present, complete the relations internal to components
   first, then reduce relations crossing components by that completed ideal.
4. Run the existing GBNP completion and final reduction on their union.

There are no algebra names, distinguished vertex labels, prescribed relation
counts, manual partitions, or expected dimensions in this optimization.
It applies to arbitrary quivers and supported coefficient fields. A path whose
endpoints belong to one strongly connected component cannot leave that
component and return. Crossing relations cannot contribute internal-corner
relations, so this is a genuine structural separation. Replacing each crossing
relation by its remainder preserves the combined ideal. All reductions stay
in GBNP's NP coordinates/order, and the final ordered reduced basis is checked
against the unpreconditioned algorithm in regression tests.

When there is no internal/crossing split, the existing completion path is used.
The generic free-algebra `SGrobner` and unchecked adapter are unchanged.

`PathAlgebra.quotient` reaches this through the normal Rust
`qpa_gbnp::quotient_by_relations` / `gbnp_groebner_basis` call chain, not through
a benchmark hook.

Implementation:

- `crates/qpa-gbnp/src/qpa_conversion.rs`
- `crates/qpa-gbnp/src/qpa_conversion/components.rs`

The earlier native sparse admissibility, sparse Hom, subrepresentation-closure,
and vector-action improvements remain in place. Those also contain no
ARC-specific rules.

## Normal user workflow

For QPA right-module matrices, ordinary calls suffice:

```python
import qpa

k = qpa.Field.rational()
R = qpa.PathAlgebra(k, quiver)
A = R.quotient(relations)  # the original, unpartitioned defining relations
M = A.module_from_matrices(dimension_vector, arrow_matrices)

cover = M.projective_cover()
kernel = cover.kernel_inclusion()
ext1_dimension = M.ext1_dimension(M)
ext2_dimension = kernel.source.ext1_dimension(M)
```

The last line uses the usual dimension shift through the first syzygy.
No supplied Gröbner basis, script-level reduction, or optimization flag is
needed. Retaining the cover/kernel in variables permits reuse of their native
weak caches; it does not change any mathematics.

The recorded example supplies **left** column-action matrices. Its replay
therefore uses the opposite presentation and transposed matrices, as required
by QPA's right-module convention. This is an input-convention conversion,
not a performance workaround: the ordinary constructor is tested on the
published orientation too.

## Measurements on the ordinary API

Same Linux/WSL2 machine and release configuration as the earlier reports;
single fresh-process observations, not medians. No profiler or supplied basis.

| Call | Observation |
|---|---:|
| Published `R.quotient(original_relations)`, before automatic component completion | stopped at 120 s |
| Published `R.quotient(original_relations)`, after | 30.244 s |
| Opposite `R.quotient(original_opposite_relations)`, after | 28.320 s |

The ordinary quotient call includes native completion and eager admissibility.
The separate later admissibility query is therefore a cache hit, not a new
near-zero-time admissibility algorithm. Both quotients have dimension and
canonical basis size **7708**; every input relation reduces to zero.

The ordinary-API Ext replay is recorded incrementally in
`arc-ext-ordinary.json` with `construction: ordinary_api`. Its inputs contain
only the presentation and witness, **not a completed-basis checkpoint**.

**Ext replay completed successfully:** both selected self-Ext dimensions
(degrees 1 and 2) are zero. Both projective-cover short exact sequences passed
the native checks. This is an engine regression result, not an all-degree
mathematical claim.

The entire ordinary-API replay, including fresh opposite construction, took
**540.65 seconds (9 minutes)** and **4,566,792 KiB (4.36 GiB) peak RSS**, within
the 900-second / 8-GiB guards. In particular, the projective-cover calls took
120.952 and 179.084 seconds; the two Ext dimension calls took 10.906 and 142.659
seconds. These costs are measured, not substituted with the earlier faster
checkpoint-based timings.

The separate published-constructor worker took **30.59 seconds end-to-end**,
with **96,416 KiB (94.16 MiB) peak RSS**.

Raw records:

- `arc-ordinary-published.stdout.log` / `.stderr.log`: 120-second baseline on
  the previously optimized engine, before automatic native component completion;
- `arc-native-published.stdout.log` / `.stderr.log`: ordinary published
  constructor after the generic native change;
- `arc-ext-ordinary.json`, `.stdout.log`, `.stderr.log`: ordinary opposite
  construction followed by native module/Ext operations.

The earlier staged Ext time must not be presented as the ordinary-API time.
In particular, the ordinary relation-list constructor retains the completed
relations as its defining list, following QPA, so subsequent module relation
checks can be more expensive than in the old original-ideal/checkpoint replay.

## Validation

New differential tests compare the **ordered full basis**, not just quotient
dimension, against the previous checked-adapter algorithm. They exercise
vertex-number swaps, arrow and relation permutations, nonuniform and
nonhomogeneous relations, and characteristic-dependent reductions over Q,
GF(2), GF(7), and GF(101). Graph tests cover nontrivial strongly connected
components, one-way connections, isolated vertices, and no-split fallbacks.

New Python tests use renamed vertices/arrows and interleaved internal/crossing
relations through `PathAlgebra.quotient`, with no caller-supplied partition.
They check dimension 27, all relation remainders, admissibility and nonzero Ext
between simples over Q, GF(2), GF(5), and GF(4), including a nonprime coefficient.
The existing frozen QPA fixtures also run unchanged.

Final validation passed: **224 core + 213 GBNP tests**, strict Clippy, and
**601 Python tests**. The ordered-basis differential check covers 40 combinations
of fields, vertex swaps, and arrow/relation permutations. Formatting, syntax,
and final artifact-consistency checks also pass.

Final validation commands:

```bash
cargo test -p qpa-core -p qpa-gbnp --lib
cargo clippy -p qpa-core -p qpa-gbnp --all-targets -- -D warnings
.venv/bin/python -m pytest python/tests -q
```

## Reproduction

From the repository root, with the read-only input documents listed in the
[earlier report](arc-optimization.md) available:

```bash
VIRTUAL_ENV="$PWD/.venv" .venv/bin/maturin develop --release

/usr/bin/time -v timeout -k 2 180 .venv/bin/python benchmarks/arc_materialize.py \
  --admissibility

/usr/bin/time -v timeout -k 2 900 .venv/bin/python benchmarks/arc_ext.py \
  --degree 2 --output benchmarks/results/arc-ext-ordinary-reproduced.json
```

The manual `--precondition` switch has been removed from the construction
runner. The Ext runner defaults to ordinary construction; its optional
`--basis` argument exists only to replay historical checkpoint measurements.
Input-size/dimension assertions in these benchmark loaders are regression
checks, not inputs to a specialized engine algorithm.

The recorded native extension SHA-256 is
`10c34ab8e5793782b184a6345d03e16e56b5367a1fe9b868327b48853a8e53e2`.
It is a modified working-tree release build, not an unchanged released binary.
Gröbner completion can still be expensive or nonterminating on general inputs;
this is not a claim that all quiver algebras have cheap Ext computations.
