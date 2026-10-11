# General native performance pass: all four priorities

This report follows the [performance plan](../../docs/development/performance-plan.md).
The dimension-only Ext work is documented separately in
[ext-dimension-only.md](ext-dimension-only.md). This pass implements priorities
2–4 while retaining that work. No version number or release was changed.

**All four priorities are implemented and validated.** The large ordinary-API
workflow takes **180.29 seconds**, versus **467.40 seconds** after priority 1
and **540.65 seconds** before these four priorities. Small sparse-query
regressions and other limits are disclosed below.

## 2. Projective construction and assembly

The reference-cover profile separated the previous public-API algorithm into
minimal generators, full projective family, generator maps, direct sum,
projections, projection composites, and addition:

| Reference stage | First, before | First, after | Second, before | Second, after |
|---|---:|---:|---:|---:|
| Full projective family | 116.017 s | 72.894 s | 114.388 s | 66.648 s |
| Generator maps | 0.034 s | 0.031 s | 0.389 s | 0.405 s |
| Direct sum | 7.074 s | 0.114 s | 61.448 s | 0.898 s |
| Projections | 0.105 s | 0.109 s | 0.664 s | 0.679 s |
| Projection composites | 0.003 s | 0.003 s | 0.041 s | 0.040 s |
| Adding composites | 0.016 s | 0.015 s | 0.250 s | 0.302 s |

Sources: `cover-profile-before.json`, `cover-profile-after.json` and matching logs. This is a diagnostic
reference translation through public APIs, not a special optimized workflow;
it omits Ext queries. It does not substitute a supplied algebra basis. The
after profile deliberately still requests a full projective family, whereas
ordinary native covers now request only selected members. Its total is **not**
the optimized end-to-end workflow's total.

Changes:

- `HomFromProjective` constructs only its selected `e_i A` instead of the entire
  family. A cover shares canonical-basis data across its ordered generators.
- Both selective and full-family requests recover each live weak-cache member
  before construction. An expired sibling no longer causes construction and
  subsequent disposal of a duplicate live member. Caches remain weak.
- Product coefficients go directly into the relevant vertex block through a
  word index, rather than repeatedly extracting a full-algebra coordinate
  vector. Indexes are lazy per projective: basis-only requests and live cache
  hits do not construct unused indexes. Quotient products must still carry
  normal form; every term is checked against the indexed canonical basis.
- Cover matrices are assembled by vertical concatenation in precisely the old
  summand/generator order. This equals the sum of projection composites without
  materializing those composites. The final map is still checked.
- Direct-sum component maps are lazy. Previously eager maps immediately expired
  because their only retained references were weak; asking for them rebuilt
  them anyway. Live requested maps retain their identity. Foreign-field NC
  cases retain their previous eager validation/error behavior.

Main implementation: `crates/qpa-gbnp/src/module/{algebra_invariants,projectives,
syzygies_and_tilting}.rs`, core projective-cache accessors, and
`crates/qpa-core/src/path_module/categorical.rs`.

## 3. Relation validation, not disabled checks

`crates/qpa-core/src/path_module/relation_checking.rs` replaces repeated complete
path products and scaled matrix copies with:

- borrowed input arrow matrices;
- reusable common prefixes, local to one validation call;
- exact zero-product descriptors that avoid allocating and multiplying zero
  matrices, while preserving field and shape checks;
- in-place scaled accumulation, avoiding a scaled temporary plus an added copy
  per term.

Prefix metadata is limited to 4,096 words of length at most 16. At most 1,024
values are retained, subject to an 8-MiB **inline coefficient-slot** budget.
Arbitrary-precision scalar payloads are additional; this is not a hard RSS cap.
There is no persistent cross-module validation cache.

A direct sum of already checked modules inherits uniform relations by block
composition. Only that private, structurally certified construction may reuse
validation; its field, label, dimension and matrix-shape checks still run.
Nonuniform relations and NC inputs with foreign runtime fields use the full
checker. In particular, QPA's physical zero-space stand-ins can make nonuniform
relation additions shape-dependent; the fallback retains those failures.
Indecomposable projectives and caller-supplied action matrices remain checked.
No caller-facing flag skips relation validation, and no supplied Gröbner basis
is treated as a new validation certificate.

## 4. Elimination after fill-in

Input-row order, increasing pivot-column sweep, pivot normalization and ordered
nullspace relations are unchanged.

- Individual elimination rows promote from ordered sparse maps to dense storage
  at 25% occupancy, subject to a shared 64-MiB inline coefficient-slot budget.
  A large sparse system never allocates its Cartesian product for this policy.
  This is one-way promotion; there is no claim of an optimal universal crossover.
- GF(2) systems can be packed directly from sparse nonzero positions, without
  first allocating a dense scalar matrix. The bound accounts for coefficient
  and relation bits when transformations are needed. Packing requires field
  **cardinality two**, not merely characteristic two.
- A measured sparse-GF(2) regression led to a density guard: inputs below 1/32
  occupancy stay sparse and may adapt later if rows fill. Packing also has a
  64-MiB bit-storage bound.
- Dense rank-only elimination no longer exports a second, unused basis matrix.
  Sparse rank-only elimination continues to omit transformation data.

The fill-in regression uses 96-by-768 exact matrices with late/early pivots and
dependent rows. Maximum input row occupancies were 78 (Q/GF(101)) and 76
(GF(2)/GF(4)); observed pre-promotion row occupancies reached 585, 583, 328 and
328, respectively. Ordered nullspaces match the original dense algorithm over
all four fields. See `elimination-fill-profile.log`.

Implementation: `crates/qpa-core/src/matrix/{sparse,adaptive_row,binary}.rs` and
`crates/qpa-core/src/matrix.rs`.

## Ordinary-API measurements

### Synthetic workloads

`benchmarks/native_performance.py` uses three samples of each case over Q,
GF(2), GF(101), and GF(4):

- Covers of repeated simples over a two-loop radical-four truncation and a
  12-vertex line. Each sample starts with a fresh algebra; top preparation is
  outside the cover timer.
- Checked reconstruction of the resulting projective action matrices, timed
  separately with action serialization outside the timer.
- Hom dimensions for sparse and dense presentations of a pair `(J,E_11)` on
  an 18-dimensional vector space over a free two-loop path algebra. Their
  common centralizer consists of scalars, so the expected Hom dimension is one
  over every field. The dense case conjugates by `I+u*v^t`, where `v^t*u=0`;
  this is invertible over every field, without random rational coefficient
  growth. GF(4) also uses a nonprime scalar.

Both runs use the identical runner. Raw before/after records, all samples,
binary hashes and cumulative process RSS are in `native-performance-before.*`
and `native-performance-after.*`. Per-case RSS is explicitly a process
high-water mark, not an independent case peak. These structured workloads do
not establish universal speedups.

Median milliseconds, three samples per case:

| Field | Operation | Before | After |
|---|---|---:|---:|
| Q | cover, loop | 13.096 | 1.147 |
| Q | checked module, loop | 11.366 | 6.227 |
| Q | cover, line | 0.589 | 0.130 |
| Q | checked module, line | 0.011 | 0.013 |
| Q | sparse Hom | 0.193 | 0.202 |
| Q | dense Hom | 1455.174 | 500.299 |
| GF(2) | cover, loop | 0.597 | 0.108 |
| GF(2) | checked module, loop | 0.557 | 0.442 |
| GF(2) | cover, line | 0.459 | 0.080 |
| GF(2) | checked module, line | 0.006 | 0.005 |
| GF(2) | sparse Hom | 0.068 | 0.102 |
| GF(2) | dense Hom | 312.477 | 0.787 |
| GF(101) | cover, loop | 0.801 | 0.113 |
| GF(101) | checked module, loop | 0.698 | 0.496 |
| GF(101) | cover, line | 0.483 | 0.083 |
| GF(101) | checked module, line | 0.007 | 0.006 |
| GF(101) | sparse Hom | 0.085 | 0.097 |
| GF(101) | dense Hom | 694.803 | 24.615 |
| GF(4) | cover, loop | 22.274 | 1.694 |
| GF(4) | checked module, loop | 17.760 | 8.203 |
| GF(4) | cover, line | 0.810 | 0.175 |
| GF(4) | checked module, line | 0.024 | 0.026 |
| GF(4) | sparse Hom | 0.463 | 0.456 |
| GF(4) | dense Hom | 555.270 | 283.002 |

All 24 outputs agree. The complete worker takes **9.37 -> 2.57 seconds**;
process peak RSS is **48,880 -> 50,680 KiB** (slightly higher, not lower).
The identical runner SHA-256 is
`64462611e852091f1389259ade7c7a66892c537541daefa1f7df333d0f771a3d`.

### Large ordinary-API workflow

The immediate baseline is `arc-ext-dimension-only.json`: 467.40 seconds,
2,860,504 KiB peak RSS, covers 132.131 / 185.108 seconds. It already includes
priority 1. The final after run is `arc-native-optimized.json`, starting from
original input documents through ordinary quotient/module APIs, with neither a
checkpoint nor an optimization switch. Selected source/target pairs and degrees
remain the existing `(Z,Z)` degrees one and two, not new research scope.

| Ordinary-API measurement | After priority 1 | After all four |
|---|---:|---:|
| Quotient construction | 30.249 s | 28.809 s |
| Checked witness module | 0.075 s | 0.017 s |
| First cover | 132.131 s | 4.631 s |
| First kernel | 6.463 s | 3.708 s |
| First Ext dimension | 4.906 s | 4.648 s |
| Second top projection | 36.350 s | 31.505 s |
| Second cover | 185.108 s | 41.288 s |
| Second kernel | 21.262 s | 13.331 s |
| Second Ext dimension | 50.118 s | 48.009 s |
| Whole worker | 467.40 s | 180.29 s |
| Whole-worker peak RSS | 2,860,504 KiB | 2,861,104 KiB |

The large run completes with both selected Ext dimensions zero and both short
exact sequences checked. Input hashes and all recorded module/projective/syzygy
dimension vectors agree. These are single large-process observations, not
medians; quotient construction was not changed in this pass.

Relative to priority 1, total time is about **2.59x faster**, but whole-process
peak RSS is essentially unchanged at **2.73 GiB**. First-cover cumulative RSS
falls from 1,689,192 to 361,968 KiB; later retained resolution data and Hom
elimination still determine the global peak. Relative to the original four-priority
baseline, total time is about **3x faster** and peak RSS falls from **4.36 to
2.73 GiB**. Do not attribute that whole peak-memory reduction to priorities 2–4.

## Regressions, rejected choices and limits

The first packed candidate sped up dense GF(2), but slowed the very sparse GF(2)
query from about 0.068 ms to 0.177 ms. It was not accepted as the final policy.
The input-density guard reduces that case to about 0.102 ms, but does **not**
eliminate all sparse-query overhead. The remaining roughly 34-microsecond
regression, and smaller Q/GF(101) sparse-query regressions, are visible in the
table. The policy trades small-call overhead for measured dense-system gains;
it is not a claim that every query became faster. Candidate artifacts are
retained as `native-performance-packed-candidate.*` and
`arc-native-optimized-candidate.*`; final claims use the final artifacts.

An exploratory rational dense fixture used 108 seeded random elementary basis
changes (`n=18`, seed 47019, row/column transvections of coefficient +/-1).
Its three-sample dense query batch had not completed when manually terminated
at 258.46 seconds. Partial JSON/stdout and signal-15 stderr are retained under
`native-performance-dense-exploratory.*`. That incomplete run is excluded from
completed comparisons. It exposes rational coefficient-growth sensitivity;
changing storage does not solve that problem, and no speedup is claimed for
this uncompleted fixture. The final dense fixture has a controlled exact basis
change rather than silently treating this aborted run as a completed baseline.

Tiny microsecond call differences are noisy. Report actual samples and absolute
times rather than extrapolating them to all operations. Memory budgets here
bound retained slots/bits, not arbitrary-precision payloads or output size.
Returning a full Hom basis can still require large dense output. These changes
do not promise that all computations over every finite algebra will be cheap.

## Validation and reproduction

Final validation passes:

- **230 core + 216 GBNP unit tests**, **4 integration tests**, **2 doctests**;
- **611 Python tests**, including the frozen compatibility coverage;
- strict Clippy, formatting, Python syntax, and artifact/input consistency checks;
- all 24 general benchmark outputs and the large ordinary-API results agree;
- all eight dimension-only Ext benchmark outputs agree with the priority-1 run.

The retained dimension-only benchmark's latest records are
`ext-dimension-only-final.*`. Relative to priority 1, its median batch times
(ms) are:

| Field | Mixed loop, before -> after | Kronecker, before -> after |
|---|---:|---:|
| Q | 4.026 -> 3.320 | 1.601 -> 1.881 |
| GF(2) | 0.308 -> 0.400 | 0.478 -> 0.721 |
| GF(101) | 0.386 -> 0.482 | 0.555 -> 0.770 |
| GF(4) | 5.770 -> 5.679 | 3.460 -> 3.820 |

These small, already-optimized queries show residual overhead too; they are not
hidden behind the large-workflow result. They remain substantially faster than
the pre-priority-1 basis/restriction implementation. No ordered basis API was
replaced by a dimension-only answer.

New regressions cover
selective cache lifetimes/identity, local projective coordinates against the old
full-coordinate construction, direct block maps against projection sums,
nonprime coefficients, invalid input relations, NC field errors, nonuniform
zero-space shapes, fill-in, density/memory guards, and ordered nullspaces.
Frozen QPA fixtures and the completed dimension-only Ext tests remain in the
full suites.

Commands from the repository root:

```bash
cargo test -p qpa-core -p qpa-gbnp
cargo clippy -p qpa-core -p qpa-gbnp --all-targets -- -D warnings
cargo fmt --all -- --check
VIRTUAL_ENV="$PWD/.venv" .venv/bin/maturin develop --release
.venv/bin/python -m pytest python/tests -q

/usr/bin/time -v timeout -k 2 300 .venv/bin/python benchmarks/native_performance.py \
  --output benchmarks/results/native-performance-reproduced.json
/usr/bin/time -v timeout -k 2 900 .venv/bin/python benchmarks/arc_ext.py --degree 2 \
  --output benchmarks/results/arc-native-optimized-reproduced.json
/usr/bin/time -v timeout -k 2 900 .venv/bin/python benchmarks/arc_ext.py --degree 2 \
  --profile-cover --output benchmarks/results/cover-profile-reproduced.json
/usr/bin/time -v timeout -k 2 300 .venv/bin/python benchmarks/ext_dimension_only.py \
  --output benchmarks/results/ext-dimension-only-reproduced.json
cargo test -p qpa-core fill_in_adaptation -- --nocapture
```

Synthetic workers have a 4-GiB address-space guard; the large worker has an
8-GiB guard. Large-case inputs, conventions and hashes are unchanged from the
previous ordinary-API reports and remain read-only. Builds are modified
qpa-rust 0.2.0 working trees at base revision
`6207cd941fefcadab5a695001be3aac573b1dacf-dirty`, on the same Linux/WSL2 machine.
The base commit alone does not reproduce these changes. Before measurements
require the pre-priorities-2–4 implementation, identified by extension SHA-256
`a84d9a1be30948cfd4b1fe7ec48141213ab38756aa7c1a44794cba98ea2430b9`.
The final extension SHA-256, matching all final after JSON files and the
installed extension, is
`9632c6ac045e363142645a2f3f5d301fe1715f1b355b1513f9490a307d1ef131`.
