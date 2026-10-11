# ARC: optimized rational construction and low-degree self-Ext

**Historical staged-workflow report.** The manual split described below has
since been replaced by automatic, general component preconditioning in the
native QPA constructor. See [ordinary API results](arc-native-api.md) for the
current workflow and its separately measured end-to-end Ext replay. The old
`--precondition` switch has been removed; command lines below record the
historical runs, not the current constructor CLI.

## Result and scope

Both the published ARC presentation and its opposite now materialize over
**Q**, with **dimension 7708**, canonical basis size 7708, and successful
admissibility checks. Fresh-process construction, including completion and
admissibility, takes about **30 seconds** and **102–104 MiB peak RSS** here.

The explicitly selected module pair and degrees are **(Z, Z), n = 1, 2**, where
Z is the recorded 110-dimensional witness. The exact results are

```text
Ext^1_A(Z,Z) = 0
Ext^2_A(Z,Z) = 0
```

Here these are left A-modules, implemented as right modules over A^op by
transposing the published column-action matrices. No claim is made about
higher degrees, Ext^n(Z,A), or all-degree ARC verification. This is a
computational check of the supplied presentation and witness, not an
independent audit of their identification with the research construction.
The neighboring research repository was read only.

## Changes and preservation of exact results

1. **Native block preconditioning for completion.**
   `benchmarks/arc_materialize.py` first completes the 62 diagonal relations,
   then reduces the 44 connecting relations by that native basis before
   completing their union. Replacing a connecting relation r by NF_D(r)
   preserves the combined ideal, since their difference is in the diagonal
   ideal. Completion and reduction still use native QPA/GBNP operations;
   there is no modular reconstruction, coefficient approximation, or new
   Python algebra engine. This is an explicit ARC workflow, not an automatic
   rewrite of every generic Gröbner-basis request.
2. **Sparse, cached admissibility.**
   `crates/qpa-core/src/quotient/admissibility.rs` streams exact echelon spans,
   caches multiplication by arrows, and reuses nontip prefixes. It retains
   QPA's descending arrow-ideal/product-space test until dimension stabilizes.
   Finite-dimensionality alone is **not** treated as admissibility. The result
   is cached on the ideal.
3. **Sparse Hom equations.**
   `crates/qpa-core/src/matrix/sparse.rs` and
   `crates/qpa-core/src/path_module/homomorphisms.rs` avoid allocating the full
   unknown-by-constraint matrix. Elimination preserves input-row order,
   increasing pivot-column sweeps, normalization, and ordered relation vectors.
   Small systems retain the existing dense/packed kernels; unconstrained
   dimensions do not allocate one empty container per unknown.
4. **Bounded subrepresentation closure.**
   `crates/qpa-core/src/path_module/submodules.rs` retains independent frontier
   vectors instead of every nonzero path image. By linearity, discarded
   dependent vectors cannot add new arrow-image spans. The final fully reduced
   canonical basis, and hence its module matrices, remain the same. A single
   coordinate solver is reused for all images in that basis.
5. **Vector-first path actions.**
   `crates/qpa-core/src/path_module/objects.rs` evaluates
   v(A1...Ak) as ((vA1)A2)...Ak, avoiding full path-action matrices for a single
   vector while retaining exact arithmetic and arrow shape/field checks.

## Before/after construction measurements

The [baseline diagnosis](arc-construction.md) records the unoptimized build,
commands, guards and native profiles. These are single fresh-process probes,
not medians or a general GAP/Rust performance comparison. The final probes
ran sequentially, without a profiler, on the same i5-13500 / WSL2 Linux
x86-64 environment, with Python 3.12.11 and a release build of qpa-rust 0.2.0.
Construction uses a 4-GiB address-space guard.

| Measurement | Published / Q, before | Published / Q, after | Opposite / Q, before | Opposite / Q, after |
|---|---:|---:|---:|---:|
| GBNP completion, seconds | >120, stopped | 14.141 | 24.177 | 11.428 |
| Diagonal completion + connecting reduction, seconds | not applied | 0.0254 | not applied | 0.0263 |
| Completed GB relations | not obtained | 1376 | 1406 | 1406 |
| Admissibility, seconds | not reached | 14.002 | >120, stopped | 15.646 |
| Dimension / canonical basis size | not obtained | 7708 / 7708 | 7708 / 7708 | 7708 / 7708 |
| Whole-process wall time, seconds | incomplete | 29.80 | incomplete | 30.01 |
| Memory | 42.6 MiB sampled at cutoff | 101.98 MiB peak RSS | 414.5 MiB sampled at cutoff | 103.60 MiB peak RSS |

The baseline memory figures are one-second samples at cutoff, **not guaranteed
peaks**, and the baseline never finished the full workflow. They must not be
read as an exact completed-run memory ratio. Stage timings use
`time.perf_counter`; whole-process wall time and peak RSS use GNU `time -v`.
Stage sums are 28.468 and 27.307 seconds, respectively; startup, bookkeeping,
and teardown are outside those timed stage calls.

Final raw records:

- `arc-final-rational.stdout.log` and `.stderr.log`;
- `arc-final-opposite.stdout.log` and `.stderr.log`.

Every original relation reduces to zero under the attached completed basis.
The quotient dimension and the enumerated canonical basis size are checked
separately, followed by the native admissibility test. The staged construction
is not a timing of an unchanged eager relation-list constructor.

## Exact Ext computation

Witness vertex order and dimension vector:

```text
L0 L1 L2 L3 U0 U1 U2 U3
26 30 40 13  0  0  0  1
```

The loader checks 43 arrow records, 439 nonzero entries, declared matrix
shapes, vertex/arrow order and exact rational coefficients. The native checked
module constructor verifies the defining relations. Opposite construction
reverses arrow directions **and words**, retaining coefficients and arrow
insertion order; transposing the witness matrices then gives QPA row actions.

The runner computes Ext^1 natively. For degree two it uses the projective
cover sequence and dimension shifting:

```text
Ext^2(Z,Z) = Ext^1(Omega Z,Z).
```

| Resolution step | Projective dimension vector | Projective dimension | Kernel dimension vector | Kernel dimension |
|---|---|---:|---|---:|
| P0 -> Z | [282,286,440,166,6,2,8,6] | 1196 | [256,256,400,153,6,2,8,5] | 1086 |
| P1 -> Omega Z | [678,610,1032,438,16,16,24,8] | 2822 | [422,354,632,285,10,14,16,3] | 1736 |

At both steps, native checks confirm surjectivity of the cover, injectivity of
the kernel inclusion, zero composition, and the dimension identity for the
short exact sequence. The source tops are [0,0,0,1,0,0,0,1] and
[0,0,0,1,0,0,1,0], respectively.

| Timed operation | Seconds |
|---|---:|
| First top projection | 0.069 |
| First projective cover | 22.331 |
| First kernel | 0.928 |
| Ext^1(Z,Z), exact dimension 0 | 10.049 |
| Second top projection | 26.515 |
| Second projective cover | 31.507 |
| Second kernel | 3.131 |
| Ext^1(Omega Z,Z), exact dimension 0 | 134.673 |

The whole Ext worker, **starting from a saved completed basis**, took
**258.51 seconds** and peaked at **4,576,808 KiB = 4.36 GiB RSS**. It includes
basis loading, relation checks, quotient construction and a fresh 15.619-second
admissibility check, but not fresh Gröbner completion. Its guards were
900 seconds and 8 GiB address space. It exited successfully.

Authoritative result: [`arc-ext-final.json`](arc-ext-final.json), with its
matching `.stdout.log` and `.stderr.log`.

### Failed intermediate runs, retained rather than reported as results

- `arc-ext1.*`: dense Hom equations requested **239,955,846,400 bytes** and
  aborted under the address-space guard. No Ext dimension was obtained.
- `arc-ext-sparse.*`: sparse Hom completed Ext^1 = 0, but the 300-second guard
  stopped the second cover. No Ext^2 result was obtained.
- `arc-ext-vector-action.*`: after vector-first actions, the second top
  construction reached the 8-GiB guard while retaining redundant path images.
  The final independent-frontier closure resolves that failure. These files
  were renamed from an unsuccessful `arc-ext-final` attempt.
- `arc-optimized-*`, `arc-precondition-baseline.*`, and
  `arc-sparse-admissibility.*` retain earlier construction experiments. The
  final tables above use the `arc-final-*` files instead.

## Inputs and build provenance

Read-only source documents, relative to the qpa-rust root:

```text
../homological-conjectures/results/auslander-reiten-counterexample/arc-expanded-bound-quiver.md
SHA-256 78a50028b95bd11ec02f226fef7560b7b03c0ac3e10b8d03e4a73c1b6776a69d

../homological-conjectures/results/auslander-reiten-counterexample/arc-expanded-witness.md
SHA-256 3631e1f3ff3188f965612f7c82d6e8f796d4c3730707a044d599374110b35e61
```

Final construction and Ext records all identify revision
`6207cd941fefcadab5a695001be3aac573b1dacf-dirty` and the same extension SHA-256:

```text
457d89fa5a056fa0cc5a1e27378253b5d5b7470088e8ee161a0c0f6d50d7585a
```

The dirty suffix matters: the modifications described here are not contained
in that base commit alone. Commands require the corresponding working-tree
changes. The final opposite checkpoint is byte-identical to the earlier
`target/arc-optimized-opposite-basis.json` used by the Ext run:

```text
target/arc-final-opposite-basis.json
SHA-256 336916a5f7ff33aef35e69f54139b560a7ae2d01d2e2ba5808f5c12e38726801

target/arc-final-rational-basis.json
SHA-256 28db4cc547eb76a69af6df3e5af7ca146e9d16a9f752c4386517a32fdea19429
```

These completed-basis checkpoints are generated under ignored `target/`; they
can be regenerated from the source documents with the commands below.

## Reproduction and validation

For the **before** measurements, use the commands in the baseline report with
an unmodified release build of commit `6207cd941fefcadab5a695001be3aac573b1dacf`
in a separate environment. Running that runner against the optimized extension
would measure the raw presentation on the new engine, not the historical
baseline.

For the **after** measurements, from the qpa-rust root, with the read-only
inputs above available (or supply `--input` and `--witness` explicitly):

```bash
VIRTUAL_ENV="$PWD/.venv" .venv/bin/maturin develop --release

/usr/bin/time -v timeout -k 2 180 .venv/bin/python benchmarks/arc_materialize.py \
  --precondition --admissibility --save-basis target/arc-final-rational-basis.json

/usr/bin/time -v timeout -k 2 180 .venv/bin/python benchmarks/arc_materialize.py \
  --opposite --precondition --admissibility \
  --save-basis target/arc-final-opposite-basis.json

/usr/bin/time -v timeout -k 2 900 .venv/bin/python benchmarks/arc_ext.py \
  --basis target/arc-final-opposite-basis.json --degree 2 \
  --output benchmarks/results/arc-ext-reproduced.json
```

Validation on the implementation used for the final runs:

```bash
cargo test -p qpa-core -p qpa-gbnp --lib
# 224 core + 210 GBNP tests passed
cargo clippy -p qpa-core -p qpa-gbnp --all-targets -- -D warnings
# passed
.venv/bin/python -m pytest python/tests -q
# 597 passed
```

Differential/regression coverage includes:

- sparse admissibility versus the previous dense test over Q, GF(2), and
  GF(101), including nonhomogeneous nilpotent examples, finite-dimensional
  **nonadmissible** examples, multiple vertices and the empty arrow ideal;
  cached products are compared with every native basis product in these cases;
- 2,430 matrices over Q, GF(2), GF(101): sparse rank and **ordered** nullspaces
  match dense results, and incremental accepted rows have the same canonical
  row space as all original rows;
- vector-first actions versus the previous full matrix-chain evaluation;
- Python sparse-Hom/Ext checks over Q, GF(2), GF(5), GF(4), including nonprime
  extension coefficients; an unconstrained Hom dimension of 10^10 without
  allocating that many rows; and two identical nilpotent arrow actions whose
  literal word-image closure would retain exponentially many duplicate vectors;
- the existing frozen QPA compatibility fixtures and the complete Python suite.
