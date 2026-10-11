# ARC presentation: guarded construction diagnosis

This is the **pre-optimization baseline**. See the
[optimization follow-up](arc-optimization.md) for subsequent implementation,
completed rational materialization, and selected low-degree Ext calculations.

Date: 2026-10-10. Scope: performance diagnosis authorized after the
[small synthetic pilot](algebra-construction.md). **No Ext computations,
module materialization, or all-degree ARC verification were run.** No algebra
backend code or neighboring research-repository files were changed.

## Main result

There are **two distinct bottlenecks**, and presentation orientation matters:

1. **The published rational presentation stalls in GBNP completion.** It did
   not finish within 120 seconds. A native profile locates the sampled work in
   strong normal-form reductions and multiprecision rational arithmetic,
   especially integer GCDs.
2. **The opposite rational presentation completes its Gröbner basis in
   24.18 seconds**, and then its quotient, dimension and canonical basis are
   cheap. The computed dimension is **7708**. But its subsequent admissibility
   check still did not finish within a further 120 seconds. A second native
   profile locates that work in construction of the initial arrow-ideal span,
   multiplying and reducing very many quotient elements.

Thus the generic constructor can be blocked *after* successful rational
Gröbner completion. Conversely, the first bottleneck in the published
orientation really is completion, not admissibility or Ext.

The opposite orientation is relevant because the recorded witness is a left
module, represented by the Python QPA right-module API over the opposite
algebra. We reverse both arrow directions and words, keep coefficients and
arrow insertion order unchanged, and do not otherwise simplify the input.
This is a different presentation/order for Gröbner purposes, not a change of
coefficient field or an approximation.

## Input and environment

Read-only input:

```text
../homological-conjectures/results/auslander-reiten-counterexample/arc-expanded-bound-quiver.md
SHA-256: 78a50028b95bd11ec02f226fef7560b7b03c0ac3e10b8d03e4a73c1b6776a69d
```

The parsed vertices, arrows and coefficient/path lists were also checked equal
to the local `data/auslander-reiten-counterexample/arc-explicit-presentation.json`.
That JSON is not required for reproduction from the tracked Markdown. No
consumer experiment was imported or executed to load the presentation.

Input size:

- 8 vertices, 43 arrows, 106 defining relations;
- 957 polynomial terms, maximum 94 terms in one relation;
- path lengths from 2 to 10;
- endpoint composability and uniform source/target of every relation checked.

Engine: local release build of qpa-rust 0.2.0, revision
`6207cd941fefcadab5a695001be3aac573b1dacf`, Python 3.12.11. Extension SHA-256:
`cbee0a2c587509df4185361112c8a47d807635bf80990741e5cc4a38007f83c8`.
Hardware: i5-13500, WSL2 Linux x86-64, affinity 0–19, unpinned. These large
probes use **Rust only**, not GAP. The neighboring repository's QPA preflight
also passed; its installed release is not the binary used for these timings.

The runner is [`../arc_construction.py`](../arc_construction.py). Its source
SHA-256 for these artifacts is
`617ffc9b8434eebb5525d7dfb221705e7624102345fc76d46f14a019b7a06e47`.

## Guarded stage measurements

Each unprofiled configuration was run once in a fresh process. Numbers below
are **seconds**, not medians. Parsing, presentation construction, ideal
construction, completion, basis attachment, quotient construction, dimension,
canonical basis and admissibility are separate calls.

Importantly, dimension and canonical basis were requested **before**
admissibility in this diagnostic, so an expensive admissibility check could
not hide their costs. The direct relation-list constructor normally performs
its eager admissibility check before returning the quotient. This staged
sequence is not a measurement of an unchanged direct constructor call.

| Stage | Published orientation / Q | Opposite / Q | Published orientation / GF(101) |
|---|---:|---:|---:|
| Parse list | 0.00183 | 0.00187 | 0.00189 |
| Construct native presentation | 0.00474 | 0.00549 | 0.00465 |
| Input ideal | 0.000057 | 0.000098 | 0.000068 |
| GBNP completion | **>120, stopped** | **24.177** | **6.787** |
| Completed GB relation count | not obtained | 1406 | 1376 |
| Completed ideal + GB attachment | not reached | 0.00467 | 0.00251 |
| Quotient with attached GB | not reached | 0.05745 | 0.04894 |
| Dimension | not reached | 0.00265 | 0.00204 |
| Canonical basis | not reached | 0.01236 | 0.01165 |
| Dimension / basis size | not obtained in this run | **7708 / 7708** | **7708 / 7708** |
| Admissibility | not reached | **>120, stopped** | **>60, stopped** |

All algebra workers had a **4-GiB address-space limit**. They stopped on wall
time, not on memory limits. One-second RSS samples at cutoff were approximately
42.6 MiB, 414.5 MiB and 183.4 MiB respectively, in table-column order. These are
sampled resident sizes, not guaranteed peak memory; startup high-water marks
are separately recorded where available. No workers were left running.

The finite-field probe isolates a cheaper arithmetic regime but changes the
mathematical problem. Its 7708-dimensional result is **not a rational
reconstruction or a finite-field ARC claim**. The rational opposite computation
is a direct exact calculation on the opposite of the supplied presentation;
it supplies the dimension of that presentation over Q, without independently
auditing its identification with the recorded research construction.

## Native profile 1: rational GBNP completion

Separate run using `py-spy --native`, 49 Hz, stopped after approximately 30
seconds in completion. There were **1478 sampled stacks** (30.16 seconds of
sample weight). Percentages below are **inclusive**, so they overlap and must
not be added:

| Frame or category present in sampled stack | Share |
|---|---:|
| `strong_normal_form_trees_skipping_secondary` | 94.52% |
| `reduce_indexed_term` | 91.27% |
| `np::add_np` | 90.46% |
| Numerica rational-arithmetic frames | 89.17% |
| Integer GCD frames | 56.83% |

The heaviest leaf frames include Malachite multiprecision half-GCD, division
and multiplication routines. This is evidence for a **rational-coefficient
arithmetic bottleneck during reductions**, consistent with coefficient swell;
we did not instrument intermediate numerator/denominator bit lengths.
The timed-out 120-second run does not supply a total completion time.

Relevant source chain:

```text
qpa_gbnp/src/qpa_conversion.rs: gbnp_groebner_basis
qpa_gbnp/src/sgrobner.rs: sgrobner_with_options
                        add_obstructions_to_heap
                        strong_normal_form_trees_skipping_secondary
qpa_gbnp/src/reduction.rs: reduce_indexed_term
qpa_gbnp/src/np.rs: add_np
Numerica rational arithmetic -> multiprecision integer arithmetic / GCD
```

These paths are under `crates/`. The profile does not identify Python parsing,
input conversion, or memory exhaustion as the cause of the stall.

## Native profile 2: rational opposite admissibility

A separate opposite-Q profiling run completed GBNP in 24.58 seconds, obtained
1406 relations and dimension/basis size 7708, then was stopped after roughly
45 seconds in admissibility. Of 3445 total sampled stacks, **2235** contain an
`is_admissible` frame (45.61 seconds of sample weight). Among those stacks:

| Frame or category | Inclusive share |
|---|---:|
| `quotient_ideal_span` | **99.96%** |
| `QuotientPathAlgebraElement::multiply` | 98.39% |
| `GroebnerBasis::completely_reduce` | 94.90% |
| Numerica rational-arithmetic frames | 13.96% |
| Integer GCD frames | 1.66% |

Here the Gröbner basis is **already computed**. `completely_reduce` is repeated
normal-form reduction of products, not another GBNP completion. Path equality,
path-algebra addition, ordering and allocation are prominent leaf costs. No
product-space-squaring or echelon-reduction frames were observed in these
admissibility samples: the sampled work is overwhelmingly the initial ideal
span construction, before those later steps.

Source inspection gives a concrete scale for that work. In
`crates/qpa-core/src/quotient.rs`, `PathAlgebraIdeal::is_admissible` calls
`quotient_ideal_span` on the quotient arrows. That helper executes:

```text
for left in algebra_basis:
    for arrow in arrow_generators:
        left_product = left * arrow
        for right in algebra_basis:
            product = left_product * right
            retain product if nonzero
```

For this input, completion of that loop entails

```text
7708 * 43 * 7708 = 2,554,770,352
```

inner multiplication calls, in addition to the left products. In particular,
the code does not skip the inner loop when `left_product` is zero. It retains
nonzero products and then passes them to a dense span calculation; the later
radical-squaring loop follows. The 2.55-billion figure is a **source-derived
loop count**, not a measured count of products actually executed before cutoff.

## Interpretation

- The original diagnosis is correct for the **published orientation over Q**:
  it stalls in completion before reaching Ext.
- But “Gröbner completion makes every rational materialization impossible” is
  too broad. The **opposite presentation already completes** and yields the
  full 7708-element rational canonical basis in about 24.3 seconds in this
  staged route.
- For the left-module/opposite-algebra workflow, the next serious obstruction
  is **eager admissibility**, specifically massive ideal-span product
  enumeration and normal-form reduction. Optimizing GBNP alone will not fix it.
- The two useful optimization targets are therefore separate: rational
  polynomial-reduction arithmetic, and the admissibility/ideal-span path.
  Any change must preserve the backend's exact semantics. Simply disabling
  admissibility everywhere is not a justified generic fix.
- No conclusions about actual Ext runtime are possible: no Ext was attempted.
  These are bounded diagnostics, not completed high-level-constructor runs or
  an all-degree mathematical verification.

## Reproduction and artifacts

From the qpa-rust root:

```sh
# Published orientation over Q: 120 seconds per stage, 600 seconds total.
.venv/bin/python benchmarks/arc_construction.py \
  --compare-json ../homological-conjectures/data/auslander-reiten-counterexample/arc-explicit-presentation.json \
  --output benchmarks/results/arc-construction-rational.json

# Arithmetic diagnostic only: not an ARC claim in finite characteristic.
.venv/bin/python benchmarks/arc_construction.py --field GF101 \
  --stage-timeout 60 --total-timeout 240 \
  --output benchmarks/results/arc-construction-mod101.json

# Opposite presentation over Q.
.venv/bin/python benchmarks/arc_construction.py --opposite \
  --stage-timeout 120 --total-timeout 600 \
  --output benchmarks/results/arc-construction-opposite-rational.json

# Native profiles; requires py-spy, no elevated privileges in these runs.
.venv/bin/python benchmarks/arc_construction.py --profile --stop-after groebner \
  --stage-timeout 30 --total-timeout 45 \
  --output benchmarks/results/arc-construction-rational-profile.json

.venv/bin/python benchmarks/arc_construction.py --opposite --profile \
  --stage-timeout 45 --total-timeout 90 \
  --output benchmarks/results/arc-construction-opposite-rational-profile.json
```

Omit `--compare-json` if only the tracked Markdown input is available. Each
output stem also has `.stdout.log` and `.stderr.log` files. Profile stems have
`.speedscope.json` files, viewable with a speedscope-compatible viewer. The
runner saves events and `/proc` resource samples incrementally, enforces
per-stage and total wall limits, and kills remaining descendants on exit.
Profiling runs are kept separate from unprofiled timing observations.

The reported inclusive percentages are computed by summing sample weights
whose stack contains the named frame/category, divided by total selected
sample weight. For the second profile only stacks containing `is_admissible`
are selected. Each category is counted at most once per stack.

Checks performed: Markdown/JSON equality; path/endpoint and relation schema
validation; a negative unknown-arrow control; runner syntax compilation;
source-hash consistency across all five reports; completed-phase dimension and
basis-size checks; process cleanup; and `git diff --check`. No claim is made
that a timeout is a mathematical failure or evidence of nontermination.
