# Small presentation-to-algebra benchmark

## Scope and conclusion

This is a bounded synthetic pilot motivated by the presentation-preparation
bottleneck recorded in `../homological-conjectures`. **No ARC presentation,
witness, materialization script, or Ext calculation was executed.** Only the
neighboring repository's explanatory documents and source were read.

All 288 samples completed: 12 presentations × 2 fields × 2 backends × 2
construction routes × 3 fresh-process repetitions. Algebra dimensions range
from **4 to 40**. Every computed dimension and basis size equals its independent
family formula. Gröbner relation counts also agree across both engines and
routes in these samples. These checks are sanity checks, not a general proof
of backend equivalence.

The main finding is that **the ordinary relation-list constructor does more
than Gröbner completion**. In this pilot, the subsequent admissibility check
often costs substantially more than completion. Dense changes of generators
also make the *same small algebra* significantly more expensive. Neither
finding establishes which stage dominates the unrun ARC presentation.

## Reproduce

From the repository root, with the compiled Python extension and GAP/QPA installed:

```sh
.venv/bin/python benchmarks/algebra_construction.py \
  --repeat 3 --route staged --route direct \
  --output benchmarks/results/algebra-construction.json
```

Smaller smoke test:

```sh
.venv/bin/python benchmarks/algebra_construction.py \
  --case quantum-2 --case dense-3 --case grid-2 --field Q \
  --repeat 1 --route staged --route direct --output /tmp/construction-smoke.json
```

The coordinator is Linux/POSIX-specific (`resource`, process groups,
`/usr/bin/time`). Each sample has a **10-second wall timeout** and **1-GiB
address-space limit**. GAP additionally receives `-m 64m -o 512m -K 768m`.
Timeouts kill the sample's entire process group; failed/timed-out configurations
are not repeated. Partial stage timings and errors are retained, and the report
is saved after every configuration. A separate 0.01-second timeout smoke test
confirmed that the timeout path is reported correctly.

## Environment and methodology

- Intel Core i5-13500, WSL2 Linux x86-64; CPU affinity 0–19, not pinned.
- Python 3.12.11.
- qpa-rust revision `6207cd941fefcadab5a695001be3aac573b1dacf`, release extension.
  The installed extension hash matches `target/release/lib_qpa.so`:
  `cbee0a2c587509df4185361112c8a47d807635bf80990741e5cc4a38007f83c8`.
- GAP 4.11.1, installed **QPA 1.36**, GBNP 1.1.0. The Rust port targets QPA
  1.37, so this is explicitly **not a pinned-version differential benchmark**.
- Exact arithmetic over Q and GF(101); no modular reconstruction.
- Serial execution, three fresh processes per configuration. These are cold
  API/cache timings, not steady-state warmed microbenchmarks. No core library
  implementation was changed for the experiment.
- Stage timers exclude interpreter/package startup and presentation-spec
  generation, but include native API/wrapper work and native relation creation.
  Process wall time, including startup, is recorded separately. Rust uses
  `perf_counter`; GAP uses `NanosecondsSinceEpoch` differences.
- Peak RSS covers the whole process, including startup and provenance hashing;
  maxima were **48.3 MiB Rust** and **95.1 MiB GAP**. This is not incremental
  Gröbner memory usage.
- The JSON retains exact input coefficient/path lists, individual samples,
  backend metadata, source/extension hashes, and medians. `median_pipeline_seconds`
  is the median of each sample's sum of stage times, not a sum of stage medians.

### Stages and the hidden admissibility cost

The **direct** route times `PathAlgebra.quotient(relations)` / `KQ / rels`,
then dimension and canonical-basis materialization.

The **staged** route follows its substantive work explicitly:

1. Construct field, quiver, path algebra, and relation elements.
2. Construct an input ideal and call `GBNPGroebnerBasis(I)`.
3. Construct the defining ideal on the **completed** relation list, matching
   the direct constructor, and attach that already-computed basis.
4. Form the quotient without another completion.
5. Check `IsAdmissibleIdeal`.
6. Request dimension and materialize the canonical basis.

The stages are API-level measurements, not instrumentation inside `SGrobner`.
“Gröbner” includes the checked conversion to/from GBNP. Admissibility can itself
construct a quotient, enumerate a basis and form product spaces; its timer is
inclusive of that work. The final basis timer therefore measures the remaining
work on the benchmark quotient, not every basis operation performed earlier.
The staged route adds some wrapper overhead and does not reproduce every
attribute/filter cache side effect of the direct constructor.

Source inspection explains why step 5 must not be silently omitted:

- `crates/qpa-gbnp/src/pathalg.rs::quotient_by_relations` computes the basis
  and calls `by_complete_groebner_relations`.
- `crates/qpa-core/src/quotient.rs::by_relators_and_complete_groebner_basis`
  additionally calls `ideal.is_admissible()`.
- `PathAlgebraIdeal::is_admissible` constructs the quotient arrow ideal and
  repeatedly squares its product space until the dimension stabilizes.
- Upstream QPA's relation-list `/` method in `lib/pathalg.gi` and
  `IsAdmissibleIdeal` in `lib/pathalgideal.gi` follow the same pattern.

An initial smoke split that omitted this check understated direct-constructor
cost. The saved sweep includes it explicitly; those preliminary timings are
not part of the result artifact.

## Artificial families

- **`monomial-m-d`**: one vertex, m loops, every length-d word killed.
  Dimension is `1 + m + ... + m^(d-1)`. Cases `(m,d)=(2,3),(2,4),(3,3),(3,4)`.
- **`quantum-n`**: one vertex, n loops y_i, relations `y_i^2=0` and
  `y_j*y_i - 2*y_i*y_j=0` for i<j. Ordered square-free words give dimension
  `2^n`; n=2,3,4.
- **`dense-n`**: exactly the previous algebra after substituting
  `y_i = x_i + sum_j x_j`. The change matrix is `I + all-ones`, with determinant
  n+1, nonzero over both tested fields. Dimensions remain `2^n`, while the
  relation polynomials and completed bases become denser.
- **`grid-n`**: an oriented n×n grid, with its elementary squares commuting.
  This is the product-of-chains incidence algebra, of dimension
  `[n(n+1)/2]^2`; n=2,3. It exercises multiple vertices without cycles.

## Rational results

All times below are **milliseconds**, medians of the staged route. Total
includes presentation, both ideal constructions, attachment, quotient,
admissibility, dimension, and canonical basis.

| Case | dim | completed GB size | Rust GB | Rust admissibility | Rust total | GAP GB | GAP admissibility | GAP total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| monomial-2-3 | 7 | 8 | 0.128 | 0.173 | 0.515 | 0.670 | 6.329 | 11.918 |
| monomial-2-4 | 15 | 16 | 0.250 | 0.769 | 1.248 | 1.610 | 24.332 | 64.897 |
| monomial-3-3 | 13 | 27 | 0.288 | 0.575 | 1.162 | 2.380 | 19.538 | 31.659 |
| monomial-3-4 | 40 | 81 | 1.412 | 5.899 | 7.978 | 14.776 | 198.048 | 251.481 |
| quantum-2 | 4 | 3 | 0.096 | 0.110 | 0.390 | 0.404 | 3.868 | 7.672 |
| dense-2 | 4 | 4 | 0.165 | 0.128 | 0.524 | 0.616 | 3.996 | 8.009 |
| quantum-3 | 8 | 6 | 0.122 | 0.592 | 0.959 | 0.569 | 16.147 | 22.946 |
| dense-3 | 8 | 9 | 0.631 | 1.599 | 2.592 | 3.288 | 28.917 | 43.674 |
| quantum-4 | 16 | 10 | 0.148 | 3.862 | 4.215 | 0.828 | 102.574 | 114.273 |
| dense-4 | 16 | 20 | 7.384 | 34.960 | 42.797 | 31.125 | 330.540 | 411.629 |
| grid-2 | 9 | 1 | 0.077 | 0.123 | 0.375 | 0.606 | 3.315 | 12.944 |
| grid-3 | 36 | 4 | 0.102 | 2.130 | 2.482 | 0.577 | 34.982 | 160.844 |

The grid-3 GAP total is not dominated by Gröbner work: presentation costs about
76 ms and final basis materialization about 48 ms. This is another reason not
to label all construction time “Gröbner time.”

### Same algebra, different presentation

At dimension 16, the dense substitution increases Rust Gröbner time from
**0.148 to 7.384 ms (~50×)** and total preparation from **4.215 to 42.797 ms
(~10×)**. Its completed basis grows from 10 to 20 relations. GAP Gröbner time
increases from **0.828 to 31.125 ms (~38×)**.

For dense-4 over Q, admissibility is about **82%** of Rust's staged total and
**80%** of GAP's. These are small-case observations, not asymptotic estimates.

### Finite-field comparison

For dense-4 over GF(101):

| Backend | GB (ms) | admissibility (ms) | total (ms) |
|---|---:|---:|---:|
| Rust | 2.846 | 21.033 | 24.562 |
| GAP | 24.905 | 313.447 | 385.686 |

The cheaper finite-field arithmetic helps the Rust dense example, but does not
remove the admissibility cost. This comparison says nothing about transfer of
ARC properties between characteristics.

## Implications and next steps

1. **Separate completion from admissibility/product-space work** in future
   profiling. Optimizing GBNP alone would leave most of dense-4's preparation
   time untouched in this pilot.
2. **Keep presentation density as a benchmark axis**, not just dimension or
   relation count. Dense coordinate changes are a controlled way to do this
   without invoking any research counterexample.
3. Rust is faster on these cases, but the dataset is tiny, cold-started, and
   the QPA versions differ. It does not establish that large presentations
   are now feasible, nor reproduce the reported large-case stall.
4. If a second pilot is wanted, extend only these synthetic families by one
   size step under the same guards, and profile admissibility's product-space
   and normal-form operations. Do not jump to the ARC algebra.
5. The witness manifest's quotient-without-Gröbner route is appropriate for
   checking supplied action matrices. It is **not** evidence that generic
   dimension, projective, or Ext computations can bypass completion. Likewise,
   avoiding the eager admissibility check via a lower-level ideal quotient
   changes what has been checked/cached; it is not a blanket drop-in fix.

No optimization or semantic change to either algebra backend is claimed here.

A subsequently authorized [large ARC presentation diagnosis](arc-construction.md)
separates rational Gröbner-completion cost from the much larger admissibility
product-enumeration problem. That follow-up is outside this synthetic pilot.
