# QPA Rust Port Roadmap

Last updated: 2026-07-15

## Mission

Build a source-faithful Rust port of QPA. Complete and audit the Rust port first.
Python bindings and a `pip` distribution are a follow-up campaign and are not a
condition for declaring this Rust-porting goal finished.

This is a porting campaign, not a clean-room reimplementation and not a new
collection of algorithms that merely compares its answers with QPA. QPA's GAP
methods are the computational specification. GBNP and GAP are additional
specifications wherever QPA delegates work to them.

The current completion target is: every user-visible computation available in
the pinned QPA source is available from Rust, subject to the same mathematical
preconditions and with documented QPA-compatible conventions. A later Python
layer will expose that already-complete Rust surface without moving algorithms
out of Rust.

fd-applet is explicitly outside the scope of this campaign.

## Non-negotiable principles

- [x] Use a GPL-compatible license. The port will use
  `GPL-2.0-or-later`, matching QPA.
- [ ] Give every ported operation an explicit source mapping to the QPA method
  or GBNP/GAP primitive from which it was translated.
- [ ] Preserve QPA's mathematical conventions, method selection, important
  intermediate constructions, termination behavior, and result semantics.
- [ ] Preserve QPA's left-to-right path notation: `a * b` means first `a`,
  then `b`.
- [ ] Preserve the right-module and row-vector conventions used by QPA.
- [ ] Preserve distinctions such as equality, object identity, and module
  isomorphism at the API level where they are observable.
- [ ] Record every intentional deviation before merging it. A deviation must
  include its reason, affected operations, compatibility impact, and tests.
- [ ] Treat differential testing against GAP/QPA as verification of a source
  translation, not as a substitute for translating the source method.
- [ ] Keep every supported coefficient computation exact. Numerical and
  error-tracking floating-point types are not coefficient fields for this port,
  even if a dependency implements a field-like trait for them.
- [ ] Keep deterministic behavior by default. Any parallel or probabilistic
  implementation must have reproducible ordering/seed controls.
- [ ] Do not claim feature completeness from examples or a partial test suite;
  compute coverage from the complete declaration and installed-method
  inventories of the pinned sources.

Rust data structures do not need to imitate GAP records internally. They may use
typed indices, sparse storage, arenas, interning, and caches. Such changes are
allowed only when the QPA computation and observable semantics remain
traceable. Performance refactoring comes after a faithful reference path exists.

## Pinned upstream sources

The directories below are read-only computational references for development.
They are not intended to be bundled in Python wheels.

| Source | Local path | Revision | Clone | License | Role |
|---|---|---:|---|---|---|
| QPA | `upstream/qpa` | `8d5f847258d4853700ce3a5d947db2a1441109f1` | full | GPL-2.0-or-later | Primary computational and API specification |
| GBNP | `upstream/gbnp` | `e0b46c1aac16f3fbabe743e3169dc6bd86bdc654` | full | LGPL-2.1-or-later | `SGrobner` and its transitive algorithmic dependencies |
| GAP | `upstream/gap` | `4efe877859ab0bc0ac041bf635a996fb897ee0c9` | shallow | GPL-2.0-or-later | Field, matrix, basis, generic algebra, and dispatch semantics used transitively |

At the start of the campaign, QPA contains 39 `.gi` implementation files,
about 730 declarations, and about 932 installed functions/methods. These counts
are inventory seeds, not the final coverage denominator: overloads, global
functions, undocumented public names, generated methods, and relevant GAP
operations must also be classified.

GBNP's only declared package dependency other than GAP is GAPDoc. GAPDoc is not
an algorithmic dependency and is therefore not cloned. Only the transitive GBNP
surface required by QPA is part of QPA feature completeness; implementing every
standalone GBNP feature is a possible later project.

### Initial Rust algebra substrate

- [x] Select `numerica` as the initial exact coefficient backend and field
  contract.
- [ ] Pin an exact compatible `numerica` release in the Cargo lockfile and
  record its version, license, enabled features, and provenance in
  `upstream-lock.toml`.
- [ ] Initially admit the rational field and the prime fields supported by the
  selected `numerica` version. This is a family `GF(p)`, not a field hard-coded
  to one example prime.
- [ ] Exercise `GF(2)` explicitly because `numerica` represents characteristic
  two through a specialized path distinct from its odd-prime implementations.
- [ ] Treat a `numerica` type as supported only after it is admitted by the QPA
  coefficient adapter and passes the exact GAP-oracle contract suite. Merely
  implementing a `numerica` algebraic trait does not expose it automatically.
- [ ] Keep `numerica` types behind the public coefficient-domain boundary
  so future extension-field or accelerated backends do not require breaking the
  public QPA API.

Finite extension fields are deliberately deferred until computations over the
initial rational/prime-field family work end to end. Reusing `numerica` supplies
scalar and linear-algebra primitives; it does not replace any QPA or GBNP
algorithm.

### Updating an upstream baseline

- [ ] Never update an upstream reference and the port in one opaque change.
- [ ] Record old and new revisions in this file and the porting manifest.
- [ ] Generate a declaration/method diff before updating implementation code.
- [ ] Import upstream bug fixes as source-mapped changes with regression tests.
- [ ] Re-run the complete compatibility and benchmark suites.

## Planned repository structure

```text
.
├── Cargo.toml
├── LICENSE
├── ROADMAP.md
├── crates/
│   ├── qpa-core/          # objects, Numerica-backed exact fields, QPA methods
│   ├── qpa-gbnp/          # GBNP SGrobner translation used by QPA
│   └── qpa-python/        # PyO3 extension module
├── python/                # Deferred follow-up after Rust completion
│   └── qpa/
├── porting/
│   ├── qpa-operations.toml
│   ├── gbnp-functions.toml
│   ├── gap-primitives.toml
│   ├── deviations.md
│   └── upstream-lock.toml
├── tests/
│   ├── fixtures/          # stable serialized inputs and expected structures
│   ├── gap-reference/     # scripts run by the GAP/QPA reference harness
│   ├── rust/
│   └── python/
├── benchmarks/
├── docs/
└── upstream/
    ├── qpa/
    ├── gbnp/
    └── gap/
```

Crate and PyPI names are placeholders until availability and naming are checked.

## Porting ledger

The ledger is the authoritative measure of progress. It must be generated from
the pinned source and reviewed manually rather than maintained as an informal
list.

Each QPA declaration and applicable installed method needs:

- source repository, revision, file, and source span;
- GAP operation name and method filters/signature;
- whether it is public, internal, display-only, or a transitive primitive;
- dependencies on other QPA, GBNP, and GAP operations;
- Rust implementation and tests;
- Rust exposure and documentation;
- supported coefficient-domain rows;
- status: `inventory`, `specified`, `ported`, `verified`, or `complete`;
- any deviation identifier.

`complete` means all of the following, not merely “implemented”:

1. The source algorithm has been read and mapped.
2. The Rust implementation follows that algorithm.
3. Source-level invariants and edge cases have tests.
4. Relevant upstream tests/examples pass after translation.
5. GAP/QPA differential tests pass on the supported field matrix.
6. The operation is available in the public Rust API.
7. Public semantics and limitations are documented.
8. No unresolved correctness or compatibility issue is attached to the entry.

## Cross-cutting compatibility specifications

These specifications must be written before large-scale implementation. They
prevent many locally reasonable choices from becoming globally incompatible.

### Object and convention specification

- [ ] Vertex, arrow, path, and zero-path identity and equality.
- [ ] Source/target rules and zero for non-composable path products.
- [ ] Generator numbering and external representations.
- [ ] Default length-left-lexicographic ordering and all alternative orderings.
- [ ] Coefficient/path term ordering and canonical sparse representation.
- [ ] Quotient element normalization and cache invalidation.
- [ ] Right modules, row vectors, and arrow-matrix orientation.
- [ ] Zero-dimensional vector-space and matrix-shape conventions.
- [ ] Morphism composition order.
- [ ] GAP `fail`, errors, infinity, unknown properties, and nontermination mapped
  to explicit Rust results.
- [ ] Cached GAP attributes/properties mapped to safe Rust cache semantics.

### Coefficient-domain matrix

The scalar layer must be extensible from the beginning. Feature coverage and
coefficient coverage are tracked separately so that neither can hide holes in
the other.

- [ ] Tier A: the rational field and every exact prime field admitted from the
  pinned `numerica` backend, exposed uniformly as `QQ` and `GF(p)`.
- [ ] Tier B: finite extension fields.
- [ ] Tier C: cyclotomic/algebraic number fields and rational function fields
  used by QPA workflows.
- [ ] Tier D: every remaining division ring or field accepted by the pinned QPA
  methods, or a source-supported explanation that QPA itself rejects it.
- [ ] Stable serialization including field construction data.
- [ ] Validate `GF(p)` construction and reject composite or unsupported moduli
  with stable Rust errors.
- [ ] Exhaustively test arithmetic in small prime fields, including `GF(2)`, and
  use deterministic property/differential tests across representative small,
  machine-word, 64-bit, and arbitrary-size primes that the backend admits.
- [ ] Characteristic-sensitive tests for every algorithm that branches on it.

The first useful release may support early tiers, but the campaign is not
coefficient-complete until the final tier is audited.

### GAP primitive boundary

We will not port the GAP interpreter or its entire algebra library. We will port
or provide equivalent exact implementations only for the transitive primitives
used by QPA/GBNP, while preserving the semantics on which their source methods
depend.

- [ ] Generate the transitive primitive inventory.
- [ ] Classify primitives as scalar, matrix, vector-space/basis, polynomial,
  generic algebra, graph/list, randomness, or method-dispatch behavior.
- [ ] Map `NullspaceMat`, `SolutionMat`, echelon/basis functions, matrix inverse,
  characteristic polynomial, radicals, ideals, and related calls explicitly.
- [ ] Record where `numerica` or another external Rust crate is used instead of
  translated GAP primitive code and verify the required semantics against the
  pinned GAP oracle.
- [ ] Do not let a Rust dependency introduce floating-point arithmetic into an
  exact computation.

## Verification strategy

Every phase uses several kinds of evidence.

### 1. Source-translation tests

Small tests target each branch, loop invariant, normalization rule, and special
case visible in the source. These are the primary evidence that the method was
translated rather than reinvented.

### 2. Imported upstream tests and examples

- [ ] Convert all commands in `upstream/qpa/tst/testall.tst`.
- [ ] Convert every file in `upstream/qpa/examples/`.
- [ ] Convert every documented example in the QPA manual.
- [ ] Import the GBNP tests needed by the `SGrobner` dependency closure.
- [ ] Preserve regression tests associated with upstream issues and bug fixes.

### 3. Differential reference harness

Run the same serialized case through pinned GAP/QPA and Rust, then compare
canonical mathematical structures rather than fragile pretty-printed text.

- [ ] Make the canonical runner load QPA and GBNP from `upstream/qpa` and
  `upstream/gbnp`; reject an accidentally loaded global package version.
- [ ] Pin and report the GAP, QPA, and GBNP versions/revisions in every oracle
  result. A fast local runner may use a compatible system GAP, but canonical CI
  fixtures must use the declared GAP baseline.
- [ ] Define a versioned, language-neutral fixture schema for fields,
  coefficients, shaped matrices, quivers, paths, sparse elements, modules, and
  morphisms.
- [ ] Serialize rationals as reduced integer pairs and prime-field elements as
  canonical residues together with their modulus; never compare decimal
  approximations.
- [ ] Check generated oracle fixtures into the repository with their input,
  upstream revision metadata, and regeneration command. Rust unit tests
  must remain runnable without GAP installed.
- [ ] Compare exact canonical forms when an operation has a canonical result.
  For basis-dependent or otherwise noncanonical output, additionally compare
  the defining mathematical structure (for example mutual ideal containment
  and agreement of normal forms), not only dimensions or ranks.
- [ ] Compare stable error categories and enforce explicit resource limits for
  cases that may fail, time out, or not terminate.

- [ ] Quivers, paths, and orderings.
- [ ] Sparse elements, relations, Gröbner bases, tips, nontips, and normal forms.
- [ ] Bases, dimensions, matrices, modules, and morphisms.
- [ ] Resolutions, complexes, invariants, and specialized outputs.
- [ ] Errors, precondition failures, and bounded nonterminating calculations.

### 4. Property and metamorphic tests

Use mathematical identities that QPA's methods must preserve: associativity,
normal-form idempotence, relation annihilation, functoriality, exactness,
kernel/image identities, chain conditions, and isomorphism invariance.

### 5. Performance tests

- [ ] Capture representative QPA workloads before optimizing.
- [ ] Measure time, peak memory, allocation counts, and output size.
- [ ] Keep small, medium, and research-scale benchmark sets.
- [ ] Require correctness parity before accepting a faster implementation.
- [ ] Record deterministic performance comparisons against pinned GAP/QPA.

## Campaign phases

Phases describe dependency order, not a reduction of the final Rust scope.
Python bindings are deliberately deferred until the Rust completion audit.

### Phase 0 — Campaign infrastructure and executable inventory

Source families: all QPA declarations/installations; QPA package loader; the
transitive GBNP/GAP call boundary.

- [x] Clone and pin QPA.
- [x] Clone and pin GBNP.
- [x] Clone and pin GAP as a shallow computational reference.
- [x] Confirm GPL licensing is acceptable.
- [ ] Add `LICENSE` and copyright/provenance notices.
- [x] Create the Cargo workspace and foundational crates.
- [ ] Generate `qpa-operations.toml` from all `.gd` and `.gi` files.
- [ ] Generate the GBNP transitive-function ledger rooted at `SGrobner`.
- [ ] Generate the GAP primitive ledger rooted at QPA plus that GBNP closure.
- [ ] Record all upstream revisions in `upstream-lock.toml`.
- [x] Create a pinned GAP/QPA reference runner.
- [x] Make the runner emit and consume the versioned structured oracle schema.
- [x] Add a guard test proving that the runner loaded the cloned QPA/GBNP
  packages rather than globally installed copies.
- [x] Run the unmodified upstream QPA test baseline under the pinned runner.
- [ ] Establish Rust CI on Linux, macOS, and Windows.
- [ ] Establish formatting, linting, documentation, coverage, and license gates.
- [x] Select and pin a Rust version compatible with Numerica.

Exit gate: every pinned QPA declaration and method has a stable ledger entry;
the reference runner and Rust packages run in CI.

### Phase 1 — Exact scalar and linear-algebra substrate

GAP boundary: fields, vectors, matrices, row spaces, bases, subspaces, quotient
spaces, linear maps, polynomial invariants, and deterministic randomness.

- [x] Pin and integrate `numerica` behind the coefficient-domain boundary.
- [ ] Define runtime descriptors and adapters for `QQ` and `GF(p)` over the
  exact `numerica` field contract.
- [ ] Implement Tier A coefficient domains, including the specialized `GF(2)`
  path and the admitted odd-prime backends.
- [ ] Use `numerica` exact matrix operations where their semantics match GAP;
  add compatibility wrappers or source-translated operations where orientation,
  normalization, shape, or error behavior differs.
- [x] Support dense exact matrices with explicit zero-dimensional shapes.
- [x] Implement the GAP-compatible row-vector convention.
- [ ] Port/equate echelonization, rank, nullspace, solution, inverse, image,
  kernel, and coordinate computations.
- [ ] Implement vector spaces, bases, subspaces, quotient spaces, and direct sums.
- [ ] Implement characteristic/minimal-polynomial operations required by QPA.
- [ ] Add deterministic RNG plumbing for probabilistic QPA methods.
- [ ] Differentially test rational arithmetic and representative prime fields,
  with exhaustive small-field tests and explicit characteristic-two cases.
- [ ] Differentially test shaped matrix arithmetic, echelon forms, rank,
  `NullspaceMat`, `SolutionMat`, inverse, kernel, and image against GAP.

Exit gate: every scalar and linear-algebra primitive needed by the next phases
has exact GAP-reference fixtures over Tier A fields.

### Phase 2 — Quivers, paths, orderings, posets, and graph algorithms

Primary QPA files: `quiver`, `groebnerordering`, `qpaposet`, and
`quivergraphalgorithms`.

- [x] Port vertices, arrows, zero paths, paths, walks, sources, and targets.
- [x] Port quiver construction from arrow lists and adjacency matrices.
- [x] Port path multiplication, equality, enumeration, and external forms.
- [x] Port incoming/outgoing arrows, adjacency, components, and graph properties.
- [x] Port opposite, double, separated, full-subquiver, and Dynkin constructions.
- [ ] Port all QPA monomial/path orderings and ordering composition.
- [x] Port quiver enumerators without materializing infinite path sets.
- [ ] Port QPA posets and remaining quiver graph algorithms.
- [ ] Provide Rust constructors that preserve every QPA ordering composition.

Exit gate: all quiver/path/ordering declarations and examples are complete over
the ledger, including cyclic quivers and zero/non-composable paths.

### Phase 3 — Sparse path algebras, ideals, and quotient objects

Primary QPA files: `pathalg` and `pathalgideal` excluding later high-level
properties that depend on modules or specialized algebra constructions.

- [ ] Implement sparse linear combinations of paths with canonical term order.
- [ ] Implement path rings/algebras, generators, zero, one, and element actions.
- [ ] Implement left/right uniformity and source/target decomposition.
- [ ] Implement ideals, relation lists, quotient construction, and natural maps.
- [ ] Implement quotient elements with lazy/eager normal-form state.
- [ ] Implement basis, coefficient, membership, and element conversion hooks.
- [ ] Implement admissibility, monomial, quadratic, and arrow-ideal operations.
- [ ] Preserve QPA external representation and serialization semantics.

Exit gate: non-Gröbner path-algebra operations are source mapped and verified;
normal-form-dependent operations are connected to an explicit pending backend.

### Phase 4 — QPA and GBNP Gröbner machinery

Primary QPA files: `gbnp`, `groebner`, `gbhighlevel`, `invlenlexgrb`,
`groebnerdyntree`, `groebnerdyndict`, and `groebnerstadict`.

GBNP closure begins with `SGrobner` and includes NP arithmetic/format,
strong-normal-form reduction, obstruction generation, occurrence trees, heaps,
options, and the Gröbner loop. Extra GBNP APIs are not automatically in scope.

- [ ] Port QPA-to-Cohen/NP conversion and uniform decomposition exactly.
- [ ] Port GBNP noncommutative polynomial representation and arithmetic.
- [ ] Port leading monomials/terms, cleaning, monic normalization, and reduction.
- [ ] Port occurrence trees and the k-ary todo heap.
- [ ] Port obstruction/overlap construction and strong normal forms.
- [ ] Port the four phases of `SGrobner`, including resumable iteration bounds.
- [ ] Port QPA conversion back from GBNP and path-algebra structure relations.
- [ ] Port QPA Gröbner objects, tip reduction, complete reduction, and properties.
- [ ] Port nontip enumeration and finite-dimensionality detection.
- [ ] Port QPA right Gröbner bases.
- [ ] Port `HighLevelGroebnerBasis` and its reduction operations.
- [ ] Port inverse-length/inverse-lexicographic Gröbner operations.
- [ ] Port static/dynamic dictionaries, suffix trees, and dynamic trees used by
  QPA's implementations.
- [ ] Define cancellation, progress, iteration, memory, and user-interrupt APIs
  for computations that may not terminate.

Exit gate: every QPA Gröbner backend gives source-compatible canonical results
on the complete relevant QPA/GBNP corpus, and bounded computations can resume.

### Phase 5 — Finite-dimensional algebras and algebra constructions

Primary QPA files: `algebra`, the remaining `pathalg`, `idempotent`,
`pathalgpredef`, `pathalgtensor`, `pathalgkoszul`, and `pathalgreptype`.

- [ ] Complete quotient bases, dimensions, multiplication, radicals, and centers.
- [ ] Port opposite algebras and homomorphisms.
- [ ] Port primitive-idempotent and lifting algorithms.
- [ ] Port Cartan/Coxeter, Nakayama, Frobenius, Loewy, global-dimension, and
  structural property computations.
- [ ] Port associated monomial/graded algebras and corner/quotient constructions.
- [ ] Port tensor products, enveloping algebras, and algebra-as-bimodule.
- [ ] Port trivial extensions and their structural maps.
- [ ] Port quadratic dual/Koszul helpers.
- [ ] Port all predefined QPA algebra families and their recognition properties.
- [ ] Port representation-type and finite-global-dimension tests in this layer.
- [ ] Port save/read using a safe, versioned, language-independent format while
  retaining a QPA-compatible import path.

Exit gate: all algebra-level QPA declarations are complete for Tier A fields;
every predefined algebra example and property test is translated.

### Phase 6 — Matrix right modules and special modules

Primary QPA files: `module` and `modulespecial`.

- [ ] Port module elements and the right action of paths/algebra elements.
- [ ] Port matrix-module construction with and without relation checking.
- [ ] Validate dimensions, fields, relation satisfaction, and zero spaces as QPA
  does.
- [ ] Port bases, canonical bases, coordinates, dimension vectors, and support.
- [ ] Port submodules, restriction, direct sums, generators, radical/socle series,
  Loewy length, annihilators, and module predicates.
- [ ] Port simples, indecomposable projectives/injectives, and zero modules.
- [ ] Port projective bases and element conversions.
- [ ] Port rigidity, tau-rigidity, exceptionality, and complexity operations.
- [ ] Add bulk NumPy-style import/export without changing exact arithmetic.

Exit gate: every module constructor and module-level operation in the ledger is
available in Rust and matches QPA on translated module examples.

### Phase 7 — Module homomorphisms and decomposition

Primary QPA files: `modulehom` and `moduledecomp`.

- [ ] Port matrix-defined module maps and commutativity validation.
- [ ] Port map equality, addition, scalar multiplication, and composition.
- [ ] Port `HomOverAlgebra` using QPA's linear-equation construction.
- [ ] Port endomorphism algebras and basis-coordinate functions.
- [ ] Port kernel, image, cokernel, induced maps, intersections, and sums.
- [ ] Port radical/top/socle maps and direct-sum inclusions/projections.
- [ ] Port liftings, split maps, inverses, minimal versions, and factoring through
  projectives.
- [ ] Port isomorphism/direct-summand/additive-closure computations.
- [ ] Port traces, rejects, submodule enumeration, and length-bounded module
  generation.
- [ ] Port endomorphism quiver-algebra and ADR-algebra constructions.
- [ ] Port deterministic and probabilistic module decomposition methods,
  including seed control and idempotent lifting.

Exit gate: the full homomorphism and decomposition ledgers are complete, and
all categorical identities used by later phases pass property tests.

### Phase 8 — Projective path-algebra modules, presentations, and resolutions

Primary QPA files: `present2`, `projpathalgmodule`, and `moduleprojres`.

- [ ] Port path-algebra vectors and projective modules.
- [ ] Port right-module Gröbner bases and module reduction/division operations.
- [ ] Port projective presentations and maps between projectives.
- [ ] Port QPA's presentation utilities from `present2`.
- [ ] Port the Green–Solberg–Zacharia projective-resolution construction.
- [ ] Port `XSet`, `OSet`, `NSet`, `FirstPart`, and all supporting algorithms.
- [ ] Port finite and bounded resolution representations, maps, and accessors.
- [ ] Test monomial and nonmonomial relations separately.

Exit gate: projective presentations and each QPA resolution implementation
match source-level intermediate data as well as final differentials.

### Phase 9 — Homological algebra, functors, and AR theory

Primary QPA files: `functors`, `modulehomalg`, and `moduleARtheory`.

- [ ] Port projective covers, injective envelopes, syzygies, and cosyzygies.
- [ ] Port `ExtOverAlgebra` and all QPA extension constructions.
- [ ] Port pushouts, pullbacks, approximations, and perpendicular-category tools.
- [ ] Port tilting, cotilting, Gorenstein, dominant-dimension, and periodicity
  computations present in QPA.
- [ ] Port dual, transpose, `DTr`, `TrD`, star, Nakayama, and tensor functors on
  objects and morphisms.
- [ ] Port almost-split sequences, irreducible morphisms, predecessors, and
  related AR computations.
- [ ] Verify functoriality and exact-sequence identities, not only dimensions.

Exit gate: all homological/functor/AR operations in the pinned QPA ledger are
complete over the supported coefficient tiers.

### Phase 10 — Infinite lists, complexes, chain maps, and derived computations

Primary QPA files: `inflist`, `complex`, and `derivedcat`.

- [ ] Port lazy bi-infinite/half-infinite/repeating list semantics and caching.
- [ ] Port finite, bounded, repeating, projective, and injective complexes.
- [ ] Port cycles, boundaries, homology, exactness, shifts, and truncations.
- [ ] Port chain maps, comparison liftings, mapping cones, and Yoneda products.
- [ ] Port projective resolutions of complexes and differential reduction.
- [ ] Port derived-category tau and projective/injective complex conversions.
- [ ] Verify all degree/sign conventions with explicit low-degree fixtures.

Exit gate: all complex and derived-category declarations are complete, including
lazy behavior and every documented sign convention.

### Phase 11 — Specialized and combinatorial representation theory

Primary QPA files: `stringalgebra`, `validity`, `combinatorialrep`,
`combinatorialmap`, and `degorderfinitetype`.

- [ ] Port string/band syntax, validity, inversion, enumeration, modules, and
  local AR operations.
- [ ] Port special-biserial, string, gentle, and domestic recognition methods
  that are part of QPA itself.
- [ ] Port unit forms, roots, Euler/Tits forms, and reflections.
- [ ] Port combinatorial maps, genus/faces/boundaries, cutting/joining, homology,
  and gentle-derived-equivalence computations.
- [ ] Port finite-type AR-quiver data and degeneration-order algorithms.
- [ ] Port all remaining utility functions required by these files.

Exit gate: the five specialized source families have no incomplete ledger
entries and all manual examples are executable from Rust.

### Phase 12 — Rust API completion and documentation

This phase closes Rust usability and distribution gaps across the whole system.

- [ ] Expose every complete ledger entry in the public Rust API.
- [ ] Define exceptions for invalid mathematical input, unsupported domains,
  cancellation, resource limits, and nontermination bounds.
- [ ] Provide generated Rust API reference documentation.
- [ ] Translate the QPA quickstart and manual examples into Rust.
- [ ] Write a GAP-to-Rust migration and convention guide.
- [ ] Build reproducible Rust source releases.
- [ ] Test Rust installation in clean environments with no GAP installation.
- [ ] Add source/license notices required by QPA, GBNP, GAP, and Rust dependencies.

Exit gate: a clean user can build the Rust crates, follow the migrated
documentation, and access the entire completed operation ledger without GAP.

### Phase 13 — Optimization without algorithm drift

- [ ] Establish the faithful implementation as the reference mode.
- [ ] Profile real QPA workloads before changing representations or control flow.
- [ ] Optimize path interning, sparse terms, normal-form indexes, matrices, and
  allocation hot spots.
- [ ] Add safe parallelism only where operation ordering does not change results,
  or canonicalize results before exposure.
- [ ] Compare each optimized path against the faithful path on randomized and
  adversarial corpora.
- [ ] Document performance deviations and resource-limit behavior.

Exit gate: accepted optimizations have equivalence evidence and no operation
loses source traceability.

### Phase 14 — Feature-completeness audit and 1.0 release

- [ ] Regenerate the operation/method inventory from the pinned QPA revision.
- [ ] Confirm every applicable entry is `complete`; no unsupported or silently
  omitted QPA operation remains.
- [ ] Confirm every transitive GBNP/GAP primitive has a verified implementation.
- [ ] Run all translated tests, examples, manual sessions, property suites, and
  differential suites over the declared coefficient matrix.
- [ ] Audit error behavior, cancellation, bounded nontermination, serialization,
  and cross-platform determinism.
- [ ] Audit Rust API/documentation coverage against the same ledger.
- [ ] Audit GPL/LGPL notices and source-distribution obligations.
- [ ] Publish release candidates and have representation theorists reproduce
  research computations independently.
- [ ] Freeze the 1.0 compatibility contract and publish Rust source artifacts.

Exit gate: the generated audit proves the mission-level definition of “everything
in QPA is available,” rather than relying on a manually curated feature list.

## Source-file coverage map

Every current QPA implementation file is assigned a primary phase below. A file
may have dependencies in earlier phases, but it must not disappear between
phase descriptions.

| Phase | QPA implementation files |
|---|---|
| 1 | `util.gi` and the transitive GAP primitive layer |
| 2 | `quiver.gi`, `groebnerordering.gi`, `qpaposet.gi`, `quivergraphalgorithms.gi` |
| 3 | `pathalg.gi`, `pathalgideal.gi` |
| 4 | `gbnp.gi`, `groebner.gi`, `gbhighlevel.gi`, `invlenlexgrb.gi`, `groebnerdyntree.gi`, `groebnerdyndict.gi`, `groebnerstadict.gi` |
| 5 | `algebra.gi`, `idempotent.gi`, `pathalgpredef.gi`, `pathalgtensor.gi`, `pathalgkoszul.gi`, `pathalgreptype.gi` |
| 6 | `module.gi`, `modulespecial.gi` |
| 7 | `modulehom.gi`, `moduledecomp.gi` |
| 8 | `present2.gi`, `projpathalgmodule.gi`, `moduleprojres.gi` |
| 9 | `functors.gi`, `modulehomalg.gi`, `moduleARtheory.gi` |
| 10 | `inflist.gi`, `complex.gi`, `derivedcat.gi` |
| 11 | `stringalgebra.gi`, `validity.gi`, `combinatorialrep.gi`, `combinatorialmap.gi`, `degorderfinitetype.gi` |
| 12 | `documentation.gi` |

The generated ledger, not this table, will handle individual methods installed
on GAP-generic operations or methods whose best implementation phase differs
from their file's primary phase.

## Release milestones

These milestones are progress labels, not changes to the final scope.

- **M0 — Porting laboratory:** executable inventories, reference runner, CI,
  exact Tier A linear algebra, and package skeleton.
- **M1 — Algebra kernel:** quivers, paths, path algebras, ideals, all QPA
  Gröbner backends, quotients, and finite-dimensional algebra basics.
- **M2 — Module kernel:** matrix modules, morphisms, projectives, presentations,
  decompositions, and resolutions.
- **M3 — Homological system:** functors, Ext/syzygies, AR theory, complexes, and
  derived operations.
- **M4 — Full specialist surface:** algebra constructions, string/gentle and
  combinatorial functionality, broad field tiers, and Rust documentation.
- **M5 / 1.0 — Audited QPA completeness:** Phase 14 passes with no missing
  ledger entries.

## Immediate next actions

- [ ] Complete the source/license/provenance files for the Rust workspace.
- [ ] Write the inventory extractor and commit the initial QPA operation ledger.
- [ ] Compute the exact GBNP call closure rooted at `SGrobner`.
- [x] Build the GAP/QPA reference test harness using the pinned package clones.
- [ ] Extend the structured oracle schema with rational/prime-field
  scalar and matrix fixtures.
- [x] Run the baseline QPA suite with zero failures under the pinned oracle.
- [ ] Write the object/convention and GAP primitive-boundary specifications.
- [x] Pin `numerica` and begin Phase 1 with `QQ`, `GF(p)`, and exact
  GAP-compatible row-reduction fixtures.
