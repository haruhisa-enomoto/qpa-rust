# Intentional compatibility deviations

This ledger records places where Rust's public boundary intentionally differs
from the pinned GAP/QPA surface. It does not authorize algorithmic deviations:
the computations remain translations of the cited source methods.

## ALG-001 — Finite opposite-algebra display repairs recursive QPA method

- Source: QPA `lib/algebra.gi`, `PrintObj` at lines 109–115.
- Affected Rust surface: `Display` for
  `OppositeStructureConstantsAlgebra`.
- Difference: the QPA method prints `A` while already printing `A`, recursively
  invoking itself until GAP stops at its recursion limit. Rust emits the finite
  form `OppositeAlgebra( <structure-constants algebra of dimension n> )`.
- Reason: the source evidently intended to print the underlying algebra, but
  its self-reference cannot produce a usable representation. Rust preserves
  the wrapper and a deterministic description without recursion.
- Compatibility impact: valid computations are unaffected. Code that captures
  the accidental GAP recursion diagnostic instead receives a finite string.
  Opposite elements retain QPA's exact `op( v.i )` form for canonical basis
  elements.
- Evidence: `tests/gap-reference/algebra_opposite_print.g` pins the source
  recursion, `tests/fixtures/algebra_opposite.schema-v1.json` pins element
  printing, and the structure-constant tests pin Rust's finite algebra display.

## ALG-002 — Explicit named lookup replaces GAP record-component dispatch

- Source: QPA `lib/algebra.gi`, the opposite-algebra `\.` method at lines
  189–199.
- Affected Rust operations: `StructureConstantsAlgebra::bind_named_element`
  and `OppositeStructureConstantsAlgebra::element_by_name`.
- Difference: GAP converts an integer record-name token with `NameRNam`, asks
  the underlying algebra for that dynamic component, and raises a record error
  when it is absent. Rust binds names explicitly as strings and returns
  `Ok(None)` for an absent name; a present element is reconstructed in the
  opposite family exactly as QPA does.
- Reason: Rust has neither GAP's global record-name integers nor dynamic record
  component syntax. An owner-checked name table makes the same lookup explicit
  and prevents components from containing an element of another algebra.
- Compatibility impact: successful named lookup has the same coordinates and
  opposite parent. Rust callers handle absence as an option instead of a GAP
  interpreter error.
- Evidence: `tests/gap-reference/algebra_opposite_dot_error.g` pins QPA's
  absent-component failure at `algebra.gi:198`, and the opposite-algebra fixture
  test covers successful wrapping plus the optional Rust absence result.

## ALG-003 — Prime-field linear-mapping radical repairs GAP basis failure

- Source: QPA `lib/algebra.gi`, `RadicalOfAlgebra` at lines 288–320.
- Affected Rust operation:
  `ModuleEndomorphismAlgebra::radical_basis_matrices`.
- Difference: QPA's transpose-to-matrix-algebra adapter succeeds over the
  rationals, but on the same dual-number mapping algebra over `GF(2)` the pinned
  GAP basis machinery passes `fail` to `IsZero` and aborts. Rust supports both
  characteristic zero and Numerica prime fields.
- Reason: Rust already translates the exact GAP structure-constant radical
  algorithms: trace-form nullspace in characteristic zero and the
  Cohen--Ivanyos--Wales ideal sequence in positive characteristic. The mapping
  adapter constructs that multiplication table from its retained matrix basis
  and lifts the resulting radical rows back, which is the computation QPA's
  transposition adapter requests without relying on the failing GAP dispatch.
- Compatibility impact: rational results and basis order agree with QPA. A
  supported prime-field call succeeds in Rust where pinned QPA/GAP aborts; no
  extension-field behavior is claimed.
- Evidence: `tests/fixtures/algebra_mapping_radical.schema-v1.json` pins the
  successful rational result,
  `tests/gap-reference/algebra_mapping_radical_f2_error.g` pins the upstream
  finite-field failure, and the path-module tests verify the Rust radical over
  both `QQ` and `GF(2)`.

## DT-001 — Explicit dynamic-tree ownership arena

- Source: `upstream/qpa/lib/groebnerdyntree.gd` and
  `upstream/qpa/lib/groebnerdyntree.gi`.
- Affected Rust operations: `DynamicForest::create` and `DynamicTree::link`.
- Difference: GAP independently allocates garbage-collected dynamic-tree
  records. Rust requires nodes that may be linked to come from one explicit
  `DynamicForest`; linking nodes from different forests returns
  `DynamicTreeError::DifferentForests`.
- Reason: QPA's mutable parent/left/right graph is cyclic. Arena indices retain
  stable object identity without `unsafe`, dangling weak references, or leaked
  `Rc` cycles.
- Compatibility impact: rotations, splaying, roots, linking, cutting, stored
  items, and identity agree within a forest. QPA's suffix-tree implementation
  owns one forest internally, so its user-visible dictionary behavior is
  unaffected. Direct Rust callers must allocate related nodes together.
- Evidence: `tests/gap-reference/dynamic_tree.g`,
  `tests/fixtures/dynamic_tree.schema-v1.json`, and the tests in
  `crates/qpa-core/src/dynamic_tree.rs`.

## SUFFIX-001 — Owner-tagged suffix-tree arena

- Source: `upstream/qpa/lib/groebnerdyndict.gd` and
  `upstream/qpa/lib/groebnerdyndict.gi`.
- Affected Rust surface: `SuffixTreeNodeId`, `SuffixTreeEdgeId`,
  `SuffixTree::create_internal_node`, `SuffixTree::create_external_node`,
  `SuffixTree::create_edge`, and the node/edge accessors and mutations.
- Difference: GAP can allocate free-standing garbage-collected node and edge
  records from an arbitrary compatible family and pattern descriptor. Rust
  allocates them through one `SuffixTree` arena, requires low-level labels to
  reference an existing retained pattern slot, and tags every handle with its
  owner. A handle used with another tree returns `SuffixTreeError::UnknownNode`
  or `SuffixTreeError::UnknownEdge`.
- Reason: QPA's parent edges, child edges, suffix links, and dynamic-tree
  handles form a mutable cyclic graph. Owner-tagged arena identities preserve
  QPA's object-identity equality without `unsafe`, dangling references, or
  accidental equality between the same arena index in independent trees.
- Compatibility impact: construction, labels, ordered edge lookup, edge
  splitting, parent identity, online insertion, deletion/compression, and both
  pattern queries agree inside one tree. Direct Rust callers create related
  low-level objects through their owning tree instead of assembling unrelated
  GAP records.
- Evidence: `tests/gap-reference/suffix_tree_surface.g`,
  `tests/fixtures/suffix_tree_surface.schema-v1.json`, and the suffix-tree tests
  in `crates/qpa-core/src/suffix_tree.rs`.

## PATHALG-001 — Natural projection application repairs stale GAP record access

- Source: QPA `lib/pathalg.gi`, `ImagesRepresentative` at lines 769–787 and
  `NaturalHomomorphismByIdeal` at lines 312–336.
- Reason: QPA 1.37 applies a path-algebra homomorphism through the private
  fields `alghom!.generators` and `alghom!.genimages`. The pinned GAP 4.17dev
  mapping representation no longer binds `generators`, so applying the natural
  projection raises a record-element error even though construction, source,
  range, generator images, and surjectivity are valid.
- Rust behavior: `NaturalQuotientMap::image` translates the intended
  `MappedExpression`: each vertex/arrow generator maps to the corresponding
  quotient generator, hence a sparse path-algebra representative maps to the
  identical sparse representative in the factor family. It preserves QPA's
  lazy state when the ideal has no attached Gröbner basis and reduces only when
  QPA would have installed a normal-form function.
- Compatibility impact: the Rust operation succeeds where the pinned
  QPA/GAP combination errors because of obsolete representation access; its
  mathematical result is the map specified by QPA's constructor.
- Tests: `tests/gap-reference/natural_map_image_error.g` pins the upstream
  failure; `pathalg_foundational.schema-v1.json` pins the source/range,
  surjectivity, lazy and attached-basis quotient states; Rust tests exercise the
  repaired projection over `QQ` and `GF(2)`.

## ORDERING-001 — Quiver-ordering ownership is checked at construction

- Source: QPA `lib/groebnerordering.gi`, all ordering constructors,
  `PathFromArrowList`, and the ordering-subset operations.
- Difference: QPA's GAP filters require an ordering or homogeneous generator
  list but do not require it to belong to the constructor's quiver. It can
  therefore create a cross-family ordering or formal path whose position tables
  later fail, or accidentally compute when unrelated families happen to have
  compatible generator positions. Rust requires the quiver, next/component
  orderings, assigned/preferred generators, paths, and subset generators to
  have one identity and returns `OrderingError::DifferentQuivers` or
  `PathError::DifferentQuivers` otherwise.
- Reason: generator positions are family-relative. Rejecting a malformed owner
  graph at its typed boundary prevents delayed out-of-bounds access without
  altering any QPA computation on well-owned data.
- Compatibility impact: every ordering composition on one quiver follows QPA;
  deliberately mixing independent quiver families is rejected earlier. The
  pinned oracle records that GAP accepts foreign next orderings and rejects
  heterogeneous lists only through method dispatch.
- Evidence: `tests/gap-reference/ordering_dispatch.g`,
  `tests/fixtures/ordering_dispatch.schema-v1.json`, and the ordering dispatch
  tests in `crates/qpa-core/src/ordering.rs`.

## INVLENLEX-001 — Typed zero sentinel for inverse tips and walks

- Source: QPA `lib/invlenlexgrb.gi`, `TipMonomialInvLenInvLex` and
  `TipWalkInvLenInvLex` at lines 41–65.
- Affected Rust operations: `tip_monomial_inv_len_inv_lex` and
  `tip_walk_inv_len_inv_lex`.
- Difference: on the zero input, QPA's monomial method returns `Zero(e)`, a
  path-algebra element rather than a path. Passing that result to `WalkOfPath`
  makes the walk method fail by GAP method dispatch. Rust returns `None` for
  both optional results; nonzero monomials and walks retain the exact path and
  arrow order.
- Reason: a zero polynomial has no leading monomial or walk. `Option` exposes
  that source sentinel without mixing incompatible Rust representations or
  reproducing an interpreter dispatch failure.
- Compatibility impact: callers match `None` instead of testing a zero
  algebra element or catching the subsequent GAP error. Every nonzero result
  is unchanged.
- Evidence: `tests/gap-reference/invlenlex.g` and
  `tests/fixtures/invlenlex.schema-v1.json` pin the zero monomial sentinel and
  zero-walk failure alongside the Rust optional results.

## GROEBNER-001 — Owned transformed values replace standalone GAP mutation

- Source: QPA `lib/groebner.gi`, `TipReduceGroebnerBasis` at lines 198–352
  and `CompletelyReduceGroebnerBasis` at lines 58–82.
- Affected Rust operations: `GroebnerBasis::tip_reduced_basis` and
  `GroebnerBasis::completely_reduced_basis`.
- Difference: QPA rewrites the relation component and cached properties of the
  identical mutable basis object. Rust returns an owned transformed clone and
  leaves the input value unchanged. The ideal-level
  `PathAlgebraIdeal::right_groebner_basis` explicitly writes the completely
  reduced value back to the ideal cache, preserving QPA's observable cached
  state at that public boundary.
- Reason: an owned value transformation avoids alias-dependent mutation while
  retaining exact relation order, property states, and source failures. Callers
  must use the returned basis, as they do for other fallible Rust transforms.
- Compatibility impact: computations and cached ideal results agree with QPA;
  a separately retained standalone Rust basis handle does not change behind an
  alias when another handle is reduced.
- Evidence: `tests/gap-reference/groebner_objects.g`,
  `tests/fixtures/groebner_objects.schema-v1.json`, and the Gröbner object tests
  in `crates/qpa-core/src/groebner.rs` pin transformations, flags, iterators,
  setters, right-basis construction, and the ideal cache after reduction.

## IDEMPOTENT-001 — Explicit reproducible randomness replaces GAP global state

- Source: QPA `lib/idempotent.gi`, `HOPF_ZeroDivisor` at lines 116–141,
  `HOPF_SingularIdempotent` at lines 198–215, and `PrimitiveIdempotents` at
  lines 21–67; GAP `lib/polyfinf.gi`, finite-field `Factors` and
  `FactorsCommonDegreePol`.
- Affected Rust operations:
  `StructureConstantsAlgebra::hopf_zero_divisor`,
  `StructureConstantsAlgebra::hopf_singular_idempotent`, and
  `StructureConstantsAlgebra::primitive_idempotents`.
- Difference: QPA and its finite-field factorization consume GAP's mutable
  global random source. Rust's `_with_rng` variants consume an explicit caller
  RNG, while the convenience entry points start `StdRng` from
  `HOPF_DEFAULT_SEED` on each uncached computation.
- Reason: the Rust API must make probabilistic state controllable and provide a
  reproducible default. The translated algorithm still samples field
  coefficients, computes the same minimum polynomial, applies GAP's
  square-free/distinct-degree/equal-degree factorization stages, and performs
  QPA's identical recursive corner construction.
- Compatibility impact: returned zero divisors, singular idempotents, and
  complete primitive systems can differ from a GAP session's particular random
  representatives. Their defining products, right-ideal dimensions,
  completeness, orthogonality, factorization stages, and failure conditions
  agree.
- Evidence: `tests/gap-reference/idempotent.g`,
  `tests/fixtures/idempotent.schema-v1.json`, and the idempotent tests in
  `crates/qpa-core/src/idempotent.rs` verify the structural result over
  `GF(2)` and odd prime fields, including `M_2` and a recursive `M_3` corner.

## IDEMPOTENT-002 — Attribute coordinates are cached behind owned Rust results

- Source: QPA `lib/idempotent.gd`, generated `PrimitiveIdempotents` attribute,
  and `lib/idempotent.gi` lines 21–67.
- Affected Rust operation:
  `StructureConstantsAlgebra::primitive_idempotents`.
- Difference: GAP returns the identical cached mutable list object. Rust caches
  the primitive coordinate rows on the shared algebra allocation and
  reconstructs an owned `Vec<StructureConstantsAlgebraElement<_>>` for each
  call; `cached_primitive_idempotents` exposes the bound cache state directly.
- Reason: callers should not be able to mutate an algebra attribute by editing
  a returned collection, while clones of the algebra must share the computed
  attribute.
- Compatibility impact: repeated calls return equal elements with the
  identical parent algebra and do not rerun the probabilistic computation, but
  the outer Rust `Vec` has a fresh allocation.
- Evidence: the pinned idempotent fixture checks GAP's identical cache, while
  Rust tests check unbound/bound transitions, stable repeated values, and
  shared algebra-parent identity.

## MODULEHOM-001 — Semisimple block counting uses the split-prime tier

- Source: QPA `lib/modulehom.gi`, `NumberOfNonIsoDirSummands` at lines
  1761–1779; GAP `CentralIdempotentsOfAlgebra` and
  `DirectSumDecomposition` dependencies.
- Affected Rust operation: `number_of_nonisomorphic_direct_summands`.
- Difference: QPA can factor the centre of `End(M)/rad End(M)` over the
  rationals and finite extension fields. Rust currently follows the same
  radical-quotient, central-idempotent, and block-decomposition route on
  Numerica prime fields; rational input reaches the exact radical quotient
  and then returns the typed `FiniteFieldRequired` central-splitting error.
- Reason: the pinned Rust coefficient tier has GAP-compatible finite-field
  polynomial factorization for prime fields, while Numerica 2.1 provides no
  rational polynomial factorizer or finite extension-field domain. Replacing
  QPA's block method with a decomposition-derived formula would violate the
  source-faithful porting rule.
- Compatibility impact: supported prime-field results retain QPA's block
  order and dimensions. Rational and extension-field calls are explicitly
  outside this operation's current coefficient tier instead of using a
  different computational method.
- Evidence: `tests/fixtures/moduledecomp_modules.schema-v1.json` pins the
  repeated and mixed results over `GF(2)` and `GF(5)`; the path-module test
  pins the rational typed boundary.

## MODULEHOM-002 — Endomorphism-quiver conversion handles the zero radical

- Source: QPA `lib/modulehom.gi`, `EndOfModuleAsQuiverAlgebra` at lines
  3634–3650, delegating to `AlgebraAsQuiverAlgebra` in `lib/algebra.gi` at
  lines 634–926.
- Affected Rust operation: `end_of_module_as_quiver_algebra` through
  `algebra_as_quiver_algebra`.
- Difference: on a one-dimensional simple module, QPA aborts while forming
  `e_i Basis(rad(A)/rad(A)^2) e_j` because GAP's list inner product has no
  summands when the radical basis is empty. Rust treats that basis as the zero
  space, obtains the zero-arrow quiver, and returns the field algebra.
- Reason: the empty list represents the mathematically well-defined zero
  radical. Handling it explicitly is the zero-dimensional branch of QPA's
  stated `rad(A)/rad(A)^2` construction and avoids an incidental GAP container
  error without selecting a different computational method.
- Compatibility impact: nonzero-radical successes and QPA's structural
  rejections are unchanged. Only this accidental source crash becomes a
  successful zero-arrow conversion.
- Evidence: `tests/gap-reference/modulehom_end_quiver_simple_error.g` pins the
  QPA error; the Rust algebra-conversion test pins the repaired field-algebra
  result. `tests/fixtures/modulehom_end_quiver.schema-v1.json` pins successful
  nonzero-radical conversions over `GF(2)` and `GF(5)`.

## MODULEPROJRES-001 — Stalled FirstPart reduction is reported

- Source: QPA `lib/moduleprojres.gi`, `FirstPart` at lines 648–682.
- Affected Rust operation: `first_part`.
- Difference: QPA repeats its outer reduction loop until the entered vector is
  zero, but it has no branch for a complete pass in which no entered `f` or
  `fprime` tip divides the current tip. Such input therefore loops forever.
  Rust returns `ProjectivePathModuleError::FirstPartNotSpanned` after that
  unchanged pass.
- Reason: the valid GSZ-resolution domain assumes the vector is spanned by the
  entered reducers. Reporting violation of that precondition prevents an
  uninterruptible library call without changing the literal reducer order,
  divisions, coefficient accumulation, or valid result.
- Compatibility impact: valid spanned inputs agree exactly. An invalid input
  that never terminates in QPA returns a typed Rust error.
- Evidence: `tests/fixtures/moduleprojres_xset_firstpart.schema-v1.json` pins
  ordinary and prime-generator reductions over `GF(2)`, `GF(5)`, and the
  rationals; the projective-path boundary test pins the stalled Rust result.
