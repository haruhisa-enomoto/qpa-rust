# Intentional compatibility deviations

This ledger records places where Rust's public boundary intentionally differs
from the pinned GAP/QPA surface. It does not authorize algorithmic deviations:
the computations remain translations of the cited source methods.

## COMBMAP-001 — Invalid winding pairs terminate with a typed error

- Source: QPA `lib/combinatorialmap.gi`, `WindingNumber` at lines 349–395.
- Affected Rust operation: `CombinatorialMap::winding_number`.
- Difference: QPA repeatedly applies the ordering permutation until the first
  half-edge of each pair reaches the second. If the two points lie in different
  ordering orbits, that loop never terminates. Rust detects the return to the
  starting point and returns
  `CombinatorialMapError::CurvePairInDifferentVertices`.
- Reason: the QPA documentation requires a closed curve in minimal position,
  so different-vertex pairs violate the method's input contract. A typed error
  makes this failure finite without changing any documented computation.
- Compatibility impact: valid curves have the exact QPA winding number and
  rotation behavior. Invalid different-vertex pairs return instead of hanging.
- Evidence: `tests/fixtures/combinatorialmap_foundation.schema-v1.json` pins
  the valid source results, and the Rust boundary test covers the finite error.

## COMBMAP-002 — Topological searches report violated preconditions

- Source: QPA `lib/combinatorialmap.gi`, `NonSeperatingCurve` at lines
  585–672 and `HomologyBasisOfCombinatorialMap` at lines 880–930.
- Affected Rust operations: `CombinatorialMap::non_separating_curve` and
  `CombinatorialMap::homology_basis`, including their cut/join helpers.
- Difference: QPA assumes a positive-genus map with canonical edge pairs and
  well-formed non-separating and boundary-joining curves. Outside that domain,
  its unbounded searches can run forever or later access an unbound local or
  list position. Rust validates the finite curve, boundary, and index inputs,
  detects a repeated search state or nondecreasing genus, and returns a typed
  `CombinatorialMapError`.
- Reason: these are the source algorithms' documented preconditions. Making a
  violated precondition finite prevents an uninterruptible library call while
  preserving the precise searches and permutation mutations on valid maps.
- Compatibility impact: positive-genus canonical maps and valid cut/join
  curves have QPA's exact permutations, boundary order, index remapping,
  homology basis order, winding numbers, and attribute caching. Malformed
  inputs return instead of hanging or raising a late GAP interpreter error.
- Evidence: `tests/fixtures/combinatorialmap_homology.schema-v1.json` pins the
  complete genus-one and genus-two source computations; the Rust topology
  boundary test covers the finite typed errors.

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

## PATHALG-002 — Path rings use the field coefficient tier

- Source: QPA `lib/pathalg.gd` and `lib/pathalg.gi`, `IsPathRing`, `PathRing`,
  and `PathAlgebra` at lines 4–14 and 5–88.
- Affected Rust surface: `PathAlgebra`, `PathAlgebraElement`,
  `PathAlgebra::new`, and `PathAlgebra::is_path_ring`.
- Difference: QPA's `PathRing(R, Q)` accepts a general coefficient ring and is
  only an algebra when GAP's coefficient-domain dispatch makes it one. Rust's
  initial Numerica contract requires `F: Field`; within that tier `PathRing`
  and `PathAlgebra` have the same typed representation. Nonfield examples such
  as `PathRing(Integers, Q)` are not admitted yet.
- Reason: the agreed first coefficient tier is the rational field and
  Numerica-supported prime fields. This keeps the foundational sparse and
  Gröbner algorithms exact without building unrelated general-ring machinery
  before the field-backed QPA port is complete.
- Compatibility impact: every supported Rust path algebra is both a QPA path
  ring and path algebra, with the identical quiver, identity, ordering, sparse
  normalization, and arithmetic. General coefficient rings require a later
  coefficient-tier extension.
- Evidence: `tests/fixtures/pathalg_surface.schema-v1.json` pins the field path
  ring/algebra categories and independent constructor identity; the existing
  `path_algebra` fixture pins `QQ`, `GF(2)`, and `GF(5)` arithmetic.

## PATHALG-003 — Broken list uniformity dispatch is a typed error

- Source: QPA `lib/pathalg.gi`, the list-plus-path `IsLeftUniform` and
  `IsRightUniform` methods at lines 574–663.
- Affected Rust operations: `PathAlgebra::is_left_uniform_list_at_path` and
  `PathAlgebra::is_right_uniform_list_at_path`.
- Difference: both installed GAP methods accept an `IsPath` second argument
  and then call `CoefficientsAndMagmaElements` on that path. No such method is
  installed, so every reachable call raises GAP's generic no-method error
  before inspecting the list. Rust returns
  `PathAlgebraError::UniformityListPathDispatch`.
- Reason: a typed library cannot reproduce an interpreter break loop as a
  value. The error preserves the exact failing branch without substituting an
  independently repaired endpoint test.
- Compatibility impact: the useful one-element uniformity methods agree
  exactly. Both defective list overloads fail deterministically in Rust just
  as they do in QPA, but callers receive a matchable error value.
- Evidence: `tests/gap-reference/pathalg_surface.g` pins both QPA no-method
  failures, and `tests/fixtures/pathalg_surface.schema-v1.json` plus the Rust
  fixture test pin their false success states.

## PATHALG-004 — Broken quotient reordering is a typed error

- Source: QPA `lib/pathalg.gi`, quotient `OrderedBy` at lines 791–802 and
  `FactorPathAlgebraByRelators` at lines 215–309.
- Affected Rust operation: `QuotientPathAlgebra::ordered_by`.
- Difference: the QPA factor constructor stores `pathAlgebra`, `ideal`, and
  `familyRing` components on the element family, but its `OrderedBy` method
  reads `fam!.relators`, which the constructor never assigns. Every ordinary
  quotient call therefore raises a GAP record-element error. Rust returns
  `QuotientError::UnboundRelatorsForOrderedBy` at the same branch.
- Reason: returning a typed error preserves the pinned source behavior without
  silently repairing the method by reading ideal generators instead.
- Compatibility impact: no normally constructed QPA quotient can reach the
  advertised reordering calculation. Rust likewise does not reorder it; the
  failure is deterministic and matchable instead of entering a GAP break loop.
- Evidence: `tests/gap-reference/pathalg_mapping_ordering.g` pins the unbound
  `relators` failure and
  `tests/fixtures/pathalg_mapping_ordering.schema-v1.json` pins its false
  success state alongside valid mapped-expression behavior.

## PATHALG-005 — Generator-variable assignment returns a namespace plan

- Source: QPA `lib/pathalg.gi`, `AssignGeneratorVariables` at lines
  2600–2633, and pinned GAP `IsValidIdentifier` in `lib/global.gi` lines 47–52.
- Affected Rust operations: `QuiverAlgebra::generator_variable_assignments`,
  `QuiverAlgebra::generator_variable_assignments_with_read_only`, and
  `GeneratorVariableAssignment`.
- Difference: QPA validates generator strings, warns for already-bound GAP
  globals, unbinds them, and mutates the GAP global namespace. Rust returns the
  ordered name/`One(A) * g` assignment plan. A binding supplies its namespace's
  read-only predicate and decides how to warn about or overwrite bound names.
- Reason: a reusable Rust algebra library has no GAP global-variable table and
  must not mutate a caller language's process-wide namespace implicitly. The
  future Python layer can install the returned values in an explicit mapping.
- Compatibility impact: GAP's exact identifier alphabet, keyword exclusion,
  all-digit rejection, pre-assignment validation order, read-only branch,
  vertex-then-arrow order, and repeated-name last-wins order are preserved.
  Warning emission and actual namespace mutation occur only when the caller
  applies the returned plan.
- Evidence: `tests/gap-reference/pathalg_mapping_ordering.g` pins successful
  assignment plus invalid and read-only failures;
  `tests/fixtures/pathalg_mapping_ordering.schema-v1.json` pins assigned names
  and algebra elements in source order.

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

## MODULE-001 — Scalar-zero module display failure is typed

- Source: QPA `lib/module.gi`, `Display` at lines 967–984, delegating the
  stored arrow entries to GAP's `PrintArray`.
- Affected Rust operation: `PathAlgebraModule::qpa_display_string_with`.
- Difference: for the zero module over a quiver with an arrow, QPA stores the
  integer zero as that arrow entry and aborts when `PrintArray` requires a
  dense list. Rust returns
  `PathAlgebraModuleError::DisplayScalarZeroArrowMatrix` at the same branch.
- Reason: Rust represents the QPA scalar-zero sentinel explicitly and cannot
  pass it to a matrix-only formatter. The typed result preserves the source
  branch without silently rendering a different zero matrix.
- Compatibility impact: all dense arrow matrices retain QPA's exact labels,
  alignment, spacing, and line breaks. Only QPA's interpreter error becomes a
  typed Rust error.
- Evidence: `tests/fixtures/module_presentation.txt` pins successful path-,
  quotient-, rational-, and prime-field presentation text;
  `tests/gap-reference/module_presentation_zero_display_error.g` pins the
  source error at `module.gi:981`.

## MODULE-002 — Zero-module canonical-basis failure is typed

- Source: QPA `lib/module.gi`, `BasisNC` at lines 244–261, the default
  `Basis` method at lines 298–306, and `CanonicalBasis` at lines 320–369.
- Affected Rust operation:
  `PathAlgebraModuleElementBasis::canonical_for_module`.
- Difference: `CanonicalBasis` contains an explicit empty-generator branch,
  but for `ZeroModule(A)` its preceding `Basis(V)` call reaches `BasisNC`,
  whose empty shallow-copy argument does not satisfy the installed
  `NewBasis` method and aborts. Rust returns
  `PathAlgebraModuleError::ZeroModuleCanonicalBasisFailure` before the same
  unreachable branch.
- Reason: Rust has no GAP method dispatcher whose empty-list family can fail
  dynamically. A typed error preserves the actual pinned QPA outcome without
  repairing the source method.
- Compatibility impact: every nonzero module retains QPA's vertex-major
  canonicalization and stored canonical property. Only the generic GAP
  `NewBasis` dispatch error becomes a typed Rust error.
- Evidence: `tests/fixtures/module_basis.txt` pins nonzero `QQ` and `GF(2)`
  bases and coefficients; `tests/gap-reference/module_basis_zero_error.g`
  pins the source failure at `module.gi:254`.

## MODULE-003 — Unchecked isolated-vertex failure is typed

- Source: QPA `lib/module.gi`, the two-argument
  `RightModuleOverPathAlgebraNC` method at lines 596–754.
- Affected Rust operation:
  `right_module_over_path_algebra_infer_dimension_nc`.
- Difference: the checked source method replaces every isolated vertex's
  dimension sentinel `-1` by zero, but the NC copy omits that loop. Its
  `CreateModuleBasis` path consequently aborts with "Dimensions of element are
  incorrect" on an ordinary nonzero example. Rust detects isolated vertices
  and returns `PathAlgebraModuleError::UncheckedInferredIsolatedVertices`.
- Reason: Rust dimension vectors use `usize` and cannot construct the invalid
  negative-dimensional element family that causes QPA's later interpreter
  failure. The error is reported at the point where the sentinel remains.
- Compatibility impact: NC inputs without isolated vertices retain QPA's
  inferred dimensions, matrix order, relation test, and deliberate omission
  of coefficient-field checking. Inputs with isolated vertices return the
  typed error; callers can use the checked overload, which gives those
  vertices QPA's intended dimension zero.
- Evidence: `tests/fixtures/module_inferred_constructors.txt` pins all three
  checked input forms and the nonisolated NC result over `QQ` and `GF(2)`;
  `tests/fixtures/module_inferred_constructor_errors.schema-v1.json` pins the
  source NC isolated-vertex failure and unchecked foreign-field acceptance.

## MODULE-004 — Matrix-submodule source failures are typed

- Source: QPA `lib/module.gi`, matrix-module `SubAlgebraModule` at lines
  1016–1060, `SubmoduleAsModule` at lines 1070–1128, and
  `NaturalHomomorphismBySubAlgebraModule` at lines 1137–1252.
- Affected Rust operations: `subalgebra_module`,
  `PathAlgebraModuleSubmodule::canonical_basis_elements`,
  `submodule_as_module`, and
  `natural_homomorphism_by_subalgebra_module`.
- Difference: QPA's nonempty quotient-module submodule acts by vertices from
  the original path-algebra family and fails method dispatch; dependent
  retained uniform generators reach a nonsquare matrix inverse; zero
  submodule conversion indexes a missing basis vector; zero endpoint
  dimensions are passed as integer pairs to a positional constructor and
  misdispatch; and the natural map's zero-submodule subspace check cannot
  obtain images. Rust returns a distinct `PathAlgebraModuleError` at each
  corresponding branch.
- Reason: Rust cannot recreate GAP method-dispatch, unbound-list, or nonsquare
  inversion interpreter failures as values. The typed errors preserve the
  pinned public outcome without repairing the source branch. GAP's generic
  higher-ranking empty-list dispatch and empty canonical basis remain
  successful in Rust.
- Compatibility impact: valid path-algebra submodules retain QPA's
  generator-major vertex splitting, canonical-basis closure order, standalone
  matrices, quotient dimensions, and basis images. Invalid cases return typed
  errors instead of the recorded GAP messages.
- Evidence: `tests/fixtures/module_submodules.txt` pins successful `QQ` and
  `GF(2)` results, including the natural map's leaked earlier-vertex
  coordinates; `tests/fixtures/module_submodule_errors.schema-v1.json` pins
  every listed dispatch and interpreter failure.

## MODULE-005 — Generic Hom/End dispatch failures are typed

- Source: QPA `lib/module.gi`, `HOPF_MatrixModuleHom` at lines 1255–1360,
  generic `Hom` at lines 1369–1411, and generic `End` at lines 1419–1441.
- Affected Rust operations: `generic_hom_over_path_ring` and
  `generic_end_over_path_ring`.
- Difference: QPA's helper indexes an unavailable ordinary basis for a zero
  module, and `NullspaceMat(TransposedMat(equations))` loses its intended
  width when the equation matrix has no rows. Its installed methods also do
  not dispatch on quotient path algebras, while a structurally equal but
  nonidentical acting algebra is rejected explicitly. Rust returns a distinct
  `PathAlgebraModuleError` for each branch.
- Reason: GAP method-selection, unbound-basis, and empty-matrix nullspace
  failures cannot be represented as Rust return values. The errors preserve
  the source outcomes without extending QPA's domain or substituting the
  separately ported `HomOverAlgebra` algorithm.
- Compatibility impact: valid modules and submodules retain QPA's exact HOPF
  equation layout, nullspace order, full matrices, canonical-basis rebasing,
  installed algebra-homomorphism flag, and End basis. Invalid cases return
  typed errors instead of the pinned GAP diagnostics.
- Evidence: `tests/fixtures/module_generic_hom.txt` pins all four ambient and
  submodule source/range combinations over `QQ` and `GF(2)`;
  `tests/fixtures/module_generic_hom_errors.schema-v1.json` pins the zero,
  arrowless, wrong-algebra, quotient-dispatch, and returned-space basis
  failures.

## MODULE-006 — General right modules use a based action contract

- Source: QPA `lib/module.gd` line 41 and
  `RightAlgebraModuleToPathAlgebraMatModule` in `lib/module.gi` at lines
  1973–2041.
- Affected Rust operations: `BasedRightAlgebraModule::new`,
  `BasedRightAlgebraModule::from_path_algebra_module`, and
  `right_algebra_module_to_path_algebra_mat_module`.
- Difference: GAP accepts its open-ended
  `IsRightAlgebraModuleElementCollection` hierarchy and then rejects
  nonfinite modules or modules over unrelated algebras at runtime. Rust uses a
  finite based-module boundary whose common square matrices give right action
  by each vertex and arrow generator; count, field, and shape violations are
  typed at construction. An invalid claimed action that leaves a target
  vertex space or violates quotient relations is also typed.
- Reason: Rust has no ambient GAP object/filter system from which to query an
  arbitrary element action. Generator matrices are the finite data read by
  QPA's method and are sufficient to execute the source algorithm without
  introducing an independent representation method.
- Compatibility impact: every representable finite path-algebra or quotient
  module retains QPA's entered-basis projection order,
  insertion-preserving `Basis(Subspace(...))` order, arrow coordinate solves,
  omitted zero-endpoint maps, and final checked constructor. Nonfinite and
  unrelated-algebra objects are excluded by the Rust input type rather than
  failing after GAP dispatch.
- Evidence: `tests/fixtures/module_right_algebra_conversion.txt` pins regular,
  changed-basis, zero, already-matrix, quotient, and arrowless inputs over
  `QQ` and `GF(2)`; the changed basis distinguishes QPA's semi-echelon
  insertion order from a canonical RREF replacement.
  `tests/fixtures/module_right_algebra_conversion_errors.schema-v1.json` pins
  QPA's nonfinite and unrelated-algebra failures, while Rust unit tests pin
  the based-action boundary failures.

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

## MODULEHOMALG-001 — Finite-field zero Iyama failure is typed

- Source: QPA `lib/modulehomalg.gi`, `IyamaGenerator` at lines 692–715.
- Affected Rust operation: `iyama_generator`.
- Difference: for a zero module over a finite field, QPA obtains `fail` from
  `DecomposeModuleWithMultiplicities` and then indexes `L[1]`, aborting with a
  generic GAP list error. Rust returns
  `IyamaGeneratorError::MissingMultiplicityDecomposition` at the same branch.
- Reason: Rust has no untyped `fail` value to index. The typed result identifies
  the precise failed QPA intermediate without repairing or bypassing it.
- Compatibility impact: every nonzero rational and finite-prime-field radical
  filtration and finite-field basicification retains QPA's summand and matrix
  order. Rational zero input still succeeds; finite-field zero input returns a
  typed error instead of the source interpreter error.
- Evidence: `tests/fixtures/modulehomalg_iyama.schema-v1.json` pins partial and
  complete generators for three nontrivial module shapes over `QQ` and
  `GF(2)`, while `tests/gap-reference/modulehomalg_iyama_zero_error.g` pins the
  source failure at `modulehomalg.gi:710`.

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

## DEGORDER-001 — Numerical AR positions are typed and zero-based

- Source: QPA `lib/degorderfinitetype.gd` and `lib/degorderfinitetype.gi`
  lines 61–1104.
- Affected Rust operations: `ArQuiverNumerical`, `ArMesh`, `ArOrbit`,
  `ArOrbitPosition`, `NumericalModule`, and `DimensionVectorInput`.
- Difference: GAP encodes meshes, translate orbits, indecomposable selectors,
  and dimension/multiplicity vectors in nested untyped lists with one-based
  positions. Rust gives each form its own type and uses zero-based positions;
  results that are vectors retain QPA's original coordinate order.
- Reason: typed inputs prevent an orbit position, dimension vector, and module
  multiplicity vector from being silently interchanged, while zero-based
  indices are the native Rust boundary convention.
- Compatibility impact: converting every GAP position by subtracting one gives
  the identical AR matrix, inverse-transpose Hom matrix, simple positions,
  enumeration order, and degeneration comparisons. Numerica's exact rational
  matrix inverse and arbitrary-precision integers preserve QPA's exact integer
  calculation without a field-dependent replacement.
- Evidence: `tests/fixtures/degorderfinitetype.schema-v1.json` pins explicit
  mesh and compact orbit constructors, all numerical operations, and nontrivial
  direct-order filtering against QPA.

## DEGORDER-002 — Printing operations return deterministic Rust data

- Source: QPA `lib/degorderfinitetype.gi` lines 997–1275.
- Affected Rust operations: `degeneration_order_leq`,
  `format_multiplicity_vector`, `format_multiplicity_vectors`, `Display for
  ArQuiverNumerical`, and `test_codimension_conjecture`.
- Difference: QPA writes a dimension-mismatch warning, multiplicity displays,
  conjecture progress, and runtime summaries to the current GAP output stream.
  Rust returns the checked comparison, formatted strings, and a structured
  `CodimensionConjectureReport`; the accepted `print_progress` flag does not
  write to global process output, and elapsed runtime is not part of the result.
- Reason: a library operation should not mutate global output state, and wall
  time is nondeterministic. The mathematical branch decisions, enumeration,
  early-break flag, counterexamples, and counters remain QPA's.
- Compatibility impact: callers print the returned formatting/report when
  interactive output is desired. The exact stable Print/View and multiplicity
  text is unchanged; only warnings, progress, and timing require caller policy.
- Evidence: `tests/fixtures/degorderfinitetype_print.txt` pins QPA's exact
  stable text, while `tests/fixtures/degorderfinitetype.schema-v1.json` pins the
  comparisons and conjecture result/counters.

## DEGORDER-003 — Predefined AR dispatch is a typed selector

- Source: QPA `lib/degorderfinitetype.gi`, `PredefARQuivers` at lines
  1303–1794.
- Affected Rust operations: `PredefinedArQuiver`,
  `predefined_ar_quiver_data`, `predefined_ar_quiver`, and
  `PREDEFINED_AR_QUIVER_CATALOG`.
- Difference: QPA dispatches on a variable-length untyped argument list and
  prints the `"what"` catalog before returning `false`. Rust uses an enum whose
  optional parameter represents the two- and three-argument BG forms, returns
  a typed error for QPA's wrong-argument branch, and exposes the stable catalog
  as a string constant.
- Reason: the enum makes QPA's finite dispatch domain explicit and prevents
  strings, list numbers, and parameters from being silently interchanged. A
  constant lets callers choose where to display the catalog.
- Compatibility impact: every supported table, generated series, and orbit
  encoding is identical after the documented one-based to zero-based position
  conversion. Missing, superfluous, zero, and unsupported parameters return a
  typed error instead of a GAP interpreter error.
- Evidence:
  `tests/fixtures/degorderfinitetype_predefined.schema-v1.json` pins all 36
  distinct fixed and representative parameterized dispatch cases, including
  both BG4 parity branches, and
  `tests/fixtures/degorderfinitetype_predefined_catalog.txt` pins the exact QPA
  catalog output.
