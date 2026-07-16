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
