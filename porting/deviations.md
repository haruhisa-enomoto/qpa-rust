# Intentional compatibility deviations

This ledger records places where Rust's public boundary intentionally differs
from the pinned GAP/QPA surface. It does not authorize algorithmic deviations:
the computations remain translations of the cited source methods.

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
