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
