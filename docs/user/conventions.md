# Conventions

Several QPA conventions differ from common Rust or linear-algebra defaults.
They are part of the compatibility contract.

## Indices and names

- `VertexIndex` and `ArrowIndex` are one-based, matching GAP/QPA. Their checked
  constructors reject zero.
- Rust slices, vectors, matrix row/column offsets, and iterator positions are
  ordinary zero-based Rust positions.
- Quiver generator order is all vertices followed by all arrows. Constructors
  preserve entered arrow order and QPA's default naming rules.
- Objects from equal-looking but separately constructed quivers may belong to
  different QPA families. Do not reconstruct a quiver merely to look up a
  path; retain or clone its handle.

## Multiplication and composition

- Path multiplication is left-to-right: `a.multiply(&b)` means traverse `a`
  first and then `b`.
- Path-algebra multiplication inherits that path order.
- `PathAlgebraModuleHomomorphism::compose(&right)` applies `self` first and
  `right` second, matching QPA's row-action convention.

## Right modules and matrices

Module elements are row vectors. For an arrow `a: i -> j`, its action matrix
has shape `dim(M_i) x dim(M_j)`, and a row at vertex `i` is multiplied on the
right by that matrix.

`DenseMatrix::from_rows` infers its shape, so an empty row list is `0 x 0`.
Use `DenseMatrix::from_flat(0, n, vec![], field)` or its `m x 0` counterpart
when a particular empty shape matters.

QPA stores physical one-row or one-column zero matrices at zero-dimensional
module vertices because GAP matrices do not retain arbitrary empty shapes.
The Rust module layer preserves those stand-ins where they are observable;
logical dimension vectors and basis coefficient lists still use dimension
zero. For the checked three-argument module constructor, omit an arrow matrix
incident with a zero-dimensional space so the source-compatible stand-in can
be created.

## Exact order and output

Basis, nullspace, nontip, Hom, and decomposition order often affects later QPA
computations. Returned lists therefore preserve source ordering even when a
different order would span the same space. Display helpers reproduce QPA
presentation where it is part of tested behavior; they are not parsers or a
substitute for structured accessors.
