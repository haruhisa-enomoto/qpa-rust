# Conventions

Several QPA conventions differ from common Python, Rust, or linear-algebra defaults.
They are part of the compatibility contract.

## Indices and names

- Mathematical vertex and arrow indices are one-based, matching GAP/QPA.
  Python methods such as `vertex_at(1)`, `arrow_at(1)`, and `factor(1)` retain
  that convention. Rust `VertexIndex` and `ArrowIndex` constructors reject
  zero.
- Python `list` positions and `obj[index]` protocols are ordinary zero-based
  positions. Rust slices, vectors, matrix row/column offsets, and iterator
  positions are likewise zero-based. Named methods whose argument is described
  as a QPA position remain one-based.
- Quiver generator order is all vertices followed by all arrows. Constructors
  preserve entered arrow order and QPA's default naming rules.
- Objects from equal-looking but separately constructed quivers may belong to
  different QPA families. Do not reconstruct a quiver merely to look up a
  path; retain or clone its handle.

## Multiplication and composition

- Path multiplication is left-to-right: Python `a * b` and Rust
  `a.multiply(&b)` mean traverse `a` first and then `b`.
- Path-algebra multiplication inherits that path order.
- `PathAlgebraModuleHomomorphism::compose(&right)` applies `self` first and
  `right` second, matching QPA's row-action convention.

## Right modules and matrices

Module elements are row vectors. For an arrow `a: i -> j`, its action matrix
has shape `dim(M_i) x dim(M_j)`, and a row at vertex `i` is multiplied on the
right by that matrix.

In Python, matrices are nested row lists. `module_from_integer_matrices`
accepts integer matrices and converts their entries exactly into the parent
field. For the more general constructors, use `field.element(...)` values.

Rust `DenseMatrix::from_rows` infers its shape, so an empty row list is `0 x 0`.
Use `DenseMatrix::from_flat(0, n, vec![], field)` or its `m x 0` counterpart
when a particular empty shape matters.

QPA stores physical one-row or one-column zero matrices at zero-dimensional
module vertices because GAP matrices do not retain arbitrary empty shapes.
The Rust and Python module layers preserve those stand-ins where they are observable;
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
