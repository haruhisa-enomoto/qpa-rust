# Right modules and homomorphisms in Python

QPA represents finite-dimensional right modules by a dimension vector and one
right-action matrix for each relevant arrow. The Python binding keeps that
model explicit.

## Construct a right matrix module

```python
import qpa

field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
algebra = qpa.PathAlgebra(field, quiver)

module = algebra.module_from_integer_matrices(
    [1, 1],
    [("a", [[1]])],
)

assert module.dimension_vector == [1, 1]
assert module.dimension == 2
assert module.arrow_matrix_shapes == [(1, 1)]
```

For an arrow `a: i -> j`, the matrix has shape
`dim(M_i) x dim(M_j)`. Module elements are row vectors and act on the right.
Integer entries are converted exactly into the algebra's field.

Use `module_from_matrices` on a path, quotient, or common `QuiverAlgebra`
handle for entries containing exact Python integers, `fractions.Fraction`,
or `qpa.Scalar`. Omit arrow actions that are zero, including arrows incident
to a zero-dimensional vertex; do not supply `[]` or `[[]]` for these actions.

For a quotient algebra, the checked constructor verifies the defining
relations and raises `qpa.RelationNotSatisfiedError` if the matrices do not
define a module. Invalid matrix shapes raise `qpa.InvalidInputError`;
unsupported coefficient types raise `TypeError`.

## Standard modules and structure

```python
simple_at_v1, simple_at_v2 = algebra.simple_modules()
projective_at_v1, projective_at_v2 = algebra.indecomposable_projectives()

assert simple_at_v1.dimension_vector == [1, 0]
assert projective_at_v1.dimension_vector == [1, 1]
assert module.radical_series() == [[1, 0], [0, 1]]
```

The order of returned simples and indecomposable projectives follows QPA's
vertex order. Do not sort these lists unless the subsequent mathematics is
independent of that order.

## Hom spaces and maps

For matrix modules, the structured Hom result exposes its dimension, basis
maps, and coordinate function:

```python
hom = simple_at_v1.hom_basis_with_coefficients_to(module)
assert hom is not None
assert hom.dimension == 0

basis = simple_at_v2.hom_basis_to(module)
assert basis is not None
```

`hom_basis_to` returns a list of `ModuleHomomorphism` objects, while
`hom_basis_with_coefficients_to` returns an object with `.dimension`, `.basis`,
and `.coefficients(map)`. The separate `PathAlgebra.hom(source, target)` route
exposes QPA's generic Hom space: its `.generators` are
`GenericModuleHomomorphism` objects, which do not provide the kernel and
composition methods of `ModuleHomomorphism`.

When only the dimension is needed, use `source.hom_dimension_to(target)`.
It returns the same dimension as `len(source.hom_basis_to(target))` without
constructing the basis maps, or `None` for different acting algebras. Both
methods reject incompatible coefficient fields.

For degree-one Ext, `source.ext1_dimension(target)` computes the dimension
without quotient representatives. To reuse the source's projective cover and
syzygy across several targets, use `source.ext1_dimensions(targets)`; the
returned dimensions follow the input order, and an empty list returns `[]`.
All targets must belong to the source's acting algebra. Use `source.ext1(target)`
when representatives or their coefficient function are needed.

For an explicit map, use `source.homomorphism_to(range, vertex_matrices)`.
It accepts the same exact entry types as `module_from_matrices`. Supply one
matrix per vertex, with source dimension as the number of rows and target
dimension as the number of columns. Internally, a zero-dimensional vertex
uses a physical dimension of one: any incident map matrix must be zero and
have shape `max(1, source_dimension) x max(1, target_dimension)`.
Its `kernel`, `image_module`, and `cokernel` methods preserve QPA's
constructions, and the corresponding inclusion/projection methods retain the
canonical maps:

```python
identity = module.identity()
assert identity.is_isomorphism
assert identity.kernel().is_zero
assert identity.cokernel().is_zero
```

Composition follows QPA's row-action order: `left.compose(right)` applies
`left` first and `right` second.

| Operation | Result |
| --- | --- |
| `f.kernel()` | Kernel module |
| `f.kernel_inclusion()` | Map from the kernel into `f.source` |
| `f.image_module()` | Image module |
| `f.image(element)` | Value of the map at an element |
| `f.cokernel()` | Cokernel module |
| `f.cokernel_projection()` | Map from `f.range` onto the cokernel |
| `MatrixModule.direct_sum(modules)` | Sum module, or `None` for an empty list |
| `sum_module.direct_sum_inclusions()` | Ordered inclusions of the summands |
| `sum_module.direct_sum_projections()` | Ordered projections onto the summands |

Dimensions, `dimension_vector`, and predicates such as `is_zero` are
properties: use `module.dimension`, not `module.dimension()`. Act on a module
element with `element.act_by(algebra_element)`; it accepts concrete and common
algebra elements while checking their parent algebra and field.

The runnable version of the main example is
[`docs/examples/python/modules.py`](../examples/python/modules.py).
