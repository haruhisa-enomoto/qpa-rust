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

For a quotient algebra, the checked constructor verifies the defining
relations and raises `qpa.RelationNotSatisfiedError` if the matrices do not
define a module. Other malformed inputs raise `qpa.InvalidInputError`.

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

The high-level Hom-space wrapper exposes its dimension and generators:

```python
hom = algebra.hom(simple_at_v1, module)
assert hom.dimension == 0

basis = simple_at_v2.hom_basis_to(module)
assert basis is not None
```

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

The runnable version of the main example is
[`docs/examples/python/modules.py`](../examples/python/modules.py).
