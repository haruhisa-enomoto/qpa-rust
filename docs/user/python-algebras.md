# Quivers, path algebras, and quotients in Python

This page follows the same construction chain as QPA: define a quiver, form
its path algebra over an exact field, express relations as algebra elements,
and pass those relations to QPA's quotient route.

## A bound quiver algebra

```python
import qpa

quiver = qpa.Quiver(
    ["v1", "v2", "v3"],
    [("v1", "v2", "a"), ("v2", "v3", "b")],
)
path_algebra = qpa.PathAlgebra(qpa.Field.rational(), quiver)

a = path_algebra.generator("a")
b = path_algebra.generator("b")
algebra = path_algebra.quotient([a * b])

assert algebra.dimension() == "5"
assert (algebra.generator("a") * algebra.generator("b")).is_zero
```

The complete checked example is
[`docs/examples/python/algebras.py`](../examples/python/algebras.py).

`PathAlgebra.quotient` follows QPA's configured default
`GBNPGroebnerBasis` call chain. `PathAlgebra.quotient_high_level` exposes the
other translated high-level QPA route explicitly; it is not used implicitly
as a replacement checker.

## Constructing quivers

The direct constructor uses explicit names:

```python
qpa.Quiver(
    ["source", "sink"],
    [("source", "sink", "arrow")],
)
```

QPA-style numbered vertices use one-based indices:

```python
indexed_quiver = qpa.Quiver.indexed(3, [(1, 2, "a"), (2, 3, "b")])
assert indexed_quiver.vertex_at(1).name == "v1"
assert indexed_quiver.arrow_at(1).name == "a"
```

Constructors also exist for adjacency matrices and oriented Dynkin quivers;
their exact signatures are listed in the [Python API reference](python-api.md).

## Paths and algebra elements

`Quiver.path` accepts arrow names in traversal order. `PathAlgebra.generator`
returns an algebra element for either a vertex or an arrow name.

```python
path = quiver.path(["a", "b"])
element = path_algebra.element(path)

assert path.source == "v1"
assert path.target == "v3"
assert element == a * b
```

Elements from separately constructed, equal-looking algebras cannot be mixed.
Keep the original algebra and quiver handles; see
[identity and caches](identity-and-randomness.md).

## Working with a common algebra handle

`PathAlgebra` and `QuotientAlgebra` provide `as_quiver_algebra()` to obtain a
`QuiverAlgebra` handle. The handle retains the original algebra's identity,
relations, and caches. It provides the same standard module constructors:

```python
common = algebra.as_quiver_algebra()
simples = common.simple_modules()
projectives = common.indecomposable_projectives()
injectives = common.indecomposable_injectives()

assert [module.dimension_vector for module in projectives] == [
    [1, 1, 0], [0, 1, 1], [0, 0, 1]
]
assert simples[0].hom_dimension_to(simples[0]) == 1
assert simples[0].ext1_dimension(simples[1]) == 1
assert common.zero_module().dimension == 0
```

These lists use vertex order. Module `dimension` and `dimension_vector` are
properties. For an explicit Hom basis use `hom_basis_to`; for Ext representatives
use `ext1`, whose result has a `dimension` property.

The common handle also provides `generator(name)`, `generators`, `zero()`, and
`one()`. These return `QuiverAlgebraElement` values tied to the same algebra.
The `generators` property lists vertices followed by arrows, matching the
concrete Python algebra. The identity is available separately through `one()`.

## Exact module matrices

All three algebra handles support `module_from_matrices(dimension_vector,
actions)`. Each action pairs an arrow name with its matrix; entries can be
Python integers, `fractions.Fraction` values, or `qpa.Scalar` values from the
algebra's coefficient field. Matrices use right-module row actions.

```python
from fractions import Fraction

module = common.module_from_matrices(
    [1, 1, 1],
    [("a", [[Fraction(1, 2)]]), ("b", [[0]])],
)
assert module.dimension == 3
assert module.right_acting_algebra.is_identical(common)
```

The constructor checks dimensions and quotient relations. A failed relation
raises `RelationNotSatisfiedError`. Floating-point entries raise `TypeError`;
an incompatible scalar field or a fraction whose denominator is zero in the
coefficient field raises `InvalidInputError`. `module_from_integer_matrices`
remains available for matrices of machine-sized integer entries.

## Tensor products and enveloping algebras

Path, quotient, and common algebra handles provide `tensor_product(other)` and
`enveloping_algebra()`. The tensor method accepts a path algebra, quotient, or
common handle as its other factor. Both methods return a `TensorProductAlgebra`
that retains factor and inclusion-map data. Use its `algebra` property for
common algebra operations. The enveloping algebra follows QPA's factor order
`A.opposite()` then `A`:

```python
enveloping = algebra.enveloping_algebra()
enveloping_projectives = enveloping.algebra.indecomposable_projectives()
assert len(enveloping_projectives) == 9
assert enveloping.factor(2).is_identical(common)
```

The `quotient` property provides the concrete quotient presentation. For
example, `enveloping.quotient.dimension()` computes its dimension.

## Dimensions and termination

`QuotientAlgebra.dimension()` returns a decimal string for a known finite
dimension and `None` for an infinite-dimensional quotient. Convert a known finite result with `int` before using it in
arithmetic. Module dimensions are Python integers.

Quotient and Gröbner calculations can be very expensive or nonterminating on
mathematically difficult input. A finite-dimensional presentation is not by
itself a practical resource bound.

Algebra elements support `x ** n` for nonnegative integers, using the same
operation as `x.pow(n)` on concrete elements. For a quotient element, the
existing `x.pow_integer(n)` handles negative powers and can return `None` when
an inverse is unavailable. `**` deliberately does not change that return
contract. Use `x.scale_left(field.element(2))` for scalar multiplication;
Python integer or scalar multiplication syntax is not defined for algebra
elements or module maps.
