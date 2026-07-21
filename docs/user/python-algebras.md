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

## Dimensions and termination

`QuotientAlgebra.dimension()` returns a decimal string for a known finite
dimension and `None` when QPA does not establish a finite value through the
applicable route. Strings preserve arbitrary-precision dimensions across the
Python boundary.

Quotient and Gröbner calculations can be very expensive or nonterminating on
mathematically difficult input. A finite-dimensional presentation is not by
itself a practical resource bound.
