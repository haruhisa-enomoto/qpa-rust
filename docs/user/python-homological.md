# Homological calculations in Python

The homological Python methods call the translated QPA module and complex
routines. Lazy resolutions preserve QPA's degree conventions and caches.

## Syzygies and projective dimension

```python
import qpa

field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
algebra = qpa.PathAlgebra(field, quiver)
simple_at_v1, simple_at_v2 = algebra.simple_modules()

syzygy = simple_at_v1.first_syzygy()
assert syzygy.dimension_vector == [0, 1]
assert simple_at_v1.projective_dimension(3) == 1
```

The bound in `projective_dimension(bound)` is a computational bound, not a
claim that every module has projective dimension at most that value. A `None`
result means that the translated method did not establish a finite answer
within the requested bound.

## Ext

`module.ext1(coefficients_in)` returns a structured result containing the QPA
basis and its dimension:

```python
ext = simple_at_v1.ext1(simple_at_v2)
assert ext.dimension == 1
assert len(ext.basis) == 1
```

The argument order matches QPA's `ExtOverAlgebra` route as exposed by the
binding. Consult the operation index before translating a convention from a
different package.

## Lazy projective resolutions

```python
resolution = simple_at_v1.projective_resolution()

assert resolution.object(0).dimension_vector == [1, 1]
assert resolution.object(1).dimension_vector == [0, 1]
assert resolution.object(2).is_zero
```

QPA includes the resolved module in degree `-1`; `object(0)` is its first
projective. Calling `object(n)` or `differential(n)` extends the lazy
resolution as needed and reuses cached terms.

For complexes, use `MatrixModule.projective_resolution_complex()` or the
explicit `ModuleComplex` constructors. Homology is obtained with
`complex.homology(degree)`.

The checked complete example is
[`docs/examples/python/homological.py`](../examples/python/homological.py).
