# Python quickstart

Install the package first with the [`uv` guide](python-installation.md). The
smallest complete computation constructs the truncated polynomial algebra
`F_5[x]/(x^3)`:

```python
import qpa

field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v"], [("v", "v", "x")])
algebra = qpa.QuotientAlgebra.truncated(field, quiver, 3)
x = algebra.generator("x")

assert algebra.dimension() == "3"
assert x.pow(3).is_zero
```

Run the checked version from this repository with:

```console
uv run python docs/examples/python/quickstart.py
```

`QuotientAlgebra.truncated` follows QPA's truncated-path-algebra construction;
it is not an independently selected quotient algorithm.

## Exact fields and scalars

The Python coefficient contract provides exact rationals, prime fields, and
canonical finite extensions:

```python
rationals = qpa.Field.rational()
two_thirds = rationals.element(2, 3)

field7 = qpa.Field.prime(7)
four = field7.element(11)  # 4 mod 7

field9 = qpa.Field.extension(3, 2)
z = field9.primitive_element()

assert two_thirds.numerator == "2"
assert two_thirds.denominator == "3"
assert four.representative == 4
assert z.coefficients == [0, 1]
assert field9.frobenius(z, 1).coefficients == [1, 2]
```

`Field.prime(p)` and `Field.extension(p, n)` verify the field descriptor.
Extension fields use GAP's Conway model, so their designated generators and
canonical subfield embeddings agree with QPA. Floating-point coefficients are
not accepted: construct exact integers, fractions, or polynomial-basis
coefficient vectors instead.

## Quivers and paths

```python
quiver = qpa.Quiver(
    ["v1", "v2", "v3"],
    [("v1", "v2", "a"), ("v2", "v3", "b")],
)

a = quiver.arrow_path("a")
b = quiver.arrow_path("b")
assert a * b == quiver.path(["a", "b"])
```

Path multiplication is written in traversal order: `a * b` means first `a`,
then `b`. This QPA convention is important when translating code from another
computer-algebra system.

## Where to continue

- [Quivers, algebras, and quotients](python-algebras.md)
- [Modules and homomorphisms](python-modules.md)
- [Homological calculations](python-homological.md)
- [Gröbner bases](python-groebner.md)
- [Python API and operation index](python-api.md)
