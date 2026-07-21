# Gröbner bases in Python

The default Python routes preserve QPA's ordering, uniformization, and GBNP
call structure. Numerica supplies exact field arithmetic; it does not replace
the QPA/GBNP algorithm.

## Compute a basis and a remainder

```python
import qpa

field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v"], [("v", "v", "x"), ("v", "v", "y")])
path_algebra = qpa.PathAlgebra(field, quiver)
x = path_algebra.generator("x")
y = path_algebra.generator("y")

relations = [x.pow(2), x * y - y * x, y.pow(2)]
basis = path_algebra.groebner_basis(relations)
remainder = path_algebra.remainder_of_division(y * x, basis.relations)

assert len(basis.relations) == 3
assert remainder == x * y
```

The runnable file is
[`docs/examples/python/groebner.py`](../examples/python/groebner.py).

`PathAlgebra.groebner_basis` is QPA's translated high-level Gröbner route.
`PathAlgebra.gbnp_groebner_basis` and
`PathAlgebra.gbnp_groebner_basis_nc` expose the corresponding checked and
unchecked GBNP bridges. The `nc` spelling means that QPA's applicability
checks are intentionally omitted; it does not mean commutative algebra.

## Inspecting a basis

`GroebnerBasis.relations` is the structured relation list. Iteration follows
the QPA enumerator order:

```python
for relation in basis:
    print(relation.print_obj())

assert basis.is_complete
assert basis.kind in {"two-sided", "right"}
```

Use structured paths and coefficients for subsequent computation. `print_obj`
and `view` reproduce source-compatible display forms where tested, but they
are not serialization formats.

## Ordering and resource behavior

Do not reorder a returned basis merely because another ordering spans the same
ideal. Tip order and output order are observable inputs to later QPA methods.

Noncommutative Gröbner completion may consume unbounded time and memory. The
binding deliberately does not hide that behavior behind a different
terminating heuristic. Start with small examples and use QPA's existing
preconditions and bounded operations where available.
