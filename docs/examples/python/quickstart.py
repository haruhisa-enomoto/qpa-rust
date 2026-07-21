"""Small first computation with the Python QPA binding."""

import qpa


field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v"], [("v", "v", "x")])

# The quotient is F_5[x]/(x^3), constructed through the translated QPA route.
algebra = qpa.QuotientAlgebra.truncated(field, quiver, 3)
x = algebra.generator("x")

assert algebra.dimension() == "3"
assert x.pow(3).is_zero

print(f"{algebra.view()}; dimension = {algebra.dimension()}")
