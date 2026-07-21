"""Quivers, paths, relations, and a quotient algebra."""

import qpa


quiver = qpa.Quiver(
    ["v1", "v2", "v3"],
    [("v1", "v2", "a"), ("v2", "v3", "b")],
)

# Path multiplication follows traversal order: first a, then b.
path_ab = quiver.path(["a", "b"])
assert quiver.arrow_path("a") * quiver.arrow_path("b") == path_ab

path_algebra = qpa.PathAlgebra(qpa.Field.rational(), quiver)
a = path_algebra.generator("a")
b = path_algebra.generator("b")

# PathAlgebra.quotient uses QPA's configured default GBNP route.
algebra = path_algebra.quotient([a * b])
assert algebra.dimension() == "5"
assert (algebra.generator("a") * algebra.generator("b")).is_zero

print(f"relation: {path_ab} = 0; dimension = {algebra.dimension()}")
