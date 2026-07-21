"""A right matrix module and a Hom-space computation."""

import qpa


field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
algebra = qpa.PathAlgebra(field, quiver)

# Right modules use row vectors. For a: v1 -> v2, this is a 1-by-1 matrix.
module = algebra.module_from_integer_matrices([1, 1], [("a", [[1]])])
simple_at_v1, simple_at_v2 = algebra.simple_modules()

assert module.dimension_vector == [1, 1]
assert module.arrow_matrix_shapes == [(1, 1)]
assert algebra.hom(simple_at_v1, module).dimension == 0
assert simple_at_v2.hom_basis_to(module) is not None

print(f"dimension vector: {module.dimension_vector}")
print(f"radical series: {module.radical_series()}")
