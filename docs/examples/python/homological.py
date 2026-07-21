"""Syzygies, Ext, and a lazy projective resolution."""

import qpa


field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
algebra = qpa.PathAlgebra(field, quiver)
simple_at_v1, simple_at_v2 = algebra.simple_modules()

syzygy = simple_at_v1.first_syzygy()
assert syzygy.dimension_vector == [0, 1]
assert simple_at_v1.projective_dimension(3) == 1
assert simple_at_v1.ext1(simple_at_v2).dimension == 1

# The QPA resolution includes the resolved module in degree -1. The first
# projective and its kernel occur in degrees 0 and 1.
resolution = simple_at_v1.projective_resolution()
assert resolution.object(0).dimension_vector == [1, 1]
assert resolution.object(1).dimension_vector == [0, 1]
assert resolution.object(2).is_zero

print(f"first syzygy: {syzygy.dimension_vector}")
print(f"Ext^1 dimension: {simple_at_v1.ext1(simple_at_v2).dimension}")
