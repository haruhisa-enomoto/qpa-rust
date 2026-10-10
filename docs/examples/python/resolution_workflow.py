"""A bounded resolution, Hom/Ext dimensions, and an exact complex."""

import qpa


field = qpa.Field.prime(2)
parent = qpa.PathAlgebra(field, qpa.Quiver.indexed(1, [(1, 1, "x")]))
x = parent.generator("x")
algebra = parent.quotient([x * x])
simple = algebra.simple_modules()[0]

# The lazy resolution includes the augmentation target in degree -1.
resolution = simple.projective_resolution()
assert resolution.object(-1).is_identical(simple)
terms = [resolution.object(degree) for degree in range(3)]
differentials = [resolution.differential(degree) for degree in range(3)]
assert [term.dimension for term in terms] == [2, 2, 2]
for degree, differential in enumerate(differentials):
    assert differential.source.is_identical(resolution.object(degree))
    assert differential.range.is_identical(resolution.object(degree - 1))
    if degree:
        # compose follows application order: first d_degree, then d_(degree-1).
        assert differential.compose(differentials[degree - 1]).is_zero

# For the dual numbers, each projective top has dimension equal to its free rank.
betti = [term.top().dimension for term in terms]
assert betti == [1, 1, 1]
print("projective dimensions:", [term.dimension for term in terms])
print("Betti numbers:", betti)
print("d_1:", differentials[1].underlying_linear_matrix())

hom = simple.hom_basis_with_coefficients_to(simple)
assert hom is not None
assert hom.dimension == simple.hom_dimension_to(simple) == 1
assert len(hom.basis) == len(simple.hom_basis_to(simple)) == 1
ext = simple.ext1(simple)
assert ext.dimension == simple.ext1_dimension(simple) == len(ext.basis) == 1
# Batches vary the coefficient module, not the cohomological degree.
assert simple.ext1_dimensions([simple, terms[0]]) == [1, 0]
print("Hom and Ext^1 dimensions:", hom.dimension, ext.dimension)

cover = differentials[0]
kernel = cover.kernel()  # A module.
inclusion = cover.kernel_inclusion()  # A map from that module into P_0.
assert inclusion.source.isomorphic_to(kernel)
assert inclusion.range.is_identical(cover.source)
assert inclusion.compose(cover).is_zero

# List entries are d_0, d_1: P_0 -> S and kernel -> P_0.
exact = simple.finite_complex(0, [cover, inclusion])
assert [exact.object(degree).dimension for degree in (-1, 0, 1)] == [1, 2, 1]
assert all(exact.homology(degree).is_zero for degree in (-1, 0, 1))
assert exact.is_short_exact_sequence
checked_exact = inclusion.short_exact_sequence_with(cover)
assert checked_exact.is_short_exact_sequence
print("short exact sequence:", checked_exact.is_short_exact_sequence)
