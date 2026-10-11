"""Public-API regression checks for the native performance pass."""
import gc

import pytest
import qpa


@pytest.mark.parametrize('field', [qpa.Field.rational(), qpa.Field.prime(2),
                                   qpa.Field.prime(101), qpa.Field.extension(2, 2)])
def test_selective_projectives_cover_blocks_and_lazy_maps(field):
    parent = qpa.PathAlgebra(field, qpa.Quiver.indexed(3, [
        (1, 1, 'r'), (1, 1, 's'), (1, 2, 't'), (2, 3, 'u')]))
    r, s, t, u = [parent.generator(name) for name in ['r', 's', 't', 'u']]
    scalar = field.primitive_element() if (field.degree or 1) > 1 else field.element(1)
    algebra = parent.quotient([r*r, s*s, s*r - (r*s).scale_left(scalar), t*u])
    simples = algebra.simple_modules()
    selected = simples[0].basis_elements()[0].hom_from_projective()
    held = selected.source
    family = algebra.indecomposable_projectives()
    assert held.is_identical(family[0])
    del family
    gc.collect()
    assert held.is_identical(algebra.indecomposable_projectives()[0])
    source = qpa.MatrixModule.direct_sum([simples[0], simples[1], simples[0], algebra.zero_module()])
    cover = source.projective_cover()
    assert cover.is_surjective and cover.range.is_identical(source)
    assert source.projective_cover().is_identical(cover)
    maps = [g.hom_from_projective() for g in source.minimal_generating_set()]
    reference_sum = qpa.MatrixModule.direct_sum([m.source for m in maps])
    composites = [p.compose(m) for p, m in zip(reference_sum.direct_sum_projections(), maps)]
    reference = composites[0]
    for other in composites[1:]:
        reference = reference.add(other)
    assert cover.vertex_matrices == reference.vertex_matrices
    assert cover.source.dimension_vector == reference_sum.dimension_vector
    assert cover.source.arrow_matrices == reference_sum.arrow_matrices
    for summand, mapping in zip(cover.source.direct_sum_summands(), maps):
        assert summand.is_identical(mapping.source)
    projections = cover.source.direct_sum_projections()
    inclusions = cover.source.direct_sum_inclusions()
    assert all(p.is_identical(q) for p, q in zip(projections, cover.source.direct_sum_projections()))
    for i, inclusion in enumerate(inclusions):
        for j, projection in enumerate(projections):
            composite = inclusion.compose(projection)
            if i == j:
                assert composite.vertex_matrices == inclusion.source.identity().vertex_matrices
            else:
                assert composite.is_zero
    kernel = cover.kernel_inclusion()
    assert kernel.is_injective and kernel.compose(cover).is_zero
    assert cover.source.dimension == source.dimension + kernel.source.dimension


@pytest.mark.parametrize('field', [qpa.Field.rational(), qpa.Field.prime(2),
                                   qpa.Field.prime(101), qpa.Field.extension(2, 2)])
def test_checked_relations_still_reject_inputs_and_validate_direct_sum(field):
    parent = qpa.PathAlgebra(field, qpa.Quiver.indexed(1, [(1, 1, 'x')]))
    x = parent.generator('x')
    algebra = parent.quotient([x*x*x])
    module = algebra.module_from_integer_matrices([3], [('x', [[0, 1, 0], [0, 0, 1], [0, 0, 0]])])
    direct_sum = qpa.MatrixModule.direct_sum([module, algebra.zero_module(), module])
    checked = algebra.module_from_matrices(direct_sum.dimension_vector, [('x', direct_sum.arrow_matrices[0])])
    assert checked.arrow_matrices == direct_sum.arrow_matrices
    with pytest.raises(qpa.RelationNotSatisfiedError):
        algebra.module_from_integer_matrices([1], [('x', [[1]])])
