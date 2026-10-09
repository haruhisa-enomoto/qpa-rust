"""Dimension-only paths agree with the QPA basis constructions across fields."""

import pytest
import qpa


FIELDS = [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5), qpa.Field.extension(2, 2)]


@pytest.mark.parametrize('field', FIELDS)
@pytest.mark.parametrize('family', ['bound_quiver', 'truncated_loop'])
def test_dimensions_and_batch_match_basis_constructions(field, family):
    if family == 'bound_quiver':
        parent = qpa.PathAlgebra(field, qpa.Quiver.indexed(3, [(1, 2, 'a'), (2, 3, 'b')]))
        algebra = parent.quotient([parent.generator('a') * parent.generator('b')])
    else:
        algebra = qpa.QuotientAlgebra.truncated(field, qpa.Quiver.indexed(1, [(1, 1, 'x')]), 3)
    simples = algebra.simple_modules()
    projectives = algebra.indecomposable_projectives()
    modules = [algebra.zero_module(), *simples, *projectives,
               qpa.MatrixModule.direct_sum([*simples, *projectives])]
    if family == 'truncated_loop':
        modules.append(algebra.module_from_integer_matrices([2], [('x', [[0, 1], [0, 0]])]))
    # Repeated targets and repeated batches must not disturb caches or order.
    targets = [*modules, modules[1], modules[0]]
    for source in modules:
        expected = [source.ext1(target).dimension for target in targets]
        assert source.ext1_dimensions([]) == []
        assert source.ext1_dimensions(targets) == expected
        assert source.ext1_dimensions(targets) == expected
        for target, dimension in zip(targets, expected):
            assert source.hom_dimension_to(target) == len(source.hom_basis_to(target))
            assert source.ext1_dimension(target) == dimension


def test_nonprime_extension_coefficients_are_not_binary():
    field = qpa.Field.extension(2, 2)
    algebra = qpa.QuiverAlgebra.from_path_algebra(qpa.PathAlgebra(
        field, qpa.Quiver.indexed(2, [(1, 2, 'a'), (1, 2, 'b')])))
    z = field.primitive_element()
    one = field.element(1)
    # Two nonisomorphic regular simples over GF(4).
    left = algebra.right_module_over_path_algebra_with_dimensions([1, 1], [('a', [[one]]), ('b', [[z]])])
    right = algebra.right_module_over_path_algebra_with_dimensions([1, 1], [('a', [[one]]), ('b', [[field.element_from_coefficients([1, 1])]])])
    for source in [left, right]:
        assert source.hom_dimension_to(source) == source.ext1_dimension(source) == 1
        other = right if source is left else left
        assert source.hom_dimension_to(other) == source.ext1_dimension(other) == 0
        assert source.ext1_dimensions([source, other]) == [1, 0]


def test_dimension_error_boundaries_match_existing_apis():
    quiver = qpa.Quiver.indexed(1, [])
    left = qpa.PathAlgebra(qpa.Field.rational(), quiver)
    right = qpa.PathAlgebra(qpa.Field.rational(), qpa.Quiver.indexed(2, []))
    binary = qpa.PathAlgebra(qpa.Field.prime(2), quiver)
    source, other = left.zero_module(), right.zero_module()
    assert source.hom_basis_to(other) is None
    assert source.hom_dimension_to(other) is None
    for query in [lambda: source.ext1_dimension(other),
                  lambda: source.ext1_dimensions([source, other]),
                  lambda: source.hom_dimension_to(binary.zero_module()),
                  lambda: source.ext1_dimension(binary.zero_module()),
                  lambda: source.ext1_dimensions([binary.zero_module()])]:
        with pytest.raises(qpa.InvalidInputError):
            query()


def test_binary_ext_grid_with_nonzero_restriction_image():
    algebra = qpa.QuotientAlgebra.truncated(
        qpa.Field.prime(2), qpa.Quiver.indexed(1, [(1, 1, 'x'), (1, 1, 'y')]), 2
    )
    simple = algebra.simple_modules()[0]
    projective = algebra.indecomposable_projectives()[0]
    source = qpa.MatrixModule.direct_sum([simple] * 4 + [projective] * 2)
    targets = [qpa.MatrixModule.direct_sum([simple] * n + [projective] * p)
               for n, p in [(4, 0), (2, 2), (0, 4), (3, 1), (1, 3), (2, 0)]]
    # dim Ext(S,S) = 2 and dim Ext(S,P) = 3 for this two-loop algebra.
    expected = [32, 40, 48, 36, 44, 16]
    assert [source.ext1(target).dimension for target in targets] == expected
    assert [source.ext1_dimension(target) for target in targets] == expected
    assert source.ext1_dimensions(targets) == expected
