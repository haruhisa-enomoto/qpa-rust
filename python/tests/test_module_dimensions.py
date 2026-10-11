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


@pytest.mark.parametrize('family', ['infinite_path', 'nonadmissible_quotient'])
def test_ext_dimension_zero_source_preserves_unsupported_algebra_boundaries(family):
    parent = qpa.PathAlgebra(qpa.Field.rational(), qpa.Quiver.indexed(1, [(1, 1, 'x')]))
    if family == 'infinite_path':
        algebra = parent
    else:
        x = parent.generator('x')
        algebra = parent.quotient([x*x - x])
        assert not algebra.is_admissible
    zero = algebra.zero_module()
    module = algebra.module_from_integer_matrices([1], [('x', [[0]])])
    assert zero.ext1_dimension(module) == zero.ext1(module).dimension == 0
    assert zero.ext1_dimensions([module, zero, module]) == [0, 0, 0]
    # Do not introduce a target-zero shortcut which suppresses cover errors.
    for query in [lambda: module.ext1_dimension(zero), lambda: module.ext1(zero),
                  lambda: module.ext1_dimensions([zero])]:
        with pytest.raises(qpa.InvalidInputError):
            query()


def test_unconstrained_hom_dimension_does_not_materialize_unknown_rows():
    algebra = qpa.QuiverAlgebra.from_path_algebra(qpa.PathAlgebra(
        qpa.Field.rational(), qpa.Quiver.indexed(1, [])))
    module = algebra.right_module_over_path_algebra_with_dimensions([100_000], [])
    assert module.hom_dimension_to(module) == 10_000_000_000


@pytest.mark.parametrize('field', FIELDS)
def test_sparse_hom_equations_preserve_regular_kronecker_ext(field):
    # 800 unknowns and 800 constraints: above the small-dense cutoff.
    # End(R_n) = k[t]/(t^n), and Ext^1(R_n,R_n) also has dimension n.
    n = 20
    algebra = qpa.QuiverAlgebra.from_path_algebra(qpa.PathAlgebra(
        field, qpa.Quiver.indexed(2, [(1, 2, 'a'), (1, 2, 'b')])))
    zero, one = field.element(0), field.element(1)
    eigenvalue = field.primitive_element() if (field.degree or 1) > 1 else zero
    identity = [[one if i == j else zero for j in range(n)] for i in range(n)]
    jordan = [[eigenvalue if i == j else one if j == i + 1 else zero
               for j in range(n)] for i in range(n)]
    module = algebra.right_module_over_path_algebra_with_dimensions(
        [n, n], [('a', identity), ('b', jordan)])
    assert module.hom_dimension_to(module) == n
    assert len(module.hom_basis_to(module)) == n
    assert module.ext1_dimension(module) == n


@pytest.mark.parametrize('field', FIELDS)
def test_subrepresentation_discards_exponentially_redundant_path_images(field):
    n = 24
    parent = qpa.PathAlgebra(field, qpa.Quiver.indexed(1, [(1, 1, 'x'), (1, 1, 'y')]))
    jordan = [[int(j == i + 1) for j in range(n)] for i in range(n)]
    module = parent.module_from_integer_matrices([n], [('x', jordan), ('y', jordan)])
    # All 2**length words give the same vector until length n. The full
    # subrepresentation is the original module in its canonical basis.
    inclusion = module.subrepresentation_inclusion([module.basis_elements()[0]])
    assert inclusion.source == module
    assert inclusion.range.is_identical(module)
    assert inclusion.vertex_matrices == module.identity().vertex_matrices


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
