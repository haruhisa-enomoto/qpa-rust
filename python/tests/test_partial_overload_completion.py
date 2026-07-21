from __future__ import annotations

from fractions import Fraction

import pytest

import qpa


def _value(scalar: qpa.Scalar) -> Fraction | int:
    if scalar.representative is not None:
        return scalar.representative
    assert scalar.numerator is not None and scalar.denominator is not None
    return Fraction(int(scalar.numerator), int(scalar.denominator))


def _matrices(
    matrices: list[list[list[qpa.Scalar]]],
) -> list[list[list[Fraction | int]]]:
    return [
        [[_value(entry) for entry in row] for row in matrix]
        for matrix in matrices
    ]


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_checked_module_constructor_overloads_follow_qpa(field: qpa.Field) -> None:
    quiver = qpa.Quiver(
        ["v1", "v2", "v3"],
        [("v1", "v2", "a"), ("v2", "v3", "b")],
    )
    path = qpa.PathAlgebra(field, quiver)
    algebra = qpa.QuiverAlgebra.from_path_algebra(path)
    a_matrix = [[field.element(1), field.element(2)]]
    b_matrix = [[field.element(3)], [field.element(4)]]

    modules = [
        algebra.right_module_over_path_algebra_with_dimensions(
            [1, 2, 1], [("a", a_matrix), ("b", b_matrix)]
        ),
        algebra.right_module_over_path_algebra([a_matrix, b_matrix]),
        algebra.right_module_over_path_algebra_named(
            [("b", b_matrix), ("a", a_matrix)]
        ),
        algebra.right_module_over_path_algebra_indexed(
            [(quiver.arrow("b"), b_matrix), (quiver.arrow("a"), a_matrix)]
        ),
    ]
    expected = [
        [[_value(field.element(1)), _value(field.element(2))]],
        [[_value(field.element(3))], [_value(field.element(4))]],
    ]
    for module in modules:
        assert module is not None
        assert module.dimension_vector == [1, 2, 1]
        assert _matrices(module.arrow_matrices) == expected

    zero_source = algebra.right_module_over_path_algebra_named(
        [("a", [0, 2]), ("b", [2, 0])]
    )
    assert zero_source is not None
    assert zero_source.dimension_vector == [0, 2, 0]

    isolated_quiver = qpa.Quiver(
        ["i1", "i2", "i3"], [("i1", "i2", "c")]
    )
    isolated = qpa.QuiverAlgebra.from_path_algebra(
        qpa.PathAlgebra(field, isolated_quiver)
    ).right_module_over_path_algebra([[[field.element(1)]]])
    assert isolated is not None
    assert isolated.dimension_vector == [1, 1, 0]


def test_checked_constructor_rejects_foreign_field_and_returns_qpa_fail() -> None:
    field3 = qpa.Field.prime(3)
    field5 = qpa.Field.prime(5)
    path = qpa.PathAlgebra(
        field3, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    algebra = qpa.QuiverAlgebra.from_path_algebra(path)
    with pytest.raises(qpa.InvalidInputError, match=r"GF\(3\)"):
        algebra.right_module_over_path_algebra([[[field5.element(1)]]])

    bound_path = qpa.PathAlgebra(
        field3,
        qpa.Quiver(
            ["q1", "q2", "q3"],
            [("q1", "q2", "x"), ("q2", "q3", "y")],
        ),
    )
    quotient = bound_path.quotient_high_level(
        [bound_path.generator("x") * bound_path.generator("y")]
    )
    bound = qpa.QuiverAlgebra.from_quotient_algebra(quotient)
    one = field3.element(1)
    assert bound.right_module_over_path_algebra([[[one]], [[one]]]) is None


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_list_intersection_and_sum_follow_qpa_folds(field: qpa.Field) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    simple1, simple2 = path.simple_modules()
    target = qpa.MatrixModule.direct_sum([simple1, simple1, simple2])
    assert target is not None
    inclusions = target.direct_sum_inclusions()

    intersection = qpa.ModuleHomomorphism.intersection_list(inclusions)
    binary_intersection = inclusions[0].intersection(inclusions[1])[2].intersection(
        inclusions[2]
    )[2]
    assert intersection == binary_intersection
    assert intersection.source.is_zero
    assert intersection.range == target
    assert qpa.ModuleHomomorphism.intersection_list(
        [inclusions[0]]
    ).is_identical(inclusions[0])

    total = qpa.ModuleHomomorphism.sum_list(inclusions)
    binary_total = inclusions[0].sum(inclusions[1])[0].sum(inclusions[2])[0]
    assert total == binary_total
    assert total.source.dimension_vector == target.dimension_vector
    assert total.range == target
    assert qpa.ModuleHomomorphism.sum_list([inclusions[0]]).is_identical(
        inclusions[0]
    )

    with pytest.raises(qpa.InvalidInputError, match="empty"):
        qpa.ModuleHomomorphism.intersection_list([])
    with pytest.raises(qpa.InvalidInputError, match="empty"):
        qpa.ModuleHomomorphism.sum_list([])


def test_product_decomposition_supports_qpa_one_based_indexing() -> None:
    left = qpa.Quiver(["u1", "u2"], [("u1", "u2", "a")])
    right = qpa.Quiver(["v1", "v2"], [("v1", "v2", "b")])
    decomposition = left.product(right).decomposition
    assert decomposition[1].is_identical(left)
    assert decomposition[2].is_identical(right)
    with pytest.raises(qpa.InvalidInputError, match="outside 1..=2"):
        _ = decomposition[0]


def test_display_overloads_preserve_qpa_rendering() -> None:
    field = qpa.Field.rational()
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    module = path.module_from_integer_matrices([1, 1], [("a", [[1]])])
    module_print = (
        "<Module over <Rationals[<quiver with 2 vertices and 1 arrows>]>"
        " with dimension vector [ 1, 1 ]>"
    )
    assert module.view() == "<[ 1, 1 ]>"
    assert module.print_obj() == module_print
    assert module.display() == (
        f"{module_print} and linear maps given by\n"
        "for arrow a:\n"
        "[ [  1 ] ]\n"
    )

    identity = module.identity()
    map_print = f"<{module_print} ---> {module_print}>"
    assert identity.view() == "<<[ 1, 1 ]> ---> <[ 1, 1 ]>>"
    assert identity.print_obj() == map_print
    assert identity.display() == (
        f"{map_print}\n"
        "with linear map for vertex number 1:\n"
        "[ [  1 ] ]\n"
        "linear map for vertex number 2:\n"
        "[ [  1 ] ]\n"
    )

    decomposition = path.quiver.product(path.quiver).decomposition
    assert decomposition.view() == "<quiver product decomposition>"
    assert decomposition.print_obj() == "<quiver product decomposition>"
    assert decomposition.display() == "<quiver product decomposition>"


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_loewy_length_algebra_overload_follows_qpa(field: qpa.Field) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    assert path.loewy_length == 2
    assert qpa.QuiverAlgebra.from_path_algebra(path).loewy_length == 2

    quotient = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["w"], [("w", "w", "x")]), 2
    )
    assert quotient.loewy_length == 2
    assert qpa.QuiverAlgebra.from_quotient_algebra(quotient).loewy_length == 2


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_finite_type_decisions_preserve_qpa_values_and_messages(
    field: qpa.Field,
) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    path_decision = path.finite_type_decision()
    assert path.is_finite_type_algebra
    assert path_decision.value is True
    assert path_decision.reason == "DynkinComponents"
    assert path_decision.messages == [
        "A_2\n",
        "Finite type!\nQuiver is a (union of) Dynkin quiver(s).\n",
    ]

    quotient = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["w"], [("w", "w", "x")]), 2
    )
    quotient_decision = quotient.finite_type_decision()
    assert quotient.is_finite_type_algebra is True
    assert quotient_decision.value is True
    assert quotient_decision.reason == "RadicalSquareZero"
    assert quotient_decision.messages == [
        "Quiver contains an (un)oriented cycle.\n",
        "A_2\n",
    ]


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_uniform_element_projective_family_preserves_qpa_behavior(
    field: qpa.Field,
) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    v1, v2, arrow = path.generators
    module = path.right_projective_module_from_uniform_elements([arrow])
    assert isinstance(module, qpa.RightUniformPathAlgebraModule)
    assert module.algebra.is_identical(path)
    assert module.elements == [arrow]
    assert module.is_vertex_projective

    generator = module.generators[0]
    zero = module.zero()
    assert isinstance(generator, qpa.RightUniformPathAlgebraVector)
    assert generator.module.is_identical(module)
    assert generator.components == [arrow]
    assert generator.leading_position == 1
    assert generator.is_normal_form and not generator.is_zero
    assert generator.leading_term == generator
    assert generator.leading_coefficient == field.element(1)
    assert generator.leading_component == arrow
    assert generator.source_vertex == v1.leading_monomial
    assert generator.target_vertex == v2.leading_monomial
    assert generator.is_left_uniform
    assert generator.is_right_uniform
    assert generator.is_uniform
    assert generator.compare(generator) == 0
    assert generator.is_left_divisible(generator)

    assert zero.leading_position == 2 and zero.is_zero
    assert generator.act(v1).is_zero
    negative = -generator
    assert generator + negative == zero
    assert generator.scale_left(field.element(0)) == zero
    assert generator.scale_right(field.element(1)) == generator
    expected_two = 0 if field.characteristic == 2 else 2
    assert generator.scale_left(field.element(2)).components == [
        arrow.scale_left(field.element(expected_two))
    ]

    with pytest.raises(qpa.InvalidInputError, match="no vertices component"):
        generator.act(v2)
    with pytest.raises(qpa.InvalidInputError, match="no vertices component"):
        module.vector([arrow])
    with pytest.raises(qpa.InvalidInputError, match="no vertices component"):
        module.canonical_basis()
    with pytest.raises(qpa.InvalidInputError, match="no vertices component"):
        module.uniform_generators()

    duplicate = qpa.RightUniformPathAlgebraModule.from_elements(path, [arrow])
    assert duplicate.elements == module.elements
    assert not duplicate.is_identical(module)
    with pytest.raises(qpa.InvalidInputError, match="not uniform"):
        path.right_projective_module_from_uniform_elements([v1 + v2])
    with pytest.raises(qpa.InvalidInputError, match="cannot extract coefficients"):
        path.is_left_uniform_list_at_path([arrow], v1.leading_monomial)
    with pytest.raises(qpa.InvalidInputError, match="cannot extract coefficients"):
        path.is_right_uniform_list_at_path([arrow], v1.leading_monomial)


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_quotient_canonical_basis_object_and_coefficients_follow_qpa(
    field: qpa.Field,
) -> None:
    quotient = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver(["w"], [("w", "w", "x")]), 2
    )
    basis = quotient.canonical_basis_object()
    assert isinstance(basis, qpa.QuotientCanonicalBasis)
    assert basis.quotient.is_identical(quotient)
    assert basis.has_basis_vectors
    assert basis.basis_vectors == quotient.canonical_basis()
    assert basis.is_identical(quotient.canonical_basis_object())

    element = quotient.one().scale_left(field.element(2)) + quotient.generator(
        "x"
    ).scale_left(field.element(3))
    assert [_value(value) for value in basis.coefficients(element)] == [
        _value(field.element(2)),
        _value(field.element(3)),
    ]
