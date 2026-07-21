from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _value(scalar: qpa.Scalar | int) -> Fraction | int:
    if isinstance(scalar, int):
        return scalar
    if scalar.representative is not None:
        return scalar.representative
    assert scalar.numerator is not None and scalar.denominator is not None
    return Fraction(int(scalar.numerator), int(scalar.denominator))


def _rows(rows: list[list[qpa.Scalar]]) -> list[list[Fraction | int]]:
    return [[_value(entry) for entry in row] for row in rows]


def _oracle_rows(rows: list[list[str]]) -> list[list[Fraction | int]]:
    return [
        [
            0
            if entry.startswith("0*")
            else 1
            if entry.startswith("Z(")
            else Fraction(entry)
            for entry in row
        ]
        for row in rows
    ]


def _matrices(
    matrices: list[list[list[qpa.Scalar]]],
) -> list[list[list[Fraction | int]]]:
    return [_rows(matrix) for matrix in matrices]


def _length_data(module: qpa.MatrixModule) -> dict[str, object]:
    matrices = module.matrices_of_path_algebra_module()
    return {
        "dimension": module.dimension_vector,
        "matrices": [
            [] if matrix == 0 else _rows(matrix)
            for matrix in matrices
        ],
    }


def _module_data(module: qpa.MatrixModule) -> dict[str, object]:
    return {
        "dimension": module.dimension_vector,
        "matrices": _matrices(module.arrow_matrices),
    }


@pytest.mark.parametrize(
    ("prime", "fixture_key"),
    [(2, "binary"), (5, "quinary")],
)
def test_length_enumeration_operations_follow_gap_qpa_oracle(
    prime: int, fixture_key: str
) -> None:
    expected = json.loads(
        (FIXTURES / "modulehom_length.schema-v1.json").read_text(encoding="utf-8")
    )[fixture_key]
    field = qpa.Field.prime(prime)
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    algebra = qpa.QuiverAlgebra.from_path_algebra(path)
    simples = path.simple_modules()

    assert [
        _length_data(module)
        for module in qpa.MatrixModule.all_modules_of_length_plus_one(simples)
    ] == expected["a2_plus"]
    assert [
        _length_data(module)
        for module in algebra.all_modules_of_length_at_most(2)
    ] == expected["a2_at_most_two"]
    assert [
        _length_data(module)
        for module in algebra.all_indec_modules_of_length_at_most(2)
    ] == expected["a2_indec_two"]
    assert algebra.all_modules_of_length_at_most(-1) == expected["a2_negative"]
    assert [
        _length_data(module)
        for module in algebra.all_modules_of_length_at_most(0)
    ] == expected["a2_zero"]

    loop = qpa.PathAlgebra(field, qpa.Quiver(["w"], [("w", "w", "x")]))
    x = loop.generator("x")
    dual = loop.quotient_high_level([x * x])
    length_two = qpa.MatrixModule.all_modules_of_length_plus_one(
        dual.simple_modules()
    )
    length_three = qpa.MatrixModule.all_modules_of_length_plus_one(length_two)
    assert [_length_data(module) for module in length_two] == expected[
        "dual_length_two"
    ]
    assert [_length_data(module) for module in length_three] == expected[
        "dual_length_three"
    ]


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2)])
def test_projective_basis_and_element_conversions_follow_gap_qpa_oracle(
    field: qpa.Field,
) -> None:
    quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    path = qpa.PathAlgebra(field, quiver)
    algebra = qpa.QuiverAlgebra.from_path_algebra(path)

    basis = algebra.basis_of_projectives()
    assert basis is not None
    assert [
        [
            [element.terms[0][0].external_representation for element in component]
            for component in projective
        ]
        for projective in basis
    ] == [[[[1]], [[3]]], [[], [[2]]]]

    vertex = algebra.include_path(quiver.vertex_path("v1"))
    arrow = algebra.include_path(quiver.arrow_path("a"))
    projective_element = algebra.element_as_indec_projective(vertex + arrow)
    assert projective_element.module.dimension_vector == [1, 1]
    assert projective_element.support_vertices == [1, 2]
    assert [_value(entry) for entry in projective_element.basis_coefficients] == [1, 1]
    recovered = projective_element.element_in_indec_projective(algebra, 1)
    assert [
        (term.external_representation, _value(coefficient))
        for term, coefficient in recovered.terms
    ] == [([1], 1), ([3], 1)]


def test_complement_in_full_row_space_follows_gap_qpa_oracle() -> None:
    expected = json.loads(
        (FIXTURES / "moduledecomp_foundation.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )["complement"]
    rationals = qpa.Field.rational()
    binary = qpa.Field.prime(2)
    quinary = qpa.Field.prime(5)

    rational_rows = [[1, 2, 0, 1], [0, 0, 1, 3]]
    rational_matrix = [
        [rationals.element(entry) for entry in row] for row in rational_rows
    ]
    binary_rows = [[1, 1, 0, 1], [0, 0, 1, 1]]
    binary_matrix = [[binary.element(entry) for entry in row] for row in binary_rows]
    quinary_identity = [
        [quinary.element(int(row == column)) for column in range(3)]
        for row in range(3)
    ]

    assert _rows(qpa.complement_in_full_row_space(rational_matrix, rationals)) == expected[
        "rational"
    ]
    assert _rows(qpa.complement_in_full_row_space(binary_matrix, binary)) == expected["f2"]
    assert _rows(qpa.complement_in_full_row_space(quinary_identity, quinary)) == expected[
        "full_f5"
    ]


@pytest.mark.parametrize(
    ("prime", "fixture_key"),
    [(2, "binary"), (5, "quinary")],
)
def test_alternative_decompositions_follow_gap_qpa_oracle(
    prime: int, fixture_key: str
) -> None:
    expected = json.loads(
        (FIXTURES / "moduledecomp_modules.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )[fixture_key]
    field = qpa.Field.prime(prime)
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    mixed = path.module_from_integer_matrices([2, 1], [("a", [[1], [0]])])
    vector = qpa.PathAlgebra(field, qpa.Quiver(["v"], []))
    repeated = vector.module_from_integer_matrices([2], [])
    simple = vector.module_from_integer_matrices([1], [])

    assert [_module_data(module) for module in mixed.decompose_via_top()] == expected[
        "via_top"
    ]["mixed"]
    assert [_module_data(module) for module in repeated.decompose_via_top()] == expected[
        "via_top"
    ]["repeated"]
    assert [
        _module_data(module) for module in repeated.decompose_via_char_poly()
    ] == expected["char_poly"]["repeated"]
    assert [_module_data(module) for module in simple.decompose_via_char_poly()] == expected[
        "char_poly"
    ]["simple"]

    repeated_basis = repeated.hom_basis_to(repeated)
    simple_basis = simple.hom_basis_to(simple)
    assert repeated_basis is not None and simple_basis is not None
    repeated_result = repeated.decompose_probabilistic(repeated_basis)
    simple_result = simple.decompose_probabilistic(simple_basis)
    assert isinstance(repeated_result, list)
    assert isinstance(simple_result, qpa.MatrixModule)
    assert [_module_data(module) for module in repeated_result] == expected[
        "probabilistic"
    ]["repeated"]["modules"]
    assert [_module_data(simple_result)] == expected["probabilistic"]["simple"][
        "modules"
    ]


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2)])
def test_based_module_conversion_and_nc_constructor_follow_qpa(
    field: qpa.Field,
) -> None:
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    module = path.module_from_integer_matrices(
        [2, 2], [("a", [[1, 1], [0, 1]])]
    )
    based = module.as_based_right_algebra_module()
    rebuilt_based = qpa.BasedRightAlgebraModule.from_actions(
        based.algebra,
        based.dimension,
        based.vertex_actions,
        based.arrow_actions,
    )
    for converted in (
        based.to_path_algebra_matrix_module(),
        rebuilt_based.to_path_algebra_matrix_module(),
    ):
        assert converted.dimension_vector == [2, 2]
        assert _matrices(converted.arrow_matrices) == [[[1, 1], [0, 1]]]

    one = field.element(1)
    unchecked = qpa.QuiverAlgebra.from_path_algebra(
        path
    ).right_module_over_path_algebra_nc([[[one]]])
    assert unchecked is not None
    assert _module_data(unchecked) == {
        "dimension": [1, 1],
        "matrices": [[[1]]],
    }

    chain_quiver = qpa.Quiver(
        ["c1", "c2", "c3"],
        [("c1", "c2", "a"), ("c2", "c3", "b")],
    )
    chain_path = qpa.PathAlgebra(field, chain_quiver)
    chain = qpa.QuiverAlgebra.from_path_algebra(chain_path)
    a_matrix = [[field.element(1), field.element(2)]]
    b_matrix = [[field.element(3)], [field.element(4)]]
    positional = chain.right_module_over_path_algebra_nc([a_matrix, b_matrix])
    named = chain.right_module_over_path_algebra_nc_named(
        [("b", b_matrix), ("a", a_matrix)]
    )
    indexed = chain.right_module_over_path_algebra_nc_indexed(
        [
            (chain_quiver.arrow("b"), b_matrix),
            (chain_quiver.arrow("a"), a_matrix),
        ]
    )
    for inferred in (positional, named, indexed):
        assert inferred is not None
        assert _module_data(inferred) == {
            "dimension": [1, 2, 1],
            "matrices": [
                [[_value(field.element(1)), _value(field.element(2))]],
                [[_value(field.element(3))], [_value(field.element(4))]],
            ],
        }

    zero_source = qpa.QuiverAlgebra.from_path_algebra(
        path
    ).right_module_over_path_algebra_nc_named([("a", [0, 2])])
    assert zero_source is not None
    assert zero_source.dimension_vector == [0, 2]

    isolated_path = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["i1", "i2", "i3"],
            [("i1", "i2", "c")],
        ),
    )
    with pytest.raises(qpa.InvalidInputError, match="isolated vertices"):
        qpa.QuiverAlgebra.from_path_algebra(
            isolated_path
        ).right_module_over_path_algebra_nc([[[one]]])


def test_nc_constructor_preserves_qpa_foreign_odd_prime_acceptance() -> None:
    field3 = qpa.Field.prime(3)
    field5 = qpa.Field.prime(5)
    path = qpa.PathAlgebra(
        field3, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    module = qpa.QuiverAlgebra.from_path_algebra(
        path
    ).right_module_over_path_algebra_nc([[[field5.element(0)]]])
    assert module is not None
    assert module.dimension_vector == [1, 1]
    assert module.arrow_matrices[0][0][0].field == field5


def test_bongartz_and_random_module_follow_gap_qpa_oracles() -> None:
    expected_bongartz = json.loads(
        (FIXTURES / "reptype_bongartz.schema-v1.json").read_text(encoding="utf-8")
    )
    rationals = qpa.Field.rational()
    loop = qpa.PathAlgebra(
        rationals, qpa.Quiver(["v"], [("v", "v", "x")])
    )
    x = loop.generator("x")
    dual = qpa.QuiverAlgebra.from_quotient_algebra(
        loop.quotient_high_level([x * x])
    )
    for bound, key in [(0, "bound_zero"), (1, "bound_one"), (-1, "negative_bound")]:
        decision = dual.bongartz_test(bound)
        assert decision.value == expected_bongartz[key]["value"]
        assert decision.messages == ["Entering Bongartz' test for infinite type.\n"]

    field2 = qpa.Field.prime(2)
    wild_parent = qpa.PathAlgebra(
        field2,
        qpa.Quiver(["w"], [("w", "w", "x"), ("w", "w", "y")]),
    )
    wx = wild_parent.generator("x")
    wy = wild_parent.generator("y")
    wild = qpa.QuiverAlgebra.from_quotient_algebra(
        wild_parent.quotient_high_level([wx * wx, wx * wy, wy * wx, wy * wy])
    )
    assert wild.bongartz_test(20).value == expected_bongartz["wild"]["value"]

    expected_random = json.loads(
        (FIXTURES / "module_random.schema-v1.json").read_text(encoding="utf-8")
    )
    random_cases = [
        (qpa.Field.prime(2), 1, 4, "binary"),
        (qpa.Field.prime(3), 1, 1, "ternary"),
        (qpa.Field.rational(), 1, 1, "rational"),
        (qpa.Field.prime(2), 100, 1, "isolated"),
    ]
    for field, vertices, seed, key in random_cases:
        path = qpa.PathAlgebra(
            field, qpa.Quiver([f"v{index}" for index in range(vertices)], [])
        )
        module = qpa.QuiverAlgebra.from_path_algebra(path).random_module(1, 1, seed)
        assert (module is None) == expected_random[key]["is_fail"]
        assert (
            None if module is None else module.dimension_vector
        ) == expected_random[key]["dimension_vector"]


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2)])
def test_tensor_product_of_modules_follows_gap_qpa_oracle(field: qpa.Field) -> None:
    expected = json.loads(
        (FIXTURES / "functors_tensor_modules.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )["binary" if field.characteristic == 2 else "rational"]

    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    simples = path.simple_modules()
    zero_right = simples[1].dual()
    zero_tensor = simples[0].tensor_product(zero_right)
    assert zero_tensor.dimension == expected["easy_zero"]["space"]["dimension"]
    assert _rows(zero_tensor.basis_vectors) == expected["easy_zero"]["space"]["basis"]
    assert [
        _value(entry)
        for entry in zero_tensor.elementary_tensor(
            simples[0].basis_elements()[0], zero_right.basis_elements()[0]
        )
    ] == [0]

    loop = qpa.PathAlgebra(field, qpa.Quiver(["w"], [("w", "w", "x")]))
    x = loop.generator("x")
    cubic = loop.quotient_high_level([x * x * x])
    left = cubic.module_from_integer_matrices([2], [("x", [[0, 1], [0, 0]])])
    right = cubic.indecomposable_projectives()[0].dual()
    tensor = left.tensor_product(right)
    assert tensor.dimension == expected["cubic_quotient"]["space"]["dimension"]
    assert tensor.ambient_dimension == 3
    assert _rows(tensor.relation_basis) == [[1, 0, 0]]
    assert _rows(tensor.basis_lifts) == [[0, 1, 0], [0, 0, 1]]
    assert _rows(tensor.projection) == [[0, 0], [1, 0], [0, 1]]
    assert _rows(tensor.basis_vectors) == _oracle_rows(
        expected["cubic_quotient"]["space"]["basis"]
    )
    assert [
        [
            [_value(entry) for entry in tensor.elementary_tensor(left_element, right_element)]
            for right_element in right.basis_elements()
        ]
        for left_element in left.basis_elements()
    ] == [
        [[0, 0], [1, 0], [0, 1]],
        [[0, 0], [0, 0], [1, 0]],
    ]
