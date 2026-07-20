from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _coefficient_data(value: qpa.Scalar) -> list[int] | int:
    if value.numerator is not None:
        return [int(value.numerator), int(value.denominator)]
    representative = value.representative
    assert representative is not None
    return representative


def _named_element_data(element: qpa.AlgebraElement) -> list[dict[str, Any]]:
    return [
        {"path": repr(path), "coefficient": _coefficient_data(coefficient)}
        for path, coefficient in element.terms
    ]


def test_null_list_and_positions_nonzero_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "util.schema-v1.json").read_text(encoding="utf-8")
    )
    rational = qpa.Field.rational()
    rational_values = [
        rational.element(0),
        rational.element(2),
        rational.element(0),
        rational.element(-1),
        rational.element(0),
    ]
    assert qpa.positions_nonzero(rational_values) == fixture["rational"]["positions"]
    rational_null = rational.null_list(4)
    assert [value.is_zero for value in rational_null] == [True] * 4

    binary = qpa.Field.prime(2)
    binary_values = [
        binary.element(0),
        binary.element(1),
        binary.element(0),
        binary.element(1),
        binary.element(0),
    ]
    assert qpa.positions_nonzero(binary_values) == fixture["binary"]["positions"]
    assert [_coefficient_data(value) for value in binary.null_list(4)] == [0, 0, 0, 0]
    assert qpa.positions_nonzero([[[0]], [[1]], [[0]], [[-2]]]) == fixture[
        "generic_matrices"
    ]["positions"]
    with pytest.raises(qpa.InvalidInputError, match="positive"):
        rational.null_list(0)


def test_paths_of_length_two_and_homogeneity_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_koszul.schema-v1.json").read_text(encoding="utf-8")
    )
    loop_quiver = qpa.Quiver(
        ["v1"], [("v1", "v1", "x"), ("v1", "v1", "y")]
    )
    assert [repr(path) for path in loop_quiver.paths_of_length_two()] == fixture[
        "paths_loop"
    ]
    multi = qpa.Quiver(
        ["v1", "v2"],
        [("v1", "v2", "a"), ("v2", "v1", "b"), ("v1", "v1", "c")],
    )
    assert [repr(path) for path in multi.paths_of_length_two()] == fixture["paths_multi"]

    ideal_fixture = json.loads(
        (FIXTURES / "pathalg_ideal.schema-v1.json").read_text(encoding="utf-8")
    )["homogeneous"]
    algebra = qpa.PathAlgebra(qpa.Field.rational(), loop_quiver)
    x = algebra.generator("x")
    assert algebra.is_homogeneous_list_of_elements([], 2) is ideal_fixture[
        "empty_degree_2"
    ]
    assert algebra.is_homogeneous_list_of_elements(
        [x.pow(2), x.pow(2).scale_left(algebra.field.element(2)), algebra.zero()], 2
    ) is ideal_fixture["degree_2"]
    assert algebra.is_homogeneous_list_of_elements([x.pow(2) + x.pow(3)], 2) is ideal_fixture[
        "mixed_degree_2"
    ]
    with pytest.raises(qpa.InvalidInputError, match="positive"):
        algebra.is_homogeneous_list_of_elements([x], 0)


@pytest.mark.parametrize(
    ("field", "fixture_key"),
    [
        (qpa.Field.rational(), "rational_perpendicular"),
        (qpa.Field.prime(2), "f2_perpendicular"),
        (qpa.Field.prime(5), "f5_perpendicular"),
    ],
)
def test_quadratic_perpendicular_follows_qpa(
    field: qpa.Field, fixture_key: str
) -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_koszul.schema-v1.json").read_text(encoding="utf-8")
    )
    algebra = qpa.PathAlgebra(
        field,
        qpa.Quiver(["v1"], [("v1", "v1", "x"), ("v1", "v1", "y")]),
    )
    x, y = (algebra.generator(name) for name in ("x", "y"))
    second_relation = x * y
    if field.characteristic != 2:
        second_relation = second_relation - (y * x).scale_left(field.element(2))
    relations = [x.pow(2) + y.pow(2), second_relation]

    assert algebra.is_quadratic_ideal([]) is fixture["quadratic"]["empty"]
    assert algebra.is_quadratic_ideal(relations) is fixture["quadratic"]["relations"]
    assert algebra.is_quadratic_ideal([algebra.zero()]) is fixture["quadratic"]["zero"]
    assert algebra.is_quadratic_ideal([x]) is fixture["quadratic"]["linear"]

    result = algebra.quadratic_perpendicular_of_path_algebra_ideal(relations)
    assert result is not None
    opposite, basis = result
    assert opposite.quiver.vertices == fixture["opposite_quiver"]["vertices"]
    assert [name for _, _, name in opposite.quiver.arrows] == fixture["opposite_quiver"][
        "arrows"
    ]
    assert [_named_element_data(value) for value in basis] == fixture[fixture_key]


def test_quadratic_perpendicular_boundary_results_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_koszul.schema-v1.json").read_text(encoding="utf-8")
    )
    algebra = qpa.PathAlgebra(
        qpa.Field.rational(),
        qpa.Quiver(["v1"], [("v1", "v1", "x"), ("v1", "v1", "y")]),
    )
    x, y = (algebra.generator(name) for name in ("x", "y"))
    zero_result = algebra.quadratic_perpendicular_of_path_algebra_ideal([algebra.zero()])
    assert zero_result is not None
    assert [_named_element_data(value) for value in zero_result[1]] == fixture[
        "zero_perpendicular"
    ]
    full_result = algebra.quadratic_perpendicular_of_path_algebra_ideal(
        [x.pow(2), x * y, y * x, y.pow(2)]
    )
    assert full_result is not None
    assert [_named_element_data(value) for value in full_result[1]] == fixture[
        "full_perpendicular"
    ]
    assert algebra.quadratic_perpendicular_of_path_algebra_ideal([x]) is None
    with pytest.raises(qpa.InvalidInputError, match="first generator"):
        algebra.quadratic_perpendicular_of_path_algebra_ideal([])

    semisimple = qpa.PathAlgebra(qpa.Field.rational(), qpa.Quiver(["v1"], []))
    with pytest.raises(qpa.InvalidInputError, match="NullspaceMat"):
        semisimple.quadratic_perpendicular_of_path_algebra_ideal([semisimple.zero()])
