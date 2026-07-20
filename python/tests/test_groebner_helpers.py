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


def _external_element_data(element: qpa.AlgebraElement) -> list[dict[str, Any]]:
    return [
        {
            "path": path.external_representation,
            "coefficient": _coefficient_data(coefficient),
        }
        for path, coefficient in element.terms
    ]


def _named_element_data(element: qpa.AlgebraElement) -> list[dict[str, Any]]:
    return [
        {"path": repr(path), "coefficient": _coefficient_data(coefficient)}
        for path, coefficient in element.terms
    ]


def _cohen_data(
    representations: list[tuple[list[list[int]], list[qpa.Scalar]]],
) -> list[list[dict[str, Any]]]:
    return [
        [
            {"monomial": monomial, "coefficient": _coefficient_data(coefficient)}
            for monomial, coefficient in zip(monomials, coefficients, strict=True)
        ]
        for monomials, coefficients in representations
    ]


def _loop_algebra(
    field: qpa.Field, arrow_names: tuple[str, ...] = ("x",)
) -> qpa.PathAlgebra:
    return qpa.PathAlgebra(
        field,
        qpa.Quiver(["v1"], [("v1", "v1", name) for name in arrow_names]),
    )


@pytest.mark.parametrize(
    ("field", "key", "middle_coefficient"),
    [
        (qpa.Field.rational(), "rational", -1),
        (qpa.Field.prime(5), "f5", -2),
        (qpa.Field.prime(2), "f2", 1),
    ],
)
def test_high_level_groebner_basis_and_division_follow_qpa(
    field: qpa.Field, key: str, middle_coefficient: int
) -> None:
    fixture = json.loads(
        (FIXTURES / "groebner_high_level.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )[key]
    algebra = _loop_algebra(field, ("x", "y"))
    x, y = (algebra.generator(name) for name in ("x", "y"))
    relations = [
        x.pow(2),
        x * y + (y * x).scale_left(field.element(middle_coefficient)),
        y.pow(2),
    ]
    basis = algebra.groebner_basis(relations).relations

    assert [_named_element_data(value) for value in basis] == fixture["basis"]
    assert _named_element_data(algebra.remainder_of_division(y * x, basis)) == fixture[
        "remainder_yx"
    ]

    if key == "rational":
        assert _named_element_data(
            algebra.remainder_of_division(y * x * y + x, basis)
        ) == fixture["remainder_yxy_plus_x"]
        assert _named_element_data(
            algebra.remainder_of_division(x.pow(2) + x, [algebra.one()])
        ) == fixture["remainder_by_vertex"]
        assert [
            _named_element_data(value)
            for value in algebra.reduced_list_qpa(
                [x.pow(2), x.pow(2) + x, algebra.zero()]
            )
        ] == fixture["reduced_list"]
        assert [
            _named_element_data(value)
            for value in algebra.tip_reduced_list(
                [algebra.one(), x, x.pow(2), x.pow(2), algebra.zero()]
            )
        ] == fixture["tip_reduced_list"]
        assert [arrow.name for arrow in (x + y * x).tip_walk()] == fixture["tip_walk"]
        assert algebra.groebner_basis([]).relations == fixture["empty_basis"]


def test_uniformization_and_generic_leftmost_occurrence_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "groebner_high_level.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    algebra = qpa.PathAlgebra(
        qpa.Field.rational(),
        qpa.Quiver(
            ["v1", "v2"],
            [
                ("v1", "v1", "a"),
                ("v1", "v2", "b"),
                ("v2", "v1", "c"),
            ],
        ),
    )
    a, b, c = (algebra.generator(name) for name in ("a", "b", "c"))

    assert [
        _named_element_data(value) for value in algebra.make_uniform(a + b + c)
    ] == fixture["uniform_components"]
    assert [
        _named_element_data(value)
        for value in algebra.make_uniform_list([a + b + c])
    ] == fixture["uniform_components"]
    assert [
        _named_element_data(value)
        for value in algebra.groebner_basis([a + b + c]).relations
    ] == fixture["split_basis"]

    expected = fixture["leftmost_occurrence"]
    # Python deliberately uses normal zero-based, half-open slice endpoints.
    assert qpa.leftmost_occurrence([1, 2, 1, 2], [1, 2]) == (0, 2)
    assert qpa.leftmost_occurrence([1, 2], []) == (0, 0)
    assert qpa.leftmost_occurrence([1, 2], [2, 1]) is expected["absent"]
    assert qpa.leftmost_occurrence(["a", "b", "a"], ["b", "a"]) == (1, 3)
    assert qpa.leftmost_occurrence(algebra.quiver.arrow_objects, [algebra.quiver.arrow("b")]) == (
        1,
        2,
    )


@pytest.mark.parametrize(
    ("field", "key"),
    [
        (qpa.Field.rational(), "rational"),
        (qpa.Field.prime(5), "f5"),
        (qpa.Field.prime(2), "f2"),
    ],
)
def test_unchecked_gbnp_bridge_follows_qpa(
    field: qpa.Field, key: str
) -> None:
    fixture = json.loads(
        (FIXTURES / "qpa_gbnp.schema-v1.json").read_text(encoding="utf-8")
    )[key]
    algebra = _loop_algebra(field, ("a", "b"))
    vertex, a, b = (algebra.generator(name) for name in ("v1", "a", "b"))
    relations = [a * a * b - vertex, a * b * b - vertex]

    assert [
        _external_element_data(value)
        for value in algebra.gbnp_groebner_basis_nc(relations)
    ] == fixture["nc_basis"]
    if key == "rational":
        assert not algebra.qpa_in_arrow_ideal([vertex])
        assert algebra.qpa_in_arrow_ideal([a, b])
        ideal = algebra.ideal([a.pow(2), a * b])
        assert [
            _external_element_data(value)
            for value in ideal.gbnp_groebner_basis()
        ] == fixture["checked_basis"]


@pytest.mark.parametrize(
    ("field", "key"),
    [
        (qpa.Field.rational(), "rational"),
        (qpa.Field.prime(5), "f5"),
        (qpa.Field.prime(2), "f2"),
    ],
)
def test_qpa_cohen_conversions_follow_qpa(field: qpa.Field, key: str) -> None:
    fixture = json.loads(
        (FIXTURES / "qpa_gbnp.schema-v1.json").read_text(encoding="utf-8")
    )[key]
    algebra = _loop_algebra(field, ("a", "b"))
    vertex, a, b = (algebra.generator(name) for name in ("v1", "a", "b"))
    relations = [a * a * b - vertex, a * b * b - vertex]

    free = algebra.qpa_path_to_cohen_free(relations)
    assert _cohen_data(free) == fixture["free"]
    assert algebra.qpa_cohen_to_path(free) == relations
    if key == "rational":
        raw = algebra.qpa_path_to_cohen(relations)
        assert _cohen_data(raw) == fixture["raw"]
        assert algebra.qpa_cohen_to_path(raw) == relations


def test_qpa_path_algebra_structure_relations_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "qpa_gbnp.schema-v1.json").read_text(encoding="utf-8")
    )["multi_vertex"]
    algebra = qpa.PathAlgebra(
        qpa.Field.rational(),
        qpa.Quiver(["v1", "v2"], [("v1", "v2", "c")]),
    )
    structure = algebra.qpa_relations_for_path_algebra()

    assert _cohen_data(structure) == fixture["structure_relations"]
    assert [
        _external_element_data(value) for value in algebra.qpa_cohen_to_path(structure)
    ] == fixture["structure_as_paths"]
    with pytest.raises(qpa.InvalidInputError, match="one-vertex"):
        algebra.qpa_path_to_cohen_free([algebra.generator("c")])
    with pytest.raises(qpa.InvalidInputError, match="only 0 coefficients"):
        algebra.qpa_cohen_to_path([([[3]], [])])
    with pytest.raises(qpa.InvalidInputError, match="belongs to a prime field"):
        algebra.qpa_cohen_to_path([([[3]], [qpa.Field.prime(2).element(1)])])


def _inverse_fixture_element(
    field: qpa.Field,
) -> tuple[qpa.PathAlgebra, qpa.AlgebraElement, qpa.AlgebraElement]:
    algebra = _loop_algebra(field)
    vertex, x = (algebra.generator(name) for name in ("v1", "x"))
    x2 = x.pow(2)
    x3 = x2 * x
    coefficients = (1, 1, 1) if field.characteristic == 2 else (2, 3, 4)
    mixed = (
        vertex
        + x.scale_left(field.element(coefficients[0]))
        + x2.scale_left(field.element(coefficients[1]))
        + x3.scale_left(field.element(coefficients[2]))
    )
    return algebra, x, mixed


@pytest.mark.parametrize(
    ("field", "key"),
    [
        (qpa.Field.rational(), "rational"),
        (qpa.Field.prime(5), "f5"),
        (qpa.Field.prime(2), "f2"),
    ],
)
def test_inverse_length_inverse_lex_element_operations_follow_qpa(
    field: qpa.Field, key: str
) -> None:
    fixture = json.loads(
        (FIXTURES / "invlenlex.schema-v1.json").read_text(encoding="utf-8")
    )[key]
    algebra, x, mixed = _inverse_fixture_element(field)
    x3 = x.pow(3)

    assert _external_element_data(mixed.tip_inv_len_inv_lex()) == fixture["tip"]
    assert _external_element_data(mixed.truncate_element(2)) == fixture["truncate_2"]
    assert [
        _external_element_data(value)
        for value in algebra.groebner_basis_inv_len_inv_lex([x3])
    ] == fixture["basis"]

    if key == "rational":
        assert _coefficient_data(mixed.tip_coefficient_inv_len_inv_lex) == fixture[
            "tip_coefficient"
        ]
        tip_monomial = mixed.tip_monomial_inv_len_inv_lex
        assert tip_monomial is not None
        assert tip_monomial.external_representation == fixture["tip_monomial"]
        assert [arrow.name for arrow in mixed.tip_walk_inv_len_inv_lex()] == fixture[
            "tip_walk"
        ]
        assert _external_element_data(mixed.truncate_element(1)) == fixture["truncate_1"]
        positive = x.scale_left(field.element(2)) + x.pow(2).scale_left(field.element(3))
        assert _external_element_data(
            positive.biggest_homogeneous_part_inv_len_inv_lex()
        ) == fixture["biggest_homogeneous"]
        assert [
            _external_element_data(value)
            for value in algebra.biggest_homogeneous_parts_inv_len_inv_lex(
                [mixed, positive]
            )
        ] == fixture["biggest_homogeneous_list"]
        assert _external_element_data(
            algebra.remainder_of_division_inv_len_inv_lex(x3 + x, [x.pow(2)])
        ) == fixture["remainder"]
        assert [
            _external_element_data(value)
            for value in algebra.reduced_list_inv_len_inv_lex([x.pow(2), x3])
        ] == fixture["reduced_list"]
        assert [
            _external_element_data(value)
            for value in algebra.tip_reduced_list_inv_len_inv_lex([x.pow(2), x3])
        ] == fixture["tip_reduced_list"]


def test_inverse_length_inverse_lex_boundaries_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "invlenlex.schema-v1.json").read_text(encoding="utf-8")
    )["rational"]
    algebra, _, mixed = _inverse_fixture_element(qpa.Field.rational())
    zero = algebra.zero()

    assert _external_element_data(zero.tip_inv_len_inv_lex()) == fixture["zero"]["tip"]
    assert _coefficient_data(zero.tip_coefficient_inv_len_inv_lex) == fixture["zero"][
        "tip_coefficient"
    ]
    assert zero.tip_monomial_inv_len_inv_lex is None
    with pytest.raises(qpa.InvalidInputError, match="undefined for zero"):
        zero.tip_walk_inv_len_inv_lex()
    with pytest.raises(qpa.InvalidInputError):
        mixed.truncate_element(0)

    assert qpa.leftmost_occurrence_inv_len_inv_lex([1, 2, 1, 2], [1, 2]) == (
        0,
        2,
    )
    assert qpa.leftmost_occurrence_inv_len_inv_lex([1, 2], []) == (0, 0)
    assert qpa.leftmost_occurrence_inv_len_inv_lex([1, 2], [2, 1]) is None

    split = qpa.PathAlgebra(
        qpa.Field.rational(),
        qpa.Quiver(
            ["v1", "v2"],
            [("v1", "v1", "a"), ("v2", "v2", "b")],
        ),
    )
    a, b = (split.generator(name) for name in ("a", "b"))
    assert [
        _external_element_data(value)
        for value in split.groebner_basis_inv_len_inv_lex([a.pow(2) + b.pow(2)])
    ] == fixture["nonuniform_basis"]
