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


def _path_element_data(element: qpa.AlgebraElement) -> list[dict[str, Any]]:
    return [
        {
            "path": path.external_representation,
            "coefficient": _coefficient_data(coefficient),
        }
        for path, coefficient in element.terms
    ]


def _quotient_element_data(element: qpa.QuotientElement) -> list[dict[str, Any]]:
    return _path_element_data(element.representative)


def test_path_and_quotient_centres_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_centre.schema-v1.json").read_text(encoding="utf-8")
    )
    rationals = qpa.Field.rational()
    cycle_quiver = qpa.Quiver(
        ["v1", "v2"], [("v1", "v2", "a"), ("v2", "v1", "b")]
    )
    parent = qpa.PathAlgebra(rationals, cycle_quiver)
    parent_centre = parent.centre()
    assert parent_centre.ambient.is_identical(parent)
    assert parent_centre.dimension is None
    assert [_path_element_data(element) for element in parent_centre.generators] == fixture[
        "parent_cycle"
    ]["generators"]
    assert parent.one() in parent_centre
    assert parent.generator("a") not in parent_centre

    acyclic = qpa.PathAlgebra(
        rationals,
        qpa.Quiver(
            ["v1", "v2", "v3"], [("v1", "v2", "u"), ("v2", "v3", "v")]
        ),
    )
    acyclic_centre = acyclic.centre()
    assert acyclic_centre.dimension == "1"
    assert [_path_element_data(element) for element in acyclic_centre.generators] == fixture[
        "parent_acyclic"
    ]["generators"]

    discrete = qpa.PathAlgebra(rationals, qpa.Quiver(["v1", "v2"], []))
    discrete_centre = discrete.centre()
    assert discrete_centre.dimension == "2"
    assert [_path_element_data(element) for element in discrete_centre.generators] == fixture[
        "parent_discrete"
    ]["generators"]

    quotient = qpa.QuotientAlgebra.truncated(rationals, cycle_quiver, 4)
    quotient_centre = quotient.centre()
    a = quotient.generator("a")
    b = quotient.generator("b")
    central = a * b + b * a
    assert quotient_centre.ambient.is_identical(quotient)
    assert quotient.dimension() == str(fixture["quotient"]["ambient_dimension"])
    assert quotient_centre.dimension == fixture["quotient"]["dimension"]
    assert [_quotient_element_data(element) for element in quotient_centre.basis] == fixture[
        "quotient"
    ]["basis"]
    assert (central in quotient_centre) is fixture["quotient"]["central_member"]
    assert (a in quotient_centre) is fixture["quotient"]["arrow_member"]

    binary = qpa.QuotientAlgebra.truncated(qpa.Field.prime(2), cycle_quiver, 4)
    binary_centre = binary.centre()
    assert binary_centre.dimension == fixture["f2"]["dimension"]
    assert [_quotient_element_data(element) for element in binary_centre.basis] == fixture[
        "f2"
    ]["basis"]


def _coordinate_data(group: qpa.QuiverAlgebraUnitGroup) -> list[list[int]]:
    return [
        [
            coefficient.representative
            if coefficient.representative is not None
            else int(coefficient.numerator)
            for coefficient in row
        ]
        for row in group.coordinates
    ]


def _assert_group(
    group: qpa.QuiverAlgebraUnitGroup, expected: dict[str, Any]
) -> None:
    assert len(group) == group.order == expected["size"]
    assert _coordinate_data(group) == expected["elements"]
    assert group.multiplication_table == expected["table"]
    assert group.identity_index == 0
    for left, row in enumerate(expected["table"]):
        for right, one_based_product in enumerate(row):
            assert group.product_index(left, right) == one_based_product - 1
    assert all(group.inverse_index(index) is not None for index in range(group.order))
    assert group.product_index(group.order, 0) is None
    assert group.inverse_index(group.order) is None


def test_finite_unit_groups_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_units.schema-v1.json").read_text(encoding="utf-8")
    )
    binary_a2 = qpa.PathAlgebra(
        qpa.Field.prime(2), qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    assert binary_a2.has_elements_of_units is fixture["binary_a2"]["had_elements"]
    with pytest.raises(qpa.InvalidInputError, match="ElementsOfUnits is not bound"):
        binary_a2.elements_of_units()
    binary_group = binary_a2.units()
    _assert_group(binary_group, fixture["binary_a2"])
    assert binary_a2.has_elements_of_units is fixture["binary_a2"]["has_elements"]
    assert binary_a2.elements_of_units() == binary_group.elements
    binary_union = qpa.QuiverAlgebra.from_path_algebra(binary_a2)
    assert binary_union.has_elements_of_units
    assert binary_union.elements_of_units() == binary_group.elements

    loop_quiver = qpa.Quiver(["v1"], [("v1", "v1", "x")])
    dual_parent = qpa.PathAlgebra(qpa.Field.prime(3), loop_quiver)
    dual = dual_parent.quotient([dual_parent.generator("x").pow(2)])
    dual_group = dual.units()
    _assert_group(dual_group, fixture["ternary_dual_numbers"])
    assert dual.elements_of_units() == dual_group.elements

    cubic_parent = qpa.PathAlgebra(qpa.Field.prime(2), loop_quiver)
    cubic = cubic_parent.quotient([cubic_parent.generator("x").pow(3)])
    _assert_group(cubic.units(), fixture["binary_cubic_numbers"])

    semisimple = qpa.PathAlgebra(qpa.Field.prime(3), qpa.Quiver(["u", "v"], []))
    _assert_group(semisimple.units(), fixture["ternary_semisimple"])

    quinary = qpa.PathAlgebra(qpa.Field.prime(5), qpa.Quiver(["v1"], []))
    _assert_group(quinary.units(), fixture["quinary_a1"])


def test_unit_group_errors_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_units.schema-v1.json").read_text(encoding="utf-8")
    )["errors"]
    infinite = qpa.PathAlgebra(
        qpa.Field.prime(2), qpa.Quiver(["v1"], [("v1", "v1", "z")])
    )
    with pytest.raises(qpa.InvalidInputError, match="not finite dimensional"):
        infinite.units()
    assert fixture["infinite_success"] is False

    rational = qpa.PathAlgebra(qpa.Field.rational(), qpa.Quiver(["v1"], []))
    with pytest.raises(qpa.InvalidInputError, match="not over a finite field"):
        rational.units()
    assert fixture["rational_success"] is False
