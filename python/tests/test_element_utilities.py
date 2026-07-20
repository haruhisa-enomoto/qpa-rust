from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _element_data(element: qpa.AlgebraElement) -> list[dict[str, Any]]:
    return [
        {
            "path": path.external_representation,
            "coefficient": (
                [int(coefficient.numerator), int(coefficient.denominator)]
                if coefficient.numerator is not None
                else coefficient.representative
            ),
        }
        for path, coefficient in element.terms
    ]


def _utility_parent(field: qpa.Field) -> qpa.PathAlgebra:
    return qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3", "v4"],
            [
                ("v1", "v2", "a"),
                ("v2", "v3", "b"),
                ("v1", "v3", "c"),
                ("v3", "v4", "d"),
                ("v2", "v4", "e"),
            ],
        ),
    )


def test_path_removal_and_support_follow_qpa() -> None:
    expected = json.loads(
        (FIXTURES / "pathalg_element_utilities.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    field = qpa.Field.rational()
    algebra = _utility_parent(field)
    v1, a, b, c, d = (algebra.generator(name) for name in ("v1", "a", "b", "c", "d"))
    ab, bd = a * b, b * d
    x = (
        ab.scale_left(field.element(2))
        + c.scale_left(field.element(3))
        + bd.scale_left(field.element(5))
        + d.scale_left(field.element(7))
    )
    m = (
        v1
        + a.scale_left(field.element(2))
        + ab.scale_left(field.element(3))
        + b.scale_left(field.element(5))
        + d.scale_left(field.element(7))
    )
    rational = expected["rational"]
    assert _element_data(x.path_removal([])) == rational["path_removal"]["none"]
    assert _element_data(x.path_removal([2])) == rational["path_removal"]["two"]
    assert _element_data(x.path_removal([4])) == rational["path_removal"]["four"]
    assert _element_data(x.path_removal([0, 5, 99, -10, 10**80])) == rational[
        "path_removal"
    ]["ignored"]
    assert algebra.left_support_of_quiver_algebra_element(m) == rational["support"]["left"]
    assert algebra.right_support_of_quiver_algebra_element(m) == rational["support"]["right"]
    assert algebra.support_of_quiver_algebra_element(m) == rational["support"]["both"]
    assert algebra.support_of_quiver_algebra_element(algebra.zero()) == rational["support"][
        "zero"
    ]

    quotient = algebra.quotient_high_level([ab - c, bd])
    qv1, qa, qb, qc, qd = (
        quotient.generator(name) for name in ("v1", "a", "b", "c", "d")
    )
    qm = (
        qv1
        + qa.scale_left(field.element(2))
        + qb.scale_left(field.element(3))
        + qd.scale_left(field.element(5))
    )
    qab = qa * qb
    assert quotient.support_of_quiver_algebra_element(qm) == rational["quotient"]["support"]
    assert quotient.support_of_quiver_algebra_element(qab + qc) == rational["quotient"][
        "mix_support"
    ]
    assert quotient.support_of_quiver_algebra_element(qab - qc) == rational["quotient"][
        "cancelled_support"
    ]
    assert quotient.left_support_of_quiver_algebra_element(qm) == rational["quotient"][
        "support"
    ][0]
    assert quotient.right_support_of_quiver_algebra_element(qm) == rational["quotient"][
        "support"
    ][1]

    field2 = qpa.Field.prime(2)
    algebra2 = _utility_parent(field2)
    a2, b2, d2 = (algebra2.generator(name) for name in ("a", "b", "d"))
    x2 = a2 + b2 * d2 + d2
    assert _element_data(x2.path_removal([2])) == expected["f2"]["removed_two"]
    assert _element_data(x2.path_removal([4])) == expected["f2"]["removed_four"]
    assert algebra2.support_of_quiver_algebra_element(x2) == expected["f2"]["support"]

    field5 = qpa.Field.prime(5)
    algebra5 = _utility_parent(field5)
    arrow5 = algebra5.generator("a")
    assert algebra5.support_of_quiver_algebra_element(arrow5) == [[1], [2]]


def test_element_utility_errors_follow_qpa() -> None:
    field = qpa.Field.rational()
    first = _utility_parent(field)
    second = _utility_parent(field)
    with pytest.raises(qpa.InvalidInputError, match="different path algebras"):
        first.left_support_of_quiver_algebra_element(second.one())

    vertex_parent = qpa.PathAlgebra(field, qpa.Quiver(["v"], []))
    with pytest.raises(qpa.InvalidInputError, match="vertex term"):
        vertex_parent.generator("v").path_removal([1])


def _foundational_parent(field: qpa.Field) -> qpa.PathAlgebra:
    return qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [
                ("v1", "v1", "a"),
                ("v1", "v2", "b"),
                ("v2", "v3", "c"),
                ("v3", "v1", "d"),
                ("v2", "v2", "e"),
            ],
        ),
    )


def test_add_nth_power_to_relations_preserves_qpa_mutation() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_foundational.schema-v1.json").read_text(encoding="utf-8")
    )
    algebra = _foundational_parent(qpa.Field.rational())
    relations = [algebra.generator("a") + algebra.generator("b")]
    returned = algebra.add_nth_power_to_relations(relations, 2)
    assert returned is relations
    assert [_element_data(relation) for relation in relations] == fixture["add_power"][
        "relations"
    ]

    foreign = _foundational_parent(qpa.Field.prime(2)).generator("a")
    with pytest.raises(qpa.InvalidInputError):
        algebra.add_nth_power_to_relations([foreign], 2)


def test_opposite_relations_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_foundational.schema-v1.json").read_text(encoding="utf-8")
    )
    algebra = _foundational_parent(qpa.Field.rational())
    a, b, c, d = (algebra.generator(name) for name in ("a", "b", "c", "d"))
    opposite = qpa.opposite_relations([a.pow(2) + b * c, c * d])
    assert [_element_data(relation) for relation in opposite] == fixture["opposite_parent"][
        "relations"
    ]
    assert qpa.opposite_relations([]) == []

    foreign = _foundational_parent(qpa.Field.prime(2)).generator("a")
    with pytest.raises(qpa.InvalidInputError):
        qpa.opposite_relations([a, foreign])


def test_tip_synonyms_follow_qpa() -> None:
    field = qpa.Field.rational()
    algebra = _utility_parent(field)
    element = algebra.generator("a") + (algebra.generator("a") * algebra.generator("b"))
    assert element.tip == element.leading_term
    assert element.tip_coefficient == element.leading_coefficient
    assert element.tip_monomial == element.leading_monomial

    vector = algebra.right_projective_module([1]).generators[0]
    assert vector.tip.components == vector.leading_term.components
    assert vector.tip_coefficient == vector.leading_coefficient

    loop = qpa.Quiver(["v"], [("v", "v", "x")])
    quotient = qpa.QuotientAlgebra.truncated(field, loop, 2)
    quotient_element = quotient.one() + quotient.generator("x")
    assert quotient_element.tip == quotient_element.leading_term
    assert quotient_element.tip_coefficient == quotient_element.leading_coefficient
    assert quotient_element.tip_monomial == quotient_element.leading_monomial

    quotient_vector = quotient.right_projective_module([1]).generators[0]
    assert quotient_vector.tip.components == quotient_vector.leading_term.components
    assert quotient_vector.tip_coefficient == quotient_vector.leading_coefficient
