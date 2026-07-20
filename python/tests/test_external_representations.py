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


def _element_data(
    element: qpa.AlgebraElement | qpa.QuotientElement,
) -> list[dict[str, Any]]:
    representative = (
        element.representative if isinstance(element, qpa.QuotientElement) else element
    )
    return [
        {
            "path": path.external_representation,
            "coefficient": _coefficient_data(coefficient),
        }
        for path, coefficient in representative.terms
    ]


def _named_element_data(element: qpa.QuotientElement) -> list[dict[str, Any]]:
    return [
        {"path": repr(path), "coefficient": _coefficient_data(coefficient)}
        for path, coefficient in element.representative.terms
    ]


def _conversion_algebra(
    field: qpa.Field,
) -> tuple[
    qpa.PathAlgebra,
    qpa.AlgebraElement,
    qpa.AlgebraElement,
    qpa.AlgebraElement,
]:
    algebra = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [("v1", "v2", "a"), ("v2", "v3", "b"), ("v1", "v3", "c")],
        ),
    )
    a, b, c = (algebra.generator(name) for name in ("a", "b", "c"))
    return algebra, a, b, c


def test_external_representation_round_trips_and_normal_forms_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_conversion.schema-v1.json").read_text(encoding="utf-8")
    )
    field = qpa.Field.rational()
    algebra, a, b, c = _conversion_algebra(field)
    relation = a * b - c
    quotient = algebra.quotient_high_level([relation])
    raw = a.scale_left(field.element(2)) + (a * b).scale_left(field.element(3))
    descriptor = raw.external_representation

    assert algebra.is_path_ring
    assert not quotient.is_path_ring
    assert quotient.is_full_fp_path_algebra
    assert algebra.element(a.terms[0][0]) == a
    assert algebra.element_from_external_representation(descriptor) == raw
    assert quotient.relations == quotient.relators
    assert [_element_data(value) for value in quotient.relations] == fixture["rational"][
        "quotient_relations"
    ]

    rebuilt = quotient.element_from_external_representation(descriptor)
    assert _element_data(rebuilt) == fixture["rational"]["rebuilt"]
    assert rebuilt.is_normal_form is fixture["rational"]["rebuilt_is_normal_form"]
    assert _element_data(quotient.element_with_normality(raw, False)) == fixture[
        "rational"
    ]["rebuilt"]
    assert _element_data(quotient.normal_form(raw)) == fixture["rational"]["rebuilt"]

    # The Boolean constructor really does trust `true`, just as QPA does.
    trusted = quotient.element_with_normality(raw, True)
    assert trusted.is_normal_form
    assert _element_data(trusted) == _element_data(raw)

    default_basis = algebra.groebner_basis_function([relation])
    explicit_basis = algebra.gbnp_groebner_basis([relation])
    assert default_basis.relations == explicit_basis.relations


def test_binary_external_representation_round_trip_follows_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_conversion.schema-v1.json").read_text(encoding="utf-8")
    )
    field = qpa.Field.prime(2)
    algebra, a, b, c = _conversion_algebra(field)
    quotient = algebra.quotient_high_level([a * b - c])
    raw = a + a * b

    rebuilt = quotient.element_from_external_representation(raw.external_representation)
    assert _element_data(rebuilt) == fixture["f2"]["rebuilt"]
    assert rebuilt.is_normal_form is fixture["f2"]["rebuilt_is_normal_form"]

    rational, rational_a, _, _ = _conversion_algebra(qpa.Field.rational())
    with pytest.raises(qpa.InvalidInputError, match="prime field"):
        rational.element_from_external_representation(raw.external_representation)
    assert rational.element_from_external_representation(
        rational_a.external_representation
    ) == rational_a


def test_lazy_quotient_preserves_qpa_normal_form_function_boundary() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_conversion.schema-v1.json").read_text(encoding="utf-8")
    )["lazy_default"]
    field = qpa.Field.rational()
    algebra = qpa.PathAlgebra(
        field, qpa.Quiver(["v1"], [("v1", "v1", "x")])
    )
    x = algebra.generator("x")
    ideal = algebra.ideal([x * x])
    quotient = ideal.quotient_without_groebner_basis()
    raw = x + x * x

    assert quotient.has_normal_form_function is fixture["before_normal_function"]
    lazy = quotient.element_with_normality(raw, False)
    assert _element_data(lazy) == fixture["before_element"]
    assert lazy.is_normal_form is fixture["before_element_is_normal"]
    rebuilt = quotient.element_from_external_representation(raw.external_representation)
    assert _element_data(rebuilt) == fixture["before_element"]
    assert not rebuilt.is_normal_form
    with pytest.raises(qpa.InvalidInputError, match="no normal-form function"):
        quotient.normal_form(raw)
    with pytest.raises(qpa.InvalidInputError, match="no attached Groebner basis"):
        rebuilt.normalized()

    assert quotient.dimension() == str(fixture["dimension"])
    assert quotient.has_normal_form_function is fixture["after_normal_function"]
    after = quotient.element_with_normality(raw, False)
    assert _element_data(after) == fixture["after_element"]
    assert after.is_normal_form is fixture["after_element_is_normal"]


def test_quotient_element_one_and_zero_operations_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_surface.schema-v1.json").read_text(encoding="utf-8")
    )
    field = qpa.Field.prime(2)
    algebra = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2"], [("v1", "v2", "a"), ("v2", "v1", "b")]
        ),
    )
    quotient = algebra.quotient_high_level([])
    element = quotient.generator("a") + quotient.generator("b")

    assert element.one_same_parent() == quotient.one()
    assert _named_element_data(element.one_same_parent()) == fixture["quotient"][
        "arithmetic"
    ]["one"]
    assert element.zero_same_parent() == quotient.zero()
    assert _named_element_data(element.zero_same_parent()) == fixture["quotient"][
        "arithmetic"
    ]["zero"]
    assert element.pow(0) == element.one_same_parent()
