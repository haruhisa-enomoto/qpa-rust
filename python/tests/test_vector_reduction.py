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


def _row_data(
    rows: list[list[qpa.AlgebraElement | qpa.QuotientElement]],
) -> list[list[list[dict[str, Any]]]]:
    return [[_element_data(element) for element in row] for row in rows]


def _input(field: qpa.Field) -> qpa.PathAlgebra:
    return qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )


def test_path_vector_reduction_follows_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_vector_reduction.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    field = qpa.Field.rational()
    algebra = _input(field)
    v1, v2, a = (algebra.generator(name) for name in ("v1", "v2", "a"))
    zero = algebra.zero()

    position, path, coefficient = algebra.tip_monomial_and_coefficient_of_vector(
        [v1 + a, a.scale_left(field.element(2))]
    )
    assert {
        "position": position,
        "path": path.external_representation,
        "coefficient": _coefficient_data(coefficient),
    } == fixture["rational"]["tip"]
    zero_position, zero_path, zero_coefficient = (
        algebra.tip_monomial_and_coefficient_of_vector([zero, zero])
    )
    assert {
        "position": zero_position,
        "path": zero_path.external_representation,
        "coefficient": _coefficient_data(zero_coefficient),
    } == fixture["rational"]["zero_tip"]

    rows = [[a, zero], [a.scale_left(field.element(2)), v2]]
    assert algebra.tip_reduce_vectors(rows) is rows
    assert _row_data(rows) == fixture["rational"]["reduced"]

    basis = [[a, zero], [zero, v2]]
    target = [a.scale_left(field.element(3)), v2.scale_left(field.element(5))]
    coefficients = algebra.coefficients_of_vectors(basis, target)
    assert coefficients.is_non_span_witness is False
    assert [_coefficient_data(value) for value in coefficients.coefficients or []] == fixture[
        "rational"
    ]["coefficients"]
    assert coefficients.non_span_basis is None
    assert coefficients.non_span_vector is None

    not_span = [v1, zero]
    witness = algebra.coefficients_of_vectors(basis, not_span)
    assert witness.is_non_span_witness is fixture["rational"]["failure_witness"]
    assert witness.coefficients is None
    assert witness.non_span_basis == basis
    assert witness.non_span_vector == not_span


def test_quotient_and_binary_vector_reduction_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_vector_reduction.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    field = qpa.Field.rational()
    parent = _input(field)
    quotient = parent.quotient([])
    qv2, qa = (quotient.generator(name) for name in ("v2", "a"))
    qzero = quotient.zero()
    quotient_result = quotient.coefficients_of_vectors(
        [[qa, qzero], [qzero, qv2]],
        [qa.scale_left(field.element(3)), qv2.scale_left(field.element(5))],
    )
    assert [_coefficient_data(value) for value in quotient_result.coefficients or []] == fixture[
        "rational"
    ]["quotient_coefficients"]
    quotient_rows = [[qa, qzero], [qa.scale_left(field.element(2)), qv2]]
    assert quotient.tip_reduce_vectors(quotient_rows) is quotient_rows

    binary = _input(qpa.Field.prime(2))
    b_v2, b_a = (binary.generator(name) for name in ("v2", "a"))
    b_zero = binary.zero()
    binary_rows = [[b_a, b_zero], [b_a, b_v2]]
    assert binary.tip_reduce_vectors(binary_rows) is binary_rows
    assert _row_data(binary_rows) == fixture["f2"]["reduced"]
    binary_result = binary.coefficients_of_vectors(
        [[b_a, b_zero], [b_zero, b_v2]], [b_a, b_v2]
    )
    assert [_coefficient_data(value) for value in binary_result.coefficients or []] == fixture[
        "f2"
    ]["coefficients"]


def test_vector_reduction_errors_and_empty_fallthrough_follow_qpa() -> None:
    algebra = _input(qpa.Field.rational())
    with pytest.raises(qpa.InvalidInputError, match="requires a nonempty vector"):
        algebra.tip_monomial_and_coefficient_of_vector([])
    empty: list[list[qpa.AlgebraElement]] = []
    assert algebra.tip_reduce_vectors(empty) is None

    a = algebra.generator("a")
    with pytest.raises(qpa.InvalidInputError, match="length 1, expected 2"):
        algebra.tip_reduce_vectors([[a, algebra.zero()], [a]])

    foreign = _input(qpa.Field.rational())
    with pytest.raises(qpa.InvalidInputError, match="different algebra"):
        algebra.coefficients_of_vectors([[a]], [foreign.generator("a")])
