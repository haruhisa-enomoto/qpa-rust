from __future__ import annotations

import itertools
import json
from pathlib import Path

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _relation_paths(algebra: qpa.QuiverAlgebra) -> list[list[dict[str, object]]]:
    quotient = algebra.quotient_algebra
    assert quotient is not None
    return [
        [
            {
                "path": path.external_representation,
                "coefficient": (
                    [int(coefficient.numerator), int(coefficient.denominator)]
                    if coefficient.numerator is not None
                    else coefficient.representative
                ),
            }
            for path, coefficient in relation.terms
        ]
        for relation in quotient.relations
    ]


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_associated_algebras_follow_qpa(field: qpa.Field) -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_associated.schema-v1.json").read_text(encoding="utf-8")
    )
    expected = fixture["rational" if field.characteristic is None else "f2"]

    finite_parent = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    finite = qpa.QuiverAlgebra.from_path_algebra(finite_parent)
    monomial_finite = finite.associated_monomial_algebra()
    graded_finite = finite.associated_graded_algebra()
    assert finite_parent.is_monomial is expected["finite_path_is_monomial"]
    assert monomial_finite.is_identical(finite) is expected["finite_path_identity"]
    assert graded_finite.is_identical(finite) is expected["finite_path_graded_identity"]

    loop_quiver = qpa.Quiver(["v"], [("v", "v", "x"), ("v", "v", "y")])
    parent = qpa.PathAlgebra(field, loop_quiver)
    x, y = parent.generator("x"), parent.generator("y")
    relations = [x.pow(2) + y.pow(3)]
    relations.extend(
        left * middle_left * middle_right * right
        for left, middle_left, middle_right, right in itertools.product(
            (x, y), repeat=4
        )
    )
    quotient = parent.quotient(relations)
    algebra = qpa.QuiverAlgebra.from_quotient_algebra(quotient)
    assert quotient.is_monomial is expected["nonmonomial"]

    monomial = algebra.associated_monomial_algebra()
    monomial_quotient = monomial.quotient_algebra
    assert monomial_quotient is not None
    assert monomial.is_identical(algebra) is expected["monomial_identity"]
    assert monomial_quotient.is_monomial is expected["monomial_result_is_monomial"]
    assert monomial_quotient.dimension() == str(expected["monomial_dimension"])
    assert _relation_paths(monomial) == expected["monomial_relations"]

    graded = algebra.associated_graded_algebra()
    graded_quotient = graded.quotient_algebra
    assert graded_quotient is not None
    assert graded.is_identical(algebra) is expected["graded_identity"]
    assert graded_quotient.dimension() == str(expected["graded_dimension"])
    assert _relation_paths(graded) == expected["graded_relations"]

    truncated = qpa.QuotientAlgebra.truncated(field, loop_quiver, 3)
    truncated_union = qpa.QuiverAlgebra.from_quotient_algebra(truncated)
    assert truncated_union.associated_monomial_algebra().is_identical(
        truncated_union
    ) is expected["already_monomial_identity"]


def test_associated_graded_algebra_rejects_qpa_invalid_domains() -> None:
    field = qpa.Field.prime(2)
    loop_quiver = qpa.Quiver(["v"], [("v", "v", "x")])
    cyclic = qpa.QuiverAlgebra.from_path_algebra(qpa.PathAlgebra(field, loop_quiver))
    with pytest.raises(qpa.InvalidInputError, match="not finite dimensional"):
        cyclic.associated_graded_algebra()

    parent = qpa.PathAlgebra(field, loop_quiver)
    nonadmissible = qpa.QuiverAlgebra.from_quotient_algebra(
        parent.quotient([parent.generator("x")])
    )
    with pytest.raises(qpa.InvalidInputError, match="not an admissible quotient"):
        nonadmissible.associated_graded_algebra()
