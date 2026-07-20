from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _assert_algebra_matches(algebra: qpa.QuiverAlgebra, expected: dict[str, Any]) -> None:
    assert algebra.is_quotient is (not expected["path_algebra"])
    assert algebra.quiver.vertices == expected["vertices"]
    assert [
        {
            "name": name,
            "source": expected["vertices"].index(source) + 1,
            "target": expected["vertices"].index(target) + 1,
        }
        for source, target, name in algebra.quiver.arrows
    ] == expected["arrows"]

    if algebra.path_algebra is not None:
        basis = algebra.path_algebra.canonical_basis()
        assert basis is not None and len(basis) == expected["dimension"]
        relations: list[qpa.AlgebraElement] = algebra.path_algebra.relations
    else:
        quotient = algebra.quotient_algebra
        assert quotient is not None
        assert quotient.dimension() == str(expected["dimension"])
        relations = quotient.relations
    assert [
        [path.external_representation for path, _coefficient in relation.terms]
        for relation in relations
    ] == [[term["path"] for term in relation] for relation in expected["relations"]]


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_nakayama_canonical_and_kronecker_follow_qpa(field: qpa.Field) -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_predefined.schema-v1.json").read_text(encoding="utf-8")
    )
    finite_expected = fixture["f2"]

    nakayama_result = qpa.nakayama_algebra(field, [3, 3, 3])
    assert nakayama_result.is_valid
    assert nakayama_result.invalid_sequence is None
    assert nakayama_result.diagnostic is None
    nakayama = nakayama_result.algebra
    assert nakayama is not None
    assert nakayama.kind == "Nakayama"
    assert nakayama.is_nakayama_algebra is True
    assert nakayama.is_finite_type_algebra is True
    assert nakayama.is_canonical_algebra is None
    _assert_algebra_matches(
        nakayama.algebra,
        fixture["rational" if field.characteristic is None else "f2"][
            "nakayama_cyclic"
        ],
    )

    if field.characteristic is None:
        canonical = qpa.canonical_algebra(
            field, [2, 2, 3, 4], [field.element(2), field.element(-3)]
        )
        canonical_expected = fixture["rational"]["canonical_four"]
    else:
        canonical = qpa.canonical_algebra(field, [2, 3, 4], [field.element(1)])
        canonical_expected = finite_expected["canonical"]
    assert canonical.kind == "Canonical"
    assert canonical.is_canonical_algebra is True
    assert canonical.is_finite_global_dimension_algebra is True
    _assert_algebra_matches(canonical.algebra, canonical_expected)

    kronecker = qpa.kronecker_algebra(field, 3)
    assert kronecker.kind == "Kronecker"
    assert kronecker.is_kronecker_algebra is True
    assert kronecker.is_finite_global_dimension_algebra is True
    assert kronecker.is_hereditary_algebra is True
    _assert_algebra_matches(kronecker.algebra, finite_expected["kronecker"])


def test_invalid_nakayama_sequence_preserves_qpa_result() -> None:
    result = qpa.nakayama_algebra(qpa.Field.rational(), [1, 1])
    assert not result.is_valid
    assert result.algebra is None
    assert result.invalid_sequence == [1, 1]
    assert result.diagnostic == "The admissible sequence is not a valid one.\n"
    with pytest.raises(
        qpa.InvalidInputError, match="admissible sequence must not be empty"
    ):
        qpa.nakayama_algebra(qpa.Field.rational(), [])


def test_admissible_sequences_and_remaining_predefined_algebras_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_predefined_remaining.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    assert qpa.admissible_sequence_generator(1, 3) == fixture["sequences"][
        "one_by_three"
    ]
    assert qpa.admissible_sequence_generator(3, 3) == fixture["sequences"][
        "three_by_three"
    ]
    assert qpa.admissible_sequence_generator(4, 2) == fixture["sequences"][
        "four_by_two"
    ]
    assert qpa.admissible_sequence_generator(3, 1) == fixture["sequences"][
        "three_by_one"
    ]

    rationals = qpa.Field.rational()
    a2 = qpa.PathAlgebra(
        rationals, qpa.Quiver(["v1", "v2"], [("v1", "v2", "x")])
    )
    _assert_algebra_matches(
        a2.preprojective_algebra(), fixture["rational"]["preprojective_a2"]
    )
    _assert_algebra_matches(
        qpa.finite_global_dim_green(rationals, 2), fixture["rational"]["green_2"]
    )
    _assert_algebra_matches(
        qpa.finite_global_dim_kirk_kuz(rationals, 1),
        fixture["rational"]["kirk_kuz_1"],
    )

    field2 = qpa.Field.prime(2)
    a3 = qpa.PathAlgebra(
        field2,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [("v1", "v2", "u"), ("v2", "v3", "v")],
        ),
    )
    _assert_algebra_matches(
        a3.preprojective_algebra(), fixture["f2"]["preprojective_a3"]
    )
    _assert_algebra_matches(
        qpa.finite_global_dim_green(field2, 2), fixture["f2"]["green_2"]
    )
    _assert_algebra_matches(
        qpa.finite_global_dim_kirk_kuz(field2, 2), fixture["f2"]["kirk_kuz_2"]
    )
