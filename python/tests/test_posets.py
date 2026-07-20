from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _matrix(poset: qpa.Poset) -> list[list[bool]]:
    order = poset.partial_order
    return [
        [order(left, right) for right in poset.underlying_set]
        for left in poset.underlying_set
    ]


def _assert_algebra_shape(algebra: qpa.QuiverAlgebra, expected: dict[str, Any]) -> None:
    assert algebra.is_quotient is (not expected["path_algebra"])
    assert algebra.quiver.vertices == expected["vertices"]
    assert [name for _source, _target, name in algebra.quiver.arrows] == [
        arrow["name"] for arrow in expected["arrows"]
    ]
    if algebra.path_algebra is not None:
        basis = algebra.path_algebra.canonical_basis()
        assert basis is not None and len(basis) == expected["dimension"]
        relations = algebra.path_algebra.relations
    else:
        quotient = algebra.quotient_algebra
        assert quotient is not None and quotient.dimension() == str(expected["dimension"])
        relations = quotient.relations
    assert [
        [path.external_representation for path, _coefficient in relation.terms]
        for relation in relations
    ] == [[term["path"] for term in relation] for relation in expected["relations"]]


def test_poset_objects_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "poset.schema-v1.json").read_text(encoding="utf-8")
    )
    example = qpa.Poset(
        ["d", "b", "a", "c"],
        [["a", "b", "c"], ["b", "d"], ["c", "d"]],
    )
    assert example.is_poset and example.is_default_rep
    assert len(example) == fixture["example"]["size"]
    assert example.underlying_set == fixture["example"]["underlying"]
    assert [list(pair) for pair in example.minimal_relations] == fixture["example"][
        "minimal"
    ]
    assert _matrix(example) == fixture["example"]["matrix"]
    assert example.less_than_or_equal("a", "d")
    assert not example.less_than_or_equal("b", "c")
    assert example.is_identical(example)
    assert not example.is_identical(qpa.Poset(["a"], []))
    assert qpa.Poset.family_name == "TheFamilyOfPosets"
    assert qpa.Poset.type_name == "TheTypePoset"

    expected_output = (FIXTURES / "poset_output.txt").read_text(encoding="utf-8")
    view, display = expected_output.split("\n", 1)
    assert repr(example) == view
    assert example.display() == display

    duplicates = qpa.Poset([3, 1, 2, 1], [[1, 2], [2, 3]])
    assert [str(point) for point in duplicates.underlying_set] == fixture["duplicates"][
        "underlying"
    ]
    assert _matrix(duplicates) == fixture["duplicates"]["matrix"]


def test_poset_errors_follow_qpa() -> None:
    with pytest.raises(qpa.InvalidInputError, match="cannot be empty"):
        qpa.Poset([], [])
    with pytest.raises(qpa.InvalidInputError, match="outside the underlying set"):
        qpa.Poset([1, 2], [[1, 3]])
    with pytest.raises(qpa.InvalidInputError, match="do not define a partial order"):
        qpa.Poset([1, 2], [[1, 2], [2, 1]])
    with pytest.raises(qpa.InvalidInputError, match="outside the underlying set"):
        qpa.Poset([1, 2], []).less_than_or_equal(1, 3)
    with pytest.raises(qpa.InvalidInputError, match="integer or string"):
        qpa.Poset([object()], [])


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_poset_algebra_follows_qpa(field: qpa.Field) -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_poset_algebra.schema-v1.json").read_text(encoding="utf-8")
    )
    diamond = qpa.Poset([4, 2, 1, 3], [[1, 2, 3], [2, 4], [3, 4]])
    result = qpa.poset_algebra(field, diamond)
    assert result.is_poset_algebra
    assert result.is_finite_global_dimension_algebra
    assert result.poset.is_identical(diamond)
    _assert_algebra_shape(result.algebra, fixture["rational"]["diamond"])

    if field.characteristic == 2:
        parallel = qpa.Poset(
            [1, 2, 3, 4, 5],
            [[1, 2, 3, 4], [2, 5], [3, 5], [4, 5]],
        )
        parallel_result = qpa.poset_algebra(field, parallel)
        assert parallel_result.poset.is_identical(parallel)
        _assert_algebra_shape(parallel_result.algebra, fixture["f2"]["parallel"])
