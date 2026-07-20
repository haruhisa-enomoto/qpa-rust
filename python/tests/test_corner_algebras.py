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


def _relations(algebra: qpa.QuiverAlgebra) -> list[qpa.AlgebraElement]:
    if algebra.path_algebra is not None:
        return algebra.path_algebra.relations
    quotient = algebra.quotient_algebra
    assert quotient is not None
    return quotient.relations


def _dimension(algebra: qpa.QuiverAlgebra) -> int:
    if algebra.path_algebra is not None:
        basis = algebra.path_algebra.canonical_basis()
        assert basis is not None
        return len(basis)
    quotient = algebra.quotient_algebra
    assert quotient is not None
    dimension = quotient.dimension()
    assert dimension is not None
    return int(dimension)


def _external_relations(algebra: qpa.QuiverAlgebra) -> list[list[dict[str, Any]]]:
    return [
        [
            {
                "path": path.external_representation,
                "coefficient": _coefficient_data(coefficient),
            }
            for path, coefficient in relation.terms
        ]
        for relation in _relations(algebra)
    ]


def _quiver_data(quiver: qpa.Quiver) -> dict[str, Any]:
    return {
        "vertices": quiver.vertices,
        "arrows": [
            {"source": source, "target": target, "name": name}
            for source, target, name in quiver.arrows
        ],
    }


def _vertex_quotient_input(field: qpa.Field) -> tuple[qpa.PathAlgebra, qpa.QuiverAlgebra]:
    parent = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3", "v4"],
            [
                ("v1", "v2", "a"),
                ("v2", "v3", "b"),
                ("v1", "v3", "c"),
                ("v3", "v4", "d"),
                ("v2", "v4", "e"),
                ("v1", "v4", "f"),
            ],
        ),
    )
    a, b, c, d, e, f = (parent.generator(name) for name in "abcdef")
    relations = [a * b - c, b * d - e, c * d - f]
    quotient = parent.quotient(relations)
    return parent, qpa.QuiverAlgebra.from_quotient_algebra(quotient)


def _assert_vertex_quotient(
    algebra: qpa.QuiverAlgebra, expected: dict[str, Any]
) -> None:
    assert _quiver_data(algebra.quiver) == expected["quiver"]
    assert _external_relations(algebra) == expected["relations"]
    assert _dimension(algebra) == expected["dimension"]


def test_vertex_quotients_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_vertex_quotient.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    parent, algebra = _vertex_quotient_input(qpa.Field.rational())
    rational = fixture["rational"]
    remove_two = algebra.quiver_algebra_of_a_mod_aea([2])
    remove_three = algebra.quiver_algebra_of_a_mod_aea([3])
    none = algebra.quiver_algebra_of_a_mod_aea([])
    ignored = algebra.quiver_algebra_of_a_mod_aea([0, 5, 99, -1, 10**80])
    _assert_vertex_quotient(remove_two, rational["remove_two"])
    _assert_vertex_quotient(remove_three, rational["remove_three"])
    _assert_vertex_quotient(none, rational["none"])
    assert _quiver_data(ignored.quiver) == _quiver_data(none.quiver)
    assert rational["ignored_same_quiver"] is True

    path_result = qpa.QuiverAlgebra.from_path_algebra(
        parent
    ).quiver_algebra_of_a_mod_aea([2])
    assert path_result.is_quotient is rational["path_algebra_input"]["is_quotient"]
    assert _external_relations(path_result) == rational["path_algebra_input"]["relations"]

    _parent2, algebra2 = _vertex_quotient_input(qpa.Field.prime(2))
    _assert_vertex_quotient(
        algebra2.quiver_algebra_of_a_mod_aea([2]), fixture["f2"]["remove_two"]
    )

    _parent5, algebra5 = _vertex_quotient_input(qpa.Field.prime(5))
    assert _dimension(algebra5.quiver_algebra_of_a_mod_aea([2])) == rational["remove_two"][
        "dimension"
    ]


def _corner_input(
    field: qpa.Field,
) -> tuple[qpa.PathAlgebra, qpa.QuotientAlgebra, qpa.QuiverAlgebra]:
    parent = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [("v1", "v2", "a"), ("v2", "v3", "b")],
        ),
    )
    quotient = parent.quotient([parent.generator("a") * parent.generator("b")])
    return parent, quotient, qpa.QuiverAlgebra.from_quotient_algebra(quotient)


def _corner_data(algebra: qpa.QuiverAlgebra) -> dict[str, Any]:
    return {
        "is_path": not algebra.is_quotient,
        "dimension": _dimension(algebra),
        "adjacency": algebra.quiver.adjacency_matrix,
        "relations": [
            [
                {
                    "path": repr(path),
                    "coefficient": _coefficient_data(coefficient),
                }
                for path, coefficient in relation.terms
            ]
            for relation in _relations(algebra)
        ],
    }


def _union_element(element: qpa.QuotientElement) -> qpa.QuiverAlgebraElement:
    return qpa.QuiverAlgebraElement.from_quotient_element(element)


def test_idempotent_corners_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_eae.schema-v1.json").read_text(encoding="utf-8")
    )
    _parent, quotient, algebra = _corner_input(qpa.Field.prime(2))
    v1, v2, v3, a = (quotient.generator(name) for name in ("v1", "v2", "v3", "a"))
    corners = {
        "corner12": algebra.quiver_algebra_of_eae(_union_element(v1 + v2)),
        "corner13": algebra.quiver_algebra_of_eae(_union_element(v1 + v3)),
        "conjugate": algebra.quiver_algebra_of_eae(_union_element(v1 + a)),
        "whole": algebra.quiver_algebra_of_eae(_union_element(quotient.one())),
    }
    assert {name: _corner_data(corner) for name, corner in corners.items()} == {
        name: fixture["binary"][name] for name in corners
    }

    minimal_parent = qpa.PathAlgebra(
        qpa.Field.prime(2),
        qpa.Quiver(
            ["v1", "v2", "v3", "v4"],
            [
                ("v1", "v2", "a"),
                ("v2", "v3", "b"),
                ("v3", "v4", "c"),
            ],
        ),
    )
    minimal_quotient = minimal_parent.quotient(
        [minimal_parent.generator("a") * minimal_parent.generator("b")]
    )
    minimal = qpa.QuiverAlgebra.from_quotient_algebra(minimal_quotient)
    minimal_corner = minimal.quiver_algebra_of_eae(_union_element(minimal_quotient.one()))
    assert _corner_data(minimal_corner) == fixture["binary"]["minimal"]

    _parent_q, quotient_q, algebra_q = _corner_input(qpa.Field.rational())
    e12_q = _union_element(quotient_q.generator("v1") + quotient_q.generator("v2"))
    assert _corner_data(algebra_q.quiver_algebra_of_eae(e12_q)) == fixture["rational"][
        "corner12"
    ]
    assert fixture["rational"]["success"] is True


def test_idempotent_corner_errors_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_eae.schema-v1.json").read_text(encoding="utf-8")
    )["errors"]
    parent, quotient, algebra = _corner_input(qpa.Field.prime(2))

    path_input = qpa.QuiverAlgebra.from_path_algebra(parent)
    path_idempotent = qpa.QuiverAlgebraElement.from_path_element(parent.generator("v1"))
    foreign_parent = qpa.PathAlgebra(qpa.Field.prime(2), qpa.Quiver(["w"], []))
    foreign = qpa.QuiverAlgebraElement.from_path_element(foreign_parent.generator("w"))
    nonidempotent = _union_element(quotient.generator("a"))
    zero = _union_element(quotient.zero())

    for key, source, idempotent in (
        ("path_success", path_input, path_idempotent),
        ("foreign_success", algebra, foreign),
        ("nonidempotent_success", algebra, nonidempotent),
        ("zero_success", algebra, zero),
    ):
        with pytest.raises(qpa.InvalidInputError):
            source.quiver_algebra_of_eae(idempotent)
        assert fixture[key] is False
