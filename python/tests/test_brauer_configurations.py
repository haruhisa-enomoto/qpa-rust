from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


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


def _documented_configuration() -> qpa.BrauerConfiguration:
    return qpa.BrauerConfiguration(
        [("v1", 1), ("v2", 1), ("v3", 2)],
        [("e1", ["v1", "v2"]), ("e2", ["v2", "v3"])],
        [["e1"], ["e2", "e1"], ["e2"]],
    )


def _repeated_configuration() -> qpa.BrauerConfiguration:
    return qpa.BrauerConfiguration(
        [("a", 1), ("b", 1), ("c", 1)],
        [
            ("1", ["a", "a"]),
            ("2", ["a", "b"]),
            ("3", ["b", "a"]),
            ("4", ["b", "c"]),
        ],
        [["1", "2", "3"], ["2", "4", "3"], ["4"]],
    )


def test_brauer_configuration_algebras_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_brauer_configuration.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )

    documented_configuration = _documented_configuration()
    documented_result = qpa.brauer_configuration_algebra(
        qpa.Field.rational(), documented_configuration
    )
    assert documented_result.is_valid
    assert documented_result.diagnostic is None
    assert documented_result.validation_error is None
    documented = documented_result.algebra
    assert documented is not None and documented.is_symmetric_algebra
    _assert_algebra_shape(documented.algebra, fixture["rational"]["documented"])
    assert documented_configuration.orientations == fixture["rational"][
        "orientations_after"
    ]

    repeated_configuration = _repeated_configuration()
    repeated_result = qpa.brauer_configuration_algebra(
        qpa.Field.prime(2), repeated_configuration
    )
    repeated = repeated_result.algebra
    assert repeated_result.is_valid and repeated is not None
    assert repeated.is_symmetric_algebra
    _assert_algebra_shape(repeated.algebra, fixture["f2"]["repeated"])
    assert repeated_configuration.orientations == fixture["f2"]["orientations_after"]

    # The translated constructor is field-generic, not specialized to the fixtures.
    trivial = qpa.BrauerConfiguration([("v", 1)], [("e", ["v"])], [["e"]])
    trivial_result = qpa.brauer_configuration_algebra(qpa.Field.prime(5), trivial)
    trivial_algebra = trivial_result.algebra
    assert trivial_algebra is not None
    _assert_algebra_shape(trivial_algebra.algebra, fixture["rational"]["trivial"])


def test_brauer_configuration_failures_and_mutation_follow_qpa() -> None:
    no_vertices = qpa.BrauerConfiguration([], [], [])
    result = qpa.brauer_configuration_algebra(qpa.Field.rational(), no_vertices)
    assert not result.is_valid and result.algebra is None
    assert result.validation_error == "no_vertices"
    assert result.diagnostic == "There must be at least 1 vertex\n"

    containment = qpa.BrauerConfiguration(
        [("a", 1), ("b", 1)],
        [("z", ["a"]), ("x", ["b"])],
        [["z"], ["z", "x"]],
    )
    containment_result = qpa.brauer_configuration_algebra(
        qpa.Field.rational(), containment
    )
    assert containment_result.validation_error == "orientation_containment_mismatch"
    assert containment_result.diagnostic == (
        "Orientations must contain exactly the edges or polygons which contain "
        "the corresponding vertex.\n"
    )
    assert containment.orientations == [["z"], ["x", "z"]]

    with pytest.raises(TypeError, match="exact Python integer"):
        qpa.BrauerConfiguration([("v", "one")], [], [])
