from __future__ import annotations

import json
from pathlib import Path

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _assert_cartan_coxeter(
    algebra: qpa.QuiverAlgebra, expected: dict[str, object]
) -> None:
    assert algebra.cartan_matrix() == expected["cartan"]
    matrix = algebra.coxeter_matrix()
    polynomial = algebra.coxeter_polynomial()
    assert matrix is not None and polynomial is not None
    assert [[repr(coefficient) for coefficient in row] for row in matrix] == expected[
        "coxeter"
    ]
    assert [repr(coefficient) for coefficient in polynomial] == expected[
        "polynomial_coefficients"
    ]


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_cartan_and_coxeter_invariants_follow_qpa(field: qpa.Field) -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_cartan_coxeter.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    key = "rational" if field.characteristic is None else f"f{field.characteristic}"

    a3_path = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [("v1", "v2", "a"), ("v2", "v3", "b")],
        ),
    )
    a3 = qpa.QuiverAlgebra.from_path_algebra(a3_path)
    _assert_cartan_coxeter(a3, fixture[key]["path_a3"])

    loop_path = qpa.PathAlgebra(field, qpa.Quiver(["v"], [("v", "v", "x")]))
    x = loop_path.generator("x")
    dual_numbers = qpa.QuiverAlgebra.from_quotient_algebra(loop_path.quotient([x * x]))
    _assert_cartan_coxeter(dual_numbers, fixture[key]["quotient_dual_numbers"])


def test_cartan_and_coxeter_fail_results_follow_qpa() -> None:
    loop_path = qpa.PathAlgebra(
        qpa.Field.rational(), qpa.Quiver(["v"], [("v", "v", "x")])
    )
    loop = qpa.QuiverAlgebra.from_path_algebra(loop_path)
    assert loop.cartan_matrix() is None
    assert loop.coxeter_matrix() is None
    assert loop.coxeter_polynomial() is None


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_symmetry_and_path_algebra_properties_follow_qpa(field: qpa.Field) -> None:
    symmetry = json.loads(
        (FIXTURES / "pathalg_symmetry_properties.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )["rational" if field.characteristic is None else "f2"]
    homological = json.loads(
        (FIXTURES / "algebra_homological_properties.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )["rational" if field.characteristic is None else "f2"]

    discrete = qpa.QuiverAlgebra.from_path_algebra(
        qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], []))
    )
    assert discrete.is_selfinjective_algebra() == homological["path_discrete"][
        "selfinjective"
    ]
    assert discrete.is_symmetric_algebra() == symmetry["path_discrete"]["symmetric"]
    assert discrete.is_weakly_symmetric_algebra() == symmetry["path_discrete"][
        "weakly_symmetric"
    ]
    assert discrete.is_schurian_algebra() == symmetry["path_discrete"]["schurian"]
    assert discrete.is_semicommutative_algebra() == symmetry["path_discrete"][
        "semicommutative"
    ]
    assert not discrete.is_distributive_algebra()

    a2 = qpa.QuiverAlgebra.from_path_algebra(
        qpa.PathAlgebra(
            field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
        )
    )
    assert a2.is_selfinjective_algebra() == homological["path_a2"]["selfinjective"]
    assert a2.is_symmetric_algebra() == symmetry["path_a2"]["symmetric"]
    assert a2.is_weakly_symmetric_algebra() == symmetry["path_a2"][
        "weakly_symmetric"
    ]
    assert a2.is_schurian_algebra() == symmetry["path_a2"]["schurian"]
    assert a2.is_semicommutative_algebra() == symmetry["path_a2"][
        "semicommutative"
    ]
    assert a2.is_distributive_algebra()

    loop_path = qpa.PathAlgebra(field, qpa.Quiver(["v"], [("v", "v", "x")]))
    x = loop_path.generator("x")
    dual_numbers = qpa.QuiverAlgebra.from_quotient_algebra(loop_path.quotient([x * x]))
    assert dual_numbers.is_selfinjective_algebra() == homological[
        "quotient_dual_numbers"
    ]["selfinjective"]
    assert dual_numbers.is_symmetric_algebra() == symmetry["quotient_dual_numbers"][
        "symmetric"
    ]
    assert dual_numbers.is_weakly_symmetric_algebra() == symmetry[
        "quotient_dual_numbers"
    ]["weakly_symmetric"]
    assert dual_numbers.is_schurian_algebra() == symmetry["quotient_dual_numbers"][
        "schurian"
    ]
    assert dual_numbers.is_semicommutative_algebra() == symmetry[
        "quotient_dual_numbers"
    ]["semicommutative"]


def test_quotient_distributivity_follows_qpa() -> None:
    field = qpa.Field.rational()
    path = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3", "v4", "v5", "v6", "v7"],
            [
                ("v1", "v2", "a"),
                ("v1", "v3", "b"),
                ("v2", "v4", "c"),
                ("v3", "v4", "d"),
                ("v5", "v6", "u"),
                ("v6", "v7", "z"),
            ],
        ),
    )
    relation = path.generator("u") * path.generator("z")
    diamond = qpa.QuiverAlgebra.from_quotient_algebra(path.quotient([relation]))
    expected = json.loads(
        (FIXTURES / "reptype_quotient.schema-v1.json").read_text(encoding="utf-8")
    )
    assert diamond.is_distributive_algebra() == expected["distributive"]["diamond"]
