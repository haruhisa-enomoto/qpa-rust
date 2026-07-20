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


def _element_data(element: qpa.QuotientElement) -> list[dict[str, Any]]:
    return [
        {
            "path": path.external_representation,
            "coefficient": _coefficient_data(coefficient),
        }
        for path, coefficient in element.representative.terms
    ]


def _three_vertex_quotient(
    field: qpa.Field, c_coefficient: qpa.Scalar
) -> qpa.QuiverAlgebra:
    algebra = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [("v1", "v2", "a"), ("v2", "v3", "b"), ("v1", "v3", "c")],
        ),
    )
    a, b, c = (algebra.generator(name) for name in ("a", "b", "c"))
    quotient = algebra.quotient([a * b + c.scale_left(c_coefficient)])
    return qpa.QuiverAlgebra.from_quotient_algebra(quotient)


def test_save_and_read_path_algebra_follow_qpa(tmp_path: Path) -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_persistence.schema-v1.json").read_text(encoding="utf-8")
    )["path"]
    vertices = [str(index) for index in range(1, 12)]
    arrows = [
        (str(index), str(index + 1), f"a_{index}") for index in range(1, 11)
    ]
    algebra = qpa.PathAlgebra(qpa.Field.rational(), qpa.Quiver(vertices, arrows))
    union = qpa.QuiverAlgebra.from_path_algebra(algebra)
    target = tmp_path / "path-algebra.txt"

    assert union.qpa_text() == fixture["text"]
    assert union.save_algebra(str(target), "delete") is fixture["saved"]
    assert target.read_text(encoding="utf-8") == fixture["text"]
    loaded = qpa.read_algebra(str(target))
    assert (loaded.path_algebra is not None) is fixture["read_is_path"]
    assert loaded.quotient_algebra is None
    assert loaded.quiver.number_of_vertices == fixture["read_vertices"]
    assert loaded.quiver.number_of_arrows == fixture["read_arrows"]


@pytest.mark.parametrize(
    ("field", "coefficient", "fixture_key"),
    [
        (qpa.Field.rational(), (-1, 2), "rational"),
        (qpa.Field.prime(2), (1, None), "f2"),
        (qpa.Field.prime(5), (2, None), "f5"),
    ],
)
def test_save_and_read_quotients_follow_qpa(
    tmp_path: Path,
    field: qpa.Field,
    coefficient: tuple[int, int | None],
    fixture_key: str,
) -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_persistence.schema-v1.json").read_text(encoding="utf-8")
    )[fixture_key]
    numerator, denominator = coefficient
    scalar = field.element(numerator, denominator)
    algebra = _three_vertex_quotient(field, scalar)
    target = tmp_path / f"{fixture_key}.txt"

    assert algebra.qpa_text() == fixture["text"]
    assert algebra.save_algebra(str(target), "delete old file")
    loaded = qpa.read_algebra(str(target))
    assert loaded.path_algebra is None
    assert (loaded.quotient_algebra is not None) is fixture["read_is_quotient"]
    quotient = loaded.quotient_algebra
    assert quotient is not None
    assert quotient.dimension() == str(fixture["read_dimension"])


def test_save_algebra_overwrite_branches_preserve_qpa_behavior(tmp_path: Path) -> None:
    algebra = qpa.QuiverAlgebra.from_path_algebra(
        qpa.PathAlgebra(qpa.Field.rational(), qpa.Quiver(["v"], []))
    )
    target = tmp_path / "algebra.txt"
    assert algebra.save_algebra(str(target), "delete")
    original = target.read_text(encoding="utf-8")
    assert not algebra.save_algebra(str(target), "keep existing")
    assert not algebra.save_algebra(str(target), "query", query_existing=False)
    assert algebra.save_algebra(str(target), "query", query_existing=True)
    assert target.read_text(encoding="utf-8") == original

    with pytest.raises(qpa.InvalidInputError, match="invalid SaveAlgebra"):
        algebra.save_algebra(str(target), "bad")
    target.unlink()
    with pytest.raises(qpa.InvalidInputError, match="leaves output unbound"):
        algebra.save_algebra(str(target), "query", query_existing=True)


def test_quotient_nice_basis_follows_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_conversion.schema-v1.json").read_text(encoding="utf-8")
    )["rational"]
    field = qpa.Field.rational()
    algebra = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [("v1", "v2", "a"), ("v2", "v3", "b"), ("v1", "v3", "c")],
        ),
    )
    a, b, c = (algebra.generator(name) for name in ("a", "b", "c"))
    quotient = algebra.quotient_high_level([a * b - c])
    v1, quotient_a, quotient_c = (
        quotient.generator(name) for name in ("v1", "a", "c")
    )
    basis = quotient.nice_basis([v1, quotient_a, quotient_c])
    element = (
        v1.scale_left(field.element(2))
        + quotient_a.scale_left(field.element(3))
        + quotient_c.scale_left(field.element(5))
    )

    assert basis.quotient.is_identical(quotient)
    assert [path.external_representation for path in basis.monomials] == [[1], [4], [6]]
    coordinates = basis.nice_vector(element)
    assert coordinates is not None
    assert [_coefficient_data(value) for value in coordinates] == fixture["nice"]
    rebuilt = basis.ugly_vector(coordinates)
    assert rebuilt is not None
    assert _element_data(rebuilt) == fixture["ugly"]
    assert basis.nice_vector(quotient.generator("b")) is None

    empty = quotient.nice_basis([])
    zero_coordinates = empty.nice_vector(quotient.zero())
    assert zero_coordinates == [field.element(0)]
    assert empty.element_from_nice_vector(zero_coordinates or []) == quotient.zero()
    assert empty.ugly_vector([field.element(1)]) is None


@pytest.mark.parametrize(
    ("field", "key"),
    [(qpa.Field.rational(), "rational"), (qpa.Field.prime(2), "f2")],
)
def test_path_algebra_triangular_reduced_property_follows_qpa(
    field: qpa.Field, key: str
) -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_triangular_reduced.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )[key]
    one_vertex = qpa.PathAlgebra(field, qpa.Quiver(["v"], []))
    arrow = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    assert one_vertex.is_triangular_reduced is fixture["one_vertex_path"]
    assert arrow.is_triangular_reduced is fixture["arrow_path"]
