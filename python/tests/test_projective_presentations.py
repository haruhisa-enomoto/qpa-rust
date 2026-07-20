from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _fixture_objects(
    field: qpa.Field,
) -> tuple[qpa.Quiver, qpa.PathAlgebra, list[list[qpa.AlgebraElement]]]:
    quiver = qpa.Quiver(
        ["v1", "v2", "v3"],
        [("v1", "v2", "a"), ("v1", "v3", "b"), ("v2", "v3", "c")],
    )
    algebra = qpa.PathAlgebra(field, quiver)
    a, b, c, v2 = (algebra.generator(name) for name in ("a", "b", "c", "v2"))
    return quiver, algebra, [[a + b, algebra.zero()], [c, v2]]


def _coefficient(value: qpa.Scalar) -> list[int]:
    if value.numerator is not None:
        return [int(value.numerator), int(value.denominator)]
    representative = value.representative
    assert representative is not None
    return [representative, 1]


def _matrix_data(matrix: list[list[qpa.AlgebraElement]]) -> list[list[list[dict[str, Any]]]]:
    return [
        [
            [
                {"path": repr(path), "coefficient": _coefficient(coefficient)}
                for path, coefficient in element.terms
            ]
            for element in row
        ]
        for row in matrix
    ]


def test_projective_presentation_representation_follows_qpa() -> None:
    expected = json.loads(
        (FIXTURES / "present2.schema-v1.json").read_text(encoding="utf-8")
    )
    quiver, algebra, matrix = _fixture_objects(qpa.Field.rational())
    presentation = qpa.ProjectivePresentation.from_attributes(
        algebra,
        matrix,
        2,
        2,
        [quiver.vertex("v1"), quiver.vertex("v2")],
    )

    assert presentation.is_projective_presentation is expected["category"]
    assert presentation.is_default_rep is expected["default_rep"]
    assert (
        qpa.ProjectivePresentation.family_name == "ProjectivePresentationFamily"
    ) is expected["global_family"]
    assert presentation.parent_algebra.is_identical(algebra) is expected[
        "parent_identical"
    ]
    assert (presentation.matrix_representation == matrix) is expected["matrix_equal"]
    assert _matrix_data(presentation.matrix_representation) == expected["matrix"]
    assert [presentation.matrix_num_rows, presentation.matrix_num_cols] == expected[
        "shape"
    ]
    assert [vertex.name for vertex in presentation.p0_vertex_list] == expected["p0"]
    assert presentation.has_p1_vertex_list is expected["has_p1"]
    assert presentation.p1_vertex_list is None


@pytest.mark.parametrize(
    "field", [qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_projective_presentation_representation_is_field_generic(field: qpa.Field) -> None:
    quiver, algebra, matrix = _fixture_objects(field)
    presentation = qpa.ProjectivePresentation.from_attributes(
        algebra,
        matrix,
        7,
        9,
        [quiver.vertex("v1")],
        [quiver.vertex("v3")],
    )
    assert [presentation.matrix_num_rows, presentation.matrix_num_cols] == [7, 9]
    assert [vertex.name for vertex in presentation.p1_vertex_list or []] == ["v3"]
    assert presentation.has_p1_vertex_list


def test_pinned_projective_presentation_constructor_failures() -> None:
    _quiver, algebra, matrix = _fixture_objects(qpa.Field.rational())
    with pytest.raises(qpa.InvalidInputError, match="no ProjectivePresentation method"):
        qpa.ProjectivePresentation(algebra, [])
    with pytest.raises(qpa.InvalidInputError, match="rows of equal length"):
        qpa.ProjectivePresentation(algebra, [[matrix[0][0]], matrix[0]])
    with pytest.raises(qpa.InvalidInputError, match="all-zero row"):
        qpa.ProjectivePresentation(algebra, [[algebra.zero()]])
    with pytest.raises(
        qpa.InvalidInputError, match="CoefficientsAndMagmaElements method"
    ):
        qpa.ProjectivePresentation(algebra, matrix)
