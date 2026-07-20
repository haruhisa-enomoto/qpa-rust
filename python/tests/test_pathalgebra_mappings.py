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


def _named_terms(element: qpa.QuiverAlgebraElement) -> list[dict[str, Any]]:
    if element.path_element is not None:
        terms = element.path_element.terms
    else:
        quotient_element = element.quotient_element
        assert quotient_element is not None
        terms = quotient_element.representative.terms
    return [
        {"path": repr(path), "coefficient": _coefficient_data(coefficient)}
        for path, coefficient in terms
    ]


def _external_terms(element: qpa.QuotientElement) -> list[dict[str, Any]]:
    return [
        {
            "path": path.external_representation,
            "coefficient": _coefficient_data(coefficient),
        }
        for path, coefficient in element.representative.terms
    ]


def _as_path_element(element: qpa.AlgebraElement) -> qpa.QuiverAlgebraElement:
    return qpa.QuiverAlgebraElement.from_path_element(element)


def _as_quotient_element(element: qpa.QuotientElement) -> qpa.QuiverAlgebraElement:
    return qpa.QuiverAlgebraElement.from_quotient_element(element)


def test_mapped_expression_and_images_representative_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_mapping_ordering.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    field = qpa.Field.prime(5)
    source = qpa.PathAlgebra(
        field, qpa.Quiver(["v1"], [("v1", "v1", "a"), ("v1", "v1", "b")])
    )
    target = qpa.PathAlgebra(
        field, qpa.Quiver(["v1"], [("v1", "v1", "x"), ("v1", "v1", "y")])
    )
    v1, a, b = source.generators
    target_v1, x, y = target.generators
    element = (
        v1.scale_left(field.element(2))
        + (a * a).scale_left(field.element(3))
        + (a * b).scale_left(field.element(4))
        + b
    )
    mapping = qpa.AlgebraHomomorphism.path_to_path(
        source,
        target,
        [target_v1, y, x],
        generators_with_one=True,
    )
    union = _as_path_element(element)
    mapped = mapping.mapped_expression(union)
    assert _named_terms(mapped) == fixture["path"]["mapped"]
    assert mapping.images_representative(union) == mapped
    assert mapping.image(union) == mapped
    assert _named_terms(mapping.mapped_expression(_as_path_element(source.zero()))) == fixture[
        "path"
    ]["mapped_zero"]

    source_quotient = source.quotient([a * a])
    target_quotient = target.quotient([y * y])
    quotient_v1, quotient_a, quotient_b = source_quotient.generators
    target_quotient_v1, target_x, target_y = target_quotient.generators
    quotient_mapping = qpa.AlgebraHomomorphism.quotient_to_quotient(
        source_quotient,
        target_quotient,
        [target_quotient_v1, target_y, target_x],
        generators_with_one=True,
    )
    quotient_element = (
        quotient_v1.scale_left(field.element(2))
        + (quotient_a * quotient_b).scale_left(field.element(3))
        + (quotient_b * quotient_b).scale_left(field.element(4))
    )
    quotient_union = _as_quotient_element(quotient_element)
    quotient_mapped = quotient_mapping.mapped_expression(quotient_union)
    assert _named_terms(quotient_mapped) == fixture["quotient"]["mapped"]
    assert quotient_mapping.images_representative(quotient_union) == quotient_mapped

    with pytest.raises(qpa.InvalidInputError, match="no bound relators"):
        source_quotient.ordered_by(source_quotient.ordering)
    assert fixture["quotient"]["ordered_success"] is False

    foreign = qpa.PathAlgebra(field, qpa.Quiver(["w"], []))
    with pytest.raises(qpa.InvalidInputError, match="source algebra"):
        mapping.images_representative(_as_path_element(foreign.one()))


def _foundational_algebra(field: qpa.Field) -> tuple[qpa.PathAlgebra, qpa.AlgebraElement]:
    algebra = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [
                ("v1", "v1", "a"),
                ("v1", "v2", "b"),
                ("v2", "v3", "c"),
                ("v3", "v1", "d"),
                ("v2", "v2", "e"),
            ],
        ),
    )
    return algebra, algebra.generator("a")


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_natural_homomorphism_by_ideal_follows_qpa(field: qpa.Field) -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_foundational.schema-v1.json").read_text(encoding="utf-8")
    )
    algebra, a = _foundational_algebra(field)
    relation = a * a
    ideal = algebra.ideal([relation])
    mapping = ideal.natural_homomorphism()
    assert mapping.source.is_identical(algebra)
    assert mapping.range.original_path_algebra.is_identical(algebra)
    assert mapping.kernel.is_identical(ideal)
    assert mapping.is_surjective
    assert len(mapping.generator_images) == len(algebra.generators)

    image = mapping.image(a + relation)
    assert mapping.mapped_expression(a + relation) == image
    assert mapping.images_representative(a + relation) == image
    assert image.is_normal_form is False

    if field.characteristic is None:
        expected = fixture["natural_map"]
    elif field.characteristic == 2:
        expected = fixture["f2"]
    else:
        expected = None
    if expected is not None:
        key = "image" if field.characteristic is None else "lazy_image"
        assert _external_terms(image) == expected[key]

    attached = ideal.with_complete_groebner_basis([relation])
    attached_mapping = attached.natural_homomorphism()
    attached_image = attached_mapping.image(a + relation)
    assert attached_image.is_normal_form is True
    assert attached_image == attached_mapping.range.element(a)

    foreign, _foreign_arrow = _foundational_algebra(field)
    with pytest.raises(qpa.InvalidInputError, match="different parent path algebra"):
        mapping.image(foreign.one())


def test_quotient_uniformization_and_path_parent_aliases_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_foundational.schema-v1.json").read_text(encoding="utf-8")
    )
    algebra, a = _foundational_algebra(qpa.Field.rational())
    b = algebra.generator("b")
    quotient = algebra.ideal([a * a]).natural_homomorphism().range
    lazy_sum = quotient.element(a + b)
    assert [
        _external_terms(element)
        for element in quotient.make_uniform_on_right([lazy_sum])
    ] == fixture["natural_map"]["right_uniform"]
    assert [
        _external_terms(element)
        for element in quotient.generators_times_arrows_on_right([lazy_sum])
    ] == fixture["natural_map"]["times_arrows"]

    assert algebra.original_path_algebra.is_identical(algebra)
    assert quotient.original_path_algebra.is_identical(algebra)
    assert algebra.quiver.is_identical(quotient.quiver)
    assert a.algebra.is_identical(algebra)
    assert lazy_sum.algebra.is_identical(quotient)
    union = qpa.QuiverAlgebra.from_quotient_algebra(quotient)
    assert union.original_path_algebra.is_identical(algebra)


def test_generator_variable_assignments_follow_qpa_namespace_plan() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_mapping_ordering.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )["assignment"]
    field = qpa.Field.prime(5)
    algebra = qpa.PathAlgebra(
        field, qpa.Quiver(["v1"], [("v1", "v1", "a"), ("v1", "v1", "b")])
    )
    assignments = algebra.generator_variable_assignments()
    assert [name for name, _element in assignments] == fixture["names"]
    assert [_named_terms(element) for _name, element in assignments] == fixture["values"]
    assert fixture["success"] is True

    quotient = algebra.quotient([])
    quotient_assignments = quotient.generator_variable_assignments()
    assert [name for name, _element in quotient_assignments] == fixture["names"]
    assert all(element.quotient_element is not None for _name, element in quotient_assignments)
    union = qpa.QuiverAlgebra.from_quotient_algebra(quotient)
    assert union.generator_variable_assignments() == quotient_assignments

    invalid = qpa.PathAlgebra(
        field, qpa.Quiver(["v1"], [("v1", "v1", "not-valid")])
    )
    with pytest.raises(qpa.InvalidInputError, match="valid GAP identifier"):
        invalid.generator_variable_assignments()
    assert fixture["invalid_success"] is False

    read_only = qpa.PathAlgebra(field, qpa.Quiver(["v1"], [("v1", "v1", "Z")]))
    with pytest.raises(qpa.InvalidInputError, match="read-only"):
        read_only.generator_variable_assignments(["Z"])
    assert fixture["read_only_success"] is False
