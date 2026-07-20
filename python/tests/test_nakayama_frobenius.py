from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _gap_scalar(value: qpa.Scalar) -> str:
    if value.numerator is not None:
        return repr(value)
    characteristic = value.field.characteristic
    representative = value.representative
    assert characteristic is not None and representative is not None
    if representative == 0:
        return f"0*Z({characteristic})"
    if characteristic == 2:
        return "Z(2)^0"
    if characteristic == 5:
        exponent = {1: 0, 2: 1, 4: 2, 3: 3}[representative]
        return "Z(5)" if exponent == 1 else f"Z(5)^{exponent}"
    raise AssertionError(f"test scalar formatter does not cover GF({characteristic})")


def _gap_matrix(matrix: list[list[qpa.Scalar]]) -> list[list[str]]:
    return [[_gap_scalar(value) for value in row] for row in matrix]


def _quotient_union(
    parent: qpa.PathAlgebra, relations: list[qpa.AlgebraElement]
) -> tuple[qpa.QuiverAlgebra, qpa.QuotientAlgebra]:
    quotient = parent.quotient(relations)
    return qpa.QuiverAlgebra.from_quotient_algebra(quotient), quotient


def _radical_square_zero_cycle(
    field: qpa.Field,
) -> tuple[qpa.QuiverAlgebra, qpa.QuotientAlgebra]:
    parent = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [
                ("v1", "v2", "x"),
                ("v2", "v3", "y"),
                ("v3", "v1", "z"),
            ],
        ),
    )
    x, y, z = (parent.generator(name) for name in ("x", "y", "z"))
    return _quotient_union(parent, [x * y, y * z, z * x])


def _quantum_complete_intersection(
    field: qpa.Field,
) -> tuple[qpa.QuiverAlgebra, qpa.QuotientAlgebra]:
    parent = qpa.PathAlgebra(
        field, qpa.Quiver(["v"], [("v", "v", "x"), ("v", "v", "y")])
    )
    x, y = parent.generator("x"), parent.generator("y")
    scalar = field.element(1 if field.characteristic == 2 else 2)
    return _quotient_union(parent, [x * x, y * y, x * y - (y * x).scale_left(scalar)])


def _permutation_cases(
    field: qpa.Field,
) -> dict[str, tuple[qpa.QuiverAlgebra, list[qpa.MatrixModule]]]:
    discrete_parent = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2", "v3"], []))
    discrete = qpa.QuiverAlgebra.from_path_algebra(discrete_parent)

    a2_parent = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    a2 = qpa.QuiverAlgebra.from_path_algebra(a2_parent)

    dual_parent = qpa.PathAlgebra(
        field, qpa.Quiver(["v"], [("v", "v", "t")])
    )
    t = dual_parent.generator("t")
    dual, dual_quotient = _quotient_union(dual_parent, [t * t])

    a3_parent = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [("v1", "v2", "u"), ("v2", "v3", "v")],
        ),
    )
    nonself, nonself_quotient = _quotient_union(
        a3_parent, [a3_parent.generator("u") * a3_parent.generator("v")]
    )
    cycle, cycle_quotient = _radical_square_zero_cycle(field)

    return {
        "path_discrete": (discrete, discrete_parent.simple_modules()),
        "path_a2": (a2, a2_parent.simple_modules()),
        "quotient_dual_numbers": (dual, dual_quotient.simple_modules()),
        "quotient_nonselfinjective": (nonself, nonself_quotient.simple_modules()),
        "quotient_radical_square_zero_cycle": (
            cycle,
            cycle_quotient.simple_modules(),
        ),
    }


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_nakayama_permutation_follows_qpa(field: qpa.Field) -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_nakayama_permutation.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    key = "rational" if field.characteristic is None else f"f{field.characteristic}"

    for name, (algebra, simples) in _permutation_cases(field).items():
        expected = fixture[key][name]
        permutation = algebra.nakayama_permutation()
        if expected.get("value") is False:
            assert permutation is None
            continue
        assert isinstance(permutation, qpa.NakayamaPermutation)
        assert permutation.algebra.is_identical(algebra)
        assert permutation.index_images == expected["index_images"]
        assert [permutation.image_of_index(i) for i in range(1, len(simples) + 1)] == expected[
            "index_images"
        ]
        images = [permutation.image_of_module(simple) for simple in simples]
        assert [image.dimension_vector for image in images] == expected["module_images"]
        assert [
            image.is_identical(simple) for image, simple in zip(images, simples, strict=True)
        ] == expected["module_identical"]


def test_nakayama_permutation_rejects_qpa_invalid_inputs() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_nakayama_permutation_dispatch.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    field = qpa.Field.prime(2)
    cycle, cycle_quotient = _radical_square_zero_cycle(field)
    permutation = cycle.nakayama_permutation()
    assert permutation is not None
    simples = cycle_quotient.simple_modules()
    nonsimple = qpa.MatrixModule.direct_sum(simples[:2])
    assert nonsimple is not None
    with pytest.raises(qpa.InvalidInputError):
        permutation.image_of_module(nonsimple)
    assert fixture["cycle_nonsimple_success"] is False
    for key, index in (
        ("cycle_zero_index_success", 0),
        ("cycle_negative_index_success", -1),
        ("cycle_large_index_success", 4),
    ):
        with pytest.raises(qpa.InvalidInputError):
            permutation.image_of_index(index)
        assert fixture[key] is False

    foreign_parent = qpa.PathAlgebra(field, qpa.Quiver(["w1", "w2", "w3", "w4"], []))
    with pytest.raises(qpa.InvalidInputError):
        permutation.image_of_module(foreign_parent.simple_modules()[3])
    assert fixture["foreign_simple_success"] is False

    discrete_parent = qpa.PathAlgebra(field, qpa.Quiver(["u1", "u2"], []))
    discrete = qpa.QuiverAlgebra.from_path_algebra(discrete_parent)
    discrete_permutation = discrete.nakayama_permutation()
    assert discrete_permutation is not None
    discrete_nonsimple = qpa.MatrixModule.direct_sum(discrete_parent.simple_modules())
    assert discrete_nonsimple is not None
    with pytest.raises(qpa.InvalidInputError):
        discrete_permutation.image_of_module(discrete_nonsimple)
    with pytest.raises(qpa.InvalidInputError):
        discrete_permutation.image_of_index(0)
    assert fixture["discrete_nonsimple_success"] is False
    assert fixture["discrete_zero_index_success"] is False

    arrow_parent = qpa.PathAlgebra(
        field, qpa.Quiver(["a", "b"], [("a", "b", "r")])
    )
    nonadmissible, _ = _quotient_union(arrow_parent, [arrow_parent.generator("r")])
    with pytest.raises(qpa.InvalidInputError):
        nonadmissible.nakayama_permutation()
    assert fixture["nonadmissible_success"] is False


def _automorphism_cases(
    field: qpa.Field,
) -> dict[str, tuple[qpa.QuiverAlgebra, qpa.QuotientAlgebra | None]]:
    discrete_parent = qpa.PathAlgebra(field, qpa.Quiver(["v1", "v2"], []))
    discrete = qpa.QuiverAlgebra.from_path_algebra(discrete_parent)
    a2_parent = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    a2 = qpa.QuiverAlgebra.from_path_algebra(a2_parent)

    dual_parent = qpa.PathAlgebra(
        field, qpa.Quiver(["v"], [("v", "v", "t")])
    )
    t = dual_parent.generator("t")
    dual, dual_quotient = _quotient_union(dual_parent, [t * t])

    a3_parent = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [("v1", "v2", "u"), ("v2", "v3", "v")],
        ),
    )
    nonself, nonself_quotient = _quotient_union(
        a3_parent, [a3_parent.generator("u") * a3_parent.generator("v")]
    )
    cycle, cycle_quotient = _radical_square_zero_cycle(field)
    quantum, quantum_quotient = _quantum_complete_intersection(field)
    return {
        "path_discrete": (discrete, None),
        "path_a2": (a2, None),
        "quotient_dual_numbers": (dual, dual_quotient),
        "quotient_nonselfinjective": (nonself, nonself_quotient),
        "quotient_radical_square_zero_cycle": (cycle, cycle_quotient),
        "quotient_quantum_complete_intersection": (quantum, quantum_quotient),
    }


def _canonical_basis(
    algebra: qpa.QuiverAlgebra, quotient: qpa.QuotientAlgebra | None
) -> list[qpa.QuiverAlgebraElement]:
    if quotient is not None:
        basis = quotient.canonical_basis()
        assert basis is not None
        return [qpa.QuiverAlgebraElement.from_quotient_element(value) for value in basis]
    path = algebra.path_algebra
    assert path is not None
    basis = path.canonical_basis()
    assert basis is not None
    return [qpa.QuiverAlgebraElement.from_path_element(value) for value in basis]


def _linear_combination(
    field: qpa.Field,
    basis: list[qpa.QuiverAlgebraElement],
    coefficients: list[qpa.Scalar],
) -> qpa.QuiverAlgebraElement:
    result = basis[0].scale_left(field.element(0))
    for basis_element, coefficient in zip(basis, coefficients, strict=True):
        result = result + basis_element.scale_left(coefficient)
    return result


@pytest.mark.parametrize(
    "field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)]
)
def test_nakayama_automorphism_and_frobenius_data_follow_qpa(
    field: qpa.Field,
) -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_nakayama_automorphism.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    key = "rational" if field.characteristic is None else f"f{field.characteristic}"

    for name, (algebra, quotient) in _automorphism_cases(field).items():
        expected: dict[str, Any] = fixture[key][name]
        automorphism = algebra.nakayama_automorphism()
        if expected["automorphism"] is False:
            assert automorphism is None
            if quotient is None:
                assert algebra.order_of_nakayama_automorphism() is None
            else:
                assert quotient.frobenius_form() is None
                assert quotient.frobenius_linear_functional() is None
                with pytest.raises(qpa.InvalidInputError):
                    algebra.order_of_nakayama_automorphism()
            continue

        assert isinstance(automorphism, qpa.NakayamaAutomorphism)
        assert automorphism.algebra.is_identical(algebra)
        assert _gap_matrix(automorphism.matrix) == expected["automorphism"]
        assert algebra.order_of_nakayama_automorphism() == expected["order"]

        basis = _canonical_basis(algebra, quotient)
        for index, row in enumerate(automorphism.matrix):
            expected_image = _linear_combination(field, basis, row)
            assert automorphism.image(basis[index]) == expected_image
            assert automorphism(basis[index]) == expected_image

        if quotient is None:
            continue
        form = quotient.frobenius_form()
        functional = quotient.frobenius_linear_functional()
        assert isinstance(form, qpa.FrobeniusForm)
        assert isinstance(functional, qpa.FrobeniusLinearFunctional)
        assert form.algebra.is_identical(quotient)
        assert functional.algebra.is_identical(quotient)
        assert _gap_matrix(form.matrix) == expected["form"]
        assert [_gap_scalar(value) for value in functional.canonical_basis_values] == expected[
            "functional"
        ]
        for row, left in enumerate(basis):
            assert _gap_scalar(functional.evaluate(left)) == expected["functional"][row]
            assert _gap_scalar(functional(left)) == expected["functional"][row]
            for column, right in enumerate(basis):
                assert _gap_scalar(form.evaluate(left, right)) == expected["form"][row][column]
                assert _gap_scalar(form(left, right)) == expected["form"][row][column]


def test_nakayama_frobenius_dispatch_errors_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "pathalg_nakayama_automorphism_dispatch.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )
    field = qpa.Field.prime(2)
    cycle, cycle_quotient = _radical_square_zero_cycle(field)
    automorphism = cycle.nakayama_automorphism()
    form = cycle_quotient.frobenius_form()
    functional = cycle_quotient.frobenius_linear_functional()
    assert automorphism is not None and form is not None and functional is not None
    cycle_one = qpa.QuiverAlgebraElement.from_quotient_element(cycle_quotient.one())
    foreign_parent = qpa.PathAlgebra(field, qpa.Quiver(["foreign"], []))
    foreign = qpa.QuiverAlgebraElement.from_path_element(foreign_parent.one())

    for key, operation in (
        ("foreign_auto_success", lambda: automorphism(foreign)),
        ("foreign_form_left_success", lambda: form(foreign, cycle_one)),
        ("foreign_form_right_success", lambda: form(cycle_one, foreign)),
        ("foreign_functional_success", lambda: functional(foreign)),
    ):
        with pytest.raises(qpa.InvalidInputError):
            operation()
        assert fixture[key] is False

    a3_parent = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [("v1", "v2", "u"), ("v2", "v3", "v")],
        ),
    )
    nonself, nonself_quotient = _quotient_union(
        a3_parent, [a3_parent.generator("u") * a3_parent.generator("v")]
    )
    assert nonself.nakayama_automorphism() is None
    assert nonself_quotient.frobenius_form() is None
    assert nonself_quotient.frobenius_linear_functional() is None
    with pytest.raises(qpa.InvalidInputError):
        nonself.order_of_nakayama_automorphism()
    assert fixture["nonself_auto_is_false"] is True
    assert fixture["nonself_form_is_false"] is True
    assert fixture["nonself_functional_is_false"] is True
    assert fixture["nonself_order_success"] is False

    arrow_parent = qpa.PathAlgebra(
        field, qpa.Quiver(["a", "b"], [("a", "b", "r")])
    )
    nonadmissible, nonadmissible_quotient = _quotient_union(
        arrow_parent, [arrow_parent.generator("r")]
    )
    for key, operation in (
        ("nonadmissible_auto_success", nonadmissible.nakayama_automorphism),
        ("nonadmissible_form_success", nonadmissible_quotient.frobenius_form),
        (
            "nonadmissible_functional_success",
            nonadmissible_quotient.frobenius_linear_functional,
        ),
        (
            "nonadmissible_order_success",
            nonadmissible.order_of_nakayama_automorphism,
        ),
    ):
        with pytest.raises(qpa.InvalidInputError):
            operation()
        assert fixture[key] is False
