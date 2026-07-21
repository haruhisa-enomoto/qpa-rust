from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _value(scalar: qpa.Scalar | int) -> Fraction | int:
    if isinstance(scalar, int):
        return scalar
    representative = scalar.representative
    if representative is not None:
        return representative
    assert scalar.numerator is not None and scalar.denominator is not None
    return Fraction(int(scalar.numerator), int(scalar.denominator))


def _matrix_values(
    matrices: list[list[list[qpa.Scalar]] | int],
) -> list[list[list[Fraction | int]] | int]:
    return [
        matrix
        if isinstance(matrix, int)
        else [[_value(entry) for entry in row] for row in matrix]
        for matrix in matrices
    ]


def _element_coefficients(
    element: qpa.QuiverAlgebraElement,
    basis_paths: list[list[int]],
    vertex_count: int,
    prime: int | None,
) -> list[Fraction | int]:
    positions = {tuple(path): index for index, path in enumerate(basis_paths)}
    coefficients: list[Fraction | int] = [0] * len(basis_paths)
    for path, coefficient in element.terms:
        coefficients[positions[tuple(path.external_representation)]] = _value(coefficient)
    canonical_vertices = coefficients[:vertex_count]
    for index, coefficient in enumerate(canonical_vertices):
        previous = canonical_vertices[index - 1] if index else 0
        value = coefficient - previous
        coefficients[index] = value if prime is None else int(value) % prime
    return coefficients


def _gap_coefficients(values: list[list[str]], prime: int | None) -> list[list[Any]]:
    if prime is None:
        return [[Fraction(value) for value in row] for row in values]
    return [
        [0 if value.startswith("0*") or value == "0" else 1 for value in row]
        for row in values
    ]


def _annihilator_coefficients(
    module: qpa.MatrixModule,
    basis_paths: list[list[int]],
    vertex_count: int,
    prime: int | None,
) -> list[list[Fraction | int]]:
    return [
        _element_coefficients(element, basis_paths, vertex_count, prime)
        for element in module.annihilator()
    ]


@pytest.mark.parametrize(
    ("field", "fixture_key", "prime"),
    [
        (qpa.Field.rational(), "rational", None),
        (qpa.Field.prime(2), "binary", 2),
    ],
)
def test_annihilator_of_module_follows_gap_qpa_oracle(
    field: qpa.Field, fixture_key: str, prime: int | None
) -> None:
    expected = json.loads(
        (FIXTURES / "module_annihilator.schema-v1.json").read_text(encoding="utf-8")
    )[fixture_key]

    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    simples = path.simple_modules()
    projectives = path.indecomposable_projectives()
    faithful = qpa.MatrixModule.direct_sum(projectives)
    assert faithful is not None
    path_modules = {
        "zero": path.zero_module(),
        "simple1": simples[0],
        "simple2": simples[1],
        "arrow_zero": path.module_from_integer_matrices([1, 1], []),
        "faithful": faithful,
    }
    for name, module in path_modules.items():
        assert _annihilator_coefficients(
            module, [[1], [2], [3]], 2, prime
        ) == _gap_coefficients(
            expected["path"][name], prime
        )

    loop_parent = qpa.PathAlgebra(
        field, qpa.Quiver(["w"], [("w", "w", "x")])
    )
    x = loop_parent.generator("x")
    quotient = loop_parent.quotient_high_level([x * x])
    quotient_modules = {
        "projective": quotient.indecomposable_projectives()[0],
        "simple": quotient.simple_modules()[0],
        "zero": quotient.zero_module(),
        "two_zero": quotient.module_from_integer_matrices([2], []),
    }
    for name, module in quotient_modules.items():
        assert _annihilator_coefficients(
            module, [[1], [2]], 1, prime
        ) == _gap_coefficients(
            expected["quotient"][name], prime
        )

    bound_parent = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["u1", "u2", "u3"],
            [("u1", "u2", "c"), ("u2", "u3", "d")],
        ),
    )
    bound = bound_parent.quotient_high_level(
        [bound_parent.generator("c") * bound_parent.generator("d")]
    )
    bound_modules = {
        "zero": bound.zero_module(),
        "simple3": bound.simple_modules()[2],
    }
    for name, module in bound_modules.items():
        assert _annihilator_coefficients(
            module, [[1], [2], [3], [4], [5]], 3, prime
        ) == _gap_coefficients(
            expected["bound_quotient"][name], prime
        )


@pytest.mark.parametrize(
    ("field", "fixture_key"),
    [
        (qpa.Field.rational(), "rational"),
        (qpa.Field.prime(2), "binary"),
    ],
)
def test_dimension_order_and_module_series_follow_gap_qpa_oracle(
    field: qpa.Field, fixture_key: str
) -> None:
    expected = json.loads(
        (FIXTURES / "module_series.schema-v1.json").read_text(encoding="utf-8")
    )[fixture_key]
    path = qpa.PathAlgebra(
        field,
        qpa.Quiver(
            ["v1", "v2", "v3"],
            [("v1", "v2", "a"), ("v2", "v3", "b")],
        ),
    )
    projectives = path.indecomposable_projectives()
    path_simple = path.simple_modules()[0]
    zero = path.zero_module()
    semisimple = path.module_from_integer_matrices([1, 1, 1], [])

    assert [module.dimension for module in [*projectives, zero, semisimple]] == expected[
        "path"
    ]["dimensions"]
    assert [
        projectives[2].dimension_vector_leq(projectives[0]),
        projectives[0].dimension_vector_leq(projectives[1]),
        projectives[0].dimension_vector_leq(projectives[0]),
        projectives[1].dimension_vector_leq(projectives[0]),
    ] == expected["path"]["partial_order"]
    assert [module.radical_series() for module in [*projectives, semisimple]] == expected[
        "path"
    ]["radical_series"]
    assert [module.socle_series() for module in [*projectives, semisimple]] == expected[
        "path"
    ]["socle_series"]
    assert zero.radical_series() == expected["path"]["zero_radical_series"]
    assert zero.socle_series() == expected["path"]["zero_socle_series"]

    loop_parent = qpa.PathAlgebra(
        field, qpa.Quiver(["w"], [("w", "w", "x")])
    )
    x = loop_parent.generator("x")
    quotient = loop_parent.quotient_high_level([x * x * x])
    quotient_projective = quotient.indecomposable_projectives()[0]
    quotient_simple = quotient.simple_modules()[0]
    quotient_zero = quotient.zero_module()

    assert [
        quotient_projective.dimension,
        quotient_simple.dimension,
        quotient_zero.dimension,
    ] == expected["quotient"]["dimensions"]
    assert [
        quotient_simple.dimension_vector_leq(quotient_projective),
        quotient_projective.dimension_vector_leq(quotient_simple),
    ] == expected["quotient"]["partial_order"]
    assert [
        quotient_simple.dimension_vector_leq(path_simple),
        path_simple.dimension_vector_leq(quotient_simple),
        quotient_simple.dimension_vector_leq(projectives[0]),
        projectives[0].dimension_vector_leq(quotient_simple),
        quotient_simple.dimension_vector_leq(projectives[1]),
        projectives[1].dimension_vector_leq(quotient_simple),
    ] == expected["quotient"]["cross_length_partial_order"]
    assert quotient_projective.radical_series() == expected["quotient"][
        "radical_projective"
    ]
    assert quotient_projective.socle_series() == expected["quotient"][
        "socle_projective"
    ]
    assert quotient_simple.radical_series() == expected["quotient"]["radical_simple"]
    assert quotient_simple.socle_series() == expected["quotient"]["socle_simple"]
    assert quotient_zero.radical_series() == expected["quotient"][
        "zero_radical_series"
    ]
    assert quotient_zero.socle_series() == expected["quotient"]["zero_socle_series"]


def test_module_raw_elements_matrices_and_map_algebra_preserve_qpa_boundaries() -> None:
    rationals = qpa.Field.rational()
    field5 = qpa.Field.prime(5)
    assert qpa.is_in_full_matrix_ring(
        [[rationals.element(1), rationals.element(1, 2)]], rationals
    )
    assert not qpa.is_in_full_matrix_ring([[field5.element(1)]], rationals)

    path = qpa.PathAlgebra(
        rationals, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    module = path.module_from_integer_matrices([2, 1], [("a", [[1], [-2]])])
    raw = module.element_from_external_representation(
        [[rationals.element(2), rationals.element(-3)], [rationals.element(5)]]
    )
    assert [[_value(entry) for entry in row] for row in raw.external_representation] == [
        [2, -3],
        [5],
    ]
    assert _matrix_values(module.matrices_of_path_algebra_module()) == [[[1], [-2]]]

    zero_endpoint = path.module_from_integer_matrices([1, 0], [])
    zero_row = zero_endpoint.element_from_external_representation(
        [[rationals.element(1)], []]
    )
    assert [[_value(entry) for entry in row] for row in zero_row.external_representation] == [
        [1],
        [],
    ]
    with pytest.raises(qpa.InvalidInputError, match="vertex rows"):
        module.element_from_external_representation([[rationals.element(1)]])
    with pytest.raises(qpa.InvalidInputError, match="row at vertex 1"):
        module.element_from_external_representation(
            [[rationals.element(1)], [rationals.element(0)]]
        )

    zero = path.zero_module()
    assert zero.matrices_of_path_algebra_module() == [0]
    identity = module.identity()
    assert identity.path_algebra.is_identical(qpa.QuiverAlgebra.from_path_algebra(path))

    # QPA declares LiftIdempotents for matrix modules but installs no method.
    with pytest.raises(qpa.InvalidInputError, match="installs no method"):
        module.lift_idempotents()

    # DimensionVectorPartialOrder has no common-field requirement in QPA.
    foreign = qpa.PathAlgebra(
        qpa.Field.prime(2), qpa.Quiver(["w1", "w2"], [])
    ).zero_module()
    assert foreign.dimension_vector_leq(module)


def test_module_foundation_operations_support_other_prime_fields() -> None:
    field = qpa.Field.prime(5)
    path = qpa.PathAlgebra(
        field, qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    )
    projective = path.indecomposable_projectives()[0]
    assert projective.radical_series() == [[1, 0], [0, 1]]
    assert projective.socle_series() == [[1, 0], [0, 1]]
    assert _matrix_values(projective.matrices_of_path_algebra_module()) == [[[1]]]
    assert path.zero_module().annihilator()
