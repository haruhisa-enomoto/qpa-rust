from __future__ import annotations

import pytest

import qpa


def test_exact_rational_and_prime_field_conversion() -> None:
    rationals = qpa.Field.rational()
    huge = 10**80 + 1
    value = rationals.element(huge, 3)
    assert value.numerator == str(huge)
    assert value.denominator == "3"
    assert repr(value) == f"{huge}/3"

    field5 = qpa.Field.prime(5)
    residue = field5.element(7, 3)
    assert residue.representative == 4
    assert repr(residue) == "4 mod 5"
    assert residue.field == field5

    with pytest.raises(qpa.InvalidInputError, match="not prime"):
        qpa.Field.prime(6)
    with pytest.raises(qpa.InvalidInputError, match=r"zero in GF\(2\)"):
        qpa.Field.prime(2).element(1, 2)


def test_quiver_identity_and_left_to_right_paths() -> None:
    quiver = qpa.Quiver(
        ["v1", "v2", "v3"],
        [("v1", "v2", "a"), ("v2", "v3", "b")],
    )
    path = quiver.path(["a", "b"])
    assert path.source == "v1"
    assert path.target == "v3"
    assert path.length == 2
    assert path.external_representation == [4, 5]
    assert quiver.arrow_path("a") * quiver.arrow_path("b") == path

    noncomposable = quiver.arrow_path("b") * quiver.arrow_path("a")
    assert noncomposable.length is None
    assert noncomposable.source is None

    reconstructed = qpa.Quiver(quiver.vertices, quiver.arrows)
    assert not quiver.is_identical(reconstructed)


def test_path_algebra_uses_exact_rust_operations() -> None:
    quiver = qpa.Quiver(
        ["v1", "v2", "v3"],
        [("v1", "v2", "a"), ("v2", "v3", "b")],
    )
    algebra = qpa.PathAlgebra(qpa.Field.rational(), quiver)
    assert algebra.quiver.is_identical(quiver)

    a = algebra.generator("a")
    b = algebra.generator("b")
    ab = a * b
    assert not ab.is_zero
    assert (b * a).is_zero

    half_a = algebra.term(quiver.arrow_path("a"), algebra.field.element(1, 2))
    assert "1/2" in repr(half_a)
    assert not (half_a + half_a).is_zero

    independent = qpa.PathAlgebra(qpa.Field.rational(), quiver)
    assert not algebra.is_identical(independent)
    with pytest.raises(qpa.InvalidInputError, match="different path algebras"):
        _ = a + independent.generator("a")


def test_quotient_runs_the_existing_high_level_groebner_route() -> None:
    quiver = qpa.Quiver(["v"], [("v", "v", "x")])
    algebra = qpa.PathAlgebra(qpa.Field.prime(5), quiver)
    x = algebra.generator("x")
    quotient = algebra.quotient([x * x])

    assert quotient.dimension() == "2"
    residue = quotient.generator("x")
    assert residue.is_normal_form
    assert (residue * residue).is_zero
    assert not quotient.one().is_zero


@pytest.mark.parametrize("field", [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5)])
def test_checked_matrix_module_calls_qpa_core(field: qpa.Field) -> None:
    quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
    algebra = qpa.PathAlgebra(field, quiver)
    module = algebra.module_from_integer_matrices([1, 1], [("a", [[1]])])

    assert module.dimension_vector == [1, 1]
    assert module.dimension == 2
    assert module.arrow_matrix_shapes == [(1, 1)]
    assert "dimension vector [ 1, 1 ]" in repr(module)
    assert module.is_identical(module)
