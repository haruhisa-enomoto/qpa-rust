from __future__ import annotations

import json
from pathlib import Path

import pytest

import qpa


FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "tests"
    / "fixtures"
    / "idempotent.schema-v1.json"
)


def _zero_tensor(field: qpa.Field, dimension: int) -> list[list[list[qpa.Scalar]]]:
    return [
        [
            [field.element(0) for _ in range(dimension)]
            for _ in range(dimension)
        ]
        for _ in range(dimension)
    ]


def _full_matrix_algebra(field: qpa.Field, size: int) -> qpa.StructureConstantsAlgebra:
    dimension = size * size
    products = _zero_tensor(field, dimension)
    for left_row in range(size):
        for inner in range(size):
            left = left_row * size + inner
            for right_column in range(size):
                right = inner * size + right_column
                output = left_row * size + right_column
                products[left][right][output] = field.element(1)
    return qpa.StructureConstantsAlgebra(field, products)


def _product_of_fields(field: qpa.Field, count: int) -> qpa.StructureConstantsAlgebra:
    products = _zero_tensor(field, count)
    for position in range(count):
        products[position][position][position] = field.element(1)
    return qpa.StructureConstantsAlgebra(field, products)


def _representatives(values: list[qpa.Scalar]) -> list[int]:
    return [value.representative for value in values]


def _rank_mod_prime(rows: list[list[int]], prime: int) -> int:
    rows = [[entry % prime for entry in row] for row in rows]
    rank = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next(
            (position for position in range(rank, len(rows)) if rows[position][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(entry * inverse) % prime for entry in rows[rank]]
        for position, row in enumerate(rows):
            if position == rank or row[column] == 0:
                continue
            multiple = row[column]
            rows[position] = [
                (entry - multiple * pivot_entry) % prime
                for entry, pivot_entry in zip(row, rows[rank], strict=True)
            ]
        rank += 1
    return rank


@pytest.mark.parametrize(
    ("prime", "size", "fixture_key"),
    [(2, 2, "f2_m2"), (3, 2, "f3_m2"), (2, 3, "f2_m3")],
)
def test_hopf_helpers_follow_pinned_qpa_invariants(
    prime: int, size: int, fixture_key: str
) -> None:
    expected = json.loads(FIXTURE.read_text())[fixture_key]
    field = qpa.Field.prime(prime)
    algebra = _full_matrix_algebra(field, size)
    basis = algebra.basis_elements

    assert _representatives(algebra.hopf_minimum_polynomial(basis[0])) == expected[
        "e11_minpoly"
    ]
    if "irreducible_minpoly" in expected:
        irreducible = algebra.element(
            [field.element(value) for value in [0, 1, 1, 1]]
        )
        assert _representatives(
            algebra.hopf_minimum_polynomial(irreducible)
        ) == expected["irreducible_minpoly"]

    left_identity = algebra.hopf_left_identity()
    assert left_identity is not None
    assert (left_identity == algebra.find_multiplicative_identity()) is expected[
        "left_identity"
    ]

    zero_divisor = algebra.hopf_zero_divisor()
    assert zero_divisor is not None and len(zero_divisor) == 2
    assert all(element != algebra.zero() for element in zero_divisor) is expected[
        "zero_divisor"
    ]["nonzero"]
    assert (zero_divisor[0] * zero_divisor[1] == algebra.zero()) is expected[
        "zero_divisor"
    ]["product_zero"]

    singular = algebra.hopf_singular_idempotent()
    assert (singular * singular == singular) is expected["singular"]["idempotent"]
    right_rows = [
        _representatives((singular * element).coordinates) for element in basis
    ]
    assert _rank_mod_prime(right_rows, prime) == expected["singular"][
        "right_dimension"
    ]


def test_hopf_helpers_cover_every_supported_prime_field_and_source_failures() -> None:
    for prime in (2, 3, 5, 7):
        field = qpa.Field.prime(prime)
        algebra = _full_matrix_algebra(field, 2)
        assert _representatives(
            algebra.hopf_minimum_polynomial(algebra.basis_elements[0])
        ) == [0, prime - 1, 1]
        pair = algebra.hopf_zero_divisor()
        assert pair is not None and pair[0] * pair[1] == algebra.zero()
        singular = algebra.hopf_singular_idempotent()
        assert singular * singular == singular

    expected = json.loads(FIXTURE.read_text())["commutative"]
    rational = qpa.Field.rational()
    rational_product = _product_of_fields(rational, 2)
    binary_product = _product_of_fields(qpa.Field.prime(2), 2)
    assert (rational_product.hopf_zero_divisor() is None) is expected[
        "rational_zero_divisor_empty"
    ]
    assert (binary_product.hopf_zero_divisor() is None) is expected[
        "binary_zero_divisor_empty"
    ]
    assert (
        rational_product.hopf_singular_idempotent()
        == rational_product.find_multiplicative_identity()
    ) is expected["rational_singular_is_one"]

    field3 = qpa.Field.prime(3)
    zero_algebra = qpa.StructureConstantsAlgebra(field3, _zero_tensor(field3, 1))
    assert (zero_algebra.hopf_left_identity() is None) is expected[
        "zero_left_identity_is_fail"
    ]

    with pytest.raises(qpa.InvalidInputError, match="supported finite field"):
        _full_matrix_algebra(rational, 2).hopf_zero_divisor()
    with pytest.raises(qpa.InvalidInputError, match="another coefficient field"):
        _full_matrix_algebra(qpa.Field.prime(2), 2).hopf_minimum_polynomial(
            _full_matrix_algebra(qpa.Field.prime(3), 2).basis_elements[0]
        )
