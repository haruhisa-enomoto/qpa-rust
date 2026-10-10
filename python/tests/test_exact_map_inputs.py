"""Ordinary exact Python coefficients construct the same checked maps."""
from fractions import Fraction

import pytest
import qpa


@pytest.mark.parametrize('field', [qpa.Field.rational(), qpa.Field.prime(2),
                                  qpa.Field.prime(5), qpa.Field.extension(2, 2)])
def test_exact_map_inputs_match_explicit_scalars_and_preserve_checks(field):
    algebra = qpa.PathAlgebra(field, qpa.Quiver.indexed(2, [(1, 2, 'a')]))
    module = algebra.module_from_integer_matrices([2, 2], [('a', [[1, 0], [0, 1]])])
    denominator = 3
    numerator = 2 ** 80 + 1
    literal = [[Fraction(numerator, denominator), 0], [1, -1]]
    scalar = [[field.element(numerator, denominator), field.element(0)],
              [field.element(1), field.element(-1)]]
    assert module.homomorphism_to(module, [literal, literal]) == module.homomorphism_to(module, [scalar, scalar])
    identity = [[1, 0], [0, 1]]
    assert module.homomorphism_to(module, [identity, identity]) == module.identity()
    with pytest.raises(qpa.InvalidInputError):
        module.homomorphism_to(module, [identity, [[0, 0], [0, 0]]])
    with pytest.raises(TypeError):
        module.homomorphism_to(module, [[[1.0, 0], [0, 1]], identity])
    other_field = qpa.Field.prime(7)
    with pytest.raises(qpa.InvalidInputError):
        module.homomorphism_to(module, [[[other_field.element(1), 0], [0, 1]], identity])


def test_fraction_denominator_must_be_invertible_in_map_field():
    algebra = qpa.PathAlgebra(qpa.Field.prime(5), qpa.Quiver.indexed(1, []))
    module = algebra.simple_modules()[0]
    with pytest.raises(qpa.InvalidInputError):
        module.homomorphism_to(module, [[[Fraction(1, 5)]]])
