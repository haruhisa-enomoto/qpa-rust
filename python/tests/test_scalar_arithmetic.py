from fractions import Fraction
import operator

import pytest

import qpa


FIELDS = [qpa.Field.rational(), qpa.Field.prime(2), qpa.Field.prime(5),
          qpa.Field.extension(2, 2), qpa.Field.extension(3, 2)]


@pytest.mark.parametrize('field', FIELDS)
def test_scalar_field_arithmetic(field):
    one, zero = field.element(1), field.element(0)
    value = field.primitive_element() if field.degree and field.degree > 1 else one
    assert value + zero == value
    assert value - value == zero
    assert value + (-value) == zero
    assert value * one == value
    assert value / value == one
    assert value ** 0 == one
    assert zero ** 0 == one
    assert value ** 3 == value * value * value
    assert value ** -2 == one / (value * value)
    assert (value + one).field == field


def test_rational_results_and_arbitrary_precision_int_operands():
    q = qpa.Field.rational()
    x = q.element(2, 3)
    assert x + q.element(1, 6) == q.element(5, 6)
    assert x - q.element(1, 6) == q.element(1, 2)
    assert x * q.element(3, 4) == q.element(1, 2)
    assert x / q.element(4, 9) == q.element(3, 2)
    huge = 10 ** 100
    assert huge + x == x + huge == q.element(3 * huge + 2, 3)
    assert huge - x == q.element(3 * huge - 2, 3)
    assert x - huge == q.element(2 - 3 * huge, 3)
    assert huge * x == x * huge == q.element(2 * huge, 3)
    assert huge / x == q.element(3 * huge, 2)
    assert x / huge == q.element(2, 3 * huge)
    assert sum([x, x, x]) == q.element(2)


def test_finite_field_results_and_embedding_stays_explicit():
    f5 = qpa.Field.prime(5)
    assert f5.element(3) + 4 == f5.element(2)
    assert 4 - f5.element(3) == f5.element(1)
    assert f5.element(3) / 2 == f5.element(4)
    assert 2 / f5.element(3) == f5.element(4)
    f4 = qpa.Field.extension(2, 2)
    z = f4.primitive_element()
    assert z * z == z + 1
    assert z ** -1 == z + 1
    assert z ** 3 == f4.element(1)
    f16 = qpa.Field.extension(2, 4)
    assert f16.embed(z) ** 3 == f16.element(1)


@pytest.mark.parametrize('op', [operator.add, operator.sub, operator.mul, operator.truediv])
@pytest.mark.parametrize('left,right', [
    (qpa.Field.rational().element(1), qpa.Field.prime(5).element(1)),
    (qpa.Field.prime(2).element(1), qpa.Field.prime(5).element(1)),
    (qpa.Field.extension(2, 2).element(1), qpa.Field.extension(2, 4).element(1)),
    (qpa.Field.prime(2).element(1), qpa.Field.extension(2, 2).element(1)),
])
def test_mixed_fields_rejected(op, left, right):
    with pytest.raises(qpa.InvalidInputError):
        op(left, right)
    with pytest.raises(qpa.InvalidInputError):
        op(right, left)


@pytest.mark.parametrize('field', FIELDS)
def test_division_by_zero_and_negative_zero_power_raise_python_exception(field):
    one, zero = field.element(1), field.element(0)
    for calculate in [lambda: one / zero, lambda: 1 / zero,
                      lambda: one / 0, lambda: zero ** -1]:
        with pytest.raises(ZeroDivisionError):
            calculate()
    if field.characteristic is not None:
        with pytest.raises(ZeroDivisionError):
            one / field.characteristic


@pytest.mark.parametrize('bad', [0.5, Fraction(1, 2), '1', True, object()])
@pytest.mark.parametrize('op', [operator.add, operator.sub, operator.mul, operator.truediv])
def test_inexact_or_implicit_operands_rejected(bad, op):
    value = qpa.Field.rational().element(1)
    with pytest.raises(TypeError):
        op(value, bad)
    with pytest.raises(TypeError):
        op(bad, value)


def test_power_exponent_contract_and_scalar_identity():
    q = qpa.Field.rational()
    one = q.element(1)
    for bad in [0.5, Fraction(1, 2), True, q.element(2)]:
        with pytest.raises(TypeError):
            one ** bad
    with pytest.raises(TypeError, match='modular'):
        pow(one, 2, 5)
    assert one != 1
    assert hash(one) == hash(q.element(1))


@pytest.mark.parametrize('field', FIELDS)
def test_power_backend_exponent_boundary_is_checked_for_all_fields(field):
    limit = 2 ** 32 - 1
    one, zero = field.element(1), field.element(0)
    assert one ** limit == one
    assert one ** -limit == one
    assert zero ** limit == zero
    assert (-one) ** limit == -one
    assert (-one) ** -limit == -one
    # Values whose rational powers would require huge allocations must fail
    # before invoking the backend, even when the exponent fits Python/Rust ints.
    for value in [zero, one, field.element(2), field.element(2, 3)
                  if field.characteristic != 3 else field.element(2)]:
        for too_large in [2 ** 32, 2 ** 63 - 1, 2 ** 63, 10 ** 100]:
            for exponent in [too_large, -too_large]:
                with pytest.raises(OverflowError):
                    value ** exponent
