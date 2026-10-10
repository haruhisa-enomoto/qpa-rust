"""Python power syntax preserves exact algebra arithmetic across fields."""

import pytest
import qpa


@pytest.fixture(params=["rational", "binary", "prime", "extension"])
def field(request):
    return {
        "rational": qpa.Field.rational,
        "binary": lambda: qpa.Field.prime(2),
        "prime": lambda: qpa.Field.prime(5),
        "extension": lambda: qpa.Field.extension(2, 2),
    }[request.param]()


@pytest.fixture(params=["path", "quotient", "union_path", "union_quotient"])
def elements(field, request):
    parent = qpa.PathAlgebra(field, qpa.Quiver.indexed(1, [(1, 1, "x")]))
    x = parent.generator("x")
    if "quotient" in request.param:
        algebra = parent.quotient([x * x * x])
        x = algebra.generator("x")
    else:
        algebra = parent
    one = algebra.one()
    zero = algebra.zero()
    # Include a coefficient outside the prime subfield in the GF(4) case.
    coefficient = field.primitive_element() if field.degree == 2 else field.element(1)
    value = one + x.scale_left(coefficient)
    if request.param.startswith("union"):
        wrap = (
            qpa.QuiverAlgebraElement.from_quotient_element
            if "quotient" in request.param
            else qpa.QuiverAlgebraElement.from_path_element
        )
        return tuple(wrap(element) for element in (value, one, zero, x))
    return value, one, zero, x


def test_powers_and_negation(elements):
    value, one, zero, x = elements
    assert value ** 0 == one
    assert zero ** 0 == one
    assert value ** 1 == value
    assert value ** 2 == value * value
    assert value ** 3 == value * value * value
    assert pow(value, 3, None) == value ** 3
    assert x ** 2 == x * x
    assert one ** ((1 << 64) - 1) == one
    assert value + (-value) == zero
    assert -(-value) == value
    if not isinstance(value, qpa.QuiverAlgebraElement):
        assert value ** 3 == value.pow(3)
    if isinstance(value, qpa.QuotientElement):
        inverse = value.pow_integer(-1)
        assert inverse is not None
        assert value * inverse == one


@pytest.mark.parametrize("exponent", [1.5, "2", None])
def test_noninteger_exponents_are_rejected(elements, exponent):
    with pytest.raises(TypeError):
        elements[0] ** exponent


@pytest.mark.parametrize("exponent", [-1, 1 << 64])
def test_out_of_range_exponents_are_rejected(elements, exponent):
    with pytest.raises(OverflowError):
        elements[0] ** exponent


def test_modular_powers_are_rejected(elements):
    with pytest.raises(TypeError, match="modular powers"):
        pow(elements[0], 2, 3)
