# Coefficient domains

The exact contract includes rational coefficients, Numerica's specialized
prime fields, and one runtime-configured implementation of canonical finite
extension fields.

| Mathematical field | Python value | Rust value | Notes |
| --- | --- | --- | --- |
| `Q` | `qpa.Field.rational()` | `numerica::domains::rational::Q` | Exact arbitrary-precision rational arithmetic |
| `F_2` | `qpa.Field.prime(2)` | `numerica::domains::finite_field::Z2` | Specialized characteristic-two backend |
| `F_p`, odd `p` | `qpa.Field.prime(p)` | `Zp::new(p)` | `p` must be prime |
| `GF(p^n)`, `n > 1` | `qpa.Field.extension(p, n)` | `GaloisField::conway(p, n)?` | GAP Conway polynomial and designated generator |

The Python constructor checks primality and selects the binary backend for
characteristic two. Rust values are re-exported by
`qpa_core::api::coefficients` and its prelude. The lower-level Rust `Zp::new`
contract trusts the caller's primality assertion; passing an odd composite
does not describe a field and is outside the supported QPA domain. Use `Z2`
rather than `Zp::new(2)`.

## Exact scalar conversion

Python constructs exact scalars through their owning field:

```python
import qpa

rationals = qpa.Field.rational()
minus_three = rationals.element(-3)
two_thirds = rationals.element(2, 3)

field5 = qpa.Field.prime(5)
four_in_f5 = field5.element(9)

field16 = qpa.Field.extension(2, 4)
z16 = field16.primitive_element()
z16_squared = field16.frobenius(z16, 1)

field4 = qpa.Field.extension(2, 2)
z4_in_f16 = field16.embed(field4.primitive_element())
```

`Scalar.numerator` and `Scalar.denominator` are decimal strings for rational
values, preserving arbitrary precision. `Scalar.representative` is the
canonical nonnegative residue for a prime-field value.
`Scalar.coefficients` gives extension elements in the polynomial basis
`1, z, ..., z^(n-1)`, constant coefficient first. A scalar retains its field,
and mixed-field operations are rejected unless `Field.embed` requests the
canonical Conway-subfield embedding explicitly.

Scalars support `+`, `-`, `*`, `/`, unary `-`, and integer powers. Binary
arithmetic accepts scalars from the same field and exact Python integers on
either side. Integers are interpreted in the scalar's field:

```python
assert two_thirds + 1 == rationals.element(5, 3)
assert 2 / two_thirds == rationals.element(3)
assert field5.element(3) / 2 == field5.element(4)
assert z16 ** -1 * z16 == field16.element(1)
```

Division by zero in the coefficient field and negative powers of zero raise
`ZeroDivisionError`. This includes division by the characteristic, such as
`field5.element(1) / 5`. Exponents must be Python integers from `-4294967295`
through `4294967295`, inclusive, for every coefficient field. Larger magnitudes
raise `OverflowError` before arithmetic, including powers of zero or one.
This bound follows the rational coefficient backend's exponent limit; powers
within it can still produce very large exact values. Zero to the zeroth power
is the field's one. Three-argument `pow` is
unsupported. Arithmetic does not implicitly convert floats, `Fraction`
objects, booleans, or scalars from another field. Construct fractional operands
with `field.element(numerator, denominator)` and use `Field.embed` for supported
field embeddings. Equality and hashing retain the owning field: a scalar does
not compare equal to a Python integer.

In Rust, import Numerica's `Ring`, `RingOps`, and `Field` traits when using
their methods. They are re-exported by the facade.

```rust
use qpa_core::api::coefficients::{GaloisField, Integer, Q, Ring, RingOps, Z2, Zp};

let minus_three = Q.nth(Integer::from(-3));
let two_thirds = Q.to_element(2.into(), 3.into(), true);
let one_in_f2 = Z2.one();
let f5 = Zp::new(5);
let four_in_f5 = f5.nth(Integer::from(4));
let f16 = GaloisField::conway(2, 4).expect("valid Conway field");
let z16 = f16.primitive_element();
```

`GaloisField::conway` first consults the complete precomputed Conway database
from the pinned GAP revision and follows GAP's deterministic construction when
that database has no entry. `GaloisField::from_primitive_polynomial` admits an
explicit primitive irreducible modulus as a noncanonical Rust construction.
The reference arithmetic stores fixed-length coefficient vectors and remains
independent of GAP at runtime.

These are exact operations; neither public interface accepts floating-point
approximations as algebra coefficients. Matrices retain a concrete field value
and operations reject matrices from a different field instance/domain.

## Cargo features

Both Rust crates default to Numerica's `numerica-no-gmp` feature. To select its
GMP backend consistently across the workspace, disable defaults and enable
`numerica-gmp`:

```console
cargo test --workspace --all-targets --no-default-features --features numerica-gmp
```

The coefficient backend supplies arithmetic and storage. QPA/GAP/GBNP source
still supplies the algorithms, ordering, dispatch, and observable edge-case
behavior.
