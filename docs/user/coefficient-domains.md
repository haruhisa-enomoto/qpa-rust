# Coefficient domains

The initial exact contract is deliberately small: rational coefficients and
prime fields supported by the selected Numerica types. Finite extension fields
are a later coefficient tier.

| Mathematical field | Python value | Rust value | Notes |
| --- | --- | --- | --- |
| `Q` | `qpa.Field.rational()` | `numerica::domains::rational::Q` | Exact arbitrary-precision rational arithmetic |
| `F_2` | `qpa.Field.prime(2)` | `numerica::domains::finite_field::Z2` | Specialized characteristic-two backend |
| `F_p`, odd `p` | `qpa.Field.prime(p)` | `Zp::new(p)` | `p` must be prime |

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
```

`Scalar.numerator` and `Scalar.denominator` are decimal strings for rational
values, preserving arbitrary precision. `Scalar.representative` is the
canonical nonnegative residue for a prime-field value. A scalar retains its
field, and mixed-field operations are rejected.

In Rust, import Numerica's `Ring`, `RingOps`, and `Field` traits when using
their methods. They are re-exported by the facade.

```rust
use qpa_core::api::coefficients::{Integer, Q, Ring, RingOps, Z2, Zp};

let minus_three = Q.nth(Integer::from(-3));
let two_thirds = Q.to_element(2.into(), 3.into(), true);
let one_in_f2 = Z2.one();
let f5 = Zp::new(5);
let four_in_f5 = f5.nth(Integer::from(4));
```

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
