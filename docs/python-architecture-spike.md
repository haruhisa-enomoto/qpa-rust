# Python architecture spike

The `qpa-python` crate is a deliberately small PyO3/maturin vertical slice. It
calls the existing `qpa-core` operations; it contains no representation-theory
algorithm of its own. The spike is non-published and its API may change before
the complete Python operation-ledger campaign.

## Decisions tested by the spike

- One Python `Field` runtime descriptor selects a monomorphized Rust backend:
  Numerica `Q`, `Z2`, or `Zp`. This avoids exposing Rust generics to Python.
- Algebra, element, quotient, and module wrappers use private Rust enums with
  the same variants. Operations dispatch only when variants match, after which
  the ordinary Rust identity and field checks still run.
- Python exact integers are converted through decimal text to Numerica's
  arbitrary-precision `Integer`. Rational numerators and denominators remain
  exact; finite-field division rejects a denominator that is zero modulo the
  characteristic.
- Python objects own cloned Rust handles. Those handles contain `Arc`-backed
  QPA identities, so cloning a wrapper preserves object identity and does not
  copy a mathematical object.
- The spike exposes `is_identical` explicitly. This keeps object identity
  separate from structural/module equality before Python `==` policy is frozen.
- Quiver paths retain QPA's left-to-right multiplication. User-visible QPA
  vertex/generator positions remain one-based; Python container positions will
  be zero-based when the full binding introduces indexed collections.
- Rust failures become a `QpaError` hierarchy. Entered-data and family errors
  use `InvalidInputError`; QPA's module-constructor `fail` result uses
  `RelationNotSatisfiedError` rather than an unrelated exception.
- The spike uses Pythonic snake-case methods while keeping QPA nouns. A future
  generated name map can add migration aliases without changing Rust symbols.

## Development build

```console
python3 -m venv .venv
.venv/bin/python -m pip install maturin pytest
.venv/bin/maturin develop
.venv/bin/python -m pytest python/tests
```

The future release may use Python's stable ABI, but the spike deliberately
builds against the active interpreter so it can validate all object behavior
before wheel policy is frozen.
