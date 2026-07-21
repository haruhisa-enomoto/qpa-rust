# Identity, caches, and randomness

QPA object identity is observable, so both public interfaces distinguish
identity from structural equality.

## Handles and cloning

Major quiver, algebra, ideal, quotient, module, and map objects are cloneable
handles backed by shared ownership. Cloning a handle preserves QPA identity and
shares its attribute caches. Calling a constructor again with equal data makes
a new family/object, as GAP does.

Use `is_identical` when the source question is GAP's `IsIdenticalObj`. In
Python, `==` follows the wrapper's documented structural or family-aware QPA
equality rule. In Rust, the corresponding behavior is implemented through
`PartialEq` where available. Neither language's `==` is a universal identity
test.

Many QPA attributes are lazy and identity-sensitive. The port uses weak links
where reciprocal algebra/module/map caches would otherwise form ownership
cycles. As a result, a cached object is identity-stable while a live handle is
retained, without intentionally leaking the entire object graph.

## Randomized algorithms

GAP normally consumes process-global randomness. The port makes this
dependency explicit or deterministic where practical:

- `_with_rng` variants accept a caller-supplied `rand::Rng`;
- convenience variants for idempotents and module decomposition use documented
  fixed seeds, making Rust runs reproducible;
- Rust `random_module` requires an explicit RNG; the Python
  `QuiverAlgebra.random_module(..., seed=...)` wrapper accepts an explicit seed.

The Rust random engine's byte stream is not promised to match GAP's global
Mersenne Twister. The translated sampling distribution, iteration order, and
QPA stopping rules are preserved. Use a seeded caller-owned RNG when an
application needs repeatable experiments, and record its algorithm and seed in
scientific output.
