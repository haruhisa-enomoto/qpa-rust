# Identity, caches, and randomness

QPA object identity is observable, so the Rust model distinguishes identity
from structural equality.

## Handles and cloning

Major quiver, algebra, ideal, quotient, module, and map objects are cloneable
handles backed by shared ownership. Cloning a handle preserves QPA identity and
shares its attribute caches. Calling a constructor again with equal data makes
a new family/object, as GAP does.

Use an `is_identical` method when the source question is GAP's `IsIdenticalObj`.
Where `PartialEq` is implemented, `==` follows the documented structural or
family-aware QPA equality rule; it must not be treated as a universal identity
test.

Many QPA attributes are lazy and identity-sensitive. The port uses weak links
where reciprocal algebra/module/map caches would otherwise form ownership
cycles. As a result, a cached object is identity-stable while a live handle is
retained, without intentionally leaking the entire object graph.

## Randomized algorithms

GAP normally consumes process-global randomness. Rust makes this dependency
explicit where practical:

- `_with_rng` variants accept a caller-supplied `rand::Rng`;
- convenience variants for idempotents and module decomposition use documented
  fixed seeds, making Rust runs reproducible;
- `random_module` requires an explicit RNG.

The Rust random engine's byte stream is not promised to match GAP's global
Mersenne Twister. The translated sampling distribution, iteration order, and
QPA stopping rules are preserved. Use a seeded caller-owned RNG when an
application needs repeatable experiments, and record its algorithm and seed in
scientific output.
