# Errors, source sentinels, and limits

The Rust API separates conditions that GAP represents with the same dynamic
value.

| Rust shape | Meaning |
| --- | --- |
| `Err(E)` | Invalid input, violated family/domain precondition, or failed exact computation |
| `Ok(None)` | The translated QPA method returned `fail`, `false`, or `TryNextMethod()` as documented for that function |
| Result enum | QPA returns heterogeneous values such as finite/infinite or value/false |
| `CachedAttribute<T>` | Distinguishes unbound, cached `fail`, cached `false`, and a cached value |

Read each operation's rustdoc: `Option` is not assigned one universal GAP
meaning. For example, the checked module constructor uses `None` for relation
failure, while some property dispatch uses it for `TryNextMethod()`.

Error enums preserve context and use transparent variants for lower-level exact
matrix, quotient, coefficient, and GBNP failures. Pattern-match a concrete
variant when recovery depends on the cause; otherwise propagate it with `?`.
The crate does not currently erase all failures into one global error type.

## Unsupported boundaries

The complete crate-root surface includes typed representations of QPA
declarations for which the pinned package installs no method, broken source
branches, and coefficient tiers outside the initial contract. Such a function
can be public for source coverage without being supported for an ordinary
successful computation. `#[doc(hidden)]` cross-crate hooks are also not part of
the supported facade.

Consult [`porting/deviations.md`](../../porting/deviations.md) for every
intentional Rust/QPA boundary and
[`porting/compatibility-classification.md`](../../porting/compatibility-classification.md)
for its classification. The code preserves or types these boundaries; it does
not silently replace a failing QPA method with another mathematical algorithm.

## Resource limits

- Exact integers and rationals are arbitrary precision, but container sizes,
  matrix dimensions, and Rust indices are bounded by `usize`.
- The supported odd-prime `Zp` backend uses a `u32` characteristic. Algorithms
  that enumerate every field element can additionally fail when the cardinality
  cannot be addressed or allocated.
- Gröbner bases, nontip bases, submodule enumeration, and decomposition can grow
  exponentially or worse, just as in QPA/GBNP. A finite mathematical result is
  not a memory or runtime guarantee.
- Source methods known to loop outside their intended domain have explicit
  precondition guards where the compatibility ledger records one.
- Finite extension fields and coefficient-dependent central splitting over
  unsupported domains remain later work.
