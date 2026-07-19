# Rust user guide

The supported application-facing entry points are `qpa_core::api` and
`qpa_gbnp::api`. The crate roots retain the complete translated QPA surface for
compatibility and source auditing; an item being public at the crate root does
not by itself place it in the supported facade.

Start with:

- [`quickstart.md`](quickstart.md) for a complete first computation;
- [`coefficient-domains.md`](coefficient-domains.md) for exact rational and
  prime-field scalars;
- [`conventions.md`](conventions.md) before translating GAP/QPA input;
- [`errors-and-limits.md`](errors-and-limits.md) for `Result`, `Option`, `fail`,
  unsupported operations, and resource limits;
- [`identity-and-randomness.md`](identity-and-randomness.md) for QPA object
  identity, caches, cloning, and reproducibility;
- [`gap-to-rust.md`](gap-to-rust.md) for operation-name and source mappings.

The implementation is source-faithful: QPA and the GAP/GBNP routines it calls
remain the computational specification. These facade modules reorganize names;
they do not select replacement algorithms.
