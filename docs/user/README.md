# User guide

The implementation is source-faithful: QPA and the GAP/GBNP routines it calls
remain the computational specification. Python and Rust reorganize QPA's
surface into typed APIs; they do not select replacement algorithms.

## Python: the primary interactive interface

- [Install with `uv`](python-installation.md).
- Follow the [Python quickstart](python-quickstart.md).
- Construct [quivers, path algebras, and quotients](python-algebras.md).
- Work with [right matrix modules and homomorphisms](python-modules.md).
- Compute [syzygies, Ext, and projective resolutions](python-homological.md).
- Use the translated [Gröbner-basis routes](python-groebner.md).
- Translate names with the [GAP/QPA to Python map](gap-to-python.md).
- Navigate the [Python API and provenance references](python-api.md).
- Read the [changelog](../../CHANGELOG.md) for release changes and compatibility notes.

## Shared mathematical contract

- [Coefficient domains](coefficient-domains.md) describes exact rationals and
  prime fields in both languages.
- [Conventions](conventions.md) covers indices, path order, right modules, and
  exact output order.
- [Errors, sentinels, and limits](errors-and-limits.md) explains Python
  exceptions, Rust result types, QPA `fail`, and expensive computations.
- [Identity, caches, and randomness](identity-and-randomness.md) records the
  observable QPA object model.
- [Pinned QPA manual chapters](qpa-source-manual.md) provide the original
  mathematical descriptions for uncommon operations.

## Rust: the lower-level library interface

The supported application-facing entry points are `qpa_core::api` and
`qpa_gbnp::api`. The crate roots retain the complete translated QPA surface for
compatibility and source auditing; an item being public at the crate root does
not by itself place it in the supported facade.

- Start with the [Rust quickstart](quickstart.md).
- Use the [GAP/QPA to Rust map](gap-to-rust.md) when porting an existing script.
- Generate local API documentation with
  `cargo doc --workspace --no-deps --open`.

Maintainers should use the separate [development documentation](../development/verification.md).
