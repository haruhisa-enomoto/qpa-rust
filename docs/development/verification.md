# Verification

Run all checks from the workspace root through:

```console
scripts/verify all
```

The command has smaller modes for development:

```console
scripts/verify rust
scripts/verify gmp
scripts/verify python
scripts/verify oracle
```

`rust` checks formatting, default-backend tests, strict Clippy, strict rustdoc,
and regeneration of the complete QPA operation ledger. `gmp` repeats tests and
linting with Numerica's GMP backend. `python` rebuilds the PyO3 architecture
spike into `.venv` and runs its Python tests. `oracle` validates the pinned upstream
versions, the complete checked-in GAP/QPA/GBNP fixture contract, and the
prepared compatibility reproducers. `all` runs all four modes.

The oracle mode needs the pinned GAP checkout to have been built and needs
`jq`. Python mode needs the `.venv` described in
`docs/python-architecture-spike.md`. The Rust and GMP modes do not execute GAP. A distributable Rust source
package must pass Rust tests in a clean environment where `upstream/` is not
present; the source-ledger regeneration and oracle modes are maintainer checks.

The frozen reference counts and revisions live in
`porting/rust-reference.toml`. Passing a narrower unit-test selection is not
evidence for preservation of the complete port.

The behavior-preserving source split and its unchanged-namespace rationale are
documented in [`refactor-layout.md`](refactor-layout.md).
The evidence for completion of the five post-port campaigns is recorded in the
[`completion audit`](completion-audit.md).
