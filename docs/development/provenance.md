# Computational provenance

## Authority order

QPA's installed GAP methods are the primary computational specification. GBNP
and GAP become specifications only where those QPA methods delegate to them or
depend on their observable primitive semantics. Numerica is an exact arithmetic
backend, not an alternative source of representation-theoretic algorithms.

For each user-visible QPA declaration or installed method, the generated ledger
in `porting/qpa-operations.toml` records its pinned source location, Rust symbol,
tests, and any intentional deviation identifier. `porting/deviations.md`
explains the cases where a direct observable translation is impossible or
deliberately repaired.

## Frozen reference

The annotated tag `rust-port-reference-0.1.0` identifies the completed initial
Rust translation. Its evidence and revisions are recorded in
`porting/rust-reference.toml`. Refactors compare against that tag and may not
silently rewrite its algorithms.

## Upstream update policy

Pinned upstream revisions are changed only in a dedicated update campaign:

1. Create a proposed lock-file change without modifying the reference tag.
2. Regenerate the complete operation ledger and review declaration/method
   additions, removals, and source movement.
3. Regenerate GAP oracle fixtures only after explaining every changed value or
   failure.
4. Update source translations method by method; do not make tests pass by
   substituting independent algorithms.
5. Re-run both Numerica backends, strict lint/documentation gates, the ledger
   check, and the full pinned oracle contract.
6. Record changes to deviations and compatibility notes before accepting the
   new baseline.

An upstream bug fix does not automatically change Rust behavior. It is adopted
through this process so that users can identify which QPA semantics a release
implements.

## Distribution boundary

The `upstream/` directory is a local development reference. Rust source
packages and future Python wheels must build without it and without a GAP
installation. The checked-in fixtures and source ledger contain the necessary
verification metadata but are not a substitute for the corresponding upstream
license notices.
