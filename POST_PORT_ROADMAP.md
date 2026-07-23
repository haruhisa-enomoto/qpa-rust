# Post-port roadmap

Last updated: 2026-07-21

Status: all six campaigns below have passed their exit gates. The verification
evidence is recorded in
[`docs/development/completion-audit.md`](docs/development/completion-audit.md).

The source-faithful Rust translation is frozen at the annotated tag
`rust-port-reference-0.1.0`. The historical porting phases and their source
coverage map remain in [ROADMAP.md](ROADMAP.md). This document is the active
record of the completed work following that Rust-only campaign.

The order is intentional: record the reference behavior, validate the shape of
the future binding, refactor behind a compatibility facade, and only then
stabilize the supported Rust API. No campaign below authorizes replacing a QPA
algorithm with a mathematically equivalent independent implementation.

## Campaign 1 — Reference baseline and repository hygiene

- [x] Tag the audited Rust implementation as `rust-port-reference-0.1.0`.
- [x] Record the reference revisions and verification evidence.
- [x] Provide one documented verification entry point.
- [x] Add root licensing, attribution, and provenance documentation.
- [x] Define the pinned-upstream update policy.

Exit gate: a maintainer can identify and reproduce the exact Rust/QPA/GAP/GBNP
baseline without reconstructing commands from commit history.

## Campaign 2 — QPA compatibility and upstream issue triage

- [x] Classify every intentional deviation by cause and compatibility policy.
- [x] Give every upstream bug candidate a minimal pinned GAP reproducer.
- [x] Separate QPA bugs from GAP-version interactions and invalid-input cases.
- [x] Record expected behavior, actual behavior, affected Rust behavior, and
      local test coverage.
- [x] Prepare reviewable upstream reports without submitting them automatically.

Exit gate: every claim that QPA is broken or fragile has a minimal reproducer
and an explicit confidence/status classification.

## Campaign 3 — Python architecture spike

- [x] Add a non-published `qpa-python` PyO3 crate buildable with maturin.
- [x] Bind rational and prime-field descriptors and exact scalar conversion.
- [x] Bind a small vertical slice through quivers, paths, path algebras,
      quotients, and matrix modules.
- [x] Establish Python exception, naming, indexing, identity, and ownership
      conventions.
- [x] Test Python behavior against the Rust calls used by the binding.

This campaign is an architectural experiment, not the full Python exposure
campaign. Its public surface may change before the later Python release.

Exit gate: the binding proves that the selected Rust ownership and coefficient
model works naturally from Python and records decisions needed by refactoring.

## Campaign 4 — Behavior-preserving Rust refactor

- [x] Split `qpa-core::path_module` by its QPA source/domain boundaries.
- [x] Split `qpa-gbnp::module` by its QPA source/domain boundaries.
- [x] Extract large inline test modules and retain common fixture helpers.
- [x] Preserve current root re-exports as a compatibility facade.
- [x] Keep source-method mappings adjacent to moved implementations.
- [x] Run the Rust, ledger, documentation, and applicable oracle gates after
      each structural slice.

Exit gate: no production Rust source file is an unstructured multi-domain
monolith, while reference behavior and existing public import paths remain
unchanged.

## Campaign 5 — Supported Rust API and user documentation

- [x] Mark the supported facade separately from internal implementation APIs.
- [x] Document coefficient domains and exact scalar conversion.
- [x] Document errors, unsupported-method boundaries, and resource limits.
- [x] Document indexing, path multiplication, right-module, and row-vector
      conventions.
- [x] Document equality, identity, cloning, caches, and deterministic randomness.
- [x] Provide a Rust quickstart and source-mapped QPA-to-Rust examples.

Exit gate: a Rust user can select a supported field, construct basic QPA
objects, understand failures and conventions, and write a small computation
without reading implementation modules.

## Campaign 6 — Python user documentation

- [x] Make Python the primary onboarding path and document local `uv`
      installation.
- [x] Provide task guides for algebras, modules, homological calculations, and
      Gröbner bases.
- [x] Add executable Python examples and run them in the Python verification
      gate.
- [x] Generate a complete Python signature reference from the shipped stub.
- [x] Generate a complete QPA-to-Python source and port-note index from the
      binding ledger.
- [x] Link the pinned QPA manual and document primary Python methods in runtime
      help and the type stub.
- [x] Add compile-checked examples to the supported Rust facade documentation.

Exit gate: a Python user can install the local package, perform representative
algebra and module computations, locate every Python signature and QPA source
mapping, and run all published examples as tests.

## Completed campaign — finite extension fields

The completed runtime `GF(p^n)` coefficient campaign is recorded in
[FINITE_EXTENSION_FIELDS.md](FINITE_EXTENSION_FIELDS.md). It keeps Numerica as
the generic field contract and provides a local GAP-compatible extension-field
domain for the Rust core and Python binding.

## Later campaigns

Wheel production, performance optimization, and adoption of newer upstream
revisions remain later campaigns.
