# Post-port completion audit

Audit date: 2026-07-19

The five campaigns in [`POST_PORT_ROADMAP.md`](../../POST_PORT_ROADMAP.md)
have passed their exit gates. This audit records the evidence available in the
repository; it does not move the frozen source-faithful reference tag.

## Reference baseline

- Rust reference commit: `6d5a973`
- Annotated tag: `rust-port-reference-0.1.0`
- QPA revision: `8d5f847258d4853700ce3a5d947db2a1441109f1`
- GBNP revision: `e0b46c1aac16f3fbabe743e3169dc6bd86bdc654`
- GAP revision: `4efe877859ab0bc0ac041bf635a996fb897ee0c9`
- Operation-ledger entries: 1,839

The machine-readable record is
[`porting/rust-reference.toml`](../../porting/rust-reference.toml), and the
upstream lock is [`porting/upstream-lock.toml`](../../porting/upstream-lock.toml).

## Verification matrix

| Gate | Result | Scope |
| --- | --- | --- |
| `scripts/verify rust` | Pass | Formatting, all default-backend targets, strict Clippy, strict rustdoc, and the complete operation ledger |
| `scripts/verify gmp` | Pass | All workspace targets, strict Clippy, and strict rustdoc with Numerica's GMP backend |
| `scripts/verify python` | Pass | maturin development build and 7 Python integration tests |
| `scripts/verify oracle` | Pass | Pinned GAP/QPA/GBNP metadata, fixture contract, and all compatibility reproducers |
| Clean source-tree audit | Pass | All default-backend workspace targets from a fresh copy with `.git`, `target`, `.venv`, and `upstream` absent |

The clean audit exposed and corrected a test-only dependency on the presence of
`upstream/qpa`. The affected documentation-plan test was then checked under the
default and GMP backends, followed by strict Clippy and rustdoc checks. No
production algorithm changed.

## Delivered boundaries

- The source-faithful root APIs remain available for compatibility and audit.
- `qpa_core::api` and `qpa_gbnp::api` are the supported Rust facades.
- The Python package is an architecture spike over a foundational vertical
  slice, not yet a complete binding of the 1,839-entry operation ledger.
- Upstream issue reports are local review drafts and have not been submitted.
- Full Python coverage, wheel production, performance work, finite extension
  fields, and upstream-version adoption remain later campaigns.
