# 0.2.1 release audit

Release preparation: 2026-10-11. Scope: ship the completed general optimization
pass as a source-compatible patch, not another optimization campaign. Publication
is a separate action: a local commit/tag does not mean PyPI or GitHub has been
updated.

## Compatibility cleanup

The source diff against `v0.2.0` was reviewed for existing public signatures,
struct fields, type aliases and error variants. No existing public signature or
field was removed or changed. Two hidden cross-crate cache helper methods on
`QuiverAlgebra` are additive; new implementation modules are private or
crate-private. Python stubs and public signatures remain unchanged.

The unpublished `ExtDimensionOverflow` and `InconsistentExtDimensions` variants
were removed. Checked arithmetic now uses the existing
`EnvelopingModuleError::Module(PathAlgebraModuleError::DimensionOverflow)` and
`InternalExtHomRejected` variants. The latter's diagnostic now covers inconsistent
Hom dimensions as well as missing Hom results. Error messages are not a stable
API contract; arithmetic guards and Python exception classes are retained.

`scripts/check-patch-api` compares all **182** `EnvelopingModuleError` variant
names/kinds with `v0.2.0` and compiles an exhaustive match as an external crate.
It passes on Rust 1.89.0. This is a focused compile regression plus source review,
not a claim that a general-purpose semver checker was run.

## Version and packaging changes

- Workspace crates, Python distribution, Cargo/uv lockfiles and the internal
  `qpa-gbnp -> qpa-core` minimum version are 0.2.1.
- `scripts/check-release-version [vX.Y.Z]` checks their agreement and rejects a
  mismatched tag. The release workflow gates artifact builds on this check.
- Linux release-wheel jobs now run the Python tests against the installed wheel
  and verify both distribution version and source commit before upload.
- Restored executable permissions on three shell test helpers and the GAP
  launcher invoked by `scripts/verify`. Their missing execute bits blocked a
  fresh full verification.
- Include frozen test fixtures in the sdist, so its Rust tests can compile and
  run without a source checkout or live upstream installation.
- Release notes and installation examples are updated. Historical benchmark
  hashes still describe the measured development builds, not the release wheel.

## Validation

Local platform: Linux x86_64, CPython 3.12; `rustc 1.89.0
(29483883e 2025-08-04)`, matching the declared minimum Rust toolchain. Default
Numerica backend and GMP are both required. Cross-platform wheel builds and
other Python versions remain the responsibility of CI; local results do not
claim those jobs have run.

All four verification modes pass on the release-preparation tree (run as
individual modes for the final checks, with the live oracle running separately):

| Gate | Result |
|---|---|
| Default Rust workspace tests | 451 passed, plus 2 doctests |
| GMP Rust workspace tests | 451 passed, plus 2 doctests |
| Strict workspace Clippy and rustdoc | Pass on both backends |
| Formatting, source ledger, Python ledger/generated docs | Pass |
| Python editable installation | 611 passed |
| Isolated release-profile candidate wheel, outside checkout | 611 passed |
| Unpacked sdist, no Git/upstream checkout | 450 Rust tests + 2 doctests passed |
| Pinned GAP/QPA/GBNP versions, full fixture contract and upstream reproducers | Pass (oracle exit 0) |
| External 0.2.0 exhaustive error match | 182 variants; pass |
| Version consistency and mismatched-tag rejection | Pass |
| Release workflow YAML and embedded shell syntax | Pass |

The workspace count includes the binding-ledger unit test; that maintainer crate
is not part of the Python sdist. The candidate wheel correctly identifies the
preparation tree as `6207cd941fefcadab5a695001be3aac573b1dacf-dirty`, not as a
clean release commit. Final clean-commit artifacts are regenerated separately.

Initial attempts exposed the helper/launcher permissions, a not-yet-created
documentation link, and the sdist's missing fixtures (87 Rust compile errors).
These were repaired and their checks rerun successfully; failed attempts are
not counted as passing gates. No test or relation validation was disabled.

Commands, from the repository root:

```sh
scripts/check-release-version v0.2.1
scripts/check-patch-api
# Equivalent to scripts/verify all; the oracle may run independently.
for mode in rust gmp python oracle; do
  PYO3_PYTHON=/usr/bin/python3 VIRTUAL_ENV="$PWD/.venv" \
    RUSTUP_TOOLCHAIN=1.89.0 scripts/verify "$mode"
done
PYO3_PYTHON=/usr/bin/python3 cargo test --workspace --doc
PYO3_PYTHON=/usr/bin/python3 cargo test --workspace --doc \
  --no-default-features --features numerica-gmp
.venv/bin/maturin build --release --locked --interpreter .venv/bin/python \
  --out target/release-0.2.1/candidate
.venv/bin/maturin sdist --out target/release-0.2.1/candidate

ROOT=$PWD
CHECK_ROOT=$(mktemp -d)
uv venv --python .venv/bin/python "$CHECK_ROOT/venv"
uv pip install --python "$CHECK_ROOT/venv/bin/python" pytest \
  target/release-0.2.1/candidate/*.whl
(cd "$CHECK_ROOT" && "$CHECK_ROOT/venv/bin/python" -I -m pytest -q "$ROOT/python/tests")
tar -xzf target/release-0.2.1/candidate/qpa_rust-0.2.1.tar.gz -C "$CHECK_ROOT"
(cd "$CHECK_ROOT/qpa_rust-0.2.1" && PYO3_PYTHON=/usr/bin/python3 \
  RUSTUP_TOOLCHAIN=1.89.0 CARGO_TARGET_DIR="$ROOT/target" cargo test --workspace)
```

The live oracle requires the pinned references in `porting/upstream-lock.toml`.
For this run, a local copy of an existing pinned reference cache was restored
under ignored `upstream/`; no neighboring research repository was modified.

## Final artifact procedure

After the release commit, regenerate the wheel and sdist from the clean checkout
into `dist/`, test the installed wheel in an isolated environment outside the
checkout, and verify `importlib.metadata.version("qpa-rust") == "0.2.1"` and
`qpa.PORT_GIT_COMMIT == git rev-parse HEAD`. Record final artifact SHA-256 hashes
outside the source commit (alongside the artifacts) to avoid self-referential
commit/hash records. Do not publish the dirty-tree candidate artifacts.

The local release commit/tag uses the repository's existing automated commit
identity (`Codex <codex@local>`) per command, without changing global Git config.
Pushing `v0.2.1` starts the release workflow; PyPI publication additionally
requires the configured trusted publisher and `PUBLISH_TO_PYPI=true` repository
variable. Neither a remote tag push nor publication is part of local preparation.

## Performance evidence and remaining limitations

See [the completed optimization report](../../benchmarks/results/native-performance.md)
and [the performance plan](performance-plan.md). The four-priority large workflow
fell from 540.65 to 180.29 seconds and from 4.36 to 2.73 GiB peak RSS. Priorities
2–4 reduce runtime but do not materially reduce the peak beyond priority 1.
Several tiny sparse queries retain overhead; rational coefficient growth and
retained resolution/Hom data remain future work, not release blockers. No further
algorithm tuning is included in release preparation.
