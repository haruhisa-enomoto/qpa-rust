# Install the Python package with `uv`

The Python distribution is named `qpa-rust`; Python code imports it as `qpa`.
Python 3.10 and newer are supported.

This guide documents the 0.2.1 API. This patch release speeds up existing
calculations without an API migration. Once published, select it with
`uv add "qpa-rust==0.2.1"`.

## Install from PyPI

Release wheels are published on PyPI for Linux (x86_64 and aarch64), macOS
(Intel and Apple Silicon), and Windows (x86_64) for Python 3.10 through 3.13.
Installing them needs no Rust toolchain:

```console
uv init --python 3.12
uv add qpa-rust
uv run python -c "import qpa; print(qpa.Field.prime(5))"
```

Pin an exact version when a computation must stay reproducible, for example
`uv add "qpa-rust==0.2.1"`. On a platform without a prebuilt wheel, `uv`
falls back to the source distribution, which needs the Rust toolchain
described below.

## Install the development source from GitHub

To use unreleased changes, add the Git repository instead. This compiles the
Rust extension on your machine, so a Rust toolchain (the version pinned in
`rust-toolchain.toml` is installed automatically by `rustup`) and a C compiler
are required:

```console
uv add "qpa-rust @ git+https://github.com/haruhisa-enomoto/qpa-rust"
```

A tag such as `@v0.2.1` after the URL selects a published release commit. The same Git
URL works with `pip install`.

## Use a local clone from another project

For a project on the same machine, an editable path dependency is the shortest
development workflow:

```console
git clone https://github.com/haruhisa-enomoto/qpa-rust
uv add --editable /path/to/qpa-rust
```

`uv add` records both the distribution dependency and its local source:

```toml
[project]
dependencies = ["qpa-rust"]

[tool.uv.sources]
qpa-rust = { path = "/path/to/qpa-rust", editable = true }
```

Editable installation means that Python source changes are visible directly.
Changes to the Rust extension still require `uv` to rebuild the package. Use a
non-editable local dependency when a consumer should not track ongoing source
changes:

```console
uv add /path/to/qpa-rust
```

The first installation performs an optimized Rust build through Maturin. It
can take several minutes on a cold cache, particularly when the repository is
on a mounted Windows filesystem. Later builds reuse Cargo and `uv` caches.

## What `uv` supplies

The repository's `pyproject.toml` declares Maturin as its PEP 517 build
backend. A consumer project does not need to add Maturin or PyO3 itself. A
source build does require Rust **1.89 or newer** and native build prerequisites
on the host machine. Git checkouts pin Rust 1.89.0; when building a source
distribution, select a compatible toolchain yourself (for example,
`RUSTUP_TOOLCHAIN=1.89.0 uv add "qpa-rust==0.2.1"`).

An absolute local path is intentionally not portable to another machine; use
the PyPI release or the Git URL in shared project files.

## Verify an environment

This smoke test checks the native extension as well as the Python package:

```console
uv run python - <<'PY'
import qpa

field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v"], [("v", "v", "x")])
algebra = qpa.QuotientAlgebra.truncated(field, quiver, 3)
assert algebra.dimension() == "3"
print("qpa is ready")
PY
```

If `import qpa` fails but installation succeeded, check that the command is
being run through the same project with `uv run`. Remember that the installed
name is `qpa-rust`, while the imported name is `qpa`.
