# Install the Python package with `uv`

The Python distribution is named `qpa-rust`; Python code imports it as `qpa`.
Python 3.10 and newer are supported.

## Use this checkout from another project

For a project on this machine, an editable path dependency is the shortest
development workflow:

```console
uv init --python 3.12
uv add --editable /mnt/c/ai/qpa_rust_python
uv run python -c "import qpa; print(qpa.Field.prime(5))"
```

`uv add` records both the distribution dependency and its local source:

```toml
[project]
dependencies = ["qpa-rust"]

[tool.uv.sources]
qpa-rust = { path = "/mnt/c/ai/qpa_rust_python", editable = true }
```

Editable installation means that Python source changes are visible directly.
Changes to the Rust extension still require `uv` to rebuild the package. Use a
non-editable local dependency when a consumer should not track ongoing source
changes:

```console
uv add /mnt/c/ai/qpa_rust_python
```

The first installation performs an optimized Rust build through Maturin. It
can take several minutes on a cold cache, particularly when the repository is
on a mounted Windows filesystem. Later builds reuse Cargo and `uv` caches.

## What `uv` supplies

The repository's `pyproject.toml` declares Maturin as its PEP 517 build
backend. A consumer project does not need to add Maturin or PyO3 itself. A
source build does require the Rust toolchain and native build prerequisites on
the host machine; this development machine already has them.

Plain `uv add qpa-rust` will become the portable installation command after
release wheels are published. Until then, use a local path or a future Git
source. An absolute local path is intentionally not portable to another
machine.

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
