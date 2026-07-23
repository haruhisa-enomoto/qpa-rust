# QPA for Rust

This repository is a source-faithful Rust port of the GAP package QPA
(Quivers and Path Algebras). QPA's GAP methods are the computational
specification; the GBNP and GAP sources are specifications wherever QPA
delegates to them.

The Rust port is complete for the initial exact coefficient contract. All
1,839 declarations and installed implementations extracted from the pinned QPA
source are mapped to verified public Rust behavior, including explicit typed
boundaries for declarations on which QPA itself installs no method. The
PyO3/maturin binding exposes the complete audited Rust surface. Its generated
ledger groups all 1,839 source entries into 888 QPA names: 883 are callable
from Python, two are represented directly by Python types, and three are
explicit unavailable boundaries for which pinned QPA itself installs no
callable method. No operation remains planned or partially exposed.

The exact coefficient contract uses Numerica for rational and prime fields and
a local Numerica-compatible `GaloisField` for `GF(p^n)`. Canonical extensions
use the pinned GAP Conway data and construction semantics; QPA remains the
algorithmic specification above the coefficient layer.

## Python quick start

The distribution is named `qpa-rust` and the import package is named `qpa`.
From another project on this machine, add the checkout with `uv`:

```console
uv add --editable /mnt/c/ai/qpa_rust_python
```

Then run an exact computation through the Rust extension:

```python
import qpa

field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v"], [("v", "v", "x")])
algebra = qpa.QuotientAlgebra.truncated(field, quiver, 3)

assert algebra.dimension() == "3"
assert algebra.generator("x").pow(3).is_zero
```

Begin with the [Python installation guide](docs/user/python-installation.md)
and [Python quickstart](docs/user/python-quickstart.md). The examples under
[`docs/examples/python`](docs/examples/python) are executable and tested.

See [ROADMAP.md](ROADMAP.md) for the completed source-port campaign and
[POST_PORT_ROADMAP.md](POST_PORT_ROADMAP.md) for the completed stabilization,
refactor, and API campaigns. Maintainer verification is documented in
[docs/development/verification.md](docs/development/verification.md), and
[NOTICE.md](NOTICE.md) records source provenance and attribution. The
[user-guide index](docs/user/README.md) separates shared mathematical
conventions, the primary Python workflow, and the supported Rust facades.
