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
Prebuilt wheels are published on PyPI for Linux, macOS, and Windows on
Python 3.10 through 3.13, so no Rust toolchain is needed:

```console
uv add qpa-rust
```

or, with pip, `pip install qpa-rust`. To build the current development
source instead, add the Git repository or a local clone:

```console
uv add "qpa-rust @ git+https://github.com/haruhisa-enomoto/qpa-rust"
uv add --editable /path/to/qpa-rust
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

Begin with the [Python installation guide](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/docs/user/python-installation.md)
and [Python quickstart](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/docs/user/python-quickstart.md). The examples under
[`docs/examples/python`](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/docs/examples/python) are executable and tested.

See [ROADMAP.md](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/ROADMAP.md) for the completed source-port campaign and
[POST_PORT_ROADMAP.md](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/POST_PORT_ROADMAP.md) for the completed stabilization,
refactor, and API campaigns. Maintainer verification is documented in
[docs/development/verification.md](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/docs/development/verification.md), and
[NOTICE.md](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/NOTICE.md) records source provenance and attribution. Every
intentional behavioral difference from pinned QPA is listed in
[porting/deviations.md](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/porting/deviations.md). The
[user-guide index](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/docs/user/README.md) separates shared mathematical
conventions, the primary Python workflow, and the supported Rust facades.
