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

The initial exact coefficient contract uses Numerica and admits the rational
field and supported prime fields. Finite extension fields are a later
coefficient tier.

See [ROADMAP.md](ROADMAP.md) for the completed source-port campaign and
[POST_PORT_ROADMAP.md](POST_PORT_ROADMAP.md) for the completed stabilization,
refactor, and API campaigns. Maintainer verification is documented in
[docs/development/verification.md](docs/development/verification.md), and
[NOTICE.md](NOTICE.md) records source provenance and attribution. Rust users
should begin with the [supported facade guide](docs/user/README.md) and its
[quickstart](docs/user/quickstart.md).
