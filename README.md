# QPA for Rust

This repository is a source-faithful Rust port of the GAP package QPA
(Quivers and Path Algebras). QPA's GAP methods are the computational
specification; the GBNP and GAP sources are specifications wherever QPA
delegates to them.

The Rust port is complete for the initial exact coefficient contract. All
1,839 declarations and installed implementations extracted from the pinned QPA
source are mapped to verified public Rust behavior, including explicit typed
boundaries for declarations on which QPA itself installs no method. Python
bindings remain a separate follow-up campaign.

The initial exact coefficient contract uses Numerica and admits the rational
field and supported prime fields. Finite extension fields are a later
coefficient tier.

See [ROADMAP.md](ROADMAP.md) for scope, provenance, and completion criteria.
