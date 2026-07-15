# QPA for Rust

This repository is a source-faithful Rust port of the GAP package QPA
(Quivers and Path Algebras). QPA's GAP methods are the computational
specification; the GBNP and GAP sources are specifications wherever QPA
delegates to them.

The current campaign targets the complete Rust API first. Python bindings are
deferred until the Rust port passes the source ledger, imported upstream tests,
and pinned GAP/QPA differential tests.

The initial exact coefficient contract uses Numerica and admits the rational
field and supported prime fields. Finite extension fields are a later
coefficient tier.

See [ROADMAP.md](ROADMAP.md) for scope, provenance, and completion criteria.

