# Rust quickstart

The repository currently provides workspace crates rather than a crates.io
release. From this checkout, run the checked examples with:

```console
cargo run -p qpa-core --example quickstart
cargo run -p qpa-gbnp --example quotient
```

For another local Cargo project, use path dependencies until the first
published release:

```toml
[dependencies]
qpa-core = { path = "/path/to/qpa-rust/crates/qpa-core" }
qpa-gbnp = { path = "/path/to/qpa-rust/crates/qpa-gbnp" }
```

## A quiver, path algebra, and right module

The compile-checked source is
[`crates/qpa-core/examples/quickstart.rs`](../../crates/qpa-core/examples/quickstart.rs).
Its essential construction is:

```rust
use qpa_core::api::prelude::*;

let v1 = VertexIndex::new(1).unwrap();
let v2 = VertexIndex::new(2).unwrap();
let quiver = Quiver::from_indexed(
    2,
    vec![ArrowSpec::named(v1, v2, "a")],
)?;
let path_algebra = PathAlgebra::new(Q, quiver)?;

let action = DenseMatrix::from_rows(vec![vec![Q.one()]], Q)?;
let module = right_module_over_path_algebra(
    QuiverAlgebra::Path(path_algebra),
    vec![1, 1],
    &[NamedArrowMatrix::new("a", action)],
)?;
```

The final value is `Option<PathAlgebraModule<_>>`: `None` is QPA's `fail` when
the entered matrices do not satisfy quotient relations, while `Err` reports
malformed input or a failed exact computation. See
[`errors-and-limits.md`](errors-and-limits.md).

## A quotient through GBNP

[`crates/qpa-gbnp/examples/quotient.rs`](../../crates/qpa-gbnp/examples/quotient.rs)
constructs `F_2[x]/(x^2)`. It calls `quotient_by_relations`, which follows
QPA's configured `GBNPGroebnerBasis` call chain; it does not use an independent
Rust Gröbner algorithm.

Use `qpa_core::api` for foundational quiver, algebra, matrix, and module work.
Use `qpa_gbnp::api` when the QPA method delegates to GBNP or needs the higher
module/homological layer. Both facades provide a small `prelude`; the named
facade submodules are preferable in larger programs.
