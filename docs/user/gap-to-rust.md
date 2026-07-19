# GAP/QPA to Rust map

The Rust names are snake case, but the operations below follow the indicated
pinned QPA methods.

| GAP/QPA | Supported Rust entry point | Primary QPA source |
| --- | --- | --- |
| `Quiver(...)` | `Quiver::from_indexed`, `from_named`, or `from_adjacency` | `lib/quiver.gi` |
| `PathAlgebra(K, Q)` | `PathAlgebra::new(K, Q)` | `lib/pathalg.gi` |
| `KQ / relators` with the default backend | `qpa_gbnp::api::quotients::quotient_by_relations` | `lib/pathalg.gi`, `lib/gbnp.gi`, then GBNP `SGrobner` |
| `RightModuleOverPathAlgebra` | `qpa_core::api::modules::right_module_over_path_algebra` | `lib/module.gi` |
| `RightModuleHomOverAlgebra` | `right_module_hom_over_algebra` | `lib/modulehom.gi` |
| `HomOverAlgebra` / `EndOverAlgebra` | `hom_over_algebra` / `end_over_algebra` | `lib/modulehom.gi` |
| `KernelInclusion`, `ImageProjection`, `CoKernelProjection` | same nouns in snake case | `lib/modulehom.gi` |
| `IndecProjectiveModules`, `SimpleModules`, `ProjectiveCover` | `qpa_gbnp::api::modules::*` | `lib/modulehomalg.gi` |
| `DualOfModule`, `TransposeOfModule`, `DTr`, `StarOfModule` | `qpa_gbnp::api::functors::*` | `lib/functors.gi` |

## Translation example

QPA-style pseudocode:

```gap
Q := Quiver(2, [[1, 2, "a"]]);;
KQ := PathAlgebra(Rationals, Q);;
M := RightModuleOverPathAlgebra(KQ, [1, 1], [["a", [[1]]]]);;
```

Rust:

```rust
use qpa_core::api::prelude::*;

let v1 = VertexIndex::new(1).unwrap();
let v2 = VertexIndex::new(2).unwrap();
let quiver = Quiver::from_indexed(
    2,
    vec![ArrowSpec::named(v1, v2, "a")],
)?;
let kq = PathAlgebra::new(Q, quiver)?;
let matrix = DenseMatrix::from_rows(vec![vec![Q.one()]], Q)?;
let module = right_module_over_path_algebra(
    QuiverAlgebra::Path(kq),
    vec![1, 1],
    &[NamedArrowMatrix::new("a", matrix)],
)?;
```

The outer `Result` replaces GAP errors with typed failures. The inner `Option`
preserves the constructor's possible QPA `fail` result. The full compile-checked
program is linked from [`quickstart.md`](quickstart.md).

For less common operations, search the crate-root rustdoc by the snake-case QPA
name and then check whether it has been admitted to `qpa_core::api` or
`qpa_gbnp::api`. Root re-exports preserve the complete audited translation even
when a convenience facade grouping has not yet been assigned.
