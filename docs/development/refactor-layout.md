# Refactor layout

The post-port refactor changes file layout only. It deliberately does not
change algorithms, item names, visibility, or public import paths from the
frozen `rust-port-reference-0.1.0` baseline.

## Why the domain files use `include!`

The original QPA translation accumulated mutually dependent private caches,
weak identity handles, public types, and functions in two Rust modules. Moving
those items into ordinary child modules would change privacy and name
resolution and would require a broad semantic rewrite. The facade files now
use private `include!` chunks instead: Rust still compiles each chunk in the
original module, so `qpa_core::path_module::*` and `qpa_gbnp::module::*` remain
unchanged while reviewers can work on one coherent source region at a time.

The chunks retain their original source order. Public item documentation names
the translated QPA operation and records non-obvious dispatch, ordering,
failure, and caching behavior next to the implementation.

## `qpa-core::path_module`

| File | Responsibility |
| --- | --- |
| `objects.rs` | Matrix-module objects, elements, bases, submodules, and generic based-module objects |
| `endomorphisms.rs` | Module-map objects and endomorphism/stable-endomorphism algebra representations |
| `homomorphisms.rs` | Checked maps, kernels, images, Hom and End computations |
| `submodules.rs` | Generic Hom conversion, generated submodules, radicals, tops, socles, and traces |
| `categorical.rs` | Direct sums, pushouts, pullbacks, induced maps, intersections, sums, and isomorphism tests |
| `constructors.rs` | Dimension-inferring and explicit-dimension QPA module constructors |
| `errors.rs` | Typed compatibility and construction failures |
| `tests.rs` | The formerly inline regression and oracle-fixture tests |

## `qpa-gbnp::module`

| File | Responsibility |
| --- | --- |
| `algebra_invariants.rs` | Enveloping modules, projective/injective families, homological dimensions, Cartan/Coxeter and Nakayama data |
| `projectives.rs` | Algebra properties, annihilators, projective maps, and projective covers |
| `syzygies_and_tilting.rs` | Syzygies, approximations, tilting/cotilting, complements, and mutations |
| `approximations.rs` | Perpendicular/list/radical approximations and simple endomorphism-algebra resolutions |
| `ext_and_iyama.rs` | Injective envelopes, Ext, Iyama generators, rigidity, and complexity |
| `auslander_reiten.rs` | Almost-split data, irreducible maps, predecessors, and finite-field module enumeration |
| `stable_and_functors.rs` | Stable maps, dual/transpose/star/AR/Nakayama functors and their map actions |
| `tensor_and_enveloping.rs` | Tensor products of modules and the dual regular bimodule |
| `errors.rs` | Typed failures shared by the included source regions |
| `tests.rs` | The formerly inline regression and oracle-fixture tests |

The include order in each facade is dependency-neutral because all items are
compiled in one Rust module, just as before the split. Every structural change
must still pass both coefficient backends, strict lint/rustdoc checks, the
operation ledger, and the pinned GAP/QPA oracle suite described in
[`verification.md`](verification.md).
