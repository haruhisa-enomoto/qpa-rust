# GAP/QPA to Python map

Python uses snake case and properties, while retaining the pinned QPA methods
as its computational specification.

| GAP/QPA | Python entry point | Primary QPA source |
| --- | --- | --- |
| `Quiver(...)` | `qpa.Quiver(...)`, `.indexed(...)`, `.from_adjacency(...)` | `lib/quiver.gi` |
| `PathAlgebra(K, Q)` | `qpa.PathAlgebra(field, quiver)` | `lib/pathalg.gi` |
| `KQ / relators` | `path_algebra.quotient(relations)` | `lib/pathalg.gi`, `lib/gbnp.gi`, GBNP `SGrobner` |
| `RightModuleOverPathAlgebra` | `algebra.module_from_integer_matrices(...)` or the `QuiverAlgebra.right_module_over_path_algebra*` family | `lib/module.gi` |
| `HomOverAlgebra` | `algebra.hom(source, range)` | `lib/modulehom.gi` |
| `RightModuleHomOverAlgebra` | `source.homomorphism_to(range, matrices)` | `lib/modulehom.gi` |
| `Kernel`, `Image`, `Cokernel` | `map.kernel()`, `.image_module()`, `.cokernel()` | `lib/modulehom.gi` |
| `IndecProjectiveModules`, `SimpleModules` | `algebra.indecomposable_projectives()`, `.simple_modules()` | `lib/modulehomalg.gi` |
| `ProjectiveCover`, `1stSyzygy` | `module.projective_cover()`, `.first_syzygy()` | `lib/modulehomalg.gi` |
| `ExtOverAlgebra` | `module.ext1(coefficients_in)` | `lib/modulehomalg.gi` |
| `ProjectiveResolution` | `module.projective_resolution()` | `lib/complex.gi` and module routines |
| `DualOfModule`, `TransposeOfModule`, `DTr` | `module.dual()`, `.transpose()`, `.d_tr()` | `lib/functors.gi` |

## Translation example

GAP/QPA:

```gap
Q := Quiver(2, [[1, 2, "a"]]);;
KQ := PathAlgebra(Rationals, Q);;
M := RightModuleOverPathAlgebra(KQ, [1, 1], [["a", [[1]]]]);;
```

Python:

```python
import qpa

quiver = qpa.Quiver.indexed(2, [(1, 2, "a")])
path_algebra = qpa.PathAlgebra(qpa.Field.rational(), quiver)
module = path_algebra.module_from_integer_matrices([1, 1], [("a", [[1]])])
```

The complete generated [QPA operation index](../python/qpa-operation-index.md)
contains all 888 grouped QPA names, their Python targets, statuses, and pinned
source files. Use it for less common operations and overloads.
