# Homological calculations in Python

The homological Python methods call the translated QPA module and complex
routines. Lazy resolutions preserve QPA's degree conventions and caches.

## Syzygies and projective dimension

```python
import qpa

field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")])
algebra = qpa.PathAlgebra(field, quiver)
simple_at_v1, simple_at_v2 = algebra.simple_modules()

syzygy = simple_at_v1.first_syzygy()
assert syzygy.dimension_vector == [0, 1]
assert simple_at_v1.projective_dimension(3) == 1
```

The bound in `projective_dimension(bound)` is a computational bound, not a
claim that every module has projective dimension at most that value. A `None`
result means that the translated method did not establish a finite answer
within the requested bound.

## Ext

For `source` and `target` over the same acting algebra, choose a result according
to the data you need:

| Calculation | Return value |
| --- | --- |
| `source.hom_dimension_to(target)` | Hom dimension as an integer |
| `source.hom_basis_to(target)` | List of `ModuleHomomorphism` basis maps |
| `source.hom_basis_with_coefficients_to(target)` | `ModuleHomomorphismBasis` with `.dimension`, `.basis`, and `.coefficients(map)` |
| `source.ext1_dimension(target)` | Ext¹ dimension as an integer |
| `source.ext1(target)` | `ExtOverAlgebraResult` with `.dimension`, `.basis`, and `.coefficients(map)` |
| `source.ext1_dimensions(targets)` | List of Ext¹ dimensions in target order |

The Hom methods can return `None` when the acting algebras differ. Ext methods
reject incompatible modules. `ext1_dimensions` batches coefficient modules;
it does not vary the degree.

`module.ext1(coefficients_in)` returns a structured result containing the QPA
basis and its dimension:

```python
ext = simple_at_v1.ext1(simple_at_v2)
assert ext.dimension == 1
assert len(ext.basis) == 1
```

Here the receiver is the first argument of Ext¹ and `coefficients_in` is the
second. Basis representatives are module maps from a syzygy of the receiver
to the coefficient module, modulo maps extending to its projective cover.
There is no general `ext(target, degree)` method. For higher-degree dimensions,
dimension shifting gives
`source.nth_syzygy(n - 1).ext1_dimension(target)` for `n >= 1`; degree zero is
`source.hom_dimension_to(target)`. Request only the degrees you need.

## Lazy projective resolutions

```python
resolution = simple_at_v1.projective_resolution()

assert resolution.object(0).dimension_vector == [1, 1]
assert resolution.object(1).dimension_vector == [0, 1]
assert resolution.object(2).is_zero
```

QPA includes the resolved module in degree `-1`; `object(0)` is its first
projective. Calling `object(n)` or `differential(n)` extends the lazy
resolution as needed and reuses cached terms.

`projective_resolution()` takes no bound argument. To inspect the first three
projective terms, use `[resolution.object(n) for n in range(3)]`.
The source of `differential(n)` is `object(n)` and its target is
`object(n - 1)`. In particular, `differential(0)` is the projective cover onto
the original module. A differential is a `ModuleHomomorphism`; access its
`vertex_matrices` property or call `underlying_linear_matrix()` for its full
linear matrix. Matrices act on row vectors, and `f.compose(g)` applies `f`
before `g`.

For a minimal resolution over these split quiver algebras,
`resolution.object(n).top().dimension_vector` records the multiplicities of
indecomposable projectives. For the dual numbers there is one indecomposable
projective of dimension two: its multiplicity is one, so projective term
dimensions `[2, 2, 2]` correspond to Betti numbers `[1, 1, 1]`.

## Finite complexes and exact sequences

Use `module.finite_complex(base_degree, differentials)` to build a finite
complex. The module selects the acting algebra. The first list entry is
`d_base_degree: C_base_degree -> C_(base_degree - 1)`; following entries
have increasing source degrees. For example, a zero map `S -> S` in degrees
1 and 0 is `S.finite_complex(1, [S.zero_map_to(S)])`.

The kernel module and its inclusion are separate objects:

```python
cover = simple_at_v1.projective_cover()
kernel = cover.kernel()                  # MatrixModule
inclusion = cover.kernel_inclusion()     # ModuleHomomorphism: kernel -> source
assert inclusion.source.isomorphic_to(kernel)
assert inclusion.compose(cover).is_zero

exact = simple_at_v1.finite_complex(0, [cover, inclusion])
assert exact.is_short_exact_sequence
assert all(exact.homology(n).is_zero for n in (-1, 0, 1))
```

This complex places the resolved module in degree −1, the projective in degree
0, and its kernel in degree 1. Alternatively,
`inclusion.short_exact_sequence_with(cover)` constructs and checks the short
exact sequence directly. `complex.homology(degree)` returns a `MatrixModule`.
Finite bounds such as `complex.lower_bound` and `complex.upper_bound` are
strings, as the same properties also represent infinite bounds.

For a complex version of a lazy resolution, use
`module.projective_resolution_complex()`. To construct a finite chain map,
use `source.finite_chain_map_to(target, base_degree, component_maps)`; here
`base_degree` is the degree of the first component map. The binding currently
does not provide a complex identity convenience method or an induced map on
homology.

The checked [resolution workflow](../examples/python/resolution_workflow.py)
computes the first three terms over the dual numbers, verifies that consecutive
differentials compose to zero, compares Hom/Ext dimensions with their bases,
and constructs the short exact sequence of a projective cover.

## Minimal Betti data from a GSZ resolution

A GSZ resolution need not be minimal. For an admissible quotient,
`GszProjectiveResolution.semisimple_homology()` reduces its differentials
modulo the arrow ideal and computes exact vertexwise homology:

```python
import qpa

field = qpa.Field.prime(5)
quiver = qpa.Quiver(["v"], [("v", "v", "x")])
algebra = qpa.QuotientAlgebra.truncated(field, quiver, 2)
x = algebra.generator("x")
gsz = algebra.projective_resolution_fp_path_algebra_module(
    algebra.ideal, [[x]], 4
)
betti = gsz.semisimple_homology()

assert betti.minimal_projective_betti_multiplicities == [[1], [1], [1], [1]]
```

The homology rows are the indecomposable-projective multiplicities in the
minimal resolution, even when the stored GSZ terms contain cancellable
summands. The result also exposes the raw projective multiplicities and the
vertexwise ranks of the reduced differentials. A nonterminal GSZ prefix
through `d_T` certifies Betti rows only through degree `T - 1`; an empty
terminal differential also certifies the preceding final projective degree.
The operation rejects nonadmissible quotients.

The checked complete example is
[`docs/examples/python/homological.py`](../examples/python/homological.py).
