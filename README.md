# qpa-rust

**qpa-rust computes with finite-dimensional algebras and their representations.**
Give it a quiver with relations over a field, and it computes the modules over
the resulting algebra: simples, projectives, injectives, radicals and socles,
Hom spaces, direct sum decompositions, projective resolutions, Ext groups,
projective and global dimensions, the Auslander–Reiten translate and almost
split sequences, tilting and τ-tilting checks, and more. All arithmetic is
exact, so every answer is a theorem about the algebra rather than a
floating-point estimate.

It is a Rust port of [QPA](https://gap-packages.github.io/qpa/), the
quivers-and-path-algebras package for the computer algebra system
[GAP](https://www.gap-system.org), with a Python interface. If you have used
QPA before, everything you know carries over; if you have not, you can start
from Python without learning GAP.

## What it does

```python
import qpa

# The algebra A = kQ / (ab) for the quiver 1 --a--> 2 --b--> 3 over GF(5).
field = qpa.Field.prime(5)
quiver = qpa.Quiver(["1", "2", "3"], [("1", "2", "a"), ("2", "3", "b")])
kQ = qpa.PathAlgebra(field, quiver)
a, b = kQ.generator("a"), kQ.generator("b")
A = kQ.quotient([a * b])

S1, S2, S3 = A.simple_modules()
P1, P2, P3 = A.indecomposable_projectives()

assert A.dimension() == "5"
assert A.global_dimension_bounded(10) == "2"
assert S1.projective_dimension(10) == 2
assert S1.ext1(S2).dimension == 1

# Projective resolution of S1: 0 -> P3 -> P2 -> P1 -> S1 -> 0
resolution = S1.projective_resolution()
assert [resolution.object(d).dimension_vector for d in range(3)] == [
    [1, 1, 0], [0, 1, 1], [0, 0, 1],
]

# Auslander-Reiten theory
assert S2.d_tr().dimension_vector == [0, 0, 1]          # tau(S2) = S3
f, g = S2.almost_split_sequence()                       # 0 -> S3 -> P2 -> S2 -> 0
assert f.range.dimension_vector == [0, 1, 1]

# Decomposition and tilting
M = qpa.MatrixModule.direct_sum([P1, P2, S1])
assert sorted(N.dimension_vector for N in M.decompose()) == [
    [0, 1, 1], [1, 0, 0], [1, 1, 0],
]
assert qpa.MatrixModule.direct_sum([P1, P2, P3]).tilting(10) is not None
assert A.is_finite_type_algebra and A.is_gentle
```

Modules are given by a dimension vector and one matrix per arrow, exactly as
in QPA, and can be built by hand or obtained from constructions such as
projective covers, syzygies, kernels, cokernels, tensor products, duals, and
Nakayama functors. Algebras can be entered as a quiver with relations or taken
from built-in families: Nakayama, canonical, Kronecker, poset incidence, Brauer
configuration, and truncated path algebras. Coefficients can be the rationals,
prime fields, or finite extension fields `GF(p^n)`.

## How it relates to QPA

Everything computed here follows QPA's own algorithms. QPA's GAP source code is
the specification: each QPA operation is translated method by method into Rust,
and where QPA delegates to GAP or to the
[GBNP](https://gap-packages.github.io/gbnp/) package for noncommutative Gröbner
bases, those routines are translated too. The port does not substitute
independent algorithms for QPA's, so results agree with QPA on the same input.
Exact arithmetic over the rationals and prime fields comes from
[Numerica](https://github.com/symbolica-dev/numerica); finite extension fields
use GAP's Conway polynomial data so that generators and subfield embeddings
match GAP.

The whole of QPA 1.37 is covered: all 1,839 declarations and installed methods
in its source are mapped to Rust, and 888 QPA operations are reachable from
Python. A generated ledger records, for every operation, its QPA source
location, its Rust and Python names, and its tests. The few places where the
port deliberately behaves differently from QPA, for example to repair a crash
in the original, are listed in
[porting/deviations.md](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/porting/deviations.md).

## Installation

The package is published on PyPI as `qpa-rust` and imported as `qpa`.
Prebuilt wheels exist for Linux, macOS, and Windows on Python 3.10 through
3.13, so no Rust toolchain is needed:

```console
uv add qpa-rust
```

or `pip install qpa-rust`. To use the development source instead, install from
the Git repository or from a local clone; that compiles the Rust extension and
needs a Rust toolchain:

```console
uv add "qpa-rust @ git+https://github.com/haruhisa-enomoto/qpa-rust"
uv add --editable /path/to/qpa-rust
```

## Learning more

The [user guide](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/docs/user/README.md)
covers installation, quivers and path algebras, modules and homomorphisms,
homological computations, Gröbner bases, and the mathematical conventions
shared with QPA. A
[GAP-to-Python name map](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/docs/user/gap-to-python.md)
helps translate existing QPA scripts, and the
[operation index](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/docs/python/qpa-operation-index.md)
links every Python function to the QPA source it implements. The examples in
[docs/examples/python](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/docs/examples/python)
run as part of the test suite.

Rust users can depend on the `qpa-core` and `qpa-gbnp` crates directly; see the
[Rust quickstart](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/docs/user/quickstart.md).

## Status

The current release is 0.2.0.
See the [changelog](CHANGELOG.md) for release changes and compatibility notes.
The Python API may still change. Performance has not
yet been compared with QPA under GAP. Computations that are expensive in QPA,
such as Gröbner bases and submodule enumeration, retain QPA's algorithms and
can still be expensive. The [performance audit](benchmarks/results/broad-audit.json)
records comparisons between Rust builds, consumer workload coverage, and
remaining limitations. Reproduction scripts are in [benchmarks](benchmarks/).

## Credits and licensing

This project would not exist without the original software, and users of the
port should cite the original authors' work.

- **QPA** (Quivers and Path Algebras) is developed by Edward Green, Øyvind
  Solberg, and the QPA team: <https://gap-packages.github.io/qpa/>,
  <https://github.com/gap-packages/qpa>. This port translates QPA version 1.37.
- **GAP** (Groups, Algorithms, Programming) is the computer algebra system QPA
  runs in and supplies the field, matrix, and algebra semantics the port
  follows: <https://www.gap-system.org>.
- **GBNP** by Arjeh M. Cohen, Dié A. H. Gijsbers, Jan Willem Knopper, and Chris
  Krook provides the noncommutative Gröbner basis algorithms that QPA uses:
  <https://gap-packages.github.io/gbnp/>.
- **Numerica** provides exact scalar and linear algebra:
  <https://github.com/symbolica-dev/numerica>.

qpa-rust is not an official release of QPA, GAP, or GBNP, and is not endorsed
by their authors. It is distributed under the GNU General Public License,
version 2 or later, the same license as QPA and GAP. See
[LICENSE](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/LICENSE) and
[NOTICE.md](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/NOTICE.md)
for the exact upstream versions, revisions, and license terms.

## Development

Maintainer checks, including the differential tests against GAP that verify
the port's results against QPA itself, are described in
[docs/development/verification.md](https://github.com/haruhisa-enomoto/qpa-rust/blob/main/docs/development/verification.md).
