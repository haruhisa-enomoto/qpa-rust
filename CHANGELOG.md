# Changelog

## 0.2.0 — 2026-10-10

This release improves exact module calculations and makes the Python algebra
interfaces more consistent. Existing supported calls retain their names and
return types; the new dimension queries and construction conveniences are
optional.

### Performance

- Reuse matrix elimination results and avoid computing transformations or
  basis vectors when only a solution, rank, or dimension is needed.
- Use packed binary elimination for matrices over `GF(2)`.
- Reuse projective families during projective-cover construction, benefiting
  syzygy and Ext calculations.
- Share immutable Gröbner-basis data and reuse tip dictionaries, nontips, and
  canonical coordinate ordering. Reuse the right Gröbner basis within each GSZ
  resolution step.
- Skip zero and unit arithmetic in matrix products, avoid cloning arrow
  matrices when comparing modules, and avoid redundant validation when
  composing already checked module maps. Public map constructors continue to
  check all arrow squares.

These improvements apply to existing calculations without an API migration.
The [performance audit](benchmarks/results/broad-audit.json) records measured
results, reproduction commands, and limitations. Comparisons are between Rust
builds, not against GAP/QPA; gains depend on the workload, and some small probes
became slower.

### Added

- `MatrixModule.hom_dimension_to(target)`, `ext1_dimension(target)`, and
  `ext1_dimensions(targets)` compute dimensions without constructing the full
  representative basis. The last method batches targets for Ext¹, not degrees.
  Corresponding Rust Hom/Ext dimension operations are also available.
- `QuiverAlgebra` now provides standard simple, projective, injective, and zero
  modules; integer-matrix module construction; named generators; the
  `generators` property; and `zero()` and `one()`.
- `PathAlgebra` and `QuotientAlgebra` gain `as_quiver_algebra()`,
  `tensor_product()`, and `enveloping_algebra()`. Tensor products accept any
  combination of concrete and common algebra handles.
- `module_from_matrices()` on all three algebra classes accepts exact
  `Scalar`, Python `int`, and `fractions.Fraction` entries.
  `MatrixModule.homomorphism_to()` now accepts these same entry types.
- Field scalars support arithmetic with scalars from the same field and Python
  integers, including reflected operations and signed integer powers.
- Concrete and common algebra elements support nonnegative integer powers
  through `**`; common elements also support unary negation.
- `ModuleElement.act_by()` accepts concrete or common algebra elements and
  checks their parent algebra and coefficient field.
- Concrete algebra classes expose Cartan/Coxeter calculations and the existing
  selfinjective, symmetric, weakly symmetric, Schurian, semicommutative, and
  distributive predicates.

See the [algebra](docs/user/python-algebras.md),
[module](docs/user/python-modules.md), and
[coefficient](docs/user/coefficient-domains.md) guides for usage.

### Fixed and documented

- Implemented nine `PathAlgebra` invariant methods previously advertised in
  the Python stubs but absent at runtime.
- Corrected matrix-utility keyword names in the stubs, added the missing
  `QuotientAlgebra.right_projective_module()` signature, and removed a duplicate
  path-algebra declaration. Runtime keyword names are unchanged.
- Corrected generated documentation that labeled property setters as methods.
- Added return-type guidance for Hom/Ext results, kernels, images, cokernels,
  direct sums, and tensor constructions, plus an executable
  [resolution workflow](docs/examples/python/resolution_workflow.py).
- Added tests comparing shipped Python signatures with the installed extension
  and executing the workflow examples. The
  [API audit](docs/development/python-api-audit.json) records the findings and
  validation: 588 Python tests and 434 Rust tests on each coefficient backend,
  169 consumer tests, and 29 broader calculations with matching mathematical
  results.
- Release builds pass the source commit for `qpa.PORT_GIT_COMMIT` build
  identification. See the Linux packaging limitation below.

### Compatibility notes

- No existing runtime methods were removed or renamed. Existing result types,
  row-action composition order, and checked module/map construction are
  preserved. No consumer migration is required.
- Exact matrix inputs still reject floats and implicit changes of coefficient
  field. Scalar arithmetic accepts `Scalar` and `int`; construct fractional
  operands with `field.element(numerator, denominator)`.
- New scalar powers require exponent magnitude at most `4294967295`; larger
  magnitudes raise `OverflowError` before reaching the coefficient backend.
  Division by zero and negative powers of zero raise `ZeroDivisionError`.
- Tensor and enveloping constructions still return `TensorProductAlgebra`;
  use `.algebra` or `.quotient` for module calculations. Module dimensions and
  predicates remain properties, and `QuotientAlgebra.dimension()` retains its
  string-or-`None` result.
- Error messages for invalid inputs are not guaranteed to remain identical.
- The published Linux 0.2.0 wheels report `qpa.PORT_GIT_COMMIT == "unknown"`
  because the build container did not receive the commit variable. Their source
  revision is `e9695a0b23c8dd4538db29ab33f6bc9e870eabeb` (`v0.2.0`). This affects
  build identification, not calculations; the installed Linux wheel passes all
  588 Python tests. The release workflow now forwards the variable explicitly
  and checks it before uploading future Linux wheels.

## 0.1.0

Initial public release of the Rust port of QPA with Python bindings, exact
rational and finite-field coefficients, user guides, and generated API and
QPA-operation references.
