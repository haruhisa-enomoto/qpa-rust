# Prepared QPA compatibility reports

These are local, reviewable issue drafts. They have not been submitted. Run a
finite reproducer from the repository root with, for example:

```console
scripts/gap-reference --alwaystrace \
  porting/upstream-issues/reproducers/opposite-display-recursion.g
```

All results below refer to QPA 1.37 at
`8d5f847258d4853700ce3a5d947db2a1441109f1`, GBNP 1.1.0 at
`e0b46c1aac16f3fbabe743e3169dc6bd86bdc654`, and GAP 4.17dev at
`4efe877859ab0bc0ac041bf635a996fb897ee0c9`. GAP/QPA interaction reports need
rechecking on a stable GAP release before submission.

## QPA-001 — Opposite structure-constants algebra recursively prints itself

- Classification: QPA defect; high confidence.
- Source: `lib/algebra.gi:109-115`, `PrintObj`.
- Reproducer: `reproducers/opposite-display-recursion.g`.
- Actual: printing the opposite algebra repeatedly invokes the same `PrintObj`
  method until GAP's recursion limit.
- Expected: print the opposite wrapper and its underlying algebra once.
- Rust: repairs the self-reference with finite deterministic display
  (`ALG-001`).

## QPA-002 — Path list-uniformity overloads dispatch on the path as an element

- Classification: QPA defect; high confidence.
- Source: `lib/pathalg.gi:574-663`, list-plus-path `IsLeftUniform` and
  `IsRightUniform`.
- Reproducer: `reproducers/list-uniformity-dispatch.g`.
- Actual: every reachable call asks `CoefficientsAndMagmaElements` for an
  `IsPath` object and raises a no-method error.
- Expected: either inspect the entered list relative to the path endpoint or
  remove/correct the installed overloads. The intended endpoint rule should be
  confirmed by QPA maintainers.
- Rust: does not guess the rule; it preserves a typed failure (`PATHALG-003`).

## QPA-003 — Quotient `OrderedBy` reads an unbound `relators` component

- Classification: QPA defect; high confidence.
- Source: `lib/pathalg.gi:791-802` and quotient construction at lines 215-309.
- Reproducer: `reproducers/quotient-ordered-by-unbound-relators.g`.
- Actual: an ordinarily constructed quotient has no `fam!.relators`, so
  `OrderedBy` raises a record-element error.
- Expected: reorder using the quotient's retained ideal/relations, or retain
  the expected component during construction.
- Rust: preserves a typed failure because choosing the repair would change
  pinned behavior (`PATHALG-004`).

## QPA-004 — Zero-module display passes scalar zero to `PrintArray`

- Classification: QPA defect; high confidence.
- Source: `lib/module.gi:967-984`, `Display`.
- Reproducer: `reproducers/zero-module-display.g`.
- Actual: an arrow entry of the zero module is stored as integer zero and
  `PrintArray` rejects it.
- Expected: render an appropriately shaped zero matrix or the scalar sentinel.
- Rust: preserves the branch as a typed display error (`MODULE-001`).

## QPA-005 — Zero-module canonical basis cannot reach its empty branch

- Classification: QPA defect; high confidence.
- Source: `lib/module.gi:244-369`, `BasisNC`, `Basis`, and `CanonicalBasis`.
- Reproducer: `reproducers/zero-module-canonical-basis.g`.
- Actual: `BasisNC` fails in `NewBasis` before `CanonicalBasis` reaches its
  explicit empty-generator branch.
- Expected: return the canonical empty basis promised by that branch.
- Rust: preserves the pinned failure as a typed error (`MODULE-002`).

## QPA-006 — NC inferred constructor omits isolated-vertex correction

- Classification: QPA defect; high confidence.
- Source: `lib/module.gi:596-754`, two-argument
  `RightModuleOverPathAlgebraNC`.
- Reproducer: `reproducers/nc-isolated-vertex.g`.
- Actual: the isolated vertex retains dimension sentinel `-1` and later module
  element construction fails.
- Expected: perform the same isolated-vertex `-1` to `0` loop as the checked
  constructor.
- Rust: reports a typed error and directs callers to the checked method
  (`MODULE-003`).

## QPA-007 — Zero submodule conversion indexes an absent basis vector

- Classification: QPA defect/empty-dispatch fragility; high confidence for the
  reproduced branch.
- Source: `lib/module.gi:1070-1128`, `SubmoduleAsModule`.
- Reproducer: `reproducers/zero-submodule-conversion.g`.
- Actual: conversion of the zero submodule reaches a missing basis vector.
- Expected: construct the zero module with the inherited vertex dimensions.
- Rust: preserves this and the related source branches as typed errors
  (`MODULE-004`). The broader consolidated evidence remains in
  `tests/gap-reference/module_submodule_errors.g`.

## QPA-008 — Generic Hom/End loses the width of an empty equation system

- Classification: QPA defect/GAP empty-matrix fragility; high confidence for
  the reproduced arrowless case.
- Source: `lib/module.gi:1255-1441`, `HOPF_MatrixModuleHom`, `Hom`, and `End`.
- Reproducer: `reproducers/generic-end-arrowless.g`.
- Actual: the empty equation matrix does not retain the intended nullspace
  width, so an ordinary arrowless module cannot use the generic End method.
- Expected: return the full matrix endomorphism space constrained by no arrow
  equations.
- Rust: preserves this family of failures rather than substituting
  `HomOverAlgebra` (`MODULE-005`). Consolidated cases are pinned by
  `tests/gap-reference/module_generic_hom_errors.g`.

## QPA-009 — Endomorphism-quiver conversion fails on an empty radical

- Classification: QPA defect; high confidence.
- Source: `lib/modulehom.gi:3634-3650`, delegating to
  `lib/algebra.gi:634-926`.
- Reproducer: `reproducers/end-quiver-empty-radical.g`.
- Actual: GAP's list inner product is asked to sum an empty radical basis.
- Expected: regard it as the zero space and return the one-vertex zero-arrow
  field algebra.
- Rust: implements that zero-dimensional branch of the same construction
  (`MODULEHOM-002`).

## QPA-010 — Finite-field zero Iyama generator indexes `fail`

- Classification: QPA defect; high confidence.
- Source: `lib/modulehomalg.gi:692-715`, `IyamaGenerator`.
- Reproducer: `reproducers/iyama-zero-module.g`.
- Actual: decomposition returns `fail`, after which the method evaluates
  `L[1]`.
- Expected: handle the zero module consistently with the rational branch or
  return a deliberate QPA-level failure.
- Rust: exposes the failed intermediate as a typed error
  (`MODULEHOMALG-001`).

## COMPAT-001 — Mapping-algebra radical fails over GF(2)

- Classification: GAP/QPA interaction defect; stable-GAP confirmation needed.
- Source: `lib/algebra.gi:288-320`, `RadicalOfAlgebra` mapping adapter.
- Reproducer: `reproducers/mapping-radical-gf2.g`.
- Actual: pinned GAP basis machinery passes `fail` to `IsZero` for the
  dual-number mapping algebra over GF(2); the rational analogue succeeds.
- Expected: transpose to a matrix algebra, compute its radical, and transpose
  the basis back.
- Rust: runs exactly that requested structure-constant route (`ALG-003`).

## COMPAT-002 — Natural projection application uses stale mapping internals

- Classification: GAP/QPA interaction defect; stable-GAP confirmation needed.
- Source: `lib/pathalg.gi:769-787`, `ImagesRepresentative`.
- Reproducer: `reproducers/natural-projection-stale-record.g`.
- Actual: the constructor succeeds, but application reads the absent private
  component `alghom!.generators`.
- Expected: use supported mapping accessors or retain the generator data in the
  mapping representation.
- Rust: evaluates the generator mapping installed by QPA's constructor
  (`PATHALG-001`).

## GUARD-001 — Winding search has no different-orbit termination branch

- Classification: documented-precondition guard proposal.
- Source: `lib/combinatorialmap.gi:349-395`, `WindingNumber`.
- Reproducer: `reproducers/winding-different-orbits.g` (intentionally hangs;
  do not include it in unattended test runners).
- Proposal: detect return to the starting orbit and issue a QPA-level input
  error. Valid search order and winding values remain unchanged.
- Rust: detects the repeated state (`COMBMAP-001`).

## GUARD-002 — Topological curve searches assume canonical positive-genus data

- Classification: documented-precondition guard proposal.
- Source: `lib/combinatorialmap.gi:585-930`.
- Reproducer/evidence: the valid algorithm is pinned by
  `tests/gap-reference/combinatorialmap_homology.g`; invalid curve helpers and
  boundary cases are pinned by the `combinatorialmap_*_error.g` scripts.
- Proposal: validate the stated input domain and stop on a repeated state or
  nondecreasing genus.
- Rust: adds only those finite guards (`COMBMAP-002`).

## GUARD-003 — `FirstPart` loops when reducers do not span the vector

- Classification: implicit-precondition guard proposal.
- Source: `lib/moduleprojres.gi:648-682`, `FirstPart`.
- Reproducer: `reproducers/first-part-not-spanned.g` (intentionally hangs; do
  not include it in unattended test runners).
- Proposal: after a complete unchanged reducer pass, report that the vector is
  not spanned. Reducer order and valid results remain unchanged.
- Rust: implements only that stalled-pass check (`MODULEPROJRES-001`).
