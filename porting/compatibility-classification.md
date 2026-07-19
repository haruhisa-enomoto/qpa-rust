# Compatibility classification

This index classifies every entry in `deviations.md` by cause. The detailed
entry remains authoritative for source locations, Rust behavior, compatibility
impact, and tests. Classification is not a claim that every difference is an
upstream bug.

## Classes

- **QPA defect**: an installed QPA method reaches recursion, an unbound field,
  an impossible dispatch, or another unintended interpreter failure on an
  ordinary input in its advertised domain.
- **GAP/QPA interaction defect**: the QPA method expresses a valid computation
  but relies on GAP internals or generic machinery that fails at the pinned
  combination of revisions.
- **Precondition guard**: the source algorithm is correct on its documented
  domain but can fail to terminate or fail late outside it.
- **Rust representation boundary**: Rust ownership, typing, deterministic I/O,
  or absence of GAP process-global facilities requires a different API shape;
  the computation is unchanged.
- **Coefficient-tier limit**: QPA accepts a coefficient domain outside the
  agreed initial Numerica contract.

## Complete index

| Deviation | Classification | Upstream candidate | Rust policy |
|---|---|---:|---|
| COMBMAP-001 | Precondition guard | guard proposal | Detect a repeated ordering orbit and return a typed error. |
| COMBMAP-002 | Precondition guard | guard proposal | Validate curve data and repeated/nondecreasing search states. |
| ALG-001 | QPA defect | yes | Repair the unambiguous self-recursive display branch. |
| ALG-002 | Rust representation boundary | no | Use an owner-checked optional name lookup. |
| ALG-003 | GAP/QPA interaction defect | yes | Execute QPA's requested exact structure-constant radical path. |
| DT-001 | Rust representation boundary | no | Use one owner-tagged safe arena. |
| SUFFIX-001 | Rust representation boundary | no | Use owner-tagged node and edge handles. |
| PATHALG-001 | GAP/QPA interaction defect | yes | Apply the generator map specified by the constructor. |
| PATHALG-002 | Coefficient-tier limit | no | Accept only exact supported fields. |
| PATHALG-003 | QPA defect | yes | Preserve both unreachable overloads as typed failures. |
| PATHALG-004 | QPA defect | yes | Preserve the unbound-relators branch as a typed failure. |
| PATHALG-005 | Rust representation boundary | no | Return an explicit namespace assignment plan. |
| PATHALG-006 | Rust representation boundary | no | Return the complete explicit finite table group. |
| ORDERING-001 | Rust representation boundary | no | Reject cross-owner orderings at construction. |
| INVLENLEX-001 | Rust representation boundary | no | Represent QPA's polymorphic zero sentinel with an enum. |
| GROEBNER-001 | Rust representation boundary | no | Return owned transformed values rather than mutate a standalone value. |
| IDEMPOTENT-001 | Rust representation boundary | no | Make RNG state explicit and reproducible. |
| IDEMPOTENT-002 | Rust representation boundary | no | Cache owned coordinates with stable parent identity. |
| MODULE-001 | QPA defect | yes | Preserve the scalar-zero display failure as a typed error. |
| MODULE-002 | QPA defect | yes | Preserve the unreachable zero-basis branch as a typed error. |
| MODULE-003 | QPA defect | yes | Report the omitted isolated-vertex sentinel repair in the NC method. |
| MODULE-004 | QPA defects/dispatch fragility | yes | Preserve each observed failing branch as a distinct typed error. |
| MODULE-005 | QPA defects/dispatch fragility | yes | Preserve zero/empty/quotient Hom and End failures without substituting another Hom algorithm. |
| MODULE-006 | Rust representation boundary | no | Require the finite based action data read by QPA's method. |
| MODULEHOM-001 | Coefficient-tier limit | no | Stop where QPA requires unsupported central polynomial splitting. |
| MODULEHOM-002 | QPA defect | yes | Treat the empty radical basis as the zero space in the same construction. |
| MODULEHOMALG-001 | QPA defect | yes | Preserve the failed intermediate as a typed error. |
| MODULEPROJRES-001 | Precondition guard | guard proposal | Detect a complete reduction pass that changes nothing. |
| DEGORDER-001 | Rust representation boundary | no | Use typed zero-based positions with explicit conversion. |
| DEGORDER-002 | Rust representation boundary | no | Return stable data instead of mutating process-global output. |
| DEGORDER-003 | Rust representation boundary | no | Use a typed selector and expose catalog text as data. |
| DOCUMENTATION-001 | Rust representation boundary | no | Return the GAPDoc build plan without filesystem side effects. |

## Reporting policy

The candidates above are prepared locally in `upstream-issues/`. Nothing is
submitted to an upstream tracker without human review. A report must include
the pinned revisions, a minimal script, actual behavior, expected behavior, and
the exact QPA source location. GAP/QPA interaction candidates should be checked
against a supported stable GAP release before being described as QPA-only bugs.
