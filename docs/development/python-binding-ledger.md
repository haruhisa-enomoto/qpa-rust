# Python binding ledger

[`porting/python-bindings.toml`](../../porting/python-bindings.toml) is the
complete operation-oriented inventory for Python exposure. It groups all 1,839
QPA declarations and installed overloads into 888 distinct QPA names. Every
source-ledger entry ID occurs in exactly one group.

The compact [`python/bindings.toml`](../../python/bindings.toml) manifest records
the exposure decision for every grouped QPA name. Regeneration combines it with
[`porting/qpa-operations.toml`](../../porting/qpa-operations.toml); an omitted
entry would remain explicitly `planned`, never silently unaccounted for.

The completed initial-field binding has 883 `exposed`, two `represented`, and
three `unavailable` groups. It has no `planned` or `partial` groups. The three
unavailable groups are pinned QPA declarations with no installed callable
method, not replacement algorithms or silently omitted behavior.

## Dispositions

- `operation`: a QPA operation, attribute, property, global function, installed
  overload, or helper callable that needs one Python exposure decision.
- `runtime_type_or_protocol`: a GAP category, representation, family, or
  filter represented by Python wrapper classes and runtime validation.
- `alias`: a QPA synonym that should resolve to the canonical Python target.
- `diagnostic_channel`: a GAP `InfoClass` that maps to documented Python
  diagnostics rather than a mathematical operation.

## Statuses

- `planned`: no Python exposure has yet been recorded.
- `partial`: at least one translated QPA overload is callable, but the grouped
  QPA name is not completely covered.
- `exposed`: all applicable translated behavior for the group is callable.
- `represented`: the behavior is supplied by a Python type/protocol or another
  named binding rather than a standalone callable.
- `unavailable`: QPA declares no applicable method or the Rust port records an
  explicit coefficient-tier boundary; the reason must be documented.

Run `scripts/check-python-ledger` to regenerate the report and prove that all
1,839 source entries occur exactly once. Exposure is complete only when no
operation remains `planned` or `partial` and each final status has Python tests
or a documented unavailability reason.
