# Python API and provenance references

There are three complementary sources of API information.

## Task-oriented guide

The pages in this directory explain mathematical workflows and the conventions
that signatures alone cannot express. Begin with the
[Python quickstart](python-quickstart.md), not the alphabetical reference.

## Complete Python signature reference

The generated [Python API reference](../python/api-reference.md) contains every
public class, function, method, property, and type annotation shipped in
`qpa._qpa.pyi`. Editors and type checkers consume the same stub through the
package's `py.typed` marker.

From Python, `help(qpa.Field)`, `help(qpa.Quiver)`, and the other core classes
show runtime PyO3 documentation. Explanations for the primary workflow methods
are also kept in the stub and task guide so they remain visible to static tools.

When a method is unexpectedly missing, first check the installed build and the
actual type of the value. A tensor construction, an algebra, a module, and a
map are separate objects even when one computation returns the next:

```python
import importlib.metadata
import qpa

print(importlib.metadata.version("qpa-rust"))
print(qpa.PORT_GIT_COMMIT)
print(qpa.__file__)
```

The published Linux 0.2.0 wheels report `PORT_GIT_COMMIT` as `"unknown"`
because of a build-container metadata omission. They were built from
`e9695a0b23c8dd4538db29ab33f6bc9e870eabeb` (`v0.2.0`); package version and file
path checks still work. The release workflow checks this field for future
Linux wheels.

Use `type(value)` and the reference for that class rather than inferring a
return type from a variable name. For example, `algebra.enveloping_algebra()`
returns a `TensorProductAlgebra`; `.algebra` gives its common algebra handle,
and `.quotient` gives its concrete quotient. Both algebra handles provide
standard projective, injective, and simple modules. To unwrap a
`QuiverAlgebra`, use `.quotient_algebra` or `.path_algebra` according to
`.is_quotient`. Its `.original_path_algebra` is the ambient algebra and would
discard the quotient relations.

The [module guide](python-modules.md) distinguishes module/map results and
the [homological guide](python-homological.md) records resolution and complex
degree conventions. Examples are executed in the test suite. A separate API
contract test checks the shipped class members, method/property distinction,
and callable parameters against the installed extension.
The [API consistency audit](../development/python-api-audit.json) records the
tested workflows, fixes, retained conventions, and remaining feature gaps.

## GAP/QPA operation and source index

The generated [QPA operation index](../python/qpa-operation-index.md) maps all
888 grouped GAP/QPA names to Python targets and pinned source files. It is the
right reference when the question is “where did this operation come from?” or
“what is the Python spelling of this QPA method?”

For longer mathematical descriptions of uncommon methods, read the
[pinned QPA manual chapters](qpa-source-manual.md) alongside the operation
index. The manual uses GAP syntax; the index supplies the corresponding Python
target without changing the documented QPA algorithm.

Both generated pages are checked into the repository for browsing without a
documentation tool. Regenerate and verify them with:

```console
python3 scripts/generate-python-docs
scripts/check-python-docs
```

The binding completeness ledger remains the authority for coverage status;
see the [maintainer ledger documentation](../development/python-binding-ledger.md).
