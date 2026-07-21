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
