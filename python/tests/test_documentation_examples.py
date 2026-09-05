from __future__ import annotations

import runpy
import re
from pathlib import Path

import pytest

import qpa


EXAMPLES = Path(__file__).resolve().parents[2] / "docs" / "examples" / "python"


@pytest.mark.parametrize(
    "example",
    ["quickstart.py", "algebras.py", "modules.py", "homological.py", "groebner.py"],
)
def test_documented_python_example_runs(example: str) -> None:
    runpy.run_path(str(EXAMPLES / example), run_name="__main__")


def test_local_documentation_links_resolve() -> None:
    root = Path(__file__).resolve().parents[2]
    markdown_files = [root / "README.md", *(root / "docs").rglob("*.md")]
    # The pinned upstream clones are a maintainer-only development reference
    # and are absent from distributed checkouts; links into them are verified
    # only when the directory is present.
    upstream = root / "upstream"
    missing: list[tuple[Path, str]] = []
    for markdown in markdown_files:
        text = markdown.read_text(encoding="utf-8")
        for raw_target in re.findall(r"\[[^]]*\]\(([^)]+)\)", text):
            target = raw_target.split("#", 1)[0]
            if not target or "://" in target or target.startswith(("mailto:", "<")):
                continue
            resolved = (markdown.parent / target).resolve()
            if not upstream.exists() and resolved.is_relative_to(upstream):
                continue
            if not resolved.exists():
                missing.append((markdown.relative_to(root), raw_target))
    assert missing == []


def test_primary_runtime_api_has_help_text() -> None:
    documented = {
        qpa.Field.prime: "prime field",
        qpa.Quiver.path: "left-to-right",
        qpa.PathAlgebra.quotient: "default GBNP route",
        qpa.PathAlgebra.groebner_basis: "Gröbner basis",
        qpa.MatrixModule.first_syzygy: "projective cover",
        qpa.MatrixModule.projective_resolution: "degree `-1`",
        qpa.ModuleHomomorphism.compose: "row-action order",
    }
    for operation, phrase in documented.items():
        assert operation.__doc__ is not None
        assert phrase in operation.__doc__
