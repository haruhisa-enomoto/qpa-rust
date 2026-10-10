"""Executable resolution documentation and finite-complex degree conventions."""

import runpy
from pathlib import Path

import pytest
import qpa


def test_documented_resolution_workflow_runs():
    example = (
        Path(__file__).resolve().parents[2]
        / "docs/examples/python/resolution_workflow.py"
    )
    runpy.run_path(str(example), run_name="__main__")


@pytest.mark.parametrize("field", [qpa.Field.prime(2), qpa.Field.prime(5)])
def test_finite_resolution_segment_has_expected_endpoint_homology(field):
    algebra = qpa.QuotientAlgebra.truncated(
        field, qpa.Quiver.indexed(1, [(1, 1, "x")]), 2
    )
    simple = algebra.simple_modules()[0]
    resolution = simple.projective_resolution()
    d1, d2 = (resolution.differential(degree) for degree in (1, 2))
    segment = simple.finite_complex(1, [d1, d2])
    assert [segment.object(degree).dimension for degree in range(3)] == [2, 2, 2]
    assert segment.differential(1) == d1
    assert segment.differential(2) == d2
    assert d2.compose(d1).is_zero
    # Truncation leaves kernel/cokernel at the ends, with exactness in degree 1.
    assert [segment.homology(degree).dimension for degree in range(3)] == [1, 0, 1]
    assert segment.object(-1).is_zero
    assert segment.object(3).is_zero
    assert simple.nth_syzygy(0).is_identical(simple)
    with pytest.raises(qpa.InvalidInputError, match="negative"):
        simple.nth_syzygy(-1)
