from __future__ import annotations

import json
from pathlib import Path

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def test_public_next_path_queue_protocol_follows_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "quiver_surface.schema-v1.json").read_text(encoding="utf-8")
    )["next_path"]
    quiver = qpa.Quiver(
        ["v1", "v2"], [("v1", "v2", "a"), ("v1", "v1", "b")]
    )

    first = quiver.next_path()
    assert first is not None
    assert repr(first[0]) == fixture["first"]
    assert [repr(path) for path in first[1]] == fixture["first_queue"]
    second = quiver.next_path(first[1])
    assert second is not None
    assert repr(second[0]) == fixture["second"]
    assert [repr(path) for path in second[1]] == fixture["second_queue"]

    enumerator = quiver.enumerator()
    assert repr(enumerator.at(5)) == fixture["enumerator_5"]
    assert enumerator[4] == enumerator.at(5)


def test_next_path_fail_and_error_boundaries_follow_qpa() -> None:
    empty = qpa.Quiver([], [])
    assert empty.next_path() is None
    quiver = qpa.Quiver(["v1"], [])
    with pytest.raises(qpa.InvalidInputError, match="zero-path"):
        quiver.next_path([quiver.zero_path()])
