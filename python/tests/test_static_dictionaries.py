from __future__ import annotations

import json
from pathlib import Path

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _matches(values: list[tuple[int, list[int]]]) -> list[list[object]]:
    return [[position, patterns] for position, patterns in values]


def _two_loop_quiver() -> tuple[qpa.Quiver, qpa.Path, qpa.Path]:
    quiver = qpa.Quiver(
        ["v1"], [("v1", "v1", "x"), ("v1", "v1", "y")]
    )
    return quiver, quiver.arrow_path("x"), quiver.arrow_path("y")


def test_generic_static_dictionary_follows_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "static_dictionary.schema-v1.json").read_text(encoding="utf-8")
    )["generic"]
    dictionary = qpa.StaticDictionary(fixture["patterns"])

    assert dictionary.is_static_dictionary and dictionary.is_default_rep
    assert dictionary.patterns == fixture["patterns"]
    assert _matches(dictionary.search("ushers")) == fixture["ushers"]

    overlap = qpa.StaticDictionary([[1, 2], [2], [1, 2]])
    assert _matches(overlap.search([1, 2, 1, 2])) == fixture["overlap"]

    class Symbol:
        def __init__(self, value: int) -> None:
            self.value = value

        def __eq__(self, other: object) -> bool:
            return isinstance(other, Symbol) and self.value == other.value

    symbolic = qpa.StaticDictionary([[Symbol(1), Symbol(2)]])
    assert symbolic.search([Symbol(1), Symbol(2)]) == [(2, [1])]


def test_quiver_static_dictionary_follows_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "static_dictionary.schema-v1.json").read_text(encoding="utf-8")
    )["quiver"]
    quiver, x, y = _two_loop_quiver()
    finite = qpa.QuiverStaticDictionary(quiver, [x * x, y * x, y * y])

    assert finite.is_quiver_static_dictionary and finite.is_default_rep
    assert [repr(path) for path in finite.patterns] == fixture["patterns"]
    assert _matches(finite.search(y * x * x)) == fixture["search"]
    assert _matches(finite.prefix_search(x * x)) == fixture["prefix"]
    assert finite.is_finite_difference is fixture["finite"]
    assert finite.difference_size == fixture["size"]
    words = finite.difference_words
    assert words is not None
    assert [repr(path) for path in words] == fixture["words"]

    infinite = qpa.QuiverStaticDictionary(quiver, [x * x])
    assert infinite.is_finite_difference is fixture["infinite"]
    assert infinite.difference_size == fixture["infinite_size"]
    assert (infinite.difference_words is None) is fixture["infinite_words_fail"]

    vertex_prefix = qpa.QuiverStaticDictionary(
        quiver, [quiver.vertex_path("v1"), x]
    )
    assert _matches(vertex_prefix.prefix_search(x)) == fixture["vertex_prefix"]

    overlap = qpa.QuiverStaticDictionary(quiver, [x, x * y])
    assert overlap.is_finite_difference is fixture["overlap_finite"]
    assert overlap.difference_size == fixture["overlap_size"]


def test_static_dictionary_error_boundaries_follow_qpa() -> None:
    failures = json.loads(
        (FIXTURES / "static_dictionary_errors.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )["failures"]
    generic = qpa.StaticDictionary([[1]])
    with pytest.raises(qpa.InvalidInputError, match="empty dense list"):
        generic.search([])
    assert failures["empty_search_success"] is False

    quiver, x, _y = _two_loop_quiver()
    foreign, foreign_x, _foreign_y = _two_loop_quiver()
    with pytest.raises(qpa.InvalidInputError, match="different quiver"):
        qpa.QuiverStaticDictionary(quiver, [foreign_x])
    with pytest.raises(qpa.InvalidInputError, match="zero path"):
        qpa.QuiverStaticDictionary(quiver, [quiver.zero_path()])
    dictionary = qpa.QuiverStaticDictionary(quiver, [x])
    with pytest.raises(qpa.InvalidInputError, match="different quiver"):
        dictionary.search(foreign.arrow_path("x"))
    with pytest.raises(qpa.InvalidInputError, match="zero path"):
        dictionary.search(quiver.zero_path())
    assert failures["foreign_pattern_success"] is False
    assert failures["zero_pattern_success"] is False
    assert failures["foreign_search_success"] is False
    assert failures["zero_search_success"] is False


def test_static_dictionary_output_and_info_follow_qpa() -> None:
    quiver, x, _y = _two_loop_quiver()
    actual = [
        f"GENERIC_0={qpa.StaticDictionary([])!r}",
        f"GENERIC_1={qpa.StaticDictionary([[1]])!r}",
        f"GENERIC_2={qpa.StaticDictionary([[1], [2]])!r}",
        f"QUIVER_0={qpa.QuiverStaticDictionary(quiver, [])!r}",
        f"QUIVER_1={qpa.QuiverStaticDictionary(quiver, [x])!r}",
        f"QUIVER_2={qpa.QuiverStaticDictionary(quiver, [x, x * x])!r}",
    ]
    expected = (FIXTURES / "static_dictionary_output.txt").read_text(
        encoding="utf-8"
    )
    assert "\n".join(actual) + "\n" == expected

    messages = qpa.static_dictionary_info_messages()
    expected_info = (FIXTURES / "static_dictionary_info.txt").read_text(
        encoding="utf-8"
    )
    assert "".join(messages) == expected_info
