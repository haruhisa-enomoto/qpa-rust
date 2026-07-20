from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _values(values: qpa.InfiniteList | qpa.HalfInfiniteList, positions: list[int]) -> list[Any]:
    return [values.value(position) for position in positions]


def test_half_infinite_lists_follow_qpa() -> None:
    expected = json.loads(
        (FIXTURES / "inflist.schema-v1.json").read_text(encoding="utf-8")
    )["half_lists"]
    positive = qpa.HalfInfiniteList(2, 1, ["repeat", [10, 20, 30]])
    negative = qpa.HalfInfiniteList(-1, -1, ["repeat", [7, 8]])

    assert positive.is_half_infinite_list and positive.is_default_rep
    assert {
        "start": positive.start_position,
        "direction": positive.direction,
        "type": positive.list_type,
        "repeating": positive.repeating_list,
        "storing": positive.is_storing_values,
        "initial_is_fail": positive.initial_value is None,
        "values": _values(positive, list(range(2, 9))),
        "lowest": positive.lowest_known_position,
        "highest": positive.highest_known_position,
    } == expected["positive_repeat"]
    assert {
        "values": _values(negative, list(range(-6, 0))),
        "lowest": negative.lowest_known_position,
        "highest": negative.highest_known_position,
    } == expected["negative_repeat"]
    assert positive.element_function is None

    next_callbacks: list[list[object]] = []
    repeat_data: list[object] = []

    def next_callback(position: int, direction: int, list_type: str) -> None:
        next_callbacks.append([position, direction, list_type])

    def repeat_callback(tail: list[object], repeating: qpa.HalfInfiniteList) -> None:
        repeat_data.extend(
            [
                tail,
                repeating.start_position,
                repeating.direction,
                repeating.repeating_list,
            ]
        )

    next_list = qpa.HalfInfiniteList(
        0,
        1,
        ["next/repeat", lambda value: (value + 1) % 4, -1],
        next_callback,
        repeat_callback,
    )
    assert {
        "value_10": next_list.value(10),
        "stored": next_list.stored_values,
        "callbacks": next_callbacks,
        "repeat_data": repeat_data,
    } == expected["next_repeat"]
    element_function = next_list.element_function
    assert element_function is not None and element_function(3) == 0
    callback = next_list.new_value_callback
    assert callback is not None

    position_callbacks: list[list[object]] = []

    def position_callback(position: int, direction: int, list_type: str) -> None:
        position_callbacks.append([position, direction, list_type])

    position_list = qpa.HalfInfiniteList(
        -2, 1, ["pos", lambda value: value**2, True], position_callback
    )
    lowest_before = position_list.lowest_known_position
    assert {
        "lowest_before": lowest_before,
        "value_2": position_list.value(2),
        "value_0_again": position_list.value(0),
        "stored": position_list.stored_values,
        "callbacks": position_callbacks,
        "lowest_after": position_list.lowest_known_position,
        "highest_after": position_list.highest_known_position,
    } == expected["position_storing"]

    nonstoring_callbacks: list[list[object]] = []
    nonstoring = qpa.HalfInfiniteList(
        5,
        -1,
        ["pos", lambda value: value + 100, False],
        lambda position, direction, list_type: nonstoring_callbacks.append(
            [position, direction, list_type]
        ),
    )
    assert {
        "value_1": nonstoring.value(1),
        "stored": nonstoring.stored_values,
        "callbacks": nonstoring_callbacks,
        "lowest": nonstoring.lowest_known_position,
    } == expected["position_nonstoring"]


def test_full_infinite_lists_and_transformations_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "inflist.schema-v1.json").read_text(encoding="utf-8")
    )
    full = qpa.InfiniteList(
        0, [100, 200], ["repeat", [1, 2, 3]], ["repeat", [7, 8]]
    )
    expected = fixture["full"]
    assert full.is_infinite_list and full.is_default_rep
    assert {
        "middle_start": full.middle_start,
        "middle_end": full.middle_end,
        "middle": full.middle_part,
        "positive_start": full.positive_part.start_position,
        "negative_start": full.negative_part.start_position,
        "values": _values(full, [-5, -2, -1, 0, 1, 2, 5]),
        "lowest_known": full.lowest_known_position,
        "highest_known": full.highest_known_position,
        "lower_bound": full.lower_bound,
        "upper_bound": full.upper_bound,
    } == expected

    transformations = fixture["transformations"]
    shifted = full.shift(3)
    assert {
        "middle_start": shifted.middle_start,
        "values": _values(shifted, [-5, -3, -2, -1, 2]),
    } == transformations["shift"]
    assert _values(full.map(lambda value: 2 * value), [-2, -1, 0, 1, 2]) == transformations[
        "map"
    ]

    positive_cut = full.positive_part.cut(3)
    assert {
        "middle_start": positive_cut.middle_start,
        "middle": positive_cut.middle_part,
        "positive_start": positive_cut.positive_part.start_position,
        "values": _values(positive_cut, list(range(3, 9))),
    } == transformations["positive_cut"]
    negative_cut = full.negative_part.cut(-2)
    assert {
        "middle_start": negative_cut.middle_start,
        "middle": negative_cut.middle_part,
        "negative_start": negative_cut.negative_part.start_position,
        "values": _values(negative_cut, list(range(-6, -1))),
    } == transformations["negative_cut"]
    assert _values(full.positive_part_from(3), list(range(3, 9))) == transformations[
        "positive_from_values"
    ]
    assert _values(full.negative_part_from(-2), list(range(-6, -1))) == transformations[
        "negative_from_values"
    ]

    right = qpa.InfiniteList(
        0, [20, 21], ["repeat", [22, 23]], ["repeat", [18, 19]]
    )
    assert _values(qpa.InfiniteList.splice(full, right, 0), list(range(-2, 4))) == transformations[
        "splice"
    ]
    concatenated = qpa.InfiniteList.concatenate([full, right])
    assert {
        "middle_start": concatenated.middle_start,
        "middle": concatenated.middle_part,
        "positive_start": concatenated.positive_part.start_position,
        "negative_start": concatenated.negative_part.start_position,
        "values": _values(concatenated, list(range(-4, 4))),
    } == transformations["concatenation"]
    half_mapped = qpa.HalfInfiniteList(
        2, 1, ["repeat", [10, 20, 30]]
    ).map(lambda value: value + 1)
    assert _values(half_mapped, list(range(2, 8))) == transformations["half_map"]


def test_infinite_list_constructors_and_repeat_conversion_follow_qpa() -> None:
    fixture = json.loads(
        (FIXTURES / "inflist.schema-v1.json").read_text(encoding="utf-8")
    )
    constructors = fixture["constructors"]
    finite = qpa.InfiniteList.finite(-1, [4, 5, 6])
    assert {
        "middle_start": finite.middle_start,
        "middle_end": finite.middle_end,
        "values": _values(finite, [-1, 0, 1]),
        "lower_bound": finite.lower_bound,
        "upper_bound": finite.upper_bound,
    } == constructors["finite"]
    assert _values(qpa.InfiniteList.constant(9), [-3, 0, 4]) == constructors["constant"]
    assert _values(qpa.integers_list(), [-3, 0, 4]) == constructors["integers"]
    function_list = qpa.InfiniteList.function(lambda position: position * position)
    assert _values(function_list, [-2, 3]) == [4, 9]

    def tail_next(value: int) -> int:
        return {-1: 9, 9: 5, 5: 6, 6: 7}.get(value, 5)

    positive = qpa.InfiniteList(
        0, [99], ["next/repeat", tail_next, -1], ["repeat", [8]]
    )
    positive_value = positive.value(10)
    expected = fixture["repeat_conversion"]["positive_parent"]
    assert {
        "value_10": positive_value,
        "middle_start": positive.middle_start,
        "middle": positive.middle_part,
        "positive_start": positive.positive_part.start_position,
        "period": positive.positive_part.repeating_list,
    } == expected
    negative = qpa.InfiniteList(
        0, [99], ["repeat", [8]], ["next/repeat", tail_next, -1]
    )
    negative_value = negative.value(-10)
    expected = fixture["repeat_conversion"]["negative_parent"]
    assert {
        "value_minus_10": negative_value,
        "middle_start": negative.middle_start,
        "middle": negative.middle_part,
        "negative_start": negative.negative_part.start_position,
        "period": negative.negative_part.repeating_list,
    } == expected


def test_infinite_list_source_failures_follow_qpa() -> None:
    expected = json.loads(
        (FIXTURES / "inflist_errors.schema-v1.json").read_text(encoding="utf-8")
    )["success"]
    finite = qpa.InfiniteList.finite(0, [1, 2])
    full = qpa.InfiniteList.constant(1)
    positive = qpa.HalfInfiniteList(1, 1, ["repeat", [1]])
    bad_positive = qpa.HalfInfiniteList(2, 1, ["repeat", [1]])
    negative = qpa.HalfInfiniteList(-1, -1, ["repeat", [1]])
    empty_repeat = qpa.HalfInfiniteList(0, 1, ["repeat", []])

    operations: dict[str, Callable[[], object]] = {
        "shift_finite": lambda: finite.shift(1),
        "map_finite": lambda: finite.map(lambda value: 2 * value),
        "empty_known": lambda: qpa.InfiniteList.finite(0, []).lowest_known_position,
        "shifted_splice": lambda: qpa.InfiniteList.shifted_splice(full, 0, full, 0, 0),
        "illegal_half_position": lambda: positive.value(0),
        "empty_repeat_value": lambda: empty_repeat.value(0),
        "cut_outside": lambda: positive.cut(0),
        "finite_outside": lambda: finite.value(2),
        "bad_positive_start": lambda: qpa.InfiniteList.from_half_lists(
            0, [1], bad_positive, None
        ),
        "bad_negative_start": lambda: qpa.InfiniteList.from_half_lists(
            1, [1], None, negative
        ),
    }
    for name, operation in operations.items():
        with pytest.raises(qpa.InvalidInputError):
            operation()
        assert expected[name] is False
