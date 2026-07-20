from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _item(node: qpa.DynamicTree | None) -> object | None:
    return None if node is None else node.item


def _snapshot(nodes: list[qpa.DynamicTree]) -> list[list[object | None]]:
    return [
        [_item(node), _item(node.parent), _item(node.left), _item(node.right)]
        for node in nodes
    ]


def _fixture() -> dict[str, Any]:
    return json.loads(
        (FIXTURES / "dynamic_tree.schema-v1.json").read_text(encoding="utf-8")
    )


def test_dynamic_tree_construction_and_identity_follow_qpa() -> None:
    expected = _fixture()["construction"]
    forest = qpa.DynamicForest()
    empty = forest.create()
    item = forest.create("item")
    other = forest.create("item")

    assert qpa.DynamicForest.family_name == "DynamicTreeFamily"
    assert len(forest) == 3 and not forest.is_empty
    assert empty.is_dynamic_tree and empty.is_default_rep
    assert empty.has_item is expected["empty_has_item"]
    assert item.item == expected["item"]
    assert (item == item) is expected["same"]
    assert item.is_identical(item)
    assert (item == other) is expected["different"]


def test_dynamic_tree_rotations_and_splice_follow_qpa() -> None:
    expected = _fixture()
    forest = qpa.DynamicForest()

    left = forest.create("left")
    middle = forest.create("middle")
    right = forest.create("right")
    middle.left = left
    middle.right = right
    left.parent = middle
    right.parent = middle
    middle.right_rotate()
    assert _snapshot([left, middle, right]) == expected["right_rotation"]

    left = forest.create("left")
    middle = forest.create("middle")
    right = forest.create("right")
    middle.left = left
    middle.right = right
    left.parent = middle
    right.parent = middle
    middle.left_rotate()
    assert _snapshot([left, middle, right]) == expected["left_rotation"]

    parent = forest.create("parent")
    child = forest.create("child")
    child.link(parent)
    child.splice()
    assert _snapshot([parent, child]) == expected["splice"]


def test_dynamic_tree_link_root_cut_and_relink_follow_qpa() -> None:
    expected = _fixture()
    forest = qpa.DynamicForest()
    a, b, c, d = [forest.create(label) for label in ("a", "b", "c", "d")]
    a.link(b)
    b.link(c)
    d.link(c)
    nodes = [a, b, c, d]

    assert [node.root().item for node in nodes] == expected["roots_before"]
    b.cut()
    assert [node.root().item for node in nodes] == expected["roots_after_cut"]
    b.link(d)
    assert [node.root().item for node in nodes] == expected["roots_after_relink"]
    a.splay()
    assert a.root().item == "c"


def test_dynamic_tree_error_boundaries_follow_qpa() -> None:
    failures = json.loads(
        (FIXTURES / "dynamic_tree_errors.schema-v1.json").read_text(
            encoding="utf-8"
        )
    )["failures"]
    forest = qpa.DynamicForest()
    with pytest.raises(TypeError):
        forest.create(1, 2)
    leaf = forest.create("leaf")
    with pytest.raises(qpa.InvalidInputError, match="requires a left child"):
        leaf.right_rotate()
    with pytest.raises(qpa.InvalidInputError, match="requires a right child"):
        leaf.left_rotate()

    parent = forest.create("parent")
    child = forest.create("child")
    parent.left = child
    child.parent = parent
    with pytest.raises(qpa.InvalidInputError, match="not a middle child"):
        child.splice()

    foreign = qpa.DynamicForest().create("foreign")
    with pytest.raises(qpa.InvalidInputError, match="different forests"):
        child.link(foreign)
    with pytest.raises(qpa.InvalidInputError, match="different forests"):
        child.parent = foreign

    assert failures == {
        "invalid_arity_success": False,
        "right_without_left_success": False,
        "left_without_right_success": False,
        "splice_auxiliary_child_success": False,
    }
