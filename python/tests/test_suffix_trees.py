from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

import qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def _occurrences(values: list[tuple[int, int]]) -> list[list[int]]:
    return [list(value) for value in values]


def _node_data(tree: qpa.SuffixTree, node: qpa.SuffixTreeNode) -> dict[str, Any]:
    label = tree.label_of_node(node)
    if node.is_external:
        return {
            "kind": "external",
            "label": label,
            "occurrences": _occurrences(tree.external_occurrences(node)),
        }
    suffix = tree.node_suffix_link(node)
    return {
        "kind": "internal",
        "label": label,
        "marked": tree.node_is_marked(node),
        "suffix": None if suffix is None else tree.label_of_node(suffix),
        "edges": [
            {
                "label": tree.label_of_edge(edge),
                "child": _node_data(tree, tree.edge_child(edge)),
            }
            for edge in tree.node_edges(node)
        ],
    }


def _queries(tree: qpa.SuffixTree) -> dict[str, list[list[int]]]:
    return {
        "all_in": _occurrences(tree.all_patterns_in_sequence([0, 1, 2, 1, 2, 1, 3])),
        "containing_12": _occurrences(tree.sequence_in_patterns([1, 2])),
        "containing_21": _occurrences(tree.sequence_in_patterns([2, 1])),
        "containing_empty": _occurrences(tree.sequence_in_patterns([])),
    }


def test_suffix_tree_insert_query_and_delete_follow_qpa() -> None:
    expected = json.loads(
        (FIXTURES / "suffix_tree.schema-v1.json").read_text(encoding="utf-8")
    )
    tree = qpa.SuffixTree()
    indices = [tree.insert_pattern(pattern) for pattern in expected["patterns"]]
    assert indices == expected["indices"]
    assert tree.patterns == expected["patterns"]
    assert [tree.pattern(index) for index in indices] == expected["patterns"]
    assert _queries(tree) == expected["before_delete"]
    assert _node_data(tree, tree.root) == expected["tree_before_delete"]

    tree.delete_pattern(2)
    assert _queries(tree) == expected["after_delete_2"]
    tree.delete_pattern(1)
    assert _queries(tree) == expected["after_delete_1"]
    tree.delete_pattern(3)
    assert _queries(tree) == expected["after_delete_3"]
    assert _node_data(tree, tree.root) == expected["tree_after_delete"]
    assert not tree.is_pattern_active(1)
    assert tree.is_pattern_active(4)


def test_suffix_tree_surface_follows_qpa() -> None:
    expected = json.loads(
        (FIXTURES / "suffix_tree_surface.schema-v1.json").read_text(encoding="utf-8")
    )
    tree = qpa.SuffixTree()
    tree.insert_pattern([1, 2, 3])
    root = tree.root
    internal = tree.create_internal_node(1, 1, 1)
    external = tree.create_external_node(1, 2)
    edge = tree.create_edge(internal, external, 1, 2, None)

    split_parent = tree.create_internal_node(1, 1, 1)
    split_child = tree.create_external_node(1, 1)
    split_edge = tree.create_edge(split_parent, split_child, 1, 1, 2)
    split_node = tree.create_internal_node(1, 1, 1)
    tree.add_internal_node(split_edge, 1, split_node)
    split_remainder = tree.node_edges(split_node)[0]

    added_external = tree.create_external_node(1, 1)
    added_edge = tree.add_external_node(root, added_external)

    other_tree = qpa.SuffixTree()
    other_tree.insert_pattern([1, 2, 3])
    other_root = other_tree.root
    other_edge = other_tree.node_edges(other_root)[0]

    assert (qpa.info_dynamic_dictionary_name() == "InfoDynamicDictionary") is expected[
        "info_class"
    ]
    assert [tree.is_suffix_tree, tree.is_default_rep] == expected["tree_filters"]
    assert [
        root.is_suffix_tree_element,
        root.is_suffix_tree_node,
        root.is_internal,
        root.is_external,
        root.is_default_rep,
        root.is_default_internal_rep,
    ] == expected["root_filters"]
    assert [
        external.is_suffix_tree_element,
        external.is_suffix_tree_node,
        external.is_internal,
        external.is_external,
        external.is_default_rep,
    ] == expected["external_filters"]
    assert [
        edge.is_suffix_tree_element,
        edge.is_suffix_tree_edge,
        edge.is_default_rep,
    ] == expected["edge_filters"]
    assert {
        "dead_length": tree.dead_length,
        "pattern_count": tree.number_of_patterns,
        "pattern_1": tree.pattern(1),
    } == expected["components"]
    assert tree[0] == expected["components"]["pattern_1"]

    assert {
        "internal_label": tree.label_of_node(internal),
        "external_label": tree.label_of_node(external),
        "edge_label": tree.label_of_edge(edge),
        "root_parent_is_zero": tree.parent_node(root) is None,
        "child_parent_identity": tree.parent_node(external) == internal,
        "edge_lookup_identity": tree.edge_starting_at(internal, 2) == edge,
        "missing_lookup_is_fail": tree.edge_starting_at(internal, 9) is None,
        "self_equal": edge == edge,
        "distinct_equal": edge == split_edge,
        "other_tree_root_equal": root == other_root,
        "other_tree_edge_equal": edge == other_edge,
    } == expected["manual_edge"]
    assert {
        "prefix": tree.label_of_edge(split_edge),
        "remainder": tree.label_of_edge(split_remainder),
        "old_child_parent_identity": tree.parent_node(split_child) == split_node,
    } == expected["split"]
    assert {
        "edge_label": tree.label_of_edge(added_edge),
        "parent_identity": tree.parent_node(added_external) == root,
    } == expected["add_external"]

    with pytest.raises(qpa.InvalidInputError, match="node does not exist"):
        tree.parent_node(other_root)
    with pytest.raises(qpa.InvalidInputError, match="edge does not exist"):
        tree.label_of_edge(other_edge)


def test_suffix_tree_source_failures_follow_qpa() -> None:
    tree = qpa.SuffixTree()
    with pytest.raises(qpa.InvalidInputError, match="cannot insert an empty"):
        tree.insert_pattern([])
    with pytest.raises(qpa.InvalidInputError, match="does not exist"):
        tree.pattern(1)

    for pattern in [[1, 2, 1], [2, 1], [1, 2], [3], [1, 2, 1]]:
        tree.insert_pattern(pattern)
    tree.delete_pattern(2)
    tree.delete_pattern(1)
    with pytest.raises(qpa.InvalidInputError, match="guaranteed suffix-tree edge"):
        tree.all_patterns_in_sequence([0, 1, 2, 1, 2, 1, 3])
    with pytest.raises(qpa.InvalidInputError, match="already been deleted"):
        tree.delete_pattern(1)

    strings = qpa.SuffixTree()
    assert strings.insert_pattern("aba") == 1
    assert strings.patterns == ["aba"]
    assert strings.all_patterns_in_sequence("caba") == [(1, 2)]
    assert strings.sequence_in_patterns("ba") == [(1, 2)]
