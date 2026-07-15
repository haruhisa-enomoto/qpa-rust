#!/bin/sh

set -eu

script_dir=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH='' cd -- "$script_dir/../.." && pwd)

actual=$(mktemp "${TMPDIR:-/tmp}/qpa-rust-quiver-actual.XXXXXX")
expected=$(mktemp "${TMPDIR:-/tmp}/qpa-rust-quiver-expected.XXXXXX")
cleanup() {
    rm -f -- "$actual" "$expected"
}
trap cleanup EXIT HUP INT TERM

compare_fixture() {
    source_script=$1
    fixture=$2

    "$repo_root/scripts/gap-reference" "$source_script" |
        jq --sort-keys . >"$actual"
    jq --sort-keys . "$fixture" >"$expected"
    diff --unified "$expected" "$actual"
}

compare_fixture \
    "$repo_root/tests/gap-reference/quiver_paths.g" \
    "$repo_root/tests/fixtures/quiver_paths.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/quiver_graph.g" \
    "$repo_root/tests/fixtures/quiver_graph.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/orderings.g" \
    "$repo_root/tests/fixtures/orderings.schema-v1.json"
