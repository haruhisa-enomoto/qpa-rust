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
compare_fixture \
    "$repo_root/tests/gap-reference/ordering_properties.g" \
    "$repo_root/tests/fixtures/ordering_properties.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/path_algebra.g" \
    "$repo_root/tests/fixtures/path_algebra.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/groebner_high_level.g" \
    "$repo_root/tests/fixtures/groebner_high_level.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/groebner_objects.g" \
    "$repo_root/tests/fixtures/groebner_objects.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/invlenlex.g" \
    "$repo_root/tests/fixtures/invlenlex.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/quotient_path_algebra.g" \
    "$repo_root/tests/fixtures/quotient_path_algebra.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/gbnp_np_arithmetic.g" \
    "$repo_root/tests/fixtures/gbnp_np_arithmetic.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/gbnp_obstructions.g" \
    "$repo_root/tests/fixtures/gbnp_obstructions.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/gbnp_preprocessing.g" \
    "$repo_root/tests/fixtures/gbnp_preprocessing.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/gbnp_sgrobner.g" \
    "$repo_root/tests/fixtures/gbnp_sgrobner.schema-v1.json"

# QPA's public unchecked GBNP adapter deliberately prints informational lines.
# Its fixture is the final line emitted by the reference script.
"$repo_root/scripts/gap-reference" \
    "$repo_root/tests/gap-reference/qpa_gbnp.g" |
    tail -n 1 | jq --sort-keys . >"$actual"
jq --sort-keys . \
    "$repo_root/tests/fixtures/qpa_gbnp.schema-v1.json" >"$expected"
diff --unified "$expected" "$actual"
