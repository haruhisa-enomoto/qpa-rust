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
    "$repo_root/tests/gap-reference/poset.g" \
    "$repo_root/tests/fixtures/poset.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/dynamic_tree.g" \
    "$repo_root/tests/fixtures/dynamic_tree.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/suffix_tree.g" \
    "$repo_root/tests/fixtures/suffix_tree.schema-v1.json"
if "$repo_root/scripts/gap-reference" \
    "$repo_root/tests/gap-reference/suffix_tree_empty_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted an empty suffix-tree pattern" >&2
    exit 1
fi
grep -Fq "Variable: 'currentEdge' must have an assigned value" "$actual"
if "$repo_root/scripts/gap-reference" \
    "$repo_root/tests/gap-reference/suffix_tree_duplicate_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly retained a valid duplicate-pattern suffix tree" >&2
    exit 1
fi
grep -Fq "LabelOfSuffixTreeEdge" "$actual"
grep -Fq "The 1st argument is 'fail'" "$actual"
compare_fixture \
    "$repo_root/tests/gap-reference/path_algebra.g" \
    "$repo_root/tests/fixtures/path_algebra.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_foundational.g" \
    "$repo_root/tests/fixtures/pathalg_foundational.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_conversion.g" \
    "$repo_root/tests/fixtures/pathalg_conversion.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_centre.g" \
    "$repo_root/tests/fixtures/pathalg_centre.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_vector_reduction.g" \
    "$repo_root/tests/fixtures/pathalg_vector_reduction.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_element_utilities.g" \
    "$repo_root/tests/fixtures/pathalg_element_utilities.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_vertex_quotient.g" \
    "$repo_root/tests/fixtures/pathalg_vertex_quotient.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_inverse.g" \
    "$repo_root/tests/fixtures/pathalg_inverse.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_triangular_reduced.g" \
    "$repo_root/tests/fixtures/pathalg_triangular_reduced.schema-v1.json"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_triangular_infinite_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted an infinite triangular-reduced input" >&2
    exit 1
fi
grep -Fq "The entered algebra is not finite dimensional" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_vector_empty_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly selected vector-tip dispatch for an empty list" >&2
    exit 1
fi
grep -Fq "TipMonomialandCoefficientOfVector" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_path_removal_vertex_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted a vertex term in PathRemoval" >&2
    exit 1
fi
grep -Fq "PathRemoval for an idempotent in a QuotientOfPathAlgebra" "$actual"
grep -Fq "no 1st choice method found" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/natural_map_image_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly applied its natural path-algebra projection" >&2
    exit 1
fi
grep -Fq "alghom!.generators" "$actual"
# `IsHomogeneousListOfElements([], n)` deliberately prints one source line.
"$repo_root/scripts/gap-reference" \
    "$repo_root/tests/gap-reference/pathalg_ideal.g" |
    tail -n 1 | jq --sort-keys . >"$actual"
jq --sort-keys . \
    "$repo_root/tests/fixtures/pathalg_ideal.schema-v1.json" >"$expected"
diff --unified "$expected" "$actual"
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
