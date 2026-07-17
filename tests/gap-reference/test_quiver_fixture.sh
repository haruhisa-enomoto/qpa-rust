#!/bin/sh

set -eu

script_dir=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)
repo_root=$(CDPATH='' cd -- "$script_dir/../.." && pwd)

actual=$(mktemp "${TMPDIR:-/tmp}/qpa-rust-quiver-actual.XXXXXX")
expected=$(mktemp "${TMPDIR:-/tmp}/qpa-rust-quiver-expected.XXXXXX")
raw=$(mktemp "${TMPDIR:-/tmp}/qpa-rust-quiver-raw.XXXXXX")
cleanup() {
    rm -f -- "$actual" "$expected" "$raw"
}

compare_last_json_fixture() {
    source_script=$1
    fixture=$2

    "$repo_root/scripts/gap-reference" "$source_script" >"$raw"
    grep -Fq "Computing step 1..." "$raw"
    sed -n '/^{"schema_version":/,$p' "$raw" |
        jq --sort-keys . >"$actual"
    jq --sort-keys . "$fixture" >"$expected"
    diff --unified "$expected" "$actual"
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

compare_catching_fixture() {
    source_script=$1
    fixture=$2

    "$repo_root/scripts/gap-reference" -T "$source_script" >"$raw"
    jq --sort-keys . "$raw" >"$actual"
    jq --sort-keys . "$fixture" >"$expected"
    diff --unified "$expected" "$actual"
}

compare_prefixed_catching_fixture() {
    source_script=$1
    fixture=$2
    prefix=$3
    expected_count=$4

    "$repo_root/scripts/gap-reference" -T "$source_script" >"$raw"
    actual_count=$(grep -Fc "$prefix" "$raw")
    test "$actual_count" -eq "$expected_count"
    sed -n '/^{"schema_version":/,$p' "$raw" |
        jq --sort-keys . >"$actual"
    jq --sort-keys . "$fixture" >"$expected"
    diff --unified "$expected" "$actual"
}

compare_cartan_coxeter_dispatch_fixture() {
    source_script=$1
    fixture=$2

    "$repo_root/scripts/gap-reference" -T "$source_script" >"$raw"
    test "$(grep -Fc 'Algebra is not finite dimensional!' "$raw")" -eq 1
    test "$(grep -Fc 'Not a bound quiver algebra!' "$raw")" -eq 1
    test "$(grep -Fc 'Unable to determine the Cartan matrix!' "$raw")" -eq 2
    test "$(grep -Fc 'The Cartan matrix is not invertible!' "$raw")" -eq 1
    sed -n '/^{"schema_version":/,$p' "$raw" |
        jq --sort-keys . >"$actual"
    jq --sort-keys . "$fixture" >"$expected"
    diff --unified "$expected" "$actual"
}

compare_reptype_quotient_fixture() {
    source_script=$1
    fixture=$2

    "$repo_root/scripts/gap-reference" -T "$source_script" >"$raw" 2>&1
    test "$(grep -Fc 'Finite type!' "$raw")" -eq 3
    test "$(grep -Fc 'Semisimple algebra (Q has no arrows).' "$raw")" -eq 1
    test "$(grep -Fc 'Quiver is a (union of) Dynkin quiver(s).' "$raw")" -eq 1
    test "$(grep -Fc 'Quiver contains an (un)oriented cycle.' "$raw")" -eq 5
    test "$(grep -Fc 'A_3' "$raw")" -eq 5
    test "$(grep -Fc 'A_2' "$raw")" -eq 4
    test "$(grep -Fc 'A_1' "$raw")" -eq 8
    test "$(grep -Fc 'Have checked representation type of the algebra modulo the square of the radical.' "$raw")" -eq 1
    test "$(grep -Fc 'Special biserial algebra with no unoriented cycles in Q.' "$raw")" -eq 1
    bongartz_count=$(grep -Fc "Entering Bongartz' test for infinite type." "$raw")
    test "$bongartz_count" -eq 1
    test "$(grep -Fc 'Can not determine the representation type.' "$raw")" -eq 1
    test "$(grep -Fc 'not an admissible quotient of a path algebra' "$raw")" -eq 1
    sed -n '/^{"schema_version":/,$p' "$raw" |
        jq --sort-keys . >"$actual"
    jq --sort-keys . "$fixture" >"$expected"
    diff --unified "$expected" "$actual"
}

compare_text_fixture() {
    source_script=$1
    fixture=$2

    "$repo_root/scripts/gap-reference" "$source_script" >"$actual"
    diff --unified "$fixture" "$actual"
}

compare_fixture \
    "$repo_root/tests/gap-reference/quiver_paths.g" \
    "$repo_root/tests/fixtures/quiver_paths.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/quiver_surface.g" \
    "$repo_root/tests/fixtures/quiver_surface.schema-v1.json"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/quiver_dispatch.g" \
    "$repo_root/tests/fixtures/quiver_dispatch.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/util.g" \
    "$repo_root/tests/fixtures/util.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/algebra_structure_constants.g" \
    "$repo_root/tests/fixtures/algebra_structure_constants.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/moduledecomp_foundation.g" \
    "$repo_root/tests/fixtures/moduledecomp_foundation.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/moduledecomp_modules.g" \
    "$repo_root/tests/fixtures/moduledecomp_modules.schema-v1.json"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/idempotent.g" \
    "$repo_root/tests/fixtures/idempotent.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/algebra_opposite.g" \
    "$repo_root/tests/fixtures/algebra_opposite.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/algebra_mapping_radical.g" \
    "$repo_root/tests/fixtures/algebra_mapping_radical.schema-v1.json"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/algebra_mapping_radical_f2_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly computed a GF(2) linear-mapping radical" >&2
    exit 1
fi
grep -Fq "no 1st choice method found for \`IsZero'" "$actual"
grep -Fq "The 1st argument is 'fail'" "$actual"
"$repo_root/scripts/gap-reference" \
    "$repo_root/tests/gap-reference/algebra_opposite_print.g" \
    >"$actual" 2>&1
grep -Fq "OppositeAlgebra( OppositeAlgebra(" "$actual"
grep -Fq "printing stopped, too many recursion levels!" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/algebra_opposite_dot_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted an absent opposite-algebra component" >&2
    exit 1
fi
grep -Fq "illegal access to record component \`obj.missing'" "$actual"
grep -Fq "algebra.gi:198" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/algebra_structure_constants_empty_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly found an identity for a zero-dimensional structure algebra" >&2
    exit 1
fi
grep -Fq "need <rpb>, <cpb>, <zero> if <blocks> is empty" "$actual"
grep -Fq "algebra.gi:622" "$actual"
compare_fixture \
    "$repo_root/tests/gap-reference/quiver_graph.g" \
    "$repo_root/tests/fixtures/quiver_graph.schema-v1.json"
compare_text_fixture \
    "$repo_root/tests/gap-reference/quiver_dynkin_output.g" \
    "$repo_root/tests/fixtures/quiver_dynkin_output.txt"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/quiver_graph_errors.g" \
    "$repo_root/tests/fixtures/quiver_graph_errors.schema-v1.json"
compare_text_fixture \
    "$repo_root/tests/gap-reference/poset_output.g" \
    "$repo_root/tests/fixtures/poset_output.txt"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/poset_errors.g" \
    "$repo_root/tests/fixtures/poset_errors.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalgtensor_quiver.g" \
    "$repo_root/tests/fixtures/pathalgtensor_quiver.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalgtensor_elements.g" \
    "$repo_root/tests/fixtures/pathalgtensor_elements.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalgtensor_simple.g" \
    "$repo_root/tests/fixtures/pathalgtensor_simple.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalgtensor_full.g" \
    "$repo_root/tests/fixtures/pathalgtensor_full.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalgtensor_inclusion.g" \
    "$repo_root/tests/fixtures/pathalgtensor_inclusion.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalgtensor_enveloping.g" \
    "$repo_root/tests/fixtures/pathalgtensor_enveloping.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalgtensor_regular_module.g" \
    "$repo_root/tests/fixtures/pathalgtensor_regular_module.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalgtensor_dual_module.g" \
    "$repo_root/tests/fixtures/pathalgtensor_dual_module.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalgtensor_trivial_extension.g" \
    "$repo_root/tests/fixtures/pathalgtensor_trivial_extension.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_foundation.g" \
    "$repo_root/tests/fixtures/modulehom_foundation.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_functors.g" \
    "$repo_root/tests/fixtures/modulehom_functors.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_categorical.g" \
    "$repo_root/tests/fixtures/modulehom_categorical.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_projectives.g" \
    "$repo_root/tests/fixtures/modulehom_projectives.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_projective_matrices.g" \
    "$repo_root/tests/fixtures/modulehom_projective_matrices.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_syzygy_injective.g" \
    "$repo_root/tests/fixtures/modulehom_syzygy_injective.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_global_dimension.g" \
    "$repo_root/tests/fixtures/modulehom_global_dimension.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/algebra_homological_properties.g" \
    "$repo_root/tests/fixtures/algebra_homological_properties.schema-v1.json"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/algebra_homological_dispatch.g" \
    "$repo_root/tests/fixtures/algebra_homological_dispatch.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_symmetry_properties.g" \
    "$repo_root/tests/fixtures/pathalg_symmetry_properties.schema-v1.json"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/pathalg_symmetry_dispatch.g" \
    "$repo_root/tests/fixtures/pathalg_symmetry_dispatch.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_cartan_coxeter.g" \
    "$repo_root/tests/fixtures/pathalg_cartan_coxeter.schema-v1.json"
compare_cartan_coxeter_dispatch_fixture \
    "$repo_root/tests/gap-reference/pathalg_cartan_coxeter_dispatch.g" \
    "$repo_root/tests/fixtures/pathalg_cartan_coxeter_dispatch.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_nakayama_permutation.g" \
    "$repo_root/tests/fixtures/pathalg_nakayama_permutation.schema-v1.json"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/pathalg_nakayama_permutation_dispatch.g" \
    "$repo_root/tests/fixtures/pathalg_nakayama_permutation_dispatch.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_nakayama_automorphism.g" \
    "$repo_root/tests/fixtures/pathalg_nakayama_automorphism.schema-v1.json"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/pathalg_nakayama_automorphism_dispatch.g" \
    "$repo_root/tests/fixtures/pathalg_nakayama_automorphism_dispatch.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/module_families_properties.g" \
    "$repo_root/tests/fixtures/module_families_properties.schema-v1.json"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/module_random.g" \
    "$repo_root/tests/fixtures/module_random.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/module_series.g" \
    "$repo_root/tests/fixtures/module_series.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/module_annihilator.g" \
    "$repo_root/tests/fixtures/module_annihilator.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_ext.g" \
    "$repo_root/tests/fixtures/modulehom_ext.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/functors_star_nakayama.g" \
    "$repo_root/tests/fixtures/functors_star_nakayama.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/functors_transpose.g" \
    "$repo_root/tests/fixtures/functors_transpose.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/functors_ar_maps.g" \
    "$repo_root/tests/fixtures/functors_ar_maps.schema-v1.json"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/functors_tensor_modules.g" \
    "$repo_root/tests/fixtures/functors_tensor_modules.schema-v1.json"
compare_last_json_fixture \
    "$repo_root/tests/gap-reference/functors_iterations.g" \
    "$repo_root/tests/fixtures/functors_iterations.schema-v1.json"
compare_last_json_fixture \
    "$repo_root/tests/gap-reference/functors_ar_map_iterations.g" \
    "$repo_root/tests/fixtures/functors_ar_map_iterations.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_stable.g" \
    "$repo_root/tests/fixtures/modulehom_stable.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_split.g" \
    "$repo_root/tests/fixtures/modulehom_split.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_minimal.g" \
    "$repo_root/tests/fixtures/modulehom_minimal.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_common_summand.g" \
    "$repo_root/tests/fixtures/modulehom_common_summand.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_map_utils.g" \
    "$repo_root/tests/fixtures/modulehom_map_utils.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_isomorphism.g" \
    "$repo_root/tests/fixtures/modulehom_isomorphism.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_end_quiver.g" \
    "$repo_root/tests/fixtures/modulehom_end_quiver.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/modulehom_submodules.g" \
    "$repo_root/tests/fixtures/modulehom_submodules.schema-v1.json"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_all_simple_rational_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly enumerated simple submodules over the rationals" >&2
    exit 1
fi
grep -Fq "not a module over an algebra over a finite field" "$actual"
grep -Fq "modulehom.gi:3754" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_all_submodules_rational_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly enumerated all submodules over the rationals" >&2
    exit 1
fi
grep -Fq "not a module over an algebra over a finite field" "$actual"
grep -Fq "modulehom.gi:3797" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_end_quiver_nonbasic_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly converted a nonbasic module endomorphism algebra" >&2
    exit 1
fi
grep -Fq "the entered algebra is not basic" "$actual"
grep -Fq "algebra.gi:691" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_end_quiver_simple_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly converted a simple module endomorphism algebra" >&2
    exit 1
fi
grep -Fq "Inner product multiplication of lists: no summands" "$actual"
grep -Fq "algebra.gi:811" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_inverse_noniso_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly inverted a non-isomorphism" >&2
    exit 1
fi
grep -Fq "entered map is not an isomorphism" "$actual"
grep -Fq "modulehom.gi:3330" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_images_support_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted a basis image with wrong support" >&2
    exit 1
fi
grep -Fq "does not give a homomorphism -- wrong support" "$actual"
grep -Fq "modulehom.gi:3376" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_projective_support_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted a projective matrix with discrepant support" >&2
    exit 1
fi
grep -Fq "Discrepancy between entered list of vertices and support" "$actual"
grep -Fq "modulehom.gi:4318" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/module_simple_infinite_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly constructed simples for an infinite algebra" >&2
    exit 1
fi
grep -Fq "argument entered is not a finite dimensional algebra" "$actual"
grep -Fq "modulespecial.gi:173" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/module_simple_nonadmissible_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly constructed simples for a nonadmissible quotient" >&2
    exit 1
fi
grep -Fq "argument entered is not a quotient of a path algebra by an admissible ideal" "$actual"
grep -Fq "modulespecial.gi:176" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_projective_zero_element_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly converted a zero algebra element to a projective element" >&2
    exit 1
fi
grep -Fq "The entered element is zero" "$actual"
grep -Fq "modulespecial.gi:403" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_projective_scaled_tip_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly converted a projective element with a scaled tip" >&2
    exit 1
fi
grep -Fq "no 1st choice method found for \`+' on 2 arguments" "$actual"
grep -Fq "The 2nd argument is 'fail'" "$actual"
grep -Fq "modulespecial.gi:428" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_minimal_zero_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted MoreLeftMinimalVersion of a zero endomorphism" >&2
    exit 1
fi
grep -Fq "no 1st choice method found for \`*' on 2 arguments" "$actual"
grep -Fq "modulehom.gi:3174" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_stable_zero_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted stable End of the zero module" >&2
    exit 1
fi
grep -Fq "zero module given as input for EndModuleProjOverAlgebra" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_trace_different_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly traced modules over different algebras" >&2
    exit 1
fi
grep -Fq "the entered modules are not modules over the same algebra" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/modulehom_end_zero_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted EndOverAlgebra of the zero module" >&2
    exit 1
fi
grep -Fq "<list> must be nonempty" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalgtensor_empty_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted an empty tensor-product presentation" >&2
    exit 1
fi
grep -Fq "first argument must be a ring" "$actual"
grep -Fq "second argument must be a list of generators" "$actual"
compare_fixture \
    "$repo_root/tests/gap-reference/orderings.g" \
    "$repo_root/tests/fixtures/orderings.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/ordering_properties.g" \
    "$repo_root/tests/fixtures/ordering_properties.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/ordering_surface.g" \
    "$repo_root/tests/fixtures/ordering_surface.schema-v1.json"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/ordering_dispatch.g" \
    "$repo_root/tests/fixtures/ordering_dispatch.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/poset.g" \
    "$repo_root/tests/fixtures/poset.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/dynamic_tree.g" \
    "$repo_root/tests/fixtures/dynamic_tree.schema-v1.json"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/dynamic_tree_errors.g" \
    "$repo_root/tests/fixtures/dynamic_tree_errors.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/static_dictionary.g" \
    "$repo_root/tests/fixtures/static_dictionary.schema-v1.json"
compare_text_fixture \
    "$repo_root/tests/gap-reference/static_dictionary_output.g" \
    "$repo_root/tests/fixtures/static_dictionary_output.txt"
compare_text_fixture \
    "$repo_root/tests/gap-reference/static_dictionary_info.g" \
    "$repo_root/tests/fixtures/static_dictionary_info.txt"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/static_dictionary_errors.g" \
    "$repo_root/tests/fixtures/static_dictionary_errors.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/suffix_tree.g" \
    "$repo_root/tests/fixtures/suffix_tree.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/suffix_tree_surface.g" \
    "$repo_root/tests/fixtures/suffix_tree_surface.schema-v1.json"
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
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_associated.g" \
    "$repo_root/tests/fixtures/pathalg_associated.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_admissible_properties.g" \
    "$repo_root/tests/fixtures/pathalg_admissible_properties.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_persistence.g" \
    "$repo_root/tests/fixtures/pathalg_persistence.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_prime_field_logs.g" \
    "$repo_root/tests/fixtures/pathalg_prime_field_logs.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_special_biserial.g" \
    "$repo_root/tests/fixtures/pathalg_special_biserial.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_predefined.g" \
    "$repo_root/tests/fixtures/pathalg_predefined.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_poset_algebra.g" \
    "$repo_root/tests/fixtures/pathalg_poset_algebra.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_brauer_configuration.g" \
    "$repo_root/tests/fixtures/pathalg_brauer_configuration.schema-v1.json"
compare_fixture \
    "$repo_root/tests/gap-reference/pathalg_predefined_remaining.g" \
    "$repo_root/tests/fixtures/pathalg_predefined_remaining.schema-v1.json"
compare_last_json_fixture \
    "$repo_root/tests/gap-reference/pathalg_module_preprojective_structure.g" \
    "$repo_root/tests/fixtures/pathalg_module_preprojective_structure.schema-v1.json"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/algebra_as_quiver_decomposable_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly converted a decomposable structure algebra" >&2
    exit 1
fi
grep -Fq "the entered algebra is not indecomposable" "$actual"
grep -Fq "algebra.gi:674" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/algebra_as_quiver_nonbasic_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly converted a nonbasic matrix algebra" >&2
    exit 1
fi
grep -Fq "the entered algebra is not basic" "$actual"
grep -Fq "algebra.gi:691" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_module_preprojective_bound_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted an insufficient module-preprojective bound" >&2
    exit 1
fi
grep -Fq "Must increase the second argument" "$actual"
grep -Fq "pathalgpredef.gi:852" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_module_preprojective_rational_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly converted a module-preprojective algebra over the rationals" >&2
    exit 1
fi
grep -Fq "the entered algebra is not an algebra over a finite field" "$actual"
grep -Fq "algebra.gi:653" "$actual"
compare_text_fixture \
    "$repo_root/tests/gap-reference/pathalgreptype_path.g" \
    "$repo_root/tests/fixtures/pathalgreptype_path.txt"
compare_prefixed_catching_fixture \
    "$repo_root/tests/gap-reference/reptype_bongartz.g" \
    "$repo_root/tests/fixtures/reptype_bongartz.schema-v1.json" \
    "Entering Bongartz' test for infinite type." \
    4
compare_reptype_quotient_fixture \
    "$repo_root/tests/gap-reference/reptype_quotient.g" \
    "$repo_root/tests/fixtures/reptype_quotient.schema-v1.json"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_preprojective_cycle_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted a cyclic preprojective quiver" >&2
    exit 1
fi
grep -Fq "The entered quiver has an oriented cycle." "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_preprojective_isolated_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted all-zero preprojective relations" >&2
    exit 1
fi
grep -Fq "no 1st choice method found for" "$actual"
grep -Fq "GroebnerBasis" "$actual"
"$repo_root/scripts/gap-reference" \
    "$repo_root/tests/gap-reference/pathalg_brauer_configuration_invalid.g" \
    >"$actual"
for diagnostic in \
    "Brauer Algebra input must contain 3 arrays: [vertices], [edges or polygons], [orientations]." \
    "There must be at least 1 vertex" \
    "Vertices should be of the form: [name, multiplicity]." \
    "Multiplicities must be positive integers." \
    "Each vertex must have a unique name." \
    "Edges or Polygons should be lists of length 2 or greater." \
    "An edge or polygon cannot have the same name as a vertex." \
    "An edge or polygon cannot contain a vertex which was not listed." \
    "Edges or polygons must have unique names." \
    "There must be an orientation corresponding to each vertex." \
    "Orientations must be lists of length 1 or greater." \
    "Orientations may only contain edges or polygons which were listed." \
    "Orientations must contain exactly the edges or polygons which contain the corresponding vertex."
do
    grep -Fq "$diagnostic" "$actual"
done
for label in \
    OUTER_SHAPE NO_VERTICES VERTEX_SHAPE MULTIPLICITY DUPLICATE_VERTEX \
    EMPTY_POLYGON POLYGON_VERTEX_NAME UNKNOWN_VERTEX DUPLICATE_POLYGON \
    ORIENTATION_COUNT EMPTY_ORIENTATION UNKNOWN_POLYGON CONTAINMENT
do
    grep -Fq "FAIL_$label=true" "$actual"
done
grep -Fq 'MUTATED_ORIENTATIONS=[ [ "z" ], [ "x", "z" ] ]' "$actual"
"$repo_root/scripts/gap-reference" \
    "$repo_root/tests/gap-reference/pathalg_predefined_invalid.g" \
    >"$actual"
grep -Fq "The admissible sequence is not a valid one." "$actual"
grep -Fq "INVALID_RETURN=[ 1, 1 ]" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_predefined_empty_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted an empty Nakayama admissible sequence" >&2
    exit 1
fi
grep -Fq "List Element: <list>[1] must have an assigned value" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_predefined_nonpositive_singleton_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted a nonpositive singleton Nakayama sequence" >&2
    exit 1
fi
grep -Fq "each summand of each element is not (only) a constant" "$actual"
grep -Fq "no 1st choice method found for" "$actual"
grep -Fq "GroebnerBasis" "$actual"
"$repo_root/scripts/gap-reference" \
    "$repo_root/tests/gap-reference/pathalg_save_keep.g" \
    >"$actual"
grep -Fq "A file with the same name exists already." "$actual"
grep -Fq "KEEP_RESULT=false" "$actual"
"$repo_root/scripts/gap-reference" \
    "$repo_root/tests/gap-reference/pathalg_read_prefix_empty.g" \
    >"$actual"
grep -Fq "PREFIX_EMPTY_IS_QUOTIENT=false" "$actual"
grep -Fq "PREFIX_EMPTY_RELATIONS=0" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_save_query_new_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly initialized the query-new output stream" >&2
    exit 1
fi
grep -Fq "Variable: 'output' must have an assigned value" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_save_invalid_option_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted an invalid SaveAlgebra option" >&2
    exit 1
fi
grep -Fq "Invalid overwrite option has been entered" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_triangular_infinite_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted an infinite triangular-reduced input" >&2
    exit 1
fi
grep -Fq "The entered algebra is not finite dimensional" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_associated_infinite_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted an infinite associated-graded input" >&2
    exit 1
fi
grep -Fq "The entered algebra is not finite dimensional" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_associated_nonadmissible_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted a nonadmissible associated-graded input" >&2
    exit 1
fi
grep -Fq "The entered algebra is not an admissible quotient of a path algebra" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_associated_cached_false.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly rejected cached-false admissibility before its inverse-order backend" >&2
    exit 1
fi
grep -Fq "IsAlgebraModule" "$actual"
grep -Fq "LoewyLength" "$actual"
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
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/groebner_high_level_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted a vertex relation in HighLevelGroebnerBasis" >&2
    exit 1
fi
grep -Fq "elements do not belong to the arrow ideal of the path algebra" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/groebner_high_level_foreign_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted a relation from another path algebra" >&2
    exit 1
fi
grep -Fq "elements do not belong to the arrow ideal of the path algebra" "$actual"
"$repo_root/scripts/gap-reference" \
    "$repo_root/tests/gap-reference/pathalg_koszul.g" >"$raw"
grep -Fq "IsQuadraticIdeal: The entered list of elements is empty, so quadratic ideal." "$raw"
tail -n 1 "$raw" | jq --sort-keys . >"$actual"
jq --sort-keys . \
    "$repo_root/tests/fixtures/pathalg_koszul.schema-v1.json" >"$expected"
diff --unified "$expected" "$actual"
"$repo_root/scripts/gap-reference" \
    "$repo_root/tests/gap-reference/pathalg_koszul_nonquadratic.g" >"$actual"
grep -Fq "BasisOfIperp(): ideal is not quadratic" "$actual"
grep -Fq "NONQUADRATIC_FAIL=true" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_koszul_empty_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly reached its documented empty quadratic-perpendicular branch" >&2
    exit 1
fi
grep -Fq "IsQuadraticIdeal: The entered list of elements is empty, so quadratic ideal." "$actual"
grep -Fq "List Element: <list>[1] must have an assigned value" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/pathalg_koszul_no_paths_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly retained the empty degree-two matrix shape" >&2
    exit 1
fi
grep -Fq "no 1st choice method found for \`NullspaceMat'" "$actual"
compare_fixture \
    "$repo_root/tests/gap-reference/present2.g" \
    "$repo_root/tests/fixtures/present2.schema-v1.json"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/present2_valid_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly constructed a legacy projective presentation" >&2
    exit 1
fi
grep -Fq "no 1st choice method found for \`CoefficientsAndMagmaElements'" "$actual"
grep -Fq "lib/pathalg.gi:588" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/present2_zero_row_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted an all-zero projective-presentation row" >&2
    exit 1
fi
grep -Fq "The 1st argument is 'fail'" "$actual"
grep -Fq "lib/present2.gi:41" "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/present2_ragged_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly accepted a ragged projective-presentation matrix" >&2
    exit 1
fi
grep -Fq "Matrix must have rows of equal length." "$actual"
if "$repo_root/scripts/gap-reference" --alwaystrace \
    "$repo_root/tests/gap-reference/present2_empty_error.g" \
    >"$actual" 2>&1
then
    echo "QPA unexpectedly dispatched ProjectivePresentation on an empty list" >&2
    exit 1
fi
grep -Fq "no 1st choice method found for \`ProjectivePresentation'" "$actual"
compare_catching_fixture \
    "$repo_root/tests/gap-reference/groebner_objects.g" \
    "$repo_root/tests/fixtures/groebner_objects.schema-v1.json"
compare_catching_fixture \
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
