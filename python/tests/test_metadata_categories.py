from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

import qpa
from qpa import _qpa


FIXTURES = Path(__file__).resolve().parents[2] / "tests" / "fixtures"


def test_build_exposes_port_git_commit() -> None:
    assert qpa.PORT_GIT_COMMIT == _qpa.PORT_GIT_COMMIT
    assert qpa.PORT_GIT_COMMIT == "unknown" or re.fullmatch(
        r"(?:[0-9a-f]{40}|[0-9a-f]{64})(?:-dirty)?",
        qpa.PORT_GIT_COMMIT,
    )


def _gap_list(values: list[object]) -> str:
    return "[ " + ", ".join(str(value) for value in values) + " ]"


def _gap_bool(value: bool) -> str:
    return str(value).lower()


def test_information_classes_and_documentation_plan_follow_qpa() -> None:
    expected = json.loads(
        (FIXTURES / "final_declared_surface.schema-v1.json").read_text()
    )["documentation"]

    assert qpa.info_groebner_basis_name() == "InfoGroebnerBasis"
    assert qpa.info_path_algebra_module_name() == "InfoPathAlgebraModule"
    assert qpa.info_path_algebra_module_default_level() == 0
    assert (
        qpa.info_projective_resolution_fp_path_algebra_module_name()
        == "InfoProjectiveResolutionFpPathAlgebraModule"
    )
    assert qpa.info_projective_resolution_fp_path_algebra_module_default_level() == 0
    assert qpa.info_quiver_name() == "qpa::quiver"

    plan = qpa.make_qpa_documentation("/opt/QPA")
    assert isinstance(plan, qpa.QpaDocumentationPlan)
    assert {
        "function_bound": True,
        "directory_count": 1,
        "main_exists": True,
        "main": plan.main_file,
        "files": plan.additional_files,
        "book_name": plan.book_name,
        "renderer": plan.renderer,
        "copy_html_style_files": plan.copy_html_style_files,
    } == expected
    # The Rust plan joins paths natively, so Windows renders the directory
    # with a backslash; compare as paths rather than as strings.
    assert Path(plan.directory) == Path("/opt/QPA") / "doc"
    shown = repr(plan)
    assert shown.startswith("<QPA documentation plan: ") and shown.endswith(">")
    assert Path(shown[len("<QPA documentation plan: ") : -1]) == Path(
        "/opt/QPA/doc/qpadocumentation.xml"
    )


def test_generic_cat_delegates_to_qpa_callback_category() -> None:
    properties: dict[str, object] = {
        "name": "toy abelian category",
        "zeroObj": 0,
        "isZeroObj": lambda value: value == 0,
        "zeroMap": lambda source, range: (source, range, 0),
        "isZeroMapping": lambda value: value[2] == 0,
        "composeMaps": lambda left, right: (right[0], left[1], left[2] * right[2]),
        "ker": lambda value: value[0] if value[2] == 0 else 0,
        "im": lambda value: 0 if value[2] == 0 else value[1],
        "isExact": lambda left, right: left[2] * right[2] == 0,
        "objStr": lambda value: f"({value})",
    }
    category = qpa.cat(["toy", 1], properties)

    assert isinstance(category, qpa.CallbackCategory)
    assert category.is_cat and category.is_default_rep
    assert category.identity == ["toy", 1]
    assert category.name == "toy abelian category"
    assert repr(category) == "<cat: toy abelian category>"
    assert category.property("zeroObj") == 0
    assert category.zero_object() == 0
    assert category.is_zero_object(0) and not category.is_zero_object(1)
    zero = category.zero_map(2, 3)
    assert zero == (2, 3, 0) and category.is_zero_mapping(zero)
    assert category.compose_maps((2, 3, 5), (1, 2, 7)) == (1, 3, 35)
    assert category.kernel(zero) == 2 and category.image(zero) == 0
    assert category.is_exact((2, 3, 0), (1, 2, 7))
    assert category.object_string(4) == "(4)"
    assert category == qpa.cat(["toy", 1], properties)
    assert category != qpa.cat(["toy", 2], properties)

    with pytest.raises(qpa.InvalidInputError, match="no zeroObj entry"):
        qpa.cat("broken", {"name": "broken"})


def test_path_module_element_spaces_and_bases_follow_qpa_oracle() -> None:
    field = qpa.Field.rational()
    algebra = qpa.PathAlgebra(
        field,
        qpa.Quiver(["v1", "v2"], [("v1", "v2", "a")]),
    )
    module = algebra.module_from_integer_matrices([2, 1], [("a", [[1], [2]])])
    standard = module.basis_elements()
    entered = [module.zero_element(), standard[1], standard[0], standard[2]]
    space = module.element_space(entered)

    new_basis = space.new_basis(entered)
    unchecked = space.basis_nc(entered)
    checked = space.basis_from_vectors(entered)
    assert checked is not None
    domain_basis = space.basis_of_domain()
    default_basis = space.basis()
    canonical = module.canonical_basis()
    line = module.element_space([standard[0]])
    value = (
        standard[0].scale_left(field.element(2))
        - standard[1].scale_left(field.element(3))
        + standard[2].scale_left(field.element(5))
    )
    coefficients = checked.coefficients(value)
    assert coefficients is not None

    prime = qpa.Field.prime(2)
    prime_algebra = qpa.PathAlgebra(
        prime,
        qpa.Quiver(["p1", "p2"], [("p1", "p2", "b")]),
    )
    prime_module = prime_algebra.module_from_integer_matrices(
        [2, 1], [("b", [[1], [1]])]
    )
    prime_standard = prime_module.basis_elements()
    prime_basis = prime_module.element_space(
        [prime_standard[2], prime_standard[0], prime_standard[1]]
    ).basis_from_vectors([prime_standard[2], prime_standard[0], prime_standard[1]])
    assert prime_basis is not None
    prime_coefficients = prime_basis.coefficients(prime_standard[0] + prime_standard[2])
    assert prime_coefficients is not None

    actual = "\n".join(
        [
            f"NEW_VECTORS={_gap_list(new_basis.basis_vectors)}",
            f"UNCHECKED_VECTORS={_gap_list(unchecked.basis_vectors)}",
            f"CHECKED_VECTORS={_gap_list(checked.basis_vectors)}",
            f"CHECKED_COEFFICIENTS={_gap_list(coefficients)}",
            f"DOMAIN_VECTORS={_gap_list(domain_basis.basis_vectors)}",
            f"DEFAULT_VECTORS={_gap_list(default_basis.basis_vectors)}",
            f"CANONICAL_VECTORS={_gap_list(canonical.basis_vectors)}",
            f"CANONICAL_PROPERTY={_gap_bool(canonical.is_canonical_basis)}",
            f"MEMBERSHIP={_gap_list([_gap_bool(standard[0] in line), _gap_bool(standard[1] in line), _gap_bool(module.zero_element() in line)])}",
            f"MODULE_IS_ZERO={_gap_list([_gap_bool(module.is_zero), _gap_bool(algebra.zero_module().is_zero)])}",
            "GF2_COEFFICIENTS="
            + _gap_list(
                [
                    "0*Z(2)" if value.representative == 0 else "Z(2)^0"
                    for value in prime_coefficients
                ]
            ),
        ]
    ) + "\n"
    assert actual == (FIXTURES / "module_basis.txt").read_text()

    assert module.is_algebra_module and module.is_path_algebra_mat_module
    assert standard[0].is_path_module_element
    assert space.is_path_module_element_collection
    assert checked.is_basis_of_path_module_elem_vector_space
    assert checked.underlying_space.module.is_identical(module)
    assert len(checked) == 3 and not checked.is_empty
    assert module.identity().is_algebra_module_homomorphism
    assert module.identity().is_path_algebra_mat_module_homomorphism_rep
    assert space.basis_from_vectors([standard[0], standard[0]]) is None
    with pytest.raises(qpa.InvalidInputError, match="delegate"):
        new_basis.coefficients(value)
