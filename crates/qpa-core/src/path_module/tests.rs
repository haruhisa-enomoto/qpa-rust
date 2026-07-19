use std::sync::Arc;

use numerica::domains::{
    Field, Ring,
    finite_field::{FiniteFieldCore, Z2, Zp},
    integer::Integer,
    rational::{Q, Rational},
};
use serde_json::{Value, json};

use crate::{
    ArrowIndex, ArrowSpec, DenseMatrix, PathAlgebra, QpaCoefficientField, Quiver, QuiverAlgebra,
    QuiverAlgebraElement, QuiverAlgebraHomomorphism, QuotientPathAlgebra, StructureConstantsError,
    VertexIndex,
};

use super::{
    BasedRightAlgebraModule, CommonDirectSummandResult, GenericPathAlgebraModuleHomomorphism,
    IndexedPathAlgebraModuleArrowAction, ModuleEndomorphismAlgebra, ModuleEndomorphismAlgebraData,
    MoreMinimalVersionResult, NamedArrowMatrix, NamedPathAlgebraModuleArrowAction,
    PathAlgebraModule, PathAlgebraModuleArrowAction, PathAlgebraModuleArrowActions,
    PathAlgebraModuleError, PathAlgebraModuleHomomorphism, PathAlgebraModuleMatrix,
    cokernel_projection, common_direct_summand, direct_sum_inclusions, direct_sum_of_qpa_modules,
    direct_sum_projections, end_over_algebra, generic_end_over_path_ring,
    generic_hom_over_path_ring, hom_over_algebra, hom_over_algebra_with_basis_function,
    homomorphism_from_images, identity_module_mapping, image_projection_inclusion,
    intersection_of_submodules, intersection_of_submodules_list, inverse_of_isomorphism,
    is_direct_summand, is_in_additive_closure, is_in_full_matrix_ring, is_left_minimal,
    is_right_minimal, is_split_epimorphism, is_split_monomorphism, isomorphic_modules,
    isomorphism_of_modules, kernel_inclusion, left_inverse_of_homomorphism, left_minimal_version,
    lift_idempotents, maximal_common_direct_summand, minimal_generating_set_of_module,
    more_left_minimal_version, more_right_minimal_version, morphism_on_cokernel, morphism_on_image,
    morphism_on_kernel, multiply_lists_of_maps, natural_endomorphism_quotient_by_ideal,
    natural_homomorphism_by_subalgebra_module, number_of_nonisomorphic_direct_summands, pull_back,
    push_out, radical_of_module_inclusion, reject_of_module, restriction_via_algebra_homomorphism,
    restriction_via_algebra_homomorphism_map, right_algebra_module_to_path_algebra_mat_module,
    right_fac_approximation, right_fac_approximation_via_all_basis_images,
    right_inverse_of_homomorphism, right_minimal_version, right_module_hom_over_algebra,
    right_module_over_path_algebra, right_module_over_path_algebra_infer_dimension,
    right_module_over_path_algebra_infer_dimension_nc, socle_of_module, socle_of_module_inclusion,
    subalgebra_module, sum_of_submodules, sum_of_submodules_list, top_of_module,
    top_of_module_homomorphism, top_of_module_projection, trace_of_module, zero_module,
    zero_module_mapping,
};

fn vertex(index: usize) -> VertexIndex {
    VertexIndex::new(index).unwrap()
}

fn q(value: i64) -> Rational {
    Q.nth(value.into())
}

fn qpa_matrix<F: Field>(
    matrix: &DenseMatrix<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> String {
    super::qpa_list(
        (0..matrix.nrows()).map(|row| super::qpa_list(matrix.row(row).iter().map(coefficient))),
    )
}

fn qpa_module_data<F: Field>(
    name: &str,
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> String {
    let matrices =
        super::qpa_list(module.matrices_of_path_algebra_module().iter().map(
            |matrix| match matrix {
                PathAlgebraModuleMatrix::ZeroScalar => "0".to_owned(),
                PathAlgebraModuleMatrix::Matrix(matrix) => qpa_matrix(matrix, coefficient),
            },
        ));
    format!(
        "{name}={}|{matrices}\n",
        super::qpa_usize_list(module.dimension_vector())
    )
}

fn qpa_element_list<F: Field>(
    elements: &[super::PathAlgebraModuleElement<F>],
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> String {
    super::qpa_list(
        elements
            .iter()
            .map(|element| element.qpa_print_string_with(coefficient)),
    )
}

fn qpa_module_matrix_list<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> String {
    super::qpa_list(
        module
            .matrices_of_path_algebra_module()
            .iter()
            .map(|matrix| match matrix {
                PathAlgebraModuleMatrix::ZeroScalar => "0".to_owned(),
                PathAlgebraModuleMatrix::Matrix(matrix) => qpa_matrix(matrix, coefficient),
            }),
    )
}

fn qpa_nonzero_positions<F: Field>(matrix: &DenseMatrix<F>) -> String {
    let mut positions = Vec::new();
    let mut position = 1_usize;
    for row in 0..matrix.nrows() {
        for entry in matrix.row(row) {
            if !matrix.field().is_zero(entry) {
                positions.push(position);
            }
            position += 1;
        }
    }
    super::qpa_usize_list(&positions)
}

fn qpa_generic_hom_supports<F: Field>(
    mappings: &[GenericPathAlgebraModuleHomomorphism<F>],
) -> String {
    super::qpa_list(
        mappings
            .iter()
            .map(|mapping| qpa_nonzero_positions(mapping.matrix())),
    )
}

#[test]
fn pinned_qpa_module_and_homomorphism_presentation_matches() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let a2 = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let source = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let range = right_module_over_path_algebra(a2.clone(), vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let zero = zero_module(&a2).unwrap();
    let map = zero_module_mapping(source.clone(), range.clone()).unwrap();
    let wide = right_module_over_path_algebra(
        a2.clone(),
        vec![2, 2],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![q(1), q(-2)], vec![q(0), q(3)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let loop_parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual = QuiverAlgebra::Quotient(
        QuotientPathAlgebra::by_relations_high_level(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
    );
    let dual_simple = right_module_over_path_algebra(
        dual,
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![Q.zero()]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let prime_field = Zp::new(5);
    let prime_a1 = QuiverAlgebra::Path(
        PathAlgebra::new(
            prime_field.clone(),
            Quiver::from_indexed(1, Vec::new()).unwrap(),
        )
        .unwrap(),
    );
    let prime_module = right_module_over_path_algebra(prime_a1, vec![1], &[])
        .unwrap()
        .unwrap();

    let actual = format!(
        "MODULE_VIEW={source}\nZERO_VIEW={zero}\nMODULE_PRINT={}\nQUOTIENT_MODULE_PRINT={}\nPRIME_MODULE_PRINT={}\nMAP_VIEW={map}\nMAP_PRINT={}\nMODULE_DISPLAY_BEGIN\n{}MODULE_DISPLAY_END\nMAP_DISPLAY_BEGIN\n{}MAP_DISPLAY_END\nWIDE_MODULE_DISPLAY_BEGIN\n{}WIDE_MODULE_DISPLAY_END\n",
        source.qpa_print_string(),
        dual_simple.qpa_print_string(),
        prime_module.qpa_print_string(),
        map.qpa_print_string(),
        source.qpa_display_string_with(ToString::to_string).unwrap(),
        map.qpa_display_string_with(ToString::to_string),
        wide.qpa_display_string_with(ToString::to_string).unwrap(),
    );
    assert_eq!(
        actual,
        include_str!("../../../../tests/fixtures/module_presentation.txt")
    );
    assert!(matches!(
        zero.qpa_display_string_with(ToString::to_string),
        Err(PathAlgebraModuleError::DisplayScalarZeroArrowMatrix)
    ));
}

#[test]
#[allow(clippy::too_many_lines)]
fn pinned_qpa_module_element_surface_matches() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let parent = PathAlgebra::new(Q, quiver).unwrap();
    let generators = parent.generators();
    let algebra = QuiverAlgebra::Path(parent);
    let module = right_module_over_path_algebra(
        algebra.clone(),
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![q(1)], vec![q(2)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let basis = module.basis_elements();
    let v = basis[0]
        .scale_left(&q(2))
        .subtract(&basis[1].scale_left(&q(3)))
        .unwrap()
        .add(&basis[2].scale_left(&q(5)))
        .unwrap();
    let w = basis[0]
        .add(&basis[1])
        .unwrap()
        .subtract(&basis[2])
        .unwrap();
    let rebuilt = module
        .element_from_external_representation(vec![vec![q(2), q(-3)], vec![q(5)]])
        .unwrap();
    let loop_parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuotientPathAlgebra::by_relations_high_level(loop_parent, &[x.pow(2).unwrap()]).unwrap();
    let dual_generators = dual.generators().unwrap();
    let dual_module = right_module_over_path_algebra(
        QuiverAlgebra::Quotient(dual),
        vec![2],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![q(0), q(1)], vec![q(0), q(0)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let quotient_action = dual_module.basis_elements()[0]
        .act(&QuiverAlgebraElement::Quotient(dual_generators[1].clone()))
        .unwrap();
    let qpa_basis = format!(
        "[ {} ]",
        basis
            .iter()
            .map(|element| element.qpa_print_string_with(ToString::to_string))
            .collect::<Vec<_>>()
            .join(", ")
    );
    let qpa = |element: &super::PathAlgebraModuleElement<Q>| {
        element.qpa_print_string_with(ToString::to_string)
    };
    let vertex_one = QuiverAlgebraElement::Path(generators[0].clone());
    let vertex_two = QuiverAlgebraElement::Path(generators[1].clone());
    let arrow = QuiverAlgebraElement::Path(generators[2].clone());
    let action_sum = QuiverAlgebraElement::Path(
        generators[0]
            .add(&generators[2].scale_left(&q(2)).unwrap())
            .unwrap(),
    );

    let zero_vertex_module = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let pathological = zero_vertex_module
        .element_from_external_representation(vec![vec![q(0)], Vec::new()])
        .unwrap();
    let zero_vertex_raw = zero_vertex_module.basis_elements()[0].clone();

    let prime_parent = PathAlgebra::new(
        Z2,
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "b")]).unwrap(),
    )
    .unwrap();
    let prime_generators = prime_parent.generators();
    let prime_module = right_module_over_path_algebra(
        QuiverAlgebra::Path(prime_parent),
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "b",
            DenseMatrix::from_rows(vec![vec![Z2.one()], vec![Z2.one()]], Z2).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let prime_basis = prime_module.basis_elements();
    let prime_v = prime_basis[0].add(&prime_basis[2]).unwrap();
    let prime_action = prime_v
        .act(&QuiverAlgebraElement::Path(prime_generators[2].clone()))
        .unwrap();
    let gap_z2 = |entry: &u8| {
        if Z2.is_zero(entry) {
            "0*Z(2)".to_owned()
        } else {
            "Z(2)^0".to_owned()
        }
    };
    let rational_matrix = DenseMatrix::from_rows(vec![vec![q(1), q(-2)]], Q).unwrap();
    let field5 = Zp::new(5);
    let field7 = Zp::new(7);
    let matrix5 = DenseMatrix::from_rows(vec![vec![field5.one()]], field5.clone()).unwrap();
    let full_matrix_qq = is_in_full_matrix_ring(&rational_matrix, &Q);
    let full_matrix_foreign = is_in_full_matrix_ring(&matrix5, &field7);

    let actual = format!(
        "FULL_MATRIX_QQ={full_matrix_qq}\nFULL_MATRIX_FOREIGN={full_matrix_foreign}\nBASIS={qpa_basis}\nV={}\nEXTREP={}\nSUM={}\nDIFFERENCE={}\nLEFT_SCALE={}\nRIGHT_SCALE={}\nNEGATIVE={}\nZERO={}\nEQUAL_SELF={}\nEQUAL_OTHER={}\nACTION_VERTEX_1={}\nACTION_VERTEX_2={}\nACTION_ARROW={}\nACTION_SUM={}\nREBUILT={}\nQUOTIENT_ACTION={}\nZERO_DIMENSION_EXTERNAL={}\nZERO_ROW_SUM_LEFT={}\nZERO_ROW_SUM_RIGHT={}\nZERO_ROW_DIFFERENCE_LEFT={}\nZERO_ROW_DIFFERENCE_RIGHT={}\nGF2_ACTION={}\n",
        qpa(&v),
        qpa(&v),
        qpa(&v.add(&w).unwrap()),
        qpa(&v.subtract(&w).unwrap()),
        qpa(&w.scale_left(&q(2))),
        qpa(&w.scale_right(&q(2))),
        qpa(&v.additive_inverse()),
        qpa(&module.zero_element()),
        v == v.clone(),
        v == w,
        qpa(&v.act(&vertex_one).unwrap()),
        qpa(&v.act(&vertex_two).unwrap()),
        qpa(&v.act(&arrow).unwrap()),
        qpa(&v.act(&action_sum).unwrap()),
        qpa(&rebuilt),
        qpa(&quotient_action),
        qpa(&pathological),
        qpa(&pathological.add(&zero_vertex_raw).unwrap()),
        qpa(&zero_vertex_raw.add(&pathological).unwrap()),
        qpa(&pathological.subtract(&zero_vertex_raw).unwrap()),
        qpa(&zero_vertex_raw.subtract(&pathological).unwrap()),
        prime_action.qpa_print_string_with(gap_z2),
    );
    assert_eq!(
        actual,
        include_str!("../../../../tests/fixtures/module_elements.txt")
    );
    assert_eq!(v, rebuilt);

    assert!(matches!(
        module.element_from_external_representation(vec![vec![q(1), q(0)]]),
        Err(PathAlgebraModuleError::ElementVertexCount { .. })
    ));
    assert!(matches!(
        module.element_from_external_representation(vec![vec![q(1)], vec![q(0)]]),
        Err(PathAlgebraModuleError::ElementCoordinateCount { vertex: 1, .. })
    ));
    assert_eq!(
        qpa(&pathological.add(&zero_vertex_raw).unwrap()),
        "[ [ 1 ], [ 0 ] ]"
    );
    assert_eq!(
        qpa(&zero_vertex_raw.add(&pathological).unwrap()),
        "[ [ 1 ], [ 0 ] ]"
    );
    assert_eq!(
        qpa(&pathological.subtract(&zero_vertex_raw).unwrap()),
        "[ [ -1 ], [ 0 ] ]"
    );
    assert_eq!(
        qpa(&zero_vertex_raw.subtract(&pathological).unwrap()),
        "[ [ 1 ], [ 0 ] ]"
    );
    let reconstructed = right_module_over_path_algebra(algebra, vec![2, 1], &[])
        .unwrap()
        .unwrap();
    assert_eq!(
        v.add(&reconstructed.basis_elements()[0]).unwrap_err(),
        PathAlgebraModuleError::DifferentElementFamilies
    );
    let foreign_parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "foreign")]).unwrap(),
    )
    .unwrap();
    let foreign_arrow = QuiverAlgebraElement::Path(foreign_parent.generators()[2].clone());
    assert_eq!(
        v.act(&foreign_arrow).unwrap_err(),
        PathAlgebraModuleError::ForeignActingElement
    );

    assert!(is_in_full_matrix_ring(&matrix5, &field5));
    assert!(!full_matrix_foreign);
}

#[test]
#[allow(clippy::too_many_lines)]
fn pinned_qpa_module_element_basis_surface_matches() {
    let parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
    )
    .unwrap();
    let module = right_module_over_path_algebra(
        QuiverAlgebra::Path(parent),
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![q(1)], vec![q(2)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let standard = module.basis_elements();
    let entered = vec![
        module.zero_element(),
        standard[1].clone(),
        standard[0].clone(),
        standard[2].clone(),
    ];
    let space =
        super::PathAlgebraModuleElementSpace::by_generators(module.clone(), entered.clone())
            .unwrap();
    let new_basis = space.new_basis(entered.clone()).unwrap();
    let unchecked = space.basis_nc(entered.clone()).unwrap();
    let checked = space.basis_from_vectors(entered.clone()).unwrap().unwrap();
    let domain_basis = space.basis_of_domain().unwrap();
    let default_basis = space.basis().unwrap();
    let canonical = super::PathAlgebraModuleElementBasis::canonical_for_module(&module).unwrap();
    let line = super::PathAlgebraModuleElementSpace::by_generators(
        module.clone(),
        vec![standard[0].clone()],
    )
    .unwrap();
    let membership = [
        line.contains(&standard[0]).unwrap(),
        line.contains(&standard[1]).unwrap(),
        line.contains(&module.zero_element()).unwrap(),
    ];
    let zero = zero_module(module.right_acting_algebra()).unwrap();
    let v = standard[0]
        .scale_left(&q(2))
        .subtract(&standard[1].scale_left(&q(3)))
        .unwrap()
        .add(&standard[2].scale_left(&q(5)))
        .unwrap();
    let coefficients = checked.coefficients(&v).unwrap().unwrap();
    let qpa_elements = |elements: &[super::PathAlgebraModuleElement<Q>]| {
        super::qpa_list(
            elements
                .iter()
                .map(|element| element.qpa_print_string_with(ToString::to_string)),
        )
    };

    let prime_parent = PathAlgebra::new(
        Z2,
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "b")]).unwrap(),
    )
    .unwrap();
    let prime_module = right_module_over_path_algebra(
        QuiverAlgebra::Path(prime_parent),
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "b",
            DenseMatrix::from_rows(vec![vec![Z2.one()], vec![Z2.one()]], Z2).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let prime_standard = prime_module.basis_elements();
    let prime_entered = vec![
        prime_standard[2].clone(),
        prime_standard[0].clone(),
        prime_standard[1].clone(),
    ];
    let prime_space =
        super::PathAlgebraModuleElementSpace::by_generators(prime_module, prime_entered.clone())
            .unwrap();
    let prime_basis = prime_space
        .basis_from_vectors(prime_entered)
        .unwrap()
        .unwrap();
    let prime_coefficients = prime_basis
        .coefficients(&prime_standard[0].add(&prime_standard[2]).unwrap())
        .unwrap()
        .unwrap();
    let gap_z2 = |entry: &u8| {
        if Z2.is_zero(entry) {
            "0*Z(2)".to_owned()
        } else {
            "Z(2)^0".to_owned()
        }
    };

    let actual = format!(
        "NEW_VECTORS={}\nUNCHECKED_VECTORS={}\nCHECKED_VECTORS={}\nCHECKED_COEFFICIENTS={}\nDOMAIN_VECTORS={}\nDEFAULT_VECTORS={}\nCANONICAL_VECTORS={}\nCANONICAL_PROPERTY={}\nMEMBERSHIP={}\nMODULE_IS_ZERO={}\nGF2_COEFFICIENTS={}\n",
        qpa_elements(new_basis.basis_vectors()),
        qpa_elements(unchecked.basis_vectors()),
        qpa_elements(checked.basis_vectors()),
        super::qpa_list(coefficients.iter().map(ToString::to_string)),
        qpa_elements(domain_basis.basis_vectors()),
        qpa_elements(default_basis.basis_vectors()),
        qpa_elements(canonical.basis_vectors()),
        canonical.is_canonical_basis(),
        super::qpa_list(membership.iter().map(ToString::to_string)),
        super::qpa_list(
            [module.is_zero(), zero.is_zero()]
                .iter()
                .map(ToString::to_string)
        ),
        super::qpa_list(prime_coefficients.iter().map(gap_z2)),
    );
    assert_eq!(
        actual,
        include_str!("../../../../tests/fixtures/module_basis.txt")
    );
    assert!(!canonical.is_empty());
    assert_eq!(
        new_basis.coefficients(&v).unwrap_err(),
        PathAlgebraModuleError::MissingElementBasisDelegate
    );
    assert!(
        space
            .basis_from_vectors(vec![standard[0].clone(), standard[0].clone()])
            .unwrap()
            .is_none()
    );
    assert!(
        line.basis_from_vectors(vec![standard[1].clone()])
            .unwrap()
            .is_none()
    );
    assert_eq!(
        super::PathAlgebraModuleElementBasis::canonical_for_module(&zero).unwrap_err(),
        PathAlgebraModuleError::ZeroModuleCanonicalBasisFailure
    );
}

#[test]
fn number_of_nonisomorphic_summands_reports_rational_central_split_tier() {
    let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let module = right_module_over_path_algebra(algebra, vec![2], &[])
        .unwrap()
        .unwrap();
    assert!(matches!(
        number_of_nonisomorphic_direct_summands(&module),
        Err(PathAlgebraModuleError::StructureConstants(
            StructureConstantsError::FiniteFieldRequired
        ))
    ));
}

#[test]
fn linear_mapping_algebra_radical_uses_qpa_transpose_conversion_semantics() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/algebra_mapping_radical.schema-v1.json"
    ))
    .unwrap();
    let identity = DenseMatrix::identity(2, Q);
    let nilpotent = DenseMatrix::from_rows(vec![vec![q(0), q(1)], vec![q(0), q(0)]], Q).unwrap();
    let flattened_basis = DenseMatrix::from_rows(
        vec![identity.as_flat().to_vec(), nilpotent.as_flat().to_vec()],
        Q,
    )
    .unwrap();
    let algebra = ModuleEndomorphismAlgebra(Arc::new(ModuleEndomorphismAlgebraData {
        field: Q,
        basis_matrices: vec![identity.clone(), nilpotent.clone()],
        flattened_basis,
        one: identity,
    }));

    let radical = algebra.radical_basis_matrices().unwrap();
    assert_eq!(radical, [nilpotent]);
    let actual = json!({
        "dimension": radical.len(),
        "basis_matrices": radical.iter().map(|matrix| {
            (0..matrix.nrows()).map(|row| {
                matrix.row(row).iter().map(|entry| {
                    entry.to_string().parse::<i64>().unwrap()
                }).collect::<Vec<_>>()
            }).collect::<Vec<_>>()
        }).collect::<Vec<_>>(),
        "all_square_zero": radical.iter().all(|matrix| {
            matrix.multiply(matrix).unwrap().is_zero()
        }),
    });
    assert_eq!(actual["dimension"], fixture["dimension"]);
    assert_eq!(actual["basis_matrices"], fixture["basis_matrices"]);
    assert_eq!(actual["all_square_zero"], fixture["all_square_zero"]);
    assert_eq!(fixture["had_basis_before"], false);
}

#[test]
fn linear_mapping_algebra_radical_supports_numerica_prime_fields() {
    let e11 = DenseMatrix::from_rows(
        vec![vec![Z2.one(), Z2.zero()], vec![Z2.zero(), Z2.zero()]],
        Z2,
    )
    .unwrap();
    let e12 = DenseMatrix::from_rows(
        vec![vec![Z2.zero(), Z2.one()], vec![Z2.zero(), Z2.zero()]],
        Z2,
    )
    .unwrap();
    let e22 = DenseMatrix::from_rows(
        vec![vec![Z2.zero(), Z2.zero()], vec![Z2.zero(), Z2.one()]],
        Z2,
    )
    .unwrap();
    let basis_matrices = vec![e11.clone(), e12.clone(), e22.clone()];
    let flattened_basis = DenseMatrix::from_rows(
        basis_matrices
            .iter()
            .map(|matrix| matrix.as_flat().to_vec())
            .collect(),
        Z2,
    )
    .unwrap();
    let algebra = ModuleEndomorphismAlgebra(Arc::new(ModuleEndomorphismAlgebraData {
        field: Z2,
        basis_matrices,
        flattened_basis,
        one: e11.add(&e22).unwrap(),
    }));

    assert_eq!(algebra.radical_basis_matrices().unwrap(), [e12]);
}

fn matrix_json<F: Field>(
    matrix: &DenseMatrix<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!({
        "shape": [matrix.nrows(), matrix.ncols()],
        "rows": (0..matrix.nrows())
            .map(|row| matrix.row(row).iter().map(coefficient).collect::<Vec<_>>())
            .collect::<Vec<_>>(),
    })
}

fn module_map_json<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!({
        "source_dimension": homomorphism.source().dimension_vector(),
        "range_dimension": homomorphism.range().dimension_vector(),
        "matrices": homomorphism.vertex_matrices().iter()
            .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
    })
}

fn restriction_module_json<F: Field>(
    module: &super::PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    json!({
        "dimension": module.dimension_vector(),
        "matrices": module.arrow_matrices().iter().map(|matrix| {
            (0..matrix.nrows()).map(|row| {
                matrix.row(row).iter().map(coefficient).collect::<Vec<_>>()
            }).collect::<Vec<_>>()
        }).collect::<Vec<_>>(),
    })
}

fn restriction_map_json<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    json!({
        "source": restriction_module_json(homomorphism.source(), coefficient),
        "range": restriction_module_json(homomorphism.range(), coefficient),
        "matrices": homomorphism.vertex_matrices().iter().map(|matrix| {
            (0..matrix.nrows()).map(|row| {
                matrix.row(row).iter().map(coefficient).collect::<Vec<_>>()
            }).collect::<Vec<_>>()
        }).collect::<Vec<_>>(),
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_restriction_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual = QuiverAlgebra::Quotient(
        QuotientPathAlgebra::by_relations_high_level(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
    );
    let dual_generators = dual.generators_with_one().unwrap();
    let zero_arrow_map = QuiverAlgebraHomomorphism::by_generators_with_one_images(
        dual.clone(),
        dual.clone(),
        vec![dual_generators[0].clone(), dual.zero()],
    )
    .unwrap();
    let radical_matrix = DenseMatrix::from_rows(
        vec![
            vec![field.zero(), field.one()],
            vec![field.zero(), field.zero()],
        ],
        field.clone(),
    )
    .unwrap();
    let regular = right_module_over_path_algebra(
        dual,
        vec![2],
        &[NamedArrowMatrix::new("x", radical_matrix.clone())],
    )
    .unwrap()
    .unwrap();
    let radical_map =
        right_module_hom_over_algebra(regular.clone(), regular.clone(), vec![radical_matrix])
            .unwrap();
    let restricted_regular =
        restriction_via_algebra_homomorphism(&zero_arrow_map, &regular).unwrap();
    let restricted_radical_map =
        restriction_via_algebra_homomorphism_map(&zero_arrow_map, &radical_map).unwrap();

    let a2_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());
    let field_algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
        )
        .unwrap(),
    );
    let field_inclusion = QuiverAlgebraHomomorphism::by_generators_with_one_images(
        field_algebra,
        a2.clone(),
        vec![a2.one().unwrap()],
    )
    .unwrap();
    let projective = right_module_over_path_algebra(
        a2,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![field.one()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let identity_map = identity_module_mapping(projective.clone()).unwrap();
    let restricted_projective =
        restriction_via_algebra_homomorphism(&field_inclusion, &projective).unwrap();
    let restricted_identity =
        restriction_via_algebra_homomorphism_map(&field_inclusion, &identity_map).unwrap();

    json!({
        "zero_arrow_module": restriction_module_json(&restricted_regular, coefficient),
        "zero_arrow_map": restriction_map_json(&restricted_radical_map, coefficient),
        "field_module": restriction_module_json(&restricted_projective, coefficient),
        "field_identity_map": restriction_map_json(&restricted_identity, coefficient),
    })
}

fn module_map_list_json<F: Field>(
    homomorphisms: &[&PathAlgebraModuleHomomorphism<F>],
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!(
        homomorphisms
            .iter()
            .map(|homomorphism| module_map_json(homomorphism, coefficient))
            .collect::<Vec<_>>()
    )
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_split_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let zero = field.zero();
    let one = field.one();
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let layer = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
        .unwrap()
        .unwrap();
    let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let simple_two = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let zero_module = right_module_over_path_algebra(algebra.clone(), vec![0, 0], &[])
        .unwrap()
        .unwrap();
    let projective_one = right_module_over_path_algebra(
        algebra,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let inclusion = right_module_hom_over_algebra(
        simple_one.clone(),
        layer.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
            DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
        ],
    )
    .unwrap();
    let projection = right_module_hom_over_algebra(
        layer.clone(),
        simple_one.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
            DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
        ],
    )
    .unwrap();
    let right_inverse = right_inverse_of_homomorphism(&inclusion).unwrap().unwrap();
    let left_inverse = left_inverse_of_homomorphism(&projection).unwrap().unwrap();
    let radical_inclusion = right_module_hom_over_algebra(
        simple_two,
        projective_one.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
            DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
        ],
    )
    .unwrap();
    let top_projection = right_module_hom_over_algebra(
        projective_one,
        simple_one.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![one]], field.clone()).unwrap(),
            DenseMatrix::from_rows(vec![vec![zero]], field.clone()).unwrap(),
        ],
    )
    .unwrap();
    let zero_inclusion = zero_module_mapping(zero_module.clone(), layer.clone()).unwrap();
    let zero_projection = zero_module_mapping(layer.clone(), zero_module).unwrap();
    let identity = identity_module_mapping(layer.clone()).unwrap();
    let identity_right = right_inverse_of_homomorphism(&identity).unwrap().unwrap();
    let identity_left = left_inverse_of_homomorphism(&identity).unwrap().unwrap();
    let property_inclusion = right_module_hom_over_algebra(
        simple_one.clone(),
        layer.clone(),
        vec![
            DenseMatrix::identity(1, field.clone()),
            DenseMatrix::zero(1, 1, field.clone()),
        ],
    )
    .unwrap();
    let property_projection = right_module_hom_over_algebra(
        layer,
        simple_one.clone(),
        vec![
            DenseMatrix::identity(1, field.clone()),
            DenseMatrix::zero(1, 1, field),
        ],
    )
    .unwrap();
    let property_values = [
        is_split_monomorphism(&property_inclusion).unwrap(),
        is_split_epimorphism(&property_projection).unwrap(),
    ];
    let property_right = right_inverse_of_homomorphism(&property_inclusion)
        .unwrap()
        .unwrap();
    let property_left = left_inverse_of_homomorphism(&property_projection)
        .unwrap()
        .unwrap();

    json!({
        "right_inverse": module_map_json(&right_inverse, coefficient),
        "left_inverse": module_map_json(&left_inverse, coefficient),
        "triangles": [
            inclusion.compose(&right_inverse).unwrap()
                == identity_module_mapping(simple_one.clone()).unwrap(),
            left_inverse.compose(&projection).unwrap()
                == identity_module_mapping(simple_one.clone()).unwrap(),
        ],
        "cache_identity": [
            right_inverse_of_homomorphism(&inclusion).unwrap().unwrap()
                .is_identical(&right_inverse),
            left_inverse_of_homomorphism(&projection).unwrap().unwrap()
                .is_identical(&left_inverse),
            left_inverse_of_homomorphism(&right_inverse).unwrap().unwrap()
                .is_identical(&inclusion),
            right_inverse_of_homomorphism(&left_inverse).unwrap().unwrap()
                .is_identical(&projection),
        ],
        "split_properties": [
            is_split_monomorphism(&inclusion).unwrap(),
            is_split_epimorphism(&right_inverse).unwrap(),
            is_split_epimorphism(&projection).unwrap(),
            is_split_monomorphism(&left_inverse).unwrap(),
        ],
        "nonsplit": [
            right_inverse_of_homomorphism(&radical_inclusion).unwrap().is_none(),
            is_split_monomorphism(&radical_inclusion).unwrap(),
            left_inverse_of_homomorphism(&top_projection).unwrap().is_none(),
            is_split_epimorphism(&top_projection).unwrap(),
        ],
        "zero_space_bug": [
            right_inverse_of_homomorphism(&zero_inclusion).unwrap().is_none(),
            is_split_monomorphism(&zero_inclusion).unwrap(),
            left_inverse_of_homomorphism(&zero_projection).unwrap().is_none(),
            is_split_epimorphism(&zero_projection).unwrap(),
        ],
        "identity_right": module_map_json(&identity_right, coefficient),
        "identity_left": module_map_json(&identity_left, coefficient),
        "identity_cache_same": identity_right.is_identical(&identity_left),
        "property_values": property_values,
        "property_right": module_map_json(&property_right, coefficient),
        "property_left": module_map_json(&property_left, coefficient),
        "property_triangles": [
            property_inclusion.compose(&property_right).unwrap()
                == identity_module_mapping(simple_one.clone()).unwrap(),
            property_left.compose(&property_projection).unwrap()
                == identity_module_mapping(simple_one).unwrap(),
        ],
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_map_utils_case<F: Field>(
    field: F,
    integer: impl Fn(i64) -> F::Element + Copy,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let two = right_module_over_path_algebra(algebra.clone(), vec![2], &[])
        .unwrap()
        .unwrap();
    let automorphism = right_module_hom_over_algebra(
        two.clone(),
        two.clone(),
        vec![
            DenseMatrix::from_rows(
                vec![vec![integer(1), integer(1)], vec![integer(0), integer(1)]],
                field.clone(),
            )
            .unwrap(),
        ],
    )
    .unwrap();
    let inverse = inverse_of_isomorphism(&automorphism).unwrap();
    let zero_module = right_module_over_path_algebra(algebra, vec![0], &[])
        .unwrap()
        .unwrap();
    let zero_inverse =
        inverse_of_isomorphism(&zero_module_mapping(zero_module.clone(), zero_module).unwrap())
            .unwrap();

    let quiver_a2 =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra_a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver_a2).unwrap());
    let representation = right_module_over_path_algebra(
        algebra_a2,
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![integer(1)], vec![integer(0)]], field.clone())
                .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let source_map = right_module_hom_over_algebra(
        representation.clone(),
        representation.clone(),
        vec![
            DenseMatrix::from_rows(
                vec![vec![integer(1), integer(1)], vec![integer(0), integer(1)]],
                field.clone(),
            )
            .unwrap(),
            DenseMatrix::identity(1, field.clone()),
        ],
    )
    .unwrap();
    let basis = representation.basis_elements();
    let images = basis
        .iter()
        .map(|element| source_map.image(element).unwrap())
        .collect::<Vec<_>>();
    let reconstructed =
        homomorphism_from_images(&representation, &representation, &images).unwrap();
    let mut images_with_extra = images;
    images_with_extra.push(representation.zero_element());
    let reconstructed_extra =
        homomorphism_from_images(&representation, &representation, &images_with_extra).unwrap();

    let make_map = |rows: [[i64; 2]; 2]| {
        right_module_hom_over_algebra(
            two.clone(),
            two.clone(),
            vec![
                DenseMatrix::from_rows(
                    rows.into_iter()
                        .map(|row| row.into_iter().map(integer).collect())
                        .collect(),
                    field.clone(),
                )
                .unwrap(),
            ],
        )
        .unwrap()
    };
    let p1 = make_map([[1, 1], [0, 1]]);
    let p2 = make_map([[1, 0], [1, 1]]);
    let map_matrix = vec![
        vec![make_map([[1, 0], [0, 0]]), make_map([[1, 0], [0, 1]])],
        vec![make_map([[0, 1], [1, 0]]), make_map([[1, 1], [1, 0]])],
    ];
    let i1 = make_map([[1, 0], [1, 1]]);
    let i2 = make_map([[0, 1], [1, 0]]);
    let product = multiply_lists_of_maps(
        &[p1.clone(), p2.clone()],
        &map_matrix,
        &[i1.clone(), i2.clone()],
    )
    .unwrap();
    let mut map_matrix_extra = map_matrix;
    map_matrix_extra.push(vec![p1.clone(), p2.clone()]);
    let product_extra = multiply_lists_of_maps(&[p1, p2], &map_matrix_extra, &[i1, i2]).unwrap();

    json!({
        "inverse": module_map_json(&inverse, coefficient),
        "inverse_triangles": [
            automorphism.compose(&inverse).unwrap()
                == identity_module_mapping(two.clone()).unwrap(),
            inverse.compose(&automorphism).unwrap()
                == identity_module_mapping(two).unwrap(),
        ],
        "zero_inverse": module_map_json(&zero_inverse, coefficient),
        "from_images": module_map_json(&reconstructed, coefficient),
        "from_images_equal": reconstructed == source_map,
        "extra_image_ignored": reconstructed_extra == source_map,
        "list_product": module_map_json(&product, coefficient),
        "extra_row_ignored": product_extra == product,
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_isomorphism_case<F: Field>(
    field: F,
    integer: impl Fn(i64) -> F::Element + Copy,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let one = right_module_over_path_algebra(algebra.clone(), vec![1], &[])
        .unwrap()
        .unwrap();
    let two = right_module_over_path_algebra(algebra.clone(), vec![2], &[])
        .unwrap()
        .unwrap();
    let equal_one = right_module_over_path_algebra(algebra.clone(), vec![2], &[])
        .unwrap()
        .unwrap();
    let equal_two = right_module_over_path_algebra(algebra.clone(), vec![2], &[])
        .unwrap()
        .unwrap();
    let zero_module = right_module_over_path_algebra(algebra, vec![0], &[])
        .unwrap()
        .unwrap();
    let zero_isomorphism = isomorphism_of_modules(&zero_module, &zero_module)
        .unwrap()
        .unwrap();
    let equal_isomorphism = isomorphism_of_modules(&equal_one, &equal_two)
        .unwrap()
        .unwrap();
    let different_dimension = isomorphism_of_modules(&one, &two).unwrap();

    let other_quiver = Quiver::from_indexed(2, Vec::<ArrowSpec>::new()).unwrap();
    let other_algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), other_quiver).unwrap());
    let other_module = right_module_over_path_algebra(other_algebra, vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let different_algebra = isomorphism_of_modules(&one, &other_module).unwrap();

    let quiver_a2 =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra_a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver_a2).unwrap());
    let make_representation = |rows: [[i64; 2]; 2]| {
        right_module_over_path_algebra(
            algebra_a2.clone(),
            vec![2, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(
                    rows.into_iter()
                        .map(|row| row.into_iter().map(integer).collect())
                        .collect(),
                    field.clone(),
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap()
    };
    let rank_one = make_representation([[1, 0], [0, 0]]);
    let conjugate_rank_one = make_representation([[1, 1], [1, 1]]);
    let rank_two = make_representation([[1, 0], [0, 1]]);
    let recursive = isomorphism_of_modules(&rank_one, &conjugate_rank_one)
        .unwrap()
        .unwrap();
    let reverse = isomorphism_of_modules(&conjugate_rank_one, &rank_one)
        .unwrap()
        .unwrap();
    let recursive_inverse = inverse_of_isomorphism(&recursive).unwrap();
    let nonisomorphic = isomorphism_of_modules(&rank_one, &rank_two).unwrap();

    json!({
        "zero": module_map_json(&zero_isomorphism, coefficient),
        "equal": module_map_json(&equal_isomorphism, coefficient),
        "equal_endpoints_identical": [
            equal_isomorphism.source().is_identical(&equal_one),
            equal_isomorphism.range().is_identical(&equal_two),
        ],
        "different_dimension": different_dimension.is_none(),
        "different_algebra": different_algebra.is_none(),
        "recursive": module_map_json(&recursive, coefficient),
        "recursive_isomorphism": recursive.is_isomorphism(),
        "recursive_triangles": [
            recursive.compose(&recursive_inverse).unwrap()
                == identity_module_mapping(rank_one.clone()).unwrap(),
            recursive_inverse.compose(&recursive).unwrap()
                == identity_module_mapping(conjugate_rank_one.clone()).unwrap(),
        ],
        "reverse": module_map_json(&reverse, coefficient),
        "nonisomorphic": nonisomorphic.is_none(),
    })
}

fn more_minimal_json<F: Field>(
    result: &MoreMinimalVersionResult<F>,
    original: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    match result {
        MoreMinimalVersionResult::Unchanged(homomorphism) => json!({
            "changed": false,
            "same_object": homomorphism.is_identical(original),
            "map": module_map_json(homomorphism, coefficient),
        }),
        MoreMinimalVersionResult::Reduced {
            minimal_part,
            removed_part,
        } => json!({
            "changed": true,
            "minimal_part": module_map_json(minimal_part, coefficient),
            "removed_part": module_map_json(removed_part, coefficient),
        }),
    }
}

fn minimal_version_json<F: Field>(
    result: &(
        PathAlgebraModuleHomomorphism<F>,
        Vec<super::PathAlgebraModule<F>>,
    ),
    original: &PathAlgebraModuleHomomorphism<F>,
    right: bool,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let property = if right {
        result.0.cached_is_right_minimal()
    } else {
        result.0.cached_is_left_minimal()
    };
    json!({
        "map": module_map_json(&result.0, coefficient),
        "removed_dimensions": result.1.iter()
            .map(super::PathAlgebraModule::dimension_vector).collect::<Vec<_>>(),
        "map_is_original": result.0.is_identical(original),
        "property_bound": property.is_some(),
        "property_value": property.unwrap_or(false),
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_minimal_case<F: Field>(
    field: F,
    integer: impl Fn(i64) -> F::Element + Copy,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let one = right_module_over_path_algebra(algebra.clone(), vec![1], &[])
        .unwrap()
        .unwrap();
    let three = right_module_over_path_algebra(algebra, vec![3], &[])
        .unwrap()
        .unwrap();
    let right_matrix = || {
        DenseMatrix::from_rows(
            vec![vec![integer(1)], vec![integer(0)], vec![integer(0)]],
            field.clone(),
        )
        .unwrap()
    };
    let left_matrix = || {
        DenseMatrix::from_rows(
            vec![vec![integer(1), integer(0), integer(0)]],
            field.clone(),
        )
        .unwrap()
    };
    let right_more_map =
        right_module_hom_over_algebra(three.clone(), one.clone(), vec![right_matrix()]).unwrap();
    let right_more = more_right_minimal_version(&right_more_map).unwrap();
    let left_more_map =
        right_module_hom_over_algebra(one.clone(), three.clone(), vec![left_matrix()]).unwrap();
    let left_more = more_left_minimal_version(&left_more_map).unwrap();

    let right_map =
        right_module_hom_over_algebra(three.clone(), one.clone(), vec![right_matrix()]).unwrap();
    let right_result = right_minimal_version(&right_map).unwrap();
    let right_again = right_minimal_version(&right_map).unwrap();
    let left_map = right_module_hom_over_algebra(one.clone(), three, vec![left_matrix()]).unwrap();
    let left_result = left_minimal_version(&left_map).unwrap();
    let left_again = left_minimal_version(&left_map).unwrap();
    let identity = identity_module_mapping(one.clone()).unwrap();
    let right_more_unchanged = more_right_minimal_version(&identity).unwrap();
    let left_more_unchanged = more_left_minimal_version(&identity).unwrap();
    let right_identity = right_minimal_version(&identity).unwrap();
    let left_identity = left_minimal_version(&identity).unwrap();
    let zero_map = zero_module_mapping(one.clone(), one).unwrap();
    let right_zero = right_minimal_version(&zero_map).unwrap();

    let quiver_a2 =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra_a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver_a2).unwrap());
    let projective_one = right_module_over_path_algebra(
        algebra_a2.clone(),
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![integer(1)]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let simple_two = right_module_over_path_algebra(algebra_a2, vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let direct_sum = direct_sum_of_qpa_modules(&[projective_one, simple_two])
        .unwrap()
        .unwrap();
    let projection = direct_sum_projections(&direct_sum).unwrap()[0].clone();
    let inclusion = direct_sum_inclusions(&direct_sum).unwrap()[0].clone();
    let a2_right = right_minimal_version(&projection).unwrap();
    let a2_left = left_minimal_version(&inclusion).unwrap();

    json!({
        "right_more": more_minimal_json(&right_more, &right_more_map, coefficient),
        "left_more": more_minimal_json(&left_more, &left_more_map, coefficient),
        "right_more_unchanged": more_minimal_json(
            &right_more_unchanged, &identity, coefficient),
        "left_more_unchanged": more_minimal_json(
            &left_more_unchanged, &identity, coefficient),
        "right_minimal": minimal_version_json(&right_result, &right_map, true, coefficient),
        "left_minimal": minimal_version_json(&left_result, &left_map, false, coefficient),
        "attribute_cache_identity": [
            right_result.0.is_identical(&right_again.0),
            right_result.1.iter().zip(&right_again.1)
                .all(|(left, right)| left.is_identical(right)),
            left_result.0.is_identical(&left_again.0),
            left_result.1.iter().zip(&left_again.1)
                .all(|(left, right)| left.is_identical(right)),
        ],
        "right_identity": minimal_version_json(&right_identity, &identity, true, coefficient),
        "left_identity": minimal_version_json(&left_identity, &identity, false, coefficient),
        "right_zero": minimal_version_json(&right_zero, &zero_map, true, coefficient),
        "a2_right": minimal_version_json(&a2_right, &projection, true, coefficient),
        "a2_left": minimal_version_json(&a2_left, &inclusion, false, coefficient),
    })
}

fn common_summand_dimensions<F: Field>(result: CommonDirectSummandResult<F>) -> Value {
    let CommonDirectSummandResult::Found(common) = result else {
        panic!("the pinned QPA case has a common direct summand")
    };
    json!([
        common.source_common().dimension_vector(),
        common.source_complement().dimension_vector(),
        common.range_common().dimension_vector(),
        common.range_complement().dimension_vector(),
    ])
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_common_summand_case<F: QpaCoefficientField>(field: F) -> Value {
    let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let zero_module = right_module_over_path_algebra(algebra.clone(), vec![0], &[])
        .unwrap()
        .unwrap();
    let one = right_module_over_path_algebra(algebra.clone(), vec![1], &[])
        .unwrap()
        .unwrap();
    let two = right_module_over_path_algebra(algebra, vec![2], &[])
        .unwrap()
        .unwrap();
    let two_copy = direct_sum_of_qpa_modules(&[one.clone(), one.clone()])
        .unwrap()
        .unwrap();
    let common = common_summand_dimensions(common_direct_summand(&two, &one).unwrap());
    let maximal = maximal_common_direct_summand(&two, &one).unwrap().unwrap();

    let quiver_a2 =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra_a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver_a2).unwrap());
    let projective_one = right_module_over_path_algebra(
        algebra_a2.clone(),
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![field.one()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let simple_one = right_module_over_path_algebra(algebra_a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let simple_two = right_module_over_path_algebra(algebra_a2, vec![0, 1], &[])
        .unwrap()
        .unwrap();
    assert!(matches!(
        common_direct_summand(&simple_one, &simple_two).unwrap(),
        CommonDirectSummandResult::NoCommonSummand
    ));
    let mixed = direct_sum_of_qpa_modules(&[simple_one.clone(), simple_two.clone()])
        .unwrap()
        .unwrap();
    let mixed_common =
        common_summand_dimensions(common_direct_summand(&mixed, &simple_one).unwrap());

    let identity = identity_module_mapping(one.clone()).unwrap();
    let projection = right_module_hom_over_algebra(
        two.clone(),
        one.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![field.one()], vec![field.zero()]], field.clone())
                .unwrap(),
        ],
    )
    .unwrap();
    let inclusion = right_module_hom_over_algebra(
        one.clone(),
        two.clone(),
        vec![DenseMatrix::from_rows(vec![vec![field.one(), field.zero()]], field.clone()).unwrap()],
    )
    .unwrap();
    let zero = zero_module_mapping(one.clone(), one.clone()).unwrap();
    let top_projection = top_of_module_projection(projective_one.clone()).unwrap();
    let radical_inclusion = radical_of_module_inclusion(projective_one).unwrap();
    let right_before = projection.cached_is_right_minimal().is_some();
    let left_before = inclusion.cached_is_left_minimal().is_some();

    json!({
        "common": common,
        "maximal": {
            "common": maximal.common_summands().iter()
                .map(super::PathAlgebraModule::dimension_vector).collect::<Vec<_>>(),
            "source_remainder": maximal.source_remainder().dimension_vector(),
            "range_remainder": maximal.range_remainder().dimension_vector(),
        },
        "disjoint_is_false": true,
        "zero_common_is_false": matches!(
            common_direct_summand(&zero_module, &zero_module).unwrap(),
            CommonDirectSummandResult::NoCommonSummand),
        "mixed_common": mixed_common,
        "decisions": {
            "isomorphic": isomorphic_modules(&two, &two_copy).unwrap(),
            "not_isomorphic": isomorphic_modules(&simple_one, &simple_two).unwrap(),
            "direct_summand": is_direct_summand(&one, &two).unwrap(),
            "not_direct_summand": is_direct_summand(&simple_one, &simple_two).unwrap(),
            "additive_vector_space": is_in_additive_closure(&two, &one).unwrap().unwrap(),
            "not_additive": is_in_additive_closure(&simple_one, &simple_two)
                .unwrap().unwrap(),
        },
        "minimality": {
            "right_before": right_before,
            "right_projection": is_right_minimal(&projection).unwrap(),
            "right_after": projection.cached_is_right_minimal().is_some(),
            "left_before": left_before,
            "left_inclusion": is_left_minimal(&inclusion).unwrap(),
            "left_after": inclusion.cached_is_left_minimal().is_some(),
            "identity": [
                is_right_minimal(&identity).unwrap(),
                is_left_minimal(&identity).unwrap(),
            ],
            "zero": [
                is_right_minimal(&zero).unwrap(),
                is_left_minimal(&zero).unwrap(),
            ],
            "projective_cover_right": is_right_minimal(&top_projection).unwrap(),
            "radical_inclusion_left": is_left_minimal(&radical_inclusion).unwrap(),
        }
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_categorical_case<F: Field>(
    field: F,
    integer: impl Fn(i64) -> F::Element + Copy,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let module = right_module_over_path_algebra(
        algebra.clone(),
        vec![2, 2],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(
                vec![vec![integer(1), integer(1)], vec![integer(0), integer(0)]],
                field.clone(),
            )
            .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let f = right_module_hom_over_algebra(
        module.clone(),
        module.clone(),
        vec![
            DenseMatrix::from_rows(
                vec![vec![integer(1), integer(0)], vec![integer(0), integer(0)]],
                field.clone(),
            )
            .unwrap(),
            DenseMatrix::identity(2, field.clone()),
        ],
    )
    .unwrap();
    let t = right_module_hom_over_algebra(
        module.clone(),
        module.clone(),
        vec![
            DenseMatrix::from_rows(
                vec![vec![integer(2), integer(0)], vec![integer(0), integer(3)]],
                field.clone(),
            )
            .unwrap(),
            DenseMatrix::identity(2, field.clone()).scale_left(&integer(2)),
        ],
    )
    .unwrap();
    let identity = identity_module_mapping(module.clone()).unwrap();
    let layer = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
        .unwrap()
        .unwrap();
    let first = right_module_hom_over_algebra(
        layer.clone(),
        module.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![integer(0), integer(1)]], field.clone()).unwrap(),
            DenseMatrix::from_rows(vec![vec![integer(1), integer(0)]], field.clone()).unwrap(),
        ],
    )
    .unwrap();
    let second = right_module_hom_over_algebra(
        layer.clone(),
        module.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![integer(0), integer(1)]], field.clone()).unwrap(),
            DenseMatrix::from_rows(vec![vec![integer(0), integer(1)]], field.clone()).unwrap(),
        ],
    )
    .unwrap();
    let kernel_module = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let kernel_to_module = right_module_hom_over_algebra(
        kernel_module.clone(),
        module.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![integer(0), integer(1)]], field.clone()).unwrap(),
            DenseMatrix::zero(1, 2, field.clone()),
        ],
    )
    .unwrap();
    let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let simple_two = right_module_over_path_algebra(algebra, vec![0, 1], &[])
        .unwrap()
        .unwrap();

    let direct_sum = direct_sum_of_qpa_modules(&[layer.clone(), kernel_module, simple_two.clone()])
        .unwrap()
        .unwrap();
    let projections = direct_sum_projections(&direct_sum).unwrap();
    let inclusions = direct_sum_inclusions(&direct_sum).unwrap();
    let singleton = direct_sum_of_qpa_modules(std::slice::from_ref(&layer))
        .unwrap()
        .unwrap();
    let induced_kernel = morphism_on_kernel(&f, &f, &t, &t).unwrap().unwrap();
    let induced_image = morphism_on_image(&f, &f, &t, &t).unwrap().unwrap();
    let induced_cokernel = morphism_on_cokernel(&f, &f, &t, &t).unwrap().unwrap();
    let intersection = intersection_of_submodules(&first, &second).unwrap();
    let intersection_list =
        intersection_of_submodules_list(&[first.clone(), second.clone(), kernel_to_module.clone()])
            .unwrap();
    let sum = sum_of_submodules(&first, &second).unwrap();
    let sum_list =
        sum_of_submodules_list(&[first.clone(), second.clone(), kernel_to_module]).unwrap();
    let pullback = pull_back(&first, &second).unwrap().unwrap();
    let pushout = push_out(&intersection.0, &intersection.1).unwrap().unwrap();
    let reject_zero = reject_of_module(&simple_one, &simple_two).unwrap();
    let reject_one = reject_of_module(&module, &simple_two).unwrap();
    let reject_many = reject_of_module(&module, &simple_one).unwrap();

    json!({
        "direct_sum": {
            "empty_is_fail": direct_sum_of_qpa_modules::<F>(&[]).unwrap().is_none(),
            "filter": direct_sum.is_direct_sum_of_modules(),
            "dimension": direct_sum.dimension_vector(),
            "matrices": direct_sum.arrow_matrices().iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
            "projections": projections.iter()
                .map(|map| module_map_json(map, coefficient)).collect::<Vec<_>>(),
            "inclusions": inclusions.iter()
                .map(|map| module_map_json(map, coefficient)).collect::<Vec<_>>(),
            "projection_cache_identity": projections[0]
                .is_identical(&direct_sum_projections(&direct_sum).unwrap()[0]),
            "inclusion_cache_identity": inclusions[0]
                .is_identical(&direct_sum_inclusions(&direct_sum).unwrap()[0]),
            "singleton_equal": singleton == layer,
            "singleton_identical": singleton.is_identical(&layer),
        },
        "induced": {
            "kernel": module_map_json(&induced_kernel, coefficient),
            "image": module_map_json(&induced_image, coefficient),
            "cokernel": module_map_json(&induced_cokernel, coefficient),
            "noncommuting_fails": [
                morphism_on_kernel(&f, &f, &identity, &t).unwrap().is_none(),
                morphism_on_image(&f, &f, &identity, &t).unwrap().is_none(),
                morphism_on_cokernel(&f, &f, &identity, &t).unwrap().is_none(),
            ],
        },
        "pullback": module_map_list_json(&[&pullback.0, &pullback.1], coefficient),
        "pushout": module_map_list_json(&[&pushout.0, &pushout.1], coefficient),
        "intersection": {
            "pair": module_map_list_json(
                &[&intersection.0, &intersection.1, &intersection.2], coefficient),
            "list": module_map_json(&intersection_list, coefficient),
            "singleton_identity": intersection_of_submodules_list(
                std::slice::from_ref(&first)).unwrap().is_identical(&first),
        },
        "sum": {
            "pair": module_map_list_json(&[&sum.0, &sum.1, &sum.2], coefficient),
            "list": module_map_json(&sum_list, coefficient),
            "singleton_identity": sum_of_submodules_list(
                std::slice::from_ref(&first)).unwrap().is_identical(&first),
        },
        "reject": {
            "hom_zero": module_map_json(&reject_zero, coefficient),
            "hom_one": module_map_json(&reject_one, coefficient),
            "hom_many": module_map_json(&reject_many, coefficient),
        },
    })
}

#[test]
fn checked_constructor_sets_omitted_zero_maps_in_arrow_order() {
    let quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "a"),
            ArrowSpec::named(vertex(2), vertex(3), "b"),
        ],
    )
    .unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let a = DenseMatrix::from_rows(vec![vec![q(1), q(2)]], Q).unwrap();
    let module = right_module_over_path_algebra(
        algebra.clone(),
        vec![1, 2, 0],
        &[NamedArrowMatrix::new("a", a.clone())],
    )
    .unwrap()
    .unwrap();

    assert_eq!(module.right_acting_algebra(), &algebra);
    assert_eq!(module.dimension_vector(), &[1, 2, 0]);
    assert_eq!(module.dimension(), 3);
    assert_eq!(module.arrow_matrices()[0], a);
    assert_eq!(module.arrow_matrices()[1], DenseMatrix::zero(2, 1, Q));
    assert_eq!(
        module.matrices_of_path_algebra_module(),
        vec![
            PathAlgebraModuleMatrix::Matrix(module.arrow_matrices()[0].clone()),
            PathAlgebraModuleMatrix::Matrix(module.arrow_matrices()[1].clone()),
        ]
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn dimension_inferring_constructors_match_qpa_oracle() {
    let chain_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "a"),
            ArrowSpec::named(vertex(2), vertex(3), "b"),
        ],
    )
    .unwrap();
    let chain = QuiverAlgebra::Path(PathAlgebra::new(Q, chain_quiver).unwrap());
    let a_matrix = DenseMatrix::from_rows(vec![vec![q(1), q(2)]], Q).unwrap();
    let b_matrix = DenseMatrix::from_rows(vec![vec![q(3)], vec![q(4)]], Q).unwrap();
    let named = right_module_over_path_algebra_infer_dimension(
        chain.clone(),
        &PathAlgebraModuleArrowActions::Named(vec![
            NamedPathAlgebraModuleArrowAction::new(
                "b",
                PathAlgebraModuleArrowAction::matrix(b_matrix.clone()),
            ),
            NamedPathAlgebraModuleArrowAction::new(
                "a",
                PathAlgebraModuleArrowAction::matrix(a_matrix.clone()),
            ),
        ]),
    )
    .unwrap()
    .unwrap();
    let positional = right_module_over_path_algebra_infer_dimension(
        chain.clone(),
        &PathAlgebraModuleArrowActions::Positional(vec![a_matrix.clone(), b_matrix.clone()]),
    )
    .unwrap()
    .unwrap();
    let indexed = right_module_over_path_algebra_infer_dimension(
        chain.clone(),
        &PathAlgebraModuleArrowActions::Indexed(vec![
            IndexedPathAlgebraModuleArrowAction::new(
                ArrowIndex::new(2).unwrap(),
                PathAlgebraModuleArrowAction::matrix(b_matrix.clone()),
            ),
            IndexedPathAlgebraModuleArrowAction::new(
                ArrowIndex::new(1).unwrap(),
                PathAlgebraModuleArrowAction::matrix(a_matrix.clone()),
            ),
        ]),
    )
    .unwrap()
    .unwrap();
    let unchecked = right_module_over_path_algebra_infer_dimension_nc(
        chain,
        &PathAlgebraModuleArrowActions::Positional(vec![a_matrix, b_matrix]),
    )
    .unwrap()
    .unwrap();

    let zero_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "z")]).unwrap();
    let zero_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, zero_quiver).unwrap());
    let zero_action = |source, target| {
        PathAlgebraModuleArrowActions::Named(vec![NamedPathAlgebraModuleArrowAction::new(
            "z",
            PathAlgebraModuleArrowAction::zero_space_dimensions(source, target),
        )])
    };
    let zero_source =
        right_module_over_path_algebra_infer_dimension(zero_algebra.clone(), &zero_action(0, 2))
            .unwrap()
            .unwrap();
    let zero_target =
        right_module_over_path_algebra_infer_dimension(zero_algebra.clone(), &zero_action(2, 0))
            .unwrap()
            .unwrap();
    let zero = right_module_over_path_algebra_infer_dimension(zero_algebra, &zero_action(0, 0))
        .unwrap()
        .unwrap();

    let isolated_algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(3, vec![ArrowSpec::named(vertex(1), vertex(2), "c")]).unwrap(),
        )
        .unwrap(),
    );
    let isolated = right_module_over_path_algebra_infer_dimension(
        isolated_algebra,
        &PathAlgebraModuleArrowActions::Positional(vec![
            DenseMatrix::from_rows(vec![vec![q(5)]], Q).unwrap(),
        ]),
    )
    .unwrap()
    .unwrap();

    let prime_field = Z2;
    let prime_algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            prime_field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "d")]).unwrap(),
        )
        .unwrap(),
    );
    let prime = right_module_over_path_algebra_infer_dimension(
        prime_algebra,
        &PathAlgebraModuleArrowActions::Positional(vec![
            DenseMatrix::from_rows(vec![vec![prime_field.one()]], prime_field.clone()).unwrap(),
        ]),
    )
    .unwrap()
    .unwrap();

    let quotient_parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(3), "y"),
            ],
        )
        .unwrap(),
    )
    .unwrap();
    let generators = quotient_parent.generators();
    let relation = generators[3].multiply(&generators[4]).unwrap();
    let quotient = QuiverAlgebra::Quotient(
        QuotientPathAlgebra::by_relations_high_level(quotient_parent, &[relation]).unwrap(),
    );
    let relation_actions = |second| {
        PathAlgebraModuleArrowActions::Positional(vec![
            DenseMatrix::from_rows(vec![vec![q(1)]], Q).unwrap(),
            DenseMatrix::from_rows(vec![vec![q(second)]], Q).unwrap(),
        ])
    };
    let relation_ok =
        right_module_over_path_algebra_infer_dimension(quotient.clone(), &relation_actions(0))
            .unwrap()
            .unwrap();
    let relation_fail =
        right_module_over_path_algebra_infer_dimension(quotient, &relation_actions(1))
            .unwrap()
            .is_none();

    let rational = |entry: &Rational| entry.to_string();
    let actual = format!(
        "{}{}{}{}{}{}{}{}{}The entered matrices for the module do not satisfy the relation(s).\n{}RELATION_FAIL={}\n",
        qpa_module_data("NAMED", &named, rational),
        qpa_module_data("POSITIONAL", &positional, rational),
        qpa_module_data("INDEXED", &indexed, rational),
        qpa_module_data("UNCHECKED", &unchecked, rational),
        qpa_module_data("ZERO_SOURCE", &zero_source, rational),
        qpa_module_data("ZERO_TARGET", &zero_target, rational),
        qpa_module_data("ZERO_MODULE", &zero, rational),
        qpa_module_data("ISOLATED_CHECKED", &isolated, rational),
        qpa_module_data("GF2", &prime, |entry| {
            if prime_field.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        }),
        qpa_module_data("RELATION_OK", &relation_ok, rational),
        relation_fail,
    );
    assert_eq!(
        actual,
        include_str!("../../../../tests/fixtures/module_inferred_constructors.txt")
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn dimension_inferring_constructors_preserve_qpa_failures() {
    let chain = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(
                3,
                vec![
                    ArrowSpec::named(vertex(1), vertex(2), "a"),
                    ArrowSpec::named(vertex(2), vertex(3), "b"),
                ],
            )
            .unwrap(),
        )
        .unwrap(),
    );
    let one = DenseMatrix::from_rows(vec![vec![q(1)]], Q).unwrap();
    assert_eq!(
        right_module_over_path_algebra_infer_dimension(
            chain.clone(),
            &PathAlgebraModuleArrowActions::Positional(vec![one.clone()]),
        )
        .unwrap_err(),
        PathAlgebraModuleError::IncompleteInferredArrowActions {
            actual: 1,
            expected: 2,
        }
    );
    assert!(matches!(
        right_module_over_path_algebra_infer_dimension(
            chain.clone(),
            &PathAlgebraModuleArrowActions::Positional(vec![
                DenseMatrix::zero(1, 2, Q),
                DenseMatrix::zero(3, 1, Q),
            ]),
        ),
        Err(
            PathAlgebraModuleError::IncompatibleInferredVertexDimension {
                vertex: 2,
                entered: 3,
                expected: 2,
                ..
            }
        )
    ));
    assert!(matches!(
        right_module_over_path_algebra_infer_dimension(
            chain.clone(),
            &PathAlgebraModuleArrowActions::Named(vec![
                NamedPathAlgebraModuleArrowAction::new(
                    "a",
                    PathAlgebraModuleArrowAction::zero_space_dimensions(2, 2),
                ),
                NamedPathAlgebraModuleArrowAction::new(
                    "b",
                    PathAlgebraModuleArrowAction::matrix(DenseMatrix::zero(2, 1, Q)),
                ),
            ]),
        ),
        Err(PathAlgebraModuleError::InvalidZeroSpaceDimensions { .. })
    ));
    assert!(matches!(
        right_module_over_path_algebra_infer_dimension(
            chain.clone(),
            &PathAlgebraModuleArrowActions::Named(vec![
                NamedPathAlgebraModuleArrowAction::new(
                    "a",
                    PathAlgebraModuleArrowAction::matrix(one.clone()),
                ),
                NamedPathAlgebraModuleArrowAction::new(
                    "a",
                    PathAlgebraModuleArrowAction::matrix(one.clone()),
                ),
            ]),
        ),
        Err(PathAlgebraModuleError::UnassignedInferredArrowActions { .. })
    ));
    assert!(matches!(
        right_module_over_path_algebra_infer_dimension(
            chain.clone(),
            &PathAlgebraModuleArrowActions::Indexed(vec![
                IndexedPathAlgebraModuleArrowAction::new(
                    ArrowIndex::new(3).unwrap(),
                    PathAlgebraModuleArrowAction::matrix(one.clone()),
                ),
                IndexedPathAlgebraModuleArrowAction::new(
                    ArrowIndex::new(1).unwrap(),
                    PathAlgebraModuleArrowAction::matrix(one.clone()),
                ),
            ]),
        ),
        Err(PathAlgebraModuleError::UnknownArrowIndex { index: 3 })
    ));
    assert!(matches!(
        right_module_over_path_algebra_infer_dimension(
            chain,
            &PathAlgebraModuleArrowActions::Positional(vec![
                DenseMatrix::from_flat(0, 1, Vec::new(), Q).unwrap(),
                one,
            ]),
        ),
        Err(PathAlgebraModuleError::ZeroDimensionMatrixInInferredConstructor { .. })
    ));

    let isolated = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(3, vec![ArrowSpec::named(vertex(1), vertex(2), "c")]).unwrap(),
        )
        .unwrap(),
    );
    assert_eq!(
        right_module_over_path_algebra_infer_dimension_nc(
            isolated,
            &PathAlgebraModuleArrowActions::Positional(vec![
                DenseMatrix::from_rows(vec![vec![q(1)]], Q).unwrap(),
            ]),
        )
        .unwrap_err(),
        PathAlgebraModuleError::UncheckedInferredIsolatedVertices { vertices: vec![3] }
    );

    let arrowless = QuiverAlgebra::Path(
        PathAlgebra::new(Q, Quiver::from_indexed(1, Vec::new()).unwrap()).unwrap(),
    );
    assert_eq!(
        right_module_over_path_algebra_infer_dimension(
            arrowless,
            &PathAlgebraModuleArrowActions::Positional(Vec::new()),
        )
        .unwrap_err(),
        PathAlgebraModuleError::EmptyInferredArrowActions
    );

    let prime_three = Zp::new(3);
    let prime_five = Zp::new(5);
    let prime_algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            prime_three,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "d")]).unwrap(),
        )
        .unwrap(),
    );
    let foreign = DenseMatrix::zero(1, 1, prime_five);
    assert!(matches!(
        right_module_over_path_algebra_infer_dimension(
            prime_algebra.clone(),
            &PathAlgebraModuleArrowActions::Positional(vec![foreign.clone()]),
        ),
        Err(PathAlgebraModuleError::DifferentCoefficientField { .. })
    ));
    let unchecked = right_module_over_path_algebra_infer_dimension_nc(
        prime_algebra,
        &PathAlgebraModuleArrowActions::Positional(vec![foreign.clone()]),
    )
    .unwrap()
    .unwrap();
    assert_eq!(unchecked.arrow_matrices()[0], foreign);
}

#[test]
fn matrix_submodules_and_natural_basis_mapping_match_qpa_oracle() {
    fn run<F: Field>(field: F, coefficient: impl Fn(&F::Element) -> String + Copy) -> String {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let module = right_module_over_path_algebra(
            algebra,
            vec![3, 3],
            &[NamedArrowMatrix::new("a", DenseMatrix::identity(3, field))],
        )
        .unwrap()
        .unwrap();
        let basis = module.basis_elements();
        let generator = basis[0].add(&basis[4]).unwrap();
        let submodule =
            subalgebra_module(module.clone(), &[generator, module.zero_element()]).unwrap();
        let canonical = submodule.canonical_basis_elements().unwrap();
        let standalone = submodule.as_module().unwrap();
        let mapping = natural_homomorphism_by_subalgebra_module(&module, &submodule).unwrap();
        assert!(mapping.is_surjective());
        for (source_basis, image) in basis.iter().zip(mapping.basis_images()) {
            assert_eq!(mapping.image(source_basis).unwrap(), *image);
        }
        super::qpa_list([
            qpa_element_list(submodule.generators(), coefficient),
            qpa_element_list(&canonical, coefficient),
            super::qpa_usize_list(standalone.dimension_vector()),
            qpa_module_matrix_list(&standalone, coefficient),
            super::qpa_usize_list(mapping.range().dimension_vector()),
            qpa_module_matrix_list(mapping.range(), coefficient),
            qpa_element_list(mapping.basis_images(), coefficient),
            mapping.is_surjective().to_string(),
        ])
    }

    let rational = run(Q, ToString::to_string);
    let binary = run(Z2, |entry| {
        if Z2.is_zero(entry) {
            "0*Z(2)".to_owned()
        } else {
            "Z(2)^0".to_owned()
        }
    });
    let actual = format!("QQ={rational}\nGF2={binary}\n");
    assert_eq!(
        actual,
        include_str!("../../../../tests/fixtures/module_submodules.txt")
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn generic_path_ring_hom_and_end_match_qpa_hopf_oracle() {
    fn run<F: Field>(field: F) -> String {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::identity(2, field.clone()),
            )],
        )
        .unwrap()
        .unwrap();
        let hom =
            generic_hom_over_path_ring(&algebra, module.clone().into(), module.clone().into())
                .unwrap();
        let end = generic_end_over_path_ring(&algebra, module.clone().into()).unwrap();
        let basis = module.basis_elements();
        let generator = basis[0].add(&basis[3]).unwrap();
        let submodule = subalgebra_module(module.clone(), &[generator]).unwrap();
        let canonical = submodule.canonical_basis_elements().unwrap();
        let sub_hom = generic_hom_over_path_ring(
            &algebra,
            submodule.clone().into(),
            submodule.clone().into(),
        )
        .unwrap();
        let sub_end = generic_end_over_path_ring(&algebra, submodule.clone().into()).unwrap();
        let to_submodule =
            generic_hom_over_path_ring(&algebra, module.clone().into(), submodule.clone().into())
                .unwrap();
        let from_submodule =
            generic_hom_over_path_ring(&algebra, submodule.into(), module.clone().into()).unwrap();

        assert!(
            hom.generators()
                .iter()
                .all(GenericPathAlgebraModuleHomomorphism::is_algebra_module_homomorphism)
        );
        assert_eq!(end.one().matrix(), &DenseMatrix::identity(4, field.clone()));
        assert_eq!(
            end.multiply(&end.generators()[1], &end.generators()[2])
                .unwrap()
                .matrix(),
            end.generators()[0].matrix()
        );
        assert_eq!(
            end.multiply(&end.generators()[2], &end.generators()[1])
                .unwrap()
                .matrix(),
            end.generators()[3].matrix()
        );
        assert_eq!(sub_end.one().matrix(), &DenseMatrix::identity(3, field));
        assert_eq!(
            sub_hom.generators()[0].image(&canonical[0]).unwrap(),
            canonical[0]
        );
        assert_eq!(
            sub_hom.generators()[0].image(&canonical[1]).unwrap(),
            canonical[1]
        );
        assert!(
            sub_hom.generators()[0]
                .image(&canonical[2])
                .unwrap()
                .is_zero()
        );

        let canonical_support = super::qpa_list(canonical.iter().map(|element| {
            let positions = element
                .basis_coefficients()
                .iter()
                .enumerate()
                .filter_map(|(position, entry)| {
                    (!module.right_acting_algebra().field().is_zero(entry)).then_some(position + 1)
                })
                .collect::<Vec<_>>();
            super::qpa_usize_list(&positions)
        }));
        super::qpa_list([
            hom.dimension().to_string(),
            qpa_generic_hom_supports(hom.generators()),
            super::qpa_list(
                hom.generators()
                    .iter()
                    .map(|mapping| mapping.is_algebra_module_homomorphism().to_string()),
            ),
            end.dimension().to_string(),
            qpa_generic_hom_supports(end.generators()),
            canonical_support,
            sub_hom.dimension().to_string(),
            qpa_generic_hom_supports(sub_hom.generators()),
            sub_end.dimension().to_string(),
            qpa_generic_hom_supports(sub_end.generators()),
            qpa_generic_hom_supports(to_submodule.generators()),
            qpa_generic_hom_supports(from_submodule.generators()),
        ])
    }

    let actual = format!("QQ={}\nGF2={}\n", run(Q), run(Z2));
    assert_eq!(
        actual,
        include_str!("../../../../tests/fixtures/module_generic_hom.txt")
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn generic_path_ring_hom_failures_match_qpa_dispatch() {
    let algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let module = right_module_over_path_algebra(
        algebra.clone(),
        vec![1, 1],
        &[NamedArrowMatrix::new("a", DenseMatrix::identity(1, Q))],
    )
    .unwrap()
    .unwrap();
    let zero = zero_module(&algebra).unwrap();
    assert_eq!(
        generic_hom_over_path_ring(&algebra, zero.clone().into(), module.clone().into(),)
            .unwrap_err(),
        PathAlgebraModuleError::ZeroGenericHomModuleFailure
    );
    assert_eq!(
        generic_hom_over_path_ring(&algebra, module.clone().into(), zero.clone().into(),)
            .unwrap_err(),
        PathAlgebraModuleError::ZeroGenericHomModuleFailure
    );
    assert_eq!(
        generic_end_over_path_ring(&algebra, zero.into()).unwrap_err(),
        PathAlgebraModuleError::ZeroGenericHomModuleFailure
    );

    let zero_submodule = subalgebra_module(module.clone(), &[]).unwrap();
    assert_eq!(
        generic_hom_over_path_ring(
            &algebra,
            zero_submodule.clone().into(),
            zero_submodule.clone().into(),
        )
        .unwrap_err(),
        PathAlgebraModuleError::ZeroSubmoduleAsModuleFailure
    );
    assert_eq!(
        generic_end_over_path_ring(&algebra, zero_submodule.into()).unwrap_err(),
        PathAlgebraModuleError::ZeroSubmoduleAsModuleFailure
    );

    let arrowless = QuiverAlgebra::Path(
        PathAlgebra::new(Q, Quiver::from_indexed(1, Vec::new()).unwrap()).unwrap(),
    );
    let arrowless_module = right_module_over_path_algebra(arrowless.clone(), vec![2], &[])
        .unwrap()
        .unwrap();
    assert_eq!(
        generic_hom_over_path_ring(
            &arrowless,
            arrowless_module.clone().into(),
            arrowless_module.clone().into(),
        )
        .unwrap_err(),
        PathAlgebraModuleError::DegenerateGenericHomNullspaceFailure
    );
    assert_eq!(
        generic_end_over_path_ring(&arrowless, arrowless_module.into()).unwrap_err(),
        PathAlgebraModuleError::DegenerateGenericHomNullspaceFailure
    );

    let other_algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let other_module = right_module_over_path_algebra(
        other_algebra,
        vec![1, 1],
        &[NamedArrowMatrix::new("a", DenseMatrix::identity(1, Q))],
    )
    .unwrap()
    .unwrap();
    assert_eq!(
        generic_hom_over_path_ring(&algebra, module.clone().into(), other_module.into())
            .unwrap_err(),
        PathAlgebraModuleError::GenericHomActingAlgebraMismatch
    );

    let loop_parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap(),
    )
    .unwrap();
    let t = loop_parent.generators()[1].clone();
    let quotient = QuiverAlgebra::Quotient(
        QuotientPathAlgebra::by_relations_high_level(loop_parent, &[t.pow(2).unwrap()]).unwrap(),
    );
    let quotient_module = right_module_over_path_algebra(
        quotient.clone(),
        vec![2],
        &[NamedArrowMatrix::new(
            "t",
            DenseMatrix::from_rows(vec![vec![q(0), q(1)], vec![q(0), q(0)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    assert_eq!(
        generic_hom_over_path_ring(
            &quotient,
            quotient_module.clone().into(),
            quotient_module.clone().into(),
        )
        .unwrap_err(),
        PathAlgebraModuleError::GenericHomPathRingDispatchFailure
    );
    assert_eq!(
        generic_end_over_path_ring(&quotient, quotient_module.into()).unwrap_err(),
        PathAlgebraModuleError::GenericHomPathRingDispatchFailure
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn based_right_algebra_module_conversion_matches_qpa_oracle() {
    fn data<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> String {
        super::qpa_list([
            super::qpa_usize_list(module.dimension_vector()),
            qpa_module_matrix_list(module, coefficient),
        ])
    }

    fn run<F: QpaCoefficientField>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> String {
        let zero = field.zero();
        let one = field.one();
        let negative_one = field.neg(&one);
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let standard = BasedRightAlgebraModule::new(
            algebra.clone(),
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![one.clone(), zero.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone(), zero.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::from_rows(
                    vec![
                        vec![zero.clone(), zero.clone(), zero.clone()],
                        vec![zero.clone(), one.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone(), one.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
            ],
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![zero.clone(), zero.clone(), one.clone()],
                        vec![zero.clone(), zero.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone(), zero.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
            ],
        )
        .unwrap();
        let converted = right_algebra_module_to_path_algebra_mat_module(&standard).unwrap();

        let custom = BasedRightAlgebraModule::new(
            algebra.clone(),
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![one.clone(), zero.clone(), negative_one],
                        vec![zero.clone(), zero.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone(), zero.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::from_rows(
                    vec![
                        vec![zero.clone(), zero.clone(), one.clone()],
                        vec![zero.clone(), one.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone(), one.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
            ],
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![zero.clone(), zero.clone(), one.clone()],
                        vec![zero.clone(), zero.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone(), zero.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
            ],
        )
        .unwrap();
        let custom_converted = right_algebra_module_to_path_algebra_mat_module(&custom).unwrap();

        let zero_based = BasedRightAlgebraModule::new(
            algebra.clone(),
            vec![
                DenseMatrix::zero(0, 0, field.clone()),
                DenseMatrix::zero(0, 0, field.clone()),
            ],
            vec![DenseMatrix::zero(0, 0, field.clone())],
        )
        .unwrap();
        let zero_converted = right_algebra_module_to_path_algebra_mat_module(&zero_based).unwrap();

        let matrix_module = right_module_over_path_algebra(
            algebra,
            vec![2, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(
                    vec![
                        vec![one.clone(), one.clone()],
                        vec![zero.clone(), one.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let exposed = BasedRightAlgebraModule::from_path_algebra_module(&matrix_module).unwrap();
        let reconverted = right_algebra_module_to_path_algebra_mat_module(&exposed).unwrap();

        let loop_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap(),
        )
        .unwrap();
        let t = loop_parent.generators()[1].clone();
        let quotient = QuiverAlgebra::Quotient(
            QuotientPathAlgebra::by_relations_high_level(loop_parent, &[t.pow(2).unwrap()])
                .unwrap(),
        );
        let quotient_based = BasedRightAlgebraModule::new(
            quotient,
            vec![DenseMatrix::identity(2, field.clone())],
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![zero.clone(), one.clone()],
                        vec![zero.clone(), zero.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
            ],
        )
        .unwrap();
        let quotient_converted =
            right_algebra_module_to_path_algebra_mat_module(&quotient_based).unwrap();

        let arrowless_algebra = QuiverAlgebra::Path(
            PathAlgebra::new(field.clone(), Quiver::from_indexed(1, Vec::new()).unwrap()).unwrap(),
        );
        let arrowless_based = BasedRightAlgebraModule::new(
            arrowless_algebra,
            vec![DenseMatrix::identity(1, field)],
            Vec::new(),
        )
        .unwrap();
        let arrowless_converted =
            right_algebra_module_to_path_algebra_mat_module(&arrowless_based).unwrap();

        super::qpa_list([
            data(&converted, coefficient),
            data(&custom_converted, coefficient),
            data(&zero_converted, coefficient),
            data(&reconverted, coefficient),
            data(&quotient_converted, coefficient),
            data(&arrowless_converted, coefficient),
        ])
    }

    let rational = run(Q, ToString::to_string);
    let binary = run(Z2, |entry| {
        if Z2.is_zero(entry) {
            "0*Z(2)".to_owned()
        } else {
            "Z(2)^0".to_owned()
        }
    });
    assert_eq!(
        format!("QQ={rational}\nGF2={binary}\n"),
        include_str!("../../../../tests/fixtures/module_right_algebra_conversion.txt")
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn based_right_algebra_module_contract_failures_are_typed() {
    let algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    assert_eq!(
        BasedRightAlgebraModule::new(
            algebra.clone(),
            vec![DenseMatrix::identity(1, Q)],
            vec![DenseMatrix::identity(1, Q)],
        )
        .unwrap_err(),
        PathAlgebraModuleError::BasedModuleVertexActionCount {
            actual: 1,
            expected: 2,
        }
    );
    assert_eq!(
        BasedRightAlgebraModule::new(
            algebra.clone(),
            vec![DenseMatrix::identity(1, Q), DenseMatrix::identity(1, Q)],
            Vec::new(),
        )
        .unwrap_err(),
        PathAlgebraModuleError::BasedModuleArrowActionCount {
            actual: 0,
            expected: 1,
        }
    );
    assert!(matches!(
        BasedRightAlgebraModule::new(
            algebra.clone(),
            vec![DenseMatrix::identity(2, Q), DenseMatrix::identity(1, Q)],
            vec![DenseMatrix::identity(2, Q)],
        ),
        Err(PathAlgebraModuleError::BasedModuleActionShape { .. })
    ));

    let invalid_target = BasedRightAlgebraModule::new(
        algebra,
        vec![
            DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
            DenseMatrix::from_rows(vec![vec![q(0), q(0)], vec![q(0), q(1)]], Q).unwrap(),
        ],
        vec![DenseMatrix::identity(2, Q)],
    )
    .unwrap();
    assert!(matches!(
        right_algebra_module_to_path_algebra_mat_module(&invalid_target),
        Err(PathAlgebraModuleError::BasedModuleArrowImageOutsideTarget { .. })
    ));

    let loop_parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap(),
    )
    .unwrap();
    let t = loop_parent.generators()[1].clone();
    let quotient = QuiverAlgebra::Quotient(
        QuotientPathAlgebra::by_relations_high_level(loop_parent, &[t.pow(2).unwrap()]).unwrap(),
    );
    let relation_violating = BasedRightAlgebraModule::new(
        quotient,
        vec![DenseMatrix::identity(2, Q)],
        vec![DenseMatrix::identity(2, Q)],
    )
    .unwrap();
    assert_eq!(
        right_algebra_module_to_path_algebra_mat_module(&relation_violating).unwrap_err(),
        PathAlgebraModuleError::GeneratedBasedModuleConversionRejectedRelations
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn matrix_submodule_failures_match_qpa_dispatch_and_source_bugs() {
    let algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let module = right_module_over_path_algebra(
        algebra.clone(),
        vec![3, 3],
        &[NamedArrowMatrix::new("a", DenseMatrix::identity(3, Q))],
    )
    .unwrap()
    .unwrap();
    let empty = subalgebra_module(module.clone(), &[]).unwrap();
    assert!(empty.generators().is_empty());
    assert!(empty.canonical_basis_elements().unwrap().is_empty());
    let zero = subalgebra_module(module.clone(), &[module.zero_element()]).unwrap();
    assert!(zero.generators().is_empty());
    assert!(zero.canonical_basis_elements().unwrap().is_empty());
    assert_eq!(
        zero.as_module().unwrap_err(),
        PathAlgebraModuleError::ZeroSubmoduleAsModuleFailure
    );
    assert_eq!(
        natural_homomorphism_by_subalgebra_module(zero.parent(), &zero).unwrap_err(),
        PathAlgebraModuleError::ZeroNaturalSubmoduleFailure
    );

    let basis = module.basis_elements();
    let dependent = subalgebra_module(
        module,
        &[basis[0].add(&basis[4]).unwrap(), basis[0].clone()],
    )
    .unwrap();
    assert_eq!(
        dependent.canonical_basis_elements().unwrap_err(),
        PathAlgebraModuleError::DependentSubalgebraModuleGenerators
    );

    let small = right_module_over_path_algebra(
        algebra.clone(),
        vec![1, 1],
        &[NamedArrowMatrix::new("a", DenseMatrix::identity(1, Q))],
    )
    .unwrap()
    .unwrap();
    let small_basis = small.basis_elements();
    let target_only = subalgebra_module(small, &[small_basis[1].clone()]).unwrap();
    assert!(matches!(
        target_only.as_module(),
        Err(PathAlgebraModuleError::ZeroSpaceSubmoduleAsModuleFailure { .. })
    ));

    let medium = right_module_over_path_algebra(
        algebra.clone(),
        vec![2, 2],
        &[NamedArrowMatrix::new("a", DenseMatrix::identity(2, Q))],
    )
    .unwrap()
    .unwrap();
    let medium_basis = medium.basis_elements();
    let zero_space = subalgebra_module(
        medium.clone(),
        &[medium_basis[0].add(&medium_basis[3]).unwrap()],
    )
    .unwrap();
    assert!(matches!(
        natural_homomorphism_by_subalgebra_module(&medium, &zero_space),
        Err(PathAlgebraModuleError::ZeroSpaceNaturalSubmoduleFailure { .. })
    ));
    let foreign_parent = right_module_over_path_algebra(
        algebra,
        vec![2, 2],
        &[NamedArrowMatrix::new("a", DenseMatrix::identity(2, Q))],
    )
    .unwrap()
    .unwrap();
    assert_eq!(
        natural_homomorphism_by_subalgebra_module(&foreign_parent, &zero_space).unwrap_err(),
        PathAlgebraModuleError::SubalgebraModuleNotInEnteredParent
    );

    let loop_parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let quotient = QuiverAlgebra::Quotient(
        QuotientPathAlgebra::by_relations_high_level(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
    );
    let quotient_module = right_module_over_path_algebra(
        quotient,
        vec![2],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![q(0), q(1)], vec![q(0), q(0)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    assert_eq!(
        subalgebra_module(
            quotient_module.clone(),
            &[quotient_module.basis_elements()[1].clone()],
        )
        .unwrap_err(),
        PathAlgebraModuleError::QuotientSubalgebraModuleFailure
    );
}

#[test]
fn zero_module_reports_gap_scalar_zero_entries() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let module = right_module_over_path_algebra(algebra, vec![0, 0], &[])
        .unwrap()
        .unwrap();

    assert_eq!(module.arrow_matrices()[0], DenseMatrix::zero(1, 1, Q));
    assert_eq!(
        module.matrices_of_path_algebra_module(),
        vec![PathAlgebraModuleMatrix::ZeroScalar]
    );
}

#[test]
fn quotient_relations_return_fail_when_not_satisfied() {
    let quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let parent = PathAlgebra::new(Q, quiver).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let quotient =
        QuotientPathAlgebra::by_relations_high_level(parent, &[x.multiply(&x).unwrap()]).unwrap();
    let algebra = QuiverAlgebra::Quotient(quotient);
    let nilpotent = DenseMatrix::from_rows(vec![vec![q(0), q(1)], vec![q(0), q(0)]], Q).unwrap();
    let identity = DenseMatrix::identity(2, Q);

    assert!(
        right_module_over_path_algebra(
            algebra.clone(),
            vec![2],
            &[NamedArrowMatrix::new("x", nilpotent)],
        )
        .unwrap()
        .is_some()
    );
    assert!(
        right_module_over_path_algebra(algebra, vec![2], &[NamedArrowMatrix::new("x", identity)],)
            .unwrap()
            .is_none()
    );
}

#[test]
fn checked_constructor_preserves_qpa_input_errors() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());

    assert_eq!(
        right_module_over_path_algebra(algebra.clone(), vec![1], &[]).unwrap_err(),
        PathAlgebraModuleError::DimensionVectorLength {
            actual: 1,
            expected: 2,
        }
    );
    assert_eq!(
        right_module_over_path_algebra(
            algebra.clone(),
            vec![1, 1],
            &[NamedArrowMatrix::new("missing", DenseMatrix::zero(1, 1, Q))],
        )
        .unwrap_err(),
        PathAlgebraModuleError::UnknownArrowLabels {
            labels: vec![Box::from("missing")],
        }
    );
    assert!(matches!(
        right_module_over_path_algebra(
            algebra.clone(),
            vec![0, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![q(1)]], Q).unwrap(),
            )],
        ),
        Err(PathAlgebraModuleError::NonzeroMapAtZeroSpace { .. })
    ));
    assert!(matches!(
        right_module_over_path_algebra(
            algebra.clone(),
            vec![0, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_flat(0, 1, Vec::new(), Q).unwrap(),
            )],
        ),
        Err(PathAlgebraModuleError::ZeroSpaceMapMustBeOmitted { .. })
    ));
    assert!(matches!(
        right_module_over_path_algebra(
            algebra.clone(),
            vec![1, 2],
            &[NamedArrowMatrix::new("a", DenseMatrix::zero(1, 1, Q))],
        ),
        Err(PathAlgebraModuleError::WrongArrowMatrixShape { .. })
    ));

    let prime_three = Zp::new(3);
    let prime_five = Zp::new(5);
    let prime_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let prime_algebra = QuiverAlgebra::Path(PathAlgebra::new(prime_three, prime_quiver).unwrap());
    assert!(matches!(
        right_module_over_path_algebra(
            prime_algebra,
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::zero(1, 1, prime_five),
            )],
        ),
        Err(PathAlgebraModuleError::DifferentCoefficientField { .. })
    ));
}

#[test]
fn clones_are_identical_but_reconstruction_is_only_equal() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/final_declared_surface.schema-v1.json"
    ))
    .unwrap();
    assert_eq!(
        fixture["module_equality"],
        json!([true, true, true, false, false, false, false])
    );

    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let first = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
        .unwrap()
        .unwrap();
    let clone = first.clone();
    let second = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
        .unwrap()
        .unwrap();
    let different_dimension = right_module_over_path_algebra(algebra.clone(), vec![2, 1], &[])
        .unwrap()
        .unwrap();
    let different_matrix = right_module_over_path_algebra(
        algebra,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![q(1)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let other_quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
    let other_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, other_quiver).unwrap());
    let different_algebra = right_module_over_path_algebra(other_algebra, vec![1], &[])
        .unwrap()
        .unwrap();

    assert!(first.is_identical(&clone));
    assert!(!first.is_identical(&second));
    assert_eq!(first, second);
    assert_ne!(first, different_dimension);
    assert_ne!(first, different_matrix);
    assert_ne!(first, different_algebra);
}

#[test]
#[allow(clippy::too_many_lines)]
fn direct_sum_uses_qpa_block_order_component_maps_and_live_identity_cache() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let first = right_module_over_path_algebra(
        algebra.clone(),
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![q(1)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let second = right_module_over_path_algebra(
        algebra.clone(),
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![q(1)], vec![q(0)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let third = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
        .unwrap()
        .unwrap();

    assert!(
        direct_sum_of_qpa_modules::<numerica::domains::rational::RationalField>(&[])
            .unwrap()
            .is_none()
    );
    let direct_sum = direct_sum_of_qpa_modules(&[first.clone(), second.clone(), third.clone()])
        .unwrap()
        .unwrap();
    assert!(direct_sum.is_direct_sum_of_modules());
    assert_eq!(direct_sum.dimension_vector(), &[3, 3]);
    assert_eq!(
        direct_sum.arrow_matrices()[0],
        DenseMatrix::from_rows(
            vec![
                vec![q(1), q(0), q(0)],
                vec![q(0), q(1), q(0)],
                vec![q(0), q(0), q(0)],
            ],
            Q,
        )
        .unwrap()
    );
    let summands = direct_sum.direct_sum_summands().unwrap();
    assert!(summands[0].is_identical(&first));
    assert!(summands[1].is_identical(&second));
    assert!(summands[2].is_identical(&third));

    let projections = direct_sum_projections(&direct_sum).unwrap();
    let inclusions = direct_sum_inclusions(&direct_sum).unwrap();
    assert_eq!(
        projections[1].vertex_matrices(),
        &[
            DenseMatrix::from_rows(
                vec![vec![q(0), q(0)], vec![q(1), q(0)], vec![q(0), q(1)]],
                Q,
            )
            .unwrap(),
            DenseMatrix::from_rows(vec![vec![q(0)], vec![q(1)], vec![q(0)]], Q).unwrap(),
        ]
    );
    assert_eq!(
        inclusions[2].vertex_matrices(),
        &[
            DenseMatrix::zero(1, 3, Q),
            DenseMatrix::from_rows(vec![vec![q(0), q(0), q(1)]], Q).unwrap(),
        ]
    );
    for index in 0..3 {
        assert!(projections[index].is_surjective());
        assert!(inclusions[index].is_injective());
        assert_eq!(
            inclusions[index].compose(&projections[index]).unwrap(),
            identity_module_mapping(summands[index].clone()).unwrap()
        );
        for (other, projection) in projections.iter().enumerate() {
            if index != other {
                assert!(inclusions[index].compose(projection).unwrap().is_zero());
            }
        }
    }
    let identity_decomposition = projections
        .iter()
        .zip(&inclusions)
        .map(|(projection, inclusion)| projection.compose(inclusion).unwrap())
        .reduce(|left, right| left.add(&right).unwrap())
        .unwrap();
    assert_eq!(
        identity_decomposition,
        identity_module_mapping(direct_sum.clone()).unwrap()
    );

    let retained_projection = projections[0].clone();
    drop(projections);
    let regenerated = direct_sum_projections(&direct_sum).unwrap();
    assert!(retained_projection.is_identical(&regenerated[0]));
    let repeated = direct_sum_projections(&direct_sum).unwrap();
    assert!(
        regenerated
            .iter()
            .zip(repeated)
            .all(|(left, right)| left.is_identical(&right))
    );

    let singleton = direct_sum_of_qpa_modules(std::slice::from_ref(&first))
        .unwrap()
        .unwrap();
    assert!(!singleton.is_identical(&first));
    assert_eq!(singleton, first);
    assert!(direct_sum_projections(&first).is_err());
}

#[test]
fn characteristic_two_direct_sum_uses_the_same_qpa_block_algorithm() {
    let e = |value: i64| Z2.nth(Integer::from(value));
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
    let first = right_module_over_path_algebra(
        algebra.clone(),
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![e(1)]], Z2).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let second = right_module_over_path_algebra(algebra, vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let direct_sum = direct_sum_of_qpa_modules(&[first, second])
        .unwrap()
        .unwrap();

    assert_eq!(direct_sum.dimension_vector(), &[2, 1]);
    assert_eq!(
        direct_sum.arrow_matrices()[0],
        DenseMatrix::from_rows(vec![vec![e(1)], vec![e(0)]], Z2).unwrap()
    );
    let projections = direct_sum_projections(&direct_sum).unwrap();
    let inclusions = direct_sum_inclusions(&direct_sum).unwrap();
    assert!(
        inclusions[0]
            .compose(&projections[0])
            .unwrap()
            .is_isomorphism()
    );
    assert!(
        inclusions[1]
            .compose(&projections[1])
            .unwrap()
            .is_isomorphism()
    );
}

#[test]
fn pinned_qpa_modulehom_categorical_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_categorical.schema-v1.json"
    ))
    .unwrap();
    let rational =
        modulehom_categorical_case(Q, |value| Q.nth(Integer::from(value)), ToString::to_string);
    let binary = modulehom_categorical_case(
        Z2,
        |value| Z2.nth(Integer::from(value)),
        |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        },
    );

    assert_eq!(rational, expected["rational"]);
    assert_eq!(binary, expected["binary"]);
}

#[test]
fn pinned_qpa_modulehom_split_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_split.schema-v1.json"
    ))
    .unwrap();
    let rational = modulehom_split_case(Q, ToString::to_string);
    let binary = modulehom_split_case(Z2, |entry| {
        if Z2.is_zero(entry) {
            "0*Z(2)".to_owned()
        } else {
            "Z(2)^0".to_owned()
        }
    });

    assert_eq!(rational, expected["rational"]);
    assert_eq!(binary, expected["binary"]);
}

#[test]
fn pinned_qpa_modulehom_map_utils_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_map_utils.schema-v1.json"
    ))
    .unwrap();
    let rational =
        modulehom_map_utils_case(Q, |value| Q.nth(Integer::from(value)), ToString::to_string);
    let binary = modulehom_map_utils_case(
        Z2,
        |value| Z2.nth(Integer::from(value)),
        |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        },
    );

    assert_eq!(rational, expected["rational"]);
    assert_eq!(binary, expected["binary"]);
}

#[test]
fn pinned_qpa_modulehom_restriction_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_restriction.schema-v1.json"
    ))
    .unwrap();
    let quinary = Zp::new(5);
    let quinary_coefficients = quinary.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehom_restriction",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehom_restriction_case(Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": modulehom_restriction_case(quinary, |entry| json!(
            quinary_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "rational": modulehom_restriction_case(Q, |entry| json!([
            entry.numerator_ref().to_string().parse::<i64>().unwrap(),
            entry.denominator_ref().to_string().parse::<i64>().unwrap(),
        ])),
    });

    assert_eq!(actual, expected);
}

#[test]
fn restriction_preserves_qpa_range_validation_boundaries() {
    let field_algebra = QuiverAlgebra::Path(
        PathAlgebra::new(Q, Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap()).unwrap(),
    );
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let inclusion = QuiverAlgebraHomomorphism::by_generators_with_one_images(
        field_algebra.clone(),
        a2.clone(),
        vec![a2.one().unwrap()],
    )
    .unwrap();
    let wrong = right_module_over_path_algebra(field_algebra, vec![1], &[])
        .unwrap()
        .unwrap();
    let wrong_identity = identity_module_mapping(wrong.clone()).unwrap();

    assert_eq!(
        restriction_via_algebra_homomorphism(&inclusion, &wrong).unwrap_err(),
        super::RestrictionViaAlgebraHomomorphismError::ModuleNotOverRangeAlgebra
    );
    assert_eq!(
        restriction_via_algebra_homomorphism_map(&inclusion, &wrong_identity).unwrap_err(),
        super::RestrictionViaAlgebraHomomorphismError::HomomorphismNotOverRangeAlgebra
    );
}

#[test]
fn pinned_qpa_modulehom_isomorphism_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_isomorphism.schema-v1.json"
    ))
    .unwrap();
    let rational =
        modulehom_isomorphism_case(Q, |value| Q.nth(Integer::from(value)), ToString::to_string);
    let binary = modulehom_isomorphism_case(
        Z2,
        |value| Z2.nth(Integer::from(value)),
        |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        },
    );

    assert_eq!(rational, expected["rational"]);
    assert_eq!(binary, expected["binary"]);
}

#[test]
fn map_utility_failures_preserve_qpa_validation_boundaries() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let module = right_module_over_path_algebra(
        algebra,
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![q(1)], vec![q(0)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let zero = zero_module_mapping(module.clone(), module.clone()).unwrap();
    assert_eq!(
        inverse_of_isomorphism(&zero).unwrap_err(),
        PathAlgebraModuleError::EnteredMapNotIsomorphism
    );

    let basis = module.basis_elements();
    let mut wrong_support = basis.clone();
    wrong_support[0] = basis[2].clone();
    assert_eq!(
        homomorphism_from_images(&module, &module, &wrong_support).unwrap_err(),
        PathAlgebraModuleError::WrongBasisImageSupport { position: 1 }
    );
    assert_eq!(
        homomorphism_from_images(&module, &module, &basis[..2]).unwrap_err(),
        PathAlgebraModuleError::BasisImageCount {
            actual: 2,
            expected: 3,
        }
    );
    assert_eq!(
        multiply_lists_of_maps::<numerica::domains::rational::Q>(&[], &[], &[]).unwrap_err(),
        PathAlgebraModuleError::EmptyMapProjectionList
    );
}

#[test]
fn pinned_qpa_modulehom_minimal_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_minimal.schema-v1.json"
    ))
    .unwrap();
    let rational =
        modulehom_minimal_case(Q, |value| Q.nth(Integer::from(value)), ToString::to_string);
    let binary = modulehom_minimal_case(
        Z2,
        |value| Z2.nth(Integer::from(value)),
        |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        },
    );

    assert_eq!(rational, expected["rational"]);
    assert_eq!(binary, expected["binary"]);
}

#[test]
fn pinned_qpa_modulehom_common_summand_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_common_summand.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehom_common_summand",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "rational": modulehom_common_summand_case(Q),
        "binary": modulehom_common_summand_case(Z2),
    });
    assert_eq!(actual, expected);
}

#[test]
fn left_minimal_zero_endomorphism_preserves_qpa_f_times_false_failure() {
    let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let one = right_module_over_path_algebra(algebra, vec![1], &[])
        .unwrap()
        .unwrap();
    let zero = zero_module_mapping(one.clone(), one).unwrap();

    assert_eq!(
        more_left_minimal_version(&zero).unwrap_err(),
        PathAlgebraModuleError::InternalMinimalKernelNotSplit
    );
    assert_eq!(
        left_minimal_version(&zero).unwrap_err(),
        PathAlgebraModuleError::InternalMinimalKernelNotSplit
    );
}

#[test]
fn categorical_validation_preserves_qpa_fail_and_error_precedence() {
    let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let one = right_module_over_path_algebra(algebra.clone(), vec![1], &[])
        .unwrap()
        .unwrap();
    let two = right_module_over_path_algebra(algebra.clone(), vec![2], &[])
        .unwrap()
        .unwrap();
    let three = right_module_over_path_algebra(algebra.clone(), vec![3], &[])
        .unwrap()
        .unwrap();
    let inclusion = right_module_hom_over_algebra(
        one.clone(),
        two.clone(),
        vec![DenseMatrix::from_rows(vec![vec![q(1), q(0)]], Q).unwrap()],
    )
    .unwrap();
    let zero = zero_module_mapping(one.clone(), two.clone()).unwrap();
    let different_range = right_module_hom_over_algebra(
        one.clone(),
        three,
        vec![DenseMatrix::from_rows(vec![vec![q(1), q(0), q(0)]], Q).unwrap()],
    )
    .unwrap();

    assert!(
        push_out(&inclusion, &identity_module_mapping(two.clone()).unwrap())
            .unwrap()
            .is_none()
    );
    assert!(pull_back(&inclusion, &different_range).unwrap().is_none());
    assert_eq!(
        intersection_of_submodules(&zero, &inclusion).unwrap_err(),
        PathAlgebraModuleError::FirstSubmoduleMapNotInjective
    );
    assert_eq!(
        sum_of_submodules(&inclusion, &zero).unwrap_err(),
        PathAlgebraModuleError::SecondSubmoduleMapNotInjective
    );
    assert_eq!(
        intersection_of_submodules(&inclusion, &different_range).unwrap_err(),
        PathAlgebraModuleError::SubmoduleRangesDiffer
    );
    assert_eq!(
        intersection_of_submodules_list::<numerica::domains::rational::RationalField>(&[])
            .unwrap_err(),
        PathAlgebraModuleError::EmptySubmoduleList
    );
    assert_eq!(
        sum_of_submodules_list::<numerica::domains::rational::RationalField>(&[]).unwrap_err(),
        PathAlgebraModuleError::EmptySubmoduleList
    );
    assert_eq!(
        intersection_of_submodules_list(&[
            inclusion.clone(),
            zero.clone(),
            different_range.clone(),
        ])
        .unwrap_err(),
        PathAlgebraModuleError::SomeSubmoduleMapNotInjective
    );
    assert_eq!(
        sum_of_submodules_list(&[inclusion.clone(), zero, different_range]).unwrap_err(),
        PathAlgebraModuleError::SubmoduleRangesDiffer
    );

    let foreign_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let foreign_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, foreign_quiver).unwrap());
    let foreign = right_module_over_path_algebra(
        foreign_algebra,
        vec![1],
        &[NamedArrowMatrix::new("x", DenseMatrix::zero(1, 1, Q))],
    )
    .unwrap()
    .unwrap();
    let foreign_identity = identity_module_mapping(foreign.clone()).unwrap();
    assert_eq!(
        direct_sum_of_qpa_modules(&[one.clone(), foreign.clone()]).unwrap_err(),
        PathAlgebraModuleError::DirectSumDifferentAlgebras
    );
    assert_eq!(
        intersection_of_submodules_list(&[inclusion, foreign_identity]).unwrap_err(),
        PathAlgebraModuleError::SubmoduleListDifferentAlgebras
    );
    assert_eq!(
        reject_of_module(&one, &foreign).unwrap_err(),
        PathAlgebraModuleError::RejectModulesOverDifferentAlgebras
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn foundational_module_map_factorizations_follow_qpa_row_bases() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let action = DenseMatrix::from_rows(vec![vec![q(1), q(1)], vec![q(0), q(0)]], Q).unwrap();
    let module =
        right_module_over_path_algebra(algebra, vec![2, 2], &[NamedArrowMatrix::new("a", action)])
            .unwrap()
            .unwrap();
    let map = right_module_hom_over_algebra(
        module.clone(),
        module.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
            DenseMatrix::identity(2, Q),
        ],
    )
    .unwrap();

    assert!(!map.is_injective());
    assert!(!map.is_surjective());
    assert!(!map.is_isomorphism());
    assert!(!map.is_zero());

    let kernel = kernel_inclusion(&map).unwrap();
    assert_eq!(kernel.source().dimension_vector(), &[1, 0]);
    assert_eq!(
        kernel.vertex_matrices(),
        &[
            DenseMatrix::from_rows(vec![vec![q(0), q(1)]], Q).unwrap(),
            DenseMatrix::zero(1, 2, Q),
        ]
    );
    assert!(kernel.is_injective());
    assert!(kernel.kernel_of_what().unwrap().is_identical(&map));
    assert!(kernel_inclusion(&map).unwrap().is_identical(&kernel));
    assert!(kernel.compose(&map).unwrap().is_zero());

    let (projection, inclusion) = image_projection_inclusion(&map).unwrap();
    assert_eq!(projection.range().dimension_vector(), &[1, 2]);
    assert_eq!(
        projection.vertex_matrices(),
        &[
            DenseMatrix::from_rows(vec![vec![q(1)], vec![q(0)]], Q).unwrap(),
            DenseMatrix::identity(2, Q),
        ]
    );
    assert_eq!(
        inclusion.vertex_matrices(),
        &[
            DenseMatrix::from_rows(vec![vec![q(1), q(0)]], Q).unwrap(),
            DenseMatrix::identity(2, Q),
        ]
    );
    assert!(projection.is_surjective());
    assert!(inclusion.is_injective());
    assert!(projection.image_of_what().unwrap().is_identical(&map));
    assert!(inclusion.image_of_what().unwrap().is_identical(&map));
    assert!(
        image_projection_inclusion(&map)
            .unwrap()
            .0
            .is_identical(&projection)
    );
    assert_eq!(projection.compose(&inclusion).unwrap(), map);

    let cokernel = cokernel_projection(&map).unwrap();
    assert_eq!(cokernel.range().dimension_vector(), &[1, 0]);
    assert_eq!(
        cokernel.vertex_matrices(),
        &[
            DenseMatrix::from_rows(vec![vec![q(0)], vec![q(1)]], Q).unwrap(),
            DenseMatrix::zero(2, 1, Q),
        ]
    );
    assert!(cokernel.is_surjective());
    assert!(cokernel.cokernel_of_what().unwrap().is_identical(&map));
    assert!(cokernel_projection(&map).unwrap().is_identical(&cokernel));
    assert!(map.compose(&cokernel).unwrap().is_zero());

    let zero = zero_module_mapping(module.clone(), module.clone()).unwrap();
    let identity = identity_module_mapping(module.clone()).unwrap();
    assert_eq!(map.add(&zero).unwrap(), map);
    assert!(map.add(&map.negate().unwrap()).unwrap().is_zero());
    assert!(identity.is_isomorphism());
    assert_eq!(identity.compose(&map).unwrap(), map);
    assert_eq!(map.compose(&identity).unwrap(), map);
    assert!(zero.zero().unwrap().is_zero());

    let top_map = top_of_module_homomorphism(&map).unwrap();
    assert_eq!(top_map.source().dimension_vector(), &[2, 1]);
    assert_eq!(top_map.range().dimension_vector(), &[2, 1]);
    assert_eq!(
        top_map.vertex_matrices(),
        &[
            DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
            DenseMatrix::identity(1, Q),
        ]
    );

    let hom_basis = hom_over_algebra(&module, &module).unwrap().unwrap();
    assert_eq!(hom_basis.len(), 5);
    let z = DenseMatrix::zero(2, 2, Q);
    let expected = vec![
        vec![
            DenseMatrix::from_rows(vec![vec![q(0), q(1)], vec![q(0), q(0)]], Q).unwrap(),
            z.clone(),
        ],
        vec![
            DenseMatrix::from_rows(vec![vec![q(0), q(0)], vec![q(0), q(1)]], Q).unwrap(),
            z.clone(),
        ],
        vec![
            DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
            DenseMatrix::from_rows(vec![vec![q(1), q(1)], vec![q(0), q(0)]], Q).unwrap(),
        ],
        vec![
            z,
            DenseMatrix::from_rows(vec![vec![q(-1), q(0)], vec![q(1), q(0)]], Q).unwrap(),
        ],
        map.vertex_matrices().to_vec(),
    ];
    assert_eq!(
        hom_basis
            .iter()
            .map(|homomorphism| homomorphism.vertex_matrices().to_vec())
            .collect::<Vec<_>>(),
        expected
    );
    let hom_with_basis = hom_over_algebra_with_basis_function(&module, &module)
        .unwrap()
        .unwrap();
    assert_eq!(hom_with_basis.basis(), hom_basis);
    assert_eq!(
        hom_with_basis.coefficients(&map).unwrap().unwrap(),
        vec![q(0), q(0), q(0), q(0), q(1)]
    );

    let endomorphism_algebra = end_over_algebra(&module).unwrap();
    assert_eq!(endomorphism_algebra.dimension(), 5);
    assert!(
        endomorphism_algebra.is_identical(&end_over_algebra(&module).unwrap()),
        "QPA declares EndOverAlgebra as a cached attribute"
    );
    assert_eq!(
        endomorphism_algebra.basis_matrices(),
        &[
            DenseMatrix::from_rows(
                vec![
                    vec![q(0), q(1), q(0), q(0)],
                    vec![q(0), q(0), q(0), q(0)],
                    vec![q(0), q(0), q(0), q(0)],
                    vec![q(0), q(0), q(0), q(0)],
                ],
                Q,
            )
            .unwrap(),
            DenseMatrix::from_rows(
                vec![
                    vec![q(0), q(0), q(0), q(0)],
                    vec![q(0), q(1), q(0), q(0)],
                    vec![q(0), q(0), q(0), q(0)],
                    vec![q(0), q(0), q(0), q(0)],
                ],
                Q,
            )
            .unwrap(),
            DenseMatrix::from_rows(
                vec![
                    vec![q(1), q(0), q(0), q(0)],
                    vec![q(0), q(0), q(0), q(0)],
                    vec![q(0), q(0), q(1), q(1)],
                    vec![q(0), q(0), q(0), q(0)],
                ],
                Q,
            )
            .unwrap(),
            DenseMatrix::from_rows(
                vec![
                    vec![q(0), q(0), q(0), q(0)],
                    vec![q(0), q(0), q(0), q(0)],
                    vec![q(0), q(0), q(1), q(0)],
                    vec![q(0), q(0), q(-1), q(0)],
                ],
                Q,
            )
            .unwrap(),
            DenseMatrix::from_rows(
                vec![
                    vec![q(0), q(0), q(0), q(0)],
                    vec![q(0), q(0), q(0), q(0)],
                    vec![q(0), q(0), q(0), q(1)],
                    vec![q(0), q(0), q(0), q(-1)],
                ],
                Q,
            )
            .unwrap(),
        ]
    );
    assert_eq!(*endomorphism_algebra.one(), DenseMatrix::identity(4, Q));
    assert_eq!(
        endomorphism_algebra
            .coefficients(endomorphism_algebra.one())
            .unwrap()
            .unwrap(),
        vec![q(0), q(1), q(1), q(0), q(-1)]
    );
    assert_eq!(
        endomorphism_algebra
            .multiply(
                &endomorphism_algebra.basis_matrices()[0],
                endomorphism_algebra.one(),
            )
            .unwrap(),
        endomorphism_algebra.basis_matrices()[0]
    );
    let cached_data = Arc::as_ptr(&endomorphism_algebra.0);
    drop(endomorphism_algebra);
    assert_eq!(
        Arc::as_ptr(&end_over_algebra(&module).unwrap().0),
        cached_data,
        "the GAP attribute cache survives after the caller drops its handle"
    );

    let basis = module.basis_elements();
    let nonzero_images = map.images_set(&basis).unwrap();
    assert_eq!(nonzero_images.len(), 3);
    assert!(nonzero_images.iter().all(|image| !image.is_zero()));
}

#[test]
fn characteristic_two_module_map_factorizations_match_the_same_qpa_fixture() {
    let e = |value: i64| Z2.nth(Integer::from(value));
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
    let module = right_module_over_path_algebra(
        algebra,
        vec![2, 2],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![e(1), e(1)], vec![e(0), e(0)]], Z2).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let map = right_module_hom_over_algebra(
        module.clone(),
        module,
        vec![
            DenseMatrix::from_rows(vec![vec![e(1), e(0)], vec![e(0), e(0)]], Z2).unwrap(),
            DenseMatrix::identity(2, Z2),
        ],
    )
    .unwrap();

    let kernel = kernel_inclusion(&map).unwrap();
    assert_eq!(kernel.source().dimension_vector(), &[1, 0]);
    assert_eq!(
        kernel.vertex_matrices()[0],
        DenseMatrix::from_rows(vec![vec![e(0), e(1)]], Z2).unwrap()
    );
    let (projection, inclusion) = image_projection_inclusion(&map).unwrap();
    assert_eq!(projection.range().dimension_vector(), &[1, 2]);
    assert_eq!(projection.compose(&inclusion).unwrap(), map);
    let cokernel = cokernel_projection(&map).unwrap();
    assert_eq!(cokernel.range().dimension_vector(), &[1, 0]);
    assert_eq!(
        cokernel.vertex_matrices()[0],
        DenseMatrix::from_rows(vec![vec![e(0)], vec![e(1)]], Z2).unwrap()
    );
    let hom_basis = hom_over_algebra(map.source(), map.range())
        .unwrap()
        .unwrap();
    assert_eq!(hom_basis.len(), 5);
    assert_eq!(
        hom_basis[3].vertex_matrices()[1],
        DenseMatrix::from_rows(vec![vec![e(1), e(0)], vec![e(1), e(0)]], Z2,).unwrap()
    );
    assert_eq!(
        hom_over_algebra_with_basis_function(map.source(), map.range())
            .unwrap()
            .unwrap()
            .coefficients(&map)
            .unwrap()
            .unwrap(),
        vec![e(0), e(0), e(0), e(0), e(1)]
    );
    let endomorphism_algebra = end_over_algebra(map.source()).unwrap();
    assert_eq!(endomorphism_algebra.dimension(), 5);
    assert_eq!(
        endomorphism_algebra.basis_matrices()[3],
        DenseMatrix::from_rows(
            vec![
                vec![e(0), e(0), e(0), e(0)],
                vec![e(0), e(0), e(0), e(0)],
                vec![e(0), e(0), e(1), e(0)],
                vec![e(0), e(0), e(1), e(0)],
            ],
            Z2,
        )
        .unwrap()
    );
    assert_eq!(
        endomorphism_algebra
            .coefficients(endomorphism_algebra.one())
            .unwrap()
            .unwrap(),
        vec![e(0), e(1), e(1), e(0), e(1)]
    );
}

#[test]
fn unconstrained_and_disjoint_hom_branches_match_qpa() {
    let free_quiver = Quiver::from_indexed(1, vec![]).unwrap();
    let free_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, free_quiver).unwrap());
    let free_source = right_module_over_path_algebra(free_algebra.clone(), vec![1], &[])
        .unwrap()
        .unwrap();
    let free_range = right_module_over_path_algebra(free_algebra, vec![2], &[])
        .unwrap()
        .unwrap();
    let free_basis = hom_over_algebra(&free_source, &free_range)
        .unwrap()
        .unwrap();
    assert_eq!(free_basis.len(), 2);
    assert_eq!(
        free_basis
            .iter()
            .map(|homomorphism| homomorphism.vertex_matrices()[0].clone())
            .collect::<Vec<_>>(),
        vec![
            DenseMatrix::from_rows(vec![vec![q(1), q(0)]], Q).unwrap(),
            DenseMatrix::from_rows(vec![vec![q(0), q(1)]], Q).unwrap(),
        ]
    );

    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let disjoint_source = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let disjoint_range = right_module_over_path_algebra(algebra, vec![0, 1], &[])
        .unwrap()
        .unwrap();
    assert!(
        hom_over_algebra(&disjoint_source, &disjoint_range)
            .unwrap()
            .unwrap()
            .is_empty()
    );
    let zero = zero_module_mapping(disjoint_source.clone(), disjoint_range.clone()).unwrap();
    assert_eq!(
        hom_over_algebra_with_basis_function(&disjoint_source, &disjoint_range)
            .unwrap()
            .unwrap()
            .coefficients(&zero)
            .unwrap()
            .unwrap(),
        vec![q(0)]
    );

    let zero_module = right_module_over_path_algebra(
        disjoint_source.right_acting_algebra().clone(),
        vec![0, 0],
        &[],
    )
    .unwrap()
    .unwrap();
    assert_eq!(
        end_over_algebra(&zero_module).unwrap_err(),
        PathAlgebraModuleError::EmptyEndomorphismDirectSum
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn rational_minimal_generators_socle_and_trace_match_qpa() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let module = right_module_over_path_algebra(
        algebra.clone(),
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![q(1)], vec![q(0)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let minimal = minimal_generating_set_of_module(module.clone()).unwrap();
    assert_eq!(minimal.len(), 2);
    assert_eq!(
        minimal[0].external_representation(),
        &[vec![q(1), q(0)], vec![q(0)]]
    );
    assert_eq!(
        minimal[1].external_representation(),
        &[vec![q(0), q(1)], vec![q(0)]]
    );
    let top_projection = top_of_module_projection(module.clone()).unwrap();
    assert!(top_projection.is_identical(&top_of_module_projection(module.clone()).unwrap()));
    let top = top_of_module(module.clone()).unwrap();
    assert!(top.is_identical(&top_of_module(module.clone()).unwrap()));
    assert!(top.is_identical(top_projection.range()));

    let socle = socle_of_module_inclusion(module.clone()).unwrap();
    assert_eq!(socle.source().dimension_vector(), &[1, 1]);
    assert_eq!(socle.range(), &module);
    assert_eq!(
        socle.vertex_matrices(),
        &[
            DenseMatrix::from_rows(vec![vec![q(0), q(1)]], Q).unwrap(),
            DenseMatrix::identity(1, Q),
        ]
    );
    assert_eq!(
        socle_of_module(module.clone()).unwrap().dimension_vector(),
        &[1, 1]
    );

    let zero_module = right_module_over_path_algebra(algebra.clone(), vec![0, 0], &[])
        .unwrap()
        .unwrap();
    let zero_socle = socle_of_module_inclusion(zero_module.clone()).unwrap();
    assert!(zero_socle.source().is_identical(&zero_module));
    assert!(zero_socle.range().is_identical(&zero_module));
    assert!(
        minimal_generating_set_of_module(zero_module.clone())
            .unwrap()
            .is_empty()
    );

    let target = right_module_over_path_algebra(
        algebra.clone(),
        vec![1, 2],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![q(1), q(1)]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let simple_two = right_module_over_path_algebra(algebra, vec![0, 1], &[])
        .unwrap()
        .unwrap();

    let zero_trace = trace_of_module(&simple_one, &target).unwrap();
    assert_eq!(zero_trace.source().dimension_vector(), &[0, 0]);
    assert_eq!(zero_trace.range(), &target);
    assert_eq!(
        zero_trace.vertex_matrices(),
        &[DenseMatrix::zero(1, 1, Q), DenseMatrix::zero(1, 2, Q),]
    );
    assert_eq!(trace_of_module(&zero_module, &target).unwrap(), zero_trace);
    let zero_target_trace = trace_of_module(&simple_two, &zero_module).unwrap();
    assert_eq!(zero_target_trace.source().dimension_vector(), &[0, 0]);
    assert!(zero_target_trace.range().is_identical(&zero_module));

    let trace = trace_of_module(&simple_two, &target).unwrap();
    assert_eq!(trace.source().dimension_vector(), &[0, 2]);
    assert_eq!(
        trace.vertex_matrices(),
        &[DenseMatrix::zero(1, 1, Q), DenseMatrix::identity(2, Q)]
    );
    assert_eq!(
        right_fac_approximation(&simple_two, &target).unwrap(),
        trace
    );
    assert_eq!(
        right_fac_approximation_via_all_basis_images(&simple_two, &target).unwrap(),
        trace
    );

    let foreign =
        QuiverAlgebra::Path(PathAlgebra::new(Q, Quiver::from_indexed(1, vec![]).unwrap()).unwrap());
    let foreign_module = right_module_over_path_algebra(foreign, vec![1], &[])
        .unwrap()
        .unwrap();
    assert_eq!(
        trace_of_module(&foreign_module, &target).unwrap_err(),
        PathAlgebraModuleError::DifferentHomomorphismAlgebras
    );
}

#[test]
fn binary_socle_and_trace_use_the_same_qpa_nullspaces() {
    let e = |value: i64| Z2.nth(Integer::from(value));
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
    let module = right_module_over_path_algebra(
        algebra.clone(),
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![e(1)], vec![e(0)]], Z2).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let socle = socle_of_module_inclusion(module).unwrap();
    assert_eq!(socle.source().dimension_vector(), &[1, 1]);
    assert_eq!(
        socle.vertex_matrices()[0],
        DenseMatrix::from_rows(vec![vec![e(0), e(1)]], Z2).unwrap()
    );

    let target = right_module_over_path_algebra(
        algebra.clone(),
        vec![1, 2],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![e(1), e(1)]], Z2).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let simple_two = right_module_over_path_algebra(algebra, vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let trace = trace_of_module(&simple_two, &target).unwrap();
    assert_eq!(trace.source().dimension_vector(), &[0, 2]);
    assert_eq!(trace.vertex_matrices()[1], DenseMatrix::identity(2, Z2));
}

#[test]
fn rational_radical_and_top_preserve_gap_quotient_lifts() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let action = DenseMatrix::from_rows(vec![vec![q(1), q(1)], vec![q(0), q(0)]], Q).unwrap();
    let module =
        right_module_over_path_algebra(algebra, vec![2, 2], &[NamedArrowMatrix::new("a", action)])
            .unwrap()
            .unwrap();

    let radical = radical_of_module_inclusion(module.clone()).unwrap();
    assert_eq!(radical.source().dimension_vector(), &[0, 1]);
    assert_eq!(radical.vertex_matrices()[0], DenseMatrix::zero(1, 2, Q));
    assert_eq!(
        radical.vertex_matrices()[1],
        DenseMatrix::from_rows(vec![vec![q(1), q(1)]], Q).unwrap()
    );

    let top = top_of_module_projection(module).unwrap();
    assert_eq!(top.range().dimension_vector(), &[2, 1]);
    assert_eq!(top.range().arrow_matrices()[0], DenseMatrix::zero(2, 1, Q));
    assert_eq!(top.vertex_matrices()[0], DenseMatrix::identity(2, Q));
    assert_eq!(
        top.vertex_matrices()[1],
        DenseMatrix::from_rows(vec![vec![q(1)], vec![q(-1)]], Q).unwrap()
    );

    let lifts = top
        .range()
        .basis_elements()
        .iter()
        .map(|basis| {
            top.preimage_representative(basis)
                .unwrap()
                .unwrap()
                .external_representation()
                .to_vec()
        })
        .collect::<Vec<_>>();
    assert_eq!(
        lifts,
        vec![
            vec![vec![q(1), q(0)], vec![q(0), q(0)]],
            vec![vec![q(0), q(1)], vec![q(0), q(0)]],
            vec![vec![q(0), q(0)], vec![q(1), q(0)]],
        ]
    );
}

#[test]
fn characteristic_two_top_uses_the_same_qpa_source_method() {
    let e = |value: i64| Z2.nth(Integer::from(value));
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
    let action = DenseMatrix::from_rows(vec![vec![e(1), e(1)], vec![e(0), e(0)]], Z2).unwrap();
    let module =
        right_module_over_path_algebra(algebra, vec![2, 2], &[NamedArrowMatrix::new("a", action)])
            .unwrap()
            .unwrap();

    let top = top_of_module_projection(module).unwrap();
    assert_eq!(top.range().dimension_vector(), &[2, 1]);
    assert_eq!(
        top.vertex_matrices()[1],
        DenseMatrix::from_rows(vec![vec![e(1)], vec![e(1)]], Z2).unwrap()
    );
}

#[test]
fn zero_top_is_qpa_zero_endomorphism_of_the_identical_module() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let module = right_module_over_path_algebra(algebra, vec![0, 0], &[])
        .unwrap()
        .unwrap();

    let top = top_of_module_projection(module.clone()).unwrap();
    assert!(top.source().is_identical(&module));
    assert!(top.range().is_identical(&module));
    assert_eq!(
        top.vertex_matrices(),
        &[DenseMatrix::zero(1, 1, Q), DenseMatrix::zero(1, 1, Q)]
    );
}

#[test]
fn endomorphism_ideal_closure_is_two_sided_before_quotienting() {
    let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let module = right_module_over_path_algebra(algebra, vec![2], &[])
        .unwrap()
        .unwrap();
    let endomorphisms = end_over_algebra(&module).unwrap();
    let matrix_unit = DenseMatrix::from_rows(vec![vec![q(0), q(1)], vec![q(0), q(0)]], Q).unwrap();
    let natural = natural_endomorphism_quotient_by_ideal(&endomorphisms, &[matrix_unit]).unwrap();

    assert_eq!(endomorphisms.dimension(), 4);
    assert_eq!(natural.kernel_basis_matrices().len(), 4);
    assert_eq!(natural.range().dimension(), 0);
    assert!(natural.range().one().unwrap().is_none());
    assert!(
        endomorphisms
            .basis_matrices()
            .iter()
            .all(|matrix| natural.image(matrix).unwrap().is_zero())
    );
}

#[test]
fn lift_idempotents_preserves_qpa_missing_method_boundary() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/final_declared_surface.schema-v1.json"
    ))
    .unwrap();
    assert_eq!(fixture["uninstalled_success"][0], false);

    let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
    let module = right_module_over_path_algebra(algebra, vec![1], &[])
        .unwrap()
        .unwrap();
    assert!(matches!(
        lift_idempotents(&module),
        Err(PathAlgebraModuleError::UninstalledLiftIdempotents)
    ));
}
