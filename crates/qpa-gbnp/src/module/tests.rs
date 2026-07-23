use numerica::domains::{
    Field, Ring,
    finite_field::{FiniteFieldCore, Z2, Zp},
    integer::Integer,
    rational::{Q, Rational, RationalField},
};
use qpa_core::path_module::right_fac_approximation_via_all_basis_images;
use qpa_core::{
    AlgebraGlobalDimension, ArrowSpec, CachedAttribute, DenseMatrix, GaloisField,
    ModuleDimensionSeries, NamedArrowMatrix, NamedArrowSpec, PathAlgebra, PathAlgebraModule,
    PathAlgebraModuleElement, PathAlgebraModuleHomomorphism, PathAlgebraModuleHomomorphismPair,
    PathAlgebraModuleMatrix, QpaCoefficientField, Quiver, QuiverAlgebra, QuiverAlgebraElement,
    StructureConstantsAlgebra, VertexIndex, dimension_vector_partial_order,
    direct_sum_of_qpa_modules, end_over_algebra, from_end_m_to_hom_mm, from_hom_mm_to_end_m,
    identity_module_mapping, include_in_path_algebra, is_indecomposable_module,
    is_semisimple_module, is_simple_qpa_module, kernel_inclusion, lifting_inclusion_morphisms,
    minimal_generating_set_of_module, radical_series, right_fac_approximation,
    right_module_hom_over_algebra, right_module_over_path_algebra, socle_of_module,
    socle_of_module_inclusion, top_of_module, top_of_module_projection, trace_of_module,
    underlying_linear_map, zero_module, zero_module_mapping,
};
use rand::{RngCore, SeedableRng, rngs::StdRng};
use serde_json::{Value, json};

use super::{
    AddMCoresolution, AddMResolution, AllComplementsResult, AlmostSplitSide, AuslanderReitenError,
    CotiltingModuleResult, DominantDimensionResult, EnvelopingModuleError,
    ExtAlgebraGeneratorsResult, ExtOverAlgebraResult, FaithfulDimensionResult,
    FiniteAddMResolutionResult, GlobalDimensionOfAlgebraResult, IyamaGeneratorError,
    ModuleLengthError, NakayamaAutomorphismOrder, NumberOfComplementsResult,
    PredecessorMultiplicity, PredecessorsOfModuleResult,
    ProjectiveResolutionOfSimpleModuleOverEndoResult,
    ProjectiveResolutionOfSimpleModuleOverEndoStatus, TensorProductOfModulesResult,
    TiltingModuleResult, algebra_as_module_of_enveloping_algebra,
    algebra_as_module_over_enveloping_algebra, all_complements_of_almost_complete_cotilting_module,
    all_complements_of_almost_complete_tilting_module, all_indec_modules_of_length_at_most,
    all_modules_of_length_at_most, all_modules_of_length_plus_one, almost_split_sequence,
    almost_split_sequence_from_side, almost_split_sequence_in_perp_t, annihilator_of_module,
    basis_of_projectives, canonical_basis_coefficients, cartan_matrix, complexity_of_algebra,
    complexity_of_module, cotilting_module, coxeter_matrix, coxeter_polynomial, d_tr,
    d_tr_homomorphism, d_tr_homomorphism_iterated, d_tr_iterated, dominant_dimension_of_algebra,
    dominant_dimension_of_module, dual_of_algebra_as_module_over_enveloping_algebra,
    dual_of_module, dual_of_module_homomorphism, element_in_indec_projective,
    element_in_v_a_as_element_in_indec_proj, end_modulo_projectives, ext_algebra_generators,
    ext_over_algebra, faithful_dimension, finite_canonical_basis, finite_qpa_default_basis,
    first_syzygy, frobenius_form, frobenius_linear_functional,
    from_identity_to_double_star_homomorphism, from_matrix_to_homomorphism_of_projectives,
    global_dimension_from_known_properties, global_dimension_of_algebra, gorenstein_dimension,
    gorenstein_dimension_of_algebra, have_finite_coresolution_in_add_m,
    have_finite_resolution_in_add_m, hom_factoring_through_projectives, hom_from_projective,
    indec_injective_modules, indec_projective_modules, inj_dimension, inj_dimension_of_module,
    injective_envelope, irreducible_morphisms_ending_in, irreducible_morphisms_starting_in,
    is_cotilting_module, is_exceptional_module, is_finite_global_dimension_from_known_properties,
    is_hereditary_algebra, is_injective_module, is_nth_syzygy, is_omega_periodic,
    is_omega_periodic_with_progress, is_projective_module, is_rigid_module, is_schurian_algebra,
    is_selfinjective_algebra, is_semicommutative_algebra, is_symmetric_algebra, is_tau_periodic,
    is_tau_rigid_module, is_tilting_module, is_weakly_symmetric_algebra, iyama_generator,
    left_approximation_by_add_m, left_approximation_by_add_t_hat, left_fac_m_approximation,
    left_mutation_of_cotilting_module_complement, left_mutation_of_tilting_module_complement,
    left_sub_m_approximation, lifting_morphism_from_projective,
    matrix_of_homomorphism_between_projectives, minimal_left_approximation,
    minimal_left_fac_m_approximation, minimal_left_sub_m_approximation,
    minimal_right_approximation, minimal_right_sub_m_approximation, n_rigid_module,
    nakayama_automorphism, nakayama_functor_of_module, nakayama_functor_of_module_homomorphism,
    nakayama_permutation, number_of_complements_of_almost_complete_cotilting_module,
    number_of_complements_of_almost_complete_tilting_module, opposite_nakayama_functor_of_module,
    opposite_nakayama_functor_of_module_homomorphism, order_of_nakayama_automorphism,
    partial_iyama_generator, predecessors_of_module, proj_dimension, proj_dimension_of_module,
    projective_cover, projective_resolution_of_simple_module_over_endo,
    radical_right_approximation_by_add_m, random_module, right_approximation_by_add_m,
    right_approximation_by_add_m_list, right_approximation_by_perp_t,
    right_mutation_of_cotilting_module_complement, right_mutation_of_tilting_module_complement,
    right_sub_m_approximation, simple_modules, socle_series, star_of_module,
    star_of_module_homomorphism, tensor_product_of_modules, tilting_module, tr_d,
    tr_d_homomorphism, tr_d_homomorphism_iterated, tr_d_iterated, transpose_of_module,
    transpose_of_module_homomorphism,
};
use crate::{admissible_quotient_is_semisimple, enveloping_algebra, quotient_by_relations};

fn vertex(index: usize) -> VertexIndex {
    VertexIndex::new(index).unwrap()
}

#[test]
fn extension_fields_run_qpa_syzygy_and_projective_dimension_algorithms() {
    for (prime, degree) in [(2, 2), (2, 3), (3, 2)] {
        let field = GaloisField::conway(prime, degree).unwrap();
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field, quiver).unwrap());
        let simple_one = right_module_over_path_algebra(algebra, vec![1, 0], &[])
            .unwrap()
            .unwrap();
        assert_eq!(
            first_syzygy(&simple_one).unwrap().dimension_vector(),
            [0, 1]
        );
        assert_eq!(
            proj_dimension_of_module(&simple_one, &Integer::from(3)).unwrap(),
            Some(1)
        );
    }
}

fn q(value: i64) -> Rational {
    Q.nth(Integer::from(value))
}

fn global_dimension_json(result: Option<GlobalDimensionOfAlgebraResult>) -> Value {
    match result {
        Some(GlobalDimensionOfAlgebraResult::Finite(dimension)) => dimension.into(),
        Some(GlobalDimensionOfAlgebraResult::Infinite) => "infinity".into(),
        Some(GlobalDimensionOfAlgebraResult::ExceedsBound) => false.into(),
        None => Value::Null,
    }
}

fn cached_global_dimension_json(result: Option<AlgebraGlobalDimension>) -> Value {
    match result {
        Some(AlgebraGlobalDimension::Finite(dimension)) => dimension.into(),
        Some(AlgebraGlobalDimension::Infinite) => "infinity".into(),
        None => Value::Null,
    }
}

fn bounded_dimension_json(result: Option<usize>) -> Value {
    result.map_or(Value::Bool(false), Into::into)
}

fn dominant_dimension_json(result: Option<DominantDimensionResult>) -> Value {
    match result {
        Some(DominantDimensionResult::Finite(dimension)) => dimension.into(),
        Some(DominantDimensionResult::Infinite) => "infinity".into(),
        Some(DominantDimensionResult::ExceedsBound) => false.into(),
        None => Value::Null,
    }
}

fn fresh_a2_module<F: Field>(field: &F, dimension_vector: Vec<usize>) -> PathAlgebraModule<F> {
    let algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    right_module_over_path_algebra(algebra, dimension_vector, &[])
        .unwrap()
        .unwrap()
}

fn fresh_a2_zero<F: Field>(field: &F) -> PathAlgebraModule<F> {
    let algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "z")]).unwrap(),
        )
        .unwrap(),
    );
    zero_module(&algebra).unwrap()
}

fn fresh_dual_simple<F: Field>(field: &F) -> PathAlgebraModule<F> {
    let algebra = fresh_dual_numbers(field);
    right_module_over_path_algebra(
        algebra,
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![field.zero()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap()
}

fn fresh_dual_numbers<F: Field>(field: &F) -> QuiverAlgebra<F> {
    let parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = parent.generators()[1].clone();
    QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.pow(2).unwrap()]).unwrap())
}

fn fresh_radical_square_zero_a3<F: Field>(field: &F) -> QuiverAlgebra<F> {
    let parent = PathAlgebra::new(
        field.clone(),
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
    let arrows = parent.quiver().arrows().collect::<Vec<_>>();
    let x = parent.element(arrows[0].path()).unwrap();
    let y = parent.element(arrows[1].path()).unwrap();
    QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&y).unwrap()]).unwrap())
}

fn fresh_radical_square_zero_cycle<F: Field>(field: &F) -> QuiverAlgebra<F> {
    let parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "p"),
                ArrowSpec::named(vertex(2), vertex(1), "q"),
            ],
        )
        .unwrap(),
    )
    .unwrap();
    let arrows = parent.quiver().arrows().collect::<Vec<_>>();
    let p = parent.element(arrows[0].path()).unwrap();
    let q = parent.element(arrows[1].path()).unwrap();
    QuiverAlgebra::Quotient(
        quotient_by_relations(parent, &[p.multiply(&q).unwrap(), q.multiply(&p).unwrap()]).unwrap(),
    )
}

fn isolated_path_algebra<F: Field>(field: F, vertex_count: usize) -> QuiverAlgebra<F> {
    let vertices = (1..=vertex_count)
        .map(|index| format!("v{index}").into_boxed_str())
        .collect();
    QuiverAlgebra::Path(
        PathAlgebra::new(field, Quiver::from_named(vertices, Vec::new()).unwrap()).unwrap(),
    )
}

struct ZeroRng;

impl RngCore for ZeroRng {
    fn next_u32(&mut self) -> u32 {
        0
    }

    fn next_u64(&mut self) -> u64 {
        0
    }

    fn fill_bytes(&mut self, destination: &mut [u8]) {
        destination.fill(0);
    }
}

#[test]
fn regular_bimodule_for_dual_numbers_has_qpa_row_actions() {
    let quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let parent = PathAlgebra::new(Q, quiver).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let algebra =
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap());

    let module = algebra_as_module_over_enveloping_algebra(&algebra)
        .unwrap()
        .unwrap();
    let again = algebra_as_module_of_enveloping_algebra(&algebra)
        .unwrap()
        .unwrap();
    let enveloping = enveloping_algebra(&algebra).unwrap();

    assert_eq!(module.dimension_vector(), &[2]);
    assert_eq!(module.arrow_matrices().len(), 2);
    assert_eq!(
        module.arrow_matrices()[0].as_flat(),
        &[q(0), q(1), q(0), q(0)]
    );
    assert_eq!(
        module.arrow_matrices()[1].as_flat(),
        &[q(0), q(1), q(0), q(0)]
    );
    assert!(module.is_identical(&again));
    assert!(
        module
            .retained_tensor_product()
            .unwrap()
            .is_identical(&enveloping)
    );

    let generic_dual = dual_of_module(&module).unwrap();
    let generic_back = dual_of_module(&generic_dual).unwrap();
    let bimodule_dual = dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
        .unwrap()
        .unwrap();
    assert_eq!(
        generic_dual.arrow_matrices()[0].as_flat(),
        &[q(0), q(0), q(1), q(0)]
    );
    assert!(generic_back.is_identical(&module));
    assert_eq!(
        bimodule_dual.arrow_matrices()[0].as_flat(),
        &[q(0), q(0), q(1), q(0)]
    );
}

#[test]
fn dual_of_module_homomorphism_transposes_reverses_and_caches_like_qpa() {
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
    let map = right_module_hom_over_algebra(
        module.clone(),
        module,
        vec![
            DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
            DenseMatrix::identity(1, Q),
        ],
    )
    .unwrap();

    let dual = dual_of_module_homomorphism(&map).unwrap();
    assert_eq!(dual.source().dimension_vector(), &[2, 1]);
    assert_eq!(dual.range().dimension_vector(), &[2, 1]);
    assert_eq!(
        dual.vertex_matrices(),
        &[
            DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
            DenseMatrix::identity(1, Q),
        ]
    );
    assert!(dual.is_identical(&dual_of_module_homomorphism(&map).unwrap()));
    let double_dual = dual_of_module_homomorphism(&dual).unwrap();
    assert_eq!(double_dual, map);
    assert!(
        !double_dual.is_identical(&map),
        "the QPA source does not install a reciprocal map attribute"
    );

    let e = |value: i64| Z2.nth(Integer::from(value));
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
    let module = right_module_over_path_algebra(
        algebra,
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![e(1)], vec![e(0)]], Z2).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let map = right_module_hom_over_algebra(
        module.clone(),
        module,
        vec![
            DenseMatrix::from_rows(vec![vec![e(1), e(0)], vec![e(0), e(0)]], Z2).unwrap(),
            DenseMatrix::identity(1, Z2),
        ],
    )
    .unwrap();
    assert_eq!(
        dual_of_module_homomorphism(&map).unwrap().vertex_matrices(),
        map.vertex_matrices()
    );
}

#[test]
fn regular_bimodule_supports_characteristic_two_and_zero_spaces() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
    let module = algebra_as_module_over_enveloping_algebra(&algebra)
        .unwrap()
        .unwrap();

    assert_eq!(module.dimension_vector(), &[1, 1, 0, 1]);
    assert_eq!(module.dimension(), 3);
    assert_eq!(module.arrow_matrices().len(), 4);
    assert!(module.arrow_matrices().iter().all(|matrix| {
        matrix
            .as_flat()
            .iter()
            .all(|entry| Z2.to_integer(entry) == 0 || Z2.to_integer(entry) == 1)
    }));

    let dual = dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
        .unwrap()
        .unwrap();
    assert_eq!(dual.dimension_vector(), &[1, 0, 1, 1]);
}

#[test]
fn cyclic_path_algebra_returns_qpa_fail_before_enveloping_construction() {
    let quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());

    assert!(
        algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .is_none()
    );
    assert!(algebra.cached_enveloping_algebra().is_none());
    assert!(
        dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .is_none()
    );
    assert!(algebra.has_dual_enveloping_module_attribute());
}

fn matrix_json<F: Field>(
    matrix: &DenseMatrix<F>,
    coefficient: impl Fn(&F::Element) -> String,
) -> Value {
    json!({
        "shape": [matrix.nrows(), matrix.ncols()],
        "rows": (0..matrix.nrows())
            .map(|row| matrix.row(row).iter().map(&coefficient).collect::<Vec<_>>())
            .collect::<Vec<_>>(),
    })
}

fn module_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!({
        "dimension_vector": module.dimension_vector(),
        "dimension": module.dimension(),
        "arrows": module
            .right_acting_algebra()
            .quiver()
            .arrows()
            .map(|arrow| arrow.name().to_owned())
            .collect::<Vec<_>>(),
        "matrices": module
            .arrow_matrices()
            .iter()
            .map(|matrix| matrix_json(matrix, coefficient))
            .collect::<Vec<_>>(),
    })
}

fn length_module_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    json!({
        "dimension": module.dimension_vector(),
        "matrices": module.matrices_of_path_algebra_module().iter().map(|matrix| {
            match matrix {
                PathAlgebraModuleMatrix::ZeroScalar => json!([]),
                PathAlgebraModuleMatrix::Matrix(matrix) => json!(
                    (0..matrix.nrows()).map(|row| matrix.row(row).iter()
                        .map(coefficient).collect::<Vec<_>>()).collect::<Vec<_>>()
                ),
            }
        }).collect::<Vec<_>>(),
    })
}

fn length_module_list_json<F: Field>(
    modules: &[PathAlgebraModule<F>],
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    json!(
        modules
            .iter()
            .map(|module| length_module_json(module, coefficient))
            .collect::<Vec<_>>()
    )
}

#[allow(clippy::needless_pass_by_value)]
fn iyama_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let a2_module = right_module_over_path_algebra(
        a2,
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![field.one()], vec![field.zero()]], field.clone())
                .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let vector_algebra = QuiverAlgebra::Path(
        PathAlgebra::new(field.clone(), Quiver::from_indexed(1, Vec::new()).unwrap()).unwrap(),
    );
    let vector = right_module_over_path_algebra(vector_algebra, vec![2], &[])
        .unwrap()
        .unwrap();
    let dual_algebra = fresh_dual_numbers(&field);
    let dual = right_module_over_path_algebra(
        dual_algebra,
        vec![2],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(
                vec![
                    vec![field.zero(), field.one()],
                    vec![field.zero(), field.zero()],
                ],
                field,
            )
            .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    json!({
        "a2_partial": length_module_json(
            &partial_iyama_generator(&a2_module).unwrap(), coefficient
        ),
        "a2_iyama": length_module_json(&iyama_generator(&a2_module).unwrap(), coefficient),
        "vector_partial": length_module_json(
            &partial_iyama_generator(&vector).unwrap(), coefficient
        ),
        "vector_iyama": length_module_json(&iyama_generator(&vector).unwrap(), coefficient),
        "dual_partial": length_module_json(
            &partial_iyama_generator(&dual).unwrap(), coefficient
        ),
        "dual_iyama": length_module_json(&iyama_generator(&dual).unwrap(), coefficient),
    })
}

#[test]
fn pinned_qpa_iyama_generator_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_iyama.schema-v1.json"
    ))
    .unwrap();
    let rational = iyama_case(Q, |entry: &Rational| {
        json!(entry.numerator_ref().to_i64().unwrap())
    });
    let binary = iyama_case(Z2, |entry| json!(Z2.to_integer(entry).to_i64().unwrap()));

    assert_eq!(rational, expected["rational"]);
    assert_eq!(binary, expected["binary"]);
}

#[test]
fn finite_field_zero_iyama_preserves_qpa_failed_list_index() {
    let algebra = QuiverAlgebra::Path(
        PathAlgebra::new(Z2, Quiver::from_indexed(1, Vec::new()).unwrap()).unwrap(),
    );
    let zero = right_module_over_path_algebra(algebra, vec![0], &[])
        .unwrap()
        .unwrap();

    assert!(matches!(
        iyama_generator(&zero),
        Err(IyamaGeneratorError::MissingMultiplicityDecomposition)
    ));
}

#[allow(clippy::needless_pass_by_value)]
fn modulehom_length_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simples = simple_modules(&a2).unwrap();
    let a2_plus = all_modules_of_length_plus_one(&simples).unwrap();
    let a2_at_most_two = all_modules_of_length_at_most(&a2, &Integer::from(2)).unwrap();
    let a2_indec_two = all_indec_modules_of_length_at_most(&a2, &Integer::from(2)).unwrap();
    let a2_negative = all_modules_of_length_at_most(&a2, &Integer::from(-1)).unwrap();
    let a2_zero = all_modules_of_length_at_most(&a2, &Integer::from(0)).unwrap();

    let loop_parent = PathAlgebra::new(
        field,
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simples = simple_modules(&dual).unwrap();
    let dual_length_two = all_modules_of_length_plus_one(&dual_simples).unwrap();
    let dual_length_three = all_modules_of_length_plus_one(&dual_length_two).unwrap();

    json!({
        "a2_plus": length_module_list_json(&a2_plus, coefficient),
        "a2_at_most_two": length_module_list_json(&a2_at_most_two, coefficient),
        "a2_indec_two": length_module_list_json(&a2_indec_two, coefficient),
        "a2_negative": length_module_list_json(&a2_negative, coefficient),
        "a2_zero": length_module_list_json(&a2_zero, coefficient),
        "dual_length_two": length_module_list_json(&dual_length_two, coefficient),
        "dual_length_three": length_module_list_json(&dual_length_three, coefficient),
    })
}

#[test]
fn pinned_qpa_modulehom_length_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_length.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let rational_a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let rational_zero = zero_module(&rational_a2).unwrap();
    let rational_projective = right_module_over_path_algebra(
        rational_a2,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![Q.one()]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehom_length",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehom_length_case(Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": modulehom_length_case(field5, |entry| json!(
            field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "rational_indecomposable": {
            "zero": is_indecomposable_module(&rational_zero).unwrap().unwrap(),
            "projective": is_indecomposable_module(&rational_projective).unwrap().unwrap(),
        },
    });

    assert_eq!(actual, expected);
}

#[test]
fn module_length_enumeration_preserves_qpa_validation_order() {
    assert!(matches!(
        all_modules_of_length_plus_one::<Z2>(&[]),
        Err(ModuleLengthError::EmptyModuleList)
    ));

    let rational = QuiverAlgebra::Path(
        PathAlgebra::new(Q, Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap()).unwrap(),
    );
    assert!(matches!(
        all_modules_of_length_at_most(&rational, &Integer::from(1)),
        Err(ModuleLengthError::FiniteFieldRequired)
    ));

    let first_algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Z2,
            Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
        )
        .unwrap(),
    );
    let second_algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Z2,
            Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
        )
        .unwrap(),
    );
    let first = simple_modules(&first_algebra).unwrap()[0].clone();
    let foreign = simple_modules(&second_algebra).unwrap()[0].clone();
    assert!(matches!(
        all_modules_of_length_plus_one(&[first.clone(), foreign]),
        Err(ModuleLengthError::ModulesOverDifferentAlgebras)
    ));
    let double = direct_sum_of_qpa_modules(&[first.clone(), first.clone()])
        .unwrap()
        .unwrap();
    assert!(matches!(
        all_modules_of_length_plus_one(&[first, double]),
        Err(ModuleLengthError::ModulesOfDifferentLengths)
    ));
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn module_complexity_case<F: QpaCoefficientField>(
    field: F,
    include_finite_only_case: bool,
) -> Value {
    let zero = field.zero();
    let one = field.one();
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let zero_module = zero_module(&a2).unwrap();
    let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let projective = right_module_over_path_algebra(
        a2.clone(),
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let semisimple_pair = right_module_over_path_algebra(
        a2,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual.clone(),
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let dual_regular = right_module_over_path_algebra(
        dual.clone(),
        vec![2],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(
                vec![vec![zero.clone(), one.clone()], vec![zero.clone(), zero]],
                field,
            )
            .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let common = json!({
        "exceptional_zero": is_exceptional_module(&zero_module).unwrap().unwrap(),
        "exceptional_simple_one": is_exceptional_module(&simple_one).unwrap().unwrap(),
        "exceptional_projective": is_exceptional_module(&projective).unwrap().unwrap(),
        "exceptional_dual_simple": is_exceptional_module(&dual_simple).unwrap().unwrap(),
        "complexity_a2_simple_n2": complexity_of_module(&simple_one, &Integer::from(2)).unwrap(),
        "complexity_dual_simple_n2": complexity_of_module(&dual_simple, &Integer::from(2)).unwrap(),
        "complexity_dual_simple_n3": complexity_of_module(&dual_simple, &Integer::from(3)).unwrap(),
        "complexity_dual_simple_n10": complexity_of_module(&dual_simple, &Integer::from(10)).unwrap(),
        "complexity_dual_regular_n2": complexity_of_module(&dual_regular, &Integer::from(2)).unwrap(),
        "complexity_dual_algebra_n2": complexity_of_algebra(&dual, &Integer::from(2)).unwrap(),
        "complexity_dual_algebra_n10": complexity_of_algebra(&dual, &Integer::from(10)).unwrap(),
    });
    if include_finite_only_case {
        json!({
            "common": common,
            "exceptional_semisimple_pair": is_exceptional_module(&semisimple_pair)
                .unwrap().unwrap(),
        })
    } else {
        common
    }
}

#[test]
fn pinned_qpa_module_complexity_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/module_complexity.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "module_complexity",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": module_complexity_case(Z2, true),
        "quinary": module_complexity_case(Zp::new(5), true),
        "rational": module_complexity_case(Q, false),
    });

    assert_eq!(actual, expected);
}

#[test]
fn module_complexity_preserves_qpa_dispatch_and_cache_order() {
    let loop_parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = simple_modules(&dual).unwrap()[0].clone();
    let dual_zero = zero_module(&dual).unwrap();
    assert_eq!(
        complexity_of_module(&dual_simple, &Integer::from(-2)).unwrap(),
        0
    );
    assert_eq!(
        complexity_of_module(&dual_zero, &Integer::from(-1)).unwrap(),
        0
    );
    assert!(matches!(
        complexity_of_module(&dual_simple, &Integer::from(1)),
        Err(EnvelopingModuleError::UndefinedComplexityEstimate { .. })
    ));

    let nonadmissible_parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "y")]).unwrap(),
    )
    .unwrap();
    let y = nonadmissible_parent.generators()[1].clone();
    let nonadmissible =
        QuiverAlgebra::Quotient(quotient_by_relations(nonadmissible_parent, &[y]).unwrap());
    assert!(matches!(
        complexity_of_algebra(&nonadmissible, &Integer::from(2)),
        Err(EnvelopingModuleError::ComplexityRequiresAdmissibleQuotient)
    ));

    let a2_parent = PathAlgebra::new(
        Q,
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
    )
    .unwrap();
    let a2_quotient = QuiverAlgebra::Quotient(quotient_by_relations(a2_parent, &[]).unwrap());
    assert!(matches!(
        global_dimension_of_algebra(&a2_quotient, 2).unwrap(),
        Some(GlobalDimensionOfAlgebraResult::Finite(1))
    ));
    assert_eq!(
        complexity_of_algebra(&a2_quotient, &Integer::from(1)).unwrap(),
        0
    );

    let rational_pair = right_module_over_path_algebra(
        QuiverAlgebra::Path(
            PathAlgebra::new(
                Q,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        ),
        vec![0, 2],
        &[],
    )
    .unwrap()
    .unwrap();
    assert_eq!(is_exceptional_module(&rational_pair).unwrap(), None);
    let rational_simple = right_module_over_path_algebra(
        rational_pair.right_acting_algebra().clone(),
        vec![1, 0],
        &[],
    )
    .unwrap()
    .unwrap();
    assert_eq!(is_exceptional_module(&rational_simple).unwrap(), Some(true));
    assert_eq!(rational_simple.cached_is_exceptional_module(), Some(true));
    assert_eq!(is_exceptional_module(&rational_simple).unwrap(), Some(true));
}

fn periodicity_json(period: Option<Integer>) -> Value {
    period.map_or(Value::Bool(false), |degree| {
        serde_json::from_str(&degree.to_string()).unwrap()
    })
}

#[allow(clippy::needless_pass_by_value)]
fn modulehomalg_periodicity_case<F: QpaCoefficientField>(field: F) -> Value {
    let zero = field.zero();
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let zero_module = zero_module(&a2).unwrap();
    let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let simple_two = right_module_over_path_algebra(a2, vec![0, 1], &[])
        .unwrap()
        .unwrap();

    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual,
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![zero]], field).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    json!({
        "omega_zero_n0": periodicity_json(is_omega_periodic(&zero_module, &Integer::from(0)).unwrap()),
        "omega_zero_negative": periodicity_json(is_omega_periodic(&zero_module, &Integer::from(-2)).unwrap()),
        "omega_zero_n3": periodicity_json(is_omega_periodic(&zero_module, &Integer::from(3)).unwrap()),
        "omega_a2_simple_one_n4": periodicity_json(is_omega_periodic(&simple_one, &Integer::from(4)).unwrap()),
        "omega_a2_simple_two_n4": periodicity_json(is_omega_periodic(&simple_two, &Integer::from(4)).unwrap()),
        "omega_dual_simple_n1": periodicity_json(is_omega_periodic(&dual_simple, &Integer::from(1)).unwrap()),
        "omega_dual_simple_n4": periodicity_json(is_omega_periodic(&dual_simple, &Integer::from(4)).unwrap()),
        "tau_zero_n0": periodicity_json(is_tau_periodic(&zero_module, &Integer::from(0)).unwrap()),
        "tau_zero_negative": periodicity_json(is_tau_periodic(&zero_module, &Integer::from(-2)).unwrap()),
        "tau_zero_n3": periodicity_json(is_tau_periodic(&zero_module, &Integer::from(3)).unwrap()),
        "tau_a2_simple_one_n4": periodicity_json(is_tau_periodic(&simple_one, &Integer::from(4)).unwrap()),
        "tau_a2_simple_two_n4": periodicity_json(is_tau_periodic(&simple_two, &Integer::from(4)).unwrap()),
        "tau_dual_simple_n1": periodicity_json(is_tau_periodic(&dual_simple, &Integer::from(1)).unwrap()),
        "tau_dual_simple_n4": periodicity_json(is_tau_periodic(&dual_simple, &Integer::from(4)).unwrap()),
    })
}

#[test]
fn pinned_qpa_modulehomalg_periodicity_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_periodicity.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_periodicity",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_periodicity_case(Z2),
        "quinary": modulehomalg_periodicity_case(Zp::new(5)),
        "rational": modulehomalg_periodicity_case(Q),
    });

    assert_eq!(actual, expected);
}

#[test]
fn omega_periodicity_exposes_qpa_progress_without_global_output() {
    let algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Z2,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simple_one = right_module_over_path_algebra(algebra, vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let mut progress = Vec::new();
    assert_eq!(
        is_omega_periodic_with_progress(&simple_one, &Integer::from(4), |degree| {
            progress.push(degree.clone());
        })
        .unwrap(),
        None
    );
    assert_eq!(progress, (1..=4).map(Integer::from).collect::<Vec<_>>());

    progress.clear();
    assert_eq!(
        is_omega_periodic_with_progress(&simple_one, &Integer::from(-1), |degree| {
            progress.push(degree.clone());
        })
        .unwrap(),
        None
    );
    assert!(progress.is_empty());
}

fn double_star_map_json<F: Field>(
    map: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    json!({
        "source_dimension": map.source().dimension_vector(),
        "range_dimension": map.range().dimension_vector(),
        "matrices": map.vertex_matrices().iter().map(|matrix| json!(
            (0..matrix.nrows()).map(|row| matrix.row(row).iter()
                .map(coefficient).collect::<Vec<_>>()).collect::<Vec<_>>()
        )).collect::<Vec<_>>(),
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehomalg_approximation_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let zero = field.zero();
    let one = field.one();
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let a2_zero = zero_module(&a2).unwrap();
    let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let projective = right_module_over_path_algebra(
        a2,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let repeated = direct_sum_of_qpa_modules(&[simple_one.clone(), simple_one.clone()])
        .unwrap()
        .unwrap();

    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual.clone(),
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let dual_regular = right_module_over_path_algebra(
        dual,
        vec![2],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(
                vec![vec![zero.clone(), one.clone()], vec![zero.clone(), zero]],
                field,
            )
            .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let evaluate = |map| double_star_map_json(&map, coefficient);
    json!({
        "right_zero": evaluate(right_approximation_by_add_m(&simple_one, &a2_zero).unwrap()),
        "right_projective_to_simple": evaluate(right_approximation_by_add_m(&projective, &simple_one).unwrap()),
        "right_simple_to_projective": evaluate(right_approximation_by_add_m(&simple_one, &projective).unwrap()),
        "right_repeated": evaluate(right_approximation_by_add_m(&repeated, &simple_one).unwrap()),
        "minimal_right_repeated": evaluate(minimal_right_approximation(&repeated, &simple_one).unwrap()),
        "left_zero": evaluate(left_approximation_by_add_m(&a2_zero, &simple_one).unwrap()),
        "left_projective_to_simple": evaluate(left_approximation_by_add_m(&projective, &simple_one).unwrap()),
        "left_simple_to_projective": evaluate(left_approximation_by_add_m(&simple_one, &projective).unwrap()),
        "left_repeated": evaluate(left_approximation_by_add_m(&simple_one, &repeated).unwrap()),
        "minimal_left_repeated": evaluate(minimal_left_approximation(&simple_one, &repeated).unwrap()),
        "right_dual_regular_to_simple": evaluate(right_approximation_by_add_m(&dual_regular, &dual_simple).unwrap()),
        "left_dual_simple_to_regular": evaluate(left_approximation_by_add_m(&dual_simple, &dual_regular).unwrap()),
    })
}

#[test]
fn pinned_qpa_modulehomalg_approximation_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_approximations.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_approximations",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_approximation_case(Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": modulehomalg_approximation_case(field5, |entry| json!(
            field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "rational": modulehomalg_approximation_case(Q, |entry| json!([
            entry.numerator_ref().to_string().parse::<i64>().unwrap(),
            entry.denominator_ref().to_string().parse::<i64>().unwrap(),
        ])),
    });

    assert_eq!(actual, expected);
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehomalg_left_fac_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let zero_scalar = field.zero();
    let one_scalar = field.one();
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let zero = zero_module(&a2).unwrap();
    let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let projective = right_module_over_path_algebra(
        a2,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one_scalar.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual.clone(),
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![zero_scalar.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let dual_regular = right_module_over_path_algebra(
        dual,
        vec![2],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(
                vec![
                    vec![zero_scalar.clone(), one_scalar],
                    vec![zero_scalar.clone(), zero_scalar],
                ],
                field,
            )
            .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let evaluate = |source: &PathAlgebraModule<F>, module: &PathAlgebraModule<F>| {
        json!({
            "left": double_star_map_json(
                &left_fac_m_approximation(source, module).unwrap(),
                coefficient,
            ),
            "minimal": double_star_map_json(
                &minimal_left_fac_m_approximation(source, module).unwrap(),
                coefficient,
            ),
        })
    };
    json!({
        "zero_to_simple": evaluate(&zero, &simple_one),
        "simple_to_zero": evaluate(&simple_one, &zero),
        "projective_to_simple": evaluate(&projective, &simple_one),
        "simple_to_projective": evaluate(&simple_one, &projective),
        "dual_regular_to_simple": evaluate(&dual_regular, &dual_simple),
    })
}

#[test]
fn pinned_qpa_modulehomalg_left_fac_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_left_fac.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_left_fac",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_left_fac_case(Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": modulehomalg_left_fac_case(field5, |entry| json!(
            field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "rational": modulehomalg_left_fac_case(Q, |entry| json!([
            entry.numerator_ref().to_string().parse::<i64>().unwrap(),
            entry.denominator_ref().to_string().parse::<i64>().unwrap(),
        ])),
    });
    assert_eq!(actual, expected);
}

#[test]
fn left_fac_approximations_reject_foreign_algebras() {
    let first = fresh_a2_zero(&Z2);
    let second = fresh_a2_zero(&Z2);
    assert!(matches!(
        left_fac_m_approximation(&first, &second),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
    assert!(matches!(
        minimal_left_fac_m_approximation(&first, &second),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehomalg_sub_approximations_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let zero_scalar = field.zero();
    let one_scalar = field.one();
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let zero = zero_module(&a2).unwrap();
    let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let projective = right_module_over_path_algebra(
        a2,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one_scalar.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual.clone(),
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![zero_scalar.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let dual_regular = right_module_over_path_algebra(
        dual,
        vec![2],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(
                vec![
                    vec![zero_scalar.clone(), one_scalar],
                    vec![zero_scalar.clone(), zero_scalar],
                ],
                field,
            )
            .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let evaluate = |source: &PathAlgebraModule<F>, module: &PathAlgebraModule<F>| {
        json!({
            "right": double_star_map_json(
                &right_sub_m_approximation(module, source).unwrap(),
                coefficient,
            ),
            "minimal_right": double_star_map_json(
                &minimal_right_sub_m_approximation(module, source).unwrap(),
                coefficient,
            ),
            "left": double_star_map_json(
                &left_sub_m_approximation(source, module).unwrap(),
                coefficient,
            ),
            "minimal_left": double_star_map_json(
                &minimal_left_sub_m_approximation(source, module).unwrap(),
                coefficient,
            ),
        })
    };
    json!({
        "zero_and_simple": evaluate(&zero, &simple_one),
        "simple_and_zero": evaluate(&simple_one, &zero),
        "projective_and_simple": evaluate(&projective, &simple_one),
        "dual_regular_and_simple": evaluate(&dual_regular, &dual_simple),
    })
}

#[test]
fn pinned_qpa_modulehomalg_sub_approximations_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_sub_approximations.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_sub_approximations",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_sub_approximations_case(Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": modulehomalg_sub_approximations_case(field5, |entry| json!(
            field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "rational": modulehomalg_sub_approximations_case(Q, |entry| json!([
            entry.numerator_ref().to_string().parse::<i64>().unwrap(),
            entry.denominator_ref().to_string().parse::<i64>().unwrap(),
        ])),
    });
    assert_eq!(actual, expected);
}

#[test]
fn sub_m_approximations_reject_foreign_algebras() {
    let first = fresh_a2_zero(&Z2);
    let second = fresh_a2_zero(&Z2);
    assert!(matches!(
        right_sub_m_approximation(&first, &second),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
    assert!(matches!(
        minimal_right_sub_m_approximation(&first, &second),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
    assert!(matches!(
        left_sub_m_approximation(&first, &second),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
    assert!(matches!(
        minimal_left_sub_m_approximation(&first, &second),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
}

fn add_m_coresolution_json<F: Field>(
    result: Option<AddMCoresolution<F>>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let Some(coresolution) = result else {
        return Value::Bool(false);
    };
    json!({
        "degrees": {
            "-3": double_star_map_json(
                &coresolution.differential(&Integer::from(-3)).unwrap(),
                coefficient,
            ),
            "-2": double_star_map_json(
                &coresolution.differential(&Integer::from(-2)).unwrap(),
                coefficient,
            ),
            "-1": double_star_map_json(
                &coresolution.differential(&Integer::from(-1)).unwrap(),
                coefficient,
            ),
            "0": double_star_map_json(
                &coresolution.differential(&Integer::from(0)).unwrap(),
                coefficient,
            ),
            "1": double_star_map_json(
                &coresolution.differential(&Integer::from(1)).unwrap(),
                coefficient,
            ),
        }
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehomalg_add_m_coresolution_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let zero_scalar = field.zero();
    let one_scalar = field.one();
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let projective = right_module_over_path_algebra(
        a2,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one_scalar.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual.clone(),
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![zero_scalar.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let dual_regular = right_module_over_path_algebra(
        dual,
        vec![2],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(
                vec![
                    vec![zero_scalar.clone(), one_scalar],
                    vec![zero_scalar.clone(), zero_scalar],
                ],
                field,
            )
            .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let evaluate = |source, module, bound| {
        add_m_coresolution_json(
            have_finite_coresolution_in_add_m(source, module, &Integer::from(bound)).unwrap(),
            coefficient,
        )
    };
    json!({
        "simple_in_simple_0": evaluate(&simple_one, &simple_one, 0),
        "projective_in_simple_0": evaluate(&projective, &simple_one, 0),
        "projective_in_simple_1": evaluate(&projective, &simple_one, 1),
        "dual_simple_in_regular_1": evaluate(&dual_simple, &dual_regular, 1),
        "dual_regular_in_simple_2": evaluate(&dual_regular, &dual_simple, 2),
    })
}

#[test]
fn pinned_qpa_modulehomalg_add_m_coresolution_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_add_m_coresolution.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_add_m_coresolution",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_add_m_coresolution_case(Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": modulehomalg_add_m_coresolution_case(field5, |entry| json!(
            field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "rational": modulehomalg_add_m_coresolution_case(Q, |entry| json!([
            entry.numerator_ref().to_string().parse::<i64>().unwrap(),
            entry.denominator_ref().to_string().parse::<i64>().unwrap(),
        ])),
    });
    assert_eq!(actual, expected);
}

#[test]
fn add_m_coresolution_preserves_negative_bound_and_foreign_algebra_failures() {
    let module = fresh_dual_simple(&Z2);
    assert!(matches!(
        have_finite_coresolution_in_add_m(&module, &module, &Integer::from(-1)),
        Err(EnvelopingModuleError::EmptyCoresolutionDifferentialList { bound })
            if bound == -1
    ));

    let first = fresh_a2_zero(&Z2);
    let second = fresh_a2_zero(&Z2);
    assert!(matches!(
        have_finite_coresolution_in_add_m(&first, &second, &Integer::from(0)),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
}

fn tilting_coresolution_shape<F: Field>(coresolution: &AddMCoresolution<F>) -> Value {
    let map_shape = |degree| {
        let differential = coresolution.differential(&Integer::from(degree)).unwrap();
        json!({
            "source": differential.source().dimension_vector(),
            "range": differential.range().dimension_vector(),
        })
    };
    json!({
        "-2": map_shape(-2),
        "-1": map_shape(-1),
        "0": map_shape(0),
    })
}

fn tilting_value_json<F: Field>(result: Option<TiltingModuleResult<F>>) -> Value {
    let Some(result) = result else {
        return Value::Bool(false);
    };
    json!({
        "projective_dimension": result.projective_dimension(),
        "coresolution_count": result.coresolutions().len(),
        "coresolutions": result
            .coresolutions()
            .iter()
            .map(tilting_coresolution_shape)
            .collect::<Vec<_>>(),
    })
}

fn modulehomalg_tilting_case<F: QpaCoefficientField>(field: &F) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let projectives = indec_projective_modules(&a2).unwrap();
    let regular = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
    let simple = simple_modules(&a2).unwrap()[0].clone();
    let pd_one = direct_sum_of_qpa_modules(&[projectives[0].clone(), simple.clone()])
        .unwrap()
        .unwrap();

    let dual = fresh_dual_numbers(field);
    let dual_regular = direct_sum_of_qpa_modules(&indec_projective_modules(&dual).unwrap())
        .unwrap()
        .unwrap();
    let dual_simple = simple_modules(&dual).unwrap()[0].clone();

    let check = |module: &PathAlgebraModule<F>, bound| {
        let has_before = is_tilting_module(module).is_some();
        let value = tilting_value_json(tilting_module(module, &Integer::from(bound)).unwrap());
        let has_after = is_tilting_module(module).is_some();
        json!({
            "has_before": has_before,
            "value": value,
            "has_after": has_after,
        })
    };
    json!({
        "regular_bound_0": check(&regular, 0),
        "regular_cached_negative": check(&regular, -1),
        "pd_one_bound_0": check(&pd_one, 0),
        "pd_one_bound_1": check(&pd_one, 1),
        "simple_bound_1": check(&simple, 1),
        "dual_regular_bound_0": check(&dual_regular, 0),
        "dual_simple_bound_1": check(&dual_simple, 1),
    })
}

#[test]
fn pinned_qpa_modulehomalg_tilting_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_tilting.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_tilting",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_tilting_case(&Z2),
        "quinary": modulehomalg_tilting_case(&Zp::new(5)),
        "rational": modulehomalg_tilting_case(&Q),
    });
    assert_eq!(actual, expected);
}

fn add_m_resolution_json<F: Field>(
    result: Option<FiniteAddMResolutionResult<F>>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let Some(result) = result else {
        return Value::Bool(false);
    };
    let FiniteAddMResolutionResult::Resolution(resolution) = result else {
        return Value::Bool(true);
    };
    let map = |degree| {
        double_star_map_json(
            &resolution.differential(&Integer::from(degree)).unwrap(),
            coefficient,
        )
    };
    json!({
        "degrees": {
            "-1": map(-1),
            "0": map(0),
            "1": map(1),
            "2": map(2),
            "3": map(3),
        }
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehomalg_add_m_resolution_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let zero_scalar = field.zero();
    let one_scalar = field.one();
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let projective = right_module_over_path_algebra(
        a2.clone(),
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one_scalar.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let regular = direct_sum_of_qpa_modules(&indec_projective_modules(&a2).unwrap())
        .unwrap()
        .unwrap();

    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual.clone(),
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![zero_scalar.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let dual_regular = right_module_over_path_algebra(
        dual,
        vec![2],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(
                vec![
                    vec![zero_scalar.clone(), one_scalar],
                    vec![zero_scalar.clone(), zero_scalar],
                ],
                field,
            )
            .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let evaluate = |target, module, bound| {
        add_m_resolution_json(
            have_finite_resolution_in_add_m(target, module, &Integer::from(bound)).unwrap(),
            coefficient,
        )
    };
    json!({
        "simple_in_simple_0": evaluate(&simple_one, &simple_one, 0),
        "regular_in_regular_0": evaluate(&regular, &regular, 0),
        "projective_in_simple_0": evaluate(&projective, &simple_one, 0),
        "projective_in_simple_1": evaluate(&projective, &simple_one, 1),
        "dual_simple_in_regular_1": evaluate(&dual_simple, &dual_regular, 1),
        "dual_regular_in_simple_2": evaluate(&dual_regular, &dual_simple, 2),
    })
}

#[test]
fn pinned_qpa_modulehomalg_add_m_resolution_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_add_m_resolution.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_add_m_resolution",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_add_m_resolution_case(Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": modulehomalg_add_m_resolution_case(field5, |entry| json!(
            field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "rational": modulehomalg_add_m_resolution_case(Q, |entry| json!([
            entry.numerator_ref().to_string().parse::<i64>().unwrap(),
            entry.denominator_ref().to_string().parse::<i64>().unwrap(),
        ])),
    });
    assert_eq!(actual, expected);
}

#[test]
fn add_m_resolution_preserves_negative_bound_and_foreign_algebra_failures() {
    let module = fresh_dual_simple(&Z2);
    assert!(matches!(
        have_finite_resolution_in_add_m(&module, &module, &Integer::from(-1)),
        Err(EnvelopingModuleError::EmptyResolutionDifferentialList { bound })
            if bound == -1
    ));

    let first = fresh_a2_zero(&Z2);
    let second = fresh_a2_zero(&Z2);
    assert!(matches!(
        have_finite_resolution_in_add_m(&first, &second, &Integer::from(0)),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
}

fn cotilting_resolution_shape<F: Field>(resolution: &FiniteAddMResolutionResult<F>) -> Value {
    let FiniteAddMResolutionResult::Resolution(resolution) = resolution else {
        return Value::Bool(true);
    };
    let shape = |degree| {
        let differential = resolution.differential(&Integer::from(degree)).unwrap();
        json!({
            "source": differential.source().dimension_vector(),
            "range": differential.range().dimension_vector(),
        })
    };
    json!({
        "0": shape(0),
        "1": shape(1),
        "2": shape(2),
        "3": shape(3),
    })
}

fn cotilting_value_json<F: Field>(result: Option<CotiltingModuleResult<F>>) -> Value {
    let Some(result) = result else {
        return Value::Bool(false);
    };
    json!({
        "injective_dimension": result.injective_dimension(),
        "resolution_count": result.resolutions().len(),
        "resolutions": result
            .resolutions()
            .iter()
            .map(cotilting_resolution_shape)
            .collect::<Vec<_>>(),
    })
}

fn modulehomalg_cotilting_case<F: QpaCoefficientField>(field: &F) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let injectives = indec_injective_modules(&a2).unwrap();
    let cogenerator = direct_sum_of_qpa_modules(&injectives).unwrap().unwrap();
    let simples = simple_modules(&a2).unwrap();
    let id_one = direct_sum_of_qpa_modules(&[injectives[1].clone(), simples[1].clone()])
        .unwrap()
        .unwrap();
    let semisimple = direct_sum_of_qpa_modules(&simples).unwrap().unwrap();

    let dual = fresh_dual_numbers(field);
    let dual_regular = direct_sum_of_qpa_modules(&indec_injective_modules(&dual).unwrap())
        .unwrap()
        .unwrap();
    let dual_simple = simple_modules(&dual).unwrap()[0].clone();

    let check = |module: &PathAlgebraModule<F>, bound| {
        let has_before = is_cotilting_module(module).is_some();
        let cotilting = cotilting_module(module, &Integer::from(bound)).unwrap();
        let has_after = is_cotilting_module(module).is_some();
        let mut value = json!({
            "has_before": has_before,
            "value": cotilting_value_json(cotilting),
            "has_after": has_after,
        });
        if let Some(cached) = is_cotilting_module(module) {
            value
                .as_object_mut()
                .unwrap()
                .insert("cached_value".to_owned(), Value::Bool(cached));
        }
        value
    };
    json!({
        "cogenerator_bound_0": check(&cogenerator, 0),
        "cogenerator_cached_negative": check(&cogenerator, -1),
        "id_one_bound_0": check(&id_one, 0),
        "id_one_bound_1": check(&id_one, 1),
        "semisimple_bound_1": check(&semisimple, 1),
        "dual_regular_bound_0": check(&dual_regular, 0),
        "dual_simple_bound_1": check(&dual_simple, 1),
    })
}

#[test]
fn pinned_qpa_modulehomalg_cotilting_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_cotilting.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_cotilting",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_cotilting_case(&Z2),
        "quinary": modulehomalg_cotilting_case(&Zp::new(5)),
        "rational": modulehomalg_cotilting_case(&Q),
    });
    assert_eq!(actual, expected);
}

#[test]
fn cotilting_preserves_qpa_unbound_injective_dimension_failure() {
    let module = fresh_a2_module(&Z2, vec![1, 0]);
    let dual = dual_of_module(&module).unwrap();
    let _ = dual.cache_projective_dimension(0);
    assert!(matches!(
        cotilting_module(&module, &Integer::from(1)),
        Err(EnvelopingModuleError::MissingCotiltingInjectiveDimensionCache)
    ));
}

fn complement_resolution_json<F: Field>(resolution: Option<&AddMResolution<F>>) -> Value {
    let Some(resolution) = resolution else {
        return Value::Array(Vec::new());
    };
    let shape = |degree| {
        let differential = resolution.differential(&Integer::from(degree)).unwrap();
        json!({
            "source": differential.source().dimension_vector(),
            "range": differential.range().dimension_vector(),
        })
    };
    json!({
        "-3": shape(-3),
        "-2": shape(-2),
        "-1": shape(-1),
        "0": shape(0),
        "1": shape(1),
        "2": shape(2),
        "3": shape(3),
        "4": shape(4),
    })
}

fn complement_coresolution_json<F: Field>(coresolution: Option<&AddMCoresolution<F>>) -> Value {
    let Some(coresolution) = coresolution else {
        return Value::Array(Vec::new());
    };
    let shape = |degree| {
        let differential = coresolution.differential(&Integer::from(degree)).unwrap();
        json!({
            "source": differential.source().dimension_vector(),
            "range": differential.range().dimension_vector(),
        })
    };
    json!({
        "-3": shape(-3),
        "-2": shape(-2),
        "-1": shape(-1),
        "0": shape(0),
        "1": shape(1),
        "2": shape(2),
        "3": shape(3),
        "4": shape(4),
    })
}

fn all_complements_json<F: Field>(result: &AllComplementsResult<F>) -> Value {
    json!({
        "resolution": complement_resolution_json(result.resolution()),
        "coresolution": complement_coresolution_json(result.coresolution()),
    })
}

fn modulehomalg_all_complements_case<F: QpaCoefficientField>(field: &F) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let projectives = indec_projective_modules(&a2).unwrap();
    let simples = simple_modules(&a2).unwrap();
    let zero = zero_module(&a2).unwrap();

    let dual = fresh_dual_numbers(field);
    let dual_regular = indec_projective_modules(&dual).unwrap()[0].clone();
    let dual_zero = zero_module(&dual).unwrap();

    let evaluate = |module, complement| {
        all_complements_json(
            &all_complements_of_almost_complete_tilting_module(module, complement).unwrap(),
        )
    };
    json!({
        "p1_complement_p2": evaluate(&projectives[0], &projectives[1]),
        "p2_complement_p1": evaluate(&projectives[1], &projectives[0]),
        "p1_complement_s1": evaluate(&projectives[0], &simples[0]),
        "zero_complement_simple": evaluate(&zero, &simples[0]),
        "dual_zero_complement_regular": evaluate(&dual_zero, &dual_regular),
    })
}

#[test]
fn pinned_qpa_modulehomalg_all_complements_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_all_complements.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_all_complements",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_all_complements_case(&Z2),
        "quinary": modulehomalg_all_complements_case(&Zp::new(5)),
        "rational": modulehomalg_all_complements_case(&Q),
    });
    assert_eq!(actual, expected);
}

#[test]
fn all_complements_synonym_and_foreign_algebra_follow_qpa() {
    let algebra = isolated_path_algebra(Z2, 1);
    let zero = zero_module(&algebra).unwrap();
    let simple = right_module_over_path_algebra(algebra, vec![1], &[])
        .unwrap()
        .unwrap();
    let tilting = all_complements_of_almost_complete_tilting_module(&zero, &simple).unwrap();
    let cotilting = all_complements_of_almost_complete_cotilting_module(&zero, &simple).unwrap();
    assert_eq!(
        all_complements_json(&tilting),
        all_complements_json(&cotilting)
    );

    let first = fresh_a2_zero(&Z2);
    let second = fresh_a2_zero(&Z2);
    assert!(matches!(
        all_complements_of_almost_complete_tilting_module(&first, &second),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
}

fn faithful_dimension_json(result: FaithfulDimensionResult) -> Value {
    match result {
        FaithfulDimensionResult::Finite(value) => value.into(),
        FaithfulDimensionResult::Infinite => "infinity".into(),
    }
}

fn modulehomalg_faithful_dimension_case<F: QpaCoefficientField>(field: &F) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let projectives = indec_projective_modules(&a2).unwrap();
    let regular = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
    let simples = simple_modules(&a2).unwrap();
    let pd_one = direct_sum_of_qpa_modules(&[projectives[0].clone(), simples[0].clone()])
        .unwrap()
        .unwrap();
    let zero = zero_module(&a2).unwrap();

    let dual = fresh_dual_numbers(field);
    let dual_regular = indec_projective_modules(&dual).unwrap()[0].clone();
    let dual_simple = simple_modules(&dual).unwrap()[0].clone();

    let check = |module: &PathAlgebraModule<F>| {
        let has_before = module.cached_faithful_dimension().is_some();
        let first = faithful_dimension(module).unwrap();
        let has_after_first = module.cached_faithful_dimension().is_some();
        let second = faithful_dimension(module).unwrap();
        let has_after_second = module.cached_faithful_dimension().is_some();
        json!({
            "has_before": has_before,
            "first": faithful_dimension_json(first),
            "has_after_first": has_after_first,
            "second": faithful_dimension_json(second),
            "has_after_second": has_after_second,
        })
    };
    json!({
        "regular": check(&regular),
        "projective_1": check(&projectives[0]),
        "simple_1": check(&simples[0]),
        "simple_2": check(&simples[1]),
        "pd_one": check(&pd_one),
        "zero": check(&zero),
        "dual_regular": check(&dual_regular),
        "dual_simple": check(&dual_simple),
    })
}

#[test]
fn pinned_qpa_modulehomalg_faithful_dimension_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_faithful_dimension.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_faithful_dimension",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_faithful_dimension_case(&Z2),
        "quinary": modulehomalg_faithful_dimension_case(&Zp::new(5)),
        "rational": modulehomalg_faithful_dimension_case(&Q),
    });
    assert_eq!(actual, expected);
}

fn number_of_complements_json(result: NumberOfComplementsResult) -> Value {
    match result {
        NumberOfComplementsResult::Finite(value) => value.into(),
        NumberOfComplementsResult::Infinite => "infinity".into(),
    }
}

fn modulehomalg_number_complements_case<F: QpaCoefficientField>(field: &F) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let projectives = indec_projective_modules(&a2).unwrap();
    let regular = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
    let simples = simple_modules(&a2).unwrap();
    let zero = zero_module(&a2).unwrap();

    let check = |module: &PathAlgebraModule<F>| {
        let has_before = module.cached_faithful_dimension().is_some();
        let tilting = number_of_complements_of_almost_complete_tilting_module(module).unwrap();
        let has_after = module.cached_faithful_dimension().is_some();
        let cotilting = number_of_complements_of_almost_complete_cotilting_module(module).unwrap();
        json!({
            "has_before": has_before,
            "tilting": number_of_complements_json(tilting),
            "has_after": has_after,
            "cotilting": number_of_complements_json(cotilting),
        })
    };
    json!({
        "regular": check(&regular),
        "projective_1": check(&projectives[0]),
        "simple_1": check(&simples[0]),
        "zero": check(&zero),
    })
}

#[test]
fn pinned_qpa_modulehomalg_number_complements_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_number_complements.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_number_complements",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_number_complements_case(&Z2),
        "quinary": modulehomalg_number_complements_case(&Zp::new(5)),
        "rational": modulehomalg_number_complements_case(&Q),
    });
    assert_eq!(actual, expected);
}

fn complement_mutation_json<F: Field>(module: Option<PathAlgebraModule<F>>) -> Value {
    module.map_or(
        Value::Bool(false),
        |module| json!({"dimension": module.dimension_vector()}),
    )
}

fn modulehomalg_complement_mutations_case<F: QpaCoefficientField>(field: &F) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let projectives = indec_projective_modules(&a2).unwrap();
    let simples = simple_modules(&a2).unwrap();
    let zero = zero_module(&a2).unwrap();

    let dual = fresh_dual_numbers(field);
    let dual_regular = indec_projective_modules(&dual).unwrap()[0].clone();
    let dual_simple = simple_modules(&dual).unwrap()[0].clone();
    let dual_zero = zero_module(&dual).unwrap();

    let check = |module: &PathAlgebraModule<F>, complement: &PathAlgebraModule<F>| {
        json!({
            "left": complement_mutation_json(
                left_mutation_of_tilting_module_complement(module, complement).unwrap()
            ),
            "right": complement_mutation_json(
                right_mutation_of_tilting_module_complement(module, complement).unwrap()
            ),
            "left_cotilting": complement_mutation_json(
                left_mutation_of_cotilting_module_complement(module, complement).unwrap()
            ),
            "right_cotilting": complement_mutation_json(
                right_mutation_of_cotilting_module_complement(module, complement).unwrap()
            ),
        })
    };
    json!({
        "p1_complement_p2": check(&projectives[0], &projectives[1]),
        "p1_complement_s1": check(&projectives[0], &simples[0]),
        "p2_complement_p1": check(&projectives[1], &projectives[0]),
        "zero_complement_s1": check(&zero, &simples[0]),
        "dual_regular_complement_simple": check(&dual_regular, &dual_simple),
        "dual_zero_complement_regular": check(&dual_zero, &dual_regular),
    })
}

#[test]
fn pinned_qpa_modulehomalg_complement_mutations_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_complement_mutations.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_complement_mutations",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_complement_mutations_case(&Z2),
        "quinary": modulehomalg_complement_mutations_case(&Zp::new(5)),
        "rational": modulehomalg_complement_mutations_case(&Q),
    });
    assert_eq!(actual, expected);
}

fn cotilting_approximation_map_json<F: Field>(
    map: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    json!({
        "source_dimension": map.source().dimension_vector(),
        "range_dimension": map.range().dimension_vector(),
        "matrices": map
            .vertex_matrices()
            .iter()
            .map(|matrix| (0..matrix.nrows())
                .map(|row| matrix.row(row).iter().map(coefficient).collect::<Vec<_>>())
                .collect::<Vec<_>>())
            .collect::<Vec<_>>(),
    })
}

fn modulehomalg_cotilting_approximations_case<F: QpaCoefficientField>(
    field: &F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let injectives = indec_injective_modules(&a2).unwrap();
    let simples = simple_modules(&a2).unwrap();
    let projectives = indec_projective_modules(&a2).unwrap();
    let cotilting = direct_sum_of_qpa_modules(&[injectives[1].clone(), simples[1].clone()])
        .unwrap()
        .unwrap();
    let zero = zero_module(&a2).unwrap();
    cotilting_module(&cotilting, &Integer::from(1))
        .unwrap()
        .expect("the pinned A2 module is cotilting");

    let check = |module: &PathAlgebraModule<F>| {
        json!({
            "right": cotilting_approximation_map_json(
                &right_approximation_by_perp_t(&cotilting, module).unwrap(),
                coefficient,
            ),
            "left": cotilting_approximation_map_json(
                &left_approximation_by_add_t_hat(&cotilting, module).unwrap(),
                coefficient,
            ),
        })
    };
    json!({
        "zero": check(&zero),
        "simple_one": check(&simples[0]),
        "simple_two": check(&simples[1]),
        "projective_one": check(&projectives[0]),
        "projective_two": check(&projectives[1]),
    })
}

#[test]
fn pinned_qpa_modulehomalg_cotilting_approximations_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_cotilting_approximations.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_cotilting_approximations",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_cotilting_approximations_case(&Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": modulehomalg_cotilting_approximations_case(&field5, |entry| json!(
            field5_coefficients
                .to_integer(entry)
                .to_string()
                .parse::<i64>()
                .unwrap()
        )),
        "rational": modulehomalg_cotilting_approximations_case(&Q, |entry| json!([
            entry.numerator_ref().to_string().parse::<i64>().unwrap(),
            entry.denominator_ref().to_string().parse::<i64>().unwrap(),
        ])),
    });
    assert_eq!(actual, expected);
}

#[test]
fn cotilting_approximations_preserve_attribute_check_order() {
    let unbound_algebra = isolated_path_algebra(Z2, 1);
    let unbound = zero_module(&unbound_algebra).unwrap();
    assert!(matches!(
        right_approximation_by_perp_t(&unbound, &unbound),
        Err(EnvelopingModuleError::CotiltingApproximationPropertyUnbound)
    ));
    assert!(matches!(
        left_approximation_by_add_t_hat(&unbound, &unbound),
        Err(EnvelopingModuleError::CotiltingApproximationPropertyUnbound)
    ));

    let false_algebra = isolated_path_algebra(Z2, 1);
    let false_module = zero_module(&false_algebra).unwrap();
    let _ = false_module.cache_is_cotilting_module(false);
    assert!(matches!(
        right_approximation_by_perp_t(&false_module, &false_module),
        Err(EnvelopingModuleError::FirstApproximationArgumentNotCotilting)
    ));
    assert!(matches!(
        left_approximation_by_add_t_hat(&false_module, &false_module),
        Err(EnvelopingModuleError::FirstApproximationArgumentNotCotilting)
    ));

    let true_algebra = isolated_path_algebra(Z2, 1);
    let true_without_dimension = zero_module(&true_algebra).unwrap();
    let _ = true_without_dimension.cache_is_cotilting_module(true);
    assert!(
        right_approximation_by_perp_t(&true_without_dimension, &true_without_dimension)
            .unwrap()
            .is_zero()
    );
    assert!(
        left_approximation_by_add_t_hat(&true_without_dimension, &true_without_dimension)
            .unwrap()
            .is_zero()
    );

    let nonzero = simple_modules(&true_algebra).unwrap()[0].clone();
    assert!(matches!(
        right_approximation_by_perp_t(&true_without_dimension, &nonzero),
        Err(EnvelopingModuleError::MissingCotiltingApproximationInjectiveDimensionCache)
    ));
    assert!(matches!(
        left_approximation_by_add_t_hat(&true_without_dimension, &nonzero),
        Err(EnvelopingModuleError::MissingCotiltingApproximationInjectiveDimensionCache)
    ));
}

#[test]
fn cotilting_approximations_preserve_dimension_zero_index_failure() {
    let algebra = isolated_path_algebra(Z2, 1);
    let cotilting = simple_modules(&algebra).unwrap()[0].clone();
    cotilting_module(&cotilting, &Integer::from(0))
        .unwrap()
        .expect("the isolated simple is injective cotilting");
    assert!(matches!(
        right_approximation_by_perp_t(&cotilting, &cotilting),
        Err(EnvelopingModuleError::EmptyCotiltingApproximationExactSequences)
    ));
    assert!(matches!(
        left_approximation_by_add_t_hat(&cotilting, &cotilting),
        Err(EnvelopingModuleError::EmptyCotiltingApproximationExactSequences)
    ));
}

fn modulehomalg_is_nth_syzygy_case<F: QpaCoefficientField>(field: &F) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simples = simple_modules(&a2).unwrap();
    let projectives = indec_projective_modules(&a2).unwrap();
    let zero = zero_module(&a2).unwrap();
    let dual = fresh_dual_numbers(field);
    let dual_simple = simple_modules(&dual).unwrap()[0].clone();
    let dual_regular = indec_projective_modules(&dual).unwrap()[0].clone();
    json!({
        "zero_0": is_nth_syzygy(&zero, &Integer::from(0)).unwrap(),
        "zero_1": is_nth_syzygy(&zero, &Integer::from(1)).unwrap(),
        "simple_one_0": is_nth_syzygy(&simples[0], &Integer::from(0)).unwrap(),
        "simple_one_1": is_nth_syzygy(&simples[0], &Integer::from(1)).unwrap(),
        "simple_one_2": is_nth_syzygy(&simples[0], &Integer::from(2)).unwrap(),
        "simple_two_1": is_nth_syzygy(&simples[1], &Integer::from(1)).unwrap(),
        "projective_one_1": is_nth_syzygy(&projectives[0], &Integer::from(1)).unwrap(),
        "projective_two_2": is_nth_syzygy(&projectives[1], &Integer::from(2)).unwrap(),
        "dual_simple_1": is_nth_syzygy(&dual_simple, &Integer::from(1)).unwrap(),
        "dual_simple_2": is_nth_syzygy(&dual_simple, &Integer::from(2)).unwrap(),
        "dual_regular_1": is_nth_syzygy(&dual_regular, &Integer::from(1)).unwrap(),
    })
}

#[test]
fn pinned_qpa_modulehomalg_is_nth_syzygy_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_is_nth_syzygy.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_is_nth_syzygy",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_is_nth_syzygy_case(&Z2),
        "quinary": modulehomalg_is_nth_syzygy_case(&Zp::new(5)),
        "rational": modulehomalg_is_nth_syzygy_case(&Q),
    });
    assert_eq!(actual, expected);
}

#[test]
fn is_nth_syzygy_propagates_negative_degree_failure() {
    let module = fresh_a2_module(&Z2, vec![1, 0]);
    assert!(matches!(
        is_nth_syzygy(&module, &Integer::from(-1)),
        Err(EnvelopingModuleError::NegativeSyzygyDegree { .. })
    ));
}

fn modulehomalg_list_approximation_case<F: QpaCoefficientField>(
    field: &F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simples = simple_modules(&a2).unwrap();
    let projectives = indec_projective_modules(&a2).unwrap();
    let zero = zero_module(&a2).unwrap();
    let dual = fresh_dual_numbers(field);
    let dual_simple = simple_modules(&dual).unwrap()[0].clone();
    let dual_regular = indec_projective_modules(&dual).unwrap()[0].clone();
    let evaluate = |modules: &[PathAlgebraModule<F>], target: &PathAlgebraModule<F>| {
        cotilting_approximation_map_json(
            &right_approximation_by_add_m_list(modules, target).unwrap(),
            coefficient,
        )
    };
    json!({
        "empty_to_simple_one": evaluate(&[], &simples[0]),
        "p1_p2_to_simple_one": evaluate(&projectives, &simples[0]),
        "s1_p1_to_simple_one": evaluate(
            &[simples[0].clone(), projectives[0].clone()],
            &simples[0],
        ),
        "p1_s1_to_simple_one": evaluate(
            &[projectives[0].clone(), simples[0].clone()],
            &simples[0],
        ),
        "repeated_s1_to_simple_one": evaluate(
            &[simples[0].clone(), simples[0].clone()],
            &simples[0],
        ),
        "modules_to_zero": evaluate(&projectives, &zero),
        "dual_simple_regular_to_simple": evaluate(
            &[dual_simple.clone(), dual_regular.clone()],
            &dual_simple,
        ),
        "dual_regular_simple_to_simple": evaluate(
            &[dual_regular.clone(), dual_simple.clone()],
            &dual_simple,
        ),
        "dual_regular_to_simple": evaluate(&[dual_regular], &dual_simple),
    })
}

#[test]
fn pinned_qpa_modulehomalg_list_approximation_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_list_approximation.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_list_approximation",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_list_approximation_case(&Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": modulehomalg_list_approximation_case(&field5, |entry| json!(
            field5_coefficients
                .to_integer(entry)
                .to_string()
                .parse::<i64>()
                .unwrap()
        )),
        "rational": modulehomalg_list_approximation_case(&Q, |entry| json!([
            entry.numerator_ref().to_string().parse::<i64>().unwrap(),
            entry.denominator_ref().to_string().parse::<i64>().unwrap(),
        ])),
    });
    assert_eq!(actual, expected);
}

#[test]
fn list_approximation_preserves_qpa_shortcut_and_empty_selection_order() {
    let first_algebra = isolated_path_algebra(Z2, 2);
    let second_algebra = isolated_path_algebra(Z2, 2);
    let target = simple_modules(&first_algebra).unwrap()[0].clone();
    let zero_target = zero_module(&first_algebra).unwrap();
    let foreign = simple_modules(&second_algebra).unwrap()[0].clone();
    assert!(
        right_approximation_by_add_m_list(std::slice::from_ref(&foreign), &zero_target)
            .unwrap()
            .is_zero()
    );
    assert!(
        right_approximation_by_add_m_list(&[], &target)
            .unwrap()
            .is_zero()
    );
    assert!(matches!(
        right_approximation_by_add_m_list(&[foreign], &target),
        Err(EnvelopingModuleError::ListApproximationModulesOverDifferentAlgebras)
    ));

    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            Z2,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simples = simple_modules(&a2).unwrap();
    assert!(matches!(
        right_approximation_by_add_m_list(&[simples[1].clone()], &simples[0]),
        Err(EnvelopingModuleError::EmptyListApproximationDirectSum)
    ));
}

fn modulehomalg_radical_list_approximation_case<F: QpaCoefficientField>(
    field: &F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simples = simple_modules(&a2).unwrap();
    let projectives = indec_projective_modules(&a2).unwrap();
    let dual = fresh_dual_numbers(field);
    let dual_simple = simple_modules(&dual).unwrap()[0].clone();
    let dual_regular = indec_projective_modules(&dual).unwrap()[0].clone();
    let evaluate = |modules: &[PathAlgebraModule<F>], index| {
        cotilting_approximation_map_json(
            &radical_right_approximation_by_add_m(modules, &Integer::from(index)).unwrap(),
            coefficient,
        )
    };
    json!({
        "projectives_at_one": evaluate(&projectives, 1),
        "s1_p1_at_one": evaluate(
            &[simples[0].clone(), projectives[0].clone()],
            1,
        ),
        "p1_s1_at_two": evaluate(
            &[projectives[0].clone(), simples[0].clone()],
            2,
        ),
        "single_s1": evaluate(std::slice::from_ref(&simples[0]), 1),
        "repeated_s1": evaluate(&[simples[0].clone(), simples[0].clone()], 1),
        "single_dual_simple": evaluate(std::slice::from_ref(&dual_simple), 1),
        "single_dual_regular": evaluate(std::slice::from_ref(&dual_regular), 1),
        "dual_regular_simple_at_one": evaluate(
            &[dual_regular.clone(), dual_simple.clone()],
            1,
        ),
        "dual_regular_simple_at_two": evaluate(
            &[dual_regular.clone(), dual_simple.clone()],
            2,
        ),
        "repeated_dual_regular": evaluate(
            &[dual_regular.clone(), dual_regular.clone()],
            1,
        ),
    })
}

#[test]
fn pinned_qpa_modulehomalg_radical_list_approximation_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_radical_list_approximation.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_radical_list_approximation",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_radical_list_approximation_case(&Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": modulehomalg_radical_list_approximation_case(&field5, |entry| json!(
            field5_coefficients
                .to_integer(entry)
                .to_string()
                .parse::<i64>()
                .unwrap()
        )),
        "rational": modulehomalg_radical_list_approximation_case(&Q, |entry| json!([
            entry.numerator_ref().to_string().parse::<i64>().unwrap(),
            entry.denominator_ref().to_string().parse::<i64>().unwrap(),
        ])),
    });
    assert_eq!(actual, expected);
}

fn simple_endo_resolution_json<F: Field>(
    result: &ProjectiveResolutionOfSimpleModuleOverEndoResult<F>,
) -> Value {
    let projective_dimension = match result.projective_dimension() {
        ProjectiveResolutionOfSimpleModuleOverEndoStatus::Finite(dimension) => {
            json!(dimension.to_string().parse::<i64>().unwrap())
        }
        ProjectiveResolutionOfSimpleModuleOverEndoStatus::GreaterThanModuleIndex(index) => {
            Value::String(format!("projdim > {index}"))
        }
    };
    json!({
        "projective_dimension": projective_dimension,
        "syzygy_dimension_vectors": result
            .syzygies()
            .iter()
            .map(PathAlgebraModule::dimension_vector)
            .collect::<Vec<_>>(),
    })
}

fn modulehomalg_projective_resolution_simple_endo_case<F: QpaCoefficientField>(field: &F) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simples = simple_modules(&a2).unwrap();
    let projectives = indec_projective_modules(&a2).unwrap();
    let dual = fresh_dual_numbers(field);
    let dual_simple = simple_modules(&dual).unwrap()[0].clone();
    let dual_regular = indec_projective_modules(&dual).unwrap()[0].clone();
    let evaluate = |modules: &[PathAlgebraModule<F>], index, length| {
        simple_endo_resolution_json(
            &projective_resolution_of_simple_module_over_endo(
                modules,
                &Integer::from(index),
                &Integer::from(length),
            )
            .unwrap(),
        )
    };
    json!({
        "single_s1_length_3": evaluate(std::slice::from_ref(&simples[0]), 1, 3),
        "projectives_at_one_length_3": evaluate(&projectives, 1, 3),
        "dual_regular_simple_at_two_length_3": evaluate(
            &[dual_regular.clone(), dual_simple.clone()],
            2,
            3,
        ),
        "dual_regular_simple_at_two_length_0": evaluate(
            &[dual_regular.clone(), dual_simple],
            2,
            0,
        ),
        "single_dual_regular_length_0": evaluate(
            std::slice::from_ref(&dual_regular),
            1,
            0,
        ),
        "single_dual_regular_length_2": evaluate(
            std::slice::from_ref(&dual_regular),
            1,
            2,
        ),
        "single_dual_regular_length_3": evaluate(
            std::slice::from_ref(&dual_regular),
            1,
            3,
        ),
        "single_dual_regular_length_4": evaluate(
            std::slice::from_ref(&dual_regular),
            1,
            4,
        ),
    })
}

#[test]
fn pinned_qpa_modulehomalg_projective_resolution_simple_endo_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_projective_resolution_simple_endo.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_projective_resolution_simple_endo",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_projective_resolution_simple_endo_case(&Z2),
        "quinary": modulehomalg_projective_resolution_simple_endo_case(&Zp::new(5)),
        "rational": modulehomalg_projective_resolution_simple_endo_case(&Q),
    });
    assert_eq!(actual, expected);
}

#[test]
fn simple_endo_resolution_preserves_negative_length_check_order() {
    assert!(matches!(
        projective_resolution_of_simple_module_over_endo::<Z2>(
            &[],
            &Integer::from(99),
            &Integer::from(-1),
        ),
        Err(EnvelopingModuleError::NegativeSimpleEndoProjectiveResolutionLength { .. })
    ));
    assert!(matches!(
        projective_resolution_of_simple_module_over_endo::<Z2>(
            &[],
            &Integer::from(99),
            &Integer::from(0),
        ),
        Err(EnvelopingModuleError::RadicalApproximationIndexOutOfRange { .. })
    ));
}

#[test]
fn radical_list_approximation_preserves_index_and_inner_failure_order() {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            Z2,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let projectives = indec_projective_modules(&a2).unwrap();
    for index in [-1, 0, 3] {
        assert!(matches!(
            radical_right_approximation_by_add_m(&projectives, &Integer::from(index)),
            Err(EnvelopingModuleError::RadicalApproximationIndexOutOfRange { .. })
        ));
    }
    assert!(matches!(
        radical_right_approximation_by_add_m::<Z2>(&[], &Integer::from(1)),
        Err(EnvelopingModuleError::RadicalApproximationIndexOutOfRange { .. })
    ));
    assert!(matches!(
        radical_right_approximation_by_add_m(&projectives, &Integer::from(2)),
        Err(EnvelopingModuleError::EmptyListApproximationDirectSum)
    ));
}

#[test]
fn complement_mutations_reject_foreign_algebras() {
    let first = fresh_a2_zero(&Z2);
    let second = fresh_a2_zero(&Z2);
    assert!(matches!(
        left_mutation_of_tilting_module_complement(&first, &second),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
    assert!(matches!(
        right_mutation_of_tilting_module_complement(&first, &second),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
}

#[test]
fn add_m_approximations_reject_foreign_algebras_before_zero_shortcuts() {
    let first = isolated_path_algebra(Z2, 1);
    let second = isolated_path_algebra(Z2, 1);
    let first_zero = zero_module(&first).unwrap();
    let second_zero = zero_module(&second).unwrap();
    assert!(matches!(
        right_approximation_by_add_m(&first_zero, &second_zero),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
    assert!(matches!(
        left_approximation_by_add_m(&first_zero, &second_zero),
        Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
    ));
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_double_star_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let zero = field.zero();
    let one = field.one();
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let a2_zero = zero_module(&a2).unwrap();
    let a2_simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let a2_simple_two = right_module_over_path_algebra(a2.clone(), vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let a2_projective = right_module_over_path_algebra(
        a2.clone(),
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let a2_semisimple_pair = right_module_over_path_algebra(
        a2,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual.clone(),
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let dual_regular = right_module_over_path_algebra(
        dual,
        vec![2],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(
                vec![vec![zero.clone(), one.clone()], vec![zero.clone(), zero]],
                field,
            )
            .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let evaluate = |module: &PathAlgebraModule<F>| {
        double_star_map_json(
            &from_identity_to_double_star_homomorphism(module).unwrap(),
            coefficient,
        )
    };
    json!({
        "a2_zero": evaluate(&a2_zero),
        "a2_simple_one": evaluate(&a2_simple_one),
        "a2_simple_two": evaluate(&a2_simple_two),
        "a2_projective": evaluate(&a2_projective),
        "a2_semisimple_pair": evaluate(&a2_semisimple_pair),
        "dual_simple": evaluate(&dual_simple),
        "dual_regular": evaluate(&dual_regular),
    })
}

#[test]
fn pinned_qpa_modulehom_double_star_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_double_star.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehom_double_star",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehom_double_star_case(Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": modulehom_double_star_case(field5, |entry| json!(
            field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "rational": modulehom_double_star_case(Q, |entry| json!([
            entry.numerator_ref().to_string().parse::<i64>().unwrap(),
            entry.denominator_ref().to_string().parse::<i64>().unwrap(),
        ])),
    });

    assert_eq!(actual, expected);
}

#[test]
fn identity_to_double_star_homomorphism_is_cached_with_star_modules() {
    let algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let module = right_module_over_path_algebra(
        algebra,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![Q.one()]], Q).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let first = from_identity_to_double_star_homomorphism(&module).unwrap();
    let second = from_identity_to_double_star_homomorphism(&module).unwrap();
    let double_star = star_of_module(&star_of_module(&module).unwrap()).unwrap();

    assert!(first.is_identical(&second));
    assert!(first.source().is_identical(&module));
    assert!(first.range().is_identical(&double_star));
}

fn almost_split_sequence_json<F: Field>(
    sequence: Option<PathAlgebraModuleHomomorphismPair<F>>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    match sequence {
        None => Value::Null,
        Some((monomorphism, epimorphism)) => json!({
            "monomorphism": double_star_map_json(&monomorphism, coefficient),
            "epimorphism": double_star_map_json(&epimorphism, coefficient),
        }),
    }
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn module_ar_almost_split_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let zero = field.zero();
    let one = field.one();
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let simple_two = right_module_over_path_algebra(a2.clone(), vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let projective = right_module_over_path_algebra(
        a2,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual,
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![zero]], field).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    json!({
        "ending_simple_one": almost_split_sequence_json(
            almost_split_sequence(&simple_one).unwrap(), coefficient),
        "ending_simple_two": almost_split_sequence_json(
            almost_split_sequence(&simple_two).unwrap(), coefficient),
        "ending_projective": almost_split_sequence_json(
            almost_split_sequence(&projective).unwrap(), coefficient),
        "ending_dual_simple": almost_split_sequence_json(
            almost_split_sequence(&dual_simple).unwrap(), coefficient),
        "starting_simple_one": almost_split_sequence_json(
            almost_split_sequence_from_side(&simple_one, AlmostSplitSide::Starting).unwrap(),
            coefficient),
        "starting_simple_two": almost_split_sequence_json(
            almost_split_sequence_from_side(&simple_two, AlmostSplitSide::Starting).unwrap(),
            coefficient),
        "starting_dual_simple": almost_split_sequence_json(
            almost_split_sequence_from_side(&dual_simple, AlmostSplitSide::Starting).unwrap(),
            coefficient),
    })
}

#[test]
fn pinned_qpa_module_ar_almost_split_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/module_ar_almost_split.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "module_ar_almost_split",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": module_ar_almost_split_case(Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": module_ar_almost_split_case(field5, |entry| json!(
            field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "rational": module_ar_almost_split_case(Q, |entry| json!([
            entry.numerator_ref().to_string().parse::<i64>().unwrap(),
            entry.denominator_ref().to_string().parse::<i64>().unwrap(),
        ])),
    });

    assert_eq!(actual, expected);
}

#[allow(clippy::needless_pass_by_value)]
fn module_ar_perp_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let injectives = indec_injective_modules(&a2).unwrap();
    let simples = simple_modules(&a2).unwrap();
    let cotilting = direct_sum_of_qpa_modules(&[injectives[1].clone(), simples[1].clone()])
        .unwrap()
        .unwrap();
    cotilting_module(&cotilting, &Integer::one())
        .unwrap()
        .expect("the pinned A2 module is cotilting");

    json!({
        "simple_one": almost_split_sequence_json(
            almost_split_sequence_in_perp_t(&cotilting, &simples[0]).unwrap(), coefficient
        ),
        "projective_simple_two": almost_split_sequence_json(
            almost_split_sequence_in_perp_t(&cotilting, &simples[1]).unwrap(), coefficient
        ),
    })
}

#[test]
fn pinned_qpa_module_ar_perp_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/module_ar_perp.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let rational = module_ar_perp_case(Q, |entry: &Rational| {
        json!(entry.numerator_ref().to_i64().unwrap())
    });
    let binary = module_ar_perp_case(Z2, |entry| json!(Z2.to_integer(entry).to_i64().unwrap()));
    let quinary = module_ar_perp_case(field5, |entry| {
        json!(field5_coefficients.to_integer(entry).to_i64().unwrap())
    });

    assert_eq!(rational, expected["rational"]);
    assert_eq!(binary, expected["binary"]);
    assert_eq!(quinary, expected["quinary"]);
}

#[test]
fn almost_split_sequence_in_perp_t_preserves_unbound_property_error() {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simples = simple_modules(&a2).unwrap();
    assert!(matches!(
        almost_split_sequence_in_perp_t(&simples[1], &simples[0]),
        Err(AuslanderReitenError::CotiltingPerpPropertyUnbound)
    ));
}

#[test]
fn almost_split_sequence_caches_values_and_projective_fail() {
    let algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Z2,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let nonprojective = simple_modules(&algebra).unwrap()[0].clone();
    let projective = simple_modules(&algebra).unwrap()[1].clone();
    let first = almost_split_sequence(&nonprojective).unwrap().unwrap();
    let second = almost_split_sequence(&nonprojective).unwrap().unwrap();

    assert!(first.0.is_identical(&second.0));
    assert!(first.1.is_identical(&second.1));
    assert!(almost_split_sequence(&projective).unwrap().is_none());
    assert!(matches!(
        projective.cached_almost_split_sequence(),
        CachedAttribute::Fail
    ));
    assert!(almost_split_sequence(&projective).unwrap().is_none());
}

fn irreducible_maps_json<F: Field>(
    maps: &[PathAlgebraModuleHomomorphism<F>],
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    json!(
        maps.iter()
            .map(|map| double_star_map_json(map, coefficient))
            .collect::<Vec<_>>()
    )
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn module_ar_irreducible_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let zero = field.zero();
    let one = field.one();
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let simple_two = right_module_over_path_algebra(a2.clone(), vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let projective = right_module_over_path_algebra(
        a2,
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual,
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![zero]], field).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    json!({
        "ending_simple_one": irreducible_maps_json(
            &irreducible_morphisms_ending_in(&simple_one).unwrap(), coefficient),
        "ending_projective": irreducible_maps_json(
            &irreducible_morphisms_ending_in(&projective).unwrap(), coefficient),
        "ending_dual_simple": irreducible_maps_json(
            &irreducible_morphisms_ending_in(&dual_simple).unwrap(), coefficient),
        "starting_simple_two": irreducible_maps_json(
            &irreducible_morphisms_starting_in(&simple_two).unwrap(), coefficient),
        "starting_dual_simple": irreducible_maps_json(
            &irreducible_morphisms_starting_in(&dual_simple).unwrap(), coefficient),
    })
}

#[test]
fn pinned_qpa_module_ar_irreducible_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/module_ar_irreducible.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "module_ar_irreducible",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": module_ar_irreducible_case(Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": module_ar_irreducible_case(field5, |entry| json!(
            field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
    });

    assert_eq!(actual, expected);
}

#[test]
fn irreducible_morphisms_preserve_qpa_field_error_and_live_cache() {
    let binary = QuiverAlgebra::Path(
        PathAlgebra::new(
            Z2,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simples = simple_modules(&binary).unwrap();
    let first = irreducible_morphisms_ending_in(&simples[0]).unwrap();
    let second = irreducible_morphisms_ending_in(&simples[0]).unwrap();
    assert_eq!(first.len(), second.len());
    assert!(
        first
            .iter()
            .zip(&second)
            .all(|(left, right)| left.is_identical(right))
    );
    assert!(matches!(
        irreducible_morphisms_ending_in(&simples[1]),
        Err(super::AuslanderReitenError::ZeroIrreducibleSourceDecomposition)
    ));

    let rational = QuiverAlgebra::Path(
        PathAlgebra::new(
            Q,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let rational_simple = simple_modules(&rational).unwrap()[0].clone();
    assert!(matches!(
        irreducible_morphisms_ending_in(&rational_simple),
        Err(super::AuslanderReitenError::FiniteFieldRequired)
    ));
    assert!(matches!(
        irreducible_morphisms_starting_in(&rational_simple),
        Err(super::AuslanderReitenError::FiniteFieldRequired)
    ));
}

fn predecessor_module_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    json!({
        "dimension": module.dimension_vector(),
        "matrices": module.arrow_matrices().iter().map(|matrix| json!(
            (0..matrix.nrows()).map(|row| matrix.row(row).iter()
                .map(coefficient).collect::<Vec<_>>()).collect::<Vec<_>>()
        )).collect::<Vec<_>>(),
    })
}

fn predecessors_json<F: Field>(
    result: PredecessorsOfModuleResult<F>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    match result {
        PredecessorsOfModuleResult::Module(module) => json!({
            "kind": "module",
            "module": predecessor_module_json(&module, coefficient),
        }),
        PredecessorsOfModuleResult::Layers(result) => json!({
            "kind": "layers",
            "layers": result.layers.iter().map(|layer| layer.iter()
                .map(|module| predecessor_module_json(module, coefficient))
                .collect::<Vec<_>>()).collect::<Vec<_>>(),
            "valuations": result.valuations.iter().map(|layer| layer.iter()
                .map(|valuation| {
                    let multiplicity = match valuation.multiplicity {
                        PredecessorMultiplicity::Pair { first, second } => json!([
                            first.map_or(Value::Bool(false), Value::from),
                            second.map_or(Value::Bool(false), Value::from),
                        ]),
                        PredecessorMultiplicity::Scalar(value) => json!(value),
                    };
                    json!([
                        valuation.predecessor_position,
                        valuation.successor_position,
                        multiplicity,
                    ])
                }).collect::<Vec<_>>()).collect::<Vec<_>>(),
        }),
    }
}

#[allow(clippy::needless_pass_by_value)]
fn module_ar_predecessors_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simple_one = right_module_over_path_algebra(a2, vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual,
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![field.zero()]], field).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    json!({
        "a2_n1": predecessors_json(
            predecessors_of_module(&simple_one, &Integer::from(1)).unwrap(), coefficient),
        "a2_n2": predecessors_json(
            predecessors_of_module(&simple_one, &Integer::from(2)).unwrap(), coefficient),
        "a2_n3": predecessors_json(
            predecessors_of_module(&simple_one, &Integer::from(3)).unwrap(), coefficient),
        "a2_n4": predecessors_json(
            predecessors_of_module(&simple_one, &Integer::from(4)).unwrap(), coefficient),
        "dual_n2": predecessors_json(
            predecessors_of_module(&dual_simple, &Integer::from(2)).unwrap(), coefficient),
        "dual_n3": predecessors_json(
            predecessors_of_module(&dual_simple, &Integer::from(3)).unwrap(), coefficient),
        "dual_n4": predecessors_json(
            predecessors_of_module(&dual_simple, &Integer::from(4)).unwrap(), coefficient),
    })
}

#[test]
fn pinned_qpa_module_ar_predecessors_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/module_ar_predecessors.schema-v1.json"
    ))
    .unwrap();
    let field5 = Zp::new(5);
    let field5_coefficients = field5.clone();
    let actual = json!({
        "schema_version": 1,
        "kind": "module_ar_predecessors",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": module_ar_predecessors_case(Z2, |entry| json!(
            Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
        "quinary": module_ar_predecessors_case(field5, |entry| json!(
            field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
        )),
    });

    assert_eq!(actual, expected);
}

#[test]
fn predecessors_preserve_projective_check_and_n_one_return() {
    let algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Z2,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simples = simple_modules(&algebra).unwrap();
    assert!(matches!(
        predecessors_of_module(&simples[1], &Integer::from(1)),
        Err(super::AuslanderReitenError::PredecessorModuleProjective)
    ));
    assert!(matches!(
        predecessors_of_module(&simples[0], &Integer::from(0)),
        Err(super::AuslanderReitenError::PredecessorBoundBelowOne { .. })
    ));
    let PredecessorsOfModuleResult::Module(returned) =
        predecessors_of_module(&simples[0], &Integer::from(1)).unwrap()
    else {
        panic!("QPA must return the module itself at n = 1");
    };
    assert!(returned.is_identical(&simples[0]));
}

fn module_map_json<F: Field>(
    map: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!({
        "source_dimension": map.source().dimension_vector(),
        "range_dimension": map.range().dimension_vector(),
        "matrices": map
            .vertex_matrices()
            .iter()
            .map(|matrix| matrix_json(matrix, coefficient))
            .collect::<Vec<_>>(),
    })
}

fn projective_module_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!({
        "dimension": module.dimension_vector(),
        "matrices": module
            .arrow_matrices()
            .iter()
            .map(|matrix| matrix_json(matrix, coefficient))
            .collect::<Vec<_>>(),
    })
}

fn module_elements_json<F: Field>(
    elements: &[PathAlgebraModuleElement<F>],
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!(
        elements
            .iter()
            .map(|element| element
                .external_representation()
                .iter()
                .map(|row| row.iter().map(coefficient).collect::<Vec<_>>())
                .collect::<Vec<_>>())
            .collect::<Vec<_>>()
    )
}

fn path_element_json<F: Field>(
    element: &qpa_core::PathAlgebraElement<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!(
        element
            .external_representation()
            .terms()
            .iter()
            .map(|(walk, scalar)| json!({
                "walk": walk,
                "coefficient": coefficient(scalar),
            }))
            .collect::<Vec<_>>()
    )
}

fn path_element_matrix_json<F: Field>(
    matrix: &[Vec<qpa_core::PathAlgebraElement<F>>],
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!(
        matrix
            .iter()
            .map(|row| row
                .iter()
                .map(|element| path_element_json(element, coefficient))
                .collect::<Vec<_>>())
            .collect::<Vec<_>>()
    )
}

fn physical_module_element_json<F: Field>(
    element: &PathAlgebraModuleElement<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!({
        "dimension": element.module().dimension_vector(),
        "support": element.support_vertices(),
        "coordinates": element
            .flat_external_representation()
            .iter()
            .map(coefficient)
            .collect::<Vec<_>>(),
    })
}

fn entered_path_matrix_equals<F: Field>(
    returned: &[Vec<qpa_core::PathAlgebraElement<F>>],
    entered: &[Vec<QuiverAlgebraElement<F>>],
) -> bool {
    returned.len() == entered.len()
        && returned
            .iter()
            .zip(entered)
            .all(|(returned_row, entered_row)| {
                returned_row.len() == entered_row.len()
                    && returned_row.iter().zip(entered_row).all(
                        |(returned_element, entered_element)| {
                            matches!(entered_element, QuiverAlgebraElement::Path(element)
                                if returned_element == element)
                        },
                    )
            })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_functor_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let zero = field.zero();
    let one = field.one();
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let module = right_module_over_path_algebra(
        algebra.clone(),
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one.clone()], vec![zero.clone()]], field.clone())
                .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let map = right_module_hom_over_algebra(
        module.clone(),
        module.clone(),
        vec![
            DenseMatrix::from_rows(
                vec![
                    vec![one.clone(), zero.clone()],
                    vec![zero.clone(), zero.clone()],
                ],
                field.clone(),
            )
            .unwrap(),
            DenseMatrix::identity(1, field.clone()),
        ],
    )
    .unwrap();
    let minimal = minimal_generating_set_of_module(module.clone()).unwrap();
    let top_projection = top_of_module_projection(module.clone()).unwrap();
    let top = top_of_module(module.clone()).unwrap();
    let socle_inclusion = socle_of_module_inclusion(module.clone()).unwrap();
    let dual_map = dual_of_module_homomorphism(&map).unwrap();
    let target = right_module_over_path_algebra(
        algebra.clone(),
        vec![1, 2],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one.clone(), one]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let simple_two = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let zero_module = right_module_over_path_algebra(algebra, vec![0, 0], &[])
        .unwrap()
        .unwrap();
    let zero_socle = socle_of_module_inclusion(zero_module.clone()).unwrap();
    let zero_trace = trace_of_module(&simple_one, &target).unwrap();
    let zero_source_trace = trace_of_module(&zero_module, &target).unwrap();
    let zero_target_trace = trace_of_module(&simple_two, &zero_module).unwrap();
    let nonzero_trace = trace_of_module(&simple_two, &target).unwrap();
    let legacy_trace = right_fac_approximation_via_all_basis_images(&simple_two, &target).unwrap();
    let ordered_trace = trace_of_module(&module, &module).unwrap();
    let ordered_legacy = right_fac_approximation_via_all_basis_images(&module, &module).unwrap();

    json!({
        "minimal_generators": {
            "count": minimal.len(),
            "elements": module_elements_json(&minimal, coefficient),
            "zero_count": minimal_generating_set_of_module(zero_module.clone()).unwrap().len(),
            "top_projection_cached_identity": top_projection.is_identical(
                &top_of_module_projection(module.clone()).unwrap()),
            "top_cached_identity": top.is_identical(&top_of_module(module).unwrap()),
            "top_projection_range_identity": top_projection.range().is_identical(&top),
        },
        "socle": {
            "inclusion": module_map_json(&socle_inclusion, coefficient),
            "object_dimension": socle_of_module(socle_inclusion.range().clone())
                .unwrap()
                .dimension_vector(),
            "zero_inclusion": module_map_json(&zero_socle, coefficient),
            "zero_source_identity": zero_socle.source().is_identical(&zero_module),
        },
        "dual_map": {
            "cached_identity": dual_map.is_identical(
                &dual_of_module_homomorphism(&map).unwrap()),
            "source_dimension": dual_map.source().dimension_vector(),
            "range_dimension": dual_map.range().dimension_vector(),
            "matrices": dual_map.vertex_matrices().iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
        },
        "trace": {
            "zero": module_map_json(&zero_trace, coefficient),
            "zero_source": module_map_json(&zero_source_trace, coefficient),
            "zero_target": module_map_json(&zero_target_trace, coefficient),
            "nonzero": module_map_json(&nonzero_trace, coefficient),
            "legacy": module_map_json(&legacy_trace, coefficient),
            "legacy_equal": legacy_trace == nonzero_trace,
            "ordered": module_map_json(&ordered_trace, coefficient),
            "ordered_legacy": module_map_json(&ordered_legacy, coefficient),
            "ordered_legacy_equal": ordered_legacy == ordered_trace,
            "right_fac_equal": right_fac_approximation(&simple_two, &target).unwrap()
                == nonzero_trace,
        },
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_projective_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let zero = field.zero();
    let one = field.one();
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let projectives = indec_projective_modules(&algebra).unwrap();
    let bases = basis_of_projectives(&algebra).unwrap().unwrap();
    let cached_projectives = indec_projective_modules(&algebra).unwrap();
    let basis_paths = bases
        .iter()
        .map(|projective| {
            projective
                .iter()
                .map(|vertex_basis| {
                    vertex_basis
                        .iter()
                        .map(|element| element.leading_monomial().to_string())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let module = right_module_over_path_algebra(
        algebra.clone(),
        vec![2, 2],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(
                vec![
                    vec![one.clone(), zero.clone()],
                    vec![zero.clone(), zero.clone()],
                ],
                field.clone(),
            )
            .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let generators = minimal_generating_set_of_module(module.clone()).unwrap();
    let cover = projective_cover(&module).unwrap();
    let from_generator = hom_from_projective(&generators[0], &module)
        .unwrap()
        .unwrap();
    let projective_sum = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
    let lifting = lifting_morphism_from_projective(&cover, &cover)
        .unwrap()
        .unwrap();
    let nondirect_lifting = lifting_morphism_from_projective(&cover, &from_generator).unwrap();
    let zero_module = right_module_over_path_algebra(algebra.clone(), vec![0, 0], &[])
        .unwrap()
        .unwrap();
    let zero_to_module = zero_module_mapping(zero_module.clone(), module.clone()).unwrap();
    let zero_lifting = lifting_morphism_from_projective(&cover, &zero_to_module)
        .unwrap()
        .unwrap();
    let factoring = hom_factoring_through_projectives(&module, &module).unwrap();
    let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let simple_two = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let layer = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
        .unwrap()
        .unwrap();
    let first = right_module_hom_over_algebra(
        layer.clone(),
        module.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone()).unwrap(),
            DenseMatrix::from_rows(vec![vec![one.clone(), zero.clone()]], field.clone()).unwrap(),
        ],
    )
    .unwrap();
    let second = right_module_hom_over_algebra(
        layer,
        module.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone()).unwrap(),
            DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone()).unwrap(),
        ],
    )
    .unwrap();
    let kernel_module = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let kernel_inclusion = right_module_hom_over_algebra(
        kernel_module,
        module.clone(),
        vec![
            DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone()).unwrap(),
            DenseMatrix::zero(1, 2, field.clone()),
        ],
    )
    .unwrap();
    let inclusion_lift = lifting_inclusion_morphisms(&first, &kernel_inclusion)
        .unwrap()
        .unwrap();
    let inclusion_fail = lifting_inclusion_morphisms(&first, &second).unwrap();

    let loop_quiver =
        Quiver::from_named(vec!["w".into()], vec![NamedArrowSpec::named("w", "w", "x")]).unwrap();
    let path_algebra = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let loop_element = path_algebra
        .element(path_algebra.quiver().arrows().next().unwrap().path())
        .unwrap();
    let dual_numbers = QuiverAlgebra::Quotient(
        quotient_by_relations(
            path_algebra,
            &[loop_element.multiply(&loop_element).unwrap()],
        )
        .unwrap(),
    );
    let dual_projective = indec_projective_modules(&dual_numbers).unwrap()[0].clone();
    let dual_basis = basis_of_projectives(&dual_numbers).unwrap().unwrap();
    let dual_basis_paths = dual_basis
        .iter()
        .map(|projective| {
            projective
                .iter()
                .map(|vertex_basis| {
                    vertex_basis
                        .iter()
                        .map(|element| element.leading_monomial().to_string())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let dual_simple = right_module_over_path_algebra(dual_numbers, vec![1], &[])
        .unwrap()
        .unwrap();
    let dual_cover = projective_cover(&dual_simple).unwrap();

    json!({
        "a2": {
            "basis_paths": basis_paths,
            "projectives": projectives
                .iter()
                .map(|projective| projective_module_json(projective, coefficient))
                .collect::<Vec<_>>(),
            "cache_identity": projectives
                .iter()
                .zip(&cached_projectives)
                .map(|(projective, cached)| projective.is_identical(cached))
                .collect::<Vec<_>>(),
            "installed_properties": [
                projectives[0].cached_is_projective_module() == Some(true),
                projectives[1].cached_is_projective_module() == Some(true),
                projectives[0].cached_is_indecomposable_module() == Some(true),
                projectives[1].cached_is_indecomposable_module() == Some(true),
            ],
            "projectivity": [
                is_projective_module(&projective_sum).unwrap(),
                is_projective_module(&simple_one).unwrap(),
                is_projective_module(&simple_two).unwrap(),
                is_projective_module(&zero_module).unwrap(),
            ],
            "minimal_generator_count": generators.len(),
            "from_generator": module_map_json(&from_generator, coefficient),
            "cover": module_map_json(&cover, coefficient),
            "cover_cache_identity": cover.is_identical(&projective_cover(&module).unwrap()),
            "lifting": module_map_json(&lifting, coefficient),
            "lifting_triangle": lifting.compose(&cover).unwrap() == cover,
            "nondirect_lifting_is_fail": nondirect_lifting.is_none(),
            "zero_lifting": module_map_json(&zero_lifting, coefficient),
            "factoring_count": factoring.len(),
            "factoring": factoring
                .iter()
                .map(|map| module_map_json(map, coefficient))
                .collect::<Vec<_>>(),
            "factoring_zero_counts": [
                hom_factoring_through_projectives(&zero_module, &module).unwrap().len(),
                hom_factoring_through_projectives(&module, &zero_module).unwrap().len(),
            ],
            "inclusion_lift": module_map_json(&inclusion_lift, coefficient),
            "inclusion_triangle": inclusion_lift.compose(&first).unwrap() == kernel_inclusion,
            "inclusion_fail": inclusion_fail.is_none(),
        },
        "dual_numbers": {
            "basis_paths": dual_basis_paths,
            "projective": projective_module_json(&dual_projective, coefficient),
            "generator_count": minimal_generating_set_of_module(dual_projective.clone())
                .unwrap().len(),
            "is_projective": is_projective_module(&dual_projective).unwrap(),
            "simple_is_projective": is_projective_module(&dual_simple).unwrap(),
            "cover": module_map_json(&dual_cover, coefficient),
        },
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_projective_matrix_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let path_algebra = PathAlgebra::new(field.clone(), quiver).unwrap();
    let algebra = QuiverAlgebra::Path(path_algebra);
    let vertices = algebra.quiver().vertices().collect::<Vec<_>>();
    let arrow = algebra.quiver().arrows().next().unwrap();
    let v1 = include_in_path_algebra(&vertices[0].path(), &algebra).unwrap();
    let v2 = include_in_path_algebra(&vertices[1].path(), &algebra).unwrap();
    let a = include_in_path_algebra(&arrow.path(), &algebra).unwrap();
    let entered_matrix = vec![vec![v1.clone(), a.clone()], vec![algebra.zero(), v2]];
    let map =
        from_matrix_to_homomorphism_of_projectives(&algebra, &entered_matrix, &[1, 2], &[1, 2])
            .unwrap();
    let returned_matrix = matrix_of_homomorphism_between_projectives(&map).unwrap();
    let uniform = v1.add(&a).unwrap();
    let projective_element = element_in_v_a_as_element_in_indec_proj(&algebra, &uniform).unwrap();
    let recovered = element_in_indec_projective(&algebra, &projective_element, 1).unwrap();
    let projectives = indec_projective_modules(&algebra).unwrap();
    let sum = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
    let zero_module = right_module_over_path_algebra(algebra.clone(), vec![0, 0], &[])
        .unwrap()
        .unwrap();
    let zero_to_sum = matrix_of_homomorphism_between_projectives(
        &zero_module_mapping(zero_module.clone(), sum.clone()).unwrap(),
    )
    .unwrap();
    let sum_to_zero = matrix_of_homomorphism_between_projectives(
        &zero_module_mapping(sum, zero_module.clone()).unwrap(),
    )
    .unwrap();
    let zero_to_zero = matrix_of_homomorphism_between_projectives(
        &zero_module_mapping(zero_module.clone(), zero_module).unwrap(),
    )
    .unwrap();
    let single =
        from_matrix_to_homomorphism_of_projectives(&algebra, &[vec![a]], &[2], &[1]).unwrap();

    let quotient_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "x"),
            ArrowSpec::named(vertex(2), vertex(3), "y"),
        ],
    )
    .unwrap();
    let quotient_parent = PathAlgebra::new(field, quotient_quiver).unwrap();
    let quotient_arrows = quotient_parent.quiver().arrows().collect::<Vec<_>>();
    let x = quotient_parent.element(quotient_arrows[0].path()).unwrap();
    let y = quotient_parent.element(quotient_arrows[1].path()).unwrap();
    let quotient = QuiverAlgebra::Quotient(
        quotient_by_relations(quotient_parent.clone(), &[x.multiply(&y).unwrap()]).unwrap(),
    );
    let quotient_vertices = quotient.quiver().vertices().collect::<Vec<_>>();
    let qv1 = include_in_path_algebra(&quotient_vertices[0].path(), &quotient).unwrap();
    let qv2 = include_in_path_algebra(&quotient_vertices[1].path(), &quotient).unwrap();
    let qx = include_in_path_algebra(&quotient_arrows[0].path(), &quotient).unwrap();
    let qy = include_in_path_algebra(&quotient_arrows[1].path(), &quotient).unwrap();
    let quotient_matrix = vec![vec![qx.clone(), quotient.zero()], vec![qv2, qy]];
    let quotient_map =
        from_matrix_to_homomorphism_of_projectives(&quotient, &quotient_matrix, &[2, 3], &[1, 2])
            .unwrap();
    let quotient_back = matrix_of_homomorphism_between_projectives(&quotient_map).unwrap();
    let quotient_uniform = qv1.add(&qx).unwrap();
    let quotient_projective_element =
        element_in_v_a_as_element_in_indec_proj(&quotient, &quotient_uniform).unwrap();
    let quotient_recovered =
        element_in_indec_projective(&quotient, &quotient_projective_element, 1).unwrap();

    json!({
        "path": {
            "map": module_map_json(&map, coefficient),
            "matrix": path_element_matrix_json(&returned_matrix, coefficient),
            "roundtrip_equal": entered_path_matrix_equals(
                &returned_matrix,
                &entered_matrix,
            ),
            "underlying": matrix_json(&underlying_linear_map(&map).unwrap(), coefficient),
            "uniform_projective_element": physical_module_element_json(
                &projective_element,
                coefficient,
            ),
            "uniform_recovered": path_element_json(&recovered, coefficient),
            "single": module_map_json(&single, coefficient),
            "single_underlying": matrix_json(
                &underlying_linear_map(&single).unwrap(),
                coefficient,
            ),
            "zero_to_sum": path_element_matrix_json(&zero_to_sum, coefficient),
            "sum_to_zero": path_element_matrix_json(&sum_to_zero, coefficient),
            "zero_to_zero": path_element_matrix_json(&zero_to_zero, coefficient),
        },
        "quotient": {
            "map": module_map_json(&quotient_map, coefficient),
            "matrix_parent_representatives": path_element_matrix_json(
                &quotient_back,
                coefficient,
            ),
            "underlying": matrix_json(
                &underlying_linear_map(&quotient_map).unwrap(),
                coefficient,
            ),
            "uniform_projective_element": physical_module_element_json(
                &quotient_projective_element,
                coefficient,
            ),
            "uniform_recovered_parent": path_element_json(
                &quotient_recovered,
                coefficient,
            ),
        },
    })
}

fn attribute_module_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!({
        "dimension": module.dimension_vector(),
        "matrices": module
            .matrices_of_path_algebra_module()
            .iter()
            .map(|entry| match entry {
                PathAlgebraModuleMatrix::ZeroScalar => json!("0"),
                PathAlgebraModuleMatrix::Matrix(matrix) => matrix_json(matrix, coefficient),
            })
            .collect::<Vec<_>>(),
    })
}

fn attribute_module_list_json<F: Field>(
    modules: &[PathAlgebraModule<F>],
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Vec<Value> {
    modules
        .iter()
        .map(|module| attribute_module_json(module, coefficient))
        .collect()
}

fn star_module_attribute_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = module.cached_star_of_module().is_some();
    let value = star_of_module(module).unwrap();
    let has_after = module.cached_star_of_module().is_some();
    let repeated = star_of_module(module).unwrap();
    json!({
        "has_before": has_before,
        "value": attribute_module_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
    })
}

fn nakayama_module_attribute_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = module.cached_nakayama_functor_of_module().is_some();
    let value = nakayama_functor_of_module(module).unwrap();
    let has_after = module.cached_nakayama_functor_of_module().is_some();
    let repeated = nakayama_functor_of_module(module).unwrap();
    json!({
        "has_before": has_before,
        "value": attribute_module_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
    })
}

fn opposite_nakayama_module_attribute_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = module
        .cached_opposite_nakayama_functor_of_module()
        .is_some();
    let value = opposite_nakayama_functor_of_module(module).unwrap();
    let has_after = module
        .cached_opposite_nakayama_functor_of_module()
        .is_some();
    let repeated = opposite_nakayama_functor_of_module(module).unwrap();
    json!({
        "has_before": has_before,
        "value": attribute_module_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
    })
}

fn module_functor_attribute_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let star = star_module_attribute_json(module, coefficient);
    let nakayama = nakayama_module_attribute_json(module, coefficient);
    let opposite_nakayama = opposite_nakayama_module_attribute_json(module, coefficient);
    json!({
        "star": star,
        "nakayama": nakayama,
        "opposite_nakayama": opposite_nakayama,
    })
}

fn star_map_attribute_json<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = homomorphism.cached_star_homomorphism().is_some();
    let value = star_of_module_homomorphism(homomorphism).unwrap();
    let has_after = homomorphism.cached_star_homomorphism().is_some();
    let repeated = star_of_module_homomorphism(homomorphism).unwrap();
    json!({
        "has_before": has_before,
        "value": module_map_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
    })
}

fn nakayama_map_attribute_json<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = homomorphism
        .cached_nakayama_functor_homomorphism()
        .is_some();
    let value = nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
    let has_after = homomorphism
        .cached_nakayama_functor_homomorphism()
        .is_some();
    let repeated = nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
    json!({
        "has_before": has_before,
        "value": module_map_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
    })
}

fn opposite_nakayama_map_attribute_json<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = homomorphism
        .cached_opposite_nakayama_functor_homomorphism()
        .is_some();
    let value = opposite_nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
    let has_after = homomorphism
        .cached_opposite_nakayama_functor_homomorphism()
        .is_some();
    let repeated = opposite_nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
    json!({
        "has_before": has_before,
        "value": module_map_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
    })
}

fn map_functor_attribute_json<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let star = star_map_attribute_json(homomorphism, coefficient);
    let nakayama = nakayama_map_attribute_json(homomorphism, coefficient);
    let opposite_nakayama = opposite_nakayama_map_attribute_json(homomorphism, coefficient);
    json!({
        "star": star,
        "nakayama": nakayama,
        "opposite_nakayama": opposite_nakayama,
    })
}

fn transpose_attribute_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = module.cached_transpose_of_module().is_some();
    let value = transpose_of_module(module).unwrap();
    let has_after = module.cached_transpose_of_module().is_some();
    let repeated = transpose_of_module(module).unwrap();
    json!({
        "has_before": has_before,
        "value": attribute_module_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
        "indecomposable_bound": value.cached_is_indecomposable_module().is_some(),
        "indecomposable_value": value.cached_is_indecomposable_module(),
    })
}

fn d_tr_attribute_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = module.cached_d_tr().is_some();
    let value = d_tr(module).unwrap();
    let has_after = module.cached_d_tr().is_some();
    let repeated = d_tr(module).unwrap();
    json!({
        "has_before": has_before,
        "value": attribute_module_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
    })
}

fn tr_d_attribute_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = module.cached_tr_d().is_some();
    let value = tr_d(module).unwrap();
    let has_after = module.cached_tr_d().is_some();
    let repeated = tr_d(module).unwrap();
    json!({
        "has_before": has_before,
        "value": attribute_module_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
    })
}

fn transpose_functor_attribute_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let transpose = transpose_attribute_json(module, coefficient);
    let dtr = d_tr_attribute_json(module, coefficient);
    let trd = tr_d_attribute_json(module, coefficient);
    let tau_has_before = module.cached_is_tau_rigid_module().is_some();
    let tau_value = is_tau_rigid_module(module).unwrap();
    let tau_has_after = module.cached_is_tau_rigid_module().is_some();
    json!({
        "transpose": transpose,
        "dtr": dtr,
        "trd": trd,
        "tau_rigid": {
            "has_before": tau_has_before,
            "value": tau_value,
            "has_after": tau_has_after,
        },
    })
}

fn functor_iterations_json<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let d_zero = d_tr_iterated(module, 0).unwrap();
    let d_one = d_tr_iterated(module, 1).unwrap();
    let d_one_generators = minimal_generating_set_of_module(d_one.clone()).unwrap();
    let d_minus_one = d_tr_iterated(module, -1).unwrap();
    let d_two = d_tr_iterated(module, 2).unwrap();
    let t_zero = tr_d_iterated(module, 0).unwrap();
    let t_one = tr_d_iterated(module, 1).unwrap();
    let t_minus_one = tr_d_iterated(module, -1).unwrap();
    let t_two = tr_d_iterated(module, 2).unwrap();
    json!({
        "dtr_zero_identity": d_zero.is_identical(module),
        "dtr_one_identity": d_one.is_identical(&d_tr(module).unwrap()),
        "dtr_one_generators": d_one_generators.iter().map(|element|
            element.external_representation().iter().map(|row|
                row.iter().map(coefficient).collect::<Vec<_>>()).collect::<Vec<_>>())
            .collect::<Vec<_>>(),
        "dtr_minus_one_identity": d_minus_one.is_identical(&tr_d(module).unwrap()),
        "dtr_two": attribute_module_json(&d_two, coefficient),
        "trd_zero_identity": t_zero.is_identical(module),
        "trd_one_identity": t_one.is_identical(&tr_d(module).unwrap()),
        "trd_minus_one_identity": t_minus_one.is_identical(&d_tr(module).unwrap()),
        "trd_two": attribute_module_json(&t_two, coefficient),
    })
}

fn transpose_map_attribute_json<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = homomorphism.cached_transpose_homomorphism().is_some();
    let value = transpose_of_module_homomorphism(homomorphism).unwrap();
    let has_after = homomorphism.cached_transpose_homomorphism().is_some();
    let repeated = transpose_of_module_homomorphism(homomorphism).unwrap();
    let source_attribute = transpose_of_module(homomorphism.range()).unwrap();
    let range_attribute = transpose_of_module(homomorphism.source()).unwrap();
    json!({
        "has_before": has_before,
        "value": module_map_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
        "source_attribute_identity": value.source().is_identical(&source_attribute),
        "range_attribute_identity": value.range().is_identical(&range_attribute),
    })
}

fn d_tr_map_attribute_json<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = homomorphism.cached_d_tr_homomorphism().is_some();
    let value = d_tr_homomorphism(homomorphism).unwrap();
    let has_after = homomorphism.cached_d_tr_homomorphism().is_some();
    let repeated = d_tr_homomorphism(homomorphism).unwrap();
    let source_attribute = d_tr(homomorphism.source()).unwrap();
    let range_attribute = d_tr(homomorphism.range()).unwrap();
    json!({
        "has_before": has_before,
        "value": module_map_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
        "source_attribute_identity": value.source().is_identical(&source_attribute),
        "range_attribute_identity": value.range().is_identical(&range_attribute),
    })
}

fn tr_d_map_attribute_json<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let has_before = homomorphism.cached_tr_d_homomorphism().is_some();
    let value = tr_d_homomorphism(homomorphism).unwrap();
    let has_after = homomorphism.cached_tr_d_homomorphism().is_some();
    let repeated = tr_d_homomorphism(homomorphism).unwrap();
    let source_attribute = tr_d(homomorphism.source()).unwrap();
    let range_attribute = tr_d(homomorphism.range()).unwrap();
    json!({
        "has_before": has_before,
        "value": module_map_json(&value, coefficient),
        "has_after": has_after,
        "repeat_identity": value.is_identical(&repeated),
        "source_attribute_identity": value.source().is_identical(&source_attribute),
        "range_attribute_identity": value.range().is_identical(&range_attribute),
    })
}

fn ar_map_functor_attribute_json<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!({
        "entered": module_map_json(homomorphism, coefficient),
        "transpose": transpose_map_attribute_json(homomorphism, coefficient),
        "dtr": d_tr_map_attribute_json(homomorphism, coefficient),
        "trd": tr_d_map_attribute_json(homomorphism, coefficient),
    })
}

fn ar_map_iterations_json<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let d_zero = d_tr_homomorphism_iterated(homomorphism, 0).unwrap();
    let d_one = d_tr_homomorphism_iterated(homomorphism, 1).unwrap();
    let d_minus_one = d_tr_homomorphism_iterated(homomorphism, -1).unwrap();
    let d_two = d_tr_homomorphism_iterated(homomorphism, 2).unwrap();
    let t_zero = tr_d_homomorphism_iterated(homomorphism, 0).unwrap();
    let t_one = tr_d_homomorphism_iterated(homomorphism, 1).unwrap();
    let t_minus_one = tr_d_homomorphism_iterated(homomorphism, -1).unwrap();
    let t_two = tr_d_homomorphism_iterated(homomorphism, 2).unwrap();
    json!({
        "dtr_zero_identity": d_zero.is_identical(homomorphism),
        "dtr_one_identity": d_one.is_identical(&d_tr_homomorphism(homomorphism).unwrap()),
        "dtr_minus_one_identity": d_minus_one
            .is_identical(&tr_d_homomorphism(homomorphism).unwrap()),
        "dtr_two": module_map_json(&d_two, coefficient),
        "trd_zero_identity": t_zero.is_identical(homomorphism),
        "trd_one_identity": t_one.is_identical(&tr_d_homomorphism(homomorphism).unwrap()),
        "trd_minus_one_identity": t_minus_one
            .is_identical(&d_tr_homomorphism(homomorphism).unwrap()),
        "trd_two": module_map_json(&t_two, coefficient),
    })
}

fn tensor_space_json<F: Field>(
    tensor: &TensorProductOfModulesResult<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!({
        "dimension": tensor.dimension(),
        "basis": tensor
            .basis_vectors()
            .iter()
            .map(|row| row.iter().map(coefficient).collect::<Vec<_>>())
            .collect::<Vec<_>>(),
    })
}

fn caught_tensor_element_json<F: Field>(
    tensor: &TensorProductOfModulesResult<F>,
    left: &PathAlgebraModuleElement<F>,
    right: &PathAlgebraModuleElement<F>,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    match tensor.elementary_tensor(left, right) {
        Ok(value) => json!({
            "success": true,
            "value": value.iter().map(coefficient).collect::<Vec<_>>(),
        }),
        Err(_) => json!({"success": false, "value": Value::Null}),
    }
}

fn property_bound<F: Field>(module: &PathAlgebraModule<F>) -> Vec<bool> {
    vec![
        module.cached_is_injective_module().is_some(),
        module.cached_is_simple_qpa_module().is_some(),
        module.cached_is_semisimple_module().is_some(),
    ]
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn module_families_properties_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let injectives = indec_injective_modules(&algebra).unwrap();
    let injectives_again = indec_injective_modules(&algebra).unwrap();
    let simples = simple_modules(&algebra).unwrap();
    let simples_again = simple_modules(&algebra).unwrap();
    let projectives = indec_projective_modules(&algebra).unwrap();
    let zero = zero_module(&algebra).unwrap();
    let semisimple_sum = direct_sum_of_qpa_modules(&simples).unwrap().unwrap();
    let injective_installed = injectives
        .iter()
        .map(|module| {
            vec![
                module.cached_is_injective_module() == Some(true),
                module.cached_is_indecomposable_module() == Some(true),
                module.cached_is_projective_module().is_some(),
            ]
        })
        .collect::<Vec<_>>();
    let simple_installed = simples
        .iter()
        .map(|module| {
            vec![
                module.cached_is_simple_qpa_module() == Some(true),
                module.cached_is_indecomposable_module() == Some(true),
                module.cached_is_semisimple_module() == Some(true),
            ]
        })
        .collect::<Vec<_>>();
    let injective_values = injectives
        .iter()
        .chain(&simples)
        .chain(&projectives)
        .chain(std::iter::once(&zero))
        .map(|module| is_injective_module(module).unwrap())
        .collect::<Vec<_>>();
    let simple_values = simples
        .iter()
        .chain(&projectives)
        .chain([&zero, &semisimple_sum])
        .map(is_simple_qpa_module)
        .collect::<Vec<_>>();
    let semisimple_values = simples
        .iter()
        .chain(&projectives)
        .chain([&zero, &semisimple_sum])
        .map(|module| is_semisimple_module(module).unwrap())
        .collect::<Vec<_>>();
    let properties_bound = injectives
        .iter()
        .chain(&simples)
        .chain(&projectives)
        .chain([&zero, &semisimple_sum])
        .map(property_bound)
        .collect::<Vec<_>>();

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let parent = PathAlgebra::new(field, loop_quiver).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let quotient =
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap());
    let quotient_injectives = indec_injective_modules(&quotient).unwrap();
    let quotient_injectives_again = indec_injective_modules(&quotient).unwrap();
    let quotient_simples = simple_modules(&quotient).unwrap();
    let quotient_simples_again = simple_modules(&quotient).unwrap();
    let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
    let quotient_zero = zero_module(&quotient).unwrap();
    let quotient_injective_installed = quotient_injectives
        .iter()
        .map(|module| {
            vec![
                module.cached_is_injective_module() == Some(true),
                module.cached_is_indecomposable_module() == Some(true),
                module.cached_is_projective_module().is_some(),
            ]
        })
        .collect::<Vec<_>>();
    let quotient_simple_installed = quotient_simples
        .iter()
        .map(|module| {
            vec![
                module.cached_is_simple_qpa_module() == Some(true),
                module.cached_is_indecomposable_module() == Some(true),
                module.cached_is_semisimple_module() == Some(true),
            ]
        })
        .collect::<Vec<_>>();
    let quotient_injective_values = quotient_injectives
        .iter()
        .chain(&quotient_simples)
        .chain([&quotient_projective, &quotient_zero])
        .map(|module| is_injective_module(module).unwrap())
        .collect::<Vec<_>>();
    let quotient_simple_values = quotient_simples
        .iter()
        .chain([&quotient_projective, &quotient_zero])
        .map(is_simple_qpa_module)
        .collect::<Vec<_>>();
    let quotient_semisimple_values = quotient_simples
        .iter()
        .chain([&quotient_projective, &quotient_zero])
        .map(|module| is_semisimple_module(module).unwrap())
        .collect::<Vec<_>>();
    let quotient_properties_bound = quotient_injectives
        .iter()
        .chain(&quotient_simples)
        .chain([&quotient_projective, &quotient_zero])
        .map(property_bound)
        .collect::<Vec<_>>();

    json!({
        "path": {
            "injectives": attribute_module_list_json(&injectives, coefficient),
            "injective_cache_identity": injectives.iter().zip(&injectives_again)
                .map(|(first, second)| first.is_identical(second)).collect::<Vec<_>>(),
            "injective_installed": injective_installed,
            "simples": attribute_module_list_json(&simples, coefficient),
            "simple_cache_identity": simples.iter().zip(&simples_again)
                .map(|(first, second)| first.is_identical(second)).collect::<Vec<_>>(),
            "simple_installed": simple_installed,
            "injective_values": injective_values,
            "simple_values": simple_values,
            "semisimple_values": semisimple_values,
            "property_bound": properties_bound,
        },
        "quotient": {
            "injectives": attribute_module_list_json(&quotient_injectives, coefficient),
            "injective_cache_identity": quotient_injectives.iter()
                .zip(&quotient_injectives_again).map(|(first, second)|
                    first.is_identical(second)).collect::<Vec<_>>(),
            "injective_installed": quotient_injective_installed,
            "simples": attribute_module_list_json(&quotient_simples, coefficient),
            "simple_cache_identity": quotient_simples.iter().zip(&quotient_simples_again)
                .map(|(first, second)| first.is_identical(second)).collect::<Vec<_>>(),
            "simple_installed": quotient_simple_installed,
            "injective_values": quotient_injective_values,
            "simple_values": quotient_simple_values,
            "semisimple_values": quotient_semisimple_values,
            "property_bound": quotient_properties_bound,
        },
    })
}

fn series_json(series: ModuleDimensionSeries) -> Value {
    match series {
        ModuleDimensionSeries::ZeroDimensionVector(vector) => json!(vector),
        ModuleDimensionSeries::Layers(layers) => json!(layers),
    }
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn module_series_case<F: Field>(field: F) -> Value {
    let quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "a"),
            ArrowSpec::named(vertex(2), vertex(3), "b"),
        ],
    )
    .unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let projectives = indec_projective_modules(&algebra).unwrap();
    let path_simple = simple_modules(&algebra).unwrap()[0].clone();
    let zero = zero_module(&algebra).unwrap();
    let semisimple = right_module_over_path_algebra(algebra, vec![1, 1, 1], &[])
        .unwrap()
        .unwrap();

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let parent = PathAlgebra::new(field, loop_quiver).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let x_squared = x.multiply(&x).unwrap();
    let quotient = QuiverAlgebra::Quotient(
        quotient_by_relations(parent, &[x_squared.multiply(&x).unwrap()]).unwrap(),
    );
    let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
    let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
    let quotient_zero = zero_module(&quotient).unwrap();

    json!({
        "path": {
            "dimensions": projectives.iter().map(PathAlgebraModule::dimension)
                .chain([zero.dimension(), semisimple.dimension()]).collect::<Vec<_>>(),
            "partial_order": [
                dimension_vector_partial_order(&projectives[2], &projectives[0]),
                dimension_vector_partial_order(&projectives[0], &projectives[1]),
                dimension_vector_partial_order(&projectives[0], &projectives[0]),
                dimension_vector_partial_order(&projectives[1], &projectives[0]),
            ],
            "radical_series": projectives.iter().chain(std::iter::once(&semisimple))
                .map(|module| series_json(radical_series(module).unwrap()))
                .collect::<Vec<_>>(),
            "socle_series": projectives.iter().chain(std::iter::once(&semisimple))
                .map(|module| series_json(socle_series(module).unwrap()))
                .collect::<Vec<_>>(),
            "zero_radical_series": series_json(radical_series(&zero).unwrap()),
            "zero_socle_series": series_json(socle_series(&zero).unwrap()),
        },
        "quotient": {
            "dimensions": [quotient_projective.dimension(), quotient_simple.dimension(),
                quotient_zero.dimension()],
            "partial_order": [
                dimension_vector_partial_order(&quotient_simple, &quotient_projective),
                dimension_vector_partial_order(&quotient_projective, &quotient_simple),
            ],
            "cross_length_partial_order": [
                dimension_vector_partial_order(&quotient_simple, &path_simple),
                dimension_vector_partial_order(&path_simple, &quotient_simple),
                dimension_vector_partial_order(&quotient_simple, &projectives[0]),
                dimension_vector_partial_order(&projectives[0], &quotient_simple),
                dimension_vector_partial_order(&quotient_simple, &projectives[1]),
                dimension_vector_partial_order(&projectives[1], &quotient_simple),
            ],
            "radical_projective": series_json(radical_series(&quotient_projective).unwrap()),
            "socle_projective": series_json(socle_series(&quotient_projective).unwrap()),
            "radical_simple": series_json(radical_series(&quotient_simple).unwrap()),
            "socle_simple": series_json(socle_series(&quotient_simple).unwrap()),
            "zero_radical_series": series_json(radical_series(&quotient_zero).unwrap()),
            "zero_socle_series": series_json(socle_series(&quotient_zero).unwrap()),
        },
    })
}

fn annihilator_json<F: Field>(
    algebra: &QuiverAlgebra<F>,
    elements: &[QuiverAlgebraElement<F>],
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let canonical_basis = finite_canonical_basis(algebra).unwrap().unwrap();
    let default_basis = finite_qpa_default_basis(algebra, &canonical_basis).unwrap();
    let basis_change = DenseMatrix::from_rows(
        default_basis
            .iter()
            .map(|element| {
                canonical_basis_coefficients(algebra, &canonical_basis, element).unwrap()
            })
            .collect(),
        algebra.field().clone(),
    )
    .unwrap();
    json!(
        elements
            .iter()
            .map(|element| {
                let canonical_coordinates =
                    canonical_basis_coefficients(algebra, &canonical_basis, element).unwrap();
                basis_change
                    .solution_row(&canonical_coordinates)
                    .unwrap()
                    .unwrap()
                    .iter()
                    .map(coefficient)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    )
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn module_annihilator_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let projectives = indec_projective_modules(&algebra).unwrap();
    let faithful = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
    let simples = simple_modules(&algebra).unwrap();
    let zero = zero_module(&algebra).unwrap();
    let arrow_zero = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
        .unwrap()
        .unwrap();
    let zero_annihilator = annihilator_of_module(&zero).unwrap();
    let simple1_annihilator = annihilator_of_module(&simples[0]).unwrap();
    let simple2_annihilator = annihilator_of_module(&simples[1]).unwrap();
    let arrow_zero_annihilator = annihilator_of_module(&arrow_zero).unwrap();
    let faithful_annihilator = annihilator_of_module(&faithful).unwrap();

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let quotient =
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap());
    let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
    let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
    let quotient_zero = zero_module(&quotient).unwrap();
    let quotient_two_zero = right_module_over_path_algebra(quotient.clone(), vec![2], &[])
        .unwrap()
        .unwrap();
    let quotient_projective_annihilator = annihilator_of_module(&quotient_projective).unwrap();
    let quotient_simple_annihilator = annihilator_of_module(&quotient_simple).unwrap();
    let quotient_zero_annihilator = annihilator_of_module(&quotient_zero).unwrap();
    let quotient_two_zero_annihilator = annihilator_of_module(&quotient_two_zero).unwrap();

    let bound_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "c"),
            ArrowSpec::named(vertex(2), vertex(3), "d"),
        ],
    )
    .unwrap();
    let bound_parent = PathAlgebra::new(field, bound_quiver).unwrap();
    let bound_arrows = bound_parent.quiver().arrows().collect::<Vec<_>>();
    let c = bound_parent.element(bound_arrows[0].path()).unwrap();
    let d = bound_parent.element(bound_arrows[1].path()).unwrap();
    let bound_quotient = QuiverAlgebra::Quotient(
        quotient_by_relations(bound_parent, &[c.multiply(&d).unwrap()]).unwrap(),
    );
    let bound_zero = zero_module(&bound_quotient).unwrap();
    let bound_simple = simple_modules(&bound_quotient).unwrap()[2].clone();
    let bound_zero_annihilator = annihilator_of_module(&bound_zero).unwrap();
    let bound_simple_annihilator = annihilator_of_module(&bound_simple).unwrap();

    json!({
        "path": {
            "basis_dimension": finite_canonical_basis(&algebra).unwrap().unwrap().len(),
            "zero": annihilator_json(&algebra, &zero_annihilator, coefficient),
            "simple1": annihilator_json(&algebra, &simple1_annihilator, coefficient),
            "simple2": annihilator_json(&algebra, &simple2_annihilator, coefficient),
            "arrow_zero": annihilator_json(&algebra, &arrow_zero_annihilator, coefficient),
            "faithful": annihilator_json(&algebra, &faithful_annihilator, coefficient),
            "cache_bound": [
                zero.cached_annihilator_of_module().is_some(),
                simples[0].cached_annihilator_of_module().is_some(),
                simples[1].cached_annihilator_of_module().is_some(),
                arrow_zero.cached_annihilator_of_module().is_some(),
                faithful.cached_annihilator_of_module().is_some(),
            ],
        },
        "quotient": {
            "basis_dimension": finite_canonical_basis(&quotient).unwrap().unwrap().len(),
            "projective": annihilator_json(
                &quotient,
                &quotient_projective_annihilator,
                coefficient,
            ),
            "simple": annihilator_json(
                &quotient,
                &quotient_simple_annihilator,
                coefficient,
            ),
            "zero": annihilator_json(&quotient, &quotient_zero_annihilator, coefficient),
            "two_zero": annihilator_json(
                &quotient,
                &quotient_two_zero_annihilator,
                coefficient,
            ),
            "cache_bound": [
                quotient_projective.cached_annihilator_of_module().is_some(),
                quotient_simple.cached_annihilator_of_module().is_some(),
                quotient_zero.cached_annihilator_of_module().is_some(),
                quotient_two_zero.cached_annihilator_of_module().is_some(),
            ],
        },
        "bound_quotient": {
            "basis_dimension": finite_canonical_basis(&bound_quotient)
                .unwrap().unwrap().len(),
            "zero": annihilator_json(
                &bound_quotient,
                &bound_zero_annihilator,
                coefficient,
            ),
            "simple3": annihilator_json(
                &bound_quotient,
                &bound_simple_annihilator,
                coefficient,
            ),
            "cache_bound": [
                bound_zero.cached_annihilator_of_module().is_some(),
                bound_simple.cached_annihilator_of_module().is_some(),
            ],
        },
    })
}

fn ext_json<F: Field>(
    ext: &ExtOverAlgebraResult<F>,
    field: &F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let basis_coordinates = if ext.basis().is_empty() {
        Value::Null
    } else {
        json!(
            ext.basis()
                .iter()
                .map(|homomorphism| ext
                    .coefficients(homomorphism)
                    .unwrap()
                    .unwrap()
                    .iter()
                    .map(coefficient)
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>()
        )
    };
    let (entered_coefficients, combo_coordinates) = if ext.basis().is_empty() {
        (Value::Null, Value::Null)
    } else {
        let entered = (0..ext.dimension())
            .map(|index| field.nth(Integer::from(i64::try_from(index + 2).unwrap())))
            .collect::<Vec<_>>();
        let mut terms = ext
            .basis()
            .iter()
            .zip(&entered)
            .map(|(homomorphism, scalar)| homomorphism.scale(scalar).unwrap());
        let mut combo = terms.next().unwrap();
        for term in terms {
            combo = combo.add(&term).unwrap();
        }
        let coordinates = ext.coefficients(&combo).unwrap().unwrap();
        (
            json!(entered.iter().map(coefficient).collect::<Vec<_>>()),
            json!(coordinates.iter().map(coefficient).collect::<Vec<_>>()),
        )
    };
    json!({
        "syzygy_inclusion": module_map_json(ext.syzygy_inclusion(), coefficient),
        "dimension": ext.dimension(),
        "basis": ext.basis().iter()
            .map(|homomorphism| module_map_json(homomorphism, coefficient))
            .collect::<Vec<_>>(),
        "basis_coordinates": basis_coordinates,
        "entered_coefficients": entered_coefficients,
        "combo_coordinates": combo_coordinates,
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn functors_star_nakayama_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let simples = simple_modules(&algebra).unwrap();
    let projectives = indec_projective_modules(&algebra).unwrap();
    let zero = zero_module(&algebra).unwrap();
    let cover = projective_cover(&simples[0]).unwrap();
    let inclusion = kernel_inclusion(&cover).unwrap();

    let mixed = right_module_over_path_algebra(
        algebra.clone(),
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![field.one()], vec![field.zero()]], field.clone())
                .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let mixed_map = right_module_hom_over_algebra(
        mixed.clone(),
        mixed.clone(),
        vec![
            DenseMatrix::from_rows(
                vec![
                    vec![field.one(), field.zero()],
                    vec![field.zero(), field.zero()],
                ],
                field.clone(),
            )
            .unwrap(),
            DenseMatrix::from_rows(vec![vec![field.one()]], field.clone()).unwrap(),
        ],
    )
    .unwrap();

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let quotient =
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap());
    let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
    let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
    let quotient_x = include_in_path_algebra(
        &quotient.quiver().arrows().next().unwrap().path(),
        &quotient,
    )
    .unwrap();
    let quotient_generator =
        minimal_generating_set_of_module(quotient_projective.clone()).unwrap()[0].clone();
    let radical_image = quotient_generator.act(&quotient_x).unwrap();
    let quotient_radical_map = hom_from_projective(&radical_image, &quotient_projective)
        .unwrap()
        .unwrap();

    json!({
        "modules": {
            "a2_s2": module_functor_attribute_json(&simples[1], coefficient),
            "a2_p1": module_functor_attribute_json(&projectives[0], coefficient),
            "a2_zero": module_functor_attribute_json(&zero, coefficient),
            "a2_mixed": module_functor_attribute_json(&mixed, coefficient),
            "dual_simple": module_functor_attribute_json(&quotient_simple, coefficient),
            "dual_projective": module_functor_attribute_json(
                &quotient_projective,
                coefficient,
            ),
        },
        "maps": {
            "a2_inclusion": map_functor_attribute_json(&inclusion, coefficient),
            "a2_endomorphism": map_functor_attribute_json(&mixed_map, coefficient),
            "dual_radical": map_functor_attribute_json(&quotient_radical_map, coefficient),
        },
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn functors_transpose_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let simples = simple_modules(&algebra).unwrap();
    let mixed = right_module_over_path_algebra(
        algebra,
        vec![2, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![field.one()], vec![field.zero()]], field.clone())
                .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();

    let kronecker_quiver = Quiver::from_indexed(
        2,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "p"),
            ArrowSpec::named(vertex(1), vertex(2), "q"),
        ],
    )
    .unwrap();
    let kronecker = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
    let kronecker_simple = simple_modules(&kronecker).unwrap()[0].clone();

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let x = loop_parent
        .element(loop_parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let dual_numbers = QuiverAlgebra::Quotient(
        quotient_by_relations(loop_parent, &[x.multiply(&x).unwrap()]).unwrap(),
    );
    let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();

    let bound_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "a"),
            ArrowSpec::named(vertex(2), vertex(3), "b"),
        ],
    )
    .unwrap();
    let bound_parent = PathAlgebra::new(field, bound_quiver).unwrap();
    let bound_arrows = bound_parent.quiver().arrows().collect::<Vec<_>>();
    let a = bound_parent.element(bound_arrows[0].path()).unwrap();
    let b = bound_parent.element(bound_arrows[1].path()).unwrap();
    let bound = QuiverAlgebra::Quotient(
        quotient_by_relations(bound_parent, &[a.multiply(&b).unwrap()]).unwrap(),
    );
    let bound_simple = simple_modules(&bound).unwrap()[0].clone();

    json!({
        "a2_s1": transpose_functor_attribute_json(&simples[0], coefficient),
        "a2_s2_projective": transpose_functor_attribute_json(&simples[1], coefficient),
        "a2_mixed": transpose_functor_attribute_json(&mixed, coefficient),
        "kronecker_s1": transpose_functor_attribute_json(&kronecker_simple, coefficient),
        "dual_simple": transpose_functor_attribute_json(&dual_simple, coefficient),
        "bound_a3_s1": transpose_functor_attribute_json(&bound_simple, coefficient),
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn functors_ar_maps_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let a2_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());
    let a2_module = right_module_over_path_algebra(a2, vec![2, 0], &[])
        .unwrap()
        .unwrap();
    let a2_map = right_module_hom_over_algebra(
        a2_module.clone(),
        a2_module,
        vec![
            DenseMatrix::from_rows(
                vec![
                    vec![field.one(), field.one()],
                    vec![field.zero(), field.one()],
                ],
                field.clone(),
            )
            .unwrap(),
            DenseMatrix::zero(1, 1, field.clone()),
        ],
    )
    .unwrap();

    let kronecker_quiver = Quiver::from_indexed(
        2,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "p"),
            ArrowSpec::named(vertex(1), vertex(2), "q"),
        ],
    )
    .unwrap();
    let kronecker = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
    let kronecker_module = right_module_over_path_algebra(kronecker, vec![2, 0], &[])
        .unwrap()
        .unwrap();
    let kronecker_map = right_module_hom_over_algebra(
        kronecker_module.clone(),
        kronecker_module,
        vec![
            DenseMatrix::from_rows(
                vec![
                    vec![field.one(), field.one()],
                    vec![field.zero(), field.one()],
                ],
                field.clone(),
            )
            .unwrap(),
            DenseMatrix::zero(1, 1, field.clone()),
        ],
    )
    .unwrap();

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let x = loop_parent
        .element(loop_parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let dual_numbers = QuiverAlgebra::Quotient(
        quotient_by_relations(loop_parent, &[x.multiply(&x).unwrap()]).unwrap(),
    );
    let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();
    let dual_identity = identity_module_mapping(dual_simple).unwrap();

    let bound_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "a"),
            ArrowSpec::named(vertex(2), vertex(3), "b"),
        ],
    )
    .unwrap();
    let bound_parent = PathAlgebra::new(field, bound_quiver).unwrap();
    let bound_arrows = bound_parent.quiver().arrows().collect::<Vec<_>>();
    let a = bound_parent.element(bound_arrows[0].path()).unwrap();
    let b = bound_parent.element(bound_arrows[1].path()).unwrap();
    let bound = QuiverAlgebra::Quotient(
        quotient_by_relations(bound_parent, &[a.multiply(&b).unwrap()]).unwrap(),
    );
    let bound_simple = simple_modules(&bound).unwrap()[0].clone();
    let bound_identity = identity_module_mapping(bound_simple).unwrap();

    json!({
        "a2_shear": ar_map_functor_attribute_json(&a2_map, coefficient),
        "kronecker_shear": ar_map_functor_attribute_json(&kronecker_map, coefficient),
        "dual_simple_identity": ar_map_functor_attribute_json(&dual_identity, coefficient),
        "bound_a3_identity": ar_map_functor_attribute_json(&bound_identity, coefficient),
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn functors_ar_map_iterations_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let a2_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());
    let a2_module = right_module_over_path_algebra(a2, vec![2, 0], &[])
        .unwrap()
        .unwrap();
    let a2_map = right_module_hom_over_algebra(
        a2_module.clone(),
        a2_module,
        vec![
            DenseMatrix::from_rows(
                vec![
                    vec![field.one(), field.one()],
                    vec![field.zero(), field.one()],
                ],
                field.clone(),
            )
            .unwrap(),
            DenseMatrix::zero(1, 1, field.clone()),
        ],
    )
    .unwrap();

    let kronecker_quiver = Quiver::from_indexed(
        2,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "p"),
            ArrowSpec::named(vertex(1), vertex(2), "q"),
        ],
    )
    .unwrap();
    let kronecker = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
    let kronecker_module = right_module_over_path_algebra(kronecker, vec![2, 0], &[])
        .unwrap()
        .unwrap();
    let kronecker_map = right_module_hom_over_algebra(
        kronecker_module.clone(),
        kronecker_module,
        vec![
            DenseMatrix::from_rows(
                vec![
                    vec![field.one(), field.one()],
                    vec![field.zero(), field.one()],
                ],
                field.clone(),
            )
            .unwrap(),
            DenseMatrix::zero(1, 1, field.clone()),
        ],
    )
    .unwrap();

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let loop_parent = PathAlgebra::new(field, loop_quiver).unwrap();
    let x = loop_parent
        .element(loop_parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let dual_numbers = QuiverAlgebra::Quotient(
        quotient_by_relations(loop_parent, &[x.multiply(&x).unwrap()]).unwrap(),
    );
    let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();
    let dual_identity = identity_module_mapping(dual_simple).unwrap();

    json!({
        "a2_shear": ar_map_iterations_json(&a2_map, coefficient),
        "kronecker_shear": ar_map_iterations_json(&kronecker_map, coefficient),
        "dual_simple_identity": ar_map_iterations_json(&dual_identity, coefficient),
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn functors_tensor_modules_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let a2_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());
    let simples = simple_modules(&a2).unwrap();
    let projectives = indec_projective_modules(&a2).unwrap();

    let zero_right = dual_of_module(&simples[1]).unwrap();
    let zero_tensor = tensor_product_of_modules(&simples[0], &zero_right).unwrap();
    let zero_element = zero_tensor
        .elementary_tensor(
            &simples[0].basis_elements()[0],
            &zero_right.basis_elements()[0],
        )
        .unwrap();

    let projective_right = dual_of_module(&simples[0]).unwrap();
    let projective_tensor = tensor_product_of_modules(&projectives[0], &projective_right).unwrap();
    let projective_element = caught_tensor_element_json(
        &projective_tensor,
        &projectives[0].basis_elements()[0],
        &projective_right.basis_elements()[0],
        coefficient,
    );
    assert_eq!(
        projective_tensor
            .elementary_tensor(
                &projectives[0].basis_elements()[0],
                &projective_right.basis_elements()[0],
            )
            .unwrap_err(),
        super::EnvelopingModuleError::TensorProductUnboundSizeMatrix
    );

    let path_main =
        tensor_product_of_modules(&simples[0], &dual_of_module(&projectives[0]).unwrap());
    assert_eq!(
        path_main.as_ref().unwrap_err(),
        &super::EnvelopingModuleError::TensorProductPathAlgebraP1Extraction
    );
    let wrong_algebras = tensor_product_of_modules(&simples[0], &simples[0]);
    assert_eq!(
        wrong_algebras.as_ref().unwrap_err(),
        &super::EnvelopingModuleError::TensorProductModulesOverWrongAlgebras
    );

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let x = loop_parent
        .element(loop_parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let cubic =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(3).unwrap()]).unwrap());
    let cubic_module = right_module_over_path_algebra(
        cubic.clone(),
        vec![2],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(
                vec![
                    vec![field.zero(), field.one()],
                    vec![field.zero(), field.zero()],
                ],
                field.clone(),
            )
            .unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let cubic_projective = indec_projective_modules(&cubic).unwrap()[0].clone();
    let cubic_right = dual_of_module(&cubic_projective).unwrap();
    let cubic_tensor = tensor_product_of_modules(&cubic_module, &cubic_right).unwrap();
    let cubic_elements = cubic_module
        .basis_elements()
        .iter()
        .map(|left| {
            cubic_right
                .basis_elements()
                .iter()
                .map(|right| caught_tensor_element_json(&cubic_tensor, left, right, coefficient))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let bound_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "a"),
            ArrowSpec::named(vertex(2), vertex(3), "b"),
        ],
    )
    .unwrap();
    let bound_parent = PathAlgebra::new(field, bound_quiver).unwrap();
    let bound_arrows = bound_parent.quiver().arrows().collect::<Vec<_>>();
    let a = bound_parent.element(bound_arrows[0].path()).unwrap();
    let b = bound_parent.element(bound_arrows[1].path()).unwrap();
    let bound = QuiverAlgebra::Quotient(
        quotient_by_relations(bound_parent, &[a.multiply(&b).unwrap()]).unwrap(),
    );
    let bound_simple = simple_modules(&bound).unwrap()[0].clone();
    let bound_projective = indec_projective_modules(&bound).unwrap()[0].clone();
    let bound_right = dual_of_module(&bound_projective).unwrap();
    let bound_tensor = tensor_product_of_modules(&bound_simple, &bound_right).unwrap();
    let bound_element = caught_tensor_element_json(
        &bound_tensor,
        &bound_simple.basis_elements()[0],
        &bound_right.basis_elements()[0],
        coefficient,
    );
    assert!(matches!(
        bound_tensor.elementary_tensor(
            &bound_simple.basis_elements()[0],
            &bound_right.basis_elements()[0]
        ),
        Err(
            super::EnvelopingModuleError::TensorProductClosureProjectionOutOfRange {
                index: 2,
                count: 1
            }
        )
    ));

    json!({
        "easy_zero": {
            "space": tensor_space_json(&zero_tensor, coefficient),
            "element": zero_element.iter().map(coefficient).collect::<Vec<_>>(),
        },
        "projective_branch": {
            "space": tensor_space_json(&projective_tensor, coefficient),
            "element": projective_element,
        },
        "path_main_constructor_success": path_main.is_ok(),
        "wrong_algebra_constructor_success": wrong_algebras.is_ok(),
        "cubic_quotient": {
            "space": tensor_space_json(&cubic_tensor, coefficient),
            "elements": cubic_elements,
        },
        "bound_capture_bug": {
            "space": tensor_space_json(&bound_tensor, coefficient),
            "element": bound_element,
        },
    })
}

#[allow(clippy::needless_pass_by_value)]
fn functors_iterations_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let a2_simple = simple_modules(&algebra).unwrap()[0].clone();

    let kronecker_quiver = Quiver::from_indexed(
        2,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "p"),
            ArrowSpec::named(vertex(1), vertex(2), "q"),
        ],
    )
    .unwrap();
    let kronecker = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
    let kronecker_simple = simple_modules(&kronecker).unwrap()[0].clone();

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let parent = PathAlgebra::new(field, loop_quiver).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let dual_numbers =
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap());
    let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();

    json!({
        "a2_s1": functor_iterations_json(&a2_simple, coefficient),
        "kronecker_s1": functor_iterations_json(&kronecker_simple, coefficient),
        "dual_simple": functor_iterations_json(&dual_simple, coefficient),
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_ext_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let simples = simple_modules(&algebra).unwrap();
    let zero = zero_module(&algebra).unwrap();
    let a2_12 = ext_over_algebra(&simples[0], &simples[1]).unwrap();
    let a2_11 = ext_over_algebra(&simples[0], &simples[0]).unwrap();
    let a2_22 = ext_over_algebra(&simples[1], &simples[1]).unwrap();
    let a2_02 = ext_over_algebra(&zero, &simples[1]).unwrap();

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let quotient =
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap());
    let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
    let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
    let dual_simple_ext = ext_over_algebra(&quotient_simple, &quotient_simple).unwrap();
    let dual_projective_ext = ext_over_algebra(&quotient_simple, &quotient_projective).unwrap();
    let quotient_sum =
        direct_sum_of_qpa_modules(&[quotient_simple.clone(), quotient_projective.clone()])
            .unwrap()
            .unwrap();
    let dual_sum = ext_over_algebra(&quotient_simple, &quotient_sum).unwrap();

    let kronecker_quiver = Quiver::from_indexed(
        2,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "p"),
            ArrowSpec::named(vertex(1), vertex(2), "q"),
        ],
    )
    .unwrap();
    let kronecker = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
    let kronecker_simples = simple_modules(&kronecker).unwrap();
    let kronecker_ext = ext_over_algebra(&kronecker_simples[0], &kronecker_simples[1]).unwrap();

    let rigid_json = |module: &PathAlgebraModule<F>| {
        let has_before = module.cached_is_rigid_module().is_some();
        let value = is_rigid_module(module).unwrap();
        let has_after = module.cached_is_rigid_module().is_some();
        json!({
            "has_before": has_before,
            "value": value,
            "has_after": has_after,
        })
    };

    json!({
        "a2_s1_s2": ext_json(&a2_12, &field, coefficient),
        "a2_s1_s1": ext_json(&a2_11, &field, coefficient),
        "a2_s2_s2": ext_json(&a2_22, &field, coefficient),
        "a2_zero_s2": ext_json(&a2_02, &field, coefficient),
        "dual_simple_simple": ext_json(&dual_simple_ext, &field, coefficient),
        "dual_simple_projective": ext_json(&dual_projective_ext, &field, coefficient),
        "dual_simple_sum": ext_json(&dual_sum, &field, coefficient),
        "kronecker_s1_s2": ext_json(&kronecker_ext, &field, coefficient),
        "rigid": {
            "a2_s1": rigid_json(&simples[0]),
            "a2_s2": rigid_json(&simples[1]),
            "a2_zero": rigid_json(&zero),
            "dual_simple": rigid_json(&quotient_simple),
            "dual_projective": rigid_json(&quotient_projective),
            "dual_sum": rigid_json(&quotient_sum),
            "kronecker_s1": rigid_json(&kronecker_simples[0]),
            "kronecker_s2": rigid_json(&kronecker_simples[1]),
        },
    })
}

fn ext_algebra_map_json<F: Field>(
    map: &PathAlgebraModuleHomomorphism<F>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    json!({
        "source_dimension": map.source().dimension_vector(),
        "range_dimension": map.range().dimension_vector(),
        "matrices": map.vertex_matrices().iter().map(|matrix| json!({
            "shape": [matrix.nrows(), matrix.ncols()],
            "rows": (0..matrix.nrows()).map(|row| matrix.row(row).iter()
                .map(coefficient).collect::<Vec<_>>()).collect::<Vec<_>>(),
        })).collect::<Vec<_>>(),
    })
}

fn ext_algebra_json<F: Field>(
    result: &ExtAlgebraGeneratorsResult<F>,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    json!({
        "dimensions": result.dimensions(),
        "generator_counts": result.generator_counts(),
        "generators": result.generators().iter().map(|degree| degree.iter()
            .map(|map| ext_algebra_map_json(map, coefficient))
            .collect::<Vec<_>>()).collect::<Vec<_>>(),
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehomalg_ext_algebra_case<F: QpaCoefficientField>(
    field: F,
    coefficient: impl Fn(&F::Element) -> Value + Copy,
) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let a2_sum = direct_sum_of_qpa_modules(&simple_modules(&a2).unwrap())
        .unwrap()
        .unwrap();

    let vector_algebra = QuiverAlgebra::Path(
        PathAlgebra::new(field.clone(), Quiver::from_indexed(1, Vec::new()).unwrap()).unwrap(),
    );
    let vector2 = right_module_over_path_algebra(vector_algebra, vec![2], &[])
        .unwrap()
        .unwrap();

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let parent = PathAlgebra::new(field, loop_quiver).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let x2 = x.multiply(&x).unwrap();
    let dual = QuiverAlgebra::Quotient(
        quotient_by_relations(parent.clone(), std::slice::from_ref(&x2)).unwrap(),
    );
    let dual_simple = simple_modules(&dual).unwrap()[0].clone();
    let dual_projective = indec_projective_modules(&dual).unwrap()[0].clone();
    let dual_sum = direct_sum_of_qpa_modules(&[dual_simple.clone(), dual_projective])
        .unwrap()
        .unwrap();

    let x3 = x2.multiply(&x).unwrap();
    let truncated =
        QuiverAlgebra::Quotient(quotient_by_relations(parent, std::slice::from_ref(&x3)).unwrap());
    let truncated_projective = indec_projective_modules(&truncated).unwrap()[0].clone();

    json!({
        "a2_sum_n2": ext_algebra_json(
            &ext_algebra_generators(&a2_sum, 2).unwrap(), coefficient
        ),
        "vector2_n2": ext_algebra_json(
            &ext_algebra_generators(&vector2, 2).unwrap(), coefficient
        ),
        "dual_simple_n3": ext_algebra_json(
            &ext_algebra_generators(&dual_simple, 3).unwrap(), coefficient
        ),
        "dual_sum_n3": ext_algebra_json(
            &ext_algebra_generators(&dual_sum, 3).unwrap(), coefficient
        ),
        "truncated_projective_n2": ext_algebra_json(
            &ext_algebra_generators(&truncated_projective, 2).unwrap(), coefficient
        ),
    })
}

#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn modulehom_syzygy_injective_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let zero = zero_module(&algebra).unwrap();
    let zero_again = zero_module(&algebra).unwrap();
    let simple1 = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let simple2 = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let syzygy1 = first_syzygy(&simple1).unwrap();
    let syzygy2 = first_syzygy(&simple2).unwrap();
    let syzygy0 = first_syzygy(&zero).unwrap();
    let envelope1 = injective_envelope(&simple1).unwrap();
    let envelope2 = injective_envelope(&simple2).unwrap();
    let envelope0 = injective_envelope(&zero).unwrap();

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let parent = PathAlgebra::new(field, loop_quiver).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let quotient =
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap());
    let quotient_zero = zero_module(&quotient).unwrap();
    let quotient_simple = right_module_over_path_algebra(quotient.clone(), vec![1], &[])
        .unwrap()
        .unwrap();
    let quotient_syzygy = first_syzygy(&quotient_simple).unwrap();
    let quotient_envelope = injective_envelope(&quotient_simple).unwrap();

    json!({
        "path": {
            "zero": attribute_module_json(&zero, coefficient),
            "zero_cache_identity": zero.is_identical(&zero_again),
            "syzygy_simple1": attribute_module_json(&syzygy1, coefficient),
            "syzygy_simple2": attribute_module_json(&syzygy2, coefficient),
            "syzygy_zero": attribute_module_json(&syzygy0, coefficient),
            "syzygy_cache_identity": [
                syzygy1.is_identical(&first_syzygy(&simple1).unwrap()),
                syzygy2.is_identical(&first_syzygy(&simple2).unwrap()),
                syzygy0.is_identical(&zero),
            ],
            "injective_simple1": module_map_json(&envelope1, coefficient),
            "injective_simple2": module_map_json(&envelope2, coefficient),
            "injective_zero": module_map_json(&envelope0, coefficient),
            "injective_cache_identity": [
                envelope1.is_identical(&injective_envelope(&simple1).unwrap()),
                envelope2.is_identical(&injective_envelope(&simple2).unwrap()),
                envelope0.is_identical(&injective_envelope(&zero).unwrap()),
            ],
            "injective_source_identity": [
                envelope1.source().is_identical(&simple1),
                envelope2.source().is_identical(&simple2),
                envelope0.source().is_identical(&zero),
            ],
        },
        "quotient": {
            "zero": attribute_module_json(&quotient_zero, coefficient),
            "syzygy_simple": attribute_module_json(&quotient_syzygy, coefficient),
            "injective_simple": module_map_json(&quotient_envelope, coefficient),
            "cache_identity": [
                quotient_zero.is_identical(&zero_module(&quotient).unwrap()),
                quotient_syzygy.is_identical(&first_syzygy(&quotient_simple).unwrap()),
                quotient_envelope
                    .is_identical(&injective_envelope(&quotient_simple).unwrap()),
            ],
        },
    })
}

fn stable_module_json<F: Field>(
    module: &PathAlgebraModule<F>,
    print_zero_vertex_map: bool,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let natural = end_modulo_projectives(module).unwrap();
    let factoring = hom_factoring_through_projectives(module, module).unwrap();
    let factor_matrices = factoring
        .iter()
        .map(|map| from_hom_mm_to_end_m(map).unwrap())
        .collect::<Vec<_>>();
    let quotient_basis = natural.range().basis_elements();
    let images = natural
        .source()
        .basis_matrices()
        .iter()
        .map(|matrix| {
            natural
                .image(matrix)
                .unwrap()
                .coordinates()
                .iter()
                .map(coefficient)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let one = natural.range().one().unwrap().map(|element| {
        element
            .coordinates()
            .iter()
            .map(coefficient)
            .collect::<Vec<_>>()
    });
    let multiplication = quotient_basis
        .iter()
        .map(|left| {
            quotient_basis
                .iter()
                .map(|right| {
                    natural
                        .range()
                        .multiply(left, right)
                        .unwrap()
                        .coordinates()
                        .iter()
                        .map(coefficient)
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let roundtrip_matrices = natural
        .source()
        .basis_matrices()
        .iter()
        .map(|matrix| from_hom_mm_to_end_m(&from_end_m_to_hom_mm(module, matrix).unwrap()).unwrap())
        .collect::<Vec<_>>();
    let mut value = json!({
        "end_basis": natural.source().basis_matrices().iter()
            .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
        "factor_matrices": factor_matrices.iter()
            .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
        "ideal_basis": natural.range().ideal_basis_matrices().iter()
            .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
        "kernel_basis": natural.kernel_basis_matrices().iter()
            .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
        "source_cache_identity": natural.source().is_identical(&end_over_algebra(module).unwrap()),
        "dimension": natural.range().dimension(),
        "images": images,
        "one": one,
        "multiplication": multiplication,
        "roundtrip_matrices": roundtrip_matrices.iter()
            .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
    });
    if print_zero_vertex_map {
        value["zero_vertex_map"] = module_map_json(
            &from_end_m_to_hom_mm(module, &natural.source().basis_matrices()[0]).unwrap(),
            coefficient,
        );
    }
    value
}

#[allow(clippy::needless_pass_by_value)]
fn modulehom_stable_case<F: Field>(
    field: F,
    coefficient: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let zero = field.zero();
    let one = field.one();
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
    let module = right_module_over_path_algebra(
        algebra.clone(),
        vec![2, 2],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(
                vec![
                    vec![one.clone(), zero.clone()],
                    vec![zero.clone(), zero.clone()],
                ],
                field.clone(),
            )
            .unwrap(),
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

    let loop_quiver =
        Quiver::from_named(vec!["w".into()], vec![NamedArrowSpec::named("w", "w", "x")]).unwrap();
    let path_algebra = PathAlgebra::new(field, loop_quiver).unwrap();
    let loop_element = path_algebra
        .element(path_algebra.quiver().arrows().next().unwrap().path())
        .unwrap();
    let dual_numbers = QuiverAlgebra::Quotient(
        quotient_by_relations(
            path_algebra,
            &[loop_element.multiply(&loop_element).unwrap()],
        )
        .unwrap(),
    );
    let dual_simple = right_module_over_path_algebra(dual_numbers.clone(), vec![1], &[])
        .unwrap()
        .unwrap();
    let dual_projective = indec_projective_modules(&dual_numbers).unwrap()[0].clone();

    json!({
        "a2_module": stable_module_json(&module, false, coefficient),
        "a2_simple_one": stable_module_json(&simple_one, true, coefficient),
        "a2_simple_two": stable_module_json(&simple_two, true, coefficient),
        "dual_simple": stable_module_json(&dual_simple, false, coefficient),
        "dual_projective": stable_module_json(&dual_projective, false, coefficient),
    })
}

#[test]
fn pinned_qpa_modulehom_projective_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_projectives.schema-v1.json"
    ))
    .unwrap();
    let rational = modulehom_projective_case(Q, ToString::to_string);
    let binary = modulehom_projective_case(Z2, |entry| {
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
fn pinned_qpa_modulehom_projective_matrix_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_projective_matrices.schema-v1.json"
    ))
    .unwrap();
    let rational = modulehom_projective_matrix_case(Q, ToString::to_string);
    let binary = modulehom_projective_matrix_case(Z2, |entry| {
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
fn pinned_qpa_modulehom_syzygy_injective_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_syzygy_injective.schema-v1.json"
    ))
    .unwrap();
    let rational = modulehom_syzygy_injective_case(Q, ToString::to_string);
    let binary = modulehom_syzygy_injective_case(Z2, |entry| {
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
fn pinned_qpa_module_families_properties_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/module_families_properties.schema-v1.json"
    ))
    .unwrap();
    let rational = module_families_properties_case(Q, ToString::to_string);
    let binary = module_families_properties_case(Z2, |entry| {
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
fn pinned_qpa_module_series_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/module_series.schema-v1.json"
    ))
    .unwrap();
    assert_eq!(module_series_case(Q), expected["rational"]);
    assert_eq!(module_series_case(Z2), expected["binary"]);
}

#[test]
fn pinned_qpa_module_annihilator_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/module_annihilator.schema-v1.json"
    ))
    .unwrap();
    let rational = module_annihilator_case(Q, ToString::to_string);
    let binary = module_annihilator_case(Z2, |entry| {
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
fn pinned_qpa_modulehom_ext_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_ext.schema-v1.json"
    ))
    .unwrap();
    let rational = modulehom_ext_case(Q, ToString::to_string);
    let binary = modulehom_ext_case(Z2, |entry| {
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
fn pinned_qpa_ext_algebra_generators_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_ext_algebra_generators.schema-v1.json"
    ))
    .unwrap();
    let rational = modulehomalg_ext_algebra_case(Q, |entry: &Rational| {
        json!(entry.numerator_ref().to_i64().unwrap())
    });
    let binary =
        modulehomalg_ext_algebra_case(Z2, |entry| json!(Z2.to_integer(entry).to_i64().unwrap()));
    assert_eq!(rational, expected["rational"]);
    assert_eq!(binary, expected["binary"]);
}

#[test]
fn pinned_qpa_star_and_nakayama_functor_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/functors_star_nakayama.schema-v1.json"
    ))
    .unwrap();
    let rational = functors_star_nakayama_case(Q, ToString::to_string);
    let binary = functors_star_nakayama_case(Z2, |entry| {
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
fn pinned_qpa_transpose_and_ar_functor_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/functors_transpose.schema-v1.json"
    ))
    .unwrap();
    let rational = functors_transpose_case(Q, ToString::to_string);
    let binary = functors_transpose_case(Z2, |entry| {
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
fn pinned_qpa_ar_functors_on_homomorphisms_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/functors_ar_maps.schema-v1.json"
    ))
    .unwrap();
    let rational = functors_ar_maps_case(Q, ToString::to_string);
    let binary = functors_ar_maps_case(Z2, |entry| {
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
fn pinned_qpa_ar_functor_homomorphism_iterations_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/functors_ar_map_iterations.schema-v1.json"
    ))
    .unwrap();
    let rational = functors_ar_map_iterations_case(Q, ToString::to_string);
    let binary = functors_ar_map_iterations_case(Z2, |entry| {
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
fn pinned_qpa_tensor_product_of_modules_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/functors_tensor_modules.schema-v1.json"
    ))
    .unwrap();
    let rational = functors_tensor_modules_case(Q, ToString::to_string);
    let binary = functors_tensor_modules_case(Z2, |entry| {
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
fn pinned_qpa_ar_functor_iterations_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/functors_iterations.schema-v1.json"
    ))
    .unwrap();
    let rational = functors_iterations_case(Q, ToString::to_string);
    let binary = functors_iterations_case(Z2, |entry| {
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
fn annihilator_of_module_reports_its_finite_basis_domain() {
    let quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let module = right_module_over_path_algebra(algebra, vec![1], &[])
        .unwrap()
        .unwrap();
    assert_eq!(
        annihilator_of_module(&module).unwrap_err(),
        super::EnvelopingModuleError::InfiniteAnnihilatorAlgebra
    );
}

#[test]
fn ext_over_algebra_preserves_domain_and_coefficient_boundaries() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let simples = simple_modules(&algebra).unwrap();
    let ext = ext_over_algebra(&simples[0], &simples[1]).unwrap();
    let outside = identity_module_mapping(simples[0].clone()).unwrap();
    assert_eq!(
        ext.coefficients(&outside).unwrap_err(),
        super::EnvelopingModuleError::ExtMapOutsideHomSpace
    );
    let zero_ext = ext_over_algebra(&simples[0], &simples[0]).unwrap();
    assert_eq!(zero_ext.coefficients(&outside).unwrap(), None);

    let other_quiver = Quiver::from_indexed(1, Vec::new()).unwrap();
    let other_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, other_quiver).unwrap());
    let other_simple = simple_modules(&other_algebra).unwrap()[0].clone();
    assert_eq!(
        ext_over_algebra(&simples[0], &other_simple).unwrap_err(),
        super::EnvelopingModuleError::ExtModulesOverDifferentAlgebras
    );
}

#[test]
fn simple_modules_preserves_qpa_source_domain_errors() {
    let quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
    let parent = PathAlgebra::new(Q, quiver).unwrap();
    assert_eq!(
        simple_modules(&QuiverAlgebra::Path(parent.clone())).unwrap_err(),
        super::EnvelopingModuleError::InfiniteSimpleModuleAlgebra
    );
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let quotient = QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x]).unwrap());
    assert_eq!(
        simple_modules(&quotient).unwrap_err(),
        super::EnvelopingModuleError::NonadmissibleSimpleModuleAlgebra
    );
}

#[test]
fn projective_matrix_conversions_report_qpa_domain_errors() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    assert_eq!(
        element_in_v_a_as_element_in_indec_proj(&algebra, &algebra.zero()),
        Err(super::EnvelopingModuleError::ZeroProjectiveAlgebraElement)
    );
    let v1 = include_in_path_algebra(
        &algebra.quiver().vertices().next().unwrap().path(),
        &algebra,
    )
    .unwrap();
    assert_eq!(
        from_matrix_to_homomorphism_of_projectives(&algebra, &[vec![v1]], &[2], &[1],),
        Err(super::EnvelopingModuleError::ProjectiveMatrixSupportDiscrepancy)
    );
    let scaled_arrow =
        include_in_path_algebra(&algebra.quiver().arrows().next().unwrap().path(), &algebra)
            .unwrap()
            .scale_left(&q(2))
            .unwrap();
    assert_eq!(
        element_in_v_a_as_element_in_indec_proj(&algebra, &scaled_arrow),
        Err(super::EnvelopingModuleError::ProjectiveLeadingTermMissing {
            left_vertex: 1,
            right_vertex: 2,
        })
    );
    assert_eq!(
        from_matrix_to_homomorphism_of_projectives(&algebra, &[], &[], &[]),
        Err(super::EnvelopingModuleError::EmptyProjectiveElementMatrix)
    );
}

#[test]
fn pinned_qpa_modulehom_stable_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_stable.schema-v1.json"
    ))
    .unwrap();
    let rational = modulehom_stable_case(Q, ToString::to_string);
    let binary = modulehom_stable_case(Z2, |entry| {
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
fn stable_endomorphism_quotient_rejects_the_zero_module_like_qpa() {
    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
    let zero = right_module_over_path_algebra(algebra, vec![0, 0], &[])
        .unwrap()
        .unwrap();

    assert_eq!(
        end_modulo_projectives(&zero).unwrap_err(),
        super::EnvelopingModuleError::ZeroStableEndomorphismModule
    );
}

#[test]
fn pinned_qpa_modulehom_functor_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_functors.schema-v1.json"
    ))
    .unwrap();
    let rational = modulehom_functor_case(Q, ToString::to_string);
    let binary = modulehom_functor_case(Z2, |entry| {
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
fn pinned_qpa_regular_enveloping_module_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/pathalgtensor_regular_module.schema-v1.json"
    ))
    .unwrap();

    let loop_quiver =
        Quiver::from_named(vec!["v".into()], vec![NamedArrowSpec::named("v", "v", "x")]).unwrap();
    let parent = PathAlgebra::new(Q, loop_quiver.clone()).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let rational_algebra =
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap());
    let rational_cached_before = rational_algebra.has_regular_enveloping_module_attribute();
    let rational_module = algebra_as_module_over_enveloping_algebra(&rational_algebra)
        .unwrap()
        .unwrap();
    let rational_again = algebra_as_module_over_enveloping_algebra(&rational_algebra)
        .unwrap()
        .unwrap();
    let rational_synonym = algebra_as_module_of_enveloping_algebra(&rational_algebra)
        .unwrap()
        .unwrap();
    let rational_enveloping = enveloping_algebra(&rational_algebra).unwrap();
    let rational_json = json!({
        "attribute_identity": [
            rational_cached_before,
            rational_algebra.has_regular_enveloping_module_attribute(),
            rational_module.is_identical(&rational_again),
            rational_module.is_identical(&rational_synonym),
            rational_module.right_acting_algebra() == rational_enveloping.algebra(),
        ],
        "module": module_json(&rational_module, ToString::to_string),
    });

    let path_quiver = Quiver::from_named(
        vec!["u1".into(), "u2".into()],
        vec![NamedArrowSpec::named("u1", "u2", "a")],
    )
    .unwrap();
    let binary_algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, path_quiver).unwrap());
    let binary_module = algebra_as_module_over_enveloping_algebra(&binary_algebra)
        .unwrap()
        .unwrap();
    let binary_again = algebra_as_module_of_enveloping_algebra(&binary_algebra)
        .unwrap()
        .unwrap();
    let binary_enveloping = enveloping_algebra(&binary_algebra).unwrap();
    let binary_json = json!({
        "attribute_identity": [
            binary_module.is_identical(&binary_again),
            binary_module.right_acting_algebra() == binary_enveloping.algebra(),
        ],
        "module": module_json(&binary_module, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        }),
    });

    let infinite_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, loop_quiver).unwrap());
    let infinite_result = algebra_as_module_over_enveloping_algebra(&infinite_algebra).unwrap();
    let infinite_json = json!({
        "is_fail": infinite_result.is_none(),
        "attribute_bound": infinite_algebra.has_regular_enveloping_module_attribute(),
        "enveloping_bound": infinite_algebra.cached_enveloping_algebra().is_some(),
    });

    assert_eq!(rational_json, expected["rational_quotient"]);
    assert_eq!(binary_json, expected["binary_path"]);
    assert_eq!(infinite_json, expected["infinite_path"]);
}

#[test]
fn pinned_qpa_dual_enveloping_module_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/pathalgtensor_dual_module.schema-v1.json"
    ))
    .unwrap();

    let loop_quiver =
        Quiver::from_named(vec!["v".into()], vec![NamedArrowSpec::named("v", "v", "x")]).unwrap();
    let parent = PathAlgebra::new(Q, loop_quiver.clone()).unwrap();
    let x = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let rational_algebra =
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap());
    let cached_before = rational_algebra.has_dual_enveloping_module_attribute();
    let bimodule_dual = dual_of_algebra_as_module_over_enveloping_algebra(&rational_algebra)
        .unwrap()
        .unwrap();
    let bimodule_again = dual_of_algebra_as_module_over_enveloping_algebra(&rational_algebra)
        .unwrap()
        .unwrap();
    let regular = algebra_as_module_over_enveloping_algebra(&rational_algebra)
        .unwrap()
        .unwrap();
    let generic_dual = dual_of_module(&regular).unwrap();
    let generic_again = dual_of_module(&regular).unwrap();
    let generic_back = dual_of_module(&generic_dual).unwrap();
    let rational_enveloping = enveloping_algebra(&rational_algebra).unwrap();
    let rational_json = json!({
        "attribute_identity": [
            cached_before,
            rational_algebra.has_dual_enveloping_module_attribute(),
            bimodule_dual.is_identical(&bimodule_again),
            bimodule_dual.right_acting_algebra() == rational_enveloping.algebra(),
        ],
        "generic_identity": [
            generic_dual.is_identical(&generic_again),
            generic_back.is_identical(&regular),
        ],
        "generic_dual": module_json(&generic_dual, ToString::to_string),
        "bimodule_dual": module_json(&bimodule_dual, ToString::to_string),
    });

    let path_quiver = Quiver::from_named(
        vec!["u1".into(), "u2".into()],
        vec![NamedArrowSpec::named("u1", "u2", "a")],
    )
    .unwrap();
    let binary_algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, path_quiver).unwrap());
    let binary_dual = dual_of_algebra_as_module_over_enveloping_algebra(&binary_algebra)
        .unwrap()
        .unwrap();
    let binary_again = dual_of_algebra_as_module_over_enveloping_algebra(&binary_algebra)
        .unwrap()
        .unwrap();
    let binary_enveloping = enveloping_algebra(&binary_algebra).unwrap();
    let binary_json = json!({
        "attribute_identity": [
            binary_dual.is_identical(&binary_again),
            binary_dual.right_acting_algebra() == binary_enveloping.algebra(),
        ],
        "bimodule_dual": module_json(&binary_dual, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        }),
    });

    let infinite_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, loop_quiver).unwrap());
    let infinite_result =
        dual_of_algebra_as_module_over_enveloping_algebra(&infinite_algebra).unwrap();
    let infinite_json = json!({
        "is_fail": infinite_result.is_none(),
        "attribute_bound": infinite_algebra.has_dual_enveloping_module_attribute(),
        "regular_bound": infinite_algebra.has_regular_enveloping_module_attribute(),
        "enveloping_bound": infinite_algebra.cached_enveloping_algebra().is_some(),
    });

    assert_eq!(rational_json, expected["rational_quotient"]);
    assert_eq!(binary_json, expected["binary_path"]);
    assert_eq!(infinite_json, expected["infinite_path"]);
}

#[test]
fn pinned_qpa_random_module_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/module_random.schema-v1.json"
    ))
    .unwrap();
    let one = std::num::NonZeroUsize::new(1).unwrap();

    let mut binary_rng = ZeroRng;
    let binary = random_module(&isolated_path_algebra(Z2, 1), one, one, &mut binary_rng).unwrap();
    let binary_json = json!({
        "success": true,
        "is_fail": binary.is_none(),
        "dimension_vector": binary.as_ref().map(PathAlgebraModule::dimension_vector),
    });

    let mut ternary_rng = StdRng::seed_from_u64(1);
    let ternary = random_module(
        &isolated_path_algebra(Zp::new(3), 1),
        one,
        one,
        &mut ternary_rng,
    )
    .unwrap();
    let ternary_json = json!({
        "success": true,
        "is_fail": ternary.is_none(),
        "dimension_vector": ternary.as_ref().map(PathAlgebraModule::dimension_vector),
    });

    let mut rational_rng = StdRng::seed_from_u64(1);
    let rational =
        random_module(&isolated_path_algebra(Q, 1), one, one, &mut rational_rng).unwrap();
    let rational_json = json!({
        "success": true,
        "is_fail": rational.is_none(),
        "dimension_vector": rational.as_ref().map(PathAlgebraModule::dimension_vector),
    });

    let mut isolated_rng = StdRng::seed_from_u64(1);
    let isolated =
        random_module(&isolated_path_algebra(Z2, 100), one, one, &mut isolated_rng).unwrap();
    let isolated_json = json!({
        "success": true,
        "is_fail": isolated.is_none(),
        "dimension_vector": isolated.as_ref().map(PathAlgebraModule::dimension_vector),
    });

    assert_eq!(binary_json, expected["binary"]);
    assert_eq!(ternary_json, expected["ternary"]);
    assert_eq!(rational_json, expected["rational"]);
    assert_eq!(isolated_json, expected["isolated"]);
}

fn modulehomalg_gorenstein_dimension_case<F: QpaCoefficientField>(field: &F) -> Value {
    let discrete = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
        )
        .unwrap(),
    );
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let radical_square = fresh_radical_square_zero_a3(field);
    let known_finite = fresh_radical_square_zero_a3(field);
    let dual = fresh_dual_numbers(field);
    let known_infinite = fresh_dual_numbers(field);
    let negative_one = Integer::from(-1);
    let zero_bound = Integer::from(0);
    let two_bound = Integer::from(2);

    let discrete_bound_0 = gorenstein_dimension_of_algebra(&discrete, &zero_bound).unwrap();
    let a2_bound_0 = gorenstein_dimension_of_algebra(&a2, &zero_bound).unwrap();
    let a2_attribute = gorenstein_dimension(&a2).unwrap();
    let a2_cached_negative = gorenstein_dimension_of_algebra(&a2, &negative_one).unwrap();

    let radical_initially_bound = gorenstein_dimension(&radical_square).is_some();
    let radical_bound_0 = gorenstein_dimension_of_algebra(&radical_square, &zero_bound).unwrap();
    let radical_bound_after_0 = gorenstein_dimension(&radical_square).is_some();
    let radical_bound_2 = gorenstein_dimension_of_algebra(&radical_square, &two_bound).unwrap();
    let radical_attribute = gorenstein_dimension(&radical_square).unwrap();
    let radical_cached_bound_0 =
        gorenstein_dimension_of_algebra(&radical_square, &zero_bound).unwrap();

    let known_finite_setup = global_dimension_of_algebra(&known_finite, 2).unwrap();
    let known_finite_negative =
        gorenstein_dimension_of_algebra(&known_finite, &negative_one).unwrap();

    let dual_bound_0 = gorenstein_dimension_of_algebra(&dual, &zero_bound).unwrap();
    let dual_attribute = gorenstein_dimension(&dual).unwrap();
    let known_infinite_setup = global_dimension_of_algebra(&known_infinite, 0).unwrap();
    let known_infinite_bound_0 =
        gorenstein_dimension_of_algebra(&known_infinite, &zero_bound).unwrap();

    json!({
        "discrete_bound_0": bounded_dimension_json(discrete_bound_0),
        "a2_bound_0": bounded_dimension_json(a2_bound_0),
        "a2_attribute": a2_attribute,
        "a2_cached_negative": bounded_dimension_json(a2_cached_negative),
        "radical_initially_bound": radical_initially_bound,
        "radical_bound_0": bounded_dimension_json(radical_bound_0),
        "radical_bound_after_0": radical_bound_after_0,
        "radical_bound_2": bounded_dimension_json(radical_bound_2),
        "radical_attribute": radical_attribute,
        "radical_cached_bound_0": bounded_dimension_json(radical_cached_bound_0),
        "known_finite_setup": global_dimension_json(known_finite_setup),
        "known_finite_negative": bounded_dimension_json(known_finite_negative),
        "dual_bound_0": bounded_dimension_json(dual_bound_0),
        "dual_attribute": dual_attribute,
        "known_infinite_setup": global_dimension_json(known_infinite_setup),
        "known_infinite_bound_0": bounded_dimension_json(known_infinite_bound_0),
    })
}

#[test]
fn pinned_qpa_modulehomalg_gorenstein_dimension_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_gorenstein_dimension.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_gorenstein_dimension",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_gorenstein_dimension_case(&Z2),
        "quinary": modulehomalg_gorenstein_dimension_case(&Zp::new(5)),
        "rational": modulehomalg_gorenstein_dimension_case(&Q),
    });
    assert_eq!(actual, expected);
}

fn modulehomalg_n_rigid_case<F: QpaCoefficientField>(field: &F) -> Value {
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let a2_simple = simple_modules(&a2).unwrap()[0].clone();
    let zero = zero_module(&a2).unwrap();
    let dual = fresh_dual_numbers(field);
    let dual_simple = simple_modules(&dual).unwrap()[0].clone();
    let dual_projective = indec_projective_modules(&dual).unwrap()[0].clone();
    let cycle = fresh_radical_square_zero_cycle(field);
    let cycle_simple = simple_modules(&cycle).unwrap()[0].clone();

    json!({
        "a2_simple_negative": n_rigid_module(&a2_simple, &Integer::from(-1)).unwrap(),
        "a2_simple_0": n_rigid_module(&a2_simple, &Integer::from(0)).unwrap(),
        "a2_simple_3": n_rigid_module(&a2_simple, &Integer::from(3)).unwrap(),
        "zero_3": n_rigid_module(&zero, &Integer::from(3)).unwrap(),
        "dual_simple_0": n_rigid_module(&dual_simple, &Integer::from(0)).unwrap(),
        "dual_simple_1": n_rigid_module(&dual_simple, &Integer::from(1)).unwrap(),
        "dual_projective_3": n_rigid_module(&dual_projective, &Integer::from(3)).unwrap(),
        "cycle_simple_1": n_rigid_module(&cycle_simple, &Integer::from(1)).unwrap(),
        "cycle_simple_2": n_rigid_module(&cycle_simple, &Integer::from(2)).unwrap(),
        "cycle_simple_3": n_rigid_module(&cycle_simple, &Integer::from(3)).unwrap(),
    })
}

#[test]
fn pinned_qpa_modulehomalg_n_rigid_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_n_rigid.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_n_rigid",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_n_rigid_case(&Z2),
        "quinary": modulehomalg_n_rigid_case(&Zp::new(5)),
        "rational": modulehomalg_n_rigid_case(&Q),
    });
    assert_eq!(actual, expected);
}

fn modulehomalg_dominant_dimension_case<F: QpaCoefficientField>(field: &F) -> Value {
    let zero_scalar = field.zero();
    let one_scalar = field.one();
    let a2 = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap(),
    );
    let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let simple_two = right_module_over_path_algebra(a2.clone(), vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let projective = right_module_over_path_algebra(
        a2.clone(),
        vec![1, 1],
        &[NamedArrowMatrix::new(
            "a",
            DenseMatrix::from_rows(vec![vec![one_scalar]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let zero = zero_module(&a2).unwrap();

    let loop_parent = PathAlgebra::new(
        field.clone(),
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
    )
    .unwrap();
    let x = loop_parent.generators()[1].clone();
    let dual =
        QuiverAlgebra::Quotient(quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap());
    let dual_simple = right_module_over_path_algebra(
        dual.clone(),
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![zero_scalar]], field.clone()).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    let discrete = QuiverAlgebra::Path(
        PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
        )
        .unwrap(),
    );
    let negative_one = Integer::from(-1);
    let zero_bound = Integer::from(0);
    let one_bound = Integer::from(1);
    let two_bound = Integer::from(2);

    json!({
        "simple_one_bound_0": dominant_dimension_json(
            dominant_dimension_of_module(&simple_one, &zero_bound).unwrap()),
        "simple_one_bound_2": dominant_dimension_json(
            dominant_dimension_of_module(&simple_one, &two_bound).unwrap()),
        "simple_one_bound_negative": dominant_dimension_json(
            dominant_dimension_of_module(&simple_one, &negative_one).unwrap()),
        "simple_two_bound_0": dominant_dimension_json(
            dominant_dimension_of_module(&simple_two, &zero_bound).unwrap()),
        "simple_two_bound_1": dominant_dimension_json(
            dominant_dimension_of_module(&simple_two, &one_bound).unwrap()),
        "simple_two_bound_2": dominant_dimension_json(
            dominant_dimension_of_module(&simple_two, &two_bound).unwrap()),
        "simple_two_bound_negative": dominant_dimension_json(
            dominant_dimension_of_module(&simple_two, &negative_one).unwrap()),
        "projective_injective": dominant_dimension_json(
            dominant_dimension_of_module(&projective, &zero_bound).unwrap()),
        "zero": dominant_dimension_json(
            dominant_dimension_of_module(&zero, &zero_bound).unwrap()),
        "dual_simple_bound_2": dominant_dimension_json(
            dominant_dimension_of_module(&dual_simple, &two_bound).unwrap()),
        "a2_bound_0": dominant_dimension_json(
            dominant_dimension_of_algebra(&a2, &zero_bound).unwrap()),
        "a2_bound_1": dominant_dimension_json(
            dominant_dimension_of_algebra(&a2, &one_bound).unwrap()),
        "a2_bound_2": dominant_dimension_json(
            dominant_dimension_of_algebra(&a2, &two_bound).unwrap()),
        "dual_algebra": dominant_dimension_json(
            dominant_dimension_of_algebra(&dual, &two_bound).unwrap()),
        "discrete_algebra": dominant_dimension_json(
            dominant_dimension_of_algebra(&discrete, &zero_bound).unwrap()),
    })
}

#[test]
fn pinned_qpa_modulehomalg_dominant_dimension_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_dominant_dimension.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_dominant_dimension",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_dominant_dimension_case(&Z2),
        "quinary": modulehomalg_dominant_dimension_case(&Zp::new(5)),
        "rational": modulehomalg_dominant_dimension_case(&Q),
    });
    assert_eq!(actual, expected);
}

#[test]
fn dominant_dimension_preserves_infinite_dimensional_try_next_method() {
    let algebra = QuiverAlgebra::Path(
        PathAlgebra::new(
            Z2,
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap(),
    );
    let module = right_module_over_path_algebra(
        algebra.clone(),
        vec![1],
        &[NamedArrowMatrix::new(
            "x",
            DenseMatrix::from_rows(vec![vec![Z2.zero()]], Z2).unwrap(),
        )],
    )
    .unwrap()
    .unwrap();
    assert_eq!(
        dominant_dimension_of_module(&module, &Integer::from(0)).unwrap(),
        None
    );
    assert_eq!(
        dominant_dimension_of_algebra(&algebra, &Integer::from(0)).unwrap(),
        None
    );
}

fn modulehomalg_dimensions_case<F: QpaCoefficientField>(field: &F) -> Value {
    let negative_one = Integer::from(-1);
    let zero_bound = Integer::from(0);
    let one_bound = Integer::from(1);
    let two_bound = Integer::from(2);

    let pd_one = fresh_a2_module(field, vec![1, 0]);
    let pd_one_initially_bound = proj_dimension(&pd_one).is_some();
    let pd_one_bound_0 = proj_dimension_of_module(&pd_one, &zero_bound).unwrap();
    let pd_one_bound_after_0 = proj_dimension(&pd_one).is_some();
    let pd_one_bound_1 = proj_dimension_of_module(&pd_one, &one_bound).unwrap();
    let pd_one_attribute = proj_dimension(&pd_one).unwrap();
    let pd_one_cached_bound_0 = proj_dimension_of_module(&pd_one, &zero_bound).unwrap();
    let pd_one_cached_bound_2 = proj_dimension_of_module(&pd_one, &two_bound).unwrap();

    let pd_two = fresh_a2_module(field, vec![0, 1]);
    let pd_two_bound_0 = proj_dimension_of_module(&pd_two, &zero_bound).unwrap();
    let pd_dual = fresh_dual_simple(field);
    let pd_dual_bound_2 = proj_dimension_of_module(&pd_dual, &two_bound).unwrap();
    let pd_negative = fresh_a2_module(field, vec![0, 1]);
    let pd_negative_bound = proj_dimension_of_module(&pd_negative, &negative_one).unwrap();
    let pd_negative_bound_after = proj_dimension(&pd_negative).is_some();
    let pd_zero = fresh_a2_zero(field);
    let pd_zero_bound_0 = proj_dimension_of_module(&pd_zero, &zero_bound).unwrap();

    let id_two = fresh_a2_module(field, vec![0, 1]);
    let id_two_initially_bound = inj_dimension(&id_two).is_some();
    let id_two_bound_0 = inj_dimension_of_module(&id_two, &zero_bound).unwrap();
    let id_two_bound_after_0 = inj_dimension(&id_two).is_some();
    let id_two_bound_1 = inj_dimension_of_module(&id_two, &one_bound).unwrap();
    let id_two_attribute = inj_dimension(&id_two).unwrap();
    let id_two_cached_bound_0 = inj_dimension_of_module(&id_two, &zero_bound).unwrap();
    let id_two_cached_bound_2 = inj_dimension_of_module(&id_two, &two_bound).unwrap();

    let id_one = fresh_a2_module(field, vec![1, 0]);
    let id_one_bound_0 = inj_dimension_of_module(&id_one, &zero_bound).unwrap();
    let id_dual = fresh_dual_simple(field);
    let id_dual_bound_2 = inj_dimension_of_module(&id_dual, &two_bound).unwrap();
    let id_negative = fresh_a2_module(field, vec![1, 0]);
    let id_negative_bound = inj_dimension_of_module(&id_negative, &negative_one).unwrap();
    let id_negative_bound_after = inj_dimension(&id_negative).is_some();
    let id_zero = fresh_a2_zero(field);
    let id_zero_bound_0 = inj_dimension_of_module(&id_zero, &zero_bound).unwrap();

    let shortcut = fresh_a2_module(field, vec![0, 1]);
    let shortcut_dual = dual_of_module(&shortcut).unwrap();
    let id_dual_cache_setup = proj_dimension_of_module(&shortcut_dual, &one_bound).unwrap();
    let id_dual_cache_shortcut_bound_0 = inj_dimension_of_module(&shortcut, &zero_bound).unwrap();
    let id_dual_cache_shortcut_bound_after = inj_dimension(&shortcut).is_some();

    json!({
        "pd_one_initially_bound": pd_one_initially_bound,
        "pd_one_bound_0": bounded_dimension_json(pd_one_bound_0),
        "pd_one_bound_after_0": pd_one_bound_after_0,
        "pd_one_bound_1": bounded_dimension_json(pd_one_bound_1),
        "pd_one_attribute": pd_one_attribute,
        "pd_one_cached_bound_0": bounded_dimension_json(pd_one_cached_bound_0),
        "pd_one_cached_bound_2": bounded_dimension_json(pd_one_cached_bound_2),
        "pd_two_bound_0": bounded_dimension_json(pd_two_bound_0),
        "pd_dual_bound_2": bounded_dimension_json(pd_dual_bound_2),
        "pd_negative_bound": bounded_dimension_json(pd_negative_bound),
        "pd_negative_bound_after": pd_negative_bound_after,
        "pd_zero_bound_0": bounded_dimension_json(pd_zero_bound_0),
        "id_two_initially_bound": id_two_initially_bound,
        "id_two_bound_0": bounded_dimension_json(id_two_bound_0),
        "id_two_bound_after_0": id_two_bound_after_0,
        "id_two_bound_1": bounded_dimension_json(id_two_bound_1),
        "id_two_attribute": id_two_attribute,
        "id_two_cached_bound_0": bounded_dimension_json(id_two_cached_bound_0),
        "id_two_cached_bound_2": bounded_dimension_json(id_two_cached_bound_2),
        "id_one_bound_0": bounded_dimension_json(id_one_bound_0),
        "id_dual_bound_2": bounded_dimension_json(id_dual_bound_2),
        "id_negative_bound": bounded_dimension_json(id_negative_bound),
        "id_negative_bound_after": id_negative_bound_after,
        "id_zero_bound_0": bounded_dimension_json(id_zero_bound_0),
        "id_dual_cache_setup": bounded_dimension_json(id_dual_cache_setup),
        "id_dual_cache_shortcut_bound_0": bounded_dimension_json(
            id_dual_cache_shortcut_bound_0),
        "id_dual_cache_shortcut_bound_after": id_dual_cache_shortcut_bound_after,
    })
}

#[test]
fn pinned_qpa_modulehomalg_dimensions_fixture_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehomalg_dimensions.schema-v1.json"
    ))
    .unwrap();
    let actual = json!({
        "schema_version": 1,
        "kind": "modulehomalg_dimensions",
        "oracle": {
            "gap_version": "4.17dev",
            "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa_version": "1.37",
            "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
            "gbnp_version": "1.1.0",
            "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
        },
        "binary": modulehomalg_dimensions_case(&Z2),
        "quinary": modulehomalg_dimensions_case(&Zp::new(5)),
        "rational": modulehomalg_dimensions_case(&Q),
    });
    assert_eq!(actual, expected);
}

#[test]
fn pinned_qpa_global_dimension_of_algebra_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/modulehom_global_dimension.schema-v1.json"
    ))
    .unwrap();

    let discrete = isolated_path_algebra(Z2, 1);
    let a2_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let a2 = QuiverAlgebra::Path(PathAlgebra::new(Z2, a2_quiver).unwrap());

    let a3_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "x"),
            ArrowSpec::named(vertex(2), vertex(3), "y"),
        ],
    )
    .unwrap();
    let a3_parent = PathAlgebra::new(Z2, a3_quiver).unwrap();
    let arrows = a3_parent.quiver().arrows().collect::<Vec<_>>();
    let x = a3_parent.element(arrows[0].path()).unwrap();
    let y = a3_parent.element(arrows[1].path()).unwrap();
    let radical_square_zero = QuiverAlgebra::Quotient(
        quotient_by_relations(a3_parent, &[x.multiply(&y).unwrap()]).unwrap(),
    );

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
    let loop_parent = PathAlgebra::new(Z2, loop_quiver).unwrap();
    let loop_arrow = loop_parent
        .element(loop_parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let dual_numbers = QuiverAlgebra::Quotient(
        quotient_by_relations(loop_parent, &[loop_arrow.multiply(&loop_arrow).unwrap()]).unwrap(),
    );

    let result = json!({
        "discrete": [global_dimension_json(global_dimension_of_algebra(&discrete, 0).unwrap())],
        "path_a2": [
            global_dimension_json(global_dimension_of_algebra(&a2, 0).unwrap()),
            global_dimension_json(global_dimension_of_algebra(&a2, 3).unwrap())
        ],
        "radical_square_zero_a3": [
            global_dimension_json(global_dimension_of_algebra(&radical_square_zero, 0).unwrap()),
            global_dimension_json(global_dimension_of_algebra(&radical_square_zero, 1).unwrap()),
            global_dimension_json(global_dimension_of_algebra(&radical_square_zero, 2).unwrap()),
            global_dimension_json(global_dimension_of_algebra(&radical_square_zero, 0).unwrap())
        ],
        "dual_numbers": [global_dimension_json(global_dimension_of_algebra(&dual_numbers, 0).unwrap())]
    });
    assert_eq!(result, expected["f2"]);
}

#[allow(clippy::too_many_lines)]
fn algebra_homological_property_cases<F: QpaCoefficientField>(field: F) -> Value {
    let zero = field.zero();
    let one = field.one();
    let semisimple_structure = StructureConstantsAlgebra::new(
        field.clone(),
        2,
        vec![
            one.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
        ],
    )
    .unwrap();
    let dual_structure = StructureConstantsAlgebra::new(
        field.clone(),
        2,
        vec![
            one.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
            zero.clone(),
            one.clone(),
            zero.clone(),
            zero.clone(),
        ],
    )
    .unwrap();

    let structure_semisimple_known_before = semisimple_structure
        .cached_is_semisimple_algebra()
        .is_some();
    let structure_semisimple_value = semisimple_structure.is_semisimple_algebra().unwrap();
    let structure_semisimple_known_after = semisimple_structure
        .cached_is_semisimple_algebra()
        .is_some();
    let structure_global_known_before = semisimple_structure.cached_global_dimension().is_some();
    let structure_global_dimension =
        cached_global_dimension_json(semisimple_structure.global_dimension_from_known_properties());
    let structure_global_known_after = semisimple_structure.cached_global_dimension().is_some();

    let structure_dual_known_before = dual_structure.cached_is_semisimple_algebra().is_some();
    let structure_dual_value = dual_structure.is_semisimple_algebra().unwrap();
    let structure_dual_known_after = dual_structure.cached_is_semisimple_algebra().is_some();

    let discrete = isolated_path_algebra(field.clone(), 1);
    let a2_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let a2_parent = PathAlgebra::new(field.clone(), a2_quiver).unwrap();
    let a2 = QuiverAlgebra::Path(a2_parent.clone());
    let a2_quotient =
        QuiverAlgebra::Quotient(quotient_by_relations(a2_parent.clone(), &[]).unwrap());

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
    let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let loop_path = QuiverAlgebra::Path(loop_parent.clone());
    let t = loop_parent
        .element(loop_parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let t_squared = t.multiply(&t).unwrap();
    let dual_numbers = QuiverAlgebra::Quotient(
        quotient_by_relations(loop_parent.clone(), std::slice::from_ref(&t_squared)).unwrap(),
    );
    let dual_numbers_reverse = QuiverAlgebra::Quotient(
        quotient_by_relations(loop_parent, std::slice::from_ref(&t_squared)).unwrap(),
    );

    let a3_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "x"),
            ArrowSpec::named(vertex(2), vertex(3), "y"),
        ],
    )
    .unwrap();
    let a3_parent = PathAlgebra::new(field, a3_quiver).unwrap();
    let arrows = a3_parent.quiver().arrows().collect::<Vec<_>>();
    let x = a3_parent.element(arrows[0].path()).unwrap();
    let y = a3_parent.element(arrows[1].path()).unwrap();
    let radical_square_zero_a3 = QuiverAlgebra::Quotient(
        quotient_by_relations(a3_parent, &[x.multiply(&y).unwrap()]).unwrap(),
    );

    let QuiverAlgebra::Quotient(a2_quotient_object) = &a2_quotient else {
        unreachable!()
    };
    let a2_quotient_semisimple = admissible_quotient_is_semisimple(a2_quotient_object).unwrap();
    let a2_quotient_hereditary = is_hereditary_algebra(&a2_quotient).unwrap().unwrap();

    let QuiverAlgebra::Quotient(dual_numbers_object) = &dual_numbers else {
        unreachable!()
    };
    let dual_global_initially = dual_numbers_object.cached_global_dimension().is_some();
    let dual_selfinjective = is_selfinjective_algebra(&dual_numbers).unwrap().unwrap();
    let dual_global_after_selfinjective = dual_numbers_object.cached_global_dimension().is_some();
    let dual_semisimple = admissible_quotient_is_semisimple(dual_numbers_object).unwrap();
    let dual_global_after_semisimple = dual_numbers_object.cached_global_dimension().is_some();
    let dual_global_dimension =
        cached_global_dimension_json(global_dimension_from_known_properties(&dual_numbers));
    let dual_global_after_global = dual_numbers_object.cached_global_dimension().is_some();
    let dual_hereditary = is_hereditary_algebra(&dual_numbers).unwrap().unwrap();

    let QuiverAlgebra::Quotient(dual_reverse_object) = &dual_numbers_reverse else {
        unreachable!()
    };
    let dual_reverse_semisimple = admissible_quotient_is_semisimple(dual_reverse_object).unwrap();
    let dual_reverse_global_after_semisimple =
        dual_reverse_object.cached_global_dimension().is_some();
    let dual_reverse_selfinjective = is_selfinjective_algebra(&dual_numbers_reverse)
        .unwrap()
        .unwrap();
    let dual_reverse_global_after_selfinjective =
        dual_reverse_object.cached_global_dimension().is_some();

    let QuiverAlgebra::Quotient(a3_object) = &radical_square_zero_a3 else {
        unreachable!()
    };
    let a3_selfinjective = is_selfinjective_algebra(&radical_square_zero_a3)
        .unwrap()
        .unwrap();
    let a3_semisimple = admissible_quotient_is_semisimple(a3_object).unwrap();
    let a3_hereditary = is_hereditary_algebra(&radical_square_zero_a3)
        .unwrap()
        .unwrap();
    let a3_global_after_hereditary = a3_object.cached_global_dimension().is_some();
    let a3_bounded =
        global_dimension_json(global_dimension_of_algebra(&radical_square_zero_a3, 2).unwrap());
    let a3_global_after_bounded = a3_object.cached_global_dimension().is_some();

    json!({
        "structure_semisimple": {
            "known_before": structure_semisimple_known_before,
            "value": structure_semisimple_value,
            "known_after": structure_semisimple_known_after,
            "global_known_before": structure_global_known_before,
            "global_dimension": structure_global_dimension,
            "global_known_after": structure_global_known_after,
            "finite_global_dimension": semisimple_structure
                .is_finite_global_dimension_from_known_properties().unwrap()
        },
        "structure_dual_numbers": {
            "known_before": structure_dual_known_before,
            "value": structure_dual_value,
            "known_after": structure_dual_known_after,
            "global_known": dual_structure.cached_global_dimension().is_some()
        },
        "path_discrete": {
            "selfinjective": is_selfinjective_algebra(&discrete).unwrap().unwrap(),
            "semisimple": discrete.quiver().number_of_arrows() == 0,
            "hereditary": is_hereditary_algebra(&discrete).unwrap().unwrap(),
            "finite_global_dimension":
                is_finite_global_dimension_from_known_properties(&discrete).unwrap(),
            "global_dimension": cached_global_dimension_json(
                global_dimension_from_known_properties(&discrete))
        },
        "path_a2": {
            "selfinjective": is_selfinjective_algebra(&a2).unwrap().unwrap(),
            "semisimple": a2.quiver().number_of_arrows() == 0,
            "hereditary": is_hereditary_algebra(&a2).unwrap().unwrap(),
            "finite_global_dimension":
                is_finite_global_dimension_from_known_properties(&a2).unwrap(),
            "global_dimension": cached_global_dimension_json(
                global_dimension_from_known_properties(&a2))
        },
        "path_loop": {
            "selfinjective": is_selfinjective_algebra(&loop_path).unwrap().unwrap(),
            "hereditary": is_hereditary_algebra(&loop_path).unwrap().unwrap(),
            "finite_global_dimension":
                is_finite_global_dimension_from_known_properties(&loop_path).unwrap(),
            "global_dimension": cached_global_dimension_json(
                global_dimension_from_known_properties(&loop_path))
        },
        "quotient_a2_zero_ideal": {
            "semisimple": a2_quotient_semisimple,
            "hereditary": a2_quotient_hereditary,
            "finite_global_dimension":
                is_finite_global_dimension_from_known_properties(&a2_quotient).unwrap(),
            "global_dimension": cached_global_dimension_json(
                global_dimension_from_known_properties(&a2_quotient))
        },
        "quotient_dual_numbers": {
            "global_known_initially": dual_global_initially,
            "selfinjective": dual_selfinjective,
            "global_known_after_selfinjective": dual_global_after_selfinjective,
            "semisimple": dual_semisimple,
            "global_known_after_semisimple": dual_global_after_semisimple,
            "global_dimension": dual_global_dimension,
            "global_known_after_global": dual_global_after_global,
            "hereditary": dual_hereditary
        },
        "quotient_dual_numbers_semisimple_first": {
            "semisimple": dual_reverse_semisimple,
            "global_known_after_semisimple": dual_reverse_global_after_semisimple,
            "selfinjective": dual_reverse_selfinjective,
            "global_known_after_selfinjective": dual_reverse_global_after_selfinjective,
            "global_dimension": cached_global_dimension_json(
                global_dimension_from_known_properties(&dual_numbers_reverse))
        },
        "quotient_radical_square_zero_a3": {
            "selfinjective": a3_selfinjective,
            "semisimple": a3_semisimple,
            "hereditary": a3_hereditary,
            "global_known_after_hereditary": a3_global_after_hereditary,
            "bounded_global_dimension": a3_bounded,
            "global_known_after_bounded": a3_global_after_bounded
        }
    })
}

#[test]
fn pinned_qpa_algebra_homological_properties_match() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/algebra_homological_properties.schema-v1.json"
    ))
    .unwrap();
    assert_eq!(algebra_homological_property_cases(Z2), expected["f2"]);
    assert_eq!(algebra_homological_property_cases(Q), expected["rational"]);

    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let parent = PathAlgebra::new(Z2, quiver).unwrap();
    let arrow = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let nonadmissible = QuiverAlgebra::Quotient(
        quotient_by_relations(parent, std::slice::from_ref(&arrow)).unwrap(),
    );
    assert_eq!(is_selfinjective_algebra(&nonadmissible).unwrap(), None);
    assert_eq!(is_hereditary_algebra(&nonadmissible).unwrap(), None);
}

fn symmetry_property_json<F: Field>(algebra: &QuiverAlgebra<F>) -> Value {
    json!({
        "symmetric": is_symmetric_algebra(algebra).unwrap().unwrap(),
        "weakly_symmetric": is_weakly_symmetric_algebra(algebra).unwrap().unwrap(),
        "schurian": is_schurian_algebra(algebra).unwrap().unwrap(),
        "semicommutative": is_semicommutative_algebra(algebra).unwrap().unwrap()
    })
}

#[allow(clippy::too_many_lines)]
fn pathalg_symmetry_property_cases<F: Field>(field: F) -> Value {
    let discrete = isolated_path_algebra(field.clone(), 1);
    let a2_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
    let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let loop_path = QuiverAlgebra::Path(loop_parent.clone());
    let loop_arrow = loop_parent
        .element(loop_parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let dual_numbers = QuiverAlgebra::Quotient(
        quotient_by_relations(loop_parent, &[loop_arrow.multiply(&loop_arrow).unwrap()]).unwrap(),
    );

    let diamond_quiver = Quiver::from_indexed(
        4,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "a"),
            ArrowSpec::named(vertex(1), vertex(3), "b"),
            ArrowSpec::named(vertex(2), vertex(4), "c"),
            ArrowSpec::named(vertex(3), vertex(4), "d"),
        ],
    )
    .unwrap();
    let diamond_parent = PathAlgebra::new(field.clone(), diamond_quiver).unwrap();
    let diamond_arrows = diamond_parent.quiver().arrows().collect::<Vec<_>>();
    let diamond_a = diamond_parent.element(diamond_arrows[0].path()).unwrap();
    let diamond_b = diamond_parent.element(diamond_arrows[1].path()).unwrap();
    let diamond_c = diamond_parent.element(diamond_arrows[2].path()).unwrap();
    let diamond_d = diamond_parent.element(diamond_arrows[3].path()).unwrap();
    let ac = diamond_a.multiply(&diamond_c).unwrap();
    let bd = diamond_b.multiply(&diamond_d).unwrap();
    let diamond_commutative = QuiverAlgebra::Quotient(
        quotient_by_relations(diamond_parent.clone(), &[ac.sub(&bd).unwrap()]).unwrap(),
    );
    let diamond_one_path_zero = QuiverAlgebra::Quotient(
        quotient_by_relations(diamond_parent.clone(), std::slice::from_ref(&ac)).unwrap(),
    );
    let diamond_path = QuiverAlgebra::Path(diamond_parent);

    let cycle_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "x"),
            ArrowSpec::named(vertex(2), vertex(3), "y"),
            ArrowSpec::named(vertex(3), vertex(1), "z"),
        ],
    )
    .unwrap();
    let cycle_parent = PathAlgebra::new(field, cycle_quiver).unwrap();
    let cycle_arrows = cycle_parent.quiver().arrows().collect::<Vec<_>>();
    let cycle_x = cycle_parent.element(cycle_arrows[0].path()).unwrap();
    let cycle_y = cycle_parent.element(cycle_arrows[1].path()).unwrap();
    let cycle_z = cycle_parent.element(cycle_arrows[2].path()).unwrap();
    let radical_square_zero_cycle = QuiverAlgebra::Quotient(
        quotient_by_relations(
            cycle_parent,
            &[
                cycle_x.multiply(&cycle_y).unwrap(),
                cycle_y.multiply(&cycle_z).unwrap(),
                cycle_z.multiply(&cycle_x).unwrap(),
            ],
        )
        .unwrap(),
    );

    json!({
        "path_discrete": symmetry_property_json(&discrete),
        "path_a2": symmetry_property_json(&a2),
        "path_loop": {
            "symmetric": is_symmetric_algebra(&loop_path).unwrap().unwrap(),
            "weakly_symmetric": is_weakly_symmetric_algebra(&loop_path).unwrap().unwrap()
        },
        "path_diamond": symmetry_property_json(&diamond_path),
        "quotient_dual_numbers": symmetry_property_json(&dual_numbers),
        "quotient_radical_square_zero_cycle":
            symmetry_property_json(&radical_square_zero_cycle),
        "quotient_diamond_commutative": symmetry_property_json(&diamond_commutative),
        "quotient_diamond_one_path_zero": symmetry_property_json(&diamond_one_path_zero)
    })
}

#[test]
fn pinned_qpa_pathalg_symmetry_properties_match() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/pathalg_symmetry_properties.schema-v1.json"
    ))
    .unwrap();
    assert_eq!(pathalg_symmetry_property_cases(Z2), expected["f2"]);
    assert_eq!(pathalg_symmetry_property_cases(Q), expected["rational"]);

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
    let loop_path = QuiverAlgebra::Path(PathAlgebra::new(Z2, loop_quiver).unwrap());
    assert_eq!(
        is_schurian_algebra(&loop_path),
        Err(EnvelopingModuleError::InfiniteSchurianAlgebra)
    );
    assert_eq!(
        is_semicommutative_algebra(&loop_path),
        Err(EnvelopingModuleError::InfiniteSchurianAlgebra)
    );

    let quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let parent = PathAlgebra::new(Z2, quiver).unwrap();
    let arrow = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let nonadmissible = QuiverAlgebra::Quotient(
        quotient_by_relations(parent, std::slice::from_ref(&arrow)).unwrap(),
    );
    assert_eq!(is_symmetric_algebra(&nonadmissible).unwrap(), None);
    assert_eq!(is_weakly_symmetric_algebra(&nonadmissible).unwrap(), None);
    assert_eq!(is_schurian_algebra(&nonadmissible).unwrap(), None);
    assert_eq!(is_semicommutative_algebra(&nonadmissible).unwrap(), None);
}

fn rational_matrix_rows_json(matrix: &DenseMatrix<RationalField>) -> Value {
    Value::Array(
        (0..matrix.nrows())
            .map(|row| {
                Value::Array(
                    matrix
                        .row(row)
                        .iter()
                        .map(|entry| Value::String(entry.to_string()))
                        .collect(),
                )
            })
            .collect(),
    )
}

fn cartan_coxeter_case_json<F: Field>(algebra: &QuiverAlgebra<F>) -> Value {
    let cartan_initially_known = algebra.cached_cartan_matrix_attribute().is_bound();
    let matrix_initially_known = algebra.cached_coxeter_matrix_attribute().is_bound();
    let polynomial_initially_known = algebra.cached_coxeter_polynomial_attribute().is_bound();
    let cartan = cartan_matrix(algebra).unwrap().unwrap();
    let matrix = coxeter_matrix(algebra).unwrap().unwrap();
    let polynomial = coxeter_polynomial(algebra).unwrap().unwrap();
    json!({
        "cartan_initially_known": cartan_initially_known,
        "matrix_initially_known": matrix_initially_known,
        "polynomial_initially_known": polynomial_initially_known,
        "cartan": cartan,
        "coxeter": rational_matrix_rows_json(&matrix),
        "polynomial_coefficients": polynomial
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>(),
        "cartan_known_after": algebra.cached_cartan_matrix_attribute().is_bound(),
        "matrix_known_after": algebra.cached_coxeter_matrix_attribute().is_bound(),
        "polynomial_known_after": algebra.cached_coxeter_polynomial_attribute().is_bound()
    })
}

#[allow(clippy::too_many_lines)]
fn pathalg_cartan_coxeter_cases<F: Field>(field: F) -> Value {
    let discrete = isolated_path_algebra(field.clone(), 2);
    let a3_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "x"),
            ArrowSpec::named(vertex(2), vertex(3), "y"),
        ],
    )
    .unwrap();
    let a3 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a3_quiver).unwrap());

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
    let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let t = loop_parent
        .element(loop_parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let dual_numbers = QuiverAlgebra::Quotient(
        quotient_by_relations(loop_parent, &[t.multiply(&t).unwrap()]).unwrap(),
    );

    let diamond_quiver = Quiver::from_indexed(
        4,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "a"),
            ArrowSpec::named(vertex(1), vertex(3), "b"),
            ArrowSpec::named(vertex(2), vertex(4), "c"),
            ArrowSpec::named(vertex(3), vertex(4), "d"),
        ],
    )
    .unwrap();
    let diamond_parent = PathAlgebra::new(field, diamond_quiver).unwrap();
    let arrows = diamond_parent.quiver().arrows().collect::<Vec<_>>();
    let diamond_a = diamond_parent.element(arrows[0].path()).unwrap();
    let diamond_b = diamond_parent.element(arrows[1].path()).unwrap();
    let diamond_c = diamond_parent.element(arrows[2].path()).unwrap();
    let diamond_d = diamond_parent.element(arrows[3].path()).unwrap();
    let relation = diamond_a
        .multiply(&diamond_c)
        .unwrap()
        .sub(&diamond_b.multiply(&diamond_d).unwrap())
        .unwrap();
    let commutative_diamond =
        QuiverAlgebra::Quotient(quotient_by_relations(diamond_parent, &[relation]).unwrap());

    json!({
        "path_discrete": cartan_coxeter_case_json(&discrete),
        "path_a3": cartan_coxeter_case_json(&a3),
        "quotient_dual_numbers": cartan_coxeter_case_json(&dual_numbers),
        "quotient_commutative_diamond":
            cartan_coxeter_case_json(&commutative_diamond)
    })
}

#[test]
fn pinned_qpa_pathalg_cartan_coxeter_match() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/pathalg_cartan_coxeter.schema-v1.json"
    ))
    .unwrap();
    assert_eq!(pathalg_cartan_coxeter_cases(Z2), expected["f2"]);
    assert_eq!(pathalg_cartan_coxeter_cases(Zp::new(5)), expected["f5"]);
    assert_eq!(pathalg_cartan_coxeter_cases(Q), expected["rational"]);

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
    let loop_path = QuiverAlgebra::Path(PathAlgebra::new(Z2, loop_quiver).unwrap());
    assert_eq!(
        loop_path.cached_cartan_matrix_attribute(),
        CachedAttribute::Unbound
    );
    assert_eq!(cartan_matrix(&loop_path).unwrap(), None);
    assert_eq!(
        loop_path.cached_cartan_matrix_attribute(),
        CachedAttribute::Fail
    );
    assert_eq!(coxeter_matrix(&loop_path).unwrap(), None);
    assert_eq!(
        loop_path.cached_coxeter_matrix_attribute(),
        CachedAttribute::Fail
    );
    assert_eq!(coxeter_polynomial(&loop_path).unwrap(), None);
    assert_eq!(
        loop_path.cached_coxeter_polynomial_attribute(),
        CachedAttribute::Fail
    );

    let arrow_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let nonadmissible_parent = PathAlgebra::new(Z2, arrow_quiver).unwrap();
    let arrow = nonadmissible_parent
        .element(
            nonadmissible_parent
                .quiver()
                .arrows()
                .next()
                .unwrap()
                .path(),
        )
        .unwrap();
    let nonadmissible =
        QuiverAlgebra::Quotient(quotient_by_relations(nonadmissible_parent, &[arrow]).unwrap());
    assert_eq!(cartan_matrix(&nonadmissible).unwrap(), None);
    assert_eq!(coxeter_matrix(&nonadmissible).unwrap(), None);
    assert_eq!(coxeter_polynomial(&nonadmissible).unwrap(), None);

    let cycle_quiver = Quiver::from_indexed(
        2,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "x"),
            ArrowSpec::named(vertex(2), vertex(1), "y"),
        ],
    )
    .unwrap();
    let cycle_parent = PathAlgebra::new(Z2, cycle_quiver).unwrap();
    let cycle_arrows = cycle_parent.quiver().arrows().collect::<Vec<_>>();
    let x = cycle_parent.element(cycle_arrows[0].path()).unwrap();
    let y = cycle_parent.element(cycle_arrows[1].path()).unwrap();
    let singular = QuiverAlgebra::Quotient(
        quotient_by_relations(
            cycle_parent,
            &[x.multiply(&y).unwrap(), y.multiply(&x).unwrap()],
        )
        .unwrap(),
    );
    assert_eq!(
        cartan_matrix(&singular).unwrap(),
        Some(vec![vec![1, 1], vec![1, 1]])
    );
    assert_eq!(coxeter_matrix(&singular).unwrap(), None);
    assert_eq!(
        singular.cached_coxeter_matrix_attribute(),
        CachedAttribute::Fail
    );
    assert_eq!(coxeter_polynomial(&singular).unwrap(), None);
    assert_eq!(
        singular.cached_coxeter_polynomial_attribute(),
        CachedAttribute::Fail
    );
}

fn nakayama_permutation_case_json<F: Field>(algebra: &QuiverAlgebra<F>) -> Value {
    let initially_known = algebra.cached_nakayama_permutation_attribute().is_bound();
    let Some(permutation) = nakayama_permutation(algebra).unwrap() else {
        return json!({
            "initially_known": initially_known,
            "value": false,
            "known_after": algebra
                .cached_nakayama_permutation_attribute()
                .is_bound()
        });
    };
    let simples = simple_modules(algebra).unwrap();
    let module_images = simples
        .iter()
        .map(|simple| {
            permutation
                .image_of_module(simple)
                .unwrap()
                .dimension_vector()
                .to_vec()
        })
        .collect::<Vec<_>>();
    let module_identical = simples
        .iter()
        .map(|simple| {
            permutation
                .image_of_module(simple)
                .unwrap()
                .is_identical(simple)
        })
        .collect::<Vec<_>>();
    let index_images = (1..=simples.len())
        .map(|index| permutation.image_of_index(&Integer::from(index)).unwrap())
        .collect::<Vec<_>>();
    json!({
        "initially_known": initially_known,
        "index_images": index_images,
        "module_images": module_images,
        "module_identical": module_identical,
        "known_after": algebra
            .cached_nakayama_permutation_attribute()
            .is_bound()
    })
}

#[allow(clippy::too_many_lines)]
fn pathalg_nakayama_permutation_cases<F: Field>(field: F) -> Value {
    let discrete = isolated_path_algebra(field.clone(), 3);
    let a2_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
    let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let loop_arrow = loop_parent
        .element(loop_parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let dual_numbers = QuiverAlgebra::Quotient(
        quotient_by_relations(loop_parent, &[loop_arrow.multiply(&loop_arrow).unwrap()]).unwrap(),
    );

    let a3_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "u"),
            ArrowSpec::named(vertex(2), vertex(3), "v"),
        ],
    )
    .unwrap();
    let a3_parent = PathAlgebra::new(field.clone(), a3_quiver).unwrap();
    let a3_arrows = a3_parent.quiver().arrows().collect::<Vec<_>>();
    let first_arrow = a3_parent.element(a3_arrows[0].path()).unwrap();
    let second_arrow = a3_parent.element(a3_arrows[1].path()).unwrap();
    let nonselfinjective = QuiverAlgebra::Quotient(
        quotient_by_relations(a3_parent, &[first_arrow.multiply(&second_arrow).unwrap()]).unwrap(),
    );

    let cycle_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "x"),
            ArrowSpec::named(vertex(2), vertex(3), "y"),
            ArrowSpec::named(vertex(3), vertex(1), "z"),
        ],
    )
    .unwrap();
    let cycle_parent = PathAlgebra::new(field, cycle_quiver).unwrap();
    let arrows = cycle_parent.quiver().arrows().collect::<Vec<_>>();
    let cycle_x = cycle_parent.element(arrows[0].path()).unwrap();
    let cycle_y = cycle_parent.element(arrows[1].path()).unwrap();
    let cycle_z = cycle_parent.element(arrows[2].path()).unwrap();
    let cycle = QuiverAlgebra::Quotient(
        quotient_by_relations(
            cycle_parent,
            &[
                cycle_x.multiply(&cycle_y).unwrap(),
                cycle_y.multiply(&cycle_z).unwrap(),
                cycle_z.multiply(&cycle_x).unwrap(),
            ],
        )
        .unwrap(),
    );

    json!({
        "path_discrete": nakayama_permutation_case_json(&discrete),
        "path_a2": nakayama_permutation_case_json(&a2),
        "quotient_dual_numbers": nakayama_permutation_case_json(&dual_numbers),
        "quotient_nonselfinjective": nakayama_permutation_case_json(&nonselfinjective),
        "quotient_radical_square_zero_cycle": nakayama_permutation_case_json(&cycle)
    })
}

#[test]
fn pinned_qpa_pathalg_nakayama_permutation_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/pathalg_nakayama_permutation.schema-v1.json"
    ))
    .unwrap();
    assert_eq!(pathalg_nakayama_permutation_cases(Z2), expected["f2"]);
    assert_eq!(
        pathalg_nakayama_permutation_cases(Zp::new(5)),
        expected["f5"]
    );
    assert_eq!(pathalg_nakayama_permutation_cases(Q), expected["rational"]);

    let cycle_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "x"),
            ArrowSpec::named(vertex(2), vertex(3), "y"),
            ArrowSpec::named(vertex(3), vertex(1), "z"),
        ],
    )
    .unwrap();
    let cycle_parent = PathAlgebra::new(Z2, cycle_quiver).unwrap();
    let arrows = cycle_parent.quiver().arrows().collect::<Vec<_>>();
    let x = cycle_parent.element(arrows[0].path()).unwrap();
    let y = cycle_parent.element(arrows[1].path()).unwrap();
    let z = cycle_parent.element(arrows[2].path()).unwrap();
    let cycle = QuiverAlgebra::Quotient(
        quotient_by_relations(
            cycle_parent,
            &[
                x.multiply(&y).unwrap(),
                y.multiply(&z).unwrap(),
                z.multiply(&x).unwrap(),
            ],
        )
        .unwrap(),
    );
    let permutation = nakayama_permutation(&cycle).unwrap().unwrap();
    let cycle_simples = simple_modules(&cycle).unwrap();
    let nonsimple = direct_sum_of_qpa_modules(&cycle_simples[..2])
        .unwrap()
        .unwrap();
    assert_eq!(
        permutation.image_of_module(&nonsimple),
        Err(EnvelopingModuleError::NakayamaPermutationModuleNotSimple)
    );
    for index in [0, -1, 4] {
        assert_eq!(
            permutation.image_of_index(&Integer::from(index)),
            Err(EnvelopingModuleError::IncorrectNakayamaPermutationIndex {
                index: Integer::from(index),
                vertex_count: 3
            })
        );
    }

    let foreign = isolated_path_algebra(Z2, 4);
    let foreign_simple = simple_modules(&foreign).unwrap()[3].clone();
    assert_eq!(
        permutation.image_of_module(&foreign_simple),
        Err(
            EnvelopingModuleError::NakayamaPermutationModulePositionOutOfRange {
                position: 4,
                vertex_count: 3
            }
        )
    );

    let arrow_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let parent = PathAlgebra::new(Z2, arrow_quiver).unwrap();
    let arrow = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let nonadmissible = QuiverAlgebra::Quotient(
        quotient_by_relations(parent, std::slice::from_ref(&arrow)).unwrap(),
    );
    assert_eq!(
        nakayama_permutation(&nonadmissible).unwrap_err(),
        EnvelopingModuleError::NakayamaPermutationSelfinjectivityUnavailable
    );
}

fn radical_square_zero_cycle<F: Field>(field: F) -> QuiverAlgebra<F> {
    let quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "x"),
            ArrowSpec::named(vertex(2), vertex(3), "y"),
            ArrowSpec::named(vertex(3), vertex(1), "z"),
        ],
    )
    .unwrap();
    let parent = PathAlgebra::new(field, quiver).unwrap();
    let arrows = parent.quiver().arrows().collect::<Vec<_>>();
    let x = parent.element(arrows[0].path()).unwrap();
    let y = parent.element(arrows[1].path()).unwrap();
    let z = parent.element(arrows[2].path()).unwrap();
    QuiverAlgebra::Quotient(
        quotient_by_relations(
            parent,
            &[
                x.multiply(&y).unwrap(),
                y.multiply(&z).unwrap(),
                z.multiply(&x).unwrap(),
            ],
        )
        .unwrap(),
    )
}

fn quantum_complete_intersection<F: Field>(field: F, scalar: &F::Element) -> QuiverAlgebra<F> {
    let quiver = Quiver::from_indexed(
        1,
        vec![
            ArrowSpec::named(vertex(1), vertex(1), "x"),
            ArrowSpec::named(vertex(1), vertex(1), "y"),
        ],
    )
    .unwrap();
    let parent = PathAlgebra::new(field, quiver).unwrap();
    let arrows = parent.quiver().arrows().collect::<Vec<_>>();
    let x = parent.element(arrows[0].path()).unwrap();
    let y = parent.element(arrows[1].path()).unwrap();
    let xy = x.multiply(&y).unwrap();
    let yx = y.multiply(&x).unwrap();
    QuiverAlgebra::Quotient(
        quotient_by_relations(
            parent,
            &[
                x.multiply(&x).unwrap(),
                y.multiply(&y).unwrap(),
                xy.sub(&yx.scale_left(scalar).unwrap()).unwrap(),
            ],
        )
        .unwrap(),
    )
}

fn gap_z2_scalar(value: <Z2 as numerica::domains::Set>::Element) -> String {
    if Z2.is_zero(&value) {
        "0*Z(2)".to_owned()
    } else {
        "Z(2)^0".to_owned()
    }
}

fn gap_zp_scalar(field: &Zp, value: <Zp as numerica::domains::Set>::Element) -> String {
    if field.is_zero(&value) {
        return format!("0*Z({})", field.get_prime());
    }
    let elements = field.qpa_elements().unwrap().unwrap();
    let position = elements
        .iter()
        .position(|candidate| candidate == &value)
        .unwrap();
    let exponent = position - 1;
    match exponent {
        1 => format!("Z({})", field.get_prime()),
        _ => format!("Z({})^{exponent}", field.get_prime()),
    }
}

fn scalar_matrix_json<F: Field>(
    matrix: &DenseMatrix<F>,
    scalar: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    json!(
        (0..matrix.nrows())
            .map(|row| matrix.row(row).iter().map(scalar).collect::<Vec<_>>())
            .collect::<Vec<_>>()
    )
}

fn nakayama_order_json(value: &NakayamaAutomorphismOrder) -> Value {
    match value {
        NakayamaAutomorphismOrder::Finite(order) => {
            serde_json::from_str(&order.to_string()).unwrap()
        }
        NakayamaAutomorphismOrder::Infinite => json!("infinity"),
    }
}

fn nakayama_path_case_json<F: QpaCoefficientField>(
    algebra: &QuiverAlgebra<F>,
    scalar: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let auto_initially_known = algebra.cached_nakayama_automorphism_attribute().is_bound();
    let order_initially_known = algebra
        .cached_order_of_nakayama_automorphism_attribute()
        .is_bound();
    let automorphism = nakayama_automorphism(algebra).unwrap();
    let order = order_of_nakayama_automorphism(algebra).unwrap();
    json!({
        "auto_initially_known": auto_initially_known,
        "order_initially_known": order_initially_known,
        "automorphism": automorphism.as_ref().map_or(Value::Bool(false),
            |value| scalar_matrix_json(value.matrix(), scalar)),
        "order": order.as_ref().map_or(Value::Bool(false), nakayama_order_json),
        "auto_known_after": algebra.cached_nakayama_automorphism_attribute().is_bound(),
        "order_known_after": algebra
            .cached_order_of_nakayama_automorphism_attribute().is_bound()
    })
}

#[derive(Clone, Copy)]
enum FrobeniusTrigger {
    Automorphism,
    Form,
    Functional,
}

fn nakayama_quotient_case_json<F: QpaCoefficientField>(
    algebra: &QuiverAlgebra<F>,
    trigger: FrobeniusTrigger,
    scalar: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let QuiverAlgebra::Quotient(quotient) = algebra else {
        unreachable!()
    };
    let auto_initially_known = algebra.cached_nakayama_automorphism_attribute().is_bound();
    let form_initially_known = quotient.cached_frobenius_form_attribute().is_bound();
    let functional_initially_known = quotient
        .cached_frobenius_linear_functional_attribute()
        .is_bound();
    let order_initially_known = algebra
        .cached_order_of_nakayama_automorphism_attribute()
        .is_bound();

    let (automorphism, mut form, mut functional) = match trigger {
        FrobeniusTrigger::Automorphism => (nakayama_automorphism(algebra).unwrap(), None, None),
        FrobeniusTrigger::Form => {
            let form = frobenius_form(quotient).unwrap();
            let automorphism = nakayama_automorphism(algebra).unwrap();
            (automorphism, form, None)
        }
        FrobeniusTrigger::Functional => {
            let functional = frobenius_linear_functional(quotient).unwrap();
            let automorphism = nakayama_automorphism(algebra).unwrap();
            (automorphism, None, functional)
        }
    };
    let auto_known_after_trigger = algebra.cached_nakayama_automorphism_attribute().is_bound();
    let form_known_after_trigger = quotient.cached_frobenius_form_attribute().is_bound();
    let functional_known_after_trigger = quotient
        .cached_frobenius_linear_functional_attribute()
        .is_bound();

    if form.is_none() {
        form = frobenius_form(quotient).unwrap();
    }
    if functional.is_none() {
        functional = frobenius_linear_functional(quotient).unwrap();
    }
    let mut result = json!({
        "auto_initially_known": auto_initially_known,
        "form_initially_known": form_initially_known,
        "functional_initially_known": functional_initially_known,
        "order_initially_known": order_initially_known,
        "auto_known_after_trigger": auto_known_after_trigger,
        "form_known_after_trigger": form_known_after_trigger,
        "functional_known_after_trigger": functional_known_after_trigger,
        "automorphism": automorphism.as_ref().map_or(Value::Bool(false),
            |value| scalar_matrix_json(value.matrix(), scalar)),
        "form": form.as_ref().map_or(Value::Bool(false),
            |value| scalar_matrix_json(value.matrix(), scalar)),
        "functional": functional.as_ref().map_or(Value::Bool(false), |value| {
            json!(value.canonical_basis_values().iter().map(scalar).collect::<Vec<_>>())
        }),
    });
    if automorphism.is_some() {
        result["order"] =
            nakayama_order_json(&order_of_nakayama_automorphism(algebra).unwrap().unwrap());
    }
    result["auto_known_after"] = json!(algebra.cached_nakayama_automorphism_attribute().is_bound());
    result["form_known_after"] = json!(quotient.cached_frobenius_form_attribute().is_bound());
    result["functional_known_after"] = json!(
        quotient
            .cached_frobenius_linear_functional_attribute()
            .is_bound()
    );
    result["order_known_after"] = json!(
        algebra
            .cached_order_of_nakayama_automorphism_attribute()
            .is_bound()
    );
    result
}

#[allow(clippy::too_many_lines)]
fn pathalg_nakayama_automorphism_cases<F: QpaCoefficientField>(
    field: F,
    quantum_scalar: &F::Element,
    scalar: impl Fn(&F::Element) -> String + Copy,
) -> Value {
    let discrete = isolated_path_algebra(field.clone(), 2);
    let a2_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());

    let loop_quiver =
        Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
    let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
    let loop_arrow = loop_parent
        .element(loop_parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let dual_numbers = QuiverAlgebra::Quotient(
        quotient_by_relations(loop_parent, &[loop_arrow.multiply(&loop_arrow).unwrap()]).unwrap(),
    );

    let a3_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "u"),
            ArrowSpec::named(vertex(2), vertex(3), "v"),
        ],
    )
    .unwrap();
    let a3_parent = PathAlgebra::new(field.clone(), a3_quiver).unwrap();
    let arrows = a3_parent.quiver().arrows().collect::<Vec<_>>();
    let u = a3_parent.element(arrows[0].path()).unwrap();
    let v = a3_parent.element(arrows[1].path()).unwrap();
    let nonselfinjective = QuiverAlgebra::Quotient(
        quotient_by_relations(a3_parent, &[u.multiply(&v).unwrap()]).unwrap(),
    );

    let cycle = radical_square_zero_cycle(field.clone());
    let quantum = quantum_complete_intersection(field, quantum_scalar);

    json!({
        "path_discrete": nakayama_path_case_json(&discrete, scalar),
        "path_a2": nakayama_path_case_json(&a2, scalar),
        "quotient_dual_numbers": nakayama_quotient_case_json(
            &dual_numbers, FrobeniusTrigger::Form, scalar),
        "quotient_nonselfinjective": nakayama_quotient_case_json(
            &nonselfinjective, FrobeniusTrigger::Automorphism, scalar),
        "quotient_radical_square_zero_cycle": nakayama_quotient_case_json(
            &cycle, FrobeniusTrigger::Automorphism, scalar),
        "quotient_quantum_complete_intersection": nakayama_quotient_case_json(
            &quantum, FrobeniusTrigger::Functional, scalar)
    })
}

#[test]
fn pinned_qpa_pathalg_nakayama_automorphism_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/pathalg_nakayama_automorphism.schema-v1.json"
    ))
    .unwrap();
    assert_eq!(
        pathalg_nakayama_automorphism_cases(Z2, &Z2.one(), |value| gap_z2_scalar(*value)),
        expected["f2"]
    );
    let field5 = Zp::new(5);
    assert_eq!(
        pathalg_nakayama_automorphism_cases(field5.clone(), &field5.to_element(2), |value| {
            gap_zp_scalar(&field5, *value)
        },),
        expected["f5"]
    );
    assert_eq!(
        pathalg_nakayama_automorphism_cases(Q, &q(2), ToString::to_string),
        expected["rational"]
    );
}

#[test]
fn pinned_qpa_pathalg_nakayama_automorphism_dispatch_matches() {
    let expected: Value = serde_json::from_str(include_str!(
        "../../../../tests/fixtures/pathalg_nakayama_automorphism_dispatch.schema-v1.json"
    ))
    .unwrap();
    let cycle = radical_square_zero_cycle(Z2);
    let QuiverAlgebra::Quotient(cycle_quotient) = &cycle else {
        unreachable!()
    };
    let automorphism = nakayama_automorphism(&cycle).unwrap().unwrap();
    let form = frobenius_form(cycle_quotient).unwrap().unwrap();
    let functional = frobenius_linear_functional(cycle_quotient)
        .unwrap()
        .unwrap();
    let foreign = isolated_path_algebra(Z2, 1).one().unwrap();
    let one = cycle.one().unwrap();

    let a3_quiver = Quiver::from_indexed(
        3,
        vec![
            ArrowSpec::named(vertex(1), vertex(2), "u"),
            ArrowSpec::named(vertex(2), vertex(3), "v"),
        ],
    )
    .unwrap();
    let a3_parent = PathAlgebra::new(Z2, a3_quiver).unwrap();
    let arrows = a3_parent.quiver().arrows().collect::<Vec<_>>();
    let u = a3_parent.element(arrows[0].path()).unwrap();
    let v = a3_parent.element(arrows[1].path()).unwrap();
    let nonselfinjective = QuiverAlgebra::Quotient(
        quotient_by_relations(a3_parent, &[u.multiply(&v).unwrap()]).unwrap(),
    );

    let arrow_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let parent = PathAlgebra::new(Z2, arrow_quiver).unwrap();
    let arrow = parent
        .element(parent.quiver().arrows().next().unwrap().path())
        .unwrap();
    let nonadmissible = QuiverAlgebra::Quotient(
        quotient_by_relations(parent, std::slice::from_ref(&arrow)).unwrap(),
    );
    let QuiverAlgebra::Quotient(nonself_quotient) = &nonselfinjective else {
        unreachable!()
    };
    let QuiverAlgebra::Quotient(nonadmissible_quotient) = &nonadmissible else {
        unreachable!()
    };

    let actual = json!({
        "schema_version": 1,
        "source": {
            "gap": "4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9",
            "qpa": "1.37@8d5f847258d4853700ce3a5d947db2a1441109f1"
        },
        "foreign_auto_success": automorphism.image(&foreign).is_ok(),
        "foreign_form_left_success": form.evaluate(&foreign, &one).is_ok(),
        "foreign_form_right_success": form.evaluate(&one, &foreign).is_ok(),
        "foreign_functional_success": functional.evaluate(&foreign).is_ok(),
        "nonself_auto_is_false": nakayama_automorphism(&nonselfinjective)
            .unwrap().is_none(),
        "nonself_form_is_false": frobenius_form(nonself_quotient).unwrap().is_none(),
        "nonself_functional_is_false": frobenius_linear_functional(nonself_quotient)
            .unwrap().is_none(),
        "nonself_order_success": order_of_nakayama_automorphism(&nonselfinjective).is_ok(),
        "nonadmissible_auto_success": nakayama_automorphism(&nonadmissible).is_ok(),
        "nonadmissible_form_success": frobenius_form(nonadmissible_quotient).is_ok(),
        "nonadmissible_functional_success":
            frobenius_linear_functional(nonadmissible_quotient).is_ok(),
        "nonadmissible_order_success":
            order_of_nakayama_automorphism(&nonadmissible).is_ok()
    });
    assert_eq!(actual, expected);
    assert!(matches!(
        automorphism.image(&foreign),
        Err(EnvelopingModuleError::NakayamaAutomorphismElementOutsideAlgebra)
    ));
    assert_eq!(
        order_of_nakayama_automorphism(&nonselfinjective),
        Err(EnvelopingModuleError::CalledFalseNakayamaAutomorphism)
    );
    assert!(matches!(
        nakayama_automorphism(&nonadmissible),
        Err(EnvelopingModuleError::NakayamaAutomorphismSelfinjectivityUnavailable)
    ));
}

#[test]
fn nakayama_automorphism_frobenius_data_and_order_follow_qpa() {
    let cycle = radical_square_zero_cycle(Q);
    let QuiverAlgebra::Quotient(cycle_quotient) = &cycle else {
        unreachable!()
    };
    assert!(!cycle.cached_nakayama_automorphism_attribute().is_bound());
    assert!(!cycle_quotient.cached_frobenius_form_attribute().is_bound());
    let form = frobenius_form(cycle_quotient).unwrap().unwrap();
    assert!(cycle.cached_nakayama_automorphism_attribute().is_bound());
    assert!(
        cycle_quotient
            .cached_frobenius_linear_functional_attribute()
            .is_bound()
    );
    let automorphism = nakayama_automorphism(&cycle).unwrap().unwrap();
    let expected_automorphism = DenseMatrix::from_rows(
        vec![
            vec![q(0), q(1), q(0), q(0), q(0), q(0)],
            vec![q(0), q(0), q(1), q(0), q(0), q(0)],
            vec![q(1), q(0), q(0), q(0), q(0), q(0)],
            vec![q(0), q(0), q(0), q(0), q(1), q(0)],
            vec![q(0), q(0), q(0), q(0), q(0), q(1)],
            vec![q(0), q(0), q(0), q(1), q(0), q(0)],
        ],
        Q,
    )
    .unwrap();
    assert_eq!(automorphism.matrix(), &expected_automorphism);
    assert_eq!(
        order_of_nakayama_automorphism(&cycle).unwrap(),
        Some(NakayamaAutomorphismOrder::Finite(Integer::from(3)))
    );
    let functional = frobenius_linear_functional(cycle_quotient)
        .unwrap()
        .unwrap();
    assert_eq!(
        functional.canonical_basis_values(),
        &[q(0), q(0), q(0), q(1), q(1), q(1)]
    );
    let basis = finite_canonical_basis(&cycle).unwrap().unwrap();
    for (row, left) in basis.iter().enumerate() {
        for (column, right) in basis.iter().enumerate() {
            assert_eq!(
                form.evaluate(left, right).unwrap(),
                form.matrix()[(row, column)]
            );
        }
        assert_eq!(
            functional.evaluate(left).unwrap(),
            functional.canonical_basis_values()[row]
        );
        assert_eq!(
            canonical_basis_coefficients(&cycle, &basis, &automorphism.image(left).unwrap())
                .unwrap(),
            expected_automorphism.row(row)
        );
    }

    let quantum_q = quantum_complete_intersection(Q, &q(2));
    let quantum_q_auto = nakayama_automorphism(&quantum_q).unwrap().unwrap();
    assert_eq!(
        quantum_q_auto.matrix(),
        &DenseMatrix::from_rows(
            vec![
                vec![q(1), q(0), q(0), q(0)],
                vec![q(0), q(2), q(0), q(0)],
                vec![q(0), q(0), Q.to_element(1.into(), 2.into(), true), q(0)],
                vec![q(0), q(0), q(0), q(1)],
            ],
            Q,
        )
        .unwrap()
    );
    assert_eq!(
        order_of_nakayama_automorphism(&quantum_q).unwrap(),
        Some(NakayamaAutomorphismOrder::Infinite)
    );

    let field5 = Zp::new(5);
    let quantum_f5 = quantum_complete_intersection(field5.clone(), &field5.to_element(2));
    assert_eq!(
        order_of_nakayama_automorphism(&quantum_f5).unwrap(),
        Some(NakayamaAutomorphismOrder::Finite(Integer::from(4)))
    );
    let quantum_f2 = quantum_complete_intersection(Z2, &Z2.one());
    assert_eq!(
        order_of_nakayama_automorphism(&quantum_f2).unwrap(),
        Some(NakayamaAutomorphismOrder::Finite(Integer::one()))
    );

    let discrete = isolated_path_algebra(Q, 2);
    assert_eq!(
        order_of_nakayama_automorphism(&discrete).unwrap(),
        Some(NakayamaAutomorphismOrder::Finite(Integer::one()))
    );
    let a2_quiver =
        Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
    let a2 = QuiverAlgebra::Path(PathAlgebra::new(Q, a2_quiver).unwrap());
    assert!(nakayama_automorphism(&a2).unwrap().is_none());
    assert!(order_of_nakayama_automorphism(&a2).unwrap().is_none());
}
