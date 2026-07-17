//! The transitive GBNP implementation surface used by QPA.
//!
//! The translation is rooted at GBNP's `SGrobner`; no independent Gröbner
//! algorithm is substituted here.

pub mod algebra_conversion;
pub mod complex;
pub mod heap;
pub mod higher_preprojective;
pub mod invlenlex;
pub mod module;
pub mod np;
pub mod obstruction;
pub mod occur_tree;
pub mod options;
pub mod pathalg;
pub mod persistence;
pub mod predefined;
pub mod preprocessing;
pub mod qpa_conversion;
pub mod reduction;
pub mod reptype;
pub mod sgrobner;
pub mod tensor;
pub mod trivial_extension;

pub use algebra_conversion::{
    AdrAlgebraError, AdrAlgebraOfAlgebraResult, AlgebraAsQuiverAlgebraError,
    AlgebraAsQuiverAlgebraResult, EndOfBasicModuleAsQuiverAlgebraError,
    EndOfBasicModuleAsQuiverAlgebraResult, EndOfModuleAsQuiverAlgebraError,
    EndOfModuleAsQuiverAlgebraResult, adr_algebra_of_algebra, adr_algebra_of_basic_module,
    algebra_as_quiver_algebra, end_of_basic_module_as_quiver_algebra,
    end_of_basic_module_as_quiver_algebra_with_bound, end_of_module_as_quiver_algebra,
};
pub use complex::{
    InjectiveResolution, ProjectiveResolution, injective_resolution, projective_resolution,
};
pub use heap::{TodoHeap, TodoHeapError};
pub use higher_preprojective::{
    ModulePreprojectiveError, ModulePreprojectiveStructure, module_preprojective_algebra,
    module_preprojective_structure,
};
pub use invlenlex::{
    InvLenLexError, biggest_homogeneous_part_inv_len_inv_lex,
    biggest_homogeneous_parts_inv_len_inv_lex, groebner_basis_inv_len_inv_lex,
    leftmost_occurrence_inv_len_inv_lex, reduced_list_inv_len_inv_lex,
    remainder_of_division_inv_len_inv_lex, tip_coefficient_inv_len_inv_lex, tip_inv_len_inv_lex,
    tip_monomial_inv_len_inv_lex, tip_reduced_list_inv_len_inv_lex, tip_walk_inv_len_inv_lex,
    truncate_element,
};
pub use module::{
    AddMCoresolution, AddMResolution, AllComplementsResult, AlmostSplitSide, AuslanderReitenError,
    CotiltingModuleResult, DominantDimensionResult, EnvelopingModuleError, ExceptionalModuleError,
    ExtOverAlgebraResult, FaithfulDimensionResult, FiniteAddMResolutionResult, FrobeniusForm,
    FrobeniusLinearFunctional, GlobalDimensionOfAlgebraResult, ModuleFunctorError,
    ModuleLengthError, NakayamaAutomorphism, NakayamaAutomorphismOrder, NakayamaPermutation,
    NumberOfComplementsResult, PredecessorLayers, PredecessorMultiplicity, PredecessorValuation,
    PredecessorsOfModuleResult, ProjectiveModuleBases, QpaRandomScalar,
    TensorProductOfModulesResult, TiltingModuleResult, algebra_as_module_of_enveloping_algebra,
    algebra_as_module_over_enveloping_algebra, all_complements_of_almost_complete_cotilting_module,
    all_complements_of_almost_complete_tilting_module, all_indec_modules_of_length_at_most,
    all_modules_of_length_at_most, all_modules_of_length_plus_one, almost_split_sequence,
    almost_split_sequence_from_side, annihilator_of_module, basis_of_projectives, cartan_matrix,
    complexity_of_algebra, complexity_of_module, cotilting_module, coxeter_matrix,
    coxeter_polynomial, d_tr, d_tr_homomorphism, d_tr_homomorphism_iterated, d_tr_iterated,
    dominant_dimension_of_algebra, dominant_dimension_of_module,
    dual_of_algebra_as_module_over_enveloping_algebra, dual_of_module, dual_of_module_homomorphism,
    element_in_indec_projective, element_in_v_a_as_element_in_indec_proj, end_modulo_projectives,
    ext_over_algebra, faithful_dimension, first_syzygy, frobenius_form,
    frobenius_linear_functional, from_identity_to_double_star_homomorphism,
    from_matrix_to_homomorphism_of_projectives, global_dimension_from_known_properties,
    global_dimension_of_algebra, gorenstein_dimension, gorenstein_dimension_of_algebra,
    gorenstein_dimension_of_algebra_with_progress, have_finite_coresolution_in_add_m,
    have_finite_resolution_in_add_m, hom_factoring_through_projectives, hom_from_projective,
    indec_injective_modules, indec_projective_modules, inj_dimension, inj_dimension_of_module,
    injective_envelope, irreducible_morphisms_ending_in, irreducible_morphisms_starting_in,
    is_cotilting_module, is_exceptional_module, is_finite_global_dimension_from_known_properties,
    is_hereditary_algebra, is_injective_module, is_omega_periodic, is_omega_periodic_with_progress,
    is_projective_module, is_rigid_module, is_schurian_algebra, is_selfinjective_algebra,
    is_semicommutative_algebra, is_symmetric_algebra, is_tau_periodic, is_tau_rigid_module,
    is_tilting_module, is_weakly_symmetric_algebra, left_approximation_by_add_m,
    left_fac_m_approximation, left_sub_m_approximation, lifting_morphism_from_projective,
    loewy_length_of_algebra, matrix_of_homomorphism_between_projectives,
    minimal_left_approximation, minimal_left_fac_m_approximation, minimal_left_sub_m_approximation,
    minimal_right_approximation, minimal_right_sub_m_approximation, n_rigid_module,
    nakayama_automorphism, nakayama_functor_of_module, nakayama_functor_of_module_homomorphism,
    nakayama_permutation, nth_syzygy, number_of_complements_of_almost_complete_cotilting_module,
    number_of_complements_of_almost_complete_tilting_module, opposite_nakayama_functor_of_module,
    opposite_nakayama_functor_of_module_homomorphism, order_of_nakayama_automorphism,
    predecessors_of_module, proj_dimension, proj_dimension_of_module, projective_cover,
    random_module, right_approximation_by_add_m, right_sub_m_approximation, simple_modules,
    socle_series, star_of_module, star_of_module_homomorphism, tensor_product_of_modules,
    tilting_module, tr_d, tr_d_homomorphism, tr_d_homomorphism_iterated, tr_d_iterated,
    transpose_of_module, transpose_of_module_homomorphism,
};
pub use module::{
    minimal_left_approximation as minimal_left_add_m_approximation,
    minimal_right_approximation as minimal_right_add_m_approximation,
};
pub use np::{
    NpError, NpMonomial, NpPolynomial, NpTerm, add_np, bimul_np, clean_np, gt_np, lmon_np, lt_np,
    mk_monic_np, mul_np,
};
pub use obstruction::{
    Obstruction, ObstructionError, TreeObstructions, left_obstructions, left_obstructions_tree,
    right_obstructions, right_obstructions_tree, self_obstruction, spolynomial,
};
pub use occur_tree::{OccurTree, OccurTreeError, TreeDirection};
pub use options::GbnpOptions;
pub use pathalg::{
    PathalgGbnpError, admissible_quotient_is_basic, admissible_quotient_is_elementary,
    admissible_quotient_is_semisimple, associated_graded_path_algebra, associated_graded_quotient,
    associated_monomial_path_algebra, associated_monomial_quotient, ensure_quotient_groebner_basis,
    gentle_algebra_is_gorenstein, gentle_algebra_is_special_biserial, ideal_minimal_generating_set,
    inverse_quotient_element, is_special_biserial_quiver, opposite_path_algebra,
    opposite_quotient_element, path_algebra_is_gentle, path_algebra_is_monomial,
    path_algebra_is_special_biserial, path_algebra_is_string, power_quotient_element,
    quiver_algebra_of_path_algebra_mod_aea, quiver_algebra_of_quotient_mod_aea,
    quotient_by_relations, quotient_centre, quotient_dimension, quotient_is_admissible,
    quotient_is_finite_dimensional, quotient_is_gentle, quotient_is_monomial,
    quotient_is_special_biserial, quotient_is_string, quotient_is_triangular_reduced,
    truncated_path_algebra,
};
pub use persistence::{
    LoadedQpaAlgebra, QpaPersistenceError, QpaTextField, parse_qpa_algebra_text, qpa_algebra_text,
    read_qpa_algebra_file, save_qpa_algebra_file,
};
pub use predefined::{
    BrauerConfiguration, BrauerConfigurationAlgebra, BrauerConfigurationAlgebraResult,
    BrauerConfigurationPolygon, BrauerConfigurationValidationError, BrauerConfigurationVertex,
    INVALID_ADMISSIBLE_SEQUENCE_DIAGNOSTIC, NakayamaAlgebraResult,
    PREPROJECTIVE_ORIENTED_CYCLE_DIAGNOSTIC, PosetAlgebra, PredefinedAlgebra,
    PredefinedAlgebraError, PredefinedAlgebraKind, admissible_sequence_generator,
    brauer_configuration_algebra, canonical_algebra, canonical_algebra_two_weights,
    finite_global_dim_green, finite_global_dim_kirk_kuz, kronecker_algebra, nakayama_algebra,
    nakayama_algebra_sequence_first, poset_algebra, poset_algebra_with_string_function,
    preprojective_algebra,
};
pub use preprocessing::{
    PreprocessingError, calculate_prefix_generators, reduce_cancellation, reduce_polynomials,
    reduce_polynomials_in_place, reduce_polynomials_with_options, strong_normal_form_except,
};
pub use qpa_conversion::{
    CohenPolynomial, CohenRepresentation, QpaGbnpError, cohen_to_path, gbnp_groebner_basis,
    gbnp_groebner_basis_nc, gbnp_groebner_basis_of_ideal, path_to_cohen, path_to_cohen_free,
    relations_for_path_algebra,
};
pub use reduction::{
    ReductionError, left_occurrence, normal_form_np, occur, occur_in_list, right_occurrence,
    strong_normal_form_np, strong_normal_form_two,
};
pub use reptype::{
    BONGARTZ_ENTERING_MESSAGE, BongartzTestDecision, DYNKIN_ORIENTED_CYCLE_MESSAGE,
    FINITE_TYPE_DYNKIN_MESSAGE, FINITE_TYPE_INCONCLUSIVE_MESSAGE, FINITE_TYPE_SEMISIMPLE_MESSAGE,
    FINITE_TYPE_SPECIAL_BISERIAL_MESSAGE, INFINITE_TYPE_NON_DYNKIN_MESSAGE,
    INFINITE_TYPE_SEPARATED_QUIVER_MESSAGE, PathAlgebraFiniteTypeDecision,
    PathAlgebraFiniteTypeReason, QuotientAlgebraFiniteTypeDecision,
    QuotientAlgebraFiniteTypeReason, RepresentationTypeError, bongartz_test,
    path_algebra_finite_type, path_algebra_is_distributive, quotient_algebra_finite_type,
    quotient_is_distributive,
};
pub use sgrobner::{
    SgrobnerError, SgrobnerResult, all_obstructions, sgrobner, sgrobner_loops,
    sgrobner_with_options, strong_normal_form_trees,
};
pub use tensor::{
    EMPTY_TENSOR_RELATIONS_ERROR, TensorProductAlgebra, TensorProductError, enveloping_algebra,
    is_enveloping_algebra, tensor_algebra_inclusion, tensor_product_of_algebras,
    tensor_product_of_path_algebras,
};
pub use trivial_extension::{
    TrivialExtensionError, enveloping_algebra_homomorphism, opposite_algebra_homomorphism,
    trivial_extension_level, trivial_extension_of_quiver_algebra,
    trivial_extension_of_quiver_algebra_projection,
};
