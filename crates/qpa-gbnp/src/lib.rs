//! The transitive GBNP implementation surface used by QPA.
//!
//! The translation is rooted at GBNP's `SGrobner`; no independent Gröbner
//! algorithm is substituted here.

pub mod heap;
pub mod invlenlex;
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
pub mod sgrobner;

pub use heap::{TodoHeap, TodoHeapError};
pub use invlenlex::{
    InvLenLexError, biggest_homogeneous_part_inv_len_inv_lex,
    biggest_homogeneous_parts_inv_len_inv_lex, groebner_basis_inv_len_inv_lex,
    reduced_list_inv_len_inv_lex, remainder_of_division_inv_len_inv_lex,
    tip_coefficient_inv_len_inv_lex, tip_inv_len_inv_lex, tip_monomial_inv_len_inv_lex,
    tip_reduced_list_inv_len_inv_lex, tip_walk_inv_len_inv_lex, truncate_element,
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
    gbnp_groebner_basis_nc, path_to_cohen, path_to_cohen_free, relations_for_path_algebra,
};
pub use reduction::{
    ReductionError, left_occurrence, normal_form_np, occur, occur_in_list, right_occurrence,
    strong_normal_form_np, strong_normal_form_two,
};
pub use sgrobner::{
    SgrobnerError, SgrobnerResult, all_obstructions, sgrobner, sgrobner_loops,
    sgrobner_with_options, strong_normal_form_trees,
};
