//! The transitive GBNP implementation surface used by QPA.
//!
//! The translation is rooted at GBNP's `SGrobner`; no independent Gröbner
//! algorithm is substituted here.

pub mod np;
pub mod obstruction;
pub mod reduction;

pub use np::{
    NpError, NpMonomial, NpPolynomial, NpTerm, add_np, bimul_np, clean_np, gt_np, lmon_np, lt_np,
    mk_monic_np, mul_np,
};
pub use obstruction::{
    Obstruction, ObstructionError, left_obstructions, right_obstructions, self_obstruction,
    spolynomial,
};
pub use reduction::{
    ReductionError, left_occurrence, normal_form_np, occur, occur_in_list, right_occurrence,
    strong_normal_form_np, strong_normal_form_two,
};
