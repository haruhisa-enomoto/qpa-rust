//! QPA modules constructed from enveloping algebras.
//!
//! This begins the source-order translation of `lib/pathalgtensor.gi`'s
//! module section. The regular bimodule method deliberately uses QPA's
//! canonical-basis partition and its checked path-module constructor.

use std::num::NonZeroUsize;

use numerica::domains::{
    Field, Ring, RingOps,
    finite_field::{Z2, Zp},
    integer::Integer,
    rational::{Q, Rational, RationalField},
};
use qpa_core::{
    AlgebraGlobalDimension, CachedAttribute, CachedNakayamaAutomorphismOrder,
    CoefficientDomainError, CommonDirectSummandResult, DenseMatrix, MatrixError,
    ModuleDecompositionError, ModuleDimensionSeries, ModuleEndomorphismQuotientHomomorphism,
    NamedArrowMatrix, PathAlgebraElement, PathAlgebraError, PathAlgebraModule,
    PathAlgebraModuleElement, PathAlgebraModuleError, PathAlgebraModuleHomomorphism,
    PathAlgebraModuleHomomorphismPair, PathAlgebraTensorError, PathEndpoint, PathError,
    QpaCoefficientField, QpaFieldCharacteristic, QuiverAlgebra, QuiverAlgebraElement,
    QuiverAlgebraElementError, QuiverProduct, QuiverProductError, QuotientError,
    QuotientPathAlgebra, StructureConstantsError, cokernel_of_module_homomorphism,
    cokernel_projection, common_direct_summand, decompose_module_with_inclusions,
    decompose_module_with_multiplicities, direct_sum_inclusions, direct_sum_of_qpa_modules,
    direct_sum_projections, end_over_algebra, from_end_m_to_hom_mm, from_hom_mm_to_end_m,
    hom_over_algebra, identity_module_mapping, image_projection, include_in_path_algebra,
    is_direct_summand, is_indecomposable_module, isomorphic_modules, isomorphism_of_modules,
    kernel_inclusion, kernel_of_module_homomorphism, left_minimal_version, loewy_length_of_module,
    minimal_generating_set_of_module, morphism_on_cokernel, morphism_on_kernel,
    natural_endomorphism_quotient_by_ideal, project_from_product_quiver, pull_back, push_out,
    radical_of_module, radical_of_module_inclusion, right_minimal_version,
    right_module_hom_over_algebra, right_module_over_path_algebra,
    right_module_over_path_algebra_with_tensor, socle_of_module, socle_of_module_inclusion,
    subrepresentation, subrepresentation_inclusion, top_of_module, zero_module,
    zero_module_mapping,
};
use rand::Rng;
use thiserror::Error;

use crate::complex::projective_resolution;
use crate::{
    PathalgGbnpError, TensorProductError, ensure_quotient_groebner_basis, enveloping_algebra,
    opposite_path_algebra, opposite_quotient_element, quotient_is_admissible,
    quotient_is_finite_dimensional,
};

// These private source chunks remain in this module through `include!`.
// That preserves every established `qpa_gbnp::module::*` path while
// making the QPA source domains reviewable in isolation.
include!("module/algebra_invariants.rs");
include!("module/projectives.rs");
include!("module/syzygies_and_tilting.rs");
include!("module/approximations.rs");
include!("module/ext_and_iyama.rs");
include!("module/auslander_reiten.rs");
include!("module/stable_and_functors.rs");
include!("module/tensor_and_enveloping.rs");
include!("module/errors.rs");

#[cfg(test)]
#[path = "module/tests.rs"]
mod tests;
