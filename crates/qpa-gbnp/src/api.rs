//! Supported user-facing facade for QPA operations that use the GBNP layer.
//!
//! The crate root remains the complete compatibility and source-audit surface.
//! This facade groups the operations most useful to application code while
//! preserving their source-faithful implementations.
//!
//! # Example
//!
//! ```
//! # use std::error::Error;
//! # fn main() -> Result<(), Box<dyn Error>> {
//! use qpa_gbnp::api::prelude::*;
//!
//! let vertex = VertexIndex::new(1).expect("one-based vertex");
//! let quiver = Quiver::from_indexed(
//!     1,
//!     vec![ArrowSpec::named(vertex, vertex, "x")],
//! )?;
//! let path_algebra = PathAlgebra::new(Z2, quiver.clone())?;
//! let x_path = quiver.arrows().next().expect("one loop").path();
//! let x = path_algebra.element(x_path)?;
//! let quotient = quotient_by_relations(path_algebra, &[x.multiply(&x)?])?;
//!
//! assert!(matches!(quotient_dimension(&quotient)?, AlgebraDimension::Finite(_)));
//! # Ok(())
//! # }
//! ```

/// The supported foundational facade from `qpa-core`.
pub use qpa_core::api as core;

/// Quotient construction, Gröbner completion, and common algebra properties.
pub mod quotients {
    pub use crate::{
        PathalgGbnpError, QpaGbnpError, ensure_quotient_groebner_basis, gbnp_groebner_basis,
        gbnp_groebner_basis_nc, gbnp_groebner_basis_of_ideal, opposite_path_algebra,
        opposite_quotient_element, quotient_by_relations, quotient_centre, quotient_dimension,
        quotient_is_admissible, quotient_is_finite_dimensional, quotient_is_gentle,
        quotient_is_monomial, quotient_is_special_biserial, quotient_is_string,
        quotient_is_triangular_reduced, truncated_path_algebra,
    };
}

/// Standard module families, projective/injective maps, and homological data.
pub mod modules {
    pub use crate::{
        DominantDimensionResult, ExtOverAlgebraResult, GlobalDimensionOfAlgebraResult,
        algebra_as_module_over_enveloping_algebra, annihilator_of_module, basis_of_projectives,
        cartan_matrix, coxeter_matrix, coxeter_polynomial, dominant_dimension_of_algebra,
        dominant_dimension_of_module, ext_over_algebra, first_syzygy, global_dimension_of_algebra,
        indec_injective_modules, indec_projective_modules, inj_dimension_of_module,
        injective_envelope, is_injective_module, is_projective_module, loewy_length_of_algebra,
        nth_syzygy, proj_dimension_of_module, projective_cover, simple_modules, socle_series,
    };
}

/// Duality, transpose, Auslander--Reiten, star, and Nakayama functors.
pub mod functors {
    pub use crate::{
        ModuleFunctorError, d_tr, d_tr_homomorphism, d_tr_homomorphism_iterated, d_tr_iterated,
        dual_of_module, dual_of_module_homomorphism, from_identity_to_double_star_homomorphism,
        is_tau_rigid_module, nakayama_functor_of_module, nakayama_functor_of_module_homomorphism,
        opposite_nakayama_functor_of_module, opposite_nakayama_functor_of_module_homomorphism,
        star_of_module, star_of_module_homomorphism, tr_d, tr_d_homomorphism,
        tr_d_homomorphism_iterated, tr_d_iterated, transpose_of_module,
        transpose_of_module_homomorphism,
    };
}

/// Complexes, projective/injective resolutions, and chain maps.
pub mod complexes {
    pub use crate::{
        ChainMapHalf, ComplexAndChainMaps, ComplexAndChainMapsHalf, ComplexError, ComplexHalf,
        InjectiveResolution, ModuleChainMap, ModuleComplex, ProjectiveResolution, QpaChainMap,
        QpaComplex, chain_map, complex, finite_chain_map, finite_complex, injective_resolution,
        mapping_cone, projective_resolution, short_exact_sequence, stalk_complex, yoneda_product,
        zero_chain_map, zero_complex,
    };
}

/// Common imports for quotient and module computations.
pub mod prelude {
    pub use qpa_core::api::prelude::*;

    pub use crate::{
        PathalgGbnpError, ensure_quotient_groebner_basis, indec_projective_modules,
        projective_cover, quotient_by_relations, quotient_dimension,
        quotient_is_finite_dimensional, simple_modules,
    };
}
