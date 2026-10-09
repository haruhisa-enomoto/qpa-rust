//! Supported user-facing facade for foundational QPA computations.
//!
//! The crate root continues to expose the complete source-translation surface
//! for compatibility and source auditing. This module is the smaller facade
//! intended for ordinary Rust programs and receives the stronger stability
//! commitment. Adding more translated QPA operations here is non-breaking;
//! removing or renaming an item follows the crate's versioning policy.
//!
//! # Example
//!
//! ```
//! # use std::error::Error;
//! # fn main() -> Result<(), Box<dyn Error>> {
//! use qpa_core::api::prelude::*;
//!
//! let v1 = VertexIndex::new(1).expect("one-based vertex");
//! let v2 = VertexIndex::new(2).expect("one-based vertex");
//! let quiver = Quiver::from_indexed(
//!     2,
//!     vec![ArrowSpec::named(v1, v2, "a")],
//! )?;
//! let path_algebra = PathAlgebra::new(Q, quiver)?;
//!
//! assert_eq!(path_algebra.quiver().number_of_vertices(), 2);
//! # Ok(())
//! # }
//! ```

/// Exact coefficient fields and scalar operations supplied by Numerica.
pub mod coefficients {
    pub use numerica::domains::finite_field::{FiniteFieldCore, Z2, Zp};
    pub use numerica::domains::integer::Integer;
    pub use numerica::domains::rational::{Q, Rational, RationalField};
    pub use numerica::domains::{Field, Ring, RingOps};

    pub use crate::{
        CoefficientDomain, CoefficientDomainError, FiniteFieldConstruction, GaloisField,
        GaloisFieldElement, GaloisFieldError, QpaCoefficientField, QpaFieldCharacteristic,
    };
}

/// Quivers, one-based indices, paths, and QPA path-ordering data.
pub mod quivers {
    pub use crate::{
        Arrow, ArrowIndex, ArrowSpec, DynkinDecision, DynkinFamily, DynkinOrientation, DynkinType,
        Generator, NamedArrowSpec, NextPathPosition, OrderedQuiver, OrderingError, Path,
        PathEndpoint, PathError, PathIterator, Quiver, QuiverCardinality, QuiverError,
        QuiverOrdering, QuiverOrderingKind, Vertex, VertexIndex,
    };
}

/// Path algebras, ideals, quotients, exact elements, and Gröbner data.
pub mod algebras {
    pub use crate::{
        AlgebraDimension, AlgebraGlobalDimension, GroebnerBasis, GroebnerBasisKind, GroebnerError,
        NaturalQuotientMap, PathAlgebra, PathAlgebraElement, PathAlgebraError,
        PathAlgebraExternalRepresentation, PathAlgebraIdeal, PathAlgebraTerm, QuiverAlgebra,
        QuiverAlgebraElement, QuiverAlgebraElementError, QuotientError, QuotientPathAlgebra,
        QuotientPathAlgebraElement, high_level_groebner_basis,
    };
}

/// Exact dense matrices using QPA's row-action convention.
pub mod linear_algebra {
    pub use crate::{DenseMatrix, MatrixError, SemiEchelonTransformation};
}

/// Right matrix modules, module maps, and commonly used constructions.
pub mod modules {
    pub use crate::{
        ModuleDimensionSeries, NamedArrowMatrix, PathAlgebraModule, PathAlgebraModuleArrowAction,
        PathAlgebraModuleArrowActions, PathAlgebraModuleElement, PathAlgebraModuleError,
        PathAlgebraModuleHomomorphism, PathAlgebraModuleMatrix, cokernel_of_module_homomorphism,
        direct_sum_inclusions, direct_sum_of_qpa_modules, direct_sum_projections, end_over_algebra,
        hom_dimension_over_algebra, hom_over_algebra, identity_module_mapping, image_inclusion,
        image_of_module_homomorphism, image_projection, kernel_inclusion,
        kernel_of_module_homomorphism, radical_of_module, right_module_hom_over_algebra,
        right_module_over_path_algebra, right_module_over_path_algebra_infer_dimension,
        socle_of_module, top_of_module, zero_module, zero_module_mapping,
    };
}

/// Convenient imports for a small foundational QPA program.
///
/// This prelude is deliberately compact. Import the domain modules above when
/// a program needs less common types or wants imports to remain explicit.
pub mod prelude {
    pub use numerica::domains::finite_field::{FiniteFieldCore, Z2, Zp};
    pub use numerica::domains::rational::{Q, Rational, RationalField};
    pub use numerica::domains::{Field, Ring, RingOps};

    pub use crate::{
        AlgebraDimension, ArrowIndex, ArrowSpec, DenseMatrix, GaloisField, GaloisFieldElement,
        NamedArrowMatrix, Path, PathAlgebra, PathAlgebraElement, PathAlgebraIdeal,
        PathAlgebraModule, PathAlgebraModuleError, PathAlgebraModuleHomomorphism, Quiver,
        QuiverAlgebra, QuotientPathAlgebra, VertexIndex, right_module_hom_over_algebra,
        right_module_over_path_algebra, zero_module,
    };
}
