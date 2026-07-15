//! Exact foundational structures for the source-faithful Rust port of QPA.
//!
//! QPA and the GAP primitives it calls define observable semantics. External
//! crates provide exact arithmetic and storage primitives, not replacement QPA
//! algorithms.

pub mod coefficient;
pub mod groebner;
pub mod matrix;
pub mod ordering;
pub mod path_algebra;
pub mod quiver;
pub mod quotient;
pub mod static_dictionary;

pub use coefficient::{CoefficientDomain, CoefficientDomainError};
pub use groebner::{
    GroebnerBasis, GroebnerBasisKind, GroebnerError, high_level_groebner_basis,
    leftmost_occurrence, make_uniform, qpa_in_arrow_ideal, reduced_list_qpa, remainder_of_division,
    tip_reduced_list,
};
pub use matrix::{DenseMatrix, MatrixError, SemiEchelonTransformation};
pub use ordering::{OrderingError, QuiverOrdering, check_quiver_subset_for_cycles};
pub use path_algebra::{PathAlgebra, PathAlgebraElement, PathAlgebraError, PathAlgebraTerm};
pub use quiver::{
    Arrow, ArrowIndex, ArrowSpec, DynkinFamily, DynkinOrientation, DynkinType, Generator,
    NamedArrowSpec, Path, PathEndpoint, PathError, PathIterator, Quiver, QuiverCardinality,
    QuiverError, Vertex, VertexIndex,
};
pub use quotient::{
    AlgebraDimension, PathAlgebraIdeal, QuotientError, QuotientPathAlgebra,
    QuotientPathAlgebraElement,
};
pub use static_dictionary::{DictionaryMatch, QuiverStaticDictionary, StaticDictionaryError};
