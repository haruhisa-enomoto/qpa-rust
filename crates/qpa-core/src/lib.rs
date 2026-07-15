//! Exact foundational structures for the source-faithful Rust port of QPA.
//!
//! QPA and the GAP primitives it calls define observable semantics. External
//! crates provide exact arithmetic and storage primitives, not replacement QPA
//! algorithms.

pub mod coefficient;
pub mod matrix;
pub mod ordering;
pub mod quiver;

pub use coefficient::{CoefficientDomain, CoefficientDomainError};
pub use matrix::{DenseMatrix, MatrixError, SemiEchelonTransformation};
pub use ordering::{OrderingError, QuiverOrdering};
pub use quiver::{
    Arrow, ArrowIndex, ArrowSpec, DynkinFamily, DynkinOrientation, DynkinType, Generator,
    NamedArrowSpec, Path, PathEndpoint, PathError, PathIterator, Quiver, QuiverCardinality,
    QuiverError, Vertex, VertexIndex,
};
