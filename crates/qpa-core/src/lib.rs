//! Exact foundational structures for the source-faithful Rust port of QPA.
//!
//! QPA and the GAP primitives it calls define observable semantics. External
//! crates provide exact arithmetic and storage primitives, not replacement QPA
//! algorithms.

pub mod centre;
pub mod coefficient;
pub mod dynamic_tree;
pub mod element_utilities;
pub mod groebner;
pub mod matrix;
pub mod ordering;
pub mod path_algebra;
pub mod poset;
pub mod quiver;
pub mod quotient;
pub mod static_dictionary;
pub mod suffix_tree;
pub mod vector_reduction;

pub use centre::{CentreError, PathAlgebraCentre, QuotientPathAlgebraCentre};
pub use coefficient::{CoefficientDomain, CoefficientDomainError};
pub use dynamic_tree::{DynamicForest, DynamicTree, DynamicTreeError};
pub use element_utilities::{QuiverAlgebra, VertexQuotientPresentation};
pub use groebner::{
    GroebnerBasis, GroebnerBasisKind, GroebnerError, high_level_groebner_basis,
    leftmost_occurrence, make_uniform, qpa_in_arrow_ideal, reduced_list_qpa, remainder_of_division,
    tip_reduced_list,
};
pub use matrix::{DenseMatrix, MatrixError, SemiEchelonTransformation};
pub use ordering::{OrderingError, QuiverOrdering, check_quiver_subset_for_cycles};
pub use path_algebra::{
    PathAlgebra, PathAlgebraElement, PathAlgebraError, PathAlgebraExternalRepresentation,
    PathAlgebraTerm, opposite_relations,
};
pub use poset::{Poset, PosetError};
pub use quiver::{
    Arrow, ArrowIndex, ArrowSpec, DynkinFamily, DynkinOrientation, DynkinType, Generator,
    NamedArrowSpec, Path, PathEndpoint, PathError, PathIterator, Quiver, QuiverCardinality,
    QuiverError, Vertex, VertexIndex,
};
pub use quotient::{
    AlgebraDimension, NaturalQuotientMap, PathAlgebraIdeal, QuotientCanonicalBasis, QuotientError,
    QuotientNiceBasis, QuotientPathAlgebra, QuotientPathAlgebraElement,
    is_homogeneous_list_of_elements,
};
pub use static_dictionary::{DictionaryMatch, QuiverStaticDictionary, StaticDictionaryError};
pub use suffix_tree::{
    SuffixTree, SuffixTreeEdgeId, SuffixTreeError, SuffixTreeNodeId, SuffixTreeNodeKind,
    SuffixTreeOccurrence, SuffixTreeToken,
};
pub use vector_reduction::{VectorCoefficients, VectorReductionError, VectorTip};
