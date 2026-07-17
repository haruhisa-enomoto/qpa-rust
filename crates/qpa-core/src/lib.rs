//! Exact foundational structures for the source-faithful Rust port of QPA.
//!
//! QPA and the GAP primitives it calls define observable semantics. External
//! crates provide exact arithmetic and storage primitives, not replacement QPA
//! algorithms.

pub mod centre;
pub mod coefficient;
pub mod dynamic_tree;
pub mod element_inverse;
pub mod element_utilities;
pub mod groebner;
pub mod homomorphism;
pub mod idempotent;
pub mod koszul;
pub mod matrix;
pub mod ordering;
pub mod path_algebra;
pub mod path_module;
pub mod poset;
pub mod projective_presentation;
pub mod quiver;
pub mod quotient;
pub mod static_dictionary;
pub mod structure_constants;
pub mod suffix_tree;
pub mod tensor;
pub mod triangular_reduced;
pub mod util;
pub mod vector_reduction;

/// Bound state of a GAP attribute whose computed value may itself be `fail`.
///
/// This distinguishes an attribute that has not been evaluated from one whose
/// source method ran and cached `fail`, without encoding three states as a
/// nested `Option`.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum CachedAttribute<T> {
    /// No attribute value has been computed or installed.
    Unbound,
    /// The source method ran and cached GAP's `fail` value.
    Fail,
    /// The source method ran and cached GAP's boolean `false` sentinel.
    False,
    /// The source method computed and cached a value.
    Value(T),
}

impl<T> CachedAttribute<T> {
    /// Return GAP's `HasAttribute` state.
    #[must_use]
    pub const fn is_bound(&self) -> bool {
        !matches!(self, Self::Unbound)
    }

    /// Return whether the bound value is GAP's `fail`.
    #[must_use]
    pub const fn is_fail(&self) -> bool {
        matches!(self, Self::Fail)
    }

    /// Return whether the bound value is GAP's boolean `false` sentinel.
    #[must_use]
    pub const fn is_false(&self) -> bool {
        matches!(self, Self::False)
    }
}

/// Cache payload for QPA's finite integer or `infinity`
/// `OrderOfNakayamaAutomorphism` value.
///
/// Boolean `false` remains represented by the enclosing optional attribute
/// slot; this enum distinguishes the two successful GAP values without a
/// nested `Option`.
#[doc(hidden)]
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum CachedNakayamaAutomorphismOrder {
    /// A finite positive matrix order.
    Finite(numerica::domains::integer::Integer),
    /// GAP's `infinity` object.
    Infinite,
}

pub use centre::{CentreError, PathAlgebraCentre, QuotientPathAlgebraCentre};
pub use coefficient::{
    CoefficientDomain, CoefficientDomainError, QpaCoefficientField, QpaFieldCharacteristic,
};
pub use dynamic_tree::{DynamicForest, DynamicTree, DynamicTreeError};
pub use element_inverse::ElementInverseError;
pub use element_utilities::{
    QuiverAlgebra, QuiverAlgebraElement, QuiverAlgebraElementError, VertexQuotientPresentation,
};
pub use groebner::{
    GroebnerBasis, GroebnerBasisKind, GroebnerError, INFO_GROEBNER_BASIS,
    high_level_groebner_basis, leftmost_occurrence, make_uniform, make_uniform_list,
    qpa_in_arrow_ideal, reduced_list_qpa, remainder_of_division, tip_reduced_list, tip_walk,
};
pub use homomorphism::{QuiverAlgebraHomomorphism, QuiverAlgebraHomomorphismError};
pub use idempotent::{HOPF_DEFAULT_SEED, HopfError};
pub use koszul::{
    KoszulError, QuadraticPerpendicular, is_quadratic_ideal, paths_of_length_two,
    quadratic_perpendicular_of_path_algebra_ideal,
};
pub use matrix::{DenseMatrix, MatrixError, SemiEchelonTransformation};
pub use ordering::{
    OrderingError, QuiverOrdering, QuiverOrderingKind, check_blocks_for_cycles,
    check_quiver_subset_for_cycles,
};
pub use path_algebra::{
    PathAlgebra, PathAlgebraElement, PathAlgebraError, PathAlgebraExternalRepresentation,
    PathAlgebraTerm, opposite_relations,
};
pub use path_module::{
    CommonDirectSummand, CommonDirectSummandResult, MaximalCommonDirectSummand,
    ModuleDimensionSeries, ModuleEndomorphismAlgebra, ModuleEndomorphismQuotientAlgebra,
    ModuleEndomorphismQuotientElement, ModuleEndomorphismQuotientHomomorphism,
    MoreMinimalVersionResult, NamedArrowMatrix, PathAlgebraModule, PathAlgebraModuleElement,
    PathAlgebraModuleError, PathAlgebraModuleHomomorphism, PathAlgebraModuleHomomorphismBasis,
    PathAlgebraModuleHomomorphismPair, PathAlgebraModuleHomomorphismTriple,
    PathAlgebraModuleMatrix, PathAlgebraModuleMinimalVersion, cokernel_of_module_homomorphism,
    cokernel_projection, common_direct_summand, dimension_vector_partial_order,
    direct_sum_inclusions, direct_sum_of_qpa_modules, direct_sum_projections, end_over_algebra,
    from_end_m_to_hom_mm, from_hom_mm_to_end_m, hom_over_algebra,
    hom_over_algebra_with_basis_function, homomorphism_from_images, identity_module_mapping,
    image_inclusion, image_of_module_homomorphism, image_projection, image_projection_inclusion,
    intersection_of_submodules, intersection_of_submodules_list, inverse_of_isomorphism,
    is_direct_summand, is_in_additive_closure, is_left_minimal, is_right_minimal,
    is_semisimple_module, is_simple_qpa_module, is_split_epimorphism, is_split_monomorphism,
    isomorphic_modules, isomorphism_of_modules, kernel_inclusion, kernel_of_module_homomorphism,
    left_inverse_of_homomorphism, left_minimal_version, lifting_inclusion_morphisms,
    loewy_length_of_module, maximal_common_direct_summand, minimal_generating_set_of_module,
    minimal_right_fac_m_approximation, more_left_minimal_version, more_right_minimal_version,
    morphism_on_cokernel, morphism_on_image, morphism_on_kernel, multiply_lists_of_maps,
    natural_endomorphism_quotient_by_ideal, pull_back, push_out, radical_of_module,
    radical_of_module_inclusion, radical_series, reject_of_module, right_fac_approximation,
    right_fac_m_approximation, right_inverse_of_homomorphism, right_minimal_version,
    right_module_hom_over_algebra, right_module_over_path_algebra,
    right_module_over_path_algebra_with_tensor, socle_of_module, socle_of_module_inclusion,
    subrepresentation, subrepresentation_inclusion, sum_of_submodules, sum_of_submodules_list,
    top_of_module, top_of_module_homomorphism, top_of_module_projection, trace_of_module,
    underlying_linear_map, zero_module, zero_module_mapping,
};
pub use poset::{Poset, PosetError};
pub use projective_presentation::{ProjectivePresentation, ProjectivePresentationError};
pub use quiver::{
    Arrow, ArrowIndex, ArrowSpec, DYNKIN_NOT_CONNECTED_MESSAGE, DYNKIN_UNORIENTED_CYCLE_MESSAGE,
    DynkinDecision, DynkinFamily, DynkinOrientation, DynkinType, Generator, INFO_QUIVER,
    NamedArrowSpec, NextPathPosition, OrderedQuiver, Path, PathEndpoint, PathError, PathIterator,
    Quiver, QuiverCardinality, QuiverError, Vertex, VertexIndex, unique_quiver_name,
};
pub use quotient::{
    AlgebraDimension, AlgebraGlobalDimension, NaturalQuotientMap, PathAlgebraIdeal,
    QuotientCanonicalBasis, QuotientError, QuotientNiceBasis, QuotientPathAlgebra,
    QuotientPathAlgebraElement, is_homogeneous_list_of_elements,
};
pub use static_dictionary::{
    DictionaryMatch, QuiverStaticDictionary, STATIC_DICTIONARY_CHECKING_CYCLE_INFO,
    STATIC_DICTIONARY_CREATING_AUTOMATA_INFO, StaticDictionary, StaticDictionaryError,
};
pub use structure_constants::{
    OppositeStructureConstantsAlgebra, OppositeStructureConstantsAlgebraElement,
    StructureConstantsAlgebra, StructureConstantsAlgebraElement, StructureConstantsError,
    StructureConstantsQuotient, StructureConstantsSubalgebra,
};
pub use suffix_tree::{
    INFO_DYNAMIC_DICTIONARY, SuffixTree, SuffixTreeEdgeId, SuffixTreeError, SuffixTreeNodeId,
    SuffixTreeNodeKind, SuffixTreeOccurrence, SuffixTreeToken,
};
pub use tensor::{
    PathAlgebraElementTermRecord, PathAlgebraTensorError, QuiverProduct,
    QuiverProductDecomposition, QuiverProductError, TensorAlgebraInclusionError,
    TensorProductAlgebra, TensorProductAlgebraError, include_in_path_algebra,
    include_in_product_quiver, path_algebra_element_terms, project_from_product_quiver,
    quiver_product, reverse_path, simple_tensor, vertices_of_path_algebra, walk_of_path_or_vertex,
};
pub use triangular_reduced::TriangularReducedError;
pub use util::{null_list, positions_nonzero, positions_nonzero_in_field};
pub use vector_reduction::{VectorCoefficients, VectorReductionError, VectorTip};
