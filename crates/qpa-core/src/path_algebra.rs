//! Exact sparse path algebras over Numerica fields.
//!
//! QPA's `PathRing`/`PathAlgebra` methods in `lib/pathalg.gi` delegate sparse
//! arithmetic to GAP's magma-ring implementation in `lib/mgmring.gi`. This
//! module translates that representation and its `ZippedSum`/`ZippedProduct`
//! semantics: terms are strictly path-sorted, equal paths are combined, zero
//! coefficients disappear, and QPA's zero path is quotiented out.

use std::{
    cmp::Ordering,
    fmt,
    sync::{Arc, OnceLock, RwLock, Weak},
};

use numerica::domains::{
    Field, InternalOrdering,
    rational::{Rational, RationalField},
};
use thiserror::Error;

use crate::{
    CachedAttribute, CachedNakayamaAutomorphismOrder, OrderedQuiver, OrderingError, Path,
    PathEndpoint, PathError, Quiver, QuiverError, QuiverOrdering, Vertex,
    combinatorial_map::CombinatorialMap,
    path_module::{PathAlgebraModule, WeakPathAlgebraModule, WeakPathAlgebraModuleListCache},
    quotient::{QuotientPathAlgebra, WeakQuotientPathAlgebra},
    tensor::WeakTensorProductAlgebra,
};
use crate::{DenseMatrix, TensorProductAlgebra};

struct PathAlgebraData<F: Field> {
    field: F,
    quiver: Quiver,
    ordering: QuiverOrdering,
    opposite_weak: RwLock<Option<Weak<PathAlgebraData<F>>>>,
    opposite_strong: OnceLock<PathAlgebra<F>>,
    enveloping_weak: RwLock<Option<WeakTensorProductAlgebra<F>>>,
    regular_enveloping_module_weak: RwLock<Option<WeakPathAlgebraModule<F>>>,
    regular_enveloping_module_bound: OnceLock<()>,
    regular_enveloping_module_fail: OnceLock<()>,
    dual_enveloping_module_weak: RwLock<Option<WeakPathAlgebraModule<F>>>,
    dual_enveloping_module_bound: OnceLock<()>,
    dual_enveloping_module_fail: OnceLock<()>,
    zero_module_weak: RwLock<Option<WeakPathAlgebraModule<F>>>,
    indec_projective_modules: WeakPathAlgebraModuleListCache<F>,
    indec_injective_modules: WeakPathAlgebraModuleListCache<F>,
    simple_modules: WeakPathAlgebraModuleListCache<F>,
    gentle_maximal_paths: OnceLock<Vec<Path>>,
    gentle_combinatorial_map: OnceLock<CombinatorialMap>,
    gorenstein_dimension: OnceLock<usize>,
    cartan_matrix: OnceLock<Option<Vec<Vec<usize>>>>,
    coxeter_matrix: OnceLock<Option<DenseMatrix<RationalField>>>,
    coxeter_polynomial: OnceLock<Option<Vec<Rational>>>,
    nakayama_permutation: OnceLock<Option<Vec<usize>>>,
    nakayama_automorphism: OnceLock<Option<DenseMatrix<F>>>,
    order_of_nakayama_automorphism: OnceLock<Option<CachedNakayamaAutomorphismOrder>>,
    elements_of_units: OnceLock<Vec<Vec<F::Element>>>,
    trivial_extension_level: OnceLock<usize>,
    trivial_extension_weak: RwLock<Option<WeakQuotientPathAlgebra<F>>>,
    trivial_extension_bound: OnceLock<()>,
    trivial_extension_fail: OnceLock<()>,
}

/// QPA's external representation of a path-algebra element.
///
/// GAP stores this as a two-entry list containing the coefficient zero and a
/// flattened alternating list of path external representations and
/// coefficients. Rust gives the alternating entries a typed pair shape.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PathAlgebraExternalRepresentation<C> {
    zero_coefficient: C,
    terms: Vec<(Vec<usize>, C)>,
}

impl<C> PathAlgebraExternalRepresentation<C> {
    /// Construct the typed form of GAP's two-entry external representation.
    #[must_use]
    pub const fn new(zero_coefficient: C, terms: Vec<(Vec<usize>, C)>) -> Self {
        Self {
            zero_coefficient,
            terms,
        }
    }

    /// Return the stored zero coefficient.
    #[must_use]
    pub const fn zero_coefficient(&self) -> &C {
        &self.zero_coefficient
    }

    /// Return path external representations paired with their coefficients.
    #[must_use]
    pub fn terms(&self) -> &[(Vec<usize>, C)] {
        &self.terms
    }

    /// Consume the descriptor into its two GAP representation components.
    #[must_use]
    pub fn into_parts(self) -> (C, Vec<(Vec<usize>, C)>) {
        (self.zero_coefficient, self.terms)
    }
}

/// A QPA path algebra over an exact Numerica field.
///
/// Cloning this handle preserves algebra/family identity. Calling a constructor
/// again, even with equal inputs, creates a distinct algebra family as in GAP.
#[derive(Clone)]
pub struct PathAlgebra<F: Field>(Arc<PathAlgebraData<F>>);

impl<F: Field> PathAlgebra<F> {
    /// Translate QPA's `PathAlgebra(field, quiver)` with its default ordering.
    ///
    /// # Errors
    ///
    /// Returns an error only if the source-translated default ordering cannot
    /// be constructed for the quiver.
    pub fn new(field: F, quiver: Quiver) -> Result<Self, PathAlgebraError> {
        let ordering = quiver.ordering()?;
        Self::with_ordering(field, quiver, ordering)
    }

    /// Translate `PathAlgebra(field, OrderedBy(Q, O))`.
    ///
    /// QPA stores the old ordering closure on a position-preserving fresh
    /// quiver family. Rust rebinds that closure to the fresh family before it
    /// enters the path-algebra object.
    ///
    /// # Errors
    ///
    /// Returns a checked ordering error if a foreign ordering with an
    /// incompatible generator count was supplied to `OrderedBy`.
    pub fn from_ordered_quiver(field: F, ordered: OrderedQuiver) -> Result<Self, PathAlgebraError> {
        let (quiver, ordering) = ordered.into_parts()?;
        Self::with_ordering(field, quiver, ordering)
    }

    /// Construct QPA's path algebra with an explicit quiver ordering.
    ///
    /// # Errors
    ///
    /// Returns [`PathAlgebraError::OrderingQuiverMismatch`] when the ordering
    /// was constructed for another quiver family.
    pub fn with_ordering(
        field: F,
        quiver: Quiver,
        ordering: QuiverOrdering,
    ) -> Result<Self, PathAlgebraError> {
        if !ordering.quiver().is_identical(&quiver) {
            return Err(PathAlgebraError::OrderingQuiverMismatch);
        }
        Ok(Self(Arc::new(PathAlgebraData {
            field,
            quiver,
            ordering,
            opposite_weak: RwLock::new(None),
            opposite_strong: OnceLock::new(),
            enveloping_weak: RwLock::new(None),
            regular_enveloping_module_weak: RwLock::new(None),
            regular_enveloping_module_bound: OnceLock::new(),
            regular_enveloping_module_fail: OnceLock::new(),
            dual_enveloping_module_weak: RwLock::new(None),
            dual_enveloping_module_bound: OnceLock::new(),
            dual_enveloping_module_fail: OnceLock::new(),
            zero_module_weak: RwLock::new(None),
            indec_projective_modules: WeakPathAlgebraModuleListCache::new(),
            indec_injective_modules: WeakPathAlgebraModuleListCache::new(),
            simple_modules: WeakPathAlgebraModuleListCache::new(),
            gentle_maximal_paths: OnceLock::new(),
            gentle_combinatorial_map: OnceLock::new(),
            gorenstein_dimension: OnceLock::new(),
            cartan_matrix: OnceLock::new(),
            coxeter_matrix: OnceLock::new(),
            coxeter_polynomial: OnceLock::new(),
            nakayama_permutation: OnceLock::new(),
            nakayama_automorphism: OnceLock::new(),
            order_of_nakayama_automorphism: OnceLock::new(),
            elements_of_units: OnceLock::new(),
            trivial_extension_level: OnceLock::new(),
            trivial_extension_weak: RwLock::new(None),
            trivial_extension_bound: OnceLock::new(),
            trivial_extension_fail: OnceLock::new(),
        })))
    }

    /// Return the exact coefficient field.
    #[must_use]
    pub fn field(&self) -> &F {
        &self.0.field
    }

    /// Return the underlying quiver.
    #[must_use]
    pub fn quiver(&self) -> &Quiver {
        &self.0.quiver
    }

    pub(crate) fn cached_elements_of_units_coordinates(&self) -> Option<&[Vec<F::Element>]> {
        self.0.elements_of_units.get().map(Vec::as_slice)
    }

    pub(crate) fn cache_elements_of_units_coordinates(
        &self,
        coordinates: Vec<Vec<F::Element>>,
    ) -> &[Vec<F::Element>] {
        self.0
            .elements_of_units
            .get_or_init(|| coordinates)
            .as_slice()
    }

    /// Return the QPA `IsPathRing` category value.
    ///
    /// Rust's initial coefficient contract is field-only, so [`PathAlgebra`]
    /// represents the intersection of QPA path rings and algebras. General
    /// coefficient rings such as `PathRing(Integers, Q)` are outside this
    /// coefficient tier.
    #[must_use]
    pub const fn is_path_ring(&self) -> bool {
        true
    }

    /// Return the QPA `IsPathAlgebra` synonym value.
    #[must_use]
    pub const fn is_path_algebra(&self) -> bool {
        true
    }

    /// Return the QPA `IsQuiverAlgebra` true-method value.
    #[must_use]
    pub const fn is_quiver_algebra(&self) -> bool {
        true
    }

    /// Translate QPA's concise path-algebra `ViewObj` text.
    ///
    /// The caller supplies GAP-facing field text because Numerica owns exact
    /// arithmetic values, not GAP presentation syntax.
    #[must_use]
    pub fn view_obj_with(&self, field_display: &str) -> String {
        format!("<{field_display}[{}]>", self.quiver().view_string())
    }

    /// Translate QPA's verbose path-algebra `PrintObj` text.
    ///
    /// QPA prints a bound quiver object name when one exists and otherwise its
    /// concise view. Rust quivers do not carry GAP's mutable object-name
    /// attribute, so callers may supply that optional presentation name.
    #[must_use]
    pub fn print_obj_with(&self, field_display: &str, quiver_object_name: Option<&str>) -> String {
        let quiver = quiver_object_name.map_or_else(|| self.quiver().view_string(), str::to_owned);
        format!("<Path algebra of the quiver {quiver} over the field {field_display}>")
    }

    /// Return QPA's cached `MaximalPathsOfGentleAlgebra` attribute, if bound.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_gentle_maximal_paths(&self) -> Option<&[Path]> {
        self.0.gentle_maximal_paths.get().map(Vec::as_slice)
    }

    /// Bind QPA's `MaximalPathsOfGentleAlgebra` attribute if it is unbound.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_gentle_maximal_paths(&self, paths: Vec<Path>) -> &[Path] {
        self.0.gentle_maximal_paths.get_or_init(|| paths)
    }

    /// Return QPA's cached `CombinatorialMapOfGentleAlgebra`, if bound.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_gentle_combinatorial_map(&self) -> Option<&CombinatorialMap> {
        self.0.gentle_combinatorial_map.get()
    }

    /// Bind QPA's `CombinatorialMapOfGentleAlgebra` attribute if it is unbound.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_gentle_combinatorial_map(&self, map: CombinatorialMap) -> &CombinatorialMap {
        self.0.gentle_combinatorial_map.get_or_init(|| map)
    }

    /// Return the path ordering used to canonicalize terms.
    #[must_use]
    pub fn ordering(&self) -> &QuiverOrdering {
        &self.0.ordering
    }

    /// Preserve the source failure of QPA's list-plus-path
    /// `IsLeftUniform` overload.
    ///
    /// The installed GAP method accepts `[IsList, IsPath]` but immediately
    /// calls `CoefficientsAndMagmaElements` on the path argument, for which no
    /// method exists. Consequently every reachable call fails before the list
    /// is inspected.
    ///
    /// # Errors
    ///
    /// Always returns [`PathAlgebraError::UniformityListPathDispatch`].
    pub fn is_left_uniform_list_at_path(
        &self,
        _elements: &[PathAlgebraElement<F>],
        _path: &Path,
    ) -> Result<bool, PathAlgebraError> {
        Err(PathAlgebraError::UniformityListPathDispatch)
    }

    /// Preserve the source failure of QPA's list-plus-path
    /// `IsRightUniform` overload.
    ///
    /// # Errors
    ///
    /// Always returns [`PathAlgebraError::UniformityListPathDispatch`], for
    /// the same source-level dispatch defect as
    /// [`Self::is_left_uniform_list_at_path`].
    pub fn is_right_uniform_list_at_path(
        &self,
        _elements: &[PathAlgebraElement<F>],
        _path: &Path,
    ) -> Result<bool, PathAlgebraError> {
        Err(PathAlgebraError::UniformityListPathDispatch)
    }

    /// Return whether two handles denote the identical QPA algebra object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    pub(crate) fn cached_enveloping_algebra(&self) -> Option<TensorProductAlgebra<F>> {
        self.0
            .enveloping_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakTensorProductAlgebra::upgrade)
    }

    pub(crate) fn cache_enveloping_algebra(
        &self,
        enveloping: &TensorProductAlgebra<F>,
    ) -> TensorProductAlgebra<F> {
        let mut slot = self
            .0
            .enveloping_weak
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(cached) = slot.as_ref().and_then(WeakTensorProductAlgebra::upgrade) {
            return cached;
        }
        *slot = Some(enveloping.downgrade());
        enveloping.clone()
    }

    pub(crate) fn cached_regular_enveloping_module(&self) -> Option<PathAlgebraModule<F>> {
        self.0
            .regular_enveloping_module_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModule::upgrade)
    }

    pub(crate) fn has_regular_enveloping_module_attribute(&self) -> bool {
        self.0.regular_enveloping_module_bound.get().is_some()
    }

    pub(crate) fn regular_enveloping_module_attribute_is_fail(&self) -> bool {
        self.0.regular_enveloping_module_fail.get().is_some()
    }

    pub(crate) fn mark_regular_enveloping_module_attribute_fail(&self) {
        let _ = self.0.regular_enveloping_module_fail.set(());
        let _ = self.0.regular_enveloping_module_bound.set(());
    }

    pub(crate) fn cache_regular_enveloping_module(
        &self,
        module: &PathAlgebraModule<F>,
    ) -> PathAlgebraModule<F> {
        let mut slot = self
            .0
            .regular_enveloping_module_weak
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(cached) = slot.as_ref().and_then(WeakPathAlgebraModule::upgrade) {
            return cached;
        }
        *slot = Some(module.downgrade());
        let _ = self.0.regular_enveloping_module_bound.set(());
        module.clone()
    }

    pub(crate) fn cached_dual_enveloping_module(&self) -> Option<PathAlgebraModule<F>> {
        self.0
            .dual_enveloping_module_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModule::upgrade)
    }

    pub(crate) fn cached_indec_projective_module(
        &self,
        index: usize,
    ) -> Option<PathAlgebraModule<F>> {
        self.0.indec_projective_modules.cached_at(index)
    }

    pub(crate) fn merge_indec_projective_module(
        &self,
        index: usize,
        generated: &PathAlgebraModule<F>,
    ) -> Option<PathAlgebraModule<F>> {
        self.0.indec_projective_modules.merge_at(
            index,
            self.quiver().number_of_vertices(),
            generated,
        )
    }

    pub(crate) fn cached_indec_projective_modules(&self) -> Option<Vec<PathAlgebraModule<F>>> {
        self.0.indec_projective_modules.cached()
    }

    pub(crate) fn cached_zero_module(&self) -> Option<PathAlgebraModule<F>> {
        self.0
            .zero_module_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakPathAlgebraModule::upgrade)
    }

    pub(crate) fn cache_zero_module(&self, module: &PathAlgebraModule<F>) -> PathAlgebraModule<F> {
        let mut slot = self
            .0
            .zero_module_weak
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(cached) = slot.as_ref().and_then(WeakPathAlgebraModule::upgrade) {
            return cached;
        }
        *slot = Some(module.downgrade());
        module.clone()
    }

    pub(crate) fn merge_indec_projective_modules(
        &self,
        generated: &[PathAlgebraModule<F>],
    ) -> Vec<PathAlgebraModule<F>> {
        self.0.indec_projective_modules.merge(generated)
    }

    pub(crate) fn cached_indec_injective_modules(&self) -> Option<Vec<PathAlgebraModule<F>>> {
        self.0.indec_injective_modules.cached()
    }

    pub(crate) fn merge_indec_injective_modules(
        &self,
        generated: &[PathAlgebraModule<F>],
    ) -> Vec<PathAlgebraModule<F>> {
        self.0.indec_injective_modules.merge(generated)
    }

    pub(crate) fn cached_simple_modules(&self) -> Option<Vec<PathAlgebraModule<F>>> {
        self.0.simple_modules.cached()
    }

    pub(crate) fn merge_simple_modules(
        &self,
        generated: &[PathAlgebraModule<F>],
    ) -> Vec<PathAlgebraModule<F>> {
        self.0.simple_modules.merge(generated)
    }

    pub(crate) fn cached_cartan_matrix_attribute(&self) -> CachedAttribute<Vec<Vec<usize>>> {
        match self.0.cartan_matrix.get() {
            None => CachedAttribute::Unbound,
            Some(None) => CachedAttribute::Fail,
            Some(Some(value)) => CachedAttribute::Value(value.clone()),
        }
    }

    pub(crate) fn cache_cartan_matrix_attribute(
        &self,
        value: Option<Vec<Vec<usize>>>,
    ) -> Option<Vec<Vec<usize>>> {
        self.0.cartan_matrix.get_or_init(|| value).clone()
    }

    pub(crate) fn cached_coxeter_matrix_attribute(
        &self,
    ) -> CachedAttribute<DenseMatrix<RationalField>> {
        match self.0.coxeter_matrix.get() {
            None => CachedAttribute::Unbound,
            Some(None) => CachedAttribute::Fail,
            Some(Some(value)) => CachedAttribute::Value(value.clone()),
        }
    }

    pub(crate) fn cache_coxeter_matrix_attribute(
        &self,
        value: Option<DenseMatrix<RationalField>>,
    ) -> Option<DenseMatrix<RationalField>> {
        self.0.coxeter_matrix.get_or_init(|| value).clone()
    }

    pub(crate) fn cached_coxeter_polynomial_attribute(&self) -> CachedAttribute<Vec<Rational>> {
        match self.0.coxeter_polynomial.get() {
            None => CachedAttribute::Unbound,
            Some(None) => CachedAttribute::Fail,
            Some(Some(value)) => CachedAttribute::Value(value.clone()),
        }
    }

    pub(crate) fn cache_coxeter_polynomial_attribute(
        &self,
        value: Option<Vec<Rational>>,
    ) -> Option<Vec<Rational>> {
        self.0.coxeter_polynomial.get_or_init(|| value).clone()
    }

    pub(crate) fn cached_nakayama_permutation_attribute(&self) -> CachedAttribute<Vec<usize>> {
        match self.0.nakayama_permutation.get() {
            None => CachedAttribute::Unbound,
            Some(None) => CachedAttribute::False,
            Some(Some(value)) => CachedAttribute::Value(value.clone()),
        }
    }

    pub(crate) fn cache_nakayama_permutation_attribute(
        &self,
        value: Option<Vec<usize>>,
    ) -> Option<Vec<usize>> {
        self.0.nakayama_permutation.get_or_init(|| value).clone()
    }

    pub(crate) fn cached_nakayama_automorphism_attribute(&self) -> CachedAttribute<DenseMatrix<F>> {
        match self.0.nakayama_automorphism.get() {
            None => CachedAttribute::Unbound,
            Some(None) => CachedAttribute::False,
            Some(Some(value)) => CachedAttribute::Value(value.clone()),
        }
    }

    pub(crate) fn cache_nakayama_automorphism_attribute(
        &self,
        value: Option<DenseMatrix<F>>,
    ) -> Option<DenseMatrix<F>> {
        self.0.nakayama_automorphism.get_or_init(|| value).clone()
    }

    pub(crate) fn cached_order_of_nakayama_automorphism_attribute(
        &self,
    ) -> CachedAttribute<CachedNakayamaAutomorphismOrder> {
        match self.0.order_of_nakayama_automorphism.get() {
            None => CachedAttribute::Unbound,
            Some(None) => CachedAttribute::False,
            Some(Some(value)) => CachedAttribute::Value(value.clone()),
        }
    }

    pub(crate) fn cache_order_of_nakayama_automorphism_attribute(
        &self,
        value: Option<CachedNakayamaAutomorphismOrder>,
    ) -> Option<CachedNakayamaAutomorphismOrder> {
        self.0
            .order_of_nakayama_automorphism
            .get_or_init(|| value)
            .clone()
    }

    pub(crate) fn has_dual_enveloping_module_attribute(&self) -> bool {
        self.0.dual_enveloping_module_bound.get().is_some()
    }

    pub(crate) fn dual_enveloping_module_attribute_is_fail(&self) -> bool {
        self.0.dual_enveloping_module_fail.get().is_some()
    }

    pub(crate) fn mark_dual_enveloping_module_attribute_fail(&self) {
        let _ = self.0.dual_enveloping_module_fail.set(());
        let _ = self.0.dual_enveloping_module_bound.set(());
    }

    pub(crate) fn cache_dual_enveloping_module(
        &self,
        module: &PathAlgebraModule<F>,
    ) -> PathAlgebraModule<F> {
        let mut slot = self
            .0
            .dual_enveloping_module_weak
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(cached) = slot.as_ref().and_then(WeakPathAlgebraModule::upgrade) {
            return cached;
        }
        *slot = Some(module.downgrade());
        let _ = self.0.dual_enveloping_module_bound.set(());
        module.clone()
    }

    pub(crate) fn trivial_extension_level(&self) -> usize {
        self.0.trivial_extension_level.get().copied().unwrap_or(0)
    }

    pub(crate) fn set_trivial_extension_level(&self, level: usize) {
        let _ = self.0.trivial_extension_level.set(level);
    }

    pub(crate) fn cached_trivial_extension(&self) -> Option<QuotientPathAlgebra<F>> {
        self.0
            .trivial_extension_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(WeakQuotientPathAlgebra::upgrade)
    }

    pub(crate) fn trivial_extension_attribute_is_fail(&self) -> bool {
        self.0.trivial_extension_fail.get().is_some()
    }

    pub(crate) fn mark_trivial_extension_attribute_fail(&self) {
        let _ = self.0.trivial_extension_fail.set(());
        let _ = self.0.trivial_extension_bound.set(());
    }

    pub(crate) fn cache_trivial_extension(
        &self,
        extension: &QuotientPathAlgebra<F>,
    ) -> QuotientPathAlgebra<F> {
        let mut slot = self
            .0
            .trivial_extension_weak
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(cached) = slot.as_ref().and_then(WeakQuotientPathAlgebra::upgrade) {
            return cached;
        }
        *slot = Some(extension.downgrade());
        let _ = self.0.trivial_extension_bound.set(());
        extension.clone()
    }

    /// Return the additive zero.
    #[must_use]
    pub fn zero(&self) -> PathAlgebraElement<F> {
        PathAlgebraElement {
            algebra: self.clone(),
            terms: Vec::new(),
        }
    }

    /// Return QPA's identity: the sum of all vertex paths.
    #[must_use]
    pub fn one(&self) -> PathAlgebraElement<F> {
        let coefficient = self.field().one();
        let terms = self
            .quiver()
            .vertices()
            .map(|vertex| PathAlgebraTerm {
                path: vertex.path(),
                coefficient: coefficient.clone(),
            })
            .collect();
        PathAlgebraElement {
            algebra: self.clone(),
            terms,
        }
    }

    /// Embed a QPA path with coefficient one.
    ///
    /// The distinguished zero path maps to the algebra zero, implementing
    /// `MagmaRingModuloSpanOfZero`.
    ///
    /// # Errors
    ///
    /// Returns an error for a path from another quiver.
    pub fn element(&self, path: Path) -> Result<PathAlgebraElement<F>, PathAlgebraError> {
        self.from_terms([(path, self.field().one())])
    }

    /// Construct and normalize a sparse element from path/coefficient pairs.
    ///
    /// This translates GAP's `ElementOfMagmaRing`: sort paths, merge repeated
    /// paths, remove zero coefficients, and remove the quiver zero path.
    ///
    /// # Errors
    ///
    /// Returns an error if any path belongs to another quiver.
    pub fn from_terms(
        &self,
        terms: impl IntoIterator<Item = (Path, F::Element)>,
    ) -> Result<PathAlgebraElement<F>, PathAlgebraError> {
        let terms = normalize_terms(self, terms)?;
        Ok(PathAlgebraElement {
            algebra: self.clone(),
            terms,
        })
    }

    /// Translate GAP's path-algebra `ObjByExtRep` method.
    ///
    /// QPA documents that the supplied descriptor must already be normalized.
    /// Rust verifies that precondition so that its canonical sparse element
    /// invariant cannot be bypassed.
    ///
    /// # Errors
    ///
    /// Returns an error for a nonzero zero-coefficient slot, an invalid path
    /// descriptor, or terms that are not already in canonical sparse order.
    pub fn element_from_external_representation(
        &self,
        representation: &PathAlgebraExternalRepresentation<F::Element>,
    ) -> Result<PathAlgebraElement<F>, PathAlgebraError> {
        if !self.field().is_zero(representation.zero_coefficient()) {
            return Err(PathAlgebraError::InvalidExternalZeroCoefficient);
        }
        let terms = representation
            .terms()
            .iter()
            .map(|(path, coefficient)| {
                Ok((
                    Path::from_external_representation(self.quiver(), path)?
                        .ok_or(PathError::EmptyProduct)?,
                    coefficient.clone(),
                ))
            })
            .collect::<Result<Vec<_>, PathAlgebraError>>()?;
        let element = self.from_terms(terms)?;
        if &element.external_representation() != representation {
            return Err(PathAlgebraError::ExternalRepresentationNotNormalized);
        }
        Ok(element)
    }

    /// Return path-algebra generators in QPA quiver generator order: vertices
    /// followed by arrows.
    #[must_use]
    pub fn generators(&self) -> Vec<PathAlgebraElement<F>> {
        self.quiver()
            .vertices()
            .map(|vertex| vertex.path())
            .chain(self.quiver().arrows().map(|arrow| arrow.path()))
            .map(|path| self.monomial_element(path))
            .collect()
    }

    /// Return QPA's `RelationsOfAlgebra` value for a path algebra.
    #[must_use]
    pub const fn relations(&self) -> Vec<PathAlgebraElement<F>> {
        Vec::new()
    }

    /// Return QPA's path-algebra global-dimension value set by the constructor.
    #[must_use]
    pub fn qpa_global_dimension(&self) -> usize {
        usize::from(self.quiver().number_of_arrows() > 0)
    }

    /// Return QPA's installed `GorensteinDimension` attribute, if bound.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_gorenstein_dimension(&self) -> Option<usize> {
        self.0.gorenstein_dimension.get().copied()
    }

    /// Install QPA's finite `GorensteinDimension` attribute.
    #[doc(hidden)]
    #[must_use]
    pub fn cache_gorenstein_dimension(&self, value: usize) -> usize {
        *self.0.gorenstein_dimension.get_or_init(|| value)
    }

    /// Return whether QPA can construct the canonical basis by enumerating the
    /// finite underlying quiver.
    #[must_use]
    pub fn has_finite_canonical_basis(&self) -> bool {
        self.quiver().is_acyclic()
    }

    /// Translate QPA's canonical basis for an acyclic path algebra.
    ///
    /// `None` corresponds to QPA's `fail` for a cyclic quiver.
    #[must_use]
    pub fn canonical_basis(&self) -> Option<Vec<PathAlgebraElement<F>>> {
        self.has_finite_canonical_basis().then(|| {
            self.quiver()
                .paths()
                .map(|path| self.monomial_element(path))
                .collect()
        })
    }

    /// Translate QPA's cached `OppositePathAlgebra` construction.
    ///
    /// The opposite uses QPA's default ordering on the cached opposite quiver.
    /// Reciprocal weak/strong caches preserve observable double-opposite
    /// identity without creating a permanent `Arc` cycle.
    ///
    /// # Errors
    ///
    /// Propagates opposite-quiver and default-ordering construction failures.
    pub fn opposite(&self) -> Result<Self, PathAlgebraError> {
        if let Some(opposite) = self.0.opposite_strong.get() {
            return Ok(opposite.clone());
        }
        if let Some(opposite) = self
            .0
            .opposite_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(Weak::upgrade)
        {
            return Ok(Self(opposite));
        }

        let opposite = Self::new(self.field().clone(), self.quiver().opposite()?)?;
        if opposite.0.opposite_strong.set(self.clone()).is_err() {
            return Err(PathAlgebraError::OppositeCacheRace);
        }
        let mut slot = self
            .0
            .opposite_weak
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(cached) = slot.as_ref().and_then(Weak::upgrade) {
            return Ok(Self(cached));
        }
        *slot = Some(Arc::downgrade(&opposite.0));
        Ok(opposite)
    }

    /// Translate QPA's `MakeUniformOnRight` for path-algebra elements.
    ///
    /// Target vertices and final elements are deduplicated in first-occurrence
    /// order, matching GAP's `Unique` rather than a sorted Rust set.
    ///
    /// # Errors
    ///
    /// All nonempty inputs must belong to this identical path algebra.
    pub fn make_uniform_on_right(
        &self,
        generators: &[PathAlgebraElement<F>],
    ) -> Result<Vec<PathAlgebraElement<F>>, PathAlgebraError> {
        ensure_algebra_elements(self, generators)?;
        let mut uniform = Vec::new();
        for generator in generators {
            let mut targets = Vec::new();
            for term in generator.terms() {
                if let PathEndpoint::Vertex(target) = term.path().target()
                    && !targets.contains(&target)
                {
                    targets.push(target);
                }
            }
            for target in targets {
                let vertex = self.element(target.path())?;
                let component = generator.multiply(&vertex)?;
                if !component.is_zero() && !uniform.contains(&component) {
                    uniform.push(component);
                }
            }
        }
        Ok(uniform)
    }

    /// Translate QPA's `GeneratorsTimesArrowsOnRight` for path-algebra
    /// elements.
    ///
    /// # Errors
    ///
    /// Propagates right-uniformization, quiver lookup, and multiplication
    /// failures.
    pub fn generators_times_arrows_on_right(
        &self,
        generators: &[PathAlgebraElement<F>],
    ) -> Result<Vec<PathAlgebraElement<F>>, PathAlgebraError> {
        let uniform = self.make_uniform_on_right(generators)?;
        let mut extended = Vec::new();
        for generator in uniform {
            let PathEndpoint::Vertex(target) = generator.leading_monomial().target() else {
                continue;
            };
            for arrow in self.quiver().outgoing_arrows(&target)? {
                let product = generator.multiply(&self.element(arrow.path())?)?;
                if !product.is_zero() {
                    extended.push(product);
                }
            }
        }
        Ok(extended)
    }

    /// Translate QPA's `NthPowerOfArrowIdeal` generator list.
    ///
    /// The pinned source returns the arrow list for every exponent at most one
    /// because GAP's range `[1..n-1]` is empty there; Rust preserves that
    /// observable boundary behavior.
    ///
    /// # Errors
    ///
    /// Propagates source-translated right-extension failures.
    pub fn nth_power_of_arrow_ideal(
        &self,
        exponent: i64,
    ) -> Result<Vec<PathAlgebraElement<F>>, PathAlgebraError> {
        let mut generators = self
            .quiver()
            .arrows()
            .map(|arrow| self.element(arrow.path()))
            .collect::<Result<Vec<_>, _>>()?;
        for _ in 1..exponent {
            generators = self.generators_times_arrows_on_right(&generators)?;
        }
        Ok(generators)
    }

    /// Translate QPA's mutating `AddNthPowerToRelations`.
    ///
    /// # Errors
    ///
    /// Existing relations must belong to this algebra; arrow-power generation
    /// failures are propagated.
    pub fn add_nth_power_to_relations(
        &self,
        relations: &mut Vec<PathAlgebraElement<F>>,
        exponent: i64,
    ) -> Result<(), PathAlgebraError> {
        ensure_algebra_elements(self, relations)?;
        relations.extend(self.nth_power_of_arrow_ideal(exponent)?);
        Ok(())
    }

    fn monomial_element(&self, path: Path) -> PathAlgebraElement<F> {
        debug_assert!(!path.is_zero());
        debug_assert!(path.quiver().is_identical(self.quiver()));
        PathAlgebraElement {
            algebra: self.clone(),
            terms: vec![PathAlgebraTerm {
                path,
                coefficient: self.field().one(),
            }],
        }
    }
}

impl<F: Field> fmt::Debug for PathAlgebra<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PathAlgebra")
            .field("quiver", self.quiver())
            .field("ordering", self.ordering())
            .finish_non_exhaustive()
    }
}

impl<F: Field> PartialEq for PathAlgebra<F> {
    fn eq(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}

impl<F: Field> Eq for PathAlgebra<F> {}

/// One canonical nonzero path/coefficient pair.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PathAlgebraTerm<F: Field> {
    path: Path,
    coefficient: F::Element,
}

impl<F: Field> PathAlgebraTerm<F> {
    /// Return the monomial path.
    #[must_use]
    pub const fn path(&self) -> &Path {
        &self.path
    }

    /// Return the nonzero coefficient.
    #[must_use]
    pub const fn coefficient(&self) -> &F::Element {
        &self.coefficient
    }
}

/// A canonical sparse element of a QPA path algebra.
#[derive(Clone, Debug)]
pub struct PathAlgebraElement<F: Field> {
    algebra: PathAlgebra<F>,
    terms: Vec<PathAlgebraTerm<F>>,
}

impl<F: Field> PathAlgebraElement<F> {
    /// Return the containing path algebra.
    #[must_use]
    pub const fn algebra(&self) -> &PathAlgebra<F> {
        &self.algebra
    }

    /// Return canonical terms in increasing QPA path order.
    #[must_use]
    pub fn terms(&self) -> &[PathAlgebraTerm<F>] {
        &self.terms
    }

    /// Return whether this is the algebra zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Translate GAP's two-entry `ExtRepOfObj` representation.
    #[must_use]
    pub fn external_representation(&self) -> PathAlgebraExternalRepresentation<F::Element> {
        PathAlgebraExternalRepresentation::new(
            self.algebra().field().zero(),
            self.terms()
                .iter()
                .map(|term| {
                    (
                        term.path().external_representation(),
                        term.coefficient().clone(),
                    )
                })
                .collect(),
        )
    }

    /// Translate GAP's path-ring element printing.
    ///
    /// The coefficient formatter supplies GAP-facing spelling such as
    /// `Z(2)^0`; Numerica deliberately does not encode GAP syntax.
    pub fn print_obj_with(&self, coefficient: impl Fn(&F::Element) -> String + Copy) -> String {
        if self.is_zero() {
            return "0".into();
        }
        self.terms()
            .iter()
            .map(|term| format!("({})*{}", coefficient(term.coefficient()), term.path()))
            .collect::<Vec<_>>()
            .join("+")
    }

    /// Translate QPA's `VertexPosition` on a trivial-path element.
    ///
    /// # Errors
    ///
    /// QPA documents a trivial-path precondition but otherwise indexes its
    /// representation unchecked. Rust reports violations explicitly.
    pub fn vertex_position(&self) -> Result<usize, PathAlgebraError> {
        let [term] = self.terms() else {
            return Err(PathAlgebraError::NotTrivialPathElement);
        };
        match (term.path().length(), term.path().source()) {
            (Some(0), PathEndpoint::Vertex(vertex)) => Ok(vertex.index().get()),
            _ => Err(PathAlgebraError::NotTrivialPathElement),
        }
    }

    /// Translate QPA's path-algebra element comparison from `pathalg.gi`.
    ///
    /// QPA compares terms from greatest monomial to least, comparing each
    /// monomial before its coefficient, and finally compares term counts.
    /// Numerica's [`InternalOrdering`] supplies the backend's deterministic
    /// coefficient ordering, including for fields without a mathematical
    /// order such as finite fields.
    ///
    /// # Errors
    ///
    /// Returns an error unless both values have identical algebra families.
    pub fn compare(&self, right: &Self) -> Result<Ordering, PathAlgebraError> {
        self.ensure_same_algebra(right)?;
        for (left_term, right_term) in self.terms.iter().rev().zip(right.terms.iter().rev()) {
            let path_ordering = self
                .algebra
                .ordering()
                .compare(&left_term.path, &right_term.path)?;
            if !path_ordering.is_eq() {
                return Ok(path_ordering);
            }
            let coefficient_ordering = left_term.coefficient.internal_cmp(&right_term.coefficient);
            if !coefficient_ordering.is_eq() {
                return Ok(coefficient_ordering);
            }
        }
        Ok(self.terms.len().cmp(&right.terms.len()))
    }

    /// Add elements using GAP's `ZippedSum` semantics.
    ///
    /// # Errors
    ///
    /// Returns an error unless both values have identical QPA algebra families.
    pub fn add(&self, right: &Self) -> Result<Self, PathAlgebraError> {
        self.ensure_same_algebra(right)?;
        self.algebra.from_terms(
            self.terms
                .iter()
                .chain(&right.terms)
                .map(|term| (term.path.clone(), term.coefficient.clone())),
        )
    }

    /// Return the additive inverse.
    ///
    /// # Errors
    ///
    /// This operation is infallible for a Numerica field; an explicit result
    /// preserves a uniform arithmetic API.
    pub fn neg(&self) -> Result<Self, PathAlgebraError> {
        self.algebra.from_terms(self.terms.iter().map(|term| {
            (
                term.path.clone(),
                self.algebra.field().neg(&term.coefficient),
            )
        }))
    }

    /// Subtract two elements.
    ///
    /// # Errors
    ///
    /// Returns an error unless both values have identical QPA algebra families.
    pub fn sub(&self, right: &Self) -> Result<Self, PathAlgebraError> {
        self.add(&right.neg()?)
    }

    /// Multiply elements using GAP's `ZippedProduct` and QPA's left-to-right
    /// path multiplication.
    ///
    /// # Errors
    ///
    /// Returns an error unless both values have identical QPA algebra families.
    pub fn multiply(&self, right: &Self) -> Result<Self, PathAlgebraError> {
        self.ensure_same_algebra(right)?;
        if self.is_zero() || right.is_zero() {
            return Ok(self.algebra.zero());
        }
        let mut products = Vec::with_capacity(self.terms.len().saturating_mul(right.terms.len()));
        for left_term in &self.terms {
            for right_term in &right.terms {
                products.push((
                    left_term.path.multiply(&right_term.path)?,
                    self.algebra
                        .field()
                        .mul(&left_term.coefficient, &right_term.coefficient),
                ));
            }
        }
        self.algebra.from_terms(products)
    }

    /// Multiply every coefficient on the right, translating GAP's
    /// `ElmTimesRingElm`.
    ///
    /// # Errors
    ///
    /// This operation is infallible for a Numerica field; an explicit result
    /// preserves a uniform arithmetic API.
    pub fn scale_right(&self, scalar: &F::Element) -> Result<Self, PathAlgebraError> {
        self.algebra.from_terms(self.terms.iter().map(|term| {
            (
                term.path.clone(),
                self.algebra.field().mul(&term.coefficient, scalar),
            )
        }))
    }

    /// Multiply every coefficient on the left, translating GAP's
    /// `RingElmTimesElm`.
    ///
    /// # Errors
    ///
    /// This operation is infallible for a Numerica field; an explicit result
    /// preserves a uniform arithmetic API.
    pub fn scale_left(&self, scalar: &F::Element) -> Result<Self, PathAlgebraError> {
        self.algebra.from_terms(self.terms.iter().map(|term| {
            (
                term.path.clone(),
                self.algebra.field().mul(scalar, &term.coefficient),
            )
        }))
    }

    /// Divide every coefficient by a nonzero field scalar.
    ///
    /// This is the scalar division used by QPA to make relations monic in
    /// `HighLevelGroebnerBasis`.
    ///
    /// # Errors
    ///
    /// Returns [`PathAlgebraError::DivisionByZero`] for the zero scalar.
    pub fn divide_scalar(&self, scalar: &F::Element) -> Result<Self, PathAlgebraError> {
        if self.algebra.field().is_zero(scalar) {
            return Err(PathAlgebraError::DivisionByZero);
        }
        self.algebra.from_terms(self.terms.iter().map(|term| {
            (
                term.path.clone(),
                self.algebra.field().div(&term.coefficient, scalar),
            )
        }))
    }

    /// Raise an element to a nonnegative power in the path algebra.
    ///
    /// Exponent zero returns QPA's sum-of-vertices identity.
    ///
    /// # Errors
    ///
    /// Propagates source-translated path-algebra arithmetic failures.
    pub fn pow(&self, mut exponent: u64) -> Result<Self, PathAlgebraError> {
        let mut result = self.algebra.one();
        let mut base = self.clone();
        while exponent > 0 {
            if exponent & 1 == 1 {
                result = result.multiply(&base)?;
            }
            exponent >>= 1;
            if exponent > 0 {
                base = base.multiply(&base)?;
            }
        }
        Ok(result)
    }

    /// Return QPA's `LeadingTerm` (tip), or zero for the zero element.
    #[must_use]
    pub fn leading_term(&self) -> Self {
        self.terms.last().map_or_else(
            || self.algebra.zero(),
            |term| Self {
                algebra: self.algebra.clone(),
                terms: vec![term.clone()],
            },
        )
    }

    /// Return QPA's `LeadingCoefficient`, including field zero for the zero
    /// algebra element.
    #[must_use]
    pub fn leading_coefficient(&self) -> F::Element {
        self.terms.last().map_or_else(
            || self.algebra.field().zero(),
            |term| term.coefficient.clone(),
        )
    }

    /// Return QPA's `LeadingMonomial`, including the quiver zero path for the
    /// zero algebra element.
    #[must_use]
    pub fn leading_monomial(&self) -> Path {
        self.terms
            .last()
            .map_or_else(|| self.algebra.quiver().zero(), |term| term.path.clone())
    }

    /// Translate QPA's `IsLeftUniform`.
    #[must_use]
    pub fn is_left_uniform(&self) -> bool {
        uniform_endpoint(&self.terms, Path::source)
    }

    /// Translate QPA's `IsRightUniform`.
    #[must_use]
    pub fn is_right_uniform(&self) -> bool {
        uniform_endpoint(&self.terms, Path::target)
    }

    /// Translate QPA's `IsUniform`.
    #[must_use]
    pub fn is_uniform(&self) -> bool {
        self.is_left_uniform() && self.is_right_uniform()
    }

    /// Translate QPA's `OppositePathAlgebraElement` for a parent-algebra
    /// element.
    ///
    /// # Errors
    ///
    /// Propagates opposite-algebra, opposite-path, and sparse construction
    /// failures.
    pub fn opposite(&self) -> Result<Self, PathAlgebraError> {
        let opposite = self.algebra().opposite()?;
        opposite.from_terms(
            self.terms()
                .iter()
                .map(|term| {
                    Ok((
                        self.algebra().quiver().opposite_path(term.path())?,
                        term.coefficient().clone(),
                    ))
                })
                .collect::<Result<Vec<_>, PathAlgebraError>>()?,
        )
    }

    fn ensure_same_algebra(&self, other: &Self) -> Result<(), PathAlgebraError> {
        if self.algebra.is_identical(&other.algebra) {
            Ok(())
        } else {
            Err(PathAlgebraError::DifferentAlgebras)
        }
    }
}

impl<F: Field> PartialEq for PathAlgebraElement<F> {
    fn eq(&self, other: &Self) -> bool {
        self.algebra == other.algebra && self.terms == other.terms
    }
}

impl<F: Field> Eq for PathAlgebraElement<F> {}

fn normalize_terms<F: Field>(
    algebra: &PathAlgebra<F>,
    terms: impl IntoIterator<Item = (Path, F::Element)>,
) -> Result<Vec<PathAlgebraTerm<F>>, PathAlgebraError> {
    let mut result: Vec<PathAlgebraTerm<F>> = Vec::new();
    for (path, coefficient) in terms {
        if !path.quiver().is_identical(algebra.quiver()) {
            return Err(PathAlgebraError::PathQuiverMismatch);
        }
        if path.is_zero() || algebra.field().is_zero(&coefficient) {
            continue;
        }
        if let Some(existing) = result.iter_mut().find(|term| term.path == path) {
            algebra
                .field()
                .add_assign(&mut existing.coefficient, &coefficient);
        } else {
            result.push(PathAlgebraTerm { path, coefficient });
        }
    }
    result.retain(|term| !algebra.field().is_zero(&term.coefficient));
    result.sort_by(|left, right| {
        algebra
            .ordering()
            .compare(&left.path, &right.path)
            .expect("normalized terms belong to the ordering's quiver")
    });
    Ok(result)
}

fn uniform_endpoint<F: Field>(
    terms: &[PathAlgebraTerm<F>],
    endpoint: impl Fn(&Path) -> PathEndpoint,
) -> bool {
    let mut endpoints = terms.iter().filter_map(|term| match endpoint(&term.path) {
        PathEndpoint::Vertex(vertex) => Some(vertex),
        PathEndpoint::Zero(_) => None,
    });
    let Some(first) = endpoints.next() else {
        return true;
    };
    endpoints.all(|candidate: Vertex| candidate == first)
}

fn ensure_algebra_elements<F: Field>(
    algebra: &PathAlgebra<F>,
    elements: &[PathAlgebraElement<F>],
) -> Result<(), PathAlgebraError> {
    if elements
        .iter()
        .all(|element| algebra.is_identical(element.algebra()))
    {
        Ok(())
    } else {
        Err(PathAlgebraError::DifferentAlgebras)
    }
}

/// Translate QPA's `OppositeRelations` on parent-algebra relations.
///
/// # Errors
///
/// Propagates each relation's opposite-element construction failure.
pub fn opposite_relations<F: Field>(
    relations: &[PathAlgebraElement<F>],
) -> Result<Vec<PathAlgebraElement<F>>, PathAlgebraError> {
    relations.iter().map(PathAlgebraElement::opposite).collect()
}

/// Failure to construct or operate on a QPA path algebra.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PathAlgebraError {
    /// An explicit ordering belongs to another quiver family.
    #[error("the ordering and path algebra use different quivers")]
    OrderingQuiverMismatch,
    /// A term path belongs to another quiver family.
    #[error("a term path belongs to a different quiver")]
    PathQuiverMismatch,
    /// GAP installs arithmetic only for identical element families.
    #[error("the elements belong to different path algebras")]
    DifferentAlgebras,
    /// Field-scalar division is undefined for zero.
    #[error("cannot divide a path-algebra element by zero")]
    DivisionByZero,
    /// The zero slot in a GAP external representation was not field zero.
    #[error("the external representation has an invalid zero coefficient")]
    InvalidExternalZeroCoefficient,
    /// GAP's documented normalized-descriptor precondition was violated.
    #[error("the path-algebra external representation is not normalized")]
    ExternalRepresentationNotNormalized,
    /// QPA's list-plus-path uniformity overload calls a missing GAP method on
    /// its path argument before inspecting the list.
    #[error("QPA's list uniformity overload cannot extract coefficients from a path")]
    UniformityListPathDispatch,
    /// `VertexPosition` requires one trivial-path monomial.
    #[error("the element is not a scalar multiple of one trivial path")]
    NotTrivialPathElement,
    /// QPA's `PathRemoval` implementation indexes the last arrow of each term.
    #[error("QPA PathRemoval is undefined on an element containing a vertex term")]
    VertexTermInPathRemoval,
    /// A canonical sparse term violated its stored path/quiver invariant.
    #[error("PathRemoval encountered an invalid stored path term")]
    InvalidPathTermInPathRemoval,
    /// A surviving walk could not be mapped by arrow name to the rebuilt quiver.
    #[error("the arrow {name} is absent from the rebuilt vertex-quotient quiver")]
    MissingArrowInVertexQuotient {
        /// The original QPA arrow name.
        name: Box<str>,
    },
    /// Concurrent opposite-algebra construction lost its reciprocal cache.
    #[error("concurrent opposite-path-algebra cache initialization failed")]
    OppositeCacheRace,
    /// A source-translated ordering operation failed.
    #[error(transparent)]
    Ordering(#[from] OrderingError),
    /// QPA path multiplication failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// A source-translated quiver lookup or construction failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::Field;
    use numerica::domains::{
        Ring,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, Rational, RationalField},
    };
    use serde_json::{Value, json};

    use crate::{ArrowIndex, ArrowSpec, Quiver, VertexIndex};

    use super::{PathAlgebra, PathAlgebraElement};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn arrow(index: usize) -> ArrowIndex {
        ArrowIndex::new(index).unwrap()
    }

    fn quiver() -> Quiver {
        Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "a"),
                ArrowSpec::named(vertex(1), vertex(2), "b"),
            ],
        )
        .unwrap()
    }

    fn q(numerator: i64, denominator: i64) -> Rational {
        Q.to_element(numerator.into(), denominator.into(), true)
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the path-algebra fixture uses only small integers")
    }

    fn element_json<F: Field>(
        element: &PathAlgebraElement<F>,
        coefficient_json: impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().to_string(),
                        "path_external_representation": term.path().external_representation(),
                        "coefficient": coefficient_json(term.coefficient())
                    })
                })
                .collect(),
        )
    }

    fn rational_element_json(element: &PathAlgebraElement<RationalField>) -> Value {
        element_json(element, |coefficient| {
            json!([
                small_integer(coefficient.numerator_ref()),
                small_integer(coefficient.denominator_ref())
            ])
        })
    }

    fn zp_element_json(element: &PathAlgebraElement<Zp>, field: &Zp) -> Value {
        element_json(element, |coefficient| {
            small_integer(&field.to_integer(coefficient)).into()
        })
    }

    fn z2_element_json(element: &PathAlgebraElement<Z2>) -> Value {
        element_json(element, |coefficient| {
            small_integer(&Z2.to_integer(coefficient)).into()
        })
    }

    #[test]
    fn one_is_the_sum_of_vertices_and_acts_locally() {
        let algebra = PathAlgebra::new(Q, quiver()).unwrap();
        let a = algebra
            .element(algebra.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        assert_eq!(algebra.one().multiply(&a).unwrap(), a);
        assert_eq!(a.multiply(&algebra.one()).unwrap(), a);
        assert_eq!(algebra.one().terms().len(), 2);
    }

    #[test]
    fn normalization_combines_terms_and_removes_zero_path() {
        let algebra = PathAlgebra::new(Q, quiver()).unwrap();
        let a = algebra.quiver().arrow(arrow(1)).unwrap().path();
        let element = algebra
            .from_terms([
                (a.clone(), q(2, 1)),
                (algebra.quiver().zero(), q(7, 1)),
                (a, q(-2, 1)),
            ])
            .unwrap();
        assert!(element.is_zero());
    }

    #[test]
    fn multiplication_distributes_and_combines_equal_paths() {
        let algebra = PathAlgebra::new(Q, quiver()).unwrap();
        let a = algebra
            .element(algebra.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let v1 = algebra
            .element(algebra.quiver().vertex(vertex(1)).unwrap().path())
            .unwrap();
        let sum = v1.add(&a).unwrap();
        let square = sum.multiply(&sum).unwrap();
        assert_eq!(square.terms().len(), 3);
        assert_eq!(square.terms()[0].path().to_string(), "v1");
        assert_eq!(square.terms()[1].path().to_string(), "a");
        assert_eq!(square.terms()[1].coefficient(), &q(2, 1));
        assert_eq!(square.terms()[2].path().to_string(), "a^2");
    }

    #[test]
    fn leading_data_uses_the_last_qpa_sorted_term() {
        let algebra = PathAlgebra::new(Q, quiver()).unwrap();
        let a = algebra.quiver().arrow(arrow(1)).unwrap().path();
        let aa = a.multiply(&a).unwrap();
        let element = algebra
            .from_terms([(aa.clone(), q(3, 2)), (a, q(5, 1))])
            .unwrap();
        assert_eq!(element.leading_monomial(), aa);
        assert_eq!(element.leading_coefficient(), q(3, 2));
        assert_eq!(element.leading_term().terms().len(), 1);
    }

    #[test]
    fn uniformity_checks_path_sources_and_targets() {
        let algebra = PathAlgebra::new(Q, quiver()).unwrap();
        let a = algebra.quiver().arrow(arrow(1)).unwrap().path();
        let b = algebra.quiver().arrow(arrow(2)).unwrap().path();
        let element = algebra.from_terms([(a, q(1, 1)), (b, q(1, 1))]).unwrap();
        assert!(element.is_left_uniform());
        assert!(!element.is_right_uniform());
        assert!(!element.is_uniform());
        assert!(algebra.zero().is_uniform());
    }

    #[test]
    fn numerica_prime_fields_and_characteristic_two_share_sparse_algorithm() {
        let odd = Zp::new(5);
        let odd_algebra = PathAlgebra::new(odd.clone(), quiver()).unwrap();
        let a = odd_algebra.quiver().arrow(arrow(1)).unwrap().path();
        let odd_element = odd_algebra
            .from_terms([(a.clone(), odd.to_element(2)), (a, odd.to_element(3))])
            .unwrap();
        assert!(odd_element.is_zero());

        let binary_algebra = PathAlgebra::new(Z2, quiver()).unwrap();
        let a = binary_algebra.quiver().arrow(arrow(1)).unwrap().path();
        let one = Z2.nth(Integer::from(1));
        let binary_element = binary_algebra
            .from_terms([(a.clone(), one), (a, one)])
            .unwrap();
        assert!(binary_element.is_zero());
    }

    #[test]
    fn canonical_basis_is_qpa_path_iterator_for_acyclic_quivers() {
        let quiver = Quiver::from_adjacency(vec![vec![0, 1], vec![0, 0]]).unwrap();
        let algebra = PathAlgebra::new(Q, quiver).unwrap();
        let basis = algebra.canonical_basis().unwrap();
        assert_eq!(
            basis
                .iter()
                .map(|element| element.terms()[0].path().to_string())
                .collect::<Vec<_>>(),
            ["v1", "v2", "a1"]
        );
        assert_eq!(basis.len(), 3);
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_path_algebra_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/path_algebra.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1"
            })
        );

        let quiver = quiver();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let v1 = algebra
            .element(quiver.vertex(vertex(1)).unwrap().path())
            .unwrap();
        let a_path = quiver.arrow(arrow(1)).unwrap().path();
        let b_path = quiver.arrow(arrow(2)).unwrap().path();
        let aa_path = a_path.multiply(&a_path).unwrap();
        let a = algebra.element(a_path.clone()).unwrap();
        let b = algebra.element(b_path).unwrap();
        let aa = algebra.element(aa_path.clone()).unwrap();
        let x = v1.add(&a).unwrap();
        let square = x.multiply(&x).unwrap();
        let mixed = a
            .scale_left(&q(5, 1))
            .unwrap()
            .add(&aa.scale_left(&q(3, 2)).unwrap())
            .unwrap();
        let uniform = a.add(&b).unwrap();
        let zero = algebra.zero();

        assert_eq!(
            json!({
                "zero": rational_element_json(&zero),
                "one": rational_element_json(&algebra.one()),
                "x": rational_element_json(&x),
                "x_squared": rational_element_json(&square),
                "mixed": rational_element_json(&mixed),
                "zero_path_embedding": rational_element_json(
                    &algebra.element(quiver.zero()).unwrap()
                ),
                "noncomposable_product": rational_element_json(&b.multiply(&a).unwrap()),
                "leading": {
                    "term": rational_element_json(&mixed.leading_term()),
                    "coefficient": [
                        small_integer(mixed.leading_coefficient().numerator_ref()),
                        small_integer(mixed.leading_coefficient().denominator_ref())
                    ],
                    "monomial": mixed.leading_monomial().to_string(),
                    "zero_coefficient": [
                        small_integer(zero.leading_coefficient().numerator_ref()),
                        small_integer(zero.leading_coefficient().denominator_ref())
                    ],
                    "zero_monomial": zero.leading_monomial().to_string()
                },
                "uniformity": {
                    "zero": [zero.is_left_uniform(), zero.is_right_uniform(), zero.is_uniform()],
                    "a_plus_b": [uniform.is_left_uniform(), uniform.is_right_uniform(), uniform.is_uniform()]
                },
                "generator_terms": algebra.generators().iter().map(rational_element_json).collect::<Vec<_>>(),
                "global_dimension": algebra.qpa_global_dimension(),
                "is_hereditary": algebra.qpa_global_dimension() == 1
            }),
            expected["rational"]
        );

        let f5 = Zp::new(5);
        let algebra5 = PathAlgebra::new(f5.clone(), quiver.clone()).unwrap();
        let v15 = algebra5
            .element(quiver.vertex(vertex(1)).unwrap().path())
            .unwrap();
        let a5 = algebra5.element(a_path.clone()).unwrap();
        let aa5 = algebra5.element(aa_path.clone()).unwrap();
        let f5_sum = v15.add(&a5).unwrap();
        let f5_square = f5_sum.multiply(&f5_sum).unwrap();
        let f5_reduced = aa5
            .scale_left(&f5.to_element(4))
            .unwrap()
            .add(&aa5.scale_left(&f5.to_element(3)).unwrap())
            .unwrap();
        let f5_cancelled = a5
            .scale_left(&f5.to_element(2))
            .unwrap()
            .add(&a5.scale_left(&f5.to_element(3)).unwrap())
            .unwrap();

        let algebra2 = PathAlgebra::new(Z2, quiver).unwrap();
        let v12 = algebra2
            .element(algebra2.quiver().vertex(vertex(1)).unwrap().path())
            .unwrap();
        let a2 = algebra2.element(a_path).unwrap();
        let f2_sum = v12.add(&a2).unwrap();
        let f2_square = f2_sum.multiply(&f2_sum).unwrap();
        let f2_cancelled = a2.add(&a2).unwrap();
        assert_eq!(
            json!({
                "f5_square": zp_element_json(&f5_square, &f5),
                "f5_reduced": zp_element_json(&f5_reduced, &f5),
                "f5_cancelled": zp_element_json(&f5_cancelled, &f5),
                "f2_square": z2_element_json(&f2_square),
                "f2_cancelled": z2_element_json(&f2_cancelled)
            }),
            expected["prime_fields"]
        );

        let line = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "c"),
                ArrowSpec::named(vertex(2), vertex(3), "d"),
            ],
        )
        .unwrap();
        let line_algebra = PathAlgebra::new(Q, line).unwrap();
        let basis = line_algebra.canonical_basis().unwrap();
        assert_eq!(
            json!({
                "paths": basis.iter().map(|element| element.leading_monomial().to_string()).collect::<Vec<_>>(),
                "dimension": basis.len(),
                "global_dimension": line_algebra.qpa_global_dimension()
            }),
            expected["acyclic_basis"]
        );
    }
}
