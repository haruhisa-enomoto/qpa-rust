//! Ideals and quotients of QPA path algebras.
//!
//! This module translates the foundational behavior in `pathalgideal.gi` and
//! `pathalg.gi`: two-sided ideals remember their path algebra and generators,
//! membership is a Groebner normal-form test, and `KQ / relations` first
//! computes a Groebner basis and then stores quotient elements in normal form.
//! Direct factoring by an ideal without an attached basis preserves QPA's lazy
//! unreduced representative state.

use std::{
    cmp::Ordering,
    fmt,
    sync::{Arc, OnceLock, RwLock, Weak},
};

use numerica::domains::{Field, integer::Integer};
use thiserror::Error;

use crate::{
    DenseMatrix, GroebnerBasis, GroebnerError, MatrixError, Path, PathAlgebra, PathAlgebraElement,
    PathAlgebraError, PathAlgebraExternalRepresentation, PathEndpoint, Quiver, QuiverError,
    TensorProductAlgebra,
    path_module::{PathAlgebraModule, WeakPathAlgebraModule},
    tensor::WeakTensorProductAlgebra,
};

struct PathAlgebraIdealData<F: Field> {
    algebra: PathAlgebra<F>,
    generators: Vec<PathAlgebraElement<F>>,
    groebner_basis: RwLock<Option<GroebnerBasis<F>>>,
}

/// A QPA two-sided ideal in a path algebra.
///
/// Cloning preserves ideal identity. Constructing another ideal with equal
/// generators creates a distinct ideal object, matching GAP object identity.
#[derive(Clone)]
pub struct PathAlgebraIdeal<F: Field>(Arc<PathAlgebraIdealData<F>>);

impl<F: Field> PathAlgebraIdeal<F> {
    /// Translate QPA's `TwoSidedIdealByGenerators`.
    ///
    /// # Errors
    ///
    /// Returns an error if a generator belongs to another path algebra.
    pub fn new(
        algebra: PathAlgebra<F>,
        generators: Vec<PathAlgebraElement<F>>,
    ) -> Result<Self, QuotientError> {
        ensure_parent_elements(&algebra, &generators)?;
        Ok(Self(Arc::new(PathAlgebraIdealData {
            algebra,
            generators,
            groebner_basis: RwLock::new(None),
        })))
    }

    fn from_generators_and_groebner_basis(
        algebra: PathAlgebra<F>,
        generators: Vec<PathAlgebraElement<F>>,
        basis: GroebnerBasis<F>,
    ) -> Self {
        Self(Arc::new(PathAlgebraIdealData {
            algebra,
            generators,
            groebner_basis: RwLock::new(Some(basis)),
        }))
    }

    /// Return the containing path algebra.
    #[must_use]
    pub fn algebra(&self) -> &PathAlgebra<F> {
        &self.0.algebra
    }

    /// Return QPA's `GeneratorsOfIdeal` list.
    #[must_use]
    pub fn generators(&self) -> &[PathAlgebraElement<F>] {
        &self.0.generators
    }

    /// Return whether two handles denote the identical QPA ideal object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Return QPA's explicitly handled trivial-ideal condition.
    #[must_use]
    pub fn is_trivial(&self) -> bool {
        self.generators().iter().all(PathAlgebraElement::is_zero)
    }

    /// Compute and attach QPA's high-level Groebner basis.
    ///
    /// The basis is attached to the same identity-preserving ideal object, as
    /// GAP's cached `GroebnerBasisOfIdeal` attribute is.
    ///
    /// # Errors
    ///
    /// Propagates QPA high-level Groebner failures.
    pub fn with_high_level_groebner_basis(&self) -> Result<Self, QuotientError> {
        if self.groebner_basis().is_some() {
            return Ok(self.clone());
        }
        let basis = GroebnerBasis::high_level(self.algebra(), self.generators())?;
        self.set_groebner_basis(basis);
        Ok(self.clone())
    }

    /// Attach QPA `GroebnerBasis(I, relations)` to this ideal.
    ///
    /// The returned handle denotes the same ideal object; all clones observe
    /// the attached complete basis while the generator list remains unchanged.
    ///
    /// # Errors
    ///
    /// Returns an error for foreign relations or QPA comparison failures.
    pub fn with_complete_groebner_basis(
        &self,
        relations: &[PathAlgebraElement<F>],
    ) -> Result<Self, QuotientError> {
        let basis = GroebnerBasis::from_complete_relations(self.algebra(), relations)?;
        self.set_groebner_basis(basis);
        Ok(self.clone())
    }

    /// Return an attached Groebner basis, if one has been computed.
    #[must_use]
    pub fn groebner_basis(&self) -> Option<GroebnerBasis<F>> {
        self.0
            .groebner_basis
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone()
    }

    fn set_groebner_basis(&self, basis: GroebnerBasis<F>) {
        *self
            .0
            .groebner_basis
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(basis);
    }

    /// Translate ideal-level QPA `RightGroebnerBasis(I)`.
    ///
    /// QPA installs this method only when `GroebnerBasisOfIdeal(I)` is already
    /// known; Rust reports that missing dispatch precondition explicitly.
    ///
    /// # Errors
    ///
    /// Returns an error if no basis is attached, the quotient has infinitely
    /// many nontips, or source-translated reduction fails.
    pub fn right_groebner_basis(&self) -> Result<GroebnerBasis<F>, QuotientError> {
        Ok(self
            .groebner_basis()
            .ok_or(QuotientError::MissingGroebnerBasis)?
            .right_basis()?)
    }

    /// Translate QPA ideal membership as complete Groebner reduction to zero.
    ///
    /// If the basis has not yet been attached, QPA's high-level basis is
    /// computed for this test, reflecting QPA's lazy membership behavior.
    ///
    /// # Errors
    ///
    /// Returns an error for a foreign element or Groebner failure.
    pub fn contains(&self, element: PathAlgebraElement<F>) -> Result<bool, QuotientError> {
        ensure_parent_element(self.algebra(), &element)?;
        if element.is_zero() {
            return Ok(true);
        }
        if self.is_trivial() {
            return Ok(false);
        }
        let basis = match self.groebner_basis() {
            Some(basis) => basis,
            None => GroebnerBasis::high_level(self.algebra(), self.generators())?,
        }
        .completely_reduced_basis()?;
        self.set_groebner_basis(basis.clone());
        Ok(basis.completely_reduce(element)?.is_zero())
    }

    /// Translate QPA's fast monomial-ideal test on ideal generators.
    #[must_use]
    pub fn generators_are_monomial(&self) -> bool {
        self.generators()
            .iter()
            .all(|generator| generator.terms().len() == 1)
    }

    /// Translate QPA's full `IsMonomialIdeal` method.
    ///
    /// QPA first checks the supplied generators, then computes a Gröbner basis
    /// only when a nonmonomial generator occurs.
    ///
    /// # Errors
    ///
    /// Propagates source-translated Gröbner failures.
    pub fn is_monomial(&self) -> Result<bool, QuotientError> {
        if self.is_trivial() || self.generators_are_monomial() {
            return Ok(true);
        }
        let basis = if let Some(basis) = self.groebner_basis() {
            basis
        } else {
            let basis = GroebnerBasis::high_level(self.algebra(), self.generators())?;
            self.set_groebner_basis(basis.clone());
            basis
        };
        Ok(basis
            .relations()
            .iter()
            .all(|relation| relation.terms().len() == 1))
    }

    /// Translate QPA's `IsAdmissibleIdeal` method.
    ///
    /// After checking `I \subseteq R^2`, this constructs `A/I`, forms the
    /// two-sided ideal generated by quotient arrows, and repeatedly replaces
    /// it by its square until its dimension stabilizes, exactly as QPA does.
    ///
    /// # Errors
    ///
    /// Propagates Gröbner, quotient, matrix, and path-arithmetic failures.
    pub fn is_admissible(&self) -> Result<bool, QuotientError> {
        if self.is_trivial() {
            return Ok(self.algebra().has_finite_canonical_basis());
        }

        let ideal = self.with_high_level_groebner_basis()?;
        if ideal.generators().iter().any(|relation| {
            relation
                .terms()
                .iter()
                .any(|term| term.path().length().is_none_or(|length| length < 2))
        }) {
            return Ok(false);
        }

        let quotient = QuotientPathAlgebra::by_ideal_high_level(&ideal)?;
        let Some(algebra_basis) = quotient.canonical_basis()? else {
            return Ok(false);
        };
        let one = quotient.one()?;
        let mut arrow_generators = Vec::new();
        for arrow in quotient.original_path_algebra().quiver().arrows() {
            let arrow = quotient.path(arrow.path())?;
            arrow_generators.push(one.multiply(&arrow)?);
        }
        let arrow_ideal = quotient_ideal_span(&quotient, &arrow_generators, &algebra_basis)?;
        let mut dimension = arrow_ideal.len();
        let mut power =
            quotient_product_space(&quotient, &arrow_ideal, &arrow_ideal, &algebra_basis)?;
        while dimension != power.len() {
            dimension = power.len();
            power = quotient_product_space(&quotient, &power, &power, &algebra_basis)?;
        }
        Ok(dimension == 0)
    }

    /// Translate QPA's `ProductOfIdeals(I, J)` with `self` as `I`.
    ///
    /// QPA requires an already attached complete basis for `J`, constructs its
    /// right Gröbner basis, and returns the ideal generated in nested
    /// `Generators(I)`/right-basis order.
    ///
    /// # Errors
    ///
    /// Returns an error for different parent algebras, a missing basis,
    /// infinitely many right-factor nontips, or arithmetic failure.
    pub fn product(&self, right: &Self) -> Result<Self, QuotientError> {
        if !self.algebra().is_identical(right.algebra()) {
            return Err(QuotientError::DifferentParentAlgebras);
        }
        let complete_basis = right
            .groebner_basis()
            .ok_or(QuotientError::MissingGroebnerBasis)?;
        if !complete_basis.admits_finitely_many_nontips()? {
            return Err(QuotientError::InfiniteRightFactorNontips);
        }
        let right_basis = complete_basis.right_basis()?;
        let mut generators = Vec::new();
        for left_generator in self.generators() {
            for right_generator in right_basis.relations() {
                generators.push(left_generator.multiply(right_generator)?);
            }
        }
        Self::new(self.algebra().clone(), generators)
    }

    /// Translate QPA's `NaturalHomomorphismByIdeal(A, I)` with
    /// `A = self.algebra()`.
    ///
    /// This does not compute a Gröbner basis. If none is already attached to
    /// the ideal, images retain QPA's lazy unreduced representatives.
    ///
    /// # Errors
    ///
    /// Propagates factor-quotient construction failures.
    pub fn natural_homomorphism(&self) -> Result<NaturalQuotientMap<F>, QuotientError> {
        NaturalQuotientMap::new(self)
    }

    /// Translate QPA's `MinimalGeneratingSetOfIdeal` for an admissible ideal.
    ///
    /// This forms `JI + IJ`, constructs `KQ/(JI+IJ)`, computes the GAP-style
    /// semi-echelon basis of the image ideal, and pulls its representatives
    /// back to `KQ`.
    ///
    /// # Errors
    ///
    /// Requires an admissible ideal and finite basis computations; all exact
    /// arithmetic and Gröbner failures are propagated.
    pub fn minimal_generating_set(&self) -> Result<Vec<PathAlgebraElement<F>>, QuotientError> {
        self.minimal_generating_set_with(|parent, relations| {
            QuotientPathAlgebra::by_relations_high_level(parent, relations)
        })
    }

    /// Execute QPA's `MinimalGeneratingSetOfIdeal` with an explicitly selected
    /// relation-list quotient backend.
    ///
    /// QPA routes the internal `KQ/(JI+IJ)` construction through the
    /// `GroebnerBasisFunction` installed on `KQ`. This hook lets the GBNP
    /// facade preserve that dispatch while retaining the backend-neutral
    /// linear-algebra part here.
    ///
    /// # Errors
    ///
    /// Requires an admissible ideal and finite basis computations; backend and
    /// exact-arithmetic failures are converted to `E`.
    pub fn minimal_generating_set_with<E, Factor>(
        &self,
        factor: Factor,
    ) -> Result<Vec<PathAlgebraElement<F>>, E>
    where
        E: From<QuotientError>,
        Factor:
            FnOnce(PathAlgebra<F>, &[PathAlgebraElement<F>]) -> Result<QuotientPathAlgebra<F>, E>,
    {
        if self.generators().is_empty() {
            return Ok(Vec::new());
        }
        if !self.is_admissible().map_err(E::from)? {
            return Err(E::from(QuotientError::IdealNotAdmissible));
        }

        let arrows = self
            .algebra()
            .quiver()
            .arrows()
            .map(|arrow| self.algebra().element(arrow.path()))
            .collect::<Result<Vec<_>, _>>()
            .map_err(QuotientError::from)
            .map_err(E::from)?;
        let mut radical_products = Vec::new();
        for arrow in &arrows {
            for generator in self.generators() {
                let product = generator
                    .multiply(arrow)
                    .map_err(QuotientError::from)
                    .map_err(E::from)?;
                if !product.is_zero() {
                    radical_products.push(product);
                }
            }
        }
        for arrow in &arrows {
            for generator in self.generators() {
                let product = arrow
                    .multiply(generator)
                    .map_err(QuotientError::from)
                    .map_err(E::from)?;
                if !product.is_zero() {
                    radical_products.push(product);
                }
            }
        }

        if radical_products.is_empty() {
            let algebra_basis = self
                .algebra()
                .canonical_basis()
                .ok_or(QuotientError::InfiniteParentBasis)
                .map_err(E::from)?;
            return path_algebra_ideal_span(self.algebra(), self.generators(), &algebra_basis)
                .map_err(E::from);
        }

        let quotient = factor(self.algebra().clone(), &radical_products)?;
        let algebra_basis = quotient
            .canonical_basis()
            .map_err(E::from)?
            .ok_or(QuotientError::InfiniteQuotientBasis)
            .map_err(E::from)?;
        let generators = self
            .generators()
            .iter()
            .cloned()
            .map(|generator| quotient.element(generator))
            .collect::<Result<Vec<_>, _>>()
            .map_err(E::from)?;
        Ok(quotient_ideal_span(&quotient, &generators, &algebra_basis)
            .map_err(E::from)?
            .into_iter()
            .map(|element| element.representative().clone())
            .collect())
    }
}

impl<F: Field> fmt::Debug for PathAlgebraIdeal<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PathAlgebraIdeal")
            .field("algebra", self.algebra())
            .field("generators", &self.generators())
            .field("has_groebner_basis", &self.groebner_basis().is_some())
            .finish()
    }
}

impl<F: Field> PartialEq for PathAlgebraIdeal<F> {
    fn eq(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}

impl<F: Field> Eq for PathAlgebraIdeal<F> {}

struct QuotientPathAlgebraData<F: Field> {
    parent: PathAlgebra<F>,
    ideal: PathAlgebraIdeal<F>,
    groebner_basis: Option<GroebnerBasis<F>>,
    admissible: OnceLock<bool>,
    canonical_basis: OnceLock<QuotientCanonicalBasisState<F>>,
    opposite_weak: RwLock<Option<Weak<QuotientPathAlgebraData<F>>>>,
    opposite_strong: OnceLock<QuotientPathAlgebra<F>>,
    enveloping_weak: RwLock<Option<WeakTensorProductAlgebra<F>>>,
    regular_enveloping_module_weak: RwLock<Option<WeakPathAlgebraModule<F>>>,
    regular_enveloping_module_bound: OnceLock<()>,
    regular_enveloping_module_fail: OnceLock<()>,
    dual_enveloping_module_weak: RwLock<Option<WeakPathAlgebraModule<F>>>,
    dual_enveloping_module_bound: OnceLock<()>,
    dual_enveloping_module_fail: OnceLock<()>,
    is_enveloping: OnceLock<bool>,
    trivial_extension_level: OnceLock<usize>,
    trivial_extension_weak: RwLock<Option<WeakQuotientPathAlgebra<F>>>,
    trivial_extension_bound: OnceLock<()>,
    trivial_extension_fail: OnceLock<()>,
}

#[derive(Debug)]
enum QuotientCanonicalBasisState<F: Field> {
    WithoutVectors,
    WithVectors {
        representatives: Vec<PathAlgebraElement<F>>,
        is_normal_form: bool,
    },
}

/// A quotient `KQ/I` with QPA's optional normal-form function.
#[derive(Clone)]
pub struct QuotientPathAlgebra<F: Field>(Arc<QuotientPathAlgebraData<F>>);

#[derive(Clone, Debug)]
pub(crate) struct WeakQuotientPathAlgebra<F: Field>(Weak<QuotientPathAlgebraData<F>>);

impl<F: Field> WeakQuotientPathAlgebra<F> {
    pub(crate) fn upgrade(&self) -> Option<QuotientPathAlgebra<F>> {
        self.0.upgrade().map(QuotientPathAlgebra)
    }
}

impl<F: Field> QuotientPathAlgebra<F> {
    /// Translate QPA's `KQ / relators` using its `HighLevelGroebnerBasis`
    /// backend.
    ///
    /// QPA replaces the supplied relation list by the computed Groebner basis
    /// when it creates the defining ideal; [`Self::relations`] therefore
    /// returns that basis.
    ///
    /// # Errors
    ///
    /// Returns an error for foreign relators or Groebner failures.
    pub fn by_relations_high_level(
        parent: PathAlgebra<F>,
        relators: &[PathAlgebraElement<F>],
    ) -> Result<Self, QuotientError> {
        ensure_parent_elements(&parent, relators)?;
        let complete = GroebnerBasis::high_level(&parent, relators)?;
        let complete_relations = complete.relations().to_vec();
        Self::by_relators_and_complete_groebner_basis(
            parent,
            &complete_relations,
            &complete_relations,
        )
    }

    /// Construct QPA's quotient after an external source-translated backend
    /// has computed a complete Gröbner relation list.
    ///
    /// This is the cross-crate boundary used by the GBNP translation. It
    /// models QPA's `/` method on a relation list, which replaces the defining
    /// ideal generators by the computed basis before factoring.
    ///
    /// # Errors
    ///
    /// Returns an error for foreign relations or basis normalization failure.
    pub fn by_complete_groebner_relations(
        parent: PathAlgebra<F>,
        complete_relations: &[PathAlgebraElement<F>],
    ) -> Result<Self, QuotientError> {
        Self::by_relators_and_complete_groebner_basis(
            parent,
            complete_relations,
            complete_relations,
        )
    }

    /// Translate `FactorPathAlgebraByRelators(A, I, OrderingOfAlgebra(A))`
    /// without forcing a Gröbner computation.
    ///
    /// If the ideal already carries a complete basis, the quotient installs
    /// QPA's normal-form function. Otherwise quotient elements retain their
    /// representatives and are not marked as normal forms.
    ///
    /// # Errors
    ///
    /// The ideal must belong to its recorded parent path algebra.
    pub fn by_ideal(ideal: &PathAlgebraIdeal<F>) -> Result<Self, QuotientError> {
        let groebner_basis = ideal
            .groebner_basis()
            .map(|basis| basis.tip_reduced_basis())
            .transpose()?;
        if let Some(basis) = &groebner_basis {
            ideal.set_groebner_basis(basis.clone());
        }
        Ok(Self(Arc::new(QuotientPathAlgebraData {
            parent: ideal.algebra().clone(),
            ideal: ideal.clone(),
            groebner_basis,
            admissible: OnceLock::new(),
            canonical_basis: OnceLock::new(),
            opposite_weak: RwLock::new(None),
            opposite_strong: OnceLock::new(),
            enveloping_weak: RwLock::new(None),
            regular_enveloping_module_weak: RwLock::new(None),
            regular_enveloping_module_bound: OnceLock::new(),
            regular_enveloping_module_fail: OnceLock::new(),
            dual_enveloping_module_weak: RwLock::new(None),
            dual_enveloping_module_bound: OnceLock::new(),
            dual_enveloping_module_fail: OnceLock::new(),
            is_enveloping: OnceLock::new(),
            trivial_extension_level: OnceLock::new(),
            trivial_extension_weak: RwLock::new(None),
            trivial_extension_bound: OnceLock::new(),
            trivial_extension_fail: OnceLock::new(),
        })))
    }

    /// Factor by an ideal, computing QPA's high-level basis when necessary.
    ///
    /// # Errors
    ///
    /// Propagates Groebner and parent-family failures.
    pub fn by_ideal_high_level(ideal: &PathAlgebraIdeal<F>) -> Result<Self, QuotientError> {
        let parent = ideal.algebra().clone();
        let ideal = ideal.with_high_level_groebner_basis()?;
        let groebner_basis = ideal
            .groebner_basis()
            .ok_or(QuotientError::MissingGroebnerBasis)?
            .tip_reduced_basis()?;
        ideal.set_groebner_basis(groebner_basis.clone());
        Ok(Self(Arc::new(QuotientPathAlgebraData {
            parent,
            ideal,
            groebner_basis: Some(groebner_basis),
            admissible: OnceLock::new(),
            canonical_basis: OnceLock::new(),
            opposite_weak: RwLock::new(None),
            opposite_strong: OnceLock::new(),
            enveloping_weak: RwLock::new(None),
            regular_enveloping_module_weak: RwLock::new(None),
            regular_enveloping_module_bound: OnceLock::new(),
            regular_enveloping_module_fail: OnceLock::new(),
            dual_enveloping_module_weak: RwLock::new(None),
            dual_enveloping_module_bound: OnceLock::new(),
            dual_enveloping_module_fail: OnceLock::new(),
            is_enveloping: OnceLock::new(),
            trivial_extension_level: OnceLock::new(),
            trivial_extension_weak: RwLock::new(None),
            trivial_extension_bound: OnceLock::new(),
            trivial_extension_fail: OnceLock::new(),
        })))
    }

    /// Return the original path algebra.
    #[must_use]
    pub fn original_path_algebra(&self) -> &PathAlgebra<F> {
        &self.0.parent
    }

    /// Return QPA's defining ideal.
    #[must_use]
    pub fn ideal(&self) -> &PathAlgebraIdeal<F> {
        &self.0.ideal
    }

    /// Return `RelationsOfAlgebra`, namely the defining ideal generators.
    #[must_use]
    pub fn relations(&self) -> &[PathAlgebraElement<F>] {
        self.ideal().generators()
    }

    /// Return the complete Groebner basis used for normal forms.
    #[must_use]
    pub fn groebner_basis(&self) -> Option<&GroebnerBasis<F>> {
        self.0.groebner_basis.as_ref()
    }

    /// Return whether QPA installed a normal-form function on the quotient
    /// element family.
    #[must_use]
    pub fn has_normal_form_function(&self) -> bool {
        self.groebner_basis().is_some()
    }

    pub(crate) fn trivial_extension_level(&self) -> usize {
        self.0.trivial_extension_level.get().copied().unwrap_or(0)
    }

    pub(crate) fn set_trivial_extension_level(&self, level: usize) {
        let _ = self.0.trivial_extension_level.set(level);
    }

    pub(crate) fn cached_trivial_extension(&self) -> Option<Self> {
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

    pub(crate) fn cache_trivial_extension(&self, extension: &Self) -> Self {
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

    pub(crate) fn downgrade(&self) -> WeakQuotientPathAlgebra<F> {
        WeakQuotientPathAlgebra(Arc::downgrade(&self.0))
    }

    /// Return whether two handles denote the identical quotient object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Return the quotient zero.
    #[must_use]
    pub fn zero(&self) -> QuotientPathAlgebraElement<F> {
        QuotientPathAlgebraElement {
            algebra: self.clone(),
            representative: self.original_path_algebra().zero(),
            is_normal_form: self.has_normal_form_function(),
        }
    }

    /// Return the quotient identity.
    ///
    /// # Errors
    ///
    /// Propagates normal-form failures.
    pub fn one(&self) -> Result<QuotientPathAlgebraElement<F>, QuotientError> {
        self.element_with_normality(self.original_path_algebra().one(), true)
    }

    /// Map a parent path-algebra element into the quotient and normalize it.
    ///
    /// # Errors
    ///
    /// Returns an error for a foreign parent algebra or reduction failure.
    pub fn element(
        &self,
        element: PathAlgebraElement<F>,
    ) -> Result<QuotientPathAlgebraElement<F>, QuotientError> {
        self.element_with_normality(element, false)
    }

    /// Embed a quiver path into the quotient and normalize it.
    ///
    /// # Errors
    ///
    /// Propagates path-family and normal-form failures.
    pub fn path(&self, path: Path) -> Result<QuotientPathAlgebraElement<F>, QuotientError> {
        self.element(self.original_path_algebra().element(path)?)
    }

    /// Return quotient generators corresponding to parent vertices and arrows.
    ///
    /// # Errors
    ///
    /// Propagates normal-form failures.
    pub fn generators(&self) -> Result<Vec<QuotientPathAlgebraElement<F>>, QuotientError> {
        self.original_path_algebra()
            .generators()
            .into_iter()
            .map(|generator| self.element(generator))
            .collect()
    }

    /// Translate QPA's `IsFiniteDimensional` using the complete tip
    /// dictionary.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary failures.
    pub fn is_finite_dimensional(&self) -> Result<bool, QuotientError> {
        let basis = self
            .ideal()
            .groebner_basis()
            .ok_or(QuotientError::MissingGroebnerBasis)?;
        if basis.is_complete() {
            let basis = basis.tip_reduced_basis()?;
            self.ideal().set_groebner_basis(basis.clone());
            Ok(basis.admits_finitely_many_nontips()?)
        } else if self.original_path_algebra().has_finite_canonical_basis() {
            Ok(true)
        } else {
            Err(GroebnerError::IncompleteBasis.into())
        }
    }

    /// Translate QPA's `Dimension`; infinite dimension is explicit.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary failures.
    pub fn dimension(&self) -> Result<AlgebraDimension, QuotientError> {
        let basis = self
            .ideal()
            .groebner_basis()
            .ok_or(QuotientError::MissingGroebnerBasis)?
            .tip_reduced_basis()?;
        self.ideal().set_groebner_basis(basis.clone());
        Ok(basis
            .nontip_size()?
            .map_or(AlgebraDimension::Infinite, AlgebraDimension::Finite))
    }

    /// Translate QPA's canonical basis of a finite-dimensional quotient.
    ///
    /// This is a convenience view of [`Self::canonical_basis_object`]. `None`
    /// means QPA's basis object has no bound `BasisVectors`: its ideal lacks a
    /// complete basis or its nontip language is infinite.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary, embedding, and normal-form failures.
    pub fn canonical_basis(
        &self,
    ) -> Result<Option<Vec<QuotientPathAlgebraElement<F>>>, QuotientError> {
        Ok(self.canonical_basis_object()?.into_basis_vectors())
    }

    /// Translate QPA's `CanonicalBasis` object, including its unbound-vector
    /// state before a complete finite Gröbner basis is available.
    ///
    /// # Errors
    ///
    /// Propagates nontip dictionary and quotient element construction errors.
    pub fn canonical_basis_object(&self) -> Result<QuotientCanonicalBasis<F>, QuotientError> {
        if self.0.canonical_basis.get().is_none() {
            let state = match self.ideal().groebner_basis() {
                Some(basis) if basis.is_complete() => {
                    let basis = basis.tip_reduced_basis()?;
                    self.ideal().set_groebner_basis(basis.clone());
                    if basis.admits_finitely_many_nontips()? {
                        let Some(mut paths) = basis.nontips()? else {
                            return Err(QuotientError::InvalidCanonicalBasis);
                        };
                        sort_paths_as_qpa_set(&mut paths)?;
                        let representatives = paths
                            .into_iter()
                            .map(|path| Ok(self.original_path_algebra().element(path)?))
                            .collect::<Result<Vec<_>, QuotientError>>()?;
                        QuotientCanonicalBasisState::WithVectors {
                            representatives,
                            is_normal_form: self.has_normal_form_function(),
                        }
                    } else {
                        QuotientCanonicalBasisState::WithoutVectors
                    }
                }
                _ => QuotientCanonicalBasisState::WithoutVectors,
            };
            let _ = self.0.canonical_basis.set(state);
        }
        Ok(QuotientCanonicalBasis {
            quotient: self.clone(),
        })
    }

    /// Translate `Coefficients` for QPA's canonical quotient basis.
    ///
    /// # Errors
    ///
    /// The element must belong to this quotient and carry QPA's normal-form
    /// filter. The defining ideal must carry a complete finite-nontip basis.
    pub fn canonical_basis_coefficients(
        &self,
        element: &QuotientPathAlgebraElement<F>,
    ) -> Result<Vec<F::Element>, QuotientError> {
        ensure_quotient_elements(self, std::slice::from_ref(element))?;
        if !element.is_normal_form() {
            return Err(QuotientError::ElementsNotInNormalForm);
        }
        let basis = self
            .ideal()
            .groebner_basis()
            .ok_or(QuotientError::MissingGroebnerBasis)?
            .tip_reduced_basis()?;
        self.ideal().set_groebner_basis(basis.clone());
        let mut paths = basis
            .nontips()?
            .ok_or(QuotientError::InfiniteQuotientBasis)?;
        sort_paths_as_qpa_set(&mut paths)?;
        element_coordinates(element.representative(), &paths)
    }

    /// Translate QPA's quotient `ObjByExtRep` conversion.
    ///
    /// The parent descriptor is decoded first and then passed to the same
    /// lazy-or-normalizing quotient constructor used by QPA.
    ///
    /// # Errors
    ///
    /// Propagates external-representation validation and normal-form failures.
    pub fn element_from_external_representation(
        &self,
        representation: &PathAlgebraExternalRepresentation<F::Element>,
    ) -> Result<QuotientPathAlgebraElement<F>, QuotientError> {
        let representative = self
            .original_path_algebra()
            .element_from_external_representation(representation)?;
        self.element(representative)
    }

    /// Translate QPA's `IsAdmissibleQuotientOfPathAlgebra`.
    ///
    /// # Errors
    ///
    /// Propagates the defining ideal's admissibility computation.
    pub fn is_admissible(&self) -> Result<bool, QuotientError> {
        if let Some(admissible) = self.0.admissible.get() {
            return Ok(*admissible);
        }
        let admissible = self.ideal().is_admissible()?;
        let _ = self.0.admissible.set(admissible);
        Ok(admissible)
    }

    /// Return QPA's cached `IsAdmissibleQuotientOfPathAlgebra` property.
    ///
    /// `None` corresponds to `HasIsAdmissibleQuotientOfPathAlgebra(A) = false`.
    #[must_use]
    pub fn cached_admissibility(&self) -> Option<bool> {
        self.0.admissible.get().copied()
    }

    /// Install QPA's true admissible-quotient filter without recomputation.
    #[doc(hidden)]
    pub fn mark_admissible(&self) {
        let _ = self.0.admissible.set(true);
    }

    /// Return QPA's `IsEnvelopingAlgebra` property value.
    #[must_use]
    pub fn is_enveloping_algebra(&self) -> bool {
        self.0.is_enveloping.get().copied().unwrap_or(false)
    }

    /// Install QPA's true enveloping-algebra property.
    #[doc(hidden)]
    pub fn mark_enveloping_algebra(&self) {
        let _ = self.0.is_enveloping.set(true);
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

    /// Translate QPA's cached `OppositePathAlgebra` for a quotient.
    ///
    /// QPA takes opposite relators, computes the configured Gröbner basis on
    /// the opposite parent, and factors by that ideal. This core entry point
    /// uses the already source-translated high-level QPA backend; the default
    /// GBNP-backed facade is supplied by the GBNP layer.
    ///
    /// # Errors
    ///
    /// Propagates opposite-parent, relation, and Gröbner failures.
    pub fn opposite_high_level(&self) -> Result<Self, QuotientError> {
        if let Some(opposite) = self.cached_opposite() {
            return Ok(opposite);
        }

        let opposite_parent = self.original_path_algebra().opposite()?;
        let opposite_relations = self.opposite_defining_relations()?;
        let complete = GroebnerBasis::high_level(&opposite_parent, &opposite_relations)?;
        self.install_opposite_from_complete_relations(&opposite_relations, complete.relations())
    }

    /// Return an already-live cached opposite quotient, if one exists.
    ///
    /// This backend-neutral hook lets the GBNP crate avoid recomputing a basis
    /// before delegating reciprocal-cache installation back to `qpa-core`.
    #[doc(hidden)]
    #[must_use]
    pub fn cached_opposite(&self) -> Option<Self> {
        self.0.opposite_strong.get().cloned().or_else(|| {
            self.0
                .opposite_weak
                .read()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .as_ref()
                .and_then(Weak::upgrade)
                .map(Self)
        })
    }

    /// Return QPA's termwise opposite defining relators.
    ///
    /// # Errors
    ///
    /// Propagates opposite path-algebra element failures.
    #[doc(hidden)]
    pub fn opposite_defining_relations(&self) -> Result<Vec<PathAlgebraElement<F>>, QuotientError> {
        self.relations()
            .iter()
            .map(PathAlgebraElement::opposite)
            .collect::<Result<Vec<_>, _>>()
            .map_err(Into::into)
    }

    /// Complete QPA's opposite-quotient construction after a backend has
    /// computed the opposite ideal's complete Gröbner basis.
    ///
    /// # Errors
    ///
    /// Returns an error for foreign relators, invalid basis data, or a lost
    /// reciprocal-cache race.
    #[doc(hidden)]
    pub fn install_opposite_from_complete_relations(
        &self,
        opposite_relators: &[PathAlgebraElement<F>],
        complete_relations: &[PathAlgebraElement<F>],
    ) -> Result<Self, QuotientError> {
        if let Some(opposite) = self.cached_opposite() {
            return Ok(opposite);
        }
        let opposite_parent = self.original_path_algebra().opposite()?;
        let opposite = Self::by_relators_and_complete_groebner_basis(
            opposite_parent,
            opposite_relators,
            complete_relations,
        )?;
        if opposite.0.opposite_strong.set(self.clone()).is_err() {
            return Err(QuotientError::OppositeCacheRace);
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

    /// Translate QPA's `MakeUniformOnRight` for quotient elements.
    ///
    /// The operation examines each stored representative, preserving QPA's
    /// behavior for quotient families without a normal-form function.
    ///
    /// # Errors
    ///
    /// All elements must belong to this identical quotient family.
    pub fn make_uniform_on_right(
        &self,
        generators: &[QuotientPathAlgebraElement<F>],
    ) -> Result<Vec<QuotientPathAlgebraElement<F>>, QuotientError> {
        ensure_quotient_elements(self, generators)?;
        let mut uniform = Vec::new();
        for generator in generators {
            let mut targets = Vec::new();
            for term in generator.representative().terms() {
                if let PathEndpoint::Vertex(target) = term.path().target()
                    && !targets.contains(&target)
                {
                    targets.push(target);
                }
            }
            for target in targets {
                let vertex = self.path(target.path())?;
                let component = generator.multiply(&vertex)?;
                if !component.is_zero() && !uniform.contains(&component) {
                    uniform.push(component);
                }
            }
        }
        Ok(uniform)
    }

    /// Translate QPA's `GeneratorsTimesArrowsOnRight` for quotient elements.
    ///
    /// # Errors
    ///
    /// Propagates uniformization, quiver lookup, and quotient arithmetic
    /// failures.
    pub fn generators_times_arrows_on_right(
        &self,
        generators: &[QuotientPathAlgebraElement<F>],
    ) -> Result<Vec<QuotientPathAlgebraElement<F>>, QuotientError> {
        let uniform = self.make_uniform_on_right(generators)?;
        let mut extended = Vec::new();
        for generator in uniform {
            let PathEndpoint::Vertex(target) = generator.leading_monomial().target() else {
                continue;
            };
            for arrow in self
                .original_path_algebra()
                .quiver()
                .outgoing_arrows(&target)?
            {
                let product = generator.multiply(&self.path(arrow.path())?)?;
                if !product.is_zero() {
                    extended.push(product);
                }
            }
        }
        Ok(extended)
    }

    /// Construct QPA's truncated path algebra with the high-level Gröbner
    /// backend available in `qpa-core`.
    ///
    /// # Errors
    ///
    /// Propagates path-algebra, arrow-power, and Gröbner construction failures.
    pub fn truncated_high_level(
        field: F,
        quiver: Quiver,
        exponent: i64,
    ) -> Result<Self, QuotientError> {
        let parent = PathAlgebra::new(field, quiver)?;
        let relations = parent.nth_power_of_arrow_ideal(exponent)?;
        Self::by_relations_high_level(parent, &relations)
    }

    /// Translate QPA's `ElementOfQuotientOfPathAlgebra(family, element,
    /// normal)` constructor.
    ///
    /// When a normal-form function is installed, `known_normal = true` trusts
    /// the caller and marks the representative; `false` reduces it. A lazy
    /// quotient ignores the flag, exactly as QPA's less-specific method does.
    ///
    /// # Errors
    ///
    /// Returns an error for a foreign parent or failed complete reduction.
    pub fn element_with_normality(
        &self,
        element: PathAlgebraElement<F>,
        known_normal: bool,
    ) -> Result<QuotientPathAlgebraElement<F>, QuotientError> {
        ensure_parent_element(self.original_path_algebra(), &element)?;
        let (representative, is_normal_form) = match self.groebner_basis() {
            Some(_) if known_normal => (element, true),
            Some(basis) => (basis.completely_reduce(element)?, true),
            None => (element, false),
        };
        Ok(QuotientPathAlgebraElement {
            algebra: self.clone(),
            representative,
            is_normal_form,
        })
    }

    fn by_relators_and_complete_groebner_basis(
        parent: PathAlgebra<F>,
        relators: &[PathAlgebraElement<F>],
        complete_relations: &[PathAlgebraElement<F>],
    ) -> Result<Self, QuotientError> {
        ensure_parent_elements(&parent, relators)?;
        ensure_parent_elements(&parent, complete_relations)?;
        let groebner_basis = GroebnerBasis::from_complete_relations(&parent, complete_relations)?
            .tip_reduced_basis()?;
        let ideal = PathAlgebraIdeal::from_generators_and_groebner_basis(
            parent.clone(),
            relators.to_vec(),
            groebner_basis.clone(),
        );
        let quotient = Self(Arc::new(QuotientPathAlgebraData {
            parent,
            ideal,
            groebner_basis: Some(groebner_basis),
            admissible: OnceLock::new(),
            canonical_basis: OnceLock::new(),
            opposite_weak: RwLock::new(None),
            opposite_strong: OnceLock::new(),
            enveloping_weak: RwLock::new(None),
            regular_enveloping_module_weak: RwLock::new(None),
            regular_enveloping_module_bound: OnceLock::new(),
            regular_enveloping_module_fail: OnceLock::new(),
            dual_enveloping_module_weak: RwLock::new(None),
            dual_enveloping_module_bound: OnceLock::new(),
            dual_enveloping_module_fail: OnceLock::new(),
            is_enveloping: OnceLock::new(),
            trivial_extension_level: OnceLock::new(),
            trivial_extension_weak: RwLock::new(None),
            trivial_extension_bound: OnceLock::new(),
            trivial_extension_fail: OnceLock::new(),
        }));
        if quotient.ideal().is_admissible()? {
            let _ = quotient.0.admissible.set(true);
        }
        Ok(quotient)
    }
}

impl<F: Field> fmt::Debug for QuotientPathAlgebra<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("QuotientPathAlgebra")
            .field("parent", self.original_path_algebra())
            .field("relations", &self.relations())
            .finish_non_exhaustive()
    }
}

impl<F: Field> PartialEq for QuotientPathAlgebra<F> {
    fn eq(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}

impl<F: Field> Eq for QuotientPathAlgebra<F> {}

/// QPA's canonical-basis object for a quotient path algebra.
///
/// QPA constructs this object unconditionally and binds its vectors only when
/// the defining ideal already carries a complete basis with finite nontips.
#[derive(Clone, Debug)]
pub struct QuotientCanonicalBasis<F: Field> {
    quotient: QuotientPathAlgebra<F>,
}

impl<F: Field> QuotientCanonicalBasis<F> {
    /// Return QPA's `UnderlyingLeftModule` value.
    #[must_use]
    pub const fn quotient(&self) -> &QuotientPathAlgebra<F> {
        &self.quotient
    }

    /// Return whether QPA bound the `BasisVectors` attribute.
    #[must_use]
    pub fn has_basis_vectors(&self) -> bool {
        matches!(
            self.quotient.0.canonical_basis.get(),
            Some(QuotientCanonicalBasisState::WithVectors { .. })
        )
    }

    /// Return the bound canonical vectors, if available.
    #[must_use]
    pub fn basis_vectors(&self) -> Option<Vec<QuotientPathAlgebraElement<F>>> {
        let Some(QuotientCanonicalBasisState::WithVectors {
            representatives,
            is_normal_form,
        }) = self.quotient.0.canonical_basis.get()
        else {
            return None;
        };
        Some(
            representatives
                .iter()
                .cloned()
                .map(|representative| QuotientPathAlgebraElement {
                    algebra: self.quotient.clone(),
                    representative,
                    is_normal_form: *is_normal_form,
                })
                .collect(),
        )
    }

    /// Consume this basis object into its optional vector list.
    #[must_use]
    pub fn into_basis_vectors(self) -> Option<Vec<QuotientPathAlgebraElement<F>>> {
        self.basis_vectors()
    }

    /// Return whether two handles denote QPA's identical cached basis object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        self.quotient.is_identical(&other.quotient)
    }

    /// Translate QPA's canonical-basis `Coefficients` method.
    ///
    /// # Errors
    ///
    /// Vectors must be bound and the element must be a normal form in the
    /// identical quotient family.
    pub fn coefficients(
        &self,
        element: &QuotientPathAlgebraElement<F>,
    ) -> Result<Vec<F::Element>, QuotientError> {
        if !self.has_basis_vectors() {
            return Err(QuotientError::CanonicalBasisVectorsUnavailable);
        }
        self.quotient().canonical_basis_coefficients(element)
    }
}

/// QPA's nice-basis conversion data for a space of quotient elements.
///
/// The monomial set is built from the stored representatives of the supplied
/// generators using GAP `AddSet` path order. It is intentionally not derived
/// from a separately chosen linear-algebra basis.
#[derive(Clone, Debug)]
pub struct QuotientNiceBasis<F: Field> {
    quotient: QuotientPathAlgebra<F>,
    monomials: Vec<Path>,
}

impl<F: Field> QuotientNiceBasis<F> {
    /// Translate QPA's `NiceFreeLeftModuleInfo` construction.
    ///
    /// # Errors
    ///
    /// Every generator must belong to the supplied identical quotient.
    pub fn new(
        quotient: QuotientPathAlgebra<F>,
        generators: &[QuotientPathAlgebraElement<F>],
    ) -> Result<Self, QuotientError> {
        ensure_quotient_elements(&quotient, generators)?;
        let mut monomials = Vec::new();
        for generator in generators {
            for term in generator.representative().terms() {
                if !monomials.contains(term.path()) {
                    monomials.push(term.path().clone());
                }
            }
        }
        sort_paths_as_qpa_set(&mut monomials)?;
        Ok(Self {
            quotient,
            monomials,
        })
    }

    /// Return the quotient family handled by this conversion data.
    #[must_use]
    pub const fn quotient(&self) -> &QuotientPathAlgebra<F> {
        &self.quotient
    }

    /// Return QPA's sorted monomial support set.
    #[must_use]
    pub fn monomials(&self) -> &[Path] {
        &self.monomials
    }

    /// Translate QPA's `NiceVector`, with `None` representing GAP `fail`.
    ///
    /// # Errors
    ///
    /// The element must belong to the quotient used to construct this data.
    pub fn nice_vector(
        &self,
        element: &QuotientPathAlgebraElement<F>,
    ) -> Result<Option<Vec<F::Element>>, QuotientError> {
        ensure_quotient_elements(self.quotient(), std::slice::from_ref(element))?;
        let field = self.quotient().original_path_algebra().field();
        let mut coordinates = if self.monomials.is_empty() {
            vec![field.zero()]
        } else {
            vec![field.zero(); self.monomials.len()]
        };
        for term in element.representative().terms() {
            let Some(position) = self.monomials.iter().position(|path| path == term.path()) else {
                return Ok(None);
            };
            coordinates[position] = term.coefficient().clone();
        }
        Ok(Some(coordinates))
    }

    /// Translate QPA's `UglyVector`, with `None` representing GAP `fail`.
    ///
    /// QPA passes `normal = true` to its quotient constructor. Consequently a
    /// quotient with a normal-form function trusts the supplied monomials,
    /// while a lazy quotient ignores that flag and leaves the result unmarked.
    ///
    /// # Errors
    ///
    /// Propagates parent sparse-element and quotient construction failures.
    pub fn element_from_nice_vector(
        &self,
        coordinates: &[F::Element],
    ) -> Result<Option<QuotientPathAlgebraElement<F>>, QuotientError> {
        let field = self.quotient().original_path_algebra().field();
        let expected_length = if self.monomials.is_empty() {
            1
        } else {
            self.monomials.len()
        };
        if coordinates.len() != expected_length {
            return Ok(None);
        }
        if self.monomials.is_empty() {
            return Ok(coordinates
                .iter()
                .all(|coefficient| field.is_zero(coefficient))
                .then(|| self.quotient().zero()));
        }
        let representative = self.quotient().original_path_algebra().from_terms(
            self.monomials
                .iter()
                .cloned()
                .zip(coordinates.iter().cloned()),
        )?;
        Ok(Some(
            self.quotient.element_with_normality(representative, true)?,
        ))
    }
}

/// QPA's natural surjection from a path algebra to its quotient by an ideal.
///
/// The pinned QPA method constructs the quotient without forcing a Gröbner
/// basis, maps the path-algebra generators to the corresponding quotient
/// generators, and marks the mapping surjective. Rust stores the mathematically
/// equivalent generator map as the canonical representative projection.
#[derive(Clone)]
pub struct NaturalQuotientMap<F: Field> {
    source: PathAlgebra<F>,
    range: QuotientPathAlgebra<F>,
    kernel: PathAlgebraIdeal<F>,
}

impl<F: Field> NaturalQuotientMap<F> {
    fn new(ideal: &PathAlgebraIdeal<F>) -> Result<Self, QuotientError> {
        Ok(Self {
            source: ideal.algebra().clone(),
            range: QuotientPathAlgebra::by_ideal(ideal)?,
            kernel: ideal.clone(),
        })
    }

    /// Return the source path algebra.
    #[must_use]
    pub const fn source(&self) -> &PathAlgebra<F> {
        &self.source
    }

    /// Return the factor path algebra.
    #[must_use]
    pub const fn range(&self) -> &QuotientPathAlgebra<F> {
        &self.range
    }

    /// Return the defining kernel ideal.
    #[must_use]
    pub const fn kernel(&self) -> &PathAlgebraIdeal<F> {
        &self.kernel
    }

    /// Return QPA's explicitly installed surjectivity property.
    #[must_use]
    pub const fn is_surjective(&self) -> bool {
        true
    }

    /// Map an element by replacing each path-algebra generator with its
    /// quotient generator.
    ///
    /// Since that is precisely the canonical projection, QPA's
    /// `MappedExpression` result has the same stored parent representative.
    /// The range performs normal-form reduction only when its defining ideal
    /// already carries a complete basis.
    ///
    /// # Errors
    ///
    /// Returns an error for a foreign source element or reduction failure.
    pub fn image(
        &self,
        element: &PathAlgebraElement<F>,
    ) -> Result<QuotientPathAlgebraElement<F>, QuotientError> {
        ensure_parent_element(self.source(), element)?;
        self.range().element(element.clone())
    }

    /// Return the images of QPA's vertex-then-arrow generator list.
    ///
    /// # Errors
    ///
    /// Propagates quotient element construction failures.
    pub fn generator_images(&self) -> Result<Vec<QuotientPathAlgebraElement<F>>, QuotientError> {
        self.range().generators()
    }
}

impl<F: Field> fmt::Debug for NaturalQuotientMap<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("NaturalQuotientMap")
            .field("source", self.source())
            .field("range", self.range())
            .field("surjective", &true)
            .finish_non_exhaustive()
    }
}

/// A possibly lazy element of a quotient of a path algebra.
#[derive(Clone, Debug)]
pub struct QuotientPathAlgebraElement<F: Field> {
    algebra: QuotientPathAlgebra<F>,
    representative: PathAlgebraElement<F>,
    is_normal_form: bool,
}

impl<F: Field> QuotientPathAlgebraElement<F> {
    /// Return the containing quotient algebra.
    #[must_use]
    pub const fn algebra(&self) -> &QuotientPathAlgebra<F> {
        &self.algebra
    }

    /// Return QPA's stored underlying path-algebra representative.
    ///
    /// This is a canonical normal form exactly when [`Self::is_normal_form`]
    /// is true. Quotients constructed from ideals without attached Gröbner
    /// bases deliberately retain lazy, unreduced representatives.
    #[must_use]
    pub const fn representative(&self) -> &PathAlgebraElement<F> {
        &self.representative
    }

    /// Return QPA's `IsNormalForm` filter state.
    #[must_use]
    pub const fn is_normal_form(&self) -> bool {
        self.is_normal_form
    }

    /// Return QPA's quotient `ExtRepOfObj`, delegated to the representative.
    #[must_use]
    pub fn external_representation(&self) -> PathAlgebraExternalRepresentation<F::Element> {
        self.representative().external_representation()
    }

    /// Translate QPA's `VertexPosition` on a residue class of a trivial path.
    ///
    /// # Errors
    ///
    /// Reports the documented source precondition explicitly instead of
    /// indexing a malformed representative unchecked.
    pub fn vertex_position(&self) -> Result<usize, QuotientError> {
        Ok(self.representative().vertex_position()?)
    }

    /// Normalize this representative using the quotient's installed basis.
    ///
    /// # Errors
    ///
    /// Returns [`QuotientError::MissingGroebnerBasis`] when QPA did not install
    /// a normal-form function on this quotient family.
    pub fn normalized(&self) -> Result<Self, QuotientError> {
        if self.is_normal_form() {
            return Ok(self.clone());
        }
        if !self.algebra().has_normal_form_function() {
            return Err(QuotientError::MissingGroebnerBasis);
        }
        self.algebra.element(self.representative().clone())
    }

    /// Return whether this is quotient zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.representative.is_zero()
    }

    /// Add and normalize quotient elements.
    ///
    /// # Errors
    ///
    /// Returns an error for different quotient families or reduction failure.
    pub fn add(&self, right: &Self) -> Result<Self, QuotientError> {
        self.ensure_same_quotient(right)?;
        self.algebra.element_with_normality(
            self.representative.add(&right.representative)?,
            self.is_normal_form && right.is_normal_form,
        )
    }

    /// Return the additive inverse in normal form.
    ///
    /// # Errors
    ///
    /// Propagates arithmetic and reduction failures.
    pub fn neg(&self) -> Result<Self, QuotientError> {
        self.algebra
            .element_with_normality(self.representative.neg()?, self.is_normal_form)
    }

    /// Subtract and normalize quotient elements.
    ///
    /// # Errors
    ///
    /// Returns an error for different quotient families or reduction failure.
    pub fn sub(&self, right: &Self) -> Result<Self, QuotientError> {
        self.ensure_same_quotient(right)?;
        self.algebra.element_with_normality(
            self.representative.sub(&right.representative)?,
            self.is_normal_form && right.is_normal_form,
        )
    }

    /// Multiply and normalize quotient elements.
    ///
    /// # Errors
    ///
    /// Returns an error for different quotient families or reduction failure.
    pub fn multiply(&self, right: &Self) -> Result<Self, QuotientError> {
        self.ensure_same_quotient(right)?;
        self.algebra
            .element_with_normality(self.representative.multiply(&right.representative)?, false)
    }

    /// Scale on the left and normalize.
    ///
    /// # Errors
    ///
    /// Propagates arithmetic and reduction failures.
    pub fn scale_left(&self, scalar: &F::Element) -> Result<Self, QuotientError> {
        self.algebra
            .element_with_normality(self.representative.scale_left(scalar)?, self.is_normal_form)
    }

    /// Raise to a nonnegative power and normalize each multiplication.
    ///
    /// # Errors
    ///
    /// Propagates quotient arithmetic failures.
    pub fn pow(&self, mut exponent: u64) -> Result<Self, QuotientError> {
        let mut result = self.algebra.one()?;
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

    /// Return QPA's leading term of the normal representative.
    #[must_use]
    pub fn leading_term(&self) -> Self {
        Self {
            algebra: self.algebra.clone(),
            representative: self.representative.leading_term(),
            is_normal_form: self.is_normal_form,
        }
    }

    /// Return the leading coefficient of the normal representative.
    #[must_use]
    pub fn leading_coefficient(&self) -> F::Element {
        self.representative.leading_coefficient()
    }

    /// Return the leading monomial of the normal representative.
    #[must_use]
    pub fn leading_monomial(&self) -> Path {
        self.representative.leading_monomial()
    }

    /// Translate QPA's left-uniform test on the representative.
    #[must_use]
    pub fn is_left_uniform(&self) -> bool {
        self.representative.is_left_uniform()
    }

    /// Translate QPA's right-uniform test on the representative.
    #[must_use]
    pub fn is_right_uniform(&self) -> bool {
        self.representative.is_right_uniform()
    }

    /// Translate QPA's uniform test on the representative.
    #[must_use]
    pub fn is_uniform(&self) -> bool {
        self.representative.is_uniform()
    }

    /// Compare normal representatives using QPA's path-algebra comparison.
    ///
    /// # Errors
    ///
    /// Returns an error for different quotient families.
    pub fn compare(&self, right: &Self) -> Result<Ordering, QuotientError> {
        self.ensure_same_quotient(right)?;
        if !self.is_normal_form() || !right.is_normal_form() {
            return Err(QuotientError::ElementsNotInNormalForm);
        }
        Ok(self.representative.compare(&right.representative)?)
    }

    /// Translate `OppositePathAlgebraElement` for a quotient element.
    ///
    /// QPA maps the stored representative termwise and then constructs an
    /// element of the opposite quotient, whose freshly computed basis performs
    /// normal-form reduction.
    ///
    /// # Errors
    ///
    /// Propagates opposite quotient, path, and reduction failures.
    pub fn opposite_high_level(&self) -> Result<Self, QuotientError> {
        let opposite = self.algebra().opposite_high_level()?;
        let representative = self.representative().opposite()?;
        opposite.element(representative)
    }

    fn ensure_same_quotient(&self, right: &Self) -> Result<(), QuotientError> {
        if self.algebra.is_identical(&right.algebra) {
            Ok(())
        } else {
            Err(QuotientError::DifferentQuotients)
        }
    }
}

impl<F: Field> PartialEq for QuotientPathAlgebraElement<F> {
    fn eq(&self, other: &Self) -> bool {
        self.algebra == other.algebra && self.representative == other.representative
    }
}

impl<F: Field> Eq for QuotientPathAlgebraElement<F> {}

/// QPA quotient-algebra dimension, including `infinity`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AlgebraDimension {
    /// An exact finite vector-space dimension.
    Finite(Integer),
    /// The nontip automaton contains a cycle.
    Infinite,
}

/// Translate QPA's `IsHomogeneousListOfElements(gens, degree)`.
///
/// Rust's element type guarantees the input consists of path-algebra
/// elements; an empty list is homogeneous, matching QPA's printed branch.
///
/// # Errors
///
/// QPA dispatch requires a positive degree. Nonempty inputs must belong to one
/// path-algebra family.
pub fn is_homogeneous_list_of_elements<F: Field>(
    elements: &[PathAlgebraElement<F>],
    degree: usize,
) -> Result<bool, QuotientError> {
    if degree == 0 {
        return Err(QuotientError::NonPositiveDegree);
    }
    let Some(first) = elements.first() else {
        return Ok(true);
    };
    ensure_parent_elements(first.algebra(), elements)?;
    Ok(elements.iter().all(|element| {
        element.is_zero()
            || element
                .terms()
                .iter()
                .all(|term| term.path().length() == Some(degree))
    }))
}

fn quotient_ideal_span<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
    generators: &[QuotientPathAlgebraElement<F>],
    algebra_basis: &[QuotientPathAlgebraElement<F>],
) -> Result<Vec<QuotientPathAlgebraElement<F>>, QuotientError> {
    let mut products = generators.to_vec();
    for left in algebra_basis {
        for generator in generators {
            let left_product = left.multiply(generator)?;
            for right in algebra_basis {
                let product = left_product.multiply(right)?;
                if !product.is_zero() {
                    products.push(product);
                }
            }
        }
    }
    quotient_span_basis(quotient, &products, algebra_basis)
}

fn quotient_product_space<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
    left: &[QuotientPathAlgebraElement<F>],
    right: &[QuotientPathAlgebraElement<F>],
    algebra_basis: &[QuotientPathAlgebraElement<F>],
) -> Result<Vec<QuotientPathAlgebraElement<F>>, QuotientError> {
    let mut products = Vec::new();
    for left in left {
        for right in right {
            let product = left.multiply(right)?;
            if !product.is_zero() {
                products.push(product);
            }
        }
    }
    quotient_span_basis(quotient, &products, algebra_basis)
}

fn quotient_span_basis<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
    elements: &[QuotientPathAlgebraElement<F>],
    algebra_basis: &[QuotientPathAlgebraElement<F>],
) -> Result<Vec<QuotientPathAlgebraElement<F>>, QuotientError> {
    let basis_paths = algebra_basis
        .iter()
        .map(|element| canonical_basis_path(element.representative()))
        .collect::<Result<Vec<_>, _>>()?;
    let rows = elements
        .iter()
        .map(|element| element_coordinates(element.representative(), &basis_paths))
        .collect::<Result<Vec<_>, _>>()?;
    let matrix = DenseMatrix::from_rows(rows, quotient.original_path_algebra().field().clone())?;
    let vectors = matrix.semi_echelon_transformation().vectors;
    let mut basis = Vec::with_capacity(vectors.nrows());
    for row in 0..vectors.nrows() {
        let element = quotient.original_path_algebra().from_terms(
            basis_paths
                .iter()
                .cloned()
                .zip(vectors.row(row).iter().cloned()),
        )?;
        basis.push(quotient.element(element)?);
    }
    Ok(basis)
}

fn path_algebra_ideal_span<F: Field>(
    algebra: &PathAlgebra<F>,
    generators: &[PathAlgebraElement<F>],
    algebra_basis: &[PathAlgebraElement<F>],
) -> Result<Vec<PathAlgebraElement<F>>, QuotientError> {
    let mut products = generators.to_vec();
    for left in algebra_basis {
        for generator in generators {
            let left_product = left.multiply(generator)?;
            for right in algebra_basis {
                let product = left_product.multiply(right)?;
                if !product.is_zero() {
                    products.push(product);
                }
            }
        }
    }
    let basis_paths = algebra_basis
        .iter()
        .map(canonical_basis_path)
        .collect::<Result<Vec<_>, _>>()?;
    let rows = products
        .iter()
        .map(|element| element_coordinates(element, &basis_paths))
        .collect::<Result<Vec<_>, _>>()?;
    let matrix = DenseMatrix::from_rows(rows, algebra.field().clone())?;
    let vectors = matrix.semi_echelon_transformation().vectors;
    (0..vectors.nrows())
        .map(|row| {
            Ok(algebra.from_terms(
                basis_paths
                    .iter()
                    .cloned()
                    .zip(vectors.row(row).iter().cloned()),
            )?)
        })
        .collect()
}

fn canonical_basis_path<F: Field>(element: &PathAlgebraElement<F>) -> Result<Path, QuotientError> {
    if element.terms().len() != 1 {
        return Err(QuotientError::InvalidCanonicalBasis);
    }
    Ok(element.terms()[0].path().clone())
}

fn element_coordinates<F: Field>(
    element: &PathAlgebraElement<F>,
    basis_paths: &[Path],
) -> Result<Vec<F::Element>, QuotientError> {
    let field = element.algebra().field();
    let mut coordinates = vec![field.zero(); basis_paths.len()];
    for term in element.terms() {
        let position = basis_paths
            .iter()
            .position(|path| path == term.path())
            .ok_or(QuotientError::ElementOutsideCanonicalBasis)?;
        coordinates[position] = field.add(&coordinates[position], term.coefficient());
    }
    Ok(coordinates)
}

fn sort_paths_as_qpa_set(paths: &mut [Path]) -> Result<(), QuotientError> {
    for index in 1..paths.len() {
        let mut position = index;
        while position > 0 {
            let ordering = paths[position - 1]
                .quiver()
                .compare_paths(&paths[position - 1], &paths[position])
                .map_err(PathAlgebraError::from)?;
            if !ordering.is_gt() {
                break;
            }
            paths.swap(position - 1, position);
            position -= 1;
        }
    }
    Ok(())
}

fn ensure_parent_element<F: Field>(
    parent: &PathAlgebra<F>,
    element: &PathAlgebraElement<F>,
) -> Result<(), QuotientError> {
    if parent.is_identical(element.algebra()) {
        Ok(())
    } else {
        Err(QuotientError::DifferentParentAlgebras)
    }
}

fn ensure_parent_elements<F: Field>(
    parent: &PathAlgebra<F>,
    elements: &[PathAlgebraElement<F>],
) -> Result<(), QuotientError> {
    for element in elements {
        ensure_parent_element(parent, element)?;
    }
    Ok(())
}

fn ensure_quotient_elements<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
    elements: &[QuotientPathAlgebraElement<F>],
) -> Result<(), QuotientError> {
    if elements
        .iter()
        .all(|element| quotient.is_identical(element.algebra()))
    {
        Ok(())
    } else {
        Err(QuotientError::DifferentQuotients)
    }
}

/// Failure to construct or operate on a QPA ideal or quotient path algebra.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum QuotientError {
    /// An ideal generator or representative belongs to another parent algebra.
    #[error("an element belongs to a different parent path algebra")]
    DifferentParentAlgebras,
    /// Quotient arithmetic requires identical quotient families.
    #[error("the elements belong to different quotient path algebras")]
    DifferentQuotients,
    /// QPA only installs quotient comparison for normal-form elements.
    #[error("quotient elements must be in normal form for this operation")]
    ElementsNotInNormalForm,
    /// An ideal-level method requiring an attached basis did not receive one.
    #[error("the path-algebra ideal has no attached Groebner basis")]
    MissingGroebnerBasis,
    /// `ProductOfIdeals` requires finite nontips for its right argument.
    #[error("the second ideal does not admit finitely many nontips")]
    InfiniteRightFactorNontips,
    /// `MinimalGeneratingSetOfIdeal` dispatches only for admissible ideals.
    #[error("the path-algebra ideal is not admissible")]
    IdealNotAdmissible,
    /// A required canonical basis of the parent path algebra is infinite.
    #[error("the parent path algebra has no finite canonical basis")]
    InfiniteParentBasis,
    /// A required quotient basis has infinitely many nontips.
    #[error("the quotient path algebra has no finite canonical basis")]
    InfiniteQuotientBasis,
    /// A canonical nontip basis vector was not represented by one path.
    #[error("the canonical quotient basis has an invalid vector")]
    InvalidCanonicalBasis,
    /// A normal form used a path absent from the canonical nontip basis.
    #[error("an element lies outside the canonical nontip basis")]
    ElementOutsideCanonicalBasis,
    /// QPA's canonical basis object has no bound vectors in this state.
    #[error("the canonical quotient basis has no available basis vectors")]
    CanonicalBasisVectorsUnavailable,
    /// QPA declares homogeneous-list checking only for positive degrees.
    #[error("the homogeneous degree must be positive")]
    NonPositiveDegree,
    /// Concurrent opposite-quotient construction lost its reciprocal cache.
    #[error("concurrent opposite-quotient cache initialization failed")]
    OppositeCacheRace,
    /// QPA Groebner computation or reduction failed.
    #[error(transparent)]
    Groebner(#[from] GroebnerError),
    /// Parent path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// A source-translated quiver lookup failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// Exact GAP-compatible row reduction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, Rational, RationalField},
    };
    use serde_json::{Value, json};

    use crate::{
        ArrowIndex, ArrowSpec, PathAlgebra, PathAlgebraElement, PathAlgebraExternalRepresentation,
        Quiver, VertexIndex, opposite_relations,
    };

    use super::{
        AlgebraDimension, PathAlgebraIdeal, QuotientError, QuotientNiceBasis, QuotientPathAlgebra,
        is_homogeneous_list_of_elements,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn arrow(index: usize) -> ArrowIndex {
        ArrowIndex::new(index).unwrap()
    }

    fn loop_algebra() -> PathAlgebra<numerica::domains::rational::RationalField> {
        let quiver = Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(1), "y"),
            ],
        )
        .unwrap();
        PathAlgebra::new(Q, quiver).unwrap()
    }

    fn loop_quiver() -> Quiver {
        Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(1), "y"),
            ],
        )
        .unwrap()
    }

    fn foundational_quiver() -> Quiver {
        Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "a"),
                ArrowSpec::named(vertex(1), vertex(2), "b"),
                ArrowSpec::named(vertex(2), vertex(3), "c"),
                ArrowSpec::named(vertex(3), vertex(1), "d"),
                ArrowSpec::named(vertex(2), vertex(2), "e"),
            ],
        )
        .unwrap()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the quotient fixture uses only small integers")
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
                        "coefficient": coefficient_json(term.coefficient())
                    })
                })
                .collect(),
        )
    }

    fn rational_json(element: &PathAlgebraElement<RationalField>) -> Value {
        element_json(element, |coefficient| {
            json!([
                small_integer(coefficient.numerator_ref()),
                small_integer(coefficient.denominator_ref())
            ])
        })
    }

    fn zp_json(element: &PathAlgebraElement<Zp>, field: &Zp) -> Value {
        element_json(element, |coefficient| {
            small_integer(&field.to_integer(coefficient)).into()
        })
    }

    fn z2_json(element: &PathAlgebraElement<Z2>) -> Value {
        element_json(element, |coefficient| {
            small_integer(&Z2.to_integer(coefficient)).into()
        })
    }

    fn rational_external_json(element: &PathAlgebraElement<RationalField>) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().external_representation(),
                        "coefficient": [
                            small_integer(term.coefficient().numerator_ref()),
                            small_integer(term.coefficient().denominator_ref())
                        ]
                    })
                })
                .collect(),
        )
    }

    fn z2_external_json(element: &PathAlgebraElement<Z2>) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().external_representation(),
                        "coefficient": small_integer(&Z2.to_integer(term.coefficient()))
                    })
                })
                .collect(),
        )
    }

    fn monomial_external_paths<F: Field>(elements: &[PathAlgebraElement<F>]) -> Value {
        Value::Array(
            elements
                .iter()
                .map(|element| {
                    assert_eq!(element.terms().len(), 1);
                    json!(element.terms()[0].path().external_representation())
                })
                .collect(),
        )
    }

    fn rational_coefficients_json(coefficients: &[Rational]) -> Value {
        Value::Array(
            coefficients
                .iter()
                .map(|coefficient| {
                    json!([
                        small_integer(coefficient.numerator_ref()),
                        small_integer(coefficient.denominator_ref())
                    ])
                })
                .collect(),
        )
    }

    fn z2_coefficients_json(coefficients: &[u8]) -> Value {
        Value::Array(
            coefficients
                .iter()
                .map(|coefficient| small_integer(&Z2.to_integer(coefficient)).into())
                .collect(),
        )
    }

    #[test]
    fn ideal_membership_uses_qpa_complete_reduction() {
        let algebra = loop_algebra();
        let x = algebra
            .element(algebra.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let ideal = PathAlgebraIdeal::new(algebra.clone(), vec![x.pow(2).unwrap()]).unwrap();
        assert!(ideal.contains(x.pow(3).unwrap()).unwrap());
        assert!(!ideal.contains(x).unwrap());
        assert!(ideal.generators_are_monomial());
    }

    #[test]
    fn ideal_level_right_basis_requires_and_uses_an_attached_complete_basis() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let algebra = PathAlgebra::new(Q, quiver).unwrap();
        let x = algebra
            .element(algebra.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let cubic = x.pow(3).unwrap();
        let ideal = PathAlgebraIdeal::new(algebra, vec![cubic.clone()]).unwrap();
        assert_eq!(
            ideal.right_groebner_basis().unwrap_err(),
            QuotientError::MissingGroebnerBasis
        );
        let attached = ideal
            .with_complete_groebner_basis(std::slice::from_ref(&cubic))
            .unwrap();
        let right = attached.right_groebner_basis().unwrap();
        assert_eq!(right.relations(), [cubic]);
    }

    #[test]
    fn quotient_arithmetic_is_normalized_after_each_operation() {
        let algebra = loop_algebra();
        let x = algebra
            .element(algebra.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y = algebra
            .element(algebra.quiver().arrow(arrow(2)).unwrap().path())
            .unwrap();
        let quotient = QuotientPathAlgebra::by_relations_high_level(
            algebra,
            &[
                x.pow(2).unwrap(),
                x.multiply(&y)
                    .unwrap()
                    .sub(&y.multiply(&x).unwrap())
                    .unwrap(),
                y.pow(2).unwrap(),
            ],
        )
        .unwrap();
        let generators = quotient.generators().unwrap();
        let x = &generators[1];
        let y = &generators[2];
        assert!(x.pow(2).unwrap().is_zero());
        assert_eq!(y.multiply(x).unwrap(), x.multiply(y).unwrap());
        assert!(y.multiply(x).unwrap().multiply(y).unwrap().is_zero());
        assert_eq!(quotient.relations().len(), 3);
    }

    #[test]
    #[allow(
        clippy::many_single_char_names,
        clippy::similar_names,
        clippy::too_many_lines
    )]
    fn pinned_qpa_foundational_pathalg_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_foundational.schema-v1.json"
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

        let quiver = foundational_quiver();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let arrows = (1..=5)
            .map(|index| {
                algebra
                    .element(quiver.arrow(arrow(index)).unwrap().path())
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let [a, b, c, d, _e] = arrows.as_slice() else {
            unreachable!();
        };

        let a_plus_b = a.add(b).unwrap();
        let uniform = algebra
            .make_uniform_on_right(&[a_plus_b.clone(), a.clone()])
            .unwrap();
        assert_eq!(
            Value::Array(uniform.iter().map(rational_external_json).collect()),
            expected["right_uniform"]
        );
        let extended = algebra
            .generators_times_arrows_on_right(std::slice::from_ref(&a_plus_b))
            .unwrap();
        assert_eq!(monomial_external_paths(&extended), expected["times_arrows"]);
        for exponent in [0_i64, 1, 2, 3] {
            assert_eq!(
                monomial_external_paths(&algebra.nth_power_of_arrow_ideal(exponent).unwrap()),
                expected["arrow_powers"][format!("n{exponent}")]
            );
        }
        let mut relations = vec![a_plus_b.clone()];
        algebra
            .add_nth_power_to_relations(&mut relations, 2)
            .unwrap();
        assert_eq!(
            json!({
                "same_object": true,
                "relations": relations.iter().map(rational_external_json).collect::<Vec<_>>()
            }),
            expected["add_power"]
        );

        let a_squared = a.pow(2).unwrap();
        let ideal = PathAlgebraIdeal::new(algebra.clone(), vec![a_squared.clone()]).unwrap();
        let projection = ideal.natural_homomorphism().unwrap();
        let lazy_image = projection.image(&a.add(&a_squared).unwrap()).unwrap();
        let lazy_sum = projection.image(&a_plus_b).unwrap();
        let lazy_uniform = projection
            .range()
            .make_uniform_on_right(std::slice::from_ref(&lazy_sum))
            .unwrap();
        let lazy_extended = projection
            .range()
            .generators_times_arrows_on_right(std::slice::from_ref(&lazy_sum))
            .unwrap();
        assert!(projection.source().is_identical(&algebra));
        assert!(
            projection
                .range()
                .original_path_algebra()
                .is_identical(&algebra)
        );
        assert!(projection.kernel().is_identical(&ideal));
        assert!(projection.is_surjective());
        assert!(!lazy_image.is_normal_form());
        assert_eq!(
            lazy_image.normalized().unwrap_err(),
            QuotientError::MissingGroebnerBasis
        );

        let attached_ideal = ideal
            .with_complete_groebner_basis(std::slice::from_ref(&a_squared))
            .unwrap();
        assert!(attached_ideal.is_identical(&ideal));
        assert_eq!(
            attached_ideal
                .groebner_basis()
                .unwrap()
                .tip_reduced_property(),
            None
        );
        let attached_projection = attached_ideal.natural_homomorphism().unwrap();
        assert_eq!(
            attached_ideal
                .groebner_basis()
                .unwrap()
                .tip_reduced_property(),
            Some(true)
        );
        let attached_image = attached_projection
            .image(&a.add(&a_squared).unwrap())
            .unwrap();
        assert_eq!(
            json!({
                "source_identical": projection.source().is_identical(&algebra),
                "range_original_identical": projection.range().original_path_algebra().is_identical(&algebra),
                "surjective": projection.is_surjective(),
                "image": rational_external_json(lazy_image.representative()),
                "image_is_normal_form": lazy_image.is_normal_form(),
                "right_uniform": lazy_uniform.iter().map(|element| rational_external_json(element.representative())).collect::<Vec<_>>(),
                "times_arrows": lazy_extended.iter().map(|element| rational_external_json(element.representative())).collect::<Vec<_>>(),
                "attached_image": rational_external_json(attached_image.representative()),
                "attached_image_is_normal_form": attached_image.is_normal_form()
            }),
            expected["natural_map"]
        );

        let opposite_parent = algebra.opposite().unwrap();
        assert!(opposite_parent.opposite().unwrap().is_identical(&algebra));
        assert!(
            opposite_parent
                .quiver()
                .opposite()
                .unwrap()
                .is_identical(&quiver)
        );
        let q2 = Q.to_element(2.into(), 1.into(), true);
        let q3 = Q.to_element(3.into(), 1.into(), true);
        let element = a
            .multiply(b)
            .unwrap()
            .scale_left(&q2)
            .unwrap()
            .add(&b.multiply(c).unwrap().scale_left(&q3).unwrap())
            .unwrap();
        let opposite_element = element.opposite().unwrap();
        let double_element = opposite_element.opposite().unwrap();
        let source_relations = vec![
            a_squared.add(&b.multiply(c).unwrap()).unwrap(),
            c.multiply(d).unwrap(),
        ];
        let opposite_source_relations = opposite_relations(&source_relations).unwrap();
        assert_eq!(
            json!({
                "double_identical": opposite_parent.opposite().unwrap().is_identical(&algebra),
                "quiver_double_identical": opposite_parent.quiver().opposite().unwrap().is_identical(&quiver),
                "element": rational_external_json(&opposite_element),
                "double_element": rational_external_json(&double_element),
                "relations": opposite_source_relations.iter().map(rational_external_json).collect::<Vec<_>>()
            }),
            expected["opposite_parent"]
        );

        let truncated = QuotientPathAlgebra::truncated_high_level(Q, quiver.clone(), 3).unwrap();
        let truncated_opposite = truncated.opposite_high_level().unwrap();
        let truncated_double = truncated_opposite.opposite_high_level().unwrap();
        let qa = truncated
            .path(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let qb = truncated
            .path(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let qc = truncated
            .path(quiver.arrow(arrow(3)).unwrap().path())
            .unwrap();
        let quotient_element = qa
            .add(&qa.pow(2).unwrap())
            .unwrap()
            .add(&qb.multiply(&qc).unwrap())
            .unwrap();
        let opposite_quotient_element = quotient_element.opposite_high_level().unwrap();
        let AlgebraDimension::Finite(dimension) = truncated.dimension().unwrap() else {
            panic!("the truncated algebra must be finite-dimensional");
        };
        let AlgebraDimension::Finite(opposite_dimension) = truncated_opposite.dimension().unwrap()
        else {
            panic!("the opposite truncated algebra must be finite-dimensional");
        };
        assert_eq!(
            json!({
                "dimension": small_integer(&dimension),
                "opposite_dimension": small_integer(&opposite_dimension),
                "double_identical": truncated_double.is_identical(&truncated),
                "opposite_element": rational_external_json(opposite_quotient_element.representative())
            }),
            expected["truncated"]
        );

        let algebra2 = PathAlgebra::new(Z2, quiver.clone()).unwrap();
        let a2 = algebra2
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let b2 = algebra2
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let c2 = algebra2
            .element(quiver.arrow(arrow(3)).unwrap().path())
            .unwrap();
        let a2_squared = a2.pow(2).unwrap();
        let ideal2 = PathAlgebraIdeal::new(algebra2, vec![a2_squared.clone()]).unwrap();
        let projection2 = ideal2.natural_homomorphism().unwrap();
        let image2 = projection2.image(&a2.add(&a2_squared).unwrap()).unwrap();
        let attached2 = ideal2
            .with_complete_groebner_basis(std::slice::from_ref(&a2_squared))
            .unwrap();
        let attached_image2 = attached2
            .natural_homomorphism()
            .unwrap()
            .image(&a2.add(&a2_squared).unwrap())
            .unwrap();
        let opposite_element2 = a2
            .multiply(&b2)
            .unwrap()
            .add(&b2.multiply(&c2).unwrap())
            .unwrap()
            .opposite()
            .unwrap();
        assert_eq!(
            json!({
                "lazy_image": z2_external_json(image2.representative()),
                "lazy_is_normal_form": image2.is_normal_form(),
                "attached_image": z2_external_json(attached_image2.representative()),
                "attached_is_normal_form": attached_image2.is_normal_form(),
                "opposite_element": z2_external_json(&opposite_element2)
            }),
            expected["f2"]
        );
    }

    #[test]
    #[allow(clippy::many_single_char_names, clippy::too_many_lines)]
    fn pinned_qpa_pathalg_conversion_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_conversion.schema-v1.json"
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

        let quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
                ArrowSpec::named(vertex(1), vertex(3), "c"),
            ],
        )
        .unwrap();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let a = algebra
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let b = algebra
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let c = algebra
            .element(quiver.arrow(arrow(3)).unwrap().path())
            .unwrap();
        let relation = a.multiply(&b).unwrap().sub(&c).unwrap();
        let quotient = QuotientPathAlgebra::by_relations_high_level(
            algebra.clone(),
            std::slice::from_ref(&relation),
        )
        .unwrap();
        let canonical_basis = quotient.canonical_basis_object().unwrap();
        let basis = canonical_basis.basis_vectors().unwrap();
        let generators = quotient.generators().unwrap();
        let q2 = Q.to_element(2.into(), 1.into(), true);
        let q3 = Q.to_element(3.into(), 1.into(), true);
        let q5 = Q.to_element(5.into(), 1.into(), true);
        let element = generators[0]
            .scale_left(&q2)
            .unwrap()
            .add(&generators[3].scale_left(&q3).unwrap())
            .unwrap()
            .add(&generators[5].scale_left(&q5).unwrap())
            .unwrap();
        let coefficients = canonical_basis.coefficients(&element).unwrap();

        let raw_representative = a
            .scale_left(&q2)
            .unwrap()
            .add(&a.multiply(&b).unwrap().scale_left(&q3).unwrap())
            .unwrap();
        let descriptor = raw_representative.external_representation();
        assert_eq!(
            algebra
                .element_from_external_representation(&descriptor)
                .unwrap(),
            raw_representative
        );
        let rebuilt = quotient
            .element_from_external_representation(&descriptor)
            .unwrap();
        let nice_basis = QuotientNiceBasis::new(
            quotient.clone(),
            &[
                generators[0].clone(),
                generators[3].clone(),
                generators[5].clone(),
            ],
        )
        .unwrap();
        let nice = nice_basis.nice_vector(&element).unwrap().unwrap();
        let ugly = nice_basis.element_from_nice_vector(&nice).unwrap().unwrap();
        assert_eq!(ugly, element);

        assert_eq!(
            json!({
                "parent_relations": algebra
                    .relations()
                    .iter()
                    .map(rational_external_json)
                    .collect::<Vec<_>>(),
                "quotient_relations": quotient
                    .relations()
                    .iter()
                    .map(rational_external_json)
                    .collect::<Vec<_>>(),
                "ideal_identical": quotient.ideal().is_identical(quotient.ideal()),
                "finite": quotient.is_finite_dimensional().unwrap(),
                "dimension": match quotient.dimension().unwrap() {
                    AlgebraDimension::Finite(value) => small_integer(&value),
                    AlgebraDimension::Infinite => panic!("oracle quotient must be finite"),
                },
                "basis": basis
                    .iter()
                    .map(|value| rational_external_json(value.representative()))
                    .collect::<Vec<_>>(),
                "element": rational_external_json(element.representative()),
                "coefficients": rational_coefficients_json(&coefficients),
                "rebuilt": rational_external_json(rebuilt.representative()),
                "rebuilt_is_normal_form": rebuilt.is_normal_form(),
                "nice": rational_coefficients_json(&nice),
                "ugly": rational_external_json(ugly.representative()),
                "vertex_positions": [
                    algebra
                        .element(quiver.vertex(vertex(2)).unwrap().path())
                        .unwrap()
                        .vertex_position()
                        .unwrap(),
                    generators[1].vertex_position().unwrap()
                ]
            }),
            expected["rational"]
        );

        let malformed = PathAlgebraExternalRepresentation::new(
            Q.zero(),
            descriptor.terms().iter().rev().cloned().collect(),
        );
        assert_eq!(
            algebra
                .element_from_external_representation(&malformed)
                .unwrap_err(),
            crate::PathAlgebraError::ExternalRepresentationNotNormalized
        );

        let algebra2 = PathAlgebra::new(Z2, quiver.clone()).unwrap();
        let a2 = algebra2
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let b2 = algebra2
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let c2 = algebra2
            .element(quiver.arrow(arrow(3)).unwrap().path())
            .unwrap();
        let relation2 = a2.multiply(&b2).unwrap().sub(&c2).unwrap();
        let quotient2 = QuotientPathAlgebra::by_relations_high_level(
            algebra2,
            std::slice::from_ref(&relation2),
        )
        .unwrap();
        let basis2 = quotient2.canonical_basis().unwrap().unwrap();
        let generators2 = quotient2.generators().unwrap();
        let element2 = generators2[0]
            .add(&generators2[3])
            .unwrap()
            .add(&generators2[5])
            .unwrap();
        let coefficients2 = quotient2.canonical_basis_coefficients(&element2).unwrap();
        let descriptor2 = a2
            .add(&a2.multiply(&b2).unwrap())
            .unwrap()
            .external_representation();
        let rebuilt2 = quotient2
            .element_from_external_representation(&descriptor2)
            .unwrap();
        assert_eq!(
            json!({
                "basis": basis2
                    .iter()
                    .map(|value| z2_external_json(value.representative()))
                    .collect::<Vec<_>>(),
                "coefficients": z2_coefficients_json(&coefficients2),
                "rebuilt": z2_external_json(rebuilt2.representative()),
                "rebuilt_is_normal_form": rebuilt2.is_normal_form()
            }),
            expected["f2"]
        );
    }

    #[test]
    #[allow(clippy::many_single_char_names, clippy::too_many_lines)]
    fn pinned_qpa_pathalg_ideal_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_ideal.schema-v1.json"
        ))
        .unwrap();

        let quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(1), vertex(2), "b"),
                ArrowSpec::named(vertex(2), vertex(2), "c"),
                ArrowSpec::named(vertex(2), vertex(3), "d"),
                ArrowSpec::named(vertex(3), vertex(1), "e"),
            ],
        )
        .unwrap();
        let algebra = PathAlgebra::new(Q, quiver).unwrap();
        let arrows = (1..=5)
            .map(|index| {
                algebra
                    .element(algebra.quiver().arrow(arrow(index)).unwrap().path())
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let [a, b, c, d, e] = arrows.as_slice() else {
            unreachable!();
        };
        let relations = vec![
            c.pow(2).unwrap(),
            d.multiply(e).unwrap(),
            e.multiply(a).unwrap(),
            a.multiply(c)
                .unwrap()
                .multiply(d)
                .unwrap()
                .sub(&b.multiply(d).unwrap())
                .unwrap(),
            e.multiply(b).unwrap().multiply(d).unwrap(),
        ];
        let ideal = PathAlgebraIdeal::new(algebra.clone(), relations.clone()).unwrap();
        let quotient = QuotientPathAlgebra::by_relations_high_level(algebra, &relations).unwrap();
        assert_eq!(
            json!({
                "admissible": ideal.is_admissible().unwrap(),
                "monomial": ideal.is_monomial().unwrap(),
                "quotient_admissible": quotient.is_admissible().unwrap(),
                "minimal_generators": ideal
                    .minimal_generating_set()
                    .unwrap()
                    .iter()
                    .map(rational_external_json)
                    .collect::<Vec<_>>(),
            }),
            expected["regression"]
        );

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_algebra = PathAlgebra::new(Q, loop_quiver.clone()).unwrap();
        let x = loop_algebra
            .element(loop_quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let trivial_loop = PathAlgebraIdeal::new(loop_algebra.clone(), Vec::new()).unwrap();
        let linear = PathAlgebraIdeal::new(loop_algebra.clone(), vec![x.clone()]).unwrap();
        let monomial =
            PathAlgebraIdeal::new(loop_algebra.clone(), vec![x.pow(2).unwrap()]).unwrap();
        let derived_monomial = PathAlgebraIdeal::new(
            loop_algebra.clone(),
            vec![
                x.pow(2).unwrap().add(&x.pow(3).unwrap()).unwrap(),
                x.pow(3).unwrap(),
            ],
        )
        .unwrap();

        let acyclic_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "u")]).unwrap();
        let acyclic_algebra = PathAlgebra::new(Q, acyclic_quiver).unwrap();
        let trivial_acyclic = PathAlgebraIdeal::new(acyclic_algebra, Vec::new()).unwrap();
        assert_eq!(
            json!({
                "trivial_loop_admissible": trivial_loop.is_admissible().unwrap(),
                "trivial_acyclic_admissible": trivial_acyclic.is_admissible().unwrap(),
                "linear_admissible": linear.is_admissible().unwrap(),
                "monomial_is_monomial": monomial.is_monomial().unwrap(),
                "derived_monomial": derived_monomial.is_monomial().unwrap(),
            }),
            expected["boundary"]
        );
        assert_eq!(
            json!({
                "empty_degree_2": is_homogeneous_list_of_elements::<RationalField>(&[], 2).unwrap(),
                "degree_2": is_homogeneous_list_of_elements(
                    &[
                        x.pow(2).unwrap(),
                        x.pow(2)
                            .unwrap()
                            .scale_left(&Q.to_element(2.into(), 1.into(), true))
                            .unwrap(),
                        loop_algebra.zero()
                    ],
                    2
                ).unwrap(),
                "mixed_degree_2": is_homogeneous_list_of_elements(
                    &[x.pow(2).unwrap().add(&x.pow(3).unwrap()).unwrap()],
                    2
                ).unwrap(),
            }),
            expected["homogeneous"]
        );

        let right_factor = PathAlgebraIdeal::new(loop_algebra.clone(), vec![x.pow(2).unwrap()])
            .unwrap()
            .with_complete_groebner_basis(&[x.pow(2).unwrap()])
            .unwrap();
        let left_factor = PathAlgebraIdeal::new(loop_algebra.clone(), vec![x.clone()]).unwrap();
        let product = left_factor.product(&right_factor).unwrap();
        assert_eq!(
            Value::Array(
                product
                    .generators()
                    .iter()
                    .map(rational_external_json)
                    .collect()
            ),
            expected["product_generators"]
        );
        let redundant =
            PathAlgebraIdeal::new(loop_algebra, vec![x.pow(2).unwrap(), x.pow(3).unwrap()])
                .unwrap();
        assert_eq!(
            Value::Array(
                redundant
                    .minimal_generating_set()
                    .unwrap()
                    .iter()
                    .map(rational_external_json)
                    .collect()
            ),
            expected["minimal_generators"]
        );

        let loop_algebra2 = PathAlgebra::new(Z2, loop_quiver).unwrap();
        let x2 = loop_algebra2
            .element(loop_algebra2.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let monomial2 =
            PathAlgebraIdeal::new(loop_algebra2.clone(), vec![x2.pow(2).unwrap()]).unwrap();
        let redundant2 =
            PathAlgebraIdeal::new(loop_algebra2, vec![x2.pow(2).unwrap(), x2.pow(3).unwrap()])
                .unwrap();
        assert_eq!(
            json!({
                "admissible": monomial2.is_admissible().unwrap(),
                "monomial": monomial2.is_monomial().unwrap(),
                "minimal_generators": redundant2
                    .minimal_generating_set()
                    .unwrap()
                    .iter()
                    .map(z2_external_json)
                    .collect::<Vec<_>>(),
            }),
            expected["f2"]
        );
    }

    #[test]
    #[allow(clippy::many_single_char_names, clippy::too_many_lines)]
    fn pinned_qpa_quotient_path_algebra_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/quotient_path_algebra.schema-v1.json"
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

        let quiver = loop_quiver();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let x = algebra
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y = algebra
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let rational_relators = vec![
            x.pow(2).unwrap(),
            x.multiply(&y)
                .unwrap()
                .sub(&y.multiply(&x).unwrap())
                .unwrap(),
            y.pow(2).unwrap(),
        ];
        let ideal = PathAlgebraIdeal::new(algebra.clone(), rational_relators.clone()).unwrap();
        let quotient =
            QuotientPathAlgebra::by_relations_high_level(algebra.clone(), &rational_relators)
                .unwrap();
        let generators = quotient.generators().unwrap();
        let x_quotient = &generators[1];
        let y_quotient = &generators[2];
        let one = quotient.one().unwrap();
        let yx = y_quotient.multiply(x_quotient).unwrap();
        let yxy_plus_x = yx.multiply(y_quotient).unwrap().add(x_quotient).unwrap();
        let product_with_inverse = one
            .add(x_quotient)
            .unwrap()
            .multiply(&one.sub(x_quotient).unwrap())
            .unwrap();
        let AlgebraDimension::Finite(rational_dimension) = quotient.dimension().unwrap() else {
            panic!("the pinned QPA quotient must be finite-dimensional");
        };
        let rational_basis = quotient.canonical_basis().unwrap().unwrap();
        assert_eq!(
            json!({
                "relations": quotient.relations().iter().map(rational_json).collect::<Vec<_>>(),
                "zero": rational_json(quotient.zero().representative()),
                "one": rational_json(one.representative()),
                "yx": rational_json(yx.representative()),
                "x_squared": rational_json(x_quotient.pow(2).unwrap().representative()),
                "yxy_plus_x": rational_json(yxy_plus_x.representative()),
                "one_plus_x_times_one_minus_x": rational_json(product_with_inverse.representative()),
                "finite_dimensional": quotient.is_finite_dimensional().unwrap(),
                "dimension": small_integer(&rational_dimension),
                "canonical_basis_paths": rational_basis.iter().map(|element| element.leading_monomial().to_string()).collect::<Vec<_>>(),
                "ideal_membership": {
                    "zero": ideal.contains(algebra.zero()).unwrap(),
                    "x": ideal.contains(x.clone()).unwrap(),
                    "x_cubed": ideal.contains(x.pow(3).unwrap()).unwrap()
                }
            }),
            expected["rational"]
        );

        let field5 = Zp::new(5);
        let algebra5 = PathAlgebra::new(field5.clone(), quiver.clone()).unwrap();
        let x5 = algebra5
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y5 = algebra5
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let f5_commutator = x5
            .multiply(&y5)
            .unwrap()
            .sub(
                &y5.multiply(&x5)
                    .unwrap()
                    .scale_left(&field5.to_element(2))
                    .unwrap(),
            )
            .unwrap();
        let quotient5 = QuotientPathAlgebra::by_relations_high_level(
            algebra5,
            &[x5.pow(2).unwrap(), f5_commutator, y5.pow(2).unwrap()],
        )
        .unwrap();
        let generators5 = quotient5.generators().unwrap();
        let x5 = &generators5[1];
        let y5 = &generators5[2];
        assert_eq!(
            json!({
                "yx": zp_json(y5.multiply(x5).unwrap().representative(), &field5),
                "x_squared": zp_json(x5.pow(2).unwrap().representative(), &field5)
            }),
            expected["f5"]
        );

        let algebra2 = PathAlgebra::new(Z2, quiver).unwrap();
        let x2 = algebra2
            .element(algebra2.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y2 = algebra2
            .element(algebra2.quiver().arrow(arrow(2)).unwrap().path())
            .unwrap();
        let f2_commutator = x2
            .multiply(&y2)
            .unwrap()
            .add(&y2.multiply(&x2).unwrap())
            .unwrap();
        let quotient2 = QuotientPathAlgebra::by_relations_high_level(
            algebra2,
            &[x2.pow(2).unwrap(), f2_commutator, y2.pow(2).unwrap()],
        )
        .unwrap();
        let generators2 = quotient2.generators().unwrap();
        let x2 = &generators2[1];
        let y2 = &generators2[2];
        assert_eq!(
            json!({
                "yx": z2_json(y2.multiply(x2).unwrap().representative()),
                "x_squared": z2_json(x2.pow(2).unwrap().representative())
            }),
            expected["f2"]
        );

        let regression_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(1), vertex(2), "b"),
                ArrowSpec::named(vertex(2), vertex(2), "c"),
                ArrowSpec::named(vertex(2), vertex(3), "d"),
                ArrowSpec::named(vertex(3), vertex(1), "e"),
            ],
        )
        .unwrap();
        let regression_parent = PathAlgebra::new(Q, regression_quiver).unwrap();
        let regression_arrows = (1..=5)
            .map(|index| {
                regression_parent
                    .element(
                        regression_parent
                            .quiver()
                            .arrow(arrow(index))
                            .unwrap()
                            .path(),
                    )
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let [a, b, c, d, e] = regression_arrows.as_slice() else {
            unreachable!();
        };
        let acd_minus_bd = a
            .multiply(c)
            .unwrap()
            .multiply(d)
            .unwrap()
            .sub(&b.multiply(d).unwrap())
            .unwrap();
        let regression = QuotientPathAlgebra::by_relations_high_level(
            regression_parent,
            &[
                d.multiply(e).unwrap(),
                c.pow(2).unwrap(),
                acd_minus_bd,
                e.multiply(a).unwrap(),
            ],
        )
        .unwrap();
        let AlgebraDimension::Finite(regression_dimension) = regression.dimension().unwrap() else {
            panic!("QPA's dimension-17 regression quotient must be finite");
        };
        assert_eq!(
            json!({
                "finite_dimensional": regression.is_finite_dimensional().unwrap(),
                "dimension": small_integer(&regression_dimension),
                "canonical_basis_paths": regression
                    .canonical_basis()
                    .unwrap()
                    .unwrap()
                    .iter()
                    .map(|element| element.leading_monomial().to_string())
                    .collect::<Vec<_>>()
            }),
            expected["qpa_regression_dimension_17"]
        );
    }
}
