//! Ideals and quotients of QPA path algebras.
//!
//! This module translates the foundational behavior in `pathalgideal.gi` and
//! `pathalg.gi`: two-sided ideals remember their path algebra and generators,
//! membership is a Groebner normal-form test, and `KQ / relations` first
//! computes a Groebner basis and then stores quotient elements in normal form.

use std::{cmp::Ordering, fmt, sync::Arc};

use numerica::domains::{Field, integer::Integer};
use thiserror::Error;

use crate::{
    DenseMatrix, GroebnerBasis, GroebnerError, MatrixError, Path, PathAlgebra, PathAlgebraElement,
    PathAlgebraError,
};

struct PathAlgebraIdealData<F: Field> {
    algebra: PathAlgebra<F>,
    generators: Vec<PathAlgebraElement<F>>,
    groebner_basis: Option<GroebnerBasis<F>>,
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
            groebner_basis: None,
        })))
    }

    fn from_groebner_basis(algebra: PathAlgebra<F>, basis: GroebnerBasis<F>) -> Self {
        let generators = basis.relations().to_vec();
        Self(Arc::new(PathAlgebraIdealData {
            algebra,
            generators,
            groebner_basis: Some(basis),
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
    /// The returned handle is a new ideal object with the original generator
    /// list and an attached complete basis.
    ///
    /// # Errors
    ///
    /// Propagates QPA high-level Groebner failures.
    pub fn with_high_level_groebner_basis(&self) -> Result<Self, QuotientError> {
        if self.0.groebner_basis.is_some() {
            return Ok(self.clone());
        }
        let basis = GroebnerBasis::high_level(self.algebra(), self.generators())?;
        Ok(Self(Arc::new(PathAlgebraIdealData {
            algebra: self.algebra().clone(),
            generators: self.generators().to_vec(),
            groebner_basis: Some(basis),
        })))
    }

    /// Attach QPA `GroebnerBasis(I, relations)` to this ideal.
    ///
    /// The returned handle is a new immutable Rust ideal object carrying the
    /// complete basis; the original generator list is preserved.
    ///
    /// # Errors
    ///
    /// Returns an error for foreign relations or QPA comparison failures.
    pub fn with_complete_groebner_basis(
        &self,
        relations: &[PathAlgebraElement<F>],
    ) -> Result<Self, QuotientError> {
        let basis = GroebnerBasis::from_complete_relations(self.algebra(), relations)?;
        Ok(Self(Arc::new(PathAlgebraIdealData {
            algebra: self.algebra().clone(),
            generators: self.generators().to_vec(),
            groebner_basis: Some(basis),
        })))
    }

    /// Return an attached Groebner basis, if one has been computed.
    #[must_use]
    pub fn groebner_basis(&self) -> Option<&GroebnerBasis<F>> {
        self.0.groebner_basis.as_ref()
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
        let computed;
        let basis = if let Some(basis) = self.groebner_basis() {
            basis
        } else {
            computed = GroebnerBasis::high_level(self.algebra(), self.generators())?;
            &computed
        };
        Ok(basis.remainder(element)?.is_zero())
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
        let computed;
        let basis = if let Some(basis) = self.groebner_basis() {
            basis
        } else {
            computed = GroebnerBasis::high_level(self.algebra(), self.generators())?;
            &computed
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
        if self.generators().is_empty() {
            return Ok(Vec::new());
        }
        if !self.is_admissible()? {
            return Err(QuotientError::IdealNotAdmissible);
        }

        let arrows = self
            .algebra()
            .quiver()
            .arrows()
            .map(|arrow| self.algebra().element(arrow.path()))
            .collect::<Result<Vec<_>, _>>()?;
        let mut radical_products = Vec::new();
        for arrow in &arrows {
            for generator in self.generators() {
                let product = generator.multiply(arrow)?;
                if !product.is_zero() {
                    radical_products.push(product);
                }
            }
        }
        for arrow in &arrows {
            for generator in self.generators() {
                let product = arrow.multiply(generator)?;
                if !product.is_zero() {
                    radical_products.push(product);
                }
            }
        }

        if radical_products.is_empty() {
            let algebra_basis = self
                .algebra()
                .canonical_basis()
                .ok_or(QuotientError::InfiniteParentBasis)?;
            return path_algebra_ideal_span(self.algebra(), self.generators(), &algebra_basis);
        }

        let quotient = QuotientPathAlgebra::by_relations_high_level(
            self.algebra().clone(),
            &radical_products,
        )?;
        let algebra_basis = quotient
            .canonical_basis()?
            .ok_or(QuotientError::InfiniteQuotientBasis)?;
        let generators = self
            .generators()
            .iter()
            .cloned()
            .map(|generator| quotient.element(generator))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(quotient_ideal_span(&quotient, &generators, &algebra_basis)?
            .into_iter()
            .map(|element| element.normal_form().clone())
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
    groebner_basis: GroebnerBasis<F>,
}

/// A quotient `KQ/I` whose elements are stored in QPA Groebner normal form.
#[derive(Clone)]
pub struct QuotientPathAlgebra<F: Field>(Arc<QuotientPathAlgebraData<F>>);

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
        let groebner_basis = GroebnerBasis::high_level(&parent, relators)?;
        let ideal = PathAlgebraIdeal::from_groebner_basis(parent.clone(), groebner_basis.clone());
        Ok(Self(Arc::new(QuotientPathAlgebraData {
            parent,
            ideal,
            groebner_basis,
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
            .cloned()
            .ok_or(QuotientError::MissingGroebnerBasis)?;
        Ok(Self(Arc::new(QuotientPathAlgebraData {
            parent,
            ideal,
            groebner_basis,
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
    pub fn groebner_basis(&self) -> &GroebnerBasis<F> {
        &self.0.groebner_basis
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
            normal_form: self.original_path_algebra().zero(),
        }
    }

    /// Return the quotient identity.
    ///
    /// # Errors
    ///
    /// Propagates normal-form failures.
    pub fn one(&self) -> Result<QuotientPathAlgebraElement<F>, QuotientError> {
        self.element(self.original_path_algebra().one())
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
        ensure_parent_element(self.original_path_algebra(), &element)?;
        let normal_form = self.groebner_basis().remainder(element)?;
        Ok(QuotientPathAlgebraElement {
            algebra: self.clone(),
            normal_form,
        })
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
        Ok(self.groebner_basis().admits_finitely_many_nontips()?)
    }

    /// Translate QPA's `Dimension`; infinite dimension is explicit.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary failures.
    pub fn dimension(&self) -> Result<AlgebraDimension, QuotientError> {
        Ok(self
            .groebner_basis()
            .nontip_size()?
            .map_or(AlgebraDimension::Infinite, AlgebraDimension::Finite))
    }

    /// Translate QPA's canonical basis of a finite-dimensional quotient.
    ///
    /// `None` corresponds to an infinite nontip language.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary, embedding, and normal-form failures.
    pub fn canonical_basis(
        &self,
    ) -> Result<Option<Vec<QuotientPathAlgebraElement<F>>>, QuotientError> {
        self.groebner_basis()
            .nontips()?
            .map(|paths| paths.into_iter().map(|path| self.path(path)).collect())
            .transpose()
    }

    /// Translate QPA's `IsAdmissibleQuotientOfPathAlgebra`.
    ///
    /// # Errors
    ///
    /// Propagates the defining ideal's admissibility computation.
    pub fn is_admissible(&self) -> Result<bool, QuotientError> {
        self.ideal().is_admissible()
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

/// A normalized element of a quotient of a path algebra.
#[derive(Clone, Debug)]
pub struct QuotientPathAlgebraElement<F: Field> {
    algebra: QuotientPathAlgebra<F>,
    normal_form: PathAlgebraElement<F>,
}

impl<F: Field> QuotientPathAlgebraElement<F> {
    /// Return the containing quotient algebra.
    #[must_use]
    pub const fn algebra(&self) -> &QuotientPathAlgebra<F> {
        &self.algebra
    }

    /// Return the normalized underlying path-algebra element.
    #[must_use]
    pub const fn normal_form(&self) -> &PathAlgebraElement<F> {
        &self.normal_form
    }

    /// Return whether this is quotient zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.normal_form.is_zero()
    }

    /// Add and normalize quotient elements.
    ///
    /// # Errors
    ///
    /// Returns an error for different quotient families or reduction failure.
    pub fn add(&self, right: &Self) -> Result<Self, QuotientError> {
        self.ensure_same_quotient(right)?;
        self.algebra
            .element(self.normal_form.add(&right.normal_form)?)
    }

    /// Return the additive inverse in normal form.
    ///
    /// # Errors
    ///
    /// Propagates arithmetic and reduction failures.
    pub fn neg(&self) -> Result<Self, QuotientError> {
        self.algebra.element(self.normal_form.neg()?)
    }

    /// Subtract and normalize quotient elements.
    ///
    /// # Errors
    ///
    /// Returns an error for different quotient families or reduction failure.
    pub fn sub(&self, right: &Self) -> Result<Self, QuotientError> {
        self.ensure_same_quotient(right)?;
        self.algebra
            .element(self.normal_form.sub(&right.normal_form)?)
    }

    /// Multiply and normalize quotient elements.
    ///
    /// # Errors
    ///
    /// Returns an error for different quotient families or reduction failure.
    pub fn multiply(&self, right: &Self) -> Result<Self, QuotientError> {
        self.ensure_same_quotient(right)?;
        self.algebra
            .element(self.normal_form.multiply(&right.normal_form)?)
    }

    /// Scale on the left and normalize.
    ///
    /// # Errors
    ///
    /// Propagates arithmetic and reduction failures.
    pub fn scale_left(&self, scalar: &F::Element) -> Result<Self, QuotientError> {
        self.algebra.element(self.normal_form.scale_left(scalar)?)
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
            normal_form: self.normal_form.leading_term(),
        }
    }

    /// Return the leading coefficient of the normal representative.
    #[must_use]
    pub fn leading_coefficient(&self) -> F::Element {
        self.normal_form.leading_coefficient()
    }

    /// Return the leading monomial of the normal representative.
    #[must_use]
    pub fn leading_monomial(&self) -> Path {
        self.normal_form.leading_monomial()
    }

    /// Translate QPA's left-uniform test on the representative.
    #[must_use]
    pub fn is_left_uniform(&self) -> bool {
        self.normal_form.is_left_uniform()
    }

    /// Translate QPA's right-uniform test on the representative.
    #[must_use]
    pub fn is_right_uniform(&self) -> bool {
        self.normal_form.is_right_uniform()
    }

    /// Translate QPA's uniform test on the representative.
    #[must_use]
    pub fn is_uniform(&self) -> bool {
        self.normal_form.is_uniform()
    }

    /// Compare normal representatives using QPA's path-algebra comparison.
    ///
    /// # Errors
    ///
    /// Returns an error for different quotient families.
    pub fn compare(&self, right: &Self) -> Result<Ordering, QuotientError> {
        self.ensure_same_quotient(right)?;
        Ok(self.normal_form.compare(&right.normal_form)?)
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
        self.algebra == other.algebra && self.normal_form == other.normal_form
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
        .map(|element| canonical_basis_path(element.normal_form()))
        .collect::<Result<Vec<_>, _>>()?;
    let rows = elements
        .iter()
        .map(|element| element_coordinates(element.normal_form(), &basis_paths))
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

/// Failure to construct or operate on a QPA ideal or quotient path algebra.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum QuotientError {
    /// An ideal generator or representative belongs to another parent algebra.
    #[error("an element belongs to a different parent path algebra")]
    DifferentParentAlgebras,
    /// Quotient arithmetic requires identical quotient families.
    #[error("the elements belong to different quotient path algebras")]
    DifferentQuotients,
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
    /// QPA declares homogeneous-list checking only for positive degrees.
    #[error("the homogeneous degree must be positive")]
    NonPositiveDegree,
    /// QPA Groebner computation or reduction failed.
    #[error(transparent)]
    Groebner(#[from] GroebnerError),
    /// Parent path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Exact GAP-compatible row reduction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, RationalField},
    };
    use serde_json::{Value, json};

    use crate::{ArrowIndex, ArrowSpec, PathAlgebra, PathAlgebraElement, Quiver, VertexIndex};

    use super::{
        AlgebraDimension, PathAlgebraIdeal, QuotientError, QuotientPathAlgebra,
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
                "zero": rational_json(quotient.zero().normal_form()),
                "one": rational_json(one.normal_form()),
                "yx": rational_json(yx.normal_form()),
                "x_squared": rational_json(x_quotient.pow(2).unwrap().normal_form()),
                "yxy_plus_x": rational_json(yxy_plus_x.normal_form()),
                "one_plus_x_times_one_minus_x": rational_json(product_with_inverse.normal_form()),
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
                "yx": zp_json(y5.multiply(x5).unwrap().normal_form(), &field5),
                "x_squared": zp_json(x5.pow(2).unwrap().normal_form(), &field5)
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
                "yx": z2_json(y2.multiply(x2).unwrap().normal_form()),
                "x_squared": z2_json(x2.pow(2).unwrap().normal_form())
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
