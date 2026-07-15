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
    GroebnerBasis, GroebnerError, Path, PathAlgebra, PathAlgebraElement, PathAlgebraError,
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

    /// Return an attached Groebner basis, if one has been computed.
    #[must_use]
    pub fn groebner_basis(&self) -> Option<&GroebnerBasis<F>> {
        self.0.groebner_basis.as_ref()
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
            .all(|generator| generator.terms().len() <= 1)
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
    /// An internal quotient construction invariant did not supply a basis.
    #[error("the quotient ideal has no Groebner basis")]
    MissingGroebnerBasis,
    /// QPA Groebner computation or reduction failed.
    #[error(transparent)]
    Groebner(#[from] GroebnerError),
    /// Parent path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
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

    use super::{AlgebraDimension, PathAlgebraIdeal, QuotientPathAlgebra};

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
