//! Exact sparse path algebras over Numerica fields.
//!
//! QPA's `PathRing`/`PathAlgebra` methods in `lib/pathalg.gi` delegate sparse
//! arithmetic to GAP's magma-ring implementation in `lib/mgmring.gi`. This
//! module translates that representation and its `ZippedSum`/`ZippedProduct`
//! semantics: terms are strictly path-sorted, equal paths are combined, zero
//! coefficients disappear, and QPA's zero path is quotiented out.

use std::{fmt, sync::Arc};

use numerica::domains::Field;
use thiserror::Error;

use crate::{OrderingError, Path, PathEndpoint, PathError, Quiver, QuiverOrdering, Vertex};

struct PathAlgebraData<F: Field> {
    field: F,
    quiver: Quiver,
    ordering: QuiverOrdering,
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
        let ordering = QuiverOrdering::default_for(&quiver)?;
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

    /// Return the path ordering used to canonicalize terms.
    #[must_use]
    pub fn ordering(&self) -> &QuiverOrdering {
        &self.0.ordering
    }

    /// Return whether two handles denote the identical QPA algebra object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
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

    /// Return QPA's path-algebra global-dimension value set by the constructor.
    #[must_use]
    pub fn qpa_global_dimension(&self) -> usize {
        usize::from(self.quiver().number_of_arrows() > 0)
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
    /// A source-translated ordering operation failed.
    #[error(transparent)]
    Ordering(#[from] OrderingError),
    /// QPA path multiplication failed.
    #[error(transparent)]
    Path(#[from] PathError),
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
