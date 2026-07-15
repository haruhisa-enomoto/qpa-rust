//! QPA's high-level noncommutative Groebner-basis routines.
//!
//! The algorithms here translate `lib/gbhighlevel.gi` from QPA 1.37:
//! `HighLevelGroebnerBasis`, `RemainderOfDivision`, `ReducedListQPA`,
//! `TipReducedList`, and `LeftmostOccurrence`. They intentionally operate on
//! QPA paths and QPA's monomial ordering rather than substituting a generic
//! polynomial implementation.

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    ArrowIndex, Path, PathAlgebra, PathAlgebraElement, PathAlgebraError, PathEndpoint,
    QuiverStaticDictionary, StaticDictionaryError,
};

/// A complete Groebner basis produced by QPA's high-level overlap algorithm.
#[derive(Clone, Debug)]
pub struct GroebnerBasis<F: Field> {
    algebra: PathAlgebra<F>,
    relations: Vec<PathAlgebraElement<F>>,
}

impl<F: Field> GroebnerBasis<F> {
    /// Run QPA's `HighLevelGroebnerBasis` implementation.
    ///
    /// # Errors
    ///
    /// Returns an error for relations from another algebra, relations outside
    /// the arrow ideal, or a failure in source-translated arithmetic.
    pub fn high_level(
        algebra: &PathAlgebra<F>,
        relations: &[PathAlgebraElement<F>],
    ) -> Result<Self, GroebnerError> {
        let relations = high_level_groebner_basis(relations, algebra)?;
        Ok(Self {
            algebra: algebra.clone(),
            relations,
        })
    }

    /// Return the containing path algebra.
    #[must_use]
    pub const fn algebra(&self) -> &PathAlgebra<F> {
        &self.algebra
    }

    /// Return QPA's final reduced relation list.
    #[must_use]
    pub fn relations(&self) -> &[PathAlgebraElement<F>] {
        &self.relations
    }

    /// Compute QPA's `RemainderOfDivision` by this basis.
    ///
    /// # Errors
    ///
    /// Propagates relation-family and arithmetic failures.
    pub fn remainder(
        &self,
        element: PathAlgebraElement<F>,
    ) -> Result<PathAlgebraElement<F>, GroebnerError> {
        remainder_of_division(element, &self.relations, &self.algebra)
    }

    /// Build QPA's `QuiverStaticDictionary` from leading monomials.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary validation failures.
    pub fn static_dictionary(&self) -> Result<QuiverStaticDictionary, GroebnerError> {
        let patterns = self
            .relations
            .iter()
            .map(PathAlgebraElement::leading_monomial)
            .collect::<Vec<_>>();
        Ok(QuiverStaticDictionary::new(
            self.algebra.quiver(),
            &patterns,
        )?)
    }

    /// Translate QPA's `AdmitsFinitelyManyNontips`.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary construction failures.
    pub fn admits_finitely_many_nontips(&self) -> Result<bool, GroebnerError> {
        Ok(self.static_dictionary()?.is_finite_difference())
    }

    /// Translate QPA's `NontipSize`; `None` represents `infinity`.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary construction failures.
    pub fn nontip_size(
        &self,
    ) -> Result<Option<numerica::domains::integer::Integer>, GroebnerError> {
        Ok(self.static_dictionary()?.difference_size())
    }

    /// Translate QPA's `Nontips`; `None` represents `fail` for an infinite
    /// nontip language.
    ///
    /// # Errors
    ///
    /// Propagates static-dictionary construction failures.
    pub fn nontips(&self) -> Result<Option<Vec<Path>>, GroebnerError> {
        Ok(self.static_dictionary()?.difference_words())
    }
}

/// Translate QPA's `HighLevelGroebnerBasis` overlap-completion algorithm.
///
/// This routine can be nonterminating when the ideal has no finite Groebner
/// basis for the chosen ordering, exactly as its QPA source can be.
///
/// # Errors
///
/// Returns an error for relations from another algebra, relations outside the
/// arrow ideal, or a failure in source-translated arithmetic.
pub fn high_level_groebner_basis<F: Field>(
    relations: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, GroebnerError> {
    if !qpa_in_arrow_ideal(relations, algebra)? {
        return Err(GroebnerError::OutsideArrowIdeal);
    }

    let mut pending = make_uniform_list(relations, algebra)?;
    let mut basis = Vec::new();

    while !pending.is_empty() {
        for relation in pending {
            let leading_coefficient = relation.leading_coefficient();
            basis.push(relation.divide_scalar(&leading_coefficient)?);
        }
        basis = reduced_list_qpa(&basis, algebra)?;

        let number_of_relations = basis.len();
        pending = Vec::new();
        for left_index in 0..number_of_relations {
            let left_walk = tip_walk(&basis[left_index])?;
            let left_length = left_walk.len();

            for right_index in 0..number_of_relations {
                let right_walk = tip_walk(&basis[right_index])?;
                let right_length = right_walk.len();
                let first_overlap = left_length.saturating_sub(right_length);

                for split in first_overlap..left_length {
                    let overlap_length = left_length - split;
                    if left_walk[split..] != right_walk[..overlap_length] {
                        continue;
                    }

                    let left_context = context_from_walk(algebra, &left_walk[..split])?;
                    let right_context = context_from_walk(algebra, &right_walk[overlap_length..])?;
                    let overlap = basis[left_index]
                        .multiply(&right_context)?
                        .sub(&left_context.multiply(&basis[right_index])?)?;
                    let remainder = remainder_of_division(overlap, &basis, algebra)?;
                    if !remainder.is_zero() {
                        add_set(&mut pending, remainder)?;
                    }
                }
            }
        }
    }

    let basis = tip_reduced_list(&basis, algebra)?;
    reduced_list_qpa(&basis, algebra)
}

/// Translate QPA's noncommutative `RemainderOfDivision`.
///
/// The first divisor whose tip occurs in the current tip is used, and the
/// leftmost occurrence is selected.
///
/// # Errors
///
/// Returns an error for elements from another algebra, zero divisors in the
/// divisor list, or source-translated arithmetic failures.
pub fn remainder_of_division<F: Field>(
    mut dividend: PathAlgebraElement<F>,
    divisors: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<PathAlgebraElement<F>, GroebnerError> {
    ensure_algebra(&dividend, algebra)?;
    for divisor in divisors {
        ensure_algebra(divisor, algebra)?;
        if divisor.is_zero() {
            return Err(GroebnerError::ZeroDivisor);
        }
    }

    let mut remainder = algebra.zero();
    while !dividend.is_zero() {
        let leading_term = dividend.leading_term();
        let dividend_walk = tip_walk(&leading_term)?;
        let mut reduction_applied = false;

        for divisor in divisors {
            let divisor_walk = tip_walk(divisor)?;
            let Some((start, end)) = leftmost_occurrence(&dividend_walk, &divisor_walk) else {
                continue;
            };

            let left = context_from_walk(algebra, &dividend_walk[..start])?;
            let right = context_from_walk(algebra, &dividend_walk[end + 1..])?;
            let coefficient = algebra.field().div(
                &leading_term.leading_coefficient(),
                &divisor.leading_coefficient(),
            );
            let reduction = left
                .multiply(divisor)?
                .multiply(&right)?
                .scale_left(&coefficient)?;
            dividend = dividend.sub(&reduction)?;
            reduction_applied = true;
            break;
        }

        if !reduction_applied {
            remainder = remainder.add(&leading_term)?;
            dividend = dividend.sub(&leading_term)?;
        }
    }

    Ok(remainder)
}

/// Translate QPA's `ReducedListQPA`.
///
/// # Errors
///
/// Returns an error for elements from another algebra or arithmetic failures.
pub fn reduced_list_qpa<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, GroebnerError> {
    ensure_elements(elements, algebra)?;
    let mut result = elements
        .iter()
        .filter(|element| !element.is_zero())
        .cloned()
        .collect::<Vec<_>>();

    let mut index = result.len();
    while index > 0 {
        let current = index - 1;
        let others = result
            .iter()
            .enumerate()
            .filter(|(candidate, _)| *candidate != current)
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>();
        let reduced = remainder_of_division(result[current].clone(), &others, algebra)?;
        if reduced.is_zero() {
            result.remove(current);
        } else {
            result[current] = reduced;
        }
        index -= 1;
    }

    Ok(result)
}

/// Translate QPA's `TipReducedList`.
///
/// # Errors
///
/// Returns an error for elements from another algebra or invalid zero tips.
pub fn tip_reduced_list<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, GroebnerError> {
    ensure_elements(elements, algebra)?;
    let mut result = Vec::new();
    for element in elements.iter().filter(|element| !element.is_zero()) {
        add_set(&mut result, element.clone())?;
    }

    let mut index = result.len();
    while index > 0 {
        let current = index - 1;
        let current_tip = tip_walk(&result[current])?;
        let mut contained = false;
        for prior in &result[..current] {
            if leftmost_occurrence(&current_tip, &tip_walk(prior)?).is_some() {
                contained = true;
                break;
            }
        }
        if contained {
            result.remove(current);
        }
        index -= 1;
    }
    Ok(result)
}

/// Translate QPA's `LeftmostOccurrence` using zero-based inclusive indices.
#[must_use]
pub fn leftmost_occurrence<T: PartialEq>(sequence: &[T], pattern: &[T]) -> Option<(usize, usize)> {
    if pattern.is_empty() || pattern.len() > sequence.len() {
        return None;
    }
    sequence
        .windows(pattern.len())
        .position(|window| window == pattern)
        .map(|start| (start, start + pattern.len() - 1))
}

/// Translate `QPA_InArrowIdeal` for a list of relations.
///
/// # Errors
///
/// Returns an error when a relation belongs to another algebra.
pub fn qpa_in_arrow_ideal<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<bool, GroebnerError> {
    ensure_elements(elements, algebra)?;
    Ok(elements
        .iter()
        .flat_map(PathAlgebraElement::terms)
        .all(|term| !term.path().is_vertex()))
}

/// Translate QPA's `MakeUniform` for one path-algebra element.
///
/// # Errors
///
/// Returns an error for an element from another algebra or arithmetic failure.
pub fn make_uniform<F: Field>(
    element: &PathAlgebraElement<F>,
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, GroebnerError> {
    ensure_algebra(element, algebra)?;
    let mut sources = Vec::new();
    let mut targets = Vec::new();
    for term in element.terms() {
        if let PathEndpoint::Vertex(source) = term.path().source()
            && !sources.contains(&source)
        {
            sources.push(source);
        }
        if let PathEndpoint::Vertex(target) = term.path().target()
            && !targets.contains(&target)
        {
            targets.push(target);
        }
    }

    let mut uniform = Vec::new();
    for source in &sources {
        let source_idempotent = algebra.element(source.path())?;
        for target in &targets {
            let target_idempotent = algebra.element(target.path())?;
            let component = source_idempotent
                .multiply(element)?
                .multiply(&target_idempotent)?;
            if !component.is_zero() {
                add_set(&mut uniform, component)?;
            }
        }
    }
    Ok(uniform)
}

fn make_uniform_list<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, GroebnerError> {
    let mut result = Vec::new();
    for element in elements {
        for component in make_uniform(element, algebra)? {
            add_set(&mut result, component)?;
        }
    }
    Ok(result)
}

fn context_from_walk<F: Field>(
    algebra: &PathAlgebra<F>,
    walk: &[ArrowIndex],
) -> Result<PathAlgebraElement<F>, GroebnerError> {
    if walk.is_empty() {
        return Ok(algebra.one());
    }
    let arrows = walk
        .iter()
        .map(|index| {
            algebra
                .quiver()
                .arrow(*index)
                .ok_or(GroebnerError::InvalidTipWalk)
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(algebra.element(Path::from_arrow_list(algebra.quiver(), &arrows)?)?)
}

fn tip_walk<F: Field>(element: &PathAlgebraElement<F>) -> Result<Vec<ArrowIndex>, GroebnerError> {
    if element.is_zero() {
        return Err(GroebnerError::ZeroDivisor);
    }
    element
        .leading_monomial()
        .walk()
        .map(<[ArrowIndex]>::to_vec)
        .ok_or(GroebnerError::InvalidTipWalk)
}

fn ensure_algebra<F: Field>(
    element: &PathAlgebraElement<F>,
    algebra: &PathAlgebra<F>,
) -> Result<(), GroebnerError> {
    if element.algebra().is_identical(algebra) {
        Ok(())
    } else {
        Err(GroebnerError::DifferentAlgebras)
    }
}

fn ensure_elements<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<(), GroebnerError> {
    for element in elements {
        ensure_algebra(element, algebra)?;
    }
    Ok(())
}

fn add_set<F: Field>(
    elements: &mut Vec<PathAlgebraElement<F>>,
    element: PathAlgebraElement<F>,
) -> Result<(), GroebnerError> {
    if elements.contains(&element) {
        return Ok(());
    }
    elements.push(element);
    let mut comparison_error = None;
    elements.sort_by(|left, right| match left.compare(right) {
        Ok(ordering) => ordering,
        Err(error) => {
            comparison_error = Some(error);
            std::cmp::Ordering::Equal
        }
    });
    comparison_error.map_or(Ok(()), |error| Err(error.into()))
}

/// Failure in QPA's source-translated Groebner-basis routines.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum GroebnerError {
    /// A relation belongs to another path-algebra family.
    #[error("a relation belongs to a different path algebra")]
    DifferentAlgebras,
    /// `HighLevelGroebnerBasis` accepts only relations in the arrow ideal.
    #[error("relations do not belong to the arrow ideal of the path algebra")]
    OutsideArrowIdeal,
    /// A zero relation cannot serve as a division rule.
    #[error("a Groebner divisor relation is zero")]
    ZeroDivisor,
    /// A tip was not a valid path walk.
    #[error("a relation tip does not have a valid arrow walk")]
    InvalidTipWalk,
    /// Source-translated path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// QPA path construction failed.
    #[error(transparent)]
    Path(#[from] crate::PathError),
    /// QPA static-dictionary construction failed.
    #[error(transparent)]
    StaticDictionary(#[from] StaticDictionaryError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, RationalField},
    };
    use serde_json::{Value, json};

    use crate::{ArrowIndex, ArrowSpec, PathAlgebra, PathAlgebraElement, Quiver, VertexIndex};

    use super::{GroebnerBasis, leftmost_occurrence, make_uniform, remainder_of_division};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn arrow(index: usize) -> ArrowIndex {
        ArrowIndex::new(index).unwrap()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the Groebner fixture uses only small integers")
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

    #[test]
    fn leftmost_occurrence_translates_qpa_inclusive_result() {
        assert_eq!(leftmost_occurrence(&[1, 2, 1, 2], &[1, 2]), Some((0, 1)));
        assert_eq!(leftmost_occurrence(&[1, 2], &[2, 1]), None);
        assert_eq!(leftmost_occurrence::<u8>(&[], &[]), None);
    }

    #[test]
    fn make_uniform_splits_source_target_components() {
        let quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let algebra = PathAlgebra::new(Q, quiver).unwrap();
        let a = algebra
            .element(algebra.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let b = algebra
            .element(algebra.quiver().arrow(arrow(2)).unwrap().path())
            .unwrap();
        let components = make_uniform(&a.add(&b).unwrap(), &algebra).unwrap();
        assert_eq!(components, vec![a, b]);
    }

    #[test]
    fn division_and_overlap_completion_follow_qpa_high_level_algorithm() {
        let quiver = loop_quiver();
        let algebra = PathAlgebra::new(Q, quiver).unwrap();
        let x = algebra
            .element(algebra.quiver().arrow(arrow(1)).unwrap().path())
            .unwrap();
        let y = algebra
            .element(algebra.quiver().arrow(arrow(2)).unwrap().path())
            .unwrap();
        let xx = x.multiply(&x).unwrap();
        let xy_minus_yx = x
            .multiply(&y)
            .unwrap()
            .sub(&y.multiply(&x).unwrap())
            .unwrap();
        let yy = y.multiply(&y).unwrap();
        let basis =
            GroebnerBasis::high_level(&algebra, &[xx.clone(), xy_minus_yx, yy.clone()]).unwrap();

        let word = y.multiply(&x).unwrap().multiply(&y).unwrap();
        assert!(basis.remainder(word).unwrap().is_zero());
        assert!(
            remainder_of_division(xx, basis.relations(), &algebra)
                .unwrap()
                .is_zero()
        );
        assert_eq!(basis.relations().len(), 3);
        assert_eq!(Q.one(), algebra.field().one());
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_high_level_groebner_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/groebner_high_level.schema-v1.json"
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
        let rational_basis = GroebnerBasis::high_level(
            &algebra,
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
        let yx = y.multiply(&x).unwrap();
        let yxy_plus_x = yx.multiply(&y).unwrap().add(&x).unwrap();
        assert_eq!(
            json!({
                "basis": rational_basis.relations().iter().map(rational_json).collect::<Vec<_>>(),
                "remainder_yx": rational_json(&rational_basis.remainder(yx).unwrap()),
                "remainder_yxy_plus_x": rational_json(&rational_basis.remainder(yxy_plus_x).unwrap())
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
        let relation5 = x5
            .multiply(&y5)
            .unwrap()
            .sub(
                &y5.multiply(&x5)
                    .unwrap()
                    .scale_left(&field5.to_element(2))
                    .unwrap(),
            )
            .unwrap();
        let basis5 = GroebnerBasis::high_level(
            &algebra5,
            &[x5.pow(2).unwrap(), relation5, y5.pow(2).unwrap()],
        )
        .unwrap();
        assert_eq!(
            json!({
                "basis": basis5.relations().iter().map(|element| zp_json(element, &field5)).collect::<Vec<_>>(),
                "remainder_yx": zp_json(&basis5.remainder(y5.multiply(&x5).unwrap()).unwrap(), &field5)
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
        let relation2 = x2
            .multiply(&y2)
            .unwrap()
            .add(&y2.multiply(&x2).unwrap())
            .unwrap();
        let basis2 = GroebnerBasis::high_level(
            &algebra2,
            &[x2.pow(2).unwrap(), relation2, y2.pow(2).unwrap()],
        )
        .unwrap();
        assert_eq!(
            json!({
                "basis": basis2.relations().iter().map(z2_json).collect::<Vec<_>>(),
                "remainder_yx": z2_json(&basis2.remainder(y2.multiply(&x2).unwrap()).unwrap())
            }),
            expected["f2"]
        );

        let uniform_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "a"),
                ArrowSpec::named(vertex(1), vertex(2), "b"),
                ArrowSpec::named(vertex(2), vertex(1), "c"),
            ],
        )
        .unwrap();
        let uniform_algebra = PathAlgebra::new(Q, uniform_quiver).unwrap();
        let mut sum = uniform_algebra.zero();
        for index in 1..=3 {
            sum = sum
                .add(
                    &uniform_algebra
                        .element(uniform_algebra.quiver().arrow(arrow(index)).unwrap().path())
                        .unwrap(),
                )
                .unwrap();
        }
        assert_eq!(
            json!(
                make_uniform(&sum, &uniform_algebra)
                    .unwrap()
                    .iter()
                    .map(rational_json)
                    .collect::<Vec<_>>()
            ),
            expected["uniform_components"]
        );
    }
}
