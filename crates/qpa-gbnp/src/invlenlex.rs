//! QPA's inverse-length/inverse-lexicographic Gröbner routines.
//!
//! This is a source translation of `lib/invlenlexgrb.gi`. In particular,
//! inverse tips are the *first* terms in QPA's canonical sparse term list,
//! and the complete constructor obtains its truncation bound through the same
//! default GBNP quotient route used by QPA.

use numerica::domains::Field;
use qpa_core::{
    ArrowIndex, DenseMatrix, GroebnerBasis, GroebnerError, MatrixError, Path, PathAlgebra,
    PathAlgebraElement, PathAlgebraError, PathError, leftmost_occurrence, make_uniform_list,
    qpa_in_arrow_ideal,
};
use thiserror::Error;

use crate::{QpaGbnpError, gbnp_groebner_basis};

/// Translate QPA `TipInvLenInvLex`.
///
/// # Errors
///
/// Propagates sparse path-algebra construction failures.
pub fn tip_inv_len_inv_lex<F: Field>(
    element: &PathAlgebraElement<F>,
) -> Result<PathAlgebraElement<F>, InvLenLexError> {
    let Some(term) = element.terms().first() else {
        return Ok(element.algebra().zero());
    };
    Ok(element
        .algebra()
        .from_terms([(term.path().clone(), term.coefficient().clone())])?)
}

/// Translate QPA `TipCoefficientInvLenInvLex`.
#[must_use]
pub fn tip_coefficient_inv_len_inv_lex<F: Field>(element: &PathAlgebraElement<F>) -> F::Element {
    element.terms().first().map_or_else(
        || element.algebra().field().zero(),
        |term| term.coefficient().clone(),
    )
}

/// Translate QPA `TipMonomialInvLenInvLex`.
///
/// `None` is Rust's typed replacement for QPA's zero-element sentinel.
#[must_use]
pub fn tip_monomial_inv_len_inv_lex<F: Field>(element: &PathAlgebraElement<F>) -> Option<Path> {
    element.terms().first().map(|term| term.path().clone())
}

/// Translate QPA `TipWalkInvLenInvLex`.
///
/// `None` represents the zero input, for which QPA's return crosses its path
/// and path-algebra-element representations.
#[must_use]
pub fn tip_walk_inv_len_inv_lex<F: Field>(
    element: &PathAlgebraElement<F>,
) -> Option<Vec<ArrowIndex>> {
    tip_monomial_inv_len_inv_lex(element)?
        .walk()
        .map(<[ArrowIndex]>::to_vec)
}

/// Translate `TruncateElement(element, degree)` modulo the arrow ideal power.
///
/// # Errors
///
/// QPA requires a positive degree. Path-algebra construction errors propagate.
pub fn truncate_element<F: Field>(
    element: &PathAlgebraElement<F>,
    degree: usize,
) -> Result<PathAlgebraElement<F>, InvLenLexError> {
    if degree == 0 {
        return Err(InvLenLexError::NonPositiveTruncationDegree);
    }
    let limit = element
        .terms()
        .iter()
        .filter(|term| term.path().length().is_some_and(|length| length < degree))
        .count();
    Ok(element.algebra().from_terms(
        element
            .terms()
            .iter()
            .take(limit)
            .map(|term| (term.path().clone(), term.coefficient().clone())),
    )?)
}

/// Translate element-level `BiggestHomogeneousPartInvLenInvLex`.
///
/// # Errors
///
/// Propagates path-algebra construction failures.
pub fn biggest_homogeneous_part_inv_len_inv_lex<F: Field>(
    element: &PathAlgebraElement<F>,
) -> Result<PathAlgebraElement<F>, InvLenLexError> {
    let Some(first) = element.terms().first() else {
        return Ok(element.clone());
    };
    let degree = first
        .path()
        .length()
        .ok_or(InvLenLexError::InvalidTipWalk)?;
    Ok(element.algebra().from_terms(
        element
            .terms()
            .iter()
            .filter(|term| term.path().length() == Some(degree))
            .map(|term| (term.path().clone(), term.coefficient().clone())),
    )?)
}

/// Translate list-level `BiggestHomogeneousPartInvLenInvLex`.
///
/// # Errors
///
/// Returns an error unless every element belongs to one path algebra.
pub fn biggest_homogeneous_parts_inv_len_inv_lex<F: Field>(
    elements: &[PathAlgebraElement<F>],
) -> Result<Vec<PathAlgebraElement<F>>, InvLenLexError> {
    if let Some(first) = elements.first()
        && elements
            .iter()
            .skip(1)
            .any(|element| !element.algebra().is_identical(first.algebra()))
    {
        return Err(InvLenLexError::DifferentAlgebras);
    }
    elements
        .iter()
        .map(biggest_homogeneous_part_inv_len_inv_lex)
        .collect()
}

/// Translate QPA's duplicate `LeftmostOccurrenceInvLenInvLex` method.
///
/// QPA returns one-based inclusive endpoint pairs. Rust returns the equivalent
/// zero-based half-open range; an empty pattern is `Some(0..0)` and absence is
/// `None`.
#[must_use]
pub fn leftmost_occurrence_inv_len_inv_lex<T: PartialEq>(
    sequence: &[T],
    pattern: &[T],
) -> Option<std::ops::Range<usize>> {
    leftmost_occurrence(sequence, pattern)
}

/// Translate `RemainderOfDivisionInvLenInvLexGroebnerBasis`.
///
/// # Errors
///
/// Returns an error for foreign/zero divisors or path arithmetic failures.
pub fn remainder_of_division_inv_len_inv_lex<F: Field>(
    mut dividend: PathAlgebraElement<F>,
    divisors: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<PathAlgebraElement<F>, InvLenLexError> {
    ensure_algebra(&dividend, algebra)?;
    for divisor in divisors {
        ensure_algebra(divisor, algebra)?;
        if divisor.is_zero() {
            return Err(InvLenLexError::ZeroDivisor);
        }
    }

    let mut remainder = algebra.zero();
    while !dividend.is_zero() {
        let tip = tip_inv_len_inv_lex(&dividend)?;
        let tip_walk = required_tip_walk(&tip)?;
        let mut reduction_applied = false;
        for divisor in divisors {
            let divisor_walk = required_tip_walk(divisor)?;
            let Some(occurrence) = leftmost_occurrence(&tip_walk, &divisor_walk) else {
                continue;
            };
            let left = context_from_walk(algebra, &tip_walk[..occurrence.start])?;
            let right = context_from_walk(algebra, &tip_walk[occurrence.end..])?;
            let coefficient = algebra.field().div(
                &tip_coefficient_inv_len_inv_lex(&tip),
                &tip_coefficient_inv_len_inv_lex(divisor),
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
            remainder = remainder.add(&tip)?;
            dividend = dividend.sub(&tip)?;
        }
    }
    Ok(remainder)
}

/// Translate `ReducedListInvLenInvLexQPA`.
///
/// # Errors
///
/// Propagates inverse-tip division failures.
pub fn reduced_list_inv_len_inv_lex<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, InvLenLexError> {
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
        let reduced =
            remainder_of_division_inv_len_inv_lex(result[current].clone(), &others, algebra)?;
        if reduced.is_zero() {
            result.remove(current);
        } else {
            result[current] = reduced;
        }
        index -= 1;
    }
    Ok(result)
}

/// Translate `TipReducedListInvLenInvLex`.
///
/// # Errors
///
/// Returns an error for foreign elements or invalid inverse tips.
pub fn tip_reduced_list_inv_len_inv_lex<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, InvLenLexError> {
    ensure_elements(elements, algebra)?;
    let mut result = Vec::new();
    for element in elements.iter().filter(|element| !element.is_zero()) {
        add_set(&mut result, element.clone())?;
    }
    let mut index = result.len();
    while index > 0 {
        let current = index - 1;
        let tip = required_tip_walk(&result[current])?;
        let contained = result[..current].iter().any(|prior| {
            required_tip_walk(prior)
                .is_ok_and(|prior_tip| leftmost_occurrence(&tip, &prior_tip).is_some())
        });
        if contained {
            result.remove(current);
        }
        index -= 1;
    }
    Ok(result)
}

/// Translate QPA `GroebnerBasisInvLenInvLex`.
///
/// The truncation degree is obtained from the finite nontips of the same
/// default GBNP quotient basis that QPA constructs before its inverse-order
/// overlap loop. For an admissible path-algebra quotient this is
/// `LoewyLength(A / relations) + 1`.
///
/// Like QPA's source routine, this computation can fail to terminate for an
/// input without a finite inverse-order basis.
///
/// # Errors
///
/// Returns an error outside the arrow ideal, when QPA's default quotient is
/// not finite dimensional, or for any source-translated arithmetic failure.
#[allow(clippy::too_many_lines)]
pub fn groebner_basis_inv_len_inv_lex<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, InvLenLexError> {
    if !qpa_in_arrow_ideal(elements, algebra)? {
        return Err(InvLenLexError::OutsideArrowIdeal);
    }

    let default_relations = gbnp_groebner_basis(elements, algebra)?;
    let default_basis = GroebnerBasis::from_complete_relations(algebra, &default_relations)?;
    let nontips = default_basis
        .nontips()?
        .ok_or(InvLenLexError::InfiniteDefaultQuotient)?;
    let loewy_length = admissible_loewy_length(&default_basis, &nontips)?;
    let truncation_degree = loewy_length
        .checked_add(1)
        .ok_or(InvLenLexError::DegreeOverflow)?;

    let uniform = make_uniform_list(elements, algebra)?;
    let mut pending = reduced_list_inv_len_inv_lex(&uniform, algebra)?;
    let mut basis = Vec::new();
    while !pending.is_empty() {
        for element in pending {
            let coefficient = tip_coefficient_inv_len_inv_lex(&element);
            basis.push(element.divide_scalar(&coefficient)?);
        }

        let basis_length = basis.len();
        pending = Vec::new();
        for left_index in 0..basis_length {
            let left_tip = required_tip_walk(&basis[left_index])?;
            let left_length = left_tip.len();
            for right_index in 0..basis_length {
                let right_tip = required_tip_walk(&basis[right_index])?;
                let right_length = right_tip.len();
                let first_overlap = left_length.saturating_sub(right_length);
                for split in first_overlap..left_length {
                    let overlap_length = left_length - split;
                    if left_tip[split..] != right_tip[..overlap_length] {
                        continue;
                    }
                    let left_context = context_from_walk(algebra, &left_tip[..split])?;
                    let right_context = context_from_walk(algebra, &right_tip[overlap_length..])?;
                    let overlap = basis[left_index]
                        .multiply(&right_context)?
                        .sub(&left_context.multiply(&basis[right_index])?)?;
                    let remainder =
                        remainder_of_division_inv_len_inv_lex(overlap, &basis, algebra)?;
                    let remainder = truncate_element(&remainder, truncation_degree)?;
                    if !remainder.is_zero() {
                        add_set(&mut pending, remainder)?;
                    }
                }
            }
        }
        basis = tip_reduced_list_inv_len_inv_lex(&basis, algebra)?;
        basis = reduced_list_inv_len_inv_lex(&basis, algebra)?;
    }
    Ok(basis)
}

fn required_tip_walk<F: Field>(
    element: &PathAlgebraElement<F>,
) -> Result<Vec<ArrowIndex>, InvLenLexError> {
    if element.is_zero() {
        return Err(InvLenLexError::ZeroDivisor);
    }
    tip_walk_inv_len_inv_lex(element).ok_or(InvLenLexError::InvalidTipWalk)
}

fn admissible_loewy_length<F: Field>(
    basis: &GroebnerBasis<F>,
    nontips: &[Path],
) -> Result<usize, InvLenLexError> {
    if basis.relations().iter().any(|relation| {
        relation
            .terms()
            .iter()
            .any(|term| term.path().length().is_none_or(|length| length < 2))
    }) {
        return Err(InvLenLexError::NonAdmissibleDefaultQuotient);
    }

    let arrows = basis
        .algebra()
        .quiver()
        .arrows()
        .map(|arrow| basis.algebra().element(arrow.path()))
        .collect::<Result<Vec<_>, _>>()?;
    let initial_rows = nontips
        .iter()
        .enumerate()
        .filter(|(_, path)| path.length().is_some_and(|length| length > 0))
        .map(|(index, _)| {
            let mut row = vec![basis.algebra().field().zero(); nontips.len()];
            row[index] = basis.algebra().field().one();
            row
        })
        .collect::<Vec<_>>();
    let mut power = row_basis(initial_rows, nontips.len(), basis.algebra().field().clone())?;
    let mut exponent = 1_usize;

    if power.nrows() == 0 {
        return Ok(1);
    }
    loop {
        let mut product_rows = Vec::new();
        for row in 0..power.nrows() {
            let representative = basis
                .algebra()
                .from_terms(nontips.iter().cloned().zip(power.row(row).iter().cloned()))?;
            for arrow in &arrows {
                product_rows.push(quotient_coordinates(
                    &representative.multiply(arrow)?,
                    basis,
                    nontips,
                )?);
            }
        }
        let next = row_basis(product_rows, nontips.len(), basis.algebra().field().clone())?;
        exponent = exponent
            .checked_add(1)
            .ok_or(InvLenLexError::DegreeOverflow)?;
        if next.nrows() == 0 {
            return Ok(exponent);
        }
        if next.nrows() == power.nrows() {
            return Err(InvLenLexError::NonAdmissibleDefaultQuotient);
        }
        power = next;
    }
}

fn quotient_coordinates<F: Field>(
    element: &PathAlgebraElement<F>,
    basis: &GroebnerBasis<F>,
    nontips: &[Path],
) -> Result<Vec<F::Element>, InvLenLexError> {
    let normal = basis.remainder(element.clone())?;
    let mut coordinates = vec![basis.algebra().field().zero(); nontips.len()];
    for term in normal.terms() {
        let index = nontips
            .iter()
            .position(|path| path == term.path())
            .ok_or(InvLenLexError::InvalidQuotientCoordinates)?;
        coordinates[index] = term.coefficient().clone();
    }
    Ok(coordinates)
}

fn row_basis<F: Field>(
    rows: Vec<Vec<F::Element>>,
    columns: usize,
    field: F,
) -> Result<DenseMatrix<F>, InvLenLexError> {
    let matrix = if rows.is_empty() {
        DenseMatrix::zero(0, columns, field)
    } else {
        DenseMatrix::from_rows(rows, field)?
    };
    Ok(matrix.semi_echelon_transformation().vectors)
}

fn context_from_walk<F: Field>(
    algebra: &PathAlgebra<F>,
    walk: &[ArrowIndex],
) -> Result<PathAlgebraElement<F>, InvLenLexError> {
    if walk.is_empty() {
        return Ok(algebra.one());
    }
    let arrows = walk
        .iter()
        .map(|&index| {
            algebra
                .quiver()
                .arrow(index)
                .ok_or(InvLenLexError::InvalidTipWalk)
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(algebra.element(Path::from_arrow_list(algebra.quiver(), &arrows)?)?)
}

fn ensure_algebra<F: Field>(
    element: &PathAlgebraElement<F>,
    algebra: &PathAlgebra<F>,
) -> Result<(), InvLenLexError> {
    if element.algebra().is_identical(algebra) {
        Ok(())
    } else {
        Err(InvLenLexError::DifferentAlgebras)
    }
}

fn ensure_elements<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<(), InvLenLexError> {
    for element in elements {
        ensure_algebra(element, algebra)?;
    }
    Ok(())
}

fn add_set<F: Field>(
    elements: &mut Vec<PathAlgebraElement<F>>,
    element: PathAlgebraElement<F>,
) -> Result<(), InvLenLexError> {
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

/// Failure in QPA's inverse-length/inverse-lexicographic routines.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum InvLenLexError {
    /// Relations or divisors come from different path-algebra families.
    #[error("inverse-length computation received different path algebras")]
    DifferentAlgebras,
    /// `GroebnerBasisInvLenInvLex` accepts only arrow-ideal relations.
    #[error("relations do not belong to the arrow ideal")]
    OutsideArrowIdeal,
    /// A zero relation cannot serve as an inverse-tip divisor.
    #[error("inverse-length division received a zero divisor")]
    ZeroDivisor,
    /// An inverse tip did not carry a valid arrow walk.
    #[error("inverse-length tip has no valid arrow walk")]
    InvalidTipWalk,
    /// QPA declares `TruncateElement` only for positive degrees.
    #[error("TruncateElement requires a positive degree")]
    NonPositiveTruncationDegree,
    /// QPA's default quotient must be finite dimensional before it computes a
    /// Loewy truncation bound.
    #[error("the default GBNP quotient has infinitely many nontips")]
    InfiniteDefaultQuotient,
    /// QPA's default quotient does not have nilpotent arrow ideal and hence is
    /// outside its admissible-quotient `LoewyLength` method.
    #[error("the default GBNP quotient is not admissible")]
    NonAdmissibleDefaultQuotient,
    /// A default-basis normal form was not expressible in its nontip basis.
    #[error("default quotient normal form is outside the nontip coordinate basis")]
    InvalidQuotientCoordinates,
    /// A path-length/Loewy truncation degree overflowed `usize`.
    #[error("inverse-length truncation degree overflow")]
    DegreeOverflow,
    /// QPA high-level/object Gröbner machinery failed.
    #[error(transparent)]
    Groebner(#[from] GroebnerError),
    /// QPA's GBNP adapter failed.
    #[error(transparent)]
    QpaGbnp(#[from] QpaGbnpError),
    /// Sparse path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// QPA path construction failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// Exact quotient-coordinate row reduction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring, Set,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, RationalField},
    };
    use qpa_core::{ArrowSpec, PathAlgebra, PathAlgebraElement, Quiver, VertexIndex};
    use serde_json::{Value, json};

    use super::{
        InvLenLexError, biggest_homogeneous_part_inv_len_inv_lex,
        biggest_homogeneous_parts_inv_len_inv_lex, groebner_basis_inv_len_inv_lex,
        leftmost_occurrence_inv_len_inv_lex, reduced_list_inv_len_inv_lex,
        remainder_of_division_inv_len_inv_lex, tip_coefficient_inv_len_inv_lex,
        tip_inv_len_inv_lex, tip_monomial_inv_len_inv_lex, tip_reduced_list_inv_len_inv_lex,
        tip_walk_inv_len_inv_lex, truncate_element,
    };

    struct LoopCase<F: Field> {
        algebra: PathAlgebra<F>,
        vertex: PathAlgebraElement<F>,
        x: PathAlgebraElement<F>,
        x_squared: PathAlgebraElement<F>,
        x_cubed: PathAlgebraElement<F>,
    }

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn loop_case<F: Field>(field: F) -> LoopCase<F> {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let algebra = PathAlgebra::new(field, quiver).unwrap();
        let generators = algebra.generators();
        let vertex = generators[0].clone();
        let x = generators[1].clone();
        let x_squared = x.multiply(&x).unwrap();
        let x_cubed = x_squared.multiply(&x).unwrap();
        LoopCase {
            algebra,
            vertex,
            x,
            x_squared,
            x_cubed,
        }
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the inverse-length fixture uses only small integers")
    }

    fn element_json<F: Field>(
        element: &PathAlgebraElement<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().external_representation(),
                        "coefficient": coefficient_json(term.coefficient())
                    })
                })
                .collect(),
        )
    }

    fn element_list_json<F: Field>(
        elements: &[PathAlgebraElement<F>],
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            elements
                .iter()
                .map(|element| element_json(element, coefficient_json))
                .collect(),
        )
    }

    fn rational_coefficient_json(coefficient: &<RationalField as Set>::Element) -> Value {
        json!([
            small_integer(coefficient.numerator_ref()),
            small_integer(coefficient.denominator_ref())
        ])
    }

    fn rational_multiple(
        element: &PathAlgebraElement<RationalField>,
        coefficient: i64,
    ) -> PathAlgebraElement<RationalField> {
        element.scale_left(&Q.nth(coefficient.into())).unwrap()
    }

    #[test]
    #[allow(clippy::similar_names, clippy::too_many_lines)]
    fn pinned_qpa_invlenlex_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/invlenlex.schema-v1.json"
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

        let rational = loop_case(Q);
        let mixed = rational
            .vertex
            .add(&rational_multiple(&rational.x, 2))
            .unwrap()
            .add(&rational_multiple(&rational.x_squared, 3))
            .unwrap()
            .add(&rational_multiple(&rational.x_cubed, 4))
            .unwrap();
        let positive = rational_multiple(&rational.x, 2)
            .add(&rational_multiple(&rational.x_squared, 3))
            .unwrap();
        let inverse_remainder = remainder_of_division_inv_len_inv_lex(
            rational.x_cubed.add(&rational.x).unwrap(),
            std::slice::from_ref(&rational.x_squared),
            &rational.algebra,
        )
        .unwrap();
        let list_input = [rational.x_squared.clone(), rational.x_cubed.clone()];
        let inverse_reduced = reduced_list_inv_len_inv_lex(&list_input, &rational.algebra).unwrap();
        let inverse_tip_reduced =
            tip_reduced_list_inv_len_inv_lex(&list_input, &rational.algebra).unwrap();
        let inverse_basis = groebner_basis_inv_len_inv_lex(
            std::slice::from_ref(&rational.x_cubed),
            &rational.algebra,
        )
        .unwrap();

        let field5 = Zp::new(5);
        let f5 = loop_case(field5.clone());
        let mixed5 = f5
            .vertex
            .add(&f5.x.scale_left(&field5.to_element(2)).unwrap())
            .unwrap()
            .add(&f5.x_squared.scale_left(&field5.to_element(3)).unwrap())
            .unwrap()
            .add(&f5.x_cubed.scale_left(&field5.to_element(4)).unwrap())
            .unwrap();
        let inverse_basis5 =
            groebner_basis_inv_len_inv_lex(std::slice::from_ref(&f5.x_cubed), &f5.algebra).unwrap();
        let f5_coefficient_json = |coefficient: &<Zp as Set>::Element| {
            small_integer(&field5.to_integer(coefficient)).into()
        };

        let f2 = loop_case(Z2);
        let mixed2 = f2
            .vertex
            .add(&f2.x)
            .unwrap()
            .add(&f2.x_squared)
            .unwrap()
            .add(&f2.x_cubed)
            .unwrap();
        let inverse_basis2 =
            groebner_basis_inv_len_inv_lex(std::slice::from_ref(&f2.x_cubed), &f2.algebra).unwrap();
        let f2_coefficient_json =
            |coefficient: &<Z2 as Set>::Element| small_integer(&Z2.to_integer(coefficient)).into();

        let tip = tip_inv_len_inv_lex(&mixed).unwrap();
        let tip_monomial = tip_monomial_inv_len_inv_lex(&mixed).unwrap();
        let tip_walk = tip_walk_inv_len_inv_lex(&mixed).unwrap();
        let occurrence = leftmost_occurrence_inv_len_inv_lex(&[1, 2, 1, 2], &[1, 2]).unwrap();
        let empty_occurrence = leftmost_occurrence_inv_len_inv_lex(&[1, 2], &[]).unwrap();
        let missing_occurrence = leftmost_occurrence_inv_len_inv_lex(&[1, 2], &[2, 1]);
        let zero = rational.algebra.zero();
        let zero_tip = tip_inv_len_inv_lex(&zero).unwrap();

        let split_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "a"),
                ArrowSpec::named(vertex(2), vertex(2), "b"),
            ],
        )
        .unwrap();
        let split_algebra = PathAlgebra::new(Q, split_quiver).unwrap();
        let split_generators = split_algebra.generators();
        let nonuniform_relation = split_generators[2]
            .multiply(&split_generators[2])
            .unwrap()
            .add(&split_generators[3].multiply(&split_generators[3]).unwrap())
            .unwrap();
        let nonuniform_basis = groebner_basis_inv_len_inv_lex(
            std::slice::from_ref(&nonuniform_relation),
            &split_algebra,
        )
        .unwrap();
        assert_eq!(
            json!({
                "schema_version": 1,
                "kind": "invlenlex",
                "oracle": expected["oracle"].clone(),
                "rational": {
                    "tip": element_json(&tip, &rational_coefficient_json),
                    "tip_coefficient": rational_coefficient_json(
                        &tip_coefficient_inv_len_inv_lex(&mixed)
                    ),
                    "tip_monomial": tip_monomial.external_representation(),
                    "tip_walk": tip_walk
                        .iter()
                        .map(|arrow| {
                            rational
                                .algebra
                                .quiver()
                                .arrow(*arrow)
                                .unwrap()
                                .generator_index()
                        })
                        .collect::<Vec<_>>(),
                    "truncate_1": element_json(
                        &truncate_element(&mixed, 1).unwrap(),
                        &rational_coefficient_json
                    ),
                    "truncate_2": element_json(
                        &truncate_element(&mixed, 2).unwrap(),
                        &rational_coefficient_json
                    ),
                    "biggest_homogeneous": element_json(
                        &biggest_homogeneous_part_inv_len_inv_lex(&positive).unwrap(),
                        &rational_coefficient_json
                    ),
                    "biggest_homogeneous_list": element_list_json(
                        &biggest_homogeneous_parts_inv_len_inv_lex(&[
                            mixed.clone(),
                            positive
                        ]).unwrap(),
                        &rational_coefficient_json
                    ),
                    "remainder": element_json(
                        &inverse_remainder,
                        &rational_coefficient_json
                    ),
                    "reduced_list": element_list_json(
                        &inverse_reduced,
                        &rational_coefficient_json
                    ),
                    "tip_reduced_list": element_list_json(
                        &inverse_tip_reduced,
                        &rational_coefficient_json
                    ),
                    "leftmost_occurrence": [occurrence.start + 1, occurrence.end],
                    "leftmost_empty": [empty_occurrence.start + 1, empty_occurrence.end],
                    "leftmost_missing": missing_occurrence.map(|range| [
                        range.start + 1,
                        range.end
                    ]),
                    "zero": {
                        "tip": element_json(&zero_tip, &rational_coefficient_json),
                        "tip_coefficient": rational_coefficient_json(
                            &tip_coefficient_inv_len_inv_lex(&zero)
                        ),
                        "tip_monomial_is_zero": tip_monomial_inv_len_inv_lex(&zero).is_none(),
                        "tip_walk_success": tip_walk_inv_len_inv_lex(&zero).is_some(),
                        "truncate_zero_success": truncate_element(&mixed, 0).is_ok()
                    },
                    "nonuniform_basis": element_list_json(
                        &nonuniform_basis,
                        &rational_coefficient_json
                    ),
                    "basis": element_list_json(&inverse_basis, &rational_coefficient_json)
                },
                "f5": {
                    "tip": element_json(
                        &tip_inv_len_inv_lex(&mixed5).unwrap(),
                        &f5_coefficient_json
                    ),
                    "truncate_2": element_json(
                        &truncate_element(&mixed5, 2).unwrap(),
                        &f5_coefficient_json
                    ),
                    "basis": element_list_json(&inverse_basis5, &f5_coefficient_json)
                },
                "f2": {
                    "tip": element_json(
                        &tip_inv_len_inv_lex(&mixed2).unwrap(),
                        &f2_coefficient_json
                    ),
                    "truncate_2": element_json(
                        &truncate_element(&mixed2, 2).unwrap(),
                        &f2_coefficient_json
                    ),
                    "basis": element_list_json(&inverse_basis2, &f2_coefficient_json)
                }
            }),
            expected
        );
    }

    #[test]
    fn finite_nontips_do_not_replace_qpa_admissibility() {
        let case = loop_case(Q);
        let nonnilpotent = case.x_squared.sub(&case.x_cubed).unwrap();
        assert_eq!(
            groebner_basis_inv_len_inv_lex(std::slice::from_ref(&nonnilpotent), &case.algebra)
                .unwrap_err(),
            InvLenLexError::NonAdmissibleDefaultQuotient
        );
    }
}
