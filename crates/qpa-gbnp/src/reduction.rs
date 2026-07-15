//! GBNP occurrence searches and NP normal-form reduction.
//!
//! The routines translate the corresponding sections of `lib/grobner.gi`:
//! `RightOccur`, `LeftOccur`, `Occur`, `OccurInLst`, `NormalForm2`,
//! `StrongNormalForm2`, and the public `StrongNormalFormNP` wrapper.

use numerica::domains::Field;
use thiserror::Error;

use crate::{NpError, NpMonomial, NpPolynomial, add_np, bimul_np, lt_np, mk_monic_np};

/// Translate `GBNP.RightOccur`; zero means no suffix occurrence and positive
/// results are one-based positions.
#[must_use]
pub fn right_occurrence(pattern: &NpMonomial, word: &NpMonomial) -> usize {
    if pattern == word {
        return 1;
    }
    if !lt_np(pattern, word) {
        return 0;
    }
    let pattern_length = pattern.degree();
    let word_length = word.degree();
    let start = word_length - pattern_length;
    usize::from(word.generators()[start..] == *pattern.generators()) * (start + 1)
}

/// Translate `GBNP.LeftOccur`; one means a prefix occurrence and zero means
/// none.
#[must_use]
pub fn left_occurrence(pattern: &NpMonomial, word: &NpMonomial) -> usize {
    if pattern == word {
        return 1;
    }
    if !lt_np(pattern, word) {
        return 0;
    }
    usize::from(word.generators()[..pattern.degree()] == *pattern.generators())
}

/// Translate `GBNP.Occur`, returning the first one-based start position or
/// zero when the pattern does not occur.
#[must_use]
pub fn occur(pattern: &NpMonomial, word: &NpMonomial) -> usize {
    if pattern == word || pattern.degree() == 0 {
        return 1;
    }
    if !lt_np(pattern, word) {
        return 0;
    }
    word.generators()
        .windows(pattern.degree())
        .position(|window| window == pattern.generators())
        .map_or(0, |position| position + 1)
}

/// Translate `GBNP.OccurInLst`, returning `(one_based_pattern_index,
/// one_based_start_position)` or `(0,0)`.
#[must_use]
pub fn occur_in_list(word: &NpMonomial, patterns: &[NpMonomial]) -> (usize, usize) {
    for (index, pattern) in patterns.iter().enumerate() {
        let position = occur(pattern, word);
        if position > 0 {
            return (index + 1, position);
        }
    }
    (0, 0)
}

/// Translate `GBNP.StrongNormalForm2` using two ordered reducer lists.
///
/// Every term is reduced; the first reducer in `primary` has priority over all
/// reducers in `secondary`, exactly as in the GAP source.
///
/// # Errors
///
/// Returns an error for mixed fields or a zero reducer.
pub fn strong_normal_form_two<F: Field>(
    polynomial: &NpPolynomial<F>,
    primary: &[NpPolynomial<F>],
    secondary: &[NpPolynomial<F>],
) -> Result<NpPolynomial<F>, ReductionError> {
    validate_reducers(polynomial, primary)?;
    validate_reducers(polynomial, secondary)?;
    let primary_tips = leading_monomials(primary);
    let secondary_tips = leading_monomials(secondary);
    let mut result = polynomial.clone();
    let mut term_index = 0;

    while term_index < result.terms().len() {
        let mut monomial = result.terms()[term_index].monomial().clone();
        let mut primary_match = occur_in_list(&monomial, &primary_tips);
        let mut secondary_match = occur_in_list(&monomial, &secondary_tips);
        while primary_match.0 + secondary_match.0 > 0 {
            let (reducer, position) = if primary_match.0 > 0 {
                (&primary[primary_match.0 - 1], primary_match.1)
            } else {
                (&secondary[secondary_match.0 - 1], secondary_match.1)
            };
            result = reduce_indexed_term(&result, term_index, reducer, position)?;
            if result.is_zero() || term_index >= result.terms().len() {
                return Ok(result);
            }
            monomial = result.terms()[term_index].monomial().clone();
            primary_match = occur_in_list(&monomial, &primary_tips);
            secondary_match = occur_in_list(&monomial, &secondary_tips);
        }
        term_index += 1;
    }
    Ok(result)
}

/// Translate GBNP's leading-term-only `NormalForm2` behavior for one reducer
/// list.
///
/// # Errors
///
/// Returns an error for mixed fields or a zero reducer.
pub fn normal_form_np<F: Field>(
    polynomial: &NpPolynomial<F>,
    reducers: &[NpPolynomial<F>],
) -> Result<NpPolynomial<F>, ReductionError> {
    validate_reducers(polynomial, reducers)?;
    let tips = leading_monomials(reducers);
    let mut result = polynomial.clone();
    while let Some(leading) = result.leading_monomial().cloned() {
        let occurrence = occur_in_list(&leading, &tips);
        if occurrence.0 == 0 {
            break;
        }
        result = reduce_indexed_term(&result, 0, &reducers[occurrence.0 - 1], occurrence.1)?;
    }
    Ok(result)
}

/// Translate public GBNP `StrongNormalFormNP`.
///
/// The input and all nonzero reducers are made monic, reducers are sorted by
/// increasing leading monomial, and the original input leading coefficient is
/// restored after strong reduction.
///
/// # Errors
///
/// Returns an error for mixed coefficient fields.
pub fn strong_normal_form_np<F: Field>(
    polynomial: &NpPolynomial<F>,
    reducers: &[NpPolynomial<F>],
) -> Result<NpPolynomial<F>, ReductionError> {
    if polynomial.is_zero() {
        return Ok(polynomial.clone());
    }
    validate_fields(polynomial, reducers)?;
    let Some(leading_coefficient) = polynomial.leading_coefficient().cloned() else {
        return Ok(polynomial.clone());
    };
    let monic = mk_monic_np(polynomial);
    let mut clean_reducers = reducers
        .iter()
        .filter(|reducer| !reducer.is_zero())
        .map(mk_monic_np)
        .collect::<Vec<_>>();
    clean_reducers.sort_by(|left, right| left.leading_monomial().cmp(&right.leading_monomial()));
    let reduced = strong_normal_form_two(&monic, &clean_reducers, &[])?;
    Ok(reduced.scale(&leading_coefficient))
}

pub(crate) fn reduce_indexed_term<F: Field>(
    polynomial: &NpPolynomial<F>,
    term_index: usize,
    reducer: &NpPolynomial<F>,
    one_based_position: usize,
) -> Result<NpPolynomial<F>, ReductionError> {
    let monomial = polynomial.terms()[term_index].monomial();
    let reducer_tip = reducer
        .leading_monomial()
        .ok_or(ReductionError::ZeroReducer)?;
    let start = one_based_position - 1;
    let end = start + reducer_tip.degree();
    let left = &monomial.generators()[..start];
    let right = &monomial.generators()[end..];
    let coefficient = polynomial.field().div(
        polynomial.terms()[term_index].coefficient(),
        reducer
            .leading_coefficient()
            .ok_or(ReductionError::ZeroReducer)?,
    );
    let negative = polynomial.field().neg(&coefficient);
    Ok(add_np(
        polynomial,
        &bimul_np(left, reducer, right),
        &polynomial.field().one(),
        &negative,
    )?)
}

fn leading_monomials<F: Field>(polynomials: &[NpPolynomial<F>]) -> Vec<NpMonomial> {
    polynomials
        .iter()
        .filter_map(|polynomial| polynomial.leading_monomial().cloned())
        .collect()
}

fn validate_reducers<F: Field>(
    polynomial: &NpPolynomial<F>,
    reducers: &[NpPolynomial<F>],
) -> Result<(), ReductionError> {
    validate_fields(polynomial, reducers)?;
    if reducers.iter().any(NpPolynomial::is_zero) {
        Err(ReductionError::ZeroReducer)
    } else {
        Ok(())
    }
}

fn validate_fields<F: Field>(
    polynomial: &NpPolynomial<F>,
    reducers: &[NpPolynomial<F>],
) -> Result<(), ReductionError> {
    if reducers
        .iter()
        .all(|reducer| reducer.field() == polynomial.field())
    {
        Ok(())
    } else {
        Err(ReductionError::DifferentFields)
    }
}

/// Failure in GBNP's source-translated NP reduction routines.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ReductionError {
    /// All NP polynomials in a reduction must share a coefficient field.
    #[error("NP reduction received different coefficient fields")]
    DifferentFields,
    /// Internal `StrongNormalForm2` reducer lists may not contain zero.
    #[error("NP reduction received the zero polynomial as a reducer")]
    ZeroReducer,
    /// NP scaled addition failed.
    #[error(transparent)]
    Np(#[from] NpError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{Ring, rational::Q};

    use crate::{NpMonomial, NpPolynomial, clean_np};

    use super::{left_occurrence, occur, right_occurrence, strong_normal_form_np};

    fn polynomial(
        monomials: &[&[i32]],
        coefficients: &[i64],
    ) -> NpPolynomial<numerica::domains::rational::RationalField> {
        clean_np(
            Q,
            monomials
                .iter()
                .map(|monomial| NpMonomial::from(*monomial))
                .collect(),
            coefficients
                .iter()
                .map(|&coefficient| Q.nth(coefficient.into()))
                .collect(),
        )
        .unwrap()
    }

    #[test]
    fn occurrence_functions_keep_gbnp_one_based_zero_sentinel_semantics() {
        let empty = NpMonomial::from(&[][..]);
        let ab = NpMonomial::from(&[1, 2][..]);
        let bab = NpMonomial::from(&[2, 1, 2][..]);
        assert_eq!(occur(&ab, &bab), 2);
        assert_eq!(left_occurrence(&ab, &bab), 0);
        assert_eq!(right_occurrence(&ab, &bab), 2);
        assert_eq!(occur(&empty, &bab), 1);
    }

    #[test]
    fn strong_normal_form_matches_gbnp_documented_example() {
        let first = polynomial(&[&[1, 1, 2], &[]], &[1, -1]);
        let second = polynomial(&[&[1, 2, 2], &[]], &[1, -1]);
        let input = polynomial(&[&[1, 1, 1, 2], &[2, 1], &[]], &[1, -1, 3]);
        let reduced = strong_normal_form_np(&input, &[first, second]).unwrap();
        assert_eq!(
            reduced
                .terms()
                .iter()
                .map(|term| term.monomial().generators())
                .collect::<Vec<_>>(),
            [&[2, 1][..], &[1], &[]]
        );
        assert_eq!(
            reduced
                .terms()
                .iter()
                .map(|term| term.coefficient().clone())
                .collect::<Vec<_>>(),
            [Q.nth((-1).into()), Q.one(), Q.nth(3.into())]
        );
    }

    #[test]
    fn identity_reducer_reduces_every_polynomial_to_zero() {
        let input = polynomial(&[&[1, 2], &[]], &[2, 3]);
        assert!(
            strong_normal_form_np(&input, &[NpPolynomial::one(Q)])
                .unwrap()
                .is_zero()
        );
        assert_eq!(
            strong_normal_form_np(&input, &[NpPolynomial::zero(Q)]).unwrap(),
            input
        );
    }
}
