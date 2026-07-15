//! GBNP basis preprocessing used by `SGrobner`.
//!
//! `reduce_polynomials` and its in-place phase translate `GBNP.ReducePol`,
//! `GBNP.ReducePol2`, `GBNP.StrongNormalForm2TS`, and
//! `GBNP.ReduceCancellation` from the pinned GBNP sources.

use std::collections::BTreeSet;

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    GbnpOptions, NpError, NpMonomial, NpPolynomial, OccurTree, OccurTreeError, ReductionError,
    TreeDirection, add_np, clean_np, mk_monic_np, occur, reduction::reduce_indexed_term,
};

/// Translate public GBNP `ReducePol` using default GBNP options.
///
/// # Errors
///
/// Returns an error for mixed coefficient fields or an occurrence-tree
/// invariant failure.
pub fn reduce_polynomials<F: Field>(
    input: &[NpPolynomial<F>],
) -> Result<Vec<NpPolynomial<F>>, PreprocessingError> {
    reduce_polynomials_with_options(input, &GbnpOptions::default())
}

/// Translate public GBNP `ReducePol` with explicit options.
///
/// # Errors
///
/// Returns an error for mixed coefficient fields or an occurrence-tree
/// invariant failure.
pub fn reduce_polynomials_with_options<F: Field>(
    input: &[NpPolynomial<F>],
    options: &GbnpOptions,
) -> Result<Vec<NpPolynomial<F>>, PreprocessingError> {
    validate_fields(input)?;
    let mut basis = input
        .iter()
        .map(mk_monic_np)
        .filter(|polynomial| !polynomial.is_zero())
        .collect::<Vec<_>>();

    loop {
        for polynomial in &mut basis {
            *polynomial = reduce_cancellation(polynomial, options)?;
        }
        sort_by_leading_monomial(&mut basis);
        let mut done = true;
        for index in (0..basis.len().saturating_sub(1)).rev() {
            if basis[index].leading_monomial() == basis[index + 1].leading_monomial() {
                done = false;
                let one = basis[index].field().one();
                let negative_one = basis[index].field().neg(&one);
                basis[index + 1] = mk_monic_np(&add_np(
                    &basis[index + 1],
                    &basis[index],
                    &one,
                    &negative_one,
                )?);
            }
        }
        if done {
            break;
        }
        basis.retain(|polynomial| !polynomial.is_zero());
        sort_by_leading_monomial(&mut basis);
    }

    reduce_polynomials_in_place(&mut basis, options)?;
    Ok(basis)
}

/// Translate `GBNP.ReducePol2`, mutating a clean monic basis and returning its
/// updated occurrence tree.
///
/// # Errors
///
/// Returns an error for mixed fields, zero basis elements, or a tree invariant
/// failure.
pub fn reduce_polynomials_in_place<F: Field>(
    basis: &mut Vec<NpPolynomial<F>>,
    options: &GbnpOptions,
) -> Result<OccurTree, PreprocessingError> {
    validate_fields(basis)?;
    if basis.iter().any(NpPolynomial::is_zero) {
        return Err(PreprocessingError::ZeroBasisPolynomial);
    }
    for polynomial in basis.iter_mut() {
        *polynomial = reduce_cancellation(polynomial, options)?;
    }
    sort_by_leading_monomial(basis);
    let leading = leading_monomials(basis);
    let prefix_generators = calculate_prefix_generators(&leading).max(options.prefix_generators());
    let mut tree = OccurTree::new(&leading, prefix_generators, TreeDirection::Forward)?;
    let mut pending = (0..basis.len()).collect::<BTreeSet<_>>();
    let mut basis_length = basis.len();

    while let Some(&first_index) = pending.first() {
        pending.remove(&first_index);
        let mut current_first = first_index;
        let mut second_index = current_first + 1;
        while second_index < basis_length {
            let first_tip = basis[current_first]
                .leading_monomial()
                .ok_or(PreprocessingError::ZeroBasisPolynomial)?;
            let second_tip = basis[second_index]
                .leading_monomial()
                .ok_or(PreprocessingError::ZeroBasisPolynomial)?;
            if occur(first_tip, second_tip) == 0 {
                second_index += 1;
                continue;
            }

            let reduced = strong_normal_form_except(basis, second_index, &tree)?;
            let old_tip = basis[second_index]
                .leading_monomial()
                .ok_or(PreprocessingError::ZeroBasisPolynomial)?
                .clone();
            if reduced.is_zero() {
                tree.remove(&old_tip, second_index, TreeDirection::Forward)?;
                basis.remove(second_index);
                pending.remove(&(basis_length - 1));
                basis_length -= 1;
                continue;
            }

            let insertion_index = insertion_position(basis, &reduced);
            tree.remove(&old_tip, second_index, TreeDirection::Forward)?;
            basis.remove(second_index);
            let monic = mk_monic_np(&reduced);
            let new_tip = monic
                .leading_monomial()
                .ok_or(PreprocessingError::ZeroBasisPolynomial)?
                .clone();
            tree.add(&new_tip, insertion_index, TreeDirection::Forward)?;
            basis.insert(insertion_index, monic);
            pending.remove(&second_index);
            pending = pending
                .into_iter()
                .map(|index| {
                    if (insertion_index..second_index).contains(&index) {
                        index + 1
                    } else {
                        index
                    }
                })
                .collect();
            if (insertion_index..second_index).contains(&current_first) {
                current_first += 1;
            }
            pending.insert(insertion_index);
            second_index += 1;
        }
    }
    Ok(tree)
}

/// Translate `GBNP.StrongNormalForm2TS`: reduce one basis element by every
/// other element selected through the occurrence tree.
///
/// # Errors
///
/// Returns an error for an invalid skip index or reduction/tree inconsistency.
pub fn strong_normal_form_except<F: Field>(
    basis: &[NpPolynomial<F>],
    skip_index: usize,
    tree: &OccurTree,
) -> Result<NpPolynomial<F>, PreprocessingError> {
    let Some(mut result) = basis.get(skip_index).cloned() else {
        return Err(PreprocessingError::BasisIndexOutOfRange {
            index: skip_index,
            basis_length: basis.len(),
        });
    };
    let mut term_index = 0;
    while term_index < result.terms().len() {
        let mut monomial = result.terms()[term_index].monomial().clone();
        let mut matches = tree
            .all_occurrences(&monomial, TreeDirection::Forward)
            .into_iter()
            .filter(|(index, _)| *index != skip_index)
            .collect::<Vec<_>>();
        while let Some(&(reducer_index, position)) = matches.first() {
            result = reduce_indexed_term(&result, term_index, &basis[reducer_index], position + 1)?;
            if result.is_zero() || term_index >= result.terms().len() {
                return Ok(result);
            }
            monomial = result.terms()[term_index].monomial().clone();
            matches = tree
                .all_occurrences(&monomial, TreeDirection::Forward)
                .into_iter()
                .filter(|(index, _)| *index != skip_index)
                .collect();
        }
        term_index += 1;
    }
    Ok(result)
}

/// Translate `GBNP.ReduceCancellation`.
///
/// # Errors
///
/// Returns an error only if rebuilding the canonical NP representation fails.
pub fn reduce_cancellation<F: Field>(
    polynomial: &NpPolynomial<F>,
    options: &GbnpOptions,
) -> Result<NpPolynomial<F>, PreprocessingError> {
    if options.cancellative_alphabet_size().is_none() || polynomial.terms().len() < 2 {
        return Ok(polynomial.clone());
    }
    let mut words = polynomial
        .terms()
        .iter()
        .map(|term| term.monomial().generators().to_vec())
        .collect::<Vec<_>>();
    let mut minimum_length = words.iter().map(Vec::len).min().unwrap_or(0);
    let mut same = true;
    while minimum_length > 0 && same {
        let Some(&symbol) = words[0].last() else {
            break;
        };
        if symbol < 0 {
            same = false;
        }
        if same {
            same = words
                .iter()
                .skip(1)
                .all(|word| word.last() == Some(&symbol));
        }
        if same {
            for word in &mut words {
                word.pop();
            }
            minimum_length -= 1;
        }
    }

    same = true;
    while minimum_length > 0 && same {
        let Some(&symbol) = words[0].first() else {
            break;
        };
        if symbol < 0 {
            same = false;
        }
        if same {
            same = words
                .iter()
                .skip(1)
                .all(|word| word.first() == Some(&symbol));
        }
        if same {
            for word in &mut words {
                word.remove(0);
            }
            minimum_length -= 1;
        }
    }

    Ok(clean_np(
        polynomial.field().clone(),
        words.into_iter().map(NpMonomial::from).collect(),
        polynomial
            .terms()
            .iter()
            .map(|term| term.coefficient().clone())
            .collect(),
    )?)
}

/// Translate `GBNP.CalculatePGlts`.
#[must_use]
pub fn calculate_prefix_generators(leading_monomials: &[NpMonomial]) -> usize {
    leading_monomials
        .iter()
        .filter_map(|monomial| monomial.generators().first().copied())
        .filter(|generator| *generator < 0)
        .map(|generator| generator.unsigned_abs() as usize)
        .max()
        .unwrap_or(0)
}

fn leading_monomials<F: Field>(basis: &[NpPolynomial<F>]) -> Vec<NpMonomial> {
    basis
        .iter()
        .filter_map(|polynomial| polynomial.leading_monomial().cloned())
        .collect()
}

fn sort_by_leading_monomial<F: Field>(basis: &mut [NpPolynomial<F>]) {
    basis.sort_by(|left, right| left.leading_monomial().cmp(&right.leading_monomial()));
}

fn insertion_position<F: Field>(basis: &[NpPolynomial<F>], new: &NpPolynomial<F>) -> usize {
    basis.partition_point(|existing| existing.leading_monomial() < new.leading_monomial())
}

fn validate_fields<F: Field>(polynomials: &[NpPolynomial<F>]) -> Result<(), PreprocessingError> {
    if let Some(first) = polynomials.first()
        && polynomials
            .iter()
            .skip(1)
            .any(|polynomial| polynomial.field() != first.field())
    {
        return Err(PreprocessingError::DifferentFields);
    }
    Ok(())
}

/// Failure in GBNP's source-translated basis preprocessing.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PreprocessingError {
    /// All basis polynomials must share a coefficient field.
    #[error("GBNP preprocessing received different coefficient fields")]
    DifferentFields,
    /// `ReducePol2` assumes a list of nonzero polynomials.
    #[error("GBNP ReducePol2 received a zero basis polynomial")]
    ZeroBasisPolynomial,
    /// A basis index is outside the supplied basis.
    #[error("basis index {index} is outside a basis of length {basis_length}")]
    BasisIndexOutOfRange {
        /// Invalid zero-based basis index.
        index: usize,
        /// Number of basis polynomials.
        basis_length: usize,
    },
    /// NP arithmetic failed.
    #[error(transparent)]
    Np(#[from] NpError),
    /// NP normal-form reduction failed.
    #[error(transparent)]
    Reduction(#[from] ReductionError),
    /// Occurrence-tree mutation failed.
    #[error(transparent)]
    OccurTree(#[from] OccurTreeError),
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

    use crate::{GbnpOptions, NpMonomial, NpPolynomial, OccurTree, TreeDirection, clean_np};

    use super::{reduce_cancellation, reduce_polynomials, reduce_polynomials_with_options};

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

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the preprocessing fixture uses only small integers")
    }

    fn polynomial_json<F: Field>(
        polynomial: &NpPolynomial<F>,
        coefficient_json: impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            polynomial
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "monomial": term.monomial().generators(),
                        "coefficient": coefficient_json(term.coefficient())
                    })
                })
                .collect(),
        )
    }

    fn rational_json(polynomial: &NpPolynomial<RationalField>) -> Value {
        polynomial_json(polynomial, |coefficient| {
            json!([
                small_integer(coefficient.numerator_ref()),
                small_integer(coefficient.denominator_ref())
            ])
        })
    }

    fn rational_list_json(polynomials: &[NpPolynomial<RationalField>]) -> Value {
        Value::Array(polynomials.iter().map(rational_json).collect())
    }

    fn zp_list_json(polynomials: &[NpPolynomial<Zp>], field: &Zp) -> Value {
        Value::Array(
            polynomials
                .iter()
                .map(|polynomial| {
                    polynomial_json(polynomial, |coefficient| {
                        small_integer(&field.to_integer(coefficient)).into()
                    })
                })
                .collect(),
        )
    }

    fn z2_list_json(polynomials: &[NpPolynomial<Z2>]) -> Value {
        Value::Array(
            polynomials
                .iter()
                .map(|polynomial| {
                    polynomial_json(polynomial, |coefficient| {
                        small_integer(&Z2.to_integer(coefficient)).into()
                    })
                })
                .collect(),
        )
    }

    fn matches_json(matches: &[(usize, usize)]) -> Value {
        json!(
            matches
                .iter()
                .map(|(index, position)| [index + 1, position + 1])
                .collect::<Vec<_>>()
        )
    }

    #[test]
    fn reduce_pol_resolves_equal_tips_then_internal_divisibility() {
        let reduced = reduce_polynomials(&[
            polynomial(&[&[1, 2], &[1]], &[1, 1]),
            polynomial(&[&[1, 2], &[2]], &[1, 1]),
        ])
        .unwrap();
        assert_eq!(reduced.len(), 2);
        assert_eq!(reduced[0].leading_monomial().unwrap().generators(), [2]);
        assert_eq!(reduced[1].leading_monomial().unwrap().generators(), [1, 1]);
    }

    #[test]
    fn cancellative_option_removes_common_suffix_then_prefix() {
        let input = polynomial(&[&[1, 2, 1], &[1, 1, 1]], &[1, -1]);
        let mut options = GbnpOptions::default();
        options.enable_cancellative_monoid(2);
        let reduced = reduce_cancellation(&input, &options).unwrap();
        assert_eq!(
            reduced
                .terms()
                .iter()
                .map(|term| term.monomial().generators())
                .collect::<Vec<_>>(),
            [&[2][..], &[1]]
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_gbnp_preprocessing_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/gbnp_preprocessing.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            })
        );

        let patterns = [
            NpMonomial::from(&[1, 2][..]),
            NpMonomial::from(&[1][..]),
            NpMonomial::from(&[2, 1][..]),
        ];
        let word = NpMonomial::from(&[2, 1, 2][..]);
        let lookup_tree = OccurTree::new(&patterns, 0, TreeDirection::Forward).unwrap();
        let first_lookup = lookup_tree
            .first_occurrence(&word, TreeDirection::Forward)
            .map(|(index, position)| [index + 1, position + 1])
            .unwrap();
        let all_lookup = lookup_tree.all_occurrences(&word, TreeDirection::Forward);

        let mut mutable_tree = OccurTree::new(
            &[NpMonomial::from(&[1][..]), NpMonomial::from(&[2][..])],
            0,
            TreeDirection::Forward,
        )
        .unwrap();
        let twelve = NpMonomial::from(&[1, 2][..]);
        mutable_tree
            .add(&twelve, 1, TreeDirection::Forward)
            .unwrap();
        let after_insert = mutable_tree.all_occurrences(&twelve, TreeDirection::Forward);
        mutable_tree
            .remove(&twelve, 1, TreeDirection::Forward)
            .unwrap();
        let after_remove = mutable_tree.all_occurrences(&twelve, TreeDirection::Forward);

        let rational_input = [
            polynomial(&[&[1, 2], &[1]], &[1, 1]),
            polynomial(&[&[1, 2], &[2]], &[1, 1]),
        ];
        let rational_reduced = reduce_polynomials(&rational_input).unwrap();
        let duplicate_reduced =
            reduce_polynomials(&[rational_input[0].clone(), rational_input[0].clone()]).unwrap();
        let cancellation_input = polynomial(&[&[1, 2, 1], &[1, 1, 1]], &[1, -1]);
        let mut cancellation_options = GbnpOptions::default();
        cancellation_options.enable_cancellative_monoid(2);
        let cancellation = reduce_cancellation(&cancellation_input, &cancellation_options).unwrap();

        let field5 = Zp::new(5);
        let f5_input = [
            clean_np(
                field5.clone(),
                vec![NpMonomial::from(&[1, 2][..]), NpMonomial::from(&[1][..])],
                vec![field5.one(), field5.one()],
            )
            .unwrap(),
            clean_np(
                field5.clone(),
                vec![NpMonomial::from(&[1, 2][..]), NpMonomial::from(&[2][..])],
                vec![field5.one(), field5.one()],
            )
            .unwrap(),
        ];
        let f5_reduced =
            reduce_polynomials_with_options(&f5_input, &GbnpOptions::default()).unwrap();

        let one2 = Z2.one();
        let f2_input = [
            clean_np(
                Z2,
                vec![NpMonomial::from(&[1, 2][..]), NpMonomial::from(&[1][..])],
                vec![one2, one2],
            )
            .unwrap(),
            clean_np(
                Z2,
                vec![NpMonomial::from(&[1, 2][..]), NpMonomial::from(&[2][..])],
                vec![one2, one2],
            )
            .unwrap(),
        ];
        let f2_reduced = reduce_polynomials(&f2_input).unwrap();

        assert_eq!(
            json!({
                "schema_version": 1,
                "kind": "gbnp_preprocessing",
                "oracle": expected["oracle"].clone(),
                "occurrence_tree": {
                    "first_lookup": first_lookup,
                    "all_lookup": matches_json(&all_lookup),
                    "after_insert": matches_json(&after_insert),
                    "after_remove": matches_json(&after_remove)
                },
                "rational": {
                    "reduced": rational_list_json(&rational_reduced),
                    "duplicate_reduced": rational_list_json(&duplicate_reduced),
                    "cancellation": rational_json(&cancellation)
                },
                "f5": {"reduced": zp_list_json(&f5_reduced, &field5)},
                "f2": {"reduced": z2_list_json(&f2_reduced)}
            }),
            expected
        );
    }
}
