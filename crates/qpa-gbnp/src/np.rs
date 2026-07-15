//! GBNP's NP noncommutative-polynomial representation and arithmetic.
//!
//! This module translates `lib/nparith.gi`. NP monomials are integer lists;
//! positive entries denote algebra generators and negative leading entries
//! denote free-module generators. Clean NP polynomials store distinct terms in
//! descending degree-lexicographic order.

use std::cmp::Ordering;

use numerica::domains::Field;
use thiserror::Error;

/// A GBNP NP monomial represented by its generator-number word.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct NpMonomial(Box<[i32]>);

impl NpMonomial {
    /// Construct an NP monomial from GBNP generator numbers.
    #[must_use]
    pub fn new(generators: impl Into<Box<[i32]>>) -> Self {
        Self(generators.into())
    }

    /// Return the generator word. The empty word is the identity monomial.
    #[must_use]
    pub fn generators(&self) -> &[i32] {
        &self.0
    }

    /// Return the total NP degree used as the first ordering key.
    #[must_use]
    pub fn degree(&self) -> usize {
        self.0.len()
    }

    fn concatenate(left: &[i32], middle: &Self, right: &[i32]) -> Self {
        let mut generators = Vec::with_capacity(left.len() + middle.degree() + right.len());
        generators.extend_from_slice(left);
        generators.extend_from_slice(middle.generators());
        generators.extend_from_slice(right);
        Self(generators.into_boxed_slice())
    }
}

impl From<Vec<i32>> for NpMonomial {
    fn from(value: Vec<i32>) -> Self {
        Self::new(value.into_boxed_slice())
    }
}

impl From<&[i32]> for NpMonomial {
    fn from(value: &[i32]) -> Self {
        Self::new(value.to_vec().into_boxed_slice())
    }
}

impl Ord for NpMonomial {
    fn cmp(&self, other: &Self) -> Ordering {
        self.degree()
            .cmp(&other.degree())
            .then_with(|| self.generators().cmp(other.generators()))
    }
}

impl PartialOrd for NpMonomial {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Translate GBNP's `LtNP` degree-then-lexicographic comparison.
#[must_use]
pub fn lt_np(left: &NpMonomial, right: &NpMonomial) -> bool {
    left < right
}

/// Translate GBNP's `GtNP` degree-then-lexicographic comparison.
#[must_use]
pub fn gt_np(left: &NpMonomial, right: &NpMonomial) -> bool {
    left > right
}

/// One nonzero term of a clean NP polynomial.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NpTerm<F: Field> {
    monomial: NpMonomial,
    coefficient: F::Element,
}

impl<F: Field> NpTerm<F> {
    /// Return the NP monomial.
    #[must_use]
    pub const fn monomial(&self) -> &NpMonomial {
        &self.monomial
    }

    /// Return its nonzero coefficient.
    #[must_use]
    pub const fn coefficient(&self) -> &F::Element {
        &self.coefficient
    }
}

/// A clean polynomial in GBNP's NP format.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NpPolynomial<F: Field> {
    field: F,
    terms: Vec<NpTerm<F>>,
}

impl<F: Field> NpPolynomial<F> {
    /// Construct the NP zero polynomial `[[],[]]`.
    #[must_use]
    pub fn zero(field: F) -> Self {
        Self {
            field,
            terms: Vec::new(),
        }
    }

    /// Construct the NP identity polynomial `[[[]],[1]]`.
    #[must_use]
    pub fn one(field: F) -> Self {
        let coefficient = field.one();
        Self {
            field,
            terms: vec![NpTerm {
                monomial: NpMonomial::new(Box::<[i32]>::default()),
                coefficient,
            }],
        }
    }

    /// Return the Numerica coefficient field.
    #[must_use]
    pub const fn field(&self) -> &F {
        &self.field
    }

    /// Return canonical terms, greatest NP monomial first.
    #[must_use]
    pub fn terms(&self) -> &[NpTerm<F>] {
        &self.terms
    }

    /// Return whether this is GBNP's NP zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Return the leading monomial, corresponding to `LMonNP`.
    #[must_use]
    pub fn leading_monomial(&self) -> Option<&NpMonomial> {
        self.terms.first().map(NpTerm::monomial)
    }

    /// Return the leading coefficient.
    #[must_use]
    pub fn leading_coefficient(&self) -> Option<&F::Element> {
        self.terms.first().map(NpTerm::coefficient)
    }

    /// Return the minimum number of positive algebra generators required by
    /// `NumAlgGensNP`.
    #[must_use]
    pub fn number_of_algebra_generators(&self) -> usize {
        self.terms
            .iter()
            .flat_map(|term| term.monomial.generators())
            .copied()
            .filter_map(|generator| usize::try_from(generator).ok())
            .max()
            .unwrap_or(0)
    }

    /// Return the minimum number of negative module generators required by
    /// `NumModGensNP`.
    #[must_use]
    pub fn number_of_module_generators(&self) -> usize {
        self.terms
            .iter()
            .flat_map(|term| term.monomial.generators())
            .copied()
            .filter(|generator| *generator < 0)
            .filter_map(|generator| usize::try_from(generator.unsigned_abs()).ok())
            .max()
            .unwrap_or(0)
    }

    /// Translate `GBNP.IsNPMNotZero`.
    #[must_use]
    pub fn is_nonzero_module_polynomial(&self) -> bool {
        self.terms
            .first()
            .and_then(|term| term.monomial.generators().first())
            .is_some_and(|generator| *generator < 0)
    }

    fn from_canonical(field: F, terms: Vec<NpTerm<F>>) -> Self {
        debug_assert!(
            terms
                .windows(2)
                .all(|pair| pair[0].monomial > pair[1].monomial)
        );
        debug_assert!(terms.iter().all(|term| !field.is_zero(&term.coefficient)));
        Self { field, terms }
    }

    pub(crate) fn scale(&self, scalar: &F::Element) -> Self {
        if self.field.is_zero(scalar) || self.is_zero() {
            return Self::zero(self.field.clone());
        }
        let terms = self
            .terms
            .iter()
            .map(|term| NpTerm {
                monomial: term.monomial.clone(),
                coefficient: self.field.mul(&term.coefficient, scalar),
            })
            .collect();
        Self::from_canonical(self.field.clone(), terms)
    }
}

/// Translate GBNP's `CleanNP`: sort greatest-first, combine equal monomials,
/// and remove zero coefficients.
///
/// # Errors
///
/// Returns an error if the two raw NP lists have different lengths.
pub fn clean_np<F: Field>(
    field: F,
    monomials: Vec<NpMonomial>,
    coefficients: Vec<F::Element>,
) -> Result<NpPolynomial<F>, NpError> {
    if monomials.len() != coefficients.len() {
        return Err(NpError::MismatchedRawLists {
            monomials: monomials.len(),
            coefficients: coefficients.len(),
        });
    }
    let mut raw = monomials.into_iter().zip(coefficients).collect::<Vec<_>>();
    raw.sort_by(|left, right| right.0.cmp(&left.0));

    let mut terms: Vec<NpTerm<F>> = Vec::new();
    for (monomial, coefficient) in raw {
        if let Some(previous) = terms.last_mut()
            && previous.monomial == monomial
        {
            field.add_assign(&mut previous.coefficient, &coefficient);
            if field.is_zero(&previous.coefficient) {
                terms.pop();
            }
        } else if !field.is_zero(&coefficient) {
            terms.push(NpTerm {
                monomial,
                coefficient,
            });
        }
    }
    Ok(NpPolynomial::from_canonical(field, terms))
}

/// Translate GBNP's `LMonNP`; `None` corresponds to GAP `fail` for zero.
#[must_use]
pub fn lmon_np<F: Field>(polynomial: &NpPolynomial<F>) -> Option<&NpMonomial> {
    polynomial.leading_monomial()
}

/// Translate GBNP's `MkMonicNP`.
#[must_use]
pub fn mk_monic_np<F: Field>(polynomial: &NpPolynomial<F>) -> NpPolynomial<F> {
    let Some(leading) = polynomial.leading_coefficient() else {
        return polynomial.clone();
    };
    if polynomial.field.is_one(leading) {
        return polynomial.clone();
    }
    let terms = polynomial
        .terms
        .iter()
        .map(|term| NpTerm {
            monomial: term.monomial.clone(),
            coefficient: polynomial.field.div(&term.coefficient, leading),
        })
        .collect();
    NpPolynomial::from_canonical(polynomial.field.clone(), terms)
}

/// Translate GBNP's `AddNP`, returning `c*left + d*right`.
///
/// # Errors
///
/// Returns an error unless both NP polynomials use equal Numerica fields.
pub fn add_np<F: Field>(
    left: &NpPolynomial<F>,
    right: &NpPolynomial<F>,
    left_scalar: &F::Element,
    right_scalar: &F::Element,
) -> Result<NpPolynomial<F>, NpError> {
    ensure_same_field(left, right)?;
    let field = &left.field;
    let mut left_position = 0;
    let mut right_position = 0;
    let mut terms = Vec::with_capacity(left.terms.len() + right.terms.len());

    while left_position < left.terms.len() || right_position < right.terms.len() {
        let (monomial, coefficient) = if left_position >= left.terms.len() {
            let term = &right.terms[right_position];
            right_position += 1;
            (
                term.monomial.clone(),
                field.mul(&term.coefficient, right_scalar),
            )
        } else if right_position >= right.terms.len() {
            let term = &left.terms[left_position];
            left_position += 1;
            (
                term.monomial.clone(),
                field.mul(&term.coefficient, left_scalar),
            )
        } else {
            let left_term = &left.terms[left_position];
            let right_term = &right.terms[right_position];
            match left_term.monomial.cmp(&right_term.monomial) {
                Ordering::Greater => {
                    left_position += 1;
                    (
                        left_term.monomial.clone(),
                        field.mul(&left_term.coefficient, left_scalar),
                    )
                }
                Ordering::Less => {
                    right_position += 1;
                    (
                        right_term.monomial.clone(),
                        field.mul(&right_term.coefficient, right_scalar),
                    )
                }
                Ordering::Equal => {
                    left_position += 1;
                    right_position += 1;
                    (
                        left_term.monomial.clone(),
                        field.add(
                            field.mul(&left_term.coefficient, left_scalar),
                            field.mul(&right_term.coefficient, right_scalar),
                        ),
                    )
                }
            }
        };
        if !field.is_zero(&coefficient) {
            terms.push(NpTerm {
                monomial,
                coefficient,
            });
        }
    }
    Ok(NpPolynomial::from_canonical(left.field.clone(), terms))
}

/// Translate GBNP's `BimulNP`, multiplying every monomial on the left and
/// right by fixed words.
#[must_use]
pub fn bimul_np<F: Field>(
    left_word: &[i32],
    polynomial: &NpPolynomial<F>,
    right_word: &[i32],
) -> NpPolynomial<F> {
    let terms = polynomial
        .terms
        .iter()
        .map(|term| NpTerm {
            monomial: NpMonomial::concatenate(left_word, &term.monomial, right_word),
            coefficient: term.coefficient.clone(),
        })
        .collect();
    NpPolynomial::from_canonical(polynomial.field.clone(), terms)
}

/// Translate GBNP's `MulNP`, followed by `CleanNP`.
///
/// # Errors
///
/// Returns an error unless both NP polynomials use equal Numerica fields.
pub fn mul_np<F: Field>(
    left: &NpPolynomial<F>,
    right: &NpPolynomial<F>,
) -> Result<NpPolynomial<F>, NpError> {
    ensure_same_field(left, right)?;
    let mut monomials = Vec::with_capacity(left.terms.len().saturating_mul(right.terms.len()));
    let mut coefficients = Vec::with_capacity(monomials.capacity());
    for left_term in &left.terms {
        for right_term in &right.terms {
            monomials.push(NpMonomial::concatenate(
                left_term.monomial.generators(),
                &right_term.monomial,
                &[],
            ));
            coefficients.push(
                left.field
                    .mul(&left_term.coefficient, &right_term.coefficient),
            );
        }
    }
    clean_np(left.field.clone(), monomials, coefficients)
}

fn ensure_same_field<F: Field>(
    left: &NpPolynomial<F>,
    right: &NpPolynomial<F>,
) -> Result<(), NpError> {
    if left.field == right.field {
        Ok(())
    } else {
        Err(NpError::DifferentFields)
    }
}

/// Failure in GBNP's source-translated NP representation.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum NpError {
    /// The raw monomial and coefficient lists must be parallel.
    #[error(
        "NP monomial and coefficient lists have different lengths ({monomials} and {coefficients})"
    )]
    MismatchedRawLists {
        /// Number of monomials.
        monomials: usize,
        /// Number of coefficients.
        coefficients: usize,
    },
    /// Arithmetic requires the same Numerica coefficient field.
    #[error("NP polynomials use different coefficient fields")]
    DifferentFields,
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

    use crate::{left_occurrence, occur, occur_in_list, right_occurrence, strong_normal_form_np};

    use super::{NpMonomial, NpPolynomial, add_np, bimul_np, clean_np, lt_np, mk_monic_np, mul_np};

    fn q(numerator: i64, denominator: i64) -> Rational {
        Q.to_element(numerator.into(), denominator.into(), true)
    }

    fn polynomial(
        monomials: &[&[i32]],
        coefficients: &[(i64, i64)],
    ) -> NpPolynomial<numerica::domains::rational::RationalField> {
        clean_np(
            Q,
            monomials
                .iter()
                .map(|monomial| NpMonomial::from(*monomial))
                .collect(),
            coefficients
                .iter()
                .map(|&(numerator, denominator)| q(numerator, denominator))
                .collect(),
        )
        .unwrap()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the NP fixture uses only small integers")
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

    fn zp_json(polynomial: &NpPolynomial<Zp>, field: &Zp) -> Value {
        polynomial_json(polynomial, |coefficient| {
            small_integer(&field.to_integer(coefficient)).into()
        })
    }

    fn z2_json(polynomial: &NpPolynomial<Z2>) -> Value {
        polynomial_json(polynomial, |coefficient| {
            small_integer(&Z2.to_integer(coefficient)).into()
        })
    }

    #[test]
    fn degree_lexicographic_order_is_gbnp_lt_np() {
        assert!(lt_np(
            &NpMonomial::from(&[2_i32][..]),
            &NpMonomial::from(&[1_i32, 1][..])
        ));
        assert!(lt_np(
            &NpMonomial::from(&[1_i32, 1, 2][..]),
            &NpMonomial::from(&[2_i32, 2, 1][..])
        ));
    }

    #[test]
    fn clean_add_bimultiply_and_multiply_match_gbnp_examples() {
        let dirty = polynomial(
            &[&[1, 1, 2], &[], &[1, 1, 2], &[]],
            &[(1, 1), (-3, 1), (-2, 1), (3, 1)],
        );
        assert_eq!(dirty.terms().len(), 1);
        assert_eq!(dirty.terms()[0].coefficient(), &q(-1, 1));

        let first = polynomial(&[&[1, 1, 2], &[]], &[(1, 1), (-3, 1)]);
        let second = polynomial(&[&[1, 2, 2], &[]], &[(1, 1), (-4, 1)]);
        let sum = add_np(&first, &second, &q(1, 1), &q(-1, 1)).unwrap();
        assert_eq!(
            sum.terms()
                .iter()
                .map(|term| term.monomial().generators())
                .collect::<Vec<_>>(),
            [&[1, 2, 2][..], &[1, 1, 2], &[]]
        );
        assert_eq!(
            bimul_np(&[], &first, &[2]).terms()[0]
                .monomial()
                .generators(),
            [1, 1, 2, 2]
        );
        assert_eq!(mul_np(&first, &second).unwrap().terms().len(), 4);
    }

    #[test]
    fn monic_and_prime_field_cancellation_use_numerica_fields() {
        let polynomial = polynomial(&[&[1, 1, 2], &[]], &[(2, 1), (-1, 1)]);
        let monic = mk_monic_np(&polynomial);
        assert_eq!(monic.leading_coefficient(), Some(&q(1, 1)));
        assert_eq!(monic.terms()[1].coefficient(), &q(-1, 2));

        let field5 = Zp::new(5);
        let p5 = clean_np(
            field5.clone(),
            vec![NpMonomial::from(&[1_i32][..]); 2],
            vec![field5.to_element(2), field5.to_element(3)],
        )
        .unwrap();
        assert!(p5.is_zero());

        let one2 = Z2.one();
        let p2 = clean_np(
            Z2,
            vec![NpMonomial::from(&[1_i32][..]); 2],
            vec![one2, one2],
        )
        .unwrap();
        assert!(p2.is_zero());
    }

    #[test]
    fn generator_counts_include_npm_negative_generators() {
        let module = polynomial(&[&[-2, 1], &[-3, 2]], &[(1, 1), (-1, 1)]);
        assert_eq!(module.number_of_algebra_generators(), 2);
        assert_eq!(module.number_of_module_generators(), 3);
        assert!(module.is_nonzero_module_polynomial());
        assert_eq!(NpPolynomial::zero(Q).number_of_algebra_generators(), 0);
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_gbnp_np_arithmetic_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/gbnp_np_arithmetic.schema-v1.json"
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

        let first = polynomial(&[&[1, 1, 2], &[]], &[(1, 1), (-3, 1)]);
        let second = polynomial(&[&[1, 2, 2], &[]], &[(1, 1), (-4, 1)]);
        let dirty = polynomial(
            &[&[1, 1, 2], &[], &[1, 1, 2], &[]],
            &[(1, 1), (-3, 1), (-2, 1), (3, 1)],
        );
        let product = mul_np(&first, &second).unwrap();
        let reverse_product = mul_np(&second, &first).unwrap();
        let commutator = add_np(&product, &reverse_product, &q(1, 1), &q(-1, 1)).unwrap();
        let monic = mk_monic_np(&polynomial(&[&[1, 1, 2], &[]], &[(2, 1), (-1, 1)]));
        let leading_term = clean_np(
            Q,
            vec![first.leading_monomial().unwrap().clone()],
            vec![first.leading_coefficient().unwrap().clone()],
        )
        .unwrap();
        let module = polynomial(&[&[-2, 1], &[-3, 2]], &[(1, 1), (-1, 1)]);
        let normal_input = polynomial(&[&[1, 1, 1, 2], &[2, 1], &[]], &[(1, 1), (-1, 1), (3, 1)]);
        let strong_normal =
            strong_normal_form_np(&normal_input, &[first.clone(), second.clone()]).unwrap();
        let ab = NpMonomial::from(&[1, 2][..]);
        let b = NpMonomial::from(&[2][..]);
        let bab = NpMonomial::from(&[2, 1, 2][..]);
        let empty = NpMonomial::from(&[][..]);
        assert_eq!(
            json!({
                "clean": rational_json(&dirty),
                "add_p1_minus_p2": rational_json(&add_np(&first, &second, &q(1, 1), &q(-1, 1)).unwrap()),
                "bimul_right_2": rational_json(&bimul_np(&[], &first, &[2])),
                "bimul_left_1": rational_json(&bimul_np(&[1], &second, &[])),
                "multiply": rational_json(&product),
                "commutator": rational_json(&commutator),
                "monic": rational_json(&monic),
                "leading_monomial": first.leading_monomial().unwrap().generators(),
                "leading_term": rational_json(&leading_term),
                "generator_counts": {
                    "algebra": module.number_of_algebra_generators(),
                    "module": module.number_of_module_generators()
                },
                "occurrences": {
                    "occur_ab_in_bab": occur(&ab, &bab),
                    "left_ab_in_bab": left_occurrence(&ab, &bab),
                    "right_ab_in_bab": right_occurrence(&ab, &bab),
                    "empty_in_bab": occur(&empty, &bab),
                    "first_in_list": occur_in_list(&bab, &[ab, b])
                },
                "strong_normal_form": rational_json(&strong_normal)
            }),
            expected["rational"]
        );

        let field5 = Zp::new(5);
        let f5_cancelled = clean_np(
            field5.clone(),
            vec![NpMonomial::from(&[1_i32][..]); 2],
            vec![field5.to_element(2), field5.to_element(3)],
        )
        .unwrap();
        let f5_monic = mk_monic_np(
            &clean_np(
                field5.clone(),
                vec![
                    NpMonomial::from(&[2_i32, 1][..]),
                    NpMonomial::from(&[1_i32][..]),
                ],
                vec![field5.to_element(2), field5.to_element(4)],
            )
            .unwrap(),
        );
        let f5_reducer = clean_np(
            field5.clone(),
            vec![
                NpMonomial::from(&[2_i32, 1][..]),
                NpMonomial::from(&[1_i32, 2][..]),
            ],
            vec![field5.one(), field5.to_element(2)],
        )
        .unwrap();
        let f5_input = clean_np(
            field5.clone(),
            vec![NpMonomial::from(&[2_i32, 1][..])],
            vec![field5.one()],
        )
        .unwrap();
        let f5_normal = strong_normal_form_np(&f5_input, &[f5_reducer]).unwrap();
        assert_eq!(
            json!({
                "cancelled": zp_json(&f5_cancelled, &field5),
                "monic": zp_json(&f5_monic, &field5),
                "strong_normal_form": zp_json(&f5_normal, &field5)
            }),
            expected["f5"]
        );

        let one2 = Z2.one();
        let f2_cancelled = clean_np(
            Z2,
            vec![NpMonomial::from(&[1_i32][..]); 2],
            vec![one2, one2],
        )
        .unwrap();
        let f2_clean = clean_np(
            Z2,
            vec![
                NpMonomial::from(&[2_i32, 1][..]),
                NpMonomial::from(&[][..]),
                NpMonomial::from(&[2_i32, 1][..]),
            ],
            vec![one2, one2, one2],
        )
        .unwrap();
        let f2_reducer = clean_np(
            Z2,
            vec![
                NpMonomial::from(&[2_i32, 1][..]),
                NpMonomial::from(&[1_i32, 2][..]),
            ],
            vec![one2, one2],
        )
        .unwrap();
        let f2_input = clean_np(Z2, vec![NpMonomial::from(&[2_i32, 1][..])], vec![one2]).unwrap();
        let f2_normal = strong_normal_form_np(&f2_input, &[f2_reducer]).unwrap();
        assert_eq!(
            json!({
                "cancelled": z2_json(&f2_cancelled),
                "clean": z2_json(&f2_clean),
                "strong_normal_form": z2_json(&f2_normal)
            }),
            expected["f2"]
        );
    }
}
