//! GBNP overlap obstructions and S-polynomials.
//!
//! This module translates `GBNP.SelfObs`, `GBNP.LeftObs`,
//! `GBNP.RightObs`, and `GBNP.Spoly` from `lib/grobner.gi` in the pinned
//! GBNP source. Basis indices are zero-based in the Rust API; GAP's serialized
//! obstruction records use one-based indices.

use std::cmp::Ordering;

use numerica::domains::Field;
use thiserror::Error;

use crate::{NpError, NpMonomial, NpPolynomial, add_np, bimul_np, left_occurrence};

/// GBNP's six-part overlap record
/// `[left_1, index_1, right_1, left_2, index_2, right_2]`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Obstruction {
    left_first: NpMonomial,
    first_basis_index: usize,
    right_first: NpMonomial,
    left_second: NpMonomial,
    second_basis_index: usize,
    right_second: NpMonomial,
}

impl Obstruction {
    /// Construct an obstruction using zero-based Rust basis indices.
    #[must_use]
    pub fn new(
        left_first: NpMonomial,
        first_basis_index: usize,
        right_first: NpMonomial,
        left_second: NpMonomial,
        second_basis_index: usize,
        right_second: NpMonomial,
    ) -> Self {
        Self {
            left_first,
            first_basis_index,
            right_first,
            left_second,
            second_basis_index,
            right_second,
        }
    }

    /// Word multiplying the first basis polynomial on the left.
    #[must_use]
    pub const fn left_first(&self) -> &NpMonomial {
        &self.left_first
    }

    /// Zero-based index of the first basis polynomial.
    #[must_use]
    pub const fn first_basis_index(&self) -> usize {
        self.first_basis_index
    }

    /// Word multiplying the first basis polynomial on the right.
    #[must_use]
    pub const fn right_first(&self) -> &NpMonomial {
        &self.right_first
    }

    /// Word multiplying the second basis polynomial on the left.
    #[must_use]
    pub const fn left_second(&self) -> &NpMonomial {
        &self.left_second
    }

    /// Zero-based index of the second basis polynomial.
    #[must_use]
    pub const fn second_basis_index(&self) -> usize {
        self.second_basis_index
    }

    /// Word multiplying the second basis polynomial on the right.
    #[must_use]
    pub const fn right_second(&self) -> &NpMonomial {
        &self.right_second
    }

    fn reversed_sides(&self) -> Self {
        Self::new(
            self.left_second.clone(),
            self.second_basis_index,
            self.right_second.clone(),
            self.left_first.clone(),
            self.first_basis_index,
            self.right_first.clone(),
        )
    }
}

/// Translate `GBNP.SelfObs`, returning its single smallest self-obstruction.
///
/// # Errors
///
/// Returns an error if `basis_index` is outside `leading_monomials`.
pub fn self_obstruction(
    basis_index: usize,
    leading_monomials: &[NpMonomial],
) -> Result<Option<Obstruction>, ObstructionError> {
    let word =
        leading_monomials
            .get(basis_index)
            .ok_or(ObstructionError::BasisIndexOutOfRange {
                index: basis_index,
                basis_length: leading_monomials.len(),
            })?;
    for shift in 1..word.degree() {
        let overlap_length = word.degree() - shift;
        if word.generators()[..overlap_length] == word.generators()[shift..] {
            return Ok(Some(Obstruction::new(
                empty_monomial(),
                basis_index,
                NpMonomial::from(&word.generators()[overlap_length..]),
                NpMonomial::from(&word.generators()[..shift]),
                basis_index,
                empty_monomial(),
            )));
        }
    }
    Ok(None)
}

/// Translate `GBNP.LeftObs`, including the supplied self-obstruction.
///
/// # Errors
///
/// Returns an error if `basis_index` is outside `leading_monomials`.
pub fn left_obstructions(
    basis_index: usize,
    leading_monomials: &[NpMonomial],
    self_obstruction: Option<&Obstruction>,
) -> Result<Vec<Obstruction>, ObstructionError> {
    let word =
        leading_monomials
            .get(basis_index)
            .ok_or(ObstructionError::BasisIndexOutOfRange {
                index: basis_index,
                basis_length: leading_monomials.len(),
            })?;
    let mut answer = self_obstruction
        .iter()
        .copied()
        .cloned()
        .collect::<Vec<_>>();

    for (other_index, other) in leading_monomials[..basis_index].iter().enumerate() {
        let minimum_length = word.degree().min(other.degree());
        for overlap_length in 1..minimum_length {
            if word.generators()[word.degree() - overlap_length..]
                != other.generators()[..overlap_length]
            {
                continue;
            }
            let left_second =
                NpMonomial::from(&word.generators()[..word.degree() - overlap_length]);
            let right_first = NpMonomial::from(&other.generators()[overlap_length..]);
            let mut keep_new = true;
            let mut eliminate = vec![false; answer.len()];
            for (position, existing) in answer.iter().enumerate() {
                if left_occurrence(existing.right_first(), &right_first) == 1 {
                    keep_new = false;
                    break;
                }
                if left_occurrence(&right_first, existing.right_first()) == 1 {
                    eliminate[position] = true;
                }
            }
            let mut position = 0;
            answer.retain(|_| {
                let retain = !eliminate[position];
                position += 1;
                retain
            });
            sort_as_gap_set(&mut answer);
            if keep_new {
                answer.push(Obstruction::new(
                    empty_monomial(),
                    basis_index,
                    right_first,
                    left_second,
                    other_index,
                    empty_monomial(),
                ));
            }
        }
    }
    Ok(answer)
}

/// Translate `GBNP.RightObs`, including the source's reversed self-obstruction
/// behavior.
///
/// The pinned GBNP implementation uses `Append(eli,h)` and
/// `Difference(ans,sobr)` where `h` and `sobr` are obstruction records rather
/// than singleton lists. For valid NP words those flattened values cannot
/// equal a complete obstruction, so neither expression removes an element.
/// This translation deliberately preserves that observable behavior.
///
/// # Errors
///
/// Returns an error if `basis_index` is outside `leading_monomials`.
pub fn right_obstructions(
    basis_index: usize,
    leading_monomials: &[NpMonomial],
    self_obstruction: Option<&Obstruction>,
) -> Result<Vec<Obstruction>, ObstructionError> {
    let word =
        leading_monomials
            .get(basis_index)
            .ok_or(ObstructionError::BasisIndexOutOfRange {
                index: basis_index,
                basis_length: leading_monomials.len(),
            })?;
    let mut answer = self_obstruction
        .iter()
        .copied()
        .map(Obstruction::reversed_sides)
        .collect::<Vec<_>>();

    for (other_index, other) in leading_monomials[..basis_index].iter().enumerate() {
        let minimum_length = word.degree().min(other.degree());
        for overlap_length in 1..minimum_length {
            if other.generators()[other.degree() - overlap_length..]
                != word.generators()[..overlap_length]
            {
                continue;
            }
            let left_first =
                NpMonomial::from(&other.generators()[..other.degree() - overlap_length]);
            let right_second = NpMonomial::from(&word.generators()[overlap_length..]);
            let mut keep_new = true;
            for existing in &answer {
                if reversed_left_occurrence(existing.left_first(), &left_first) == 1 {
                    keep_new = false;
                    break;
                }
                if reversed_left_occurrence(&left_first, existing.left_first()) == 1 {
                    // Preserve the pinned source's flattened `Append(eli,h)`:
                    // it does not remove this well-formed obstruction.
                }
            }
            if keep_new {
                answer.push(Obstruction::new(
                    left_first,
                    basis_index,
                    empty_monomial(),
                    empty_monomial(),
                    other_index,
                    right_second,
                ));
            }
            sort_as_gap_set(&mut answer);
        }
    }
    sort_as_gap_set(&mut answer);
    Ok(answer)
}

/// Translate `GBNP.Spoly` for one obstruction and polynomial basis.
///
/// # Errors
///
/// Returns an error for an invalid basis index or mixed coefficient fields.
pub fn spolynomial<F: Field>(
    obstruction: &Obstruction,
    basis: &[NpPolynomial<F>],
) -> Result<NpPolynomial<F>, ObstructionError> {
    let first = basis.get(obstruction.first_basis_index()).ok_or(
        ObstructionError::BasisIndexOutOfRange {
            index: obstruction.first_basis_index(),
            basis_length: basis.len(),
        },
    )?;
    let second = basis.get(obstruction.second_basis_index()).ok_or(
        ObstructionError::BasisIndexOutOfRange {
            index: obstruction.second_basis_index(),
            basis_length: basis.len(),
        },
    )?;
    let first_multiple = bimul_np(
        obstruction.left_first().generators(),
        first,
        obstruction.right_first().generators(),
    );
    let second_multiple = bimul_np(
        obstruction.left_second().generators(),
        second,
        obstruction.right_second().generators(),
    );
    Ok(add_np(
        &first_multiple,
        &second_multiple,
        &first.field().one(),
        &first.field().neg(&first.field().one()),
    )?)
}

fn empty_monomial() -> NpMonomial {
    NpMonomial::from(&[][..])
}

fn reversed_left_occurrence(pattern: &NpMonomial, word: &NpMonomial) -> usize {
    let reversed_pattern = NpMonomial::from(
        pattern
            .generators()
            .iter()
            .rev()
            .copied()
            .collect::<Vec<_>>(),
    );
    let reversed_word =
        NpMonomial::from(word.generators().iter().rev().copied().collect::<Vec<_>>());
    left_occurrence(&reversed_pattern, &reversed_word)
}

fn sort_as_gap_set(obstructions: &mut Vec<Obstruction>) {
    obstructions.sort_by(compare_gap_obstructions);
    obstructions.dedup();
}

fn compare_gap_obstructions(left: &Obstruction, right: &Obstruction) -> Ordering {
    left.left_first
        .generators()
        .cmp(right.left_first.generators())
        .then_with(|| left.first_basis_index.cmp(&right.first_basis_index))
        .then_with(|| {
            left.right_first
                .generators()
                .cmp(right.right_first.generators())
        })
        .then_with(|| {
            left.left_second
                .generators()
                .cmp(right.left_second.generators())
        })
        .then_with(|| left.second_basis_index.cmp(&right.second_basis_index))
        .then_with(|| {
            left.right_second
                .generators()
                .cmp(right.right_second.generators())
        })
}

/// Failure in GBNP's source-translated obstruction routines.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ObstructionError {
    /// An obstruction referred to a polynomial outside the supplied basis.
    #[error("basis index {index} is outside a basis of length {basis_length}")]
    BasisIndexOutOfRange {
        /// Invalid zero-based index.
        index: usize,
        /// Number of basis polynomials or leading monomials.
        basis_length: usize,
    },
    /// NP scaled addition failed, normally because fields differ.
    #[error(transparent)]
    Np(#[from] NpError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring, RingOps,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, RationalField},
    };
    use serde_json::{Value, json};

    use crate::{NpMonomial, NpPolynomial, clean_np};

    use super::{
        Obstruction, left_obstructions, right_obstructions, self_obstruction, spolynomial,
    };

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

    fn obstruction_json(obstruction: &Obstruction) -> Value {
        json!([
            obstruction.left_first().generators(),
            obstruction.first_basis_index() + 1,
            obstruction.right_first().generators(),
            obstruction.left_second().generators(),
            obstruction.second_basis_index() + 1,
            obstruction.right_second().generators()
        ])
    }

    fn obstructions_json(obstructions: &[Obstruction]) -> Value {
        Value::Array(obstructions.iter().map(obstruction_json).collect())
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the obstruction fixture uses only small integers")
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
    fn aba_bab_obstructions_follow_pinned_source_pruning() {
        let leading = vec![
            NpMonomial::from(&[1, 2, 1][..]),
            NpMonomial::from(&[2, 1, 2][..]),
        ];
        let self_second = self_obstruction(1, &leading).unwrap().unwrap();
        let left = left_obstructions(1, &leading, Some(&self_second)).unwrap();
        let right = right_obstructions(1, &leading, Some(&self_second)).unwrap();
        assert_eq!(left.len(), 1);
        assert_eq!(left[0].right_first().generators(), [1]);
        assert_eq!(left[0].left_second().generators(), [2]);
        assert_eq!(right.len(), 2);
        assert_eq!(right[0].left_first().generators(), [1]);
        assert_eq!(right[1], self_second.reversed_sides());
    }

    #[test]
    fn spoly_cancels_the_shared_overlap_word() {
        let leading = vec![
            NpMonomial::from(&[1, 2, 1][..]),
            NpMonomial::from(&[2, 1, 2][..]),
        ];
        let self_second = self_obstruction(1, &leading).unwrap().unwrap();
        let obstruction = left_obstructions(1, &leading, Some(&self_second)).unwrap();
        let basis = [
            polynomial(&[&[1, 2, 1], &[1]], &[1, -1]),
            polynomial(&[&[2, 1, 2], &[2]], &[1, 2]),
        ];
        let result = spolynomial(&obstruction[0], &basis).unwrap();
        assert_eq!(result.terms().len(), 1);
        assert_eq!(result.terms()[0].monomial().generators(), [2, 1]);
        assert_eq!(result.terms()[0].coefficient(), &Q.nth(3.into()));
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_gbnp_obstruction_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/gbnp_obstructions.schema-v1.json"
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

        let leading = vec![
            NpMonomial::from(&[1, 2, 1][..]),
            NpMonomial::from(&[2, 1, 2][..]),
        ];
        let self_first = self_obstruction(0, &leading).unwrap().unwrap();
        let self_second = self_obstruction(1, &leading).unwrap().unwrap();
        let left_second = left_obstructions(1, &leading, Some(&self_second)).unwrap();
        let right_second = right_obstructions(1, &leading, Some(&self_second)).unwrap();

        let rational_basis = [
            polynomial(&[&[1, 2, 1], &[1]], &[1, -1]),
            polynomial(&[&[2, 1, 2], &[2]], &[1, 2]),
        ];
        let rational_spolynomial = spolynomial(&left_second[0], &rational_basis).unwrap();

        let field5 = Zp::new(5);
        let f5_basis = [
            clean_np(
                field5.clone(),
                vec![NpMonomial::from(&[1, 2, 1][..]), NpMonomial::from(&[1][..])],
                vec![field5.one(), field5.neg(&field5.one())],
            )
            .unwrap(),
            clean_np(
                field5.clone(),
                vec![NpMonomial::from(&[2, 1, 2][..]), NpMonomial::from(&[2][..])],
                vec![field5.one(), field5.to_element(2)],
            )
            .unwrap(),
        ];
        let f5_spolynomial = spolynomial(&left_second[0], &f5_basis).unwrap();

        let one2 = Z2.one();
        let f2_basis = [
            clean_np(
                Z2,
                vec![NpMonomial::from(&[1, 2, 1][..]), NpMonomial::from(&[1][..])],
                vec![one2, one2],
            )
            .unwrap(),
            clean_np(
                Z2,
                vec![NpMonomial::from(&[2, 1, 2][..]), NpMonomial::from(&[2][..])],
                vec![one2, Z2.zero()],
            )
            .unwrap(),
        ];
        let f2_spolynomial = spolynomial(&left_second[0], &f2_basis).unwrap();

        assert_eq!(
            json!({
                "schema_version": 1,
                "kind": "gbnp_obstructions",
                "oracle": expected["oracle"].clone(),
                "leading_monomials": leading
                    .iter()
                    .map(NpMonomial::generators)
                    .collect::<Vec<_>>(),
                "self_first": obstructions_json(std::slice::from_ref(&self_first)),
                "self_second": obstructions_json(std::slice::from_ref(&self_second)),
                "left_second": obstructions_json(&left_second),
                "right_second": obstructions_json(&right_second),
                "rational_spolynomial": rational_json(&rational_spolynomial),
                "f5_spolynomial": zp_json(&f5_spolynomial, &field5),
                "f2_spolynomial": z2_json(&f2_spolynomial)
            }),
            expected
        );
    }
}
