//! QPA's path-vector tip reduction and coordinate extraction helpers.
//!
//! These are direct translations of the three operations in `pathalg.gi`.
//! In particular, the source's documented non-span bug is represented as an
//! explicit failure witness containing the original basis and vector.

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    Path, PathAlgebra, PathAlgebraElement, PathAlgebraError, PathError, Quiver, QuotientError,
    QuotientPathAlgebra, QuotientPathAlgebraElement,
};

/// QPA's one-based leading component, path, and coefficient triple.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VectorTip<C> {
    position: usize,
    monomial: Path,
    coefficient: C,
}

impl<C> VectorTip<C> {
    /// Return the one-based first position carrying the greatest tip path.
    #[must_use]
    pub const fn position(&self) -> usize {
        self.position
    }

    /// Return the greatest tip monomial among the vector components.
    #[must_use]
    pub const fn monomial(&self) -> &Path {
        &self.monomial
    }

    /// Return the tip coefficient in the selected component.
    #[must_use]
    pub const fn coefficient(&self) -> &C {
        &self.coefficient
    }
}

/// Result of QPA's bug-compatible `CoefficientsOfVectors` operation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum VectorCoefficients<E, C> {
    /// The vector reduced to zero and these are its coefficients.
    Coefficients(Vec<C>),
    /// QPA returns `[basis, vector]` when no matching tip exists.
    NonSpanWitness {
        /// The original candidate basis.
        basis: Vec<Vec<E>>,
        /// The original vector that failed to reduce.
        vector: Vec<E>,
    },
}

trait VectorEntry<F: Field>: Clone + Eq {
    fn same_family(&self, other: &Self) -> bool;
    fn leading_monomial_entry(&self) -> Path;
    fn leading_coefficient_entry(&self) -> F::Element;
    fn scale_entry(&self, scalar: &F::Element) -> Result<Self, VectorReductionError>;
    fn subtract_entry(&self, right: &Self) -> Result<Self, VectorReductionError>;
}

impl<F: Field> VectorEntry<F> for PathAlgebraElement<F> {
    fn same_family(&self, other: &Self) -> bool {
        self.algebra().is_identical(other.algebra())
    }

    fn leading_monomial_entry(&self) -> Path {
        self.leading_monomial()
    }

    fn leading_coefficient_entry(&self) -> F::Element {
        self.leading_coefficient()
    }

    fn scale_entry(&self, scalar: &F::Element) -> Result<Self, VectorReductionError> {
        Ok(self.scale_left(scalar)?)
    }

    fn subtract_entry(&self, right: &Self) -> Result<Self, VectorReductionError> {
        Ok(self.sub(right)?)
    }
}

impl<F: Field> VectorEntry<F> for QuotientPathAlgebraElement<F> {
    fn same_family(&self, other: &Self) -> bool {
        self.algebra().is_identical(other.algebra())
    }

    fn leading_monomial_entry(&self) -> Path {
        self.leading_monomial()
    }

    fn leading_coefficient_entry(&self) -> F::Element {
        self.leading_coefficient()
    }

    fn scale_entry(&self, scalar: &F::Element) -> Result<Self, VectorReductionError> {
        Ok(self.scale_left(scalar)?)
    }

    fn subtract_entry(&self, right: &Self) -> Result<Self, VectorReductionError> {
        Ok(self.sub(right)?)
    }
}

impl<F: Field> PathAlgebra<F> {
    /// Translate `TipMonomialandCoefficientOfVector` for parent-algebra
    /// elements.
    ///
    /// # Errors
    ///
    /// QPA indexes an empty vector out of range; Rust returns a typed error.
    /// All components must belong to this identical algebra.
    pub fn tip_monomial_and_coefficient_of_vector(
        &self,
        vector: &[PathAlgebraElement<F>],
    ) -> Result<VectorTip<F::Element>, VectorReductionError> {
        vector_tip(self.field(), self.quiver(), &self.zero(), vector)
    }

    /// Translate QPA's in-place `TipReduceVectors` loop.
    ///
    /// The Boolean is false exactly for QPA's empty-input fallthrough and true
    /// when QPA returns the mutated list.
    ///
    /// # Errors
    ///
    /// Rows must have one common length and all entries must belong to this
    /// identical path algebra.
    pub fn tip_reduce_vectors(
        &self,
        vectors: &mut [Vec<PathAlgebraElement<F>>],
    ) -> Result<bool, VectorReductionError> {
        tip_reduce_vectors(self.field(), self.quiver(), &self.zero(), vectors)
    }

    /// Translate QPA's bug-compatible `CoefficientsOfVectors` loop.
    ///
    /// # Errors
    ///
    /// Basis rows must match the vector length and all entries must belong to
    /// this identical path algebra.
    pub fn coefficients_of_vectors(
        &self,
        basis: &[Vec<PathAlgebraElement<F>>],
        vector: &[PathAlgebraElement<F>],
    ) -> Result<VectorCoefficients<PathAlgebraElement<F>, F::Element>, VectorReductionError> {
        coefficients_of_vectors(self.field(), self.quiver(), &self.zero(), basis, vector)
    }
}

impl<F: Field> QuotientPathAlgebra<F> {
    /// Translate `TipMonomialandCoefficientOfVector` for quotient elements.
    ///
    /// # Errors
    ///
    /// QPA indexes an empty vector out of range; Rust returns a typed error.
    /// All components must belong to this identical quotient.
    pub fn tip_monomial_and_coefficient_of_vector(
        &self,
        vector: &[QuotientPathAlgebraElement<F>],
    ) -> Result<VectorTip<F::Element>, VectorReductionError> {
        vector_tip(
            self.original_path_algebra().field(),
            self.original_path_algebra().quiver(),
            &self.zero(),
            vector,
        )
    }

    /// Translate QPA's in-place `TipReduceVectors` loop for quotient vectors.
    ///
    /// # Errors
    ///
    /// Rows must have one common length and all entries must belong to this
    /// identical quotient.
    pub fn tip_reduce_vectors(
        &self,
        vectors: &mut [Vec<QuotientPathAlgebraElement<F>>],
    ) -> Result<bool, VectorReductionError> {
        tip_reduce_vectors(
            self.original_path_algebra().field(),
            self.original_path_algebra().quiver(),
            &self.zero(),
            vectors,
        )
    }

    /// Translate QPA's bug-compatible `CoefficientsOfVectors` loop for
    /// quotient vectors.
    ///
    /// # Errors
    ///
    /// Basis rows must match the vector length and all entries must belong to
    /// this identical quotient.
    pub fn coefficients_of_vectors(
        &self,
        basis: &[Vec<QuotientPathAlgebraElement<F>>],
        vector: &[QuotientPathAlgebraElement<F>],
    ) -> Result<VectorCoefficients<QuotientPathAlgebraElement<F>, F::Element>, VectorReductionError>
    {
        coefficients_of_vectors(
            self.original_path_algebra().field(),
            self.original_path_algebra().quiver(),
            &self.zero(),
            basis,
            vector,
        )
    }
}

fn vector_tip<F: Field, E: VectorEntry<F>>(
    _field: &F,
    quiver: &Quiver,
    zero: &E,
    vector: &[E],
) -> Result<VectorTip<F::Element>, VectorReductionError> {
    let Some(first) = vector.first() else {
        return Err(VectorReductionError::EmptyVector);
    };
    ensure_family(zero, vector)?;
    let mut position = 0;
    let mut monomial = first.leading_monomial_entry();
    for (candidate_position, candidate) in vector.iter().enumerate().skip(1) {
        let candidate_monomial = candidate.leading_monomial_entry();
        if quiver
            .compare_paths(&monomial, &candidate_monomial)?
            .is_lt()
        {
            position = candidate_position;
            monomial = candidate_monomial;
        }
    }
    Ok(VectorTip {
        position: position + 1,
        monomial,
        coefficient: vector[position].leading_coefficient_entry(),
    })
}

fn tip_reduce_vectors<F: Field, E: VectorEntry<F>>(
    field: &F,
    quiver: &Quiver,
    zero: &E,
    vectors: &mut [Vec<E>],
) -> Result<bool, VectorReductionError> {
    let Some(width) = vectors.first().map(Vec::len) else {
        return Ok(false);
    };
    ensure_rows(zero, vectors, width)?;
    for reducer in 0..vectors.len() {
        for target in 0..vectors.len() {
            if reducer == target {
                continue;
            }
            let source_tip = vector_tip(field, quiver, zero, &vectors[reducer])?;
            let target_tip = vector_tip(field, quiver, zero, &vectors[target])?;
            if source_tip.position == target_tip.position
                && source_tip.monomial == target_tip.monomial
                && !field.is_zero(&source_tip.coefficient)
            {
                let scale = field.div(&target_tip.coefficient, &source_tip.coefficient);
                let reducer_row = vectors[reducer].clone();
                for (entry, source) in vectors[target].iter_mut().zip(reducer_row) {
                    *entry = entry.subtract_entry(&source.scale_entry(&scale)?)?;
                }
            }
        }
    }
    Ok(true)
}

fn coefficients_of_vectors<F: Field, E: VectorEntry<F>>(
    field: &F,
    quiver: &Quiver,
    zero: &E,
    basis: &[Vec<E>],
    vector: &[E],
) -> Result<VectorCoefficients<E, F::Element>, VectorReductionError> {
    ensure_family(zero, vector)?;
    ensure_rows(zero, basis, vector.len())?;
    let tip_list = basis
        .iter()
        .map(|row| vector_tip(field, quiver, zero, row))
        .collect::<Result<Vec<_>, _>>()?;
    let mut coefficients = vec![field.zero(); basis.len()];
    let mut remainder = vector.to_vec();
    if remainder.iter().all(|entry| entry == zero) {
        return Ok(VectorCoefficients::Coefficients(coefficients));
    }

    loop {
        let remainder_tip = vector_tip(field, quiver, zero, &remainder)?;
        let Some(position) = tip_list.iter().position(|tip| {
            tip.position == remainder_tip.position && tip.monomial == remainder_tip.monomial
        }) else {
            return Ok(VectorCoefficients::NonSpanWitness {
                basis: basis.to_vec(),
                vector: vector.to_vec(),
            });
        };
        let basis_tip = vector_tip(field, quiver, zero, &basis[position])?;
        if field.is_zero(&basis_tip.coefficient) {
            return Err(VectorReductionError::ZeroTipCoefficient);
        }
        let scale = field.div(&remainder_tip.coefficient, &basis_tip.coefficient);
        for (entry, source) in remainder.iter_mut().zip(&basis[position]) {
            *entry = entry.subtract_entry(&source.scale_entry(&scale)?)?;
        }
        field.add_assign(&mut coefficients[position], &scale);
        if remainder.iter().all(|entry| entry == zero) {
            return Ok(VectorCoefficients::Coefficients(coefficients));
        }
    }
}

fn ensure_family<F: Field, E: VectorEntry<F>>(
    zero: &E,
    entries: &[E],
) -> Result<(), VectorReductionError> {
    if entries.iter().all(|entry| entry.same_family(zero)) {
        Ok(())
    } else {
        Err(VectorReductionError::DifferentAlgebra)
    }
}

fn ensure_rows<F: Field, E: VectorEntry<F>>(
    zero: &E,
    rows: &[Vec<E>],
    width: usize,
) -> Result<(), VectorReductionError> {
    for (row, entries) in rows.iter().enumerate() {
        if entries.len() != width {
            return Err(VectorReductionError::RaggedVectors {
                row,
                actual: entries.len(),
                expected: width,
            });
        }
        ensure_family(zero, entries)?;
    }
    Ok(())
}

/// Failure in QPA's path-vector reduction helpers.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum VectorReductionError {
    /// QPA's tip helper indexes the first entry of an empty vector.
    #[error("QPA's vector-tip operation requires a nonempty vector")]
    EmptyVector,
    /// A component belongs to another path or quotient algebra family.
    #[error("a vector component belongs to a different algebra")]
    DifferentAlgebra,
    /// GAP list subtraction requires all vectors to have one length.
    #[error("vector row {row} has length {actual}, expected {expected}")]
    RaggedVectors {
        /// Zero-based row position.
        row: usize,
        /// Received component count.
        actual: usize,
        /// Required component count.
        expected: usize,
    },
    /// The selected basis row has zero tip coefficient.
    #[error("cannot divide by a zero vector-tip coefficient")]
    ZeroTipCoefficient,
    /// Default QPA path comparison failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// Parent path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Quotient arithmetic failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Ring,
        finite_field::Z2,
        integer::Integer,
        rational::{Q, RationalField},
    };
    use serde_json::{Value, json};

    use crate::{
        ArrowSpec, PathAlgebra, PathAlgebraElement, Quiver, QuotientPathAlgebra, VertexIndex,
    };

    use super::{VectorCoefficients, VectorReductionError};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the vector-reduction fixture uses only small integers")
    }

    fn rational_element_json(element: &PathAlgebraElement<RationalField>) -> Value {
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

    fn rational_rows_json(rows: &[Vec<PathAlgebraElement<RationalField>>]) -> Value {
        Value::Array(
            rows.iter()
                .map(|row| Value::Array(row.iter().map(rational_element_json).collect()))
                .collect(),
        )
    }

    fn z2_element_json(element: &PathAlgebraElement<Z2>) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().external_representation(),
                        "coefficient": term.coefficient()
                    })
                })
                .collect(),
        )
    }

    #[test]
    #[allow(clippy::many_single_char_names, clippy::too_many_lines)]
    fn pinned_qpa_vector_reduction_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_vector_reduction.schema-v1.json"
        ))
        .unwrap();
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let generators = algebra.generators();
        let v1 = &generators[0];
        let v2 = &generators[1];
        let a = &generators[2];
        let q2 = Q.nth(Integer::from(2));
        let q3 = Q.nth(Integer::from(3));
        let q5 = Q.nth(Integer::from(5));

        let tip = algebra
            .tip_monomial_and_coefficient_of_vector(&[
                v1.add(a).unwrap(),
                a.scale_left(&q2).unwrap(),
            ])
            .unwrap();
        let zero_tip = algebra
            .tip_monomial_and_coefficient_of_vector(&[algebra.zero(), algebra.zero()])
            .unwrap();
        let mut rows = vec![
            vec![a.clone(), algebra.zero()],
            vec![a.scale_left(&q2).unwrap(), v2.clone()],
        ];
        assert!(algebra.tip_reduce_vectors(&mut rows).unwrap());
        let basis = vec![
            vec![a.clone(), algebra.zero()],
            vec![algebra.zero(), v2.clone()],
        ];
        let target = vec![a.scale_left(&q3).unwrap(), v2.scale_left(&q5).unwrap()];
        let VectorCoefficients::Coefficients(coefficients) =
            algebra.coefficients_of_vectors(&basis, &target).unwrap()
        else {
            panic!("oracle target must lie in the vector span");
        };
        let not_span = vec![v1.clone(), algebra.zero()];
        let failure_witness = matches!(
            algebra.coefficients_of_vectors(&basis, &not_span).unwrap(),
            VectorCoefficients::NonSpanWitness {
                basis: returned_basis,
                vector: returned_vector
            } if returned_basis == basis && returned_vector == not_span
        );

        let quotient = QuotientPathAlgebra::by_relations_high_level(algebra.clone(), &[]).unwrap();
        let quotient_generators = quotient.generators().unwrap();
        let quotient_basis = vec![
            vec![quotient_generators[2].clone(), quotient.zero()],
            vec![quotient.zero(), quotient_generators[1].clone()],
        ];
        let quotient_target = vec![
            quotient_generators[2].scale_left(&q3).unwrap(),
            quotient_generators[1].scale_left(&q5).unwrap(),
        ];
        let VectorCoefficients::Coefficients(quotient_coefficients) = quotient
            .coefficients_of_vectors(&quotient_basis, &quotient_target)
            .unwrap()
        else {
            panic!("oracle quotient target must lie in the vector span");
        };

        assert_eq!(
            json!({
                "tip": {
                    "position": tip.position(),
                    "path": tip.monomial().external_representation(),
                    "coefficient": [
                        small_integer(tip.coefficient().numerator_ref()),
                        small_integer(tip.coefficient().denominator_ref())
                    ]
                },
                "zero_tip": {
                    "position": zero_tip.position(),
                    "path": zero_tip.monomial().external_representation(),
                    "coefficient": [
                        small_integer(zero_tip.coefficient().numerator_ref()),
                        small_integer(zero_tip.coefficient().denominator_ref())
                    ]
                },
                "same_rows": true,
                "reduced": rational_rows_json(&rows),
                "coefficients": coefficients
                    .iter()
                    .map(|coefficient| json!([
                        small_integer(coefficient.numerator_ref()),
                        small_integer(coefficient.denominator_ref())
                    ]))
                    .collect::<Vec<_>>(),
                "failure_witness": failure_witness,
                "quotient_coefficients": quotient_coefficients
                    .iter()
                    .map(|coefficient| json!([
                        small_integer(coefficient.numerator_ref()),
                        small_integer(coefficient.denominator_ref())
                    ]))
                    .collect::<Vec<_>>()
            }),
            expected["rational"]
        );

        let algebra2 = PathAlgebra::new(Z2, quiver).unwrap();
        let generators2 = algebra2.generators();
        let mut rows2 = vec![
            vec![generators2[2].clone(), algebra2.zero()],
            vec![generators2[2].clone(), generators2[1].clone()],
        ];
        algebra2.tip_reduce_vectors(&mut rows2).unwrap();
        let basis2 = vec![
            vec![generators2[2].clone(), algebra2.zero()],
            vec![algebra2.zero(), generators2[1].clone()],
        ];
        let VectorCoefficients::Coefficients(coefficients2) = algebra2
            .coefficients_of_vectors(&basis2, &[generators2[2].clone(), generators2[1].clone()])
            .unwrap()
        else {
            panic!("binary oracle target must lie in the vector span");
        };
        assert_eq!(
            json!({
                "reduced": rows2
                    .iter()
                    .map(|row| Value::Array(row.iter().map(z2_element_json).collect()))
                    .collect::<Vec<_>>(),
                "coefficients": coefficients2
            }),
            expected["f2"]
        );

        assert_eq!(
            algebra
                .tip_monomial_and_coefficient_of_vector(&[])
                .unwrap_err(),
            VectorReductionError::EmptyVector
        );
        assert!(!algebra.tip_reduce_vectors(&mut []).unwrap());
    }
}
