//! QPA's algebra-derived integral forms from `lib/combinatorialrep.gi`.

use numerica::{
    domains::{
        Field, Ring,
        integer::{Integer, IntegerRing, Z},
        rational::Q,
    },
    tensors::matrix::Matrix,
};
use qpa_core::{
    DenseMatrix, MatrixError, PathEndpoint, QuiverAlgebra, QuotientPathAlgebra, UnitForm,
    UnitFormError, tits_unit_form_of_path_algebra,
};
use thiserror::Error;

use crate::module::{EnvelopingModuleError, cartan_matrix};
use crate::pathalg::{
    PathalgGbnpError, ideal_minimal_generating_set, quotient_is_finite_dimensional,
};

/// Translate QPA's complete `TitsUnitFormOfAlgebra` dispatch.
///
/// # Errors
///
/// Propagates the source path-algebra or admissible-quotient method selected
/// by the entered algebra.
pub fn tits_unit_form_of_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<UnitForm, CombinatorialRepresentationError> {
    match algebra {
        QuiverAlgebra::Path(algebra) => Ok(tits_unit_form_of_path_algebra(algebra)?),
        QuiverAlgebra::Quotient(algebra) => tits_unit_form_of_quotient(algebra),
    }
}

/// Translate QPA's admissible-quotient method for `TitsUnitFormOfAlgebra`.
///
/// The minimal-relation space is obtained from QPA's translated
/// `MinimalGeneratingSetOfIdeal`, which computes the same `I/(JI+IJ)` basis
/// constructed inline by the installed Tits method.  Its source/target block
/// dimensions form the `Ext^2` matrix used by QPA.
///
/// # Errors
///
/// Rejects an infinite quotient, a loop, or a minimal relation with equal
/// source and target, and propagates QPA's exact ideal and unit-form backend.
pub fn tits_unit_form_of_quotient<F: Field>(
    algebra: &QuotientPathAlgebra<F>,
) -> Result<UnitForm, CombinatorialRepresentationError> {
    if !quotient_is_finite_dimensional(algebra)? {
        return Err(CombinatorialRepresentationError::InfiniteTitsAlgebra);
    }
    let quiver = algebra.original_path_algebra().quiver();
    let adjacency = quiver.adjacency_matrix();
    let dimension = quiver.number_of_vertices();
    if (0..dimension).any(|vertex| adjacency[vertex][vertex] != 0) {
        return Err(CombinatorialRepresentationError::TitsQuiverHasLoop);
    }

    let minimal_relations = ideal_minimal_generating_set(algebra.ideal())?;
    let mut ext_two = vec![vec![0_usize; dimension]; dimension];
    for relation in minimal_relations {
        let Some(first) = relation.terms().first() else {
            continue;
        };
        let source = first.path().source();
        let target = first.path().target();
        if relation
            .terms()
            .iter()
            .any(|term| term.path().source() != source || term.path().target() != target)
        {
            return Err(CombinatorialRepresentationError::NonUniformMinimalRelation);
        }
        let (PathEndpoint::Vertex(source), PathEndpoint::Vertex(target)) = (source, target) else {
            return Err(CombinatorialRepresentationError::InvalidMinimalRelationEndpoint);
        };
        let source = source.index().get() - 1;
        let target = target.index().get() - 1;
        ext_two[source][target] = ext_two[source][target]
            .checked_add(1)
            .ok_or(CombinatorialRepresentationError::ExtTwoDimensionOverflow)?;
    }
    if (0..dimension).any(|vertex| ext_two[vertex][vertex] != 0) {
        return Err(CombinatorialRepresentationError::CyclicMinimalRelation);
    }

    let temporary = (0..dimension)
        .map(|row| {
            (0..dimension)
                .map(|column| {
                    Integer::from(usize::from(row == column))
                        - Integer::from(adjacency[row][column])
                        + Integer::from(ext_two[row][column])
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    Ok(UnitForm::new(
        (0..dimension)
            .map(|row| {
                (0..dimension)
                    .map(|column| &temporary[row][column] + &temporary[column][row])
                    .collect()
            })
            .collect(),
    )?)
}

/// QPA's captured Euler bilinear form with its integral inverse Cartan matrix.
///
/// QPA only constructs this closure when the Cartan determinant is `1` or
/// `-1`, so the inverse matrix has integral entries.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct EulerBilinearForm {
    inverse_cartan: Matrix<IntegerRing>,
}

impl EulerBilinearForm {
    /// Return the integral inverse Cartan matrix captured by QPA's closure.
    #[must_use]
    pub const fn inverse_cartan_matrix(&self) -> &Matrix<IntegerRing> {
        &self.inverse_cartan
    }

    /// Return the dimension of the form.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.inverse_cartan.nrows()
    }

    /// Evaluate QPA's Euler closure `x * CartanMatrix(A)^-1 * y`.
    ///
    /// # Errors
    ///
    /// Rejects vectors of unequal length or a length different from the
    /// Cartan dimension.
    pub fn evaluate<R: Ring>(
        &self,
        ring: &R,
        left: &[R::Element],
        right: &[R::Element],
    ) -> Result<R::Element, CombinatorialRepresentationError> {
        if left.len() != right.len() {
            return Err(CombinatorialRepresentationError::VectorLengthsDiffer {
                left: left.len(),
                right: right.len(),
            });
        }
        if left.len() != self.dimension() {
            return Err(CombinatorialRepresentationError::WrongVectorSize {
                actual: left.len(),
                expected: self.dimension(),
            });
        }
        let mut value = ring.zero();
        for (row, left_entry) in left.iter().enumerate() {
            for (column, right_entry) in right.iter().enumerate() {
                let row = u32::try_from(row).map_err(|_| {
                    CombinatorialRepresentationError::MatrixDimensionTooLarge {
                        dimension: self.dimension(),
                    }
                })?;
                let column = u32::try_from(column).map_err(|_| {
                    CombinatorialRepresentationError::MatrixDimensionTooLarge {
                        dimension: self.dimension(),
                    }
                })?;
                let matrix_entry = ring.nth(self.inverse_cartan[(row, column)].clone());
                let left_product = ring.mul(left_entry, &matrix_entry);
                ring.add_mul_assign(&mut value, &left_product, right_entry);
            }
        }
        Ok(value)
    }
}

/// Translate QPA's `EulerBilinearFormOfAlgebra(A)`.
///
/// `Ok(None)` is QPA's `fail`: either `CartanMatrix(A)` was unavailable or
/// its determinant was not `1` or `-1`.  The successful branch uses the same
/// source-translated exact Cartan determinant and inverse routines as the
/// rest of this port.
///
/// # Errors
///
/// Propagates Cartan construction and exact matrix failures, and reports a
/// violated integrality invariant after a unimodular determinant.
pub fn euler_bilinear_form_of_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<EulerBilinearForm>, CombinatorialRepresentationError> {
    let Some(cartan_rows) = cartan_matrix(algebra)? else {
        return Ok(None);
    };
    let cartan = DenseMatrix::from_rows(
        cartan_rows
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|entry| Q.nth(Integer::from(entry)))
                    .collect()
            })
            .collect(),
        Q,
    )?;
    let determinant = cartan.determinant()?;
    if determinant != 1 && determinant != -1 {
        return Ok(None);
    }
    let inverse = cartan
        .inverse()?
        .ok_or(CombinatorialRepresentationError::UnimodularCartanNotInvertible)?;
    let entries = inverse
        .as_flat()
        .iter()
        .map(|entry| {
            if entry.denominator_ref() == &Integer::one() {
                Ok(entry.numerator_ref().clone())
            } else {
                Err(CombinatorialRepresentationError::NonIntegralCartanInverse)
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let nrows = u32::try_from(inverse.nrows()).map_err(|_| {
        CombinatorialRepresentationError::MatrixDimensionTooLarge {
            dimension: inverse.nrows(),
        }
    })?;
    let ncols = u32::try_from(inverse.ncols()).map_err(|_| {
        CombinatorialRepresentationError::MatrixDimensionTooLarge {
            dimension: inverse.ncols(),
        }
    })?;
    let inverse_cartan = Matrix::from_linear(entries, nrows, ncols, Z)
        .map_err(|_| CombinatorialRepresentationError::InvalidNumericaMatrix)?;
    Ok(Some(EulerBilinearForm { inverse_cartan }))
}

/// Failure in QPA's algebra-derived combinatorial representation forms.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum CombinatorialRepresentationError {
    /// QPA's Cartan and projective-module computation failed.
    #[error(transparent)]
    Enveloping(#[from] EnvelopingModuleError),
    /// An exact matrix operation failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// QPA's configured ideal and quotient backend failed.
    #[error(transparent)]
    Pathalg(#[from] PathalgGbnpError),
    /// QPA's unit-form validation failed.
    #[error(transparent)]
    UnitForm(#[from] UnitFormError),
    /// The quotient Tits method requires a finite-dimensional algebra.
    #[error("the entered algebra is not finite dimensional")]
    InfiniteTitsAlgebra,
    /// QPA rejects loops before constructing the quotient Tits form.
    #[error("the quiver of the algebra has at least one loop")]
    TitsQuiverHasLoop,
    /// QPA's minimal-relation basis is split by source and target vertices.
    #[error("a minimal relation is not source-target uniform")]
    NonUniformMinimalRelation,
    /// A nonzero path relation must have vertex endpoints.
    #[error("a minimal relation does not have vertex endpoints")]
    InvalidMinimalRelationEndpoint,
    /// A minimal-relation block dimension exceeded `usize`.
    #[error("an Ext^2 block dimension overflows usize")]
    ExtTwoDimensionOverflow,
    /// QPA rejects a minimal relation beginning and ending at one vertex.
    #[error("the algebra has a minimal relation starting and ending at one vertex")]
    CyclicMinimalRelation,
    /// A nonzero determinant unexpectedly failed to yield an inverse.
    #[error("a unimodular Cartan matrix was not invertible")]
    UnimodularCartanNotInvertible,
    /// A unimodular integer matrix unexpectedly had a nonintegral inverse.
    #[error("a unimodular Cartan matrix had a nonintegral inverse")]
    NonIntegralCartanInverse,
    /// Numerica matrix dimensions use `u32`.
    #[error("Euler-form dimension {dimension} exceeds Numerica's matrix limit")]
    MatrixDimensionTooLarge {
        /// Rejected dimension.
        dimension: usize,
    },
    /// A validated matrix could not be materialized by Numerica.
    #[error("Numerica rejected a validated Euler-form matrix")]
    InvalidNumericaMatrix,
    /// QPA requires equal vector lengths.
    #[error("Euler-form vector lengths {left} and {right} differ")]
    VectorLengthsDiffer {
        /// Left vector length.
        left: usize,
        /// Right vector length.
        right: usize,
    },
    /// QPA requires vectors of the Cartan dimension.
    #[error("Euler-form vector has length {actual}, expected {expected}")]
    WrongVectorSize {
        /// Entered vector length.
        actual: usize,
        /// Cartan dimension.
        expected: usize,
    },
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        finite_field::Z2,
        integer::{Integer, Z},
        rational::Q,
    };
    use qpa_core::{ArrowSpec, PathAlgebra, Quiver, QuiverAlgebra, UnitForm, VertexIndex};
    use serde_json::{Value, json};

    use super::{
        CombinatorialRepresentationError, EulerBilinearForm, euler_bilinear_form_of_algebra,
        tits_unit_form_of_algebra,
    };
    use crate::pathalg::quotient_by_relations;
    use crate::predefined::nakayama_algebra;

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn a3<F: numerica::domains::Field>(field: F) -> QuiverAlgebra<F> {
        QuiverAlgebra::Path(
            PathAlgebra::new(
                field,
                Quiver::from_indexed(
                    3,
                    vec![
                        ArrowSpec::named(vertex(1), vertex(2), "a"),
                        ArrowSpec::named(vertex(2), vertex(3), "b"),
                    ],
                )
                .unwrap(),
            )
            .unwrap(),
        )
    }

    fn radical_square_zero_a3<F: numerica::domains::Field>(field: F) -> QuiverAlgebra<F> {
        let parent = match a3(field) {
            QuiverAlgebra::Path(parent) => parent,
            QuiverAlgebra::Quotient(_) => unreachable!(),
        };
        let generators = parent.generators();
        let relation = generators[3].multiply(&generators[4]).unwrap();
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[relation]).unwrap())
    }

    fn dual_numbers<F: numerica::domains::Field>(field: F) -> QuiverAlgebra<F> {
        let parent = PathAlgebra::new(
            field,
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = parent.generators()[1].clone();
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap())
    }

    fn mixed_endpoint_relations<F: numerica::domains::Field>(field: F) -> QuiverAlgebra<F> {
        let parent = PathAlgebra::new(
            field,
            Quiver::from_indexed(
                6,
                vec![
                    ArrowSpec::named(vertex(1), vertex(2), "a1"),
                    ArrowSpec::named(vertex(2), vertex(3), "a2"),
                    ArrowSpec::named(vertex(4), vertex(5), "b1"),
                    ArrowSpec::named(vertex(5), vertex(6), "b2"),
                ],
            )
            .unwrap(),
        )
        .unwrap();
        let generators = parent.generators();
        let left = generators[6].multiply(&generators[7]).unwrap();
        let right = generators[8].multiply(&generators[9]).unwrap();
        let relation = left.add(&right).unwrap();
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[relation]).unwrap())
    }

    fn cyclic_relations<F: numerica::domains::Field>(field: F) -> QuiverAlgebra<F> {
        let parent = PathAlgebra::new(
            field,
            Quiver::from_indexed(
                2,
                vec![
                    ArrowSpec::named(vertex(1), vertex(2), "cycle_a"),
                    ArrowSpec::named(vertex(2), vertex(1), "cycle_b"),
                ],
            )
            .unwrap(),
        )
        .unwrap();
        let generators = parent.generators();
        let ab = generators[2].multiply(&generators[3]).unwrap();
        let ba = generators[3].multiply(&generators[2]).unwrap();
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[ab, ba]).unwrap())
    }

    fn matrix_rows(form: &EulerBilinearForm) -> Vec<Vec<i64>> {
        form.inverse_cartan_matrix()
            .row_iter()
            .map(|row| row.iter().map(|entry| entry.to_i64().unwrap()).collect())
            .collect()
    }

    fn unit_form_rows(form: &UnitForm) -> Vec<Vec<i64>> {
        form.symmetric_matrix()
            .row_iter()
            .map(|row| row.iter().map(|entry| entry.to_i64().unwrap()).collect())
            .collect()
    }

    #[test]
    fn pinned_qpa_quotient_tits_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/combinatorialrep_tits_quotient.schema-v1.json"
        ))
        .unwrap();
        let source_nakayama =
            nakayama_algebra(Q, &[Integer::from(3), Integer::from(2), Integer::from(1)])
                .unwrap()
                .as_algebra()
                .unwrap()
                .algebra()
                .clone();
        let radical_q = radical_square_zero_a3(Q);
        let radical_f2 = radical_square_zero_a3(Z2);
        let mixed = mixed_endpoint_relations(Q);
        assert_eq!(
            json!({
                "source_nakayama": unit_form_rows(
                    &tits_unit_form_of_algebra(&source_nakayama).unwrap()
                ),
                "radical_square_zero_a3": unit_form_rows(
                    &tits_unit_form_of_algebra(&radical_q).unwrap()
                ),
                "radical_square_zero_a3_f2": unit_form_rows(
                    &tits_unit_form_of_algebra(&radical_f2).unwrap()
                ),
                "mixed_endpoint_relations": unit_form_rows(
                    &tits_unit_form_of_algebra(&mixed).unwrap()
                ),
            }),
            expected["matrices"]
        );

        assert!(matches!(
            tits_unit_form_of_algebra(&dual_numbers(Q)),
            Err(CombinatorialRepresentationError::TitsQuiverHasLoop)
        ));
        assert!(matches!(
            tits_unit_form_of_algebra(&cyclic_relations(Q)),
            Err(CombinatorialRepresentationError::CyclicMinimalRelation)
        ));
    }

    #[test]
    fn pinned_qpa_euler_form_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/combinatorialrep_euler.schema-v1.json"
        ))
        .unwrap();
        let path_q = a3(Q);
        let quotient_q = radical_square_zero_a3(Q);
        let path_f2 = a3(Z2);
        let left = [Integer::from(1), Integer::from(2), Integer::from(-1)];
        let right = [Integer::from(3), Integer::from(0), Integer::from(2)];

        let path_form = euler_bilinear_form_of_algebra(&path_q).unwrap().unwrap();
        let quotient_form = euler_bilinear_form_of_algebra(&quotient_q)
            .unwrap()
            .unwrap();
        let finite_form = euler_bilinear_form_of_algebra(&path_f2).unwrap().unwrap();
        assert_eq!(
            json!({
                "path_a3": {
                    "inverse_cartan": matrix_rows(&path_form),
                    "value": path_form.evaluate(&Z, &left, &right).unwrap().to_i64().unwrap(),
                },
                "radical_square_zero_a3": {
                    "inverse_cartan": matrix_rows(&quotient_form),
                    "value": quotient_form.evaluate(&Z, &left, &right).unwrap().to_i64().unwrap(),
                },
                "path_a3_f2": {
                    "inverse_cartan": matrix_rows(&finite_form),
                    "value": finite_form.evaluate(&Z, &left, &right).unwrap().to_i64().unwrap(),
                },
            }),
            expected["forms"]
        );

        assert!(
            euler_bilinear_form_of_algebra(&dual_numbers(Q))
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn euler_form_argument_failures_are_typed() {
        let form = euler_bilinear_form_of_algebra(&a3(Q)).unwrap().unwrap();
        assert!(matches!(
            form.evaluate(&Z, &[Integer::one()], &[Integer::one(), Integer::zero()]),
            Err(CombinatorialRepresentationError::VectorLengthsDiffer { .. })
        ));
        assert!(matches!(
            form.evaluate(&Z, &[Integer::one()], &[Integer::one()]),
            Err(CombinatorialRepresentationError::WrongVectorSize { .. })
        ));
    }
}
