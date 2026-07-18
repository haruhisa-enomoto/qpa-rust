//! QPA's integral unit forms and source algorithms.
//!
//! The implementations in this module translate `lib/combinatorialrep.gi`.
//! Numerica supplies arbitrary-precision integers and the stored matrix, while
//! QPA's reflection and positivity procedures remain the defining algorithms.

use std::sync::OnceLock;

use numerica::{
    domains::{
        Field, Ring,
        integer::{Integer, IntegerRing, Z},
    },
    tensors::matrix::Matrix,
};
use thiserror::Error;

use crate::path_algebra::PathAlgebra;

/// Translate QPA's path-algebra method for `TitsUnitFormOfAlgebra`.
///
/// QPA forms `I - adjacency(Q)` and adds its transpose before passing the
/// result through [`UnitForm::new`].  In particular, a loop reaches the same
/// diagonal validation failure rather than being prefiltered here.
///
/// # Errors
///
/// Propagates QPA's unit-form validation of the resulting symmetric matrix.
pub fn tits_unit_form_of_path_algebra<F: Field>(
    algebra: &PathAlgebra<F>,
) -> Result<UnitForm, UnitFormError> {
    let adjacency = algebra.quiver().adjacency_matrix();
    let dimension = algebra.quiver().number_of_vertices();
    let temporary = (0..dimension)
        .map(|row| {
            (0..dimension)
                .map(|column| {
                    Integer::from(usize::from(row == column))
                        - Integer::from(adjacency[row][column])
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    UnitForm::new(
        (0..dimension)
            .map(|row| {
                (0..dimension)
                    .map(|column| &temporary[row][column] + &temporary[column][row])
                    .collect()
            })
            .collect(),
    )
}

/// An integral symmetric matrix with every diagonal entry equal to two.
///
/// This is QPA's `IsUnitForm` representation.  The positive-root attribute is
/// populated when [`Self::is_weakly_positive`] succeeds, just as in QPA.
#[derive(Debug)]
pub struct UnitForm {
    matrix: Matrix<IntegerRing>,
    weakly_nonnegative: OnceLock<bool>,
    weakly_positive: OnceLock<bool>,
    positive_roots: OnceLock<Vec<Vec<Integer>>>,
}

impl UnitForm {
    /// Construct QPA's `UnitForm(B)` from an integer row matrix.
    ///
    /// # Errors
    ///
    /// Rejects an empty or ragged row list, a nonsquare matrix, a diagonal
    /// entry other than two, a nonsymmetric matrix, or a dimension not
    /// representable by Numerica.
    pub fn new(rows: Vec<Vec<Integer>>) -> Result<Self, UnitFormError> {
        let nrows = rows.len();
        if nrows == 0 {
            return Err(UnitFormError::EmptyMatrix);
        }
        let ncols = rows.first().map_or(0, Vec::len);
        if let Some((row, actual)) = rows
            .iter()
            .enumerate()
            .find_map(|(row, entries)| (entries.len() != ncols).then_some((row, entries.len())))
        {
            return Err(UnitFormError::RaggedRows {
                row,
                actual,
                expected: ncols,
            });
        }
        if nrows != ncols {
            return Err(UnitFormError::NonSquare { nrows, ncols });
        }
        let dimension = u32::try_from(nrows)
            .map_err(|_| UnitFormError::DimensionTooLarge { dimension: nrows })?;
        for (index, row) in rows.iter().enumerate() {
            if row[index] != 2 {
                return Err(UnitFormError::DiagonalNotTwo { index });
            }
        }
        for row in 0..nrows {
            for column in row + 1..ncols {
                if rows[row][column] != rows[column][row] {
                    return Err(UnitFormError::NotSymmetric { row, column });
                }
            }
        }

        let matrix = Matrix::from_linear(
            rows.into_iter().flatten().collect(),
            dimension,
            dimension,
            Z,
        )
        .map_err(|_| UnitFormError::InvalidNumericaMatrix)?;
        Ok(Self {
            matrix,
            weakly_nonnegative: OnceLock::new(),
            weakly_positive: OnceLock::new(),
            positive_roots: OnceLock::new(),
        })
    }

    /// Return QPA's `SymmetricMatrixOfUnitForm` attribute.
    #[must_use]
    pub const fn symmetric_matrix(&self) -> &Matrix<IntegerRing> {
        &self.matrix
    }

    /// Return the dimension of the unit form.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.matrix.nrows()
    }

    /// Evaluate QPA's `BilinearFormOfUnitForm(B)(x, y)` over a Numerica ring.
    ///
    /// The integer entries of the form are mapped through [`Ring::nth`], so
    /// the same unit form can act on integer, rational, or supported finite
    /// field vectors.
    ///
    /// # Errors
    ///
    /// Rejects unequal vector lengths or vectors of the wrong dimension.
    pub fn bilinear_form<R: Ring>(
        &self,
        ring: &R,
        left: &[R::Element],
        right: &[R::Element],
    ) -> Result<R::Element, UnitFormError> {
        if left.len() != right.len() {
            return Err(UnitFormError::VectorLengthsDiffer {
                left: left.len(),
                right: right.len(),
            });
        }
        self.ensure_vector_size(left.len())?;
        Ok(self.bilinear_unchecked(ring, left, right))
    }

    /// Evaluate QPA's `QuadraticFormOfUnitForm(B)(x)` over a Numerica ring.
    ///
    /// # Errors
    ///
    /// Rejects a vector of the wrong dimension, or a coefficient ring in
    /// which division of the computed value by two is unavailable.
    pub fn quadratic_form<R: Ring>(
        &self,
        ring: &R,
        vector: &[R::Element],
    ) -> Result<R::Element, UnitFormError> {
        self.ensure_vector_size(vector.len())?;
        let doubled = self.bilinear_unchecked(ring, vector, vector);
        ring.try_div(&doubled, &ring.nth(Integer::from(2)))
            .ok_or(UnitFormError::QuadraticDivisionUnavailable)
    }

    /// Translate `ReflectionByBilinearForm(B, i, z)`.
    ///
    /// `coordinate` follows GAP's one-based convention.
    ///
    /// # Errors
    ///
    /// Rejects a coordinate outside `1..=dimension` or a vector of the wrong
    /// size.
    pub fn reflection_by_bilinear_form<R: Ring>(
        &self,
        ring: &R,
        coordinate: usize,
        vector: &[R::Element],
    ) -> Result<Vec<R::Element>, UnitFormError> {
        self.ensure_vector_size(vector.len())?;
        if coordinate == 0 || coordinate > self.dimension() {
            return Err(UnitFormError::CoordinateOutsideRange {
                coordinate,
                dimension: self.dimension(),
            });
        }
        let mut simple = vec![ring.zero(); self.dimension()];
        simple[coordinate - 1] = ring.one();
        let coefficient = self.bilinear_unchecked(ring, vector, &simple);
        let mut reflection = vector.to_vec();
        reflection[coordinate - 1] = ring.sub(&reflection[coordinate - 1], &coefficient);
        Ok(reflection)
    }

    /// Translate QPA's Dean--de la Peña weak-nonnegativity algorithm.
    #[must_use]
    pub fn is_weakly_nonnegative(&self) -> bool {
        *self
            .weakly_nonnegative
            .get_or_init(|| self.compute_weakly_nonnegative())
    }

    /// Translate QPA's von Höhne edge-reduction weak-positivity algorithm.
    ///
    /// A successful computation also binds [`Self::positive_roots`].
    #[must_use]
    pub fn is_weakly_positive(&self) -> bool {
        *self.weakly_positive.get_or_init(|| {
            let Some(roots) = self.compute_weakly_positive() else {
                return false;
            };
            let _ = self.positive_roots.set(roots);
            true
        })
    }

    /// Return QPA's `PositiveRootsOfUnitForm` attribute when it is bound.
    ///
    /// The attribute is bound only after [`Self::is_weakly_positive`] returns
    /// true.
    #[must_use]
    pub fn positive_roots(&self) -> Option<&[Vec<Integer>]> {
        self.positive_roots.get().map(Vec::as_slice)
    }

    fn ensure_vector_size(&self, actual: usize) -> Result<(), UnitFormError> {
        let expected = self.dimension();
        if actual == expected {
            Ok(())
        } else {
            Err(UnitFormError::WrongVectorSize { actual, expected })
        }
    }

    fn bilinear_unchecked<R: Ring>(
        &self,
        ring: &R,
        left: &[R::Element],
        right: &[R::Element],
    ) -> R::Element {
        let mut value = ring.zero();
        for (row, left_entry) in left.iter().enumerate() {
            for (column, right_entry) in right.iter().enumerate() {
                let matrix_entry = ring.nth(self.entry(row, column).clone());
                let left_product = ring.mul(left_entry, &matrix_entry);
                ring.add_mul_assign(&mut value, &left_product, right_entry);
            }
        }
        value
    }

    fn entry(&self, row: usize, column: usize) -> &Integer {
        let row = u32::try_from(row).expect("unit-form rows fit Numerica's index type");
        let column = u32::try_from(column).expect("unit-form columns fit Numerica's index type");
        &self.matrix[(row, column)]
    }

    fn matrix_rows(&self) -> Vec<Vec<Integer>> {
        self.matrix.row_iter().map(<[_]>::to_vec).collect()
    }

    fn compute_weakly_nonnegative(&self) -> bool {
        let simples = identity_vectors(self.dimension());
        let mut candidates = simples.clone();
        loop {
            if candidates.iter().any(|vector| {
                simples
                    .iter()
                    .any(|simple| self.integer_bilinear(vector, simple) <= -3)
            }) {
                return false;
            }
            if candidates.iter().any(|vector| {
                simples.iter().any(|simple| {
                    self.integer_bilinear(vector, simple) == -2
                        && row_times_matrix(&add_vectors(vector, simple), &self.matrix)
                            .iter()
                            .any(|entry| entry < &Integer::zero())
                })
            }) {
                return false;
            }

            let reducible = candidates
                .into_iter()
                .filter(|vector| {
                    vector.iter().all(|entry| entry <= &Integer::from(6))
                        && simples
                            .iter()
                            .any(|simple| self.integer_bilinear(vector, simple) == -1)
                })
                .collect::<Vec<_>>();
            if reducible.is_empty() {
                return true;
            }
            candidates = reducible
                .into_iter()
                .map(|vector| {
                    let coordinate = simples
                        .iter()
                        .position(|simple| self.integer_bilinear(&vector, simple) == -1)
                        .expect("QPA selected a vector with a reflecting coordinate");
                    self.reflection_by_bilinear_form(&Z, coordinate + 1, &vector)
                        .expect("QPA reflection inputs have the unit-form dimension")
                })
                .collect();
        }
    }

    fn integer_bilinear(&self, left: &[Integer], right: &[Integer]) -> Integer {
        self.bilinear_unchecked(&Z, left, right)
    }

    fn compute_weakly_positive(&self) -> Option<Vec<Vec<Integer>>> {
        let mut matrix = self.matrix_rows();
        let mut reduction_pairs = Vec::<(usize, usize)>::new();
        let mut root_count = matrix.len();
        loop {
            let size = matrix.len();
            if upper_triangle(&matrix).all(|entry| entry >= &Integer::zero()) {
                let mut roots = identity_vectors(root_count);
                let reduction_count = reduction_pairs.len();
                for reduction in (0..reduction_count).rev() {
                    let (left, right) = reduction_pairs[reduction];
                    let prefix_length = root_count + reduction - reduction_count;
                    let combined = roots[left][..prefix_length]
                        .iter()
                        .zip(&roots[right][..prefix_length])
                        .map(|(left, right)| left + right)
                        .collect::<Vec<_>>();
                    for root in &mut roots {
                        let coefficient = root[prefix_length].clone();
                        let mut shortened = root[..prefix_length].to_vec();
                        for (entry, combined_entry) in shortened.iter_mut().zip(&combined) {
                            *entry += &coefficient * combined_entry;
                        }
                        *root = shortened;
                    }
                }
                return Some(roots);
            }
            if upper_triangle(&matrix).any(|entry| entry <= &Integer::from(-2)) {
                return None;
            }

            let pair = (0..size).find_map(|column| {
                (0..=column)
                    .find(|&row| matrix[row][column] == -1)
                    .map(|row| (row, column))
            })?;
            reduction_pairs.push(pair);
            root_count += 1;

            let (left, right) = pair;
            let last_column = (0..size)
                .map(|row| &matrix[row][left] + &matrix[row][right])
                .collect::<Vec<_>>();
            let mut last_row = (0..size)
                .map(|column| &matrix[left][column] + &matrix[right][column])
                .collect::<Vec<_>>();
            for (row, entry) in matrix.iter_mut().zip(last_column) {
                row.push(entry);
            }
            last_row.push(Integer::from(2));
            matrix.push(last_row);
            matrix[left][right] = Integer::zero();
            matrix[right][left] = Integer::zero();
        }
    }
}

/// Failure to construct or evaluate a QPA unit form.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum UnitFormError {
    /// GAP's empty list does not satisfy the `IsMatrix` method filter.
    #[error("the entered value is not a matrix")]
    EmptyMatrix,
    /// A row-list input is not a matrix.
    #[error("unit-form row {row} has length {actual}, expected {expected}")]
    RaggedRows {
        /// Zero-based row index.
        row: usize,
        /// Actual row length.
        actual: usize,
        /// Expected row length.
        expected: usize,
    },
    /// QPA requires a square matrix.
    #[error("the entered matrix is not square: shape {nrows} x {ncols}")]
    NonSquare {
        /// Number of rows.
        nrows: usize,
        /// Number of columns.
        ncols: usize,
    },
    /// QPA requires all diagonal entries to equal two.
    #[error("unit-form diagonal entry {index} is not two")]
    DiagonalNotTwo {
        /// Zero-based diagonal index.
        index: usize,
    },
    /// QPA requires a symmetric matrix.
    #[error("unit-form entries ({row}, {column}) and ({column}, {row}) differ")]
    NotSymmetric {
        /// Zero-based row index.
        row: usize,
        /// Zero-based column index.
        column: usize,
    },
    /// Numerica matrix dimensions use `u32`.
    #[error("unit-form dimension {dimension} exceeds Numerica's matrix limit")]
    DimensionTooLarge {
        /// Rejected dimension.
        dimension: usize,
    },
    /// A validated matrix could not be materialized by Numerica.
    #[error("Numerica rejected a validated unit-form matrix")]
    InvalidNumericaMatrix,
    /// QPA requires equal bilinear-form vector lengths.
    #[error("bilinear-form vector lengths {left} and {right} differ")]
    VectorLengthsDiffer {
        /// Left vector length.
        left: usize,
        /// Right vector length.
        right: usize,
    },
    /// A form argument has the wrong dimension.
    #[error("unit-form vector has length {actual}, expected {expected}")]
    WrongVectorSize {
        /// Entered vector length.
        actual: usize,
        /// Unit-form dimension.
        expected: usize,
    },
    /// QPA reflection coordinates are one-based.
    #[error("reflection coordinate {coordinate} is outside 1..={dimension}")]
    CoordinateOutsideRange {
        /// Entered one-based coordinate.
        coordinate: usize,
        /// Unit-form dimension.
        dimension: usize,
    },
    /// The target ring cannot perform QPA's division by two.
    #[error("the quadratic-form value cannot be divided by two in the coefficient ring")]
    QuadraticDivisionUnavailable,
}

fn identity_vectors(size: usize) -> Vec<Vec<Integer>> {
    (0..size)
        .map(|row| {
            (0..size)
                .map(|column| Integer::from(usize::from(row == column)))
                .collect()
        })
        .collect()
}

fn add_vectors(left: &[Integer], right: &[Integer]) -> Vec<Integer> {
    left.iter()
        .zip(right)
        .map(|(left, right)| left + right)
        .collect()
}

fn row_times_matrix(vector: &[Integer], matrix: &Matrix<IntegerRing>) -> Vec<Integer> {
    (0..matrix.ncols())
        .map(|column| {
            let column = u32::try_from(column).expect("Numerica column count fits its index type");
            vector
                .iter()
                .enumerate()
                .fold(Integer::zero(), |sum, (row, entry)| {
                    let row = u32::try_from(row).expect("Numerica row count fits its index type");
                    sum + entry * &matrix[(row, column)]
                })
        })
        .collect()
}

fn upper_triangle(matrix: &[Vec<Integer>]) -> impl Iterator<Item = &Integer> {
    matrix
        .iter()
        .enumerate()
        .flat_map(|(row, entries)| entries.iter().skip(row + 1))
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Ring,
        finite_field::Z2,
        integer::{Integer, Z},
        rational::Q,
    };
    use serde_json::{Value, json};

    use crate::{ArrowSpec, PathAlgebra, Quiver, VertexIndex};

    use super::{UnitForm, UnitFormError, tits_unit_form_of_path_algebra};

    fn form(rows: &[&[i64]]) -> UnitForm {
        UnitForm::new(
            rows.iter()
                .map(|row| row.iter().copied().map(Integer::from).collect())
                .collect(),
        )
        .unwrap()
    }

    fn rows(unit_form: &UnitForm) -> Vec<Vec<i64>> {
        unit_form
            .symmetric_matrix()
            .row_iter()
            .map(|row| row.iter().map(|entry| entry.to_i64().unwrap()).collect())
            .collect()
    }

    fn integer_vector(values: &[i64]) -> Vec<Integer> {
        values.iter().copied().map(Integer::from).collect()
    }

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn path_algebra<F: numerica::domains::Field>(
        field: F,
        arrows: Vec<ArrowSpec>,
    ) -> PathAlgebra<F> {
        PathAlgebra::new(field, Quiver::from_indexed(3, arrows).unwrap()).unwrap()
    }

    #[test]
    fn pinned_qpa_unit_form_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/combinatorialrep_unit_form.schema-v1.json"
        ))
        .unwrap();
        let a2 = form(&[&[2, -1], &[-1, 2]]);
        let a3 = form(&[&[2, -1, 0], &[-1, 2, -1], &[0, -1, 2]]);
        let affine_a1 = form(&[&[2, -2], &[-2, 2]]);
        let indefinite = form(&[&[2, -3], &[-3, 2]]);
        let nonnegative_failure = form(&[&[2, -2, -1], &[-2, 2, 0], &[-1, 0, 2]]);
        let positive_off_diagonal = form(&[&[2, 1], &[1, 2]]);

        let left = integer_vector(&[1, 2]);
        let right = integer_vector(&[-1, 3]);
        assert_eq!(
            json!({
                "matrix": rows(&a2),
                "bilinear": a2.bilinear_form(&Z, &left, &right).unwrap().to_i64().unwrap(),
                "quadratic": a2.quadratic_form(&Z, &left).unwrap().to_i64().unwrap(),
                "reflection": a2.reflection_by_bilinear_form(&Z, 2, &left).unwrap()
                    .iter().map(|entry| entry.to_i64().unwrap()).collect::<Vec<_>>(),
                "weakly_nonnegative": a2.is_weakly_nonnegative(),
                "weakly_positive": a2.is_weakly_positive(),
                "positive_roots": a2.positive_roots().unwrap().iter().map(|root|
                    root.iter().map(|entry| entry.to_i64().unwrap()).collect::<Vec<_>>()
                ).collect::<Vec<_>>(),
            }),
            expected["a2"]
        );
        assert!(a3.is_weakly_nonnegative());
        assert!(a3.is_weakly_positive());
        assert_eq!(
            json!({
                "weakly_nonnegative": true,
                "weakly_positive": true,
                "positive_roots": a3.positive_roots().unwrap().iter().map(|root|
                    root.iter().map(|entry| entry.to_i64().unwrap()).collect::<Vec<_>>()
                ).collect::<Vec<_>>(),
            }),
            expected["a3"]
        );
        assert_eq!(
            json!({
                "weakly_nonnegative": affine_a1.is_weakly_nonnegative(),
                "weakly_positive": affine_a1.is_weakly_positive(),
            }),
            expected["affine_a1"]
        );
        assert!(affine_a1.positive_roots().is_none());
        assert_eq!(
            json!({
                "weakly_nonnegative": indefinite.is_weakly_nonnegative(),
                "weakly_positive": indefinite.is_weakly_positive(),
            }),
            expected["indefinite"]
        );
        assert_eq!(
            json!({
                "weakly_nonnegative": nonnegative_failure.is_weakly_nonnegative(),
            }),
            expected["nonnegative_failure"]
        );
        assert!(positive_off_diagonal.is_weakly_positive());
        assert_eq!(
            json!({
                "weakly_positive": true,
                "positive_roots": positive_off_diagonal.positive_roots().unwrap().iter().map(
                    |root| root.iter().map(|entry| entry.to_i64().unwrap()).collect::<Vec<_>>()
                ).collect::<Vec<_>>(),
            }),
            expected["positive_off_diagonal"]
        );
    }

    #[test]
    fn forms_evaluate_over_supported_numerica_fields() {
        let a2 = form(&[&[2, -1], &[-1, 2]]);
        let rational = [Q.nth(Integer::from(1)), Q.nth(Integer::from(2))];
        assert_eq!(
            a2.quadratic_form(&Q, &rational).unwrap(),
            Q.nth(Integer::from(3))
        );
        let binary = [Z2.nth(Integer::from(1)), Z2.nth(Integer::from(0))];
        assert_eq!(a2.bilinear_form(&Z2, &binary, &binary).unwrap(), Z2.zero());
        assert_eq!(
            a2.quadratic_form(&Z2, &binary),
            Err(UnitFormError::QuadraticDivisionUnavailable)
        );
    }

    #[test]
    fn pinned_qpa_path_algebra_tits_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/combinatorialrep_tits_path.schema-v1.json"
        ))
        .unwrap();
        let a3 = path_algebra(
            Q,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        );
        let kronecker = path_algebra(
            Z2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(1), vertex(2), "y"),
            ],
        );
        let opposite_orientation = path_algebra(
            Q,
            vec![
                ArrowSpec::named(vertex(2), vertex(1), "c"),
                ArrowSpec::named(vertex(3), vertex(2), "d"),
            ],
        );

        assert_eq!(
            json!({
                "a3": rows(&tits_unit_form_of_path_algebra(&a3).unwrap()),
                "kronecker_f2": rows(&tits_unit_form_of_path_algebra(&kronecker).unwrap()),
                "opposite_orientation": rows(
                    &tits_unit_form_of_path_algebra(&opposite_orientation).unwrap()
                ),
            }),
            expected["matrices"]
        );

        let loop_algebra = path_algebra(
            Q,
            vec![ArrowSpec::named(vertex(1), vertex(1), "loop_arrow")],
        );
        assert!(matches!(
            tits_unit_form_of_path_algebra(&loop_algebra),
            Err(UnitFormError::DiagonalNotTwo { index: 0 })
        ));
    }

    #[test]
    fn constructor_and_argument_failures_are_typed() {
        assert_eq!(
            UnitForm::new(Vec::new()).unwrap_err(),
            UnitFormError::EmptyMatrix
        );
        assert!(matches!(
            UnitForm::new(vec![vec![Integer::from(2), Integer::zero()]]),
            Err(UnitFormError::NonSquare { .. })
        ));
        assert!(matches!(
            UnitForm::new(vec![vec![Integer::from(2)], vec![]]),
            Err(UnitFormError::RaggedRows { .. })
        ));
        assert!(matches!(
            UnitForm::new(vec![vec![Integer::from(1)]]),
            Err(UnitFormError::DiagonalNotTwo { index: 0 })
        ));
        assert!(matches!(
            UnitForm::new(vec![
                vec![Integer::from(2), Integer::from(1)],
                vec![Integer::zero(), Integer::from(2)]
            ]),
            Err(UnitFormError::NotSymmetric { row: 0, column: 1 })
        ));

        let a2 = form(&[&[2, -1], &[-1, 2]]);
        assert!(matches!(
            a2.bilinear_form(&Z, &integer_vector(&[1]), &integer_vector(&[1, 2])),
            Err(UnitFormError::VectorLengthsDiffer { .. })
        ));
        assert!(matches!(
            a2.quadratic_form(&Z, &integer_vector(&[1])),
            Err(UnitFormError::WrongVectorSize { .. })
        ));
        assert!(matches!(
            a2.reflection_by_bilinear_form(&Z, 0, &integer_vector(&[1, 2])),
            Err(UnitFormError::CoordinateOutsideRange { .. })
        ));
    }
}
