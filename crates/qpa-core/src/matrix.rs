//! Exact, explicitly shaped matrices using GAP's row-vector conventions.
//!
//! The elimination routines in this module translate the generic GAP methods
//! `DoSemiEchelonMatTransformationDestructive`, `NullspaceMat`, and
//! `SolutionMatDestructive` from `lib/matrix.gi` in the pinned GAP source.

use std::ops::{Index, IndexMut};

use numerica::domains::Field;
use thiserror::Error;

/// A dense row-major matrix over an exact Numerica field.
///
/// Row and column counts are stored independently of the data, so shapes such
/// as `0 x n` and `n x 0` are not collapsed.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct DenseMatrix<F: Field> {
    nrows: usize,
    ncols: usize,
    data: Vec<F::Element>,
    field: F,
}

impl<F: Field> DenseMatrix<F> {
    /// Construct a zero matrix with the requested explicit shape.
    #[must_use]
    pub fn zero(nrows: usize, ncols: usize, field: F) -> Self {
        Self {
            nrows,
            ncols,
            data: vec![field.zero(); nrows.saturating_mul(ncols)],
            field,
        }
    }

    /// Construct a row-major matrix while preserving its explicit shape.
    ///
    /// # Errors
    ///
    /// Returns [`MatrixError::ShapeOverflow`] if the dimensions overflow or
    /// [`MatrixError::DataLength`] if `data` does not match the shape.
    pub fn from_flat(
        nrows: usize,
        ncols: usize,
        data: Vec<F::Element>,
        field: F,
    ) -> Result<Self, MatrixError> {
        let expected = nrows
            .checked_mul(ncols)
            .ok_or(MatrixError::ShapeOverflow { nrows, ncols })?;
        if data.len() != expected {
            return Err(MatrixError::DataLength {
                actual: data.len(),
                expected,
                nrows,
                ncols,
            });
        }

        Ok(Self {
            nrows,
            ncols,
            data,
            field,
        })
    }

    /// Construct a matrix from rows. An empty row list has shape `0 x 0`; use
    /// [`Self::from_flat`] when an explicit `0 x n` shape is required.
    ///
    /// # Errors
    ///
    /// Returns [`MatrixError::RaggedRows`] if the row lengths differ.
    pub fn from_rows(rows: Vec<Vec<F::Element>>, field: F) -> Result<Self, MatrixError> {
        let nrows = rows.len();
        let ncols = rows.first().map_or(0, Vec::len);
        if let Some((row, actual)) = rows
            .iter()
            .enumerate()
            .find_map(|(index, row)| (row.len() != ncols).then_some((index, row.len())))
        {
            return Err(MatrixError::RaggedRows {
                row,
                actual,
                expected: ncols,
            });
        }

        Self::from_flat(nrows, ncols, rows.into_iter().flatten().collect(), field)
    }

    /// Construct an identity matrix.
    #[must_use]
    pub fn identity(size: usize, field: F) -> Self {
        let mut result = Self::zero(size, size, field);
        for index in 0..size {
            result[(index, index)] = result.field.one();
        }
        result
    }

    /// Return the number of rows.
    #[must_use]
    pub const fn nrows(&self) -> usize {
        self.nrows
    }

    /// Return the number of columns.
    #[must_use]
    pub const fn ncols(&self) -> usize {
        self.ncols
    }

    /// Return the exact coefficient field.
    #[must_use]
    pub const fn field(&self) -> &F {
        &self.field
    }

    /// Return the row-major entries.
    #[must_use]
    pub fn as_flat(&self) -> &[F::Element] {
        &self.data
    }

    /// Consume the matrix and return its row-major entries.
    #[must_use]
    pub fn into_flat(self) -> Vec<F::Element> {
        self.data
    }

    /// Borrow one row, including rows of a matrix with zero columns.
    ///
    /// # Panics
    ///
    /// Panics when `row` is outside the stored shape.
    #[must_use]
    pub fn row(&self, row: usize) -> &[F::Element] {
        assert!(row < self.nrows, "matrix row index out of bounds");
        let start = row * self.ncols;
        &self.data[start..start + self.ncols]
    }

    fn row_mut(&mut self, row: usize) -> &mut [F::Element] {
        assert!(row < self.nrows, "matrix row index out of bounds");
        let start = row * self.ncols;
        &mut self.data[start..start + self.ncols]
    }

    /// Return whether every matrix entry is zero.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.data.iter().all(|entry| self.field.is_zero(entry))
    }

    /// Translate GAP's semi-echelon transformation algorithm.
    ///
    /// The returned `relations` rows are exactly the coefficient relations
    /// among the input rows and therefore form GAP's left nullspace basis.
    #[must_use]
    pub fn semi_echelon_transformation(&self) -> SemiEchelonTransformation<F> {
        let mut matrix = self.clone();
        let mut transform = Self::identity(self.nrows, self.field.clone());
        let mut heads: Vec<Option<usize>> = vec![None; self.ncols];
        let mut vectors: Vec<Vec<F::Element>> = Vec::new();
        let mut coefficients: Vec<Vec<F::Element>> = Vec::new();
        let mut relations: Vec<Vec<F::Element>> = Vec::new();

        for row_index in 0..self.nrows {
            for column in 0..self.ncols {
                let Some(head) = heads[column] else {
                    continue;
                };
                let scale = self.field.neg(&matrix[(row_index, column)]);
                if !self.field.is_zero(&scale) {
                    add_scaled_row(
                        &self.field,
                        transform.row_mut(row_index),
                        coefficients[head].as_slice(),
                        &scale,
                    );
                    add_scaled_row(
                        &self.field,
                        matrix.row_mut(row_index),
                        vectors[head].as_slice(),
                        &scale,
                    );
                }
            }

            let pivot = matrix
                .row(row_index)
                .iter()
                .position(|entry| !self.field.is_zero(entry));

            if let Some(column) = pivot {
                let inverse = self.field.inv(&matrix[(row_index, column)]);
                scale_row(&self.field, transform.row_mut(row_index), &inverse);
                scale_row(&self.field, matrix.row_mut(row_index), &inverse);

                coefficients.push(transform.row(row_index).to_vec());
                vectors.push(matrix.row(row_index).to_vec());
                heads[column] = Some(vectors.len() - 1);
            } else {
                relations.push(transform.row(row_index).to_vec());
            }
        }

        SemiEchelonTransformation {
            heads,
            vectors: rows_with_shape(vectors, self.ncols, self.field.clone()),
            coefficients: rows_with_shape(coefficients, self.nrows, self.field.clone()),
            relations: rows_with_shape(relations, self.nrows, self.field.clone()),
        }
    }

    /// Return GAP's left nullspace: rows `x` satisfying `x * self = 0`.
    #[must_use]
    pub fn nullspace_rows(&self) -> Self {
        self.semi_echelon_transformation().relations
    }

    /// Return one GAP-compatible row solution `x` of `x * self = target`.
    ///
    /// `Ok(None)` corresponds to GAP's `fail` result.
    ///
    /// # Errors
    ///
    /// Returns [`MatrixError::TargetLength`] when `target` does not have one
    /// entry per matrix column.
    pub fn solution_row(
        &self,
        target: &[F::Element],
    ) -> Result<Option<Vec<F::Element>>, MatrixError> {
        if target.len() != self.ncols {
            return Err(MatrixError::TargetLength {
                actual: target.len(),
                expected: self.ncols,
            });
        }

        let semi = self.semi_echelon_transformation();
        let mut residual = target.to_vec();
        let mut solution = vec![self.field.zero(); self.nrows];

        for column in 0..self.ncols {
            let Some(vector_number) = semi.heads[column] else {
                continue;
            };
            let scale = residual[column].clone();
            if !self.field.is_zero(&scale) {
                let negative_scale = self.field.neg(&scale);
                add_scaled_row(
                    &self.field,
                    &mut residual,
                    semi.vectors.row(vector_number),
                    &negative_scale,
                );
                add_scaled_row(
                    &self.field,
                    &mut solution,
                    semi.coefficients.row(vector_number),
                    &scale,
                );
            }
        }

        Ok(residual
            .iter()
            .all(|entry| self.field.is_zero(entry))
            .then_some(solution))
    }

    /// Return the row rank computed by GAP's generic semi-echelon algorithm.
    #[must_use]
    pub fn rank(&self) -> usize {
        self.semi_echelon_transformation().vectors.nrows
    }
}

impl<F: Field> Index<(usize, usize)> for DenseMatrix<F> {
    type Output = F::Element;

    fn index(&self, (row, column): (usize, usize)) -> &Self::Output {
        assert!(row < self.nrows, "matrix row index out of bounds");
        assert!(column < self.ncols, "matrix column index out of bounds");
        &self.data[row * self.ncols + column]
    }
}

impl<F: Field> IndexMut<(usize, usize)> for DenseMatrix<F> {
    fn index_mut(&mut self, (row, column): (usize, usize)) -> &mut Self::Output {
        assert!(row < self.nrows, "matrix row index out of bounds");
        assert!(column < self.ncols, "matrix column index out of bounds");
        &mut self.data[row * self.ncols + column]
    }
}

/// Output of GAP's semi-echelon transformation algorithm.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct SemiEchelonTransformation<F: Field> {
    /// For each input column, the corresponding basis-vector index, if any.
    pub heads: Vec<Option<usize>>,
    /// Normalized semi-echelon basis rows.
    pub vectors: DenseMatrix<F>,
    /// Coefficients expressing each basis row in terms of the input rows.
    pub coefficients: DenseMatrix<F>,
    /// Coefficient rows expressing dependencies among the input rows.
    pub relations: DenseMatrix<F>,
}

/// Failure to construct or operate on an explicitly shaped matrix.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum MatrixError {
    /// The requested dimensions overflow `usize` storage.
    #[error("matrix shape {nrows} x {ncols} overflows addressable storage")]
    ShapeOverflow {
        /// Requested row count.
        nrows: usize,
        /// Requested column count.
        ncols: usize,
    },
    /// Flat data does not match the requested shape.
    #[error("matrix data has length {actual}, expected {expected} for shape {nrows} x {ncols}")]
    DataLength {
        /// Actual number of entries.
        actual: usize,
        /// Required number of entries.
        expected: usize,
        /// Requested row count.
        nrows: usize,
        /// Requested column count.
        ncols: usize,
    },
    /// A row-list constructor received inconsistent row lengths.
    #[error("matrix row {row} has length {actual}, expected {expected}")]
    RaggedRows {
        /// Zero-based row index.
        row: usize,
        /// Actual row length.
        actual: usize,
        /// Expected row length.
        expected: usize,
    },
    /// A target for `x * matrix = target` has the wrong length.
    #[error("solution target has length {actual}, expected {expected}")]
    TargetLength {
        /// Actual target length.
        actual: usize,
        /// Required target length.
        expected: usize,
    },
}

fn add_scaled_row<F: Field>(
    field: &F,
    target: &mut [F::Element],
    source: &[F::Element],
    scale: &F::Element,
) {
    debug_assert_eq!(target.len(), source.len());
    for (target_entry, source_entry) in target.iter_mut().zip(source) {
        let product = field.mul(source_entry, scale);
        field.add_assign(target_entry, &product);
    }
}

fn scale_row<F: Field>(field: &F, row: &mut [F::Element], scale: &F::Element) {
    for entry in row {
        field.mul_assign(entry, scale);
    }
}

fn rows_with_shape<F: Field>(rows: Vec<Vec<F::Element>>, ncols: usize, field: F) -> DenseMatrix<F> {
    let nrows = rows.len();
    DenseMatrix::from_flat(nrows, ncols, rows.into_iter().flatten().collect(), field)
        .expect("internally generated matrix rows have a consistent shape")
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Ring,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, Rational},
    };

    use super::DenseMatrix;

    fn q(numerator: i64, denominator: i64) -> Rational {
        Q.to_element(numerator.into(), denominator.into(), true)
    }

    #[test]
    fn gap_rational_nullspace_example() {
        let matrix = DenseMatrix::from_rows(
            vec![
                vec![q(1, 1), q(2, 1), q(3, 1)],
                vec![q(4, 1), q(5, 1), q(6, 1)],
                vec![q(7, 1), q(8, 1), q(9, 1)],
            ],
            Q,
        )
        .unwrap();

        let nullspace = matrix.nullspace_rows();
        assert_eq!(nullspace.nrows(), 1);
        assert_eq!(nullspace.ncols(), 3);
        assert_eq!(nullspace.row(0), &[q(1, 1), q(-2, 1), q(1, 1)]);
    }

    #[test]
    fn gap_rational_solution_example() {
        let matrix = DenseMatrix::from_rows(
            vec![
                vec![q(1, 1), q(2, 1), q(3, 1)],
                vec![q(4, 1), q(5, 1), q(6, 1)],
                vec![q(7, 1), q(8, 1), q(9, 1)],
            ],
            Q,
        )
        .unwrap();

        let solution = matrix
            .solution_row(&[q(3, 1), q(5, 1), q(7, 1)])
            .unwrap()
            .unwrap();
        assert_eq!(solution, vec![q(5, 3), q(1, 3), q(0, 1)]);
    }

    #[test]
    fn prime_field_three_uses_the_same_source_algorithm() {
        let field = Zp::new(3);
        let e = |value| field.to_element(value);
        let matrix =
            DenseMatrix::from_rows(vec![vec![e(1), e(1)], vec![e(2), e(2)]], field.clone())
                .unwrap();

        let nullspace = matrix.nullspace_rows();
        assert_eq!(nullspace.row(0), &[e(1), e(1)]);
        assert_eq!(matrix.rank(), 1);
    }

    #[test]
    fn characteristic_two_special_backend_is_supported() {
        let e = |value: i64| Z2.nth(Integer::from(value));
        let matrix = DenseMatrix::from_rows(vec![vec![e(1), e(1)], vec![e(1), e(1)]], Z2).unwrap();

        let nullspace = matrix.nullspace_rows();
        assert_eq!(nullspace.row(0), &[e(1), e(1)]);
        assert_eq!(matrix.rank(), 1);
    }

    #[test]
    fn zero_column_shapes_preserve_gap_left_nullspace_semantics() {
        let matrix = DenseMatrix::from_flat(3, 0, Vec::new(), Q).unwrap();
        let nullspace = matrix.nullspace_rows();

        assert_eq!((nullspace.nrows(), nullspace.ncols()), (3, 3));
        assert_eq!(nullspace, DenseMatrix::identity(3, Q));
        assert_eq!(matrix.solution_row(&[]).unwrap(), Some(vec![q(0, 1); 3]));
    }

    #[test]
    fn zero_row_shapes_remain_explicit() {
        let matrix = DenseMatrix::from_flat(0, 4, Vec::new(), Q).unwrap();
        let nullspace = matrix.nullspace_rows();

        assert_eq!((matrix.nrows(), matrix.ncols()), (0, 4));
        assert_eq!((nullspace.nrows(), nullspace.ncols()), (0, 0));
        assert_eq!(
            matrix.solution_row(&vec![q(0, 1); 4]).unwrap(),
            Some(vec![])
        );
        assert_eq!(matrix.solution_row(&vec![q(1, 1); 4]).unwrap(), None);
    }
}
