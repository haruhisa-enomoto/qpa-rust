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

    /// Add equally shaped matrices over the same exact field.
    ///
    /// # Errors
    ///
    /// Rejects different fields or shapes.
    pub fn add(&self, right: &Self) -> Result<Self, MatrixError> {
        self.ensure_same_field(right)?;
        self.ensure_same_shape(right)?;
        Self::from_flat(
            self.nrows,
            self.ncols,
            self.data
                .iter()
                .zip(&right.data)
                .map(|(left, right)| self.field.add(left, right))
                .collect(),
            self.field.clone(),
        )
    }

    /// Multiply matrices using GAP's row-vector convention.
    ///
    /// # Errors
    ///
    /// Rejects different fields or incompatible inner dimensions.
    pub fn multiply(&self, right: &Self) -> Result<Self, MatrixError> {
        self.ensure_same_field(right)?;
        if self.ncols != right.nrows {
            return Err(MatrixError::IncompatibleProduct {
                left: (self.nrows, self.ncols),
                right: (right.nrows, right.ncols),
            });
        }
        let mut result = Self::zero(self.nrows, right.ncols, self.field.clone());
        for row in 0..self.nrows {
            for inner in 0..self.ncols {
                let left = &self[(row, inner)];
                if self.field.is_zero(left) {
                    continue;
                }
                for column in 0..right.ncols {
                    let product = self.field.mul(left, &right[(inner, column)]);
                    self.field.add_assign(&mut result[(row, column)], &product);
                }
            }
        }
        Ok(result)
    }

    /// Multiply every entry by a field scalar on the left.
    #[must_use]
    pub fn scale_left(&self, scalar: &F::Element) -> Self {
        Self {
            nrows: self.nrows,
            ncols: self.ncols,
            data: self
                .data
                .iter()
                .map(|entry| self.field.mul(scalar, entry))
                .collect(),
            field: self.field.clone(),
        }
    }

    /// Return the matrix transpose with its explicit shape reversed.
    #[must_use]
    pub fn transpose(&self) -> Self {
        let mut result = Self::zero(self.ncols, self.nrows, self.field.clone());
        for row in 0..self.nrows {
            for column in 0..self.ncols {
                result[(column, row)] = self[(row, column)].clone();
            }
        }
        result
    }

    fn ensure_same_field(&self, right: &Self) -> Result<(), MatrixError> {
        if self.field == right.field {
            Ok(())
        } else {
            Err(MatrixError::DifferentFields)
        }
    }

    fn ensure_same_shape(&self, right: &Self) -> Result<(), MatrixError> {
        if (self.nrows, self.ncols) == (right.nrows, right.ncols) {
            Ok(())
        } else {
            Err(MatrixError::DifferentShapes {
                left: (self.nrows, self.ncols),
                right: (right.nrows, right.ncols),
            })
        }
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

    /// Return GAP's canonical basis for the row space of this matrix.
    ///
    /// This translates the generic `TriangulizeMat` method used by
    /// `CanonicalBasis` for a Gaussian row space.  Unlike
    /// [`Self::semi_echelon_transformation`], every pivot column is cleared
    /// both above and below its pivot, and the surviving rows are ordered by
    /// increasing pivot column.
    #[must_use]
    pub fn canonical_row_space_basis(&self) -> Self {
        let mut matrix = self.clone();
        let mut rank = 0;

        for column in 0..self.ncols {
            let Some(pivot_row) =
                (rank..self.nrows).find(|&row| !self.field.is_zero(&matrix[(row, column)]))
            else {
                continue;
            };

            if pivot_row != rank {
                for entry_column in 0..self.ncols {
                    let pivot_index = pivot_row * self.ncols + entry_column;
                    let rank_index = rank * self.ncols + entry_column;
                    matrix.data.swap(pivot_index, rank_index);
                }
            }

            let inverse = self.field.inv(&matrix[(rank, column)]);
            scale_row(&self.field, matrix.row_mut(rank), &inverse);
            let pivot = matrix.row(rank).to_vec();
            for row in 0..self.nrows {
                if row == rank {
                    continue;
                }
                let scale = self.field.neg(&matrix[(row, column)]);
                if !self.field.is_zero(&scale) {
                    add_scaled_row(&self.field, matrix.row_mut(row), &pivot, &scale);
                }
            }
            rank += 1;
            if rank == self.nrows {
                break;
            }
        }

        Self {
            nrows: rank,
            ncols: self.ncols,
            data: matrix.data[..rank * self.ncols].to_vec(),
            field: self.field.clone(),
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

    /// Return the exact two-sided inverse of a square matrix, if it exists.
    ///
    /// Rows are obtained by GAP-compatible solutions of `x * self = e_i`, so
    /// the result uses the same row-action convention as QPA module maps.
    /// `Ok(None)` represents a singular square matrix.
    ///
    /// # Errors
    ///
    /// Returns [`MatrixError::NonSquareInverse`] for a nonsquare matrix and
    /// propagates checked matrix construction failures.
    pub fn inverse(&self) -> Result<Option<Self>, MatrixError> {
        if self.nrows != self.ncols {
            return Err(MatrixError::NonSquareInverse {
                shape: (self.nrows, self.ncols),
            });
        }
        let mut entries = Vec::with_capacity(self.data.len());
        for row in 0..self.nrows {
            let mut target = vec![self.field.zero(); self.ncols];
            target[row] = self.field.one();
            let Some(solution) = self.solution_row(&target)? else {
                return Ok(None);
            };
            entries.extend(solution);
        }
        Ok(Some(Self::from_flat(
            self.nrows,
            self.ncols,
            entries,
            self.field.clone(),
        )?))
    }

    /// Translate GAP's `Matrix_CharacteristicPolynomialSameField` spinning
    /// algorithm.
    ///
    /// Coefficients are returned in ascending degree order, exactly as
    /// `CoefficientsOfUnivariatePolynomial(CharacteristicPolynomial(...))`.
    /// The source algorithm builds cyclic subspaces from successive standard
    /// rows and multiplies their order polynomials; it does not use traces or
    /// division by coefficient-field integers.
    ///
    /// # Errors
    ///
    /// Rejects a nonsquare matrix.
    pub fn characteristic_polynomial_coefficients(&self) -> Result<Vec<F::Element>, MatrixError> {
        if self.nrows != self.ncols {
            return Err(MatrixError::NonSquareCharacteristicPolynomial {
                shape: (self.nrows, self.ncols),
            });
        }
        if self.nrows == 0 {
            return Ok(vec![self.field.one()]);
        }

        let mut basis = vec![None; self.nrows];
        let mut vector = vec![self.field.zero(); self.nrows];
        let mut characteristic = vec![self.field.one()];
        for coordinate in 0..self.nrows {
            if basis[coordinate].is_some() {
                continue;
            }
            vector[coordinate] = self.field.one();
            let order = self.order_polynomial_inner(&vector, &mut basis);
            characteristic = multiply_polynomial_coefficients(&self.field, &characteristic, &order);
            vector[coordinate] = self.field.zero();
        }
        Ok(characteristic)
    }

    fn order_polynomial_inner(
        &self,
        initial: &[F::Element],
        spanned: &mut [Option<Vec<F::Element>>],
    ) -> Vec<F::Element> {
        let dimension = initial.len();
        let mut vector = initial.to_vec();
        let mut polynomials: Vec<Option<Vec<F::Element>>> = vec![None; dimension];
        let mut zeroes = Vec::new();

        loop {
            let mut remainder = vector.clone();
            let mut polynomial = zeroes.clone();
            polynomial.push(self.field.one());
            let mut pivot = first_nonzero_from(&self.field, &remainder, 0);

            while let Some(position) = pivot {
                let Some(base) = &spanned[position] else {
                    break;
                };
                let scalar = self.field.neg(&remainder[position]);
                if let Some(reducer) = &polynomials[position] {
                    add_scaled_polynomial(&self.field, &mut polynomial, reducer, &scalar);
                }
                add_scaled_row(&self.field, &mut remainder, base, &scalar);
                pivot = first_nonzero_from(&self.field, &remainder, position + 1);
            }

            let Some(position) = pivot else {
                return polynomial;
            };
            let inverse = self.field.inv(&remainder[position]);
            for coefficient in &mut polynomial {
                self.field.mul_assign(coefficient, &inverse);
            }
            for coefficient in &mut remainder {
                self.field.mul_assign(coefficient, &inverse);
            }
            polynomials[position] = Some(polynomial);
            spanned[position] = Some(remainder);
            vector = multiply_row_by_matrix(&self.field, &vector, self);
            zeroes.push(self.field.zero());
        }
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
    /// Matrix arithmetic requires one exact coefficient field.
    #[error("matrices are over different coefficient fields")]
    DifferentFields,
    /// Matrix addition requires equal shapes.
    #[error("matrix shapes {left:?} and {right:?} are different")]
    DifferentShapes {
        /// Left matrix shape.
        left: (usize, usize),
        /// Right matrix shape.
        right: (usize, usize),
    },
    /// Matrix multiplication requires matching inner dimensions.
    #[error("cannot multiply matrices with shapes {left:?} and {right:?}")]
    IncompatibleProduct {
        /// Left matrix shape.
        left: (usize, usize),
        /// Right matrix shape.
        right: (usize, usize),
    },
    /// Only square matrices can have a two-sided inverse.
    #[error("cannot invert a nonsquare matrix of shape {shape:?}")]
    NonSquareInverse {
        /// Entered row and column counts.
        shape: (usize, usize),
    },
    /// Characteristic polynomials require square matrices.
    #[error("cannot compute a characteristic polynomial for a {shape:?} matrix")]
    NonSquareCharacteristicPolynomial {
        /// Rejected matrix shape.
        shape: (usize, usize),
    },
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

fn first_nonzero_from<F: Field>(field: &F, row: &[F::Element], start: usize) -> Option<usize> {
    row.iter()
        .enumerate()
        .skip(start)
        .find_map(|(index, entry)| (!field.is_zero(entry)).then_some(index))
}

fn add_scaled_polynomial<F: Field>(
    field: &F,
    target: &mut Vec<F::Element>,
    source: &[F::Element],
    scale: &F::Element,
) {
    if target.len() < source.len() {
        target.resize(source.len(), field.zero());
    }
    for (target_coefficient, source_coefficient) in target.iter_mut().zip(source) {
        field.add_mul_assign(target_coefficient, source_coefficient, scale);
    }
}

fn multiply_polynomial_coefficients<F: Field>(
    field: &F,
    left: &[F::Element],
    right: &[F::Element],
) -> Vec<F::Element> {
    if left.is_empty() || right.is_empty() {
        return Vec::new();
    }
    let mut result = vec![field.zero(); left.len() + right.len() - 1];
    for (left_degree, left_coefficient) in left.iter().enumerate() {
        for (right_degree, right_coefficient) in right.iter().enumerate() {
            let product = field.mul(left_coefficient, right_coefficient);
            field.add_assign(&mut result[left_degree + right_degree], &product);
        }
    }
    result
}

fn multiply_row_by_matrix<F: Field>(
    field: &F,
    row: &[F::Element],
    matrix: &DenseMatrix<F>,
) -> Vec<F::Element> {
    debug_assert_eq!(row.len(), matrix.nrows());
    let mut result = vec![field.zero(); matrix.ncols()];
    for (inner, coefficient) in row.iter().enumerate() {
        if field.is_zero(coefficient) {
            continue;
        }
        for column in 0..matrix.ncols() {
            field.add_mul_assign(&mut result[column], coefficient, &matrix[(inner, column)]);
        }
    }
    result
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

    use super::{DenseMatrix, MatrixError};

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
    fn canonical_row_basis_translates_gap_full_triangulation() {
        let matrix = DenseMatrix::from_rows(
            vec![
                vec![q(1, 1), q(1, 1), q(0, 1), q(1, 1)],
                vec![q(0, 1), q(1, 1), q(1, 1), q(0, 1)],
            ],
            Q,
        )
        .unwrap();

        assert_eq!(
            matrix.canonical_row_space_basis(),
            DenseMatrix::from_rows(
                vec![
                    vec![q(1, 1), q(0, 1), q(-1, 1), q(1, 1)],
                    vec![q(0, 1), q(1, 1), q(1, 1), q(0, 1)],
                ],
                Q,
            )
            .unwrap()
        );
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

    #[test]
    fn arithmetic_uses_gap_row_action_order() {
        let first = DenseMatrix::from_rows(vec![vec![q(1, 1), q(2, 1)], vec![q(0, 1), q(1, 1)]], Q)
            .unwrap();
        let second =
            DenseMatrix::from_rows(vec![vec![q(3, 1), q(0, 1)], vec![q(4, 1), q(5, 1)]], Q)
                .unwrap();

        assert_eq!(
            first.multiply(&second).unwrap(),
            DenseMatrix::from_rows(vec![vec![q(11, 1), q(10, 1)], vec![q(4, 1), q(5, 1)]], Q,)
                .unwrap()
        );
        assert_eq!(
            first.add(&second).unwrap(),
            DenseMatrix::from_rows(vec![vec![q(4, 1), q(2, 1)], vec![q(4, 1), q(6, 1)]], Q,)
                .unwrap()
        );
        assert_eq!(
            first.scale_left(&q(-2, 1)),
            DenseMatrix::from_rows(vec![vec![q(-2, 1), q(-4, 1)], vec![q(0, 1), q(-2, 1)]], Q,)
                .unwrap()
        );
        assert_eq!(
            first.transpose(),
            DenseMatrix::from_rows(vec![vec![q(1, 1), q(0, 1)], vec![q(2, 1), q(1, 1)]], Q,)
                .unwrap()
        );
    }

    #[test]
    fn exact_inverse_uses_gap_row_solutions() {
        let matrix =
            DenseMatrix::from_rows(vec![vec![q(1, 1), q(2, 1)], vec![q(3, 1), q(5, 1)]], Q)
                .unwrap();
        let inverse = matrix.inverse().unwrap().unwrap();

        assert_eq!(
            inverse,
            DenseMatrix::from_rows(vec![vec![q(-5, 1), q(2, 1)], vec![q(3, 1), q(-1, 1)]], Q,)
                .unwrap()
        );
        assert_eq!(
            matrix.multiply(&inverse).unwrap(),
            DenseMatrix::identity(2, Q)
        );
        assert!(DenseMatrix::zero(2, 2, Q).inverse().unwrap().is_none());
        assert_eq!(
            DenseMatrix::zero(1, 2, Q).inverse().unwrap_err(),
            MatrixError::NonSquareInverse { shape: (1, 2) }
        );
    }

    #[test]
    fn characteristic_polynomial_uses_gap_spinning_order() {
        let jordan =
            DenseMatrix::from_rows(vec![vec![q(1, 1), q(1, 1)], vec![q(0, 1), q(1, 1)]], Q)
                .unwrap();
        assert_eq!(
            jordan.characteristic_polynomial_coefficients().unwrap(),
            vec![q(1, 1), q(-2, 1), q(1, 1)]
        );

        let diagonal = DenseMatrix::from_rows(
            vec![
                vec![q(1, 1), q(0, 1), q(0, 1)],
                vec![q(0, 1), q(2, 1), q(0, 1)],
                vec![q(0, 1), q(0, 1), q(3, 1)],
            ],
            Q,
        )
        .unwrap();
        assert_eq!(
            diagonal.characteristic_polynomial_coefficients().unwrap(),
            vec![q(-6, 1), q(11, 1), q(-6, 1), q(1, 1)]
        );

        let binary = DenseMatrix::from_rows(
            vec![vec![Z2.zero(), Z2.one()], vec![Z2.one(), Z2.one()]],
            Z2,
        )
        .unwrap();
        assert_eq!(
            binary.characteristic_polynomial_coefficients().unwrap(),
            vec![Z2.one(), Z2.one(), Z2.one()]
        );
    }

    #[test]
    fn arithmetic_preserves_explicit_empty_shapes() {
        let zero_by_three = DenseMatrix::from_flat(0, 3, Vec::new(), Q).unwrap();
        let three_by_two = DenseMatrix::zero(3, 2, Q);
        let product = zero_by_three.multiply(&three_by_two).unwrap();

        assert_eq!((product.nrows(), product.ncols()), (0, 2));
        assert_eq!(
            zero_by_three.transpose(),
            DenseMatrix::from_flat(3, 0, Vec::new(), Q).unwrap()
        );
        assert_eq!(
            DenseMatrix::zero(2, 3, Q)
                .multiply(&DenseMatrix::zero(4, 1, Q))
                .unwrap_err(),
            MatrixError::IncompatibleProduct {
                left: (2, 3),
                right: (4, 1),
            }
        );
    }

    #[test]
    fn arithmetic_rejects_different_prime_fields() {
        let left = DenseMatrix::zero(1, 1, Zp::new(3));
        let right = DenseMatrix::zero(1, 1, Zp::new(5));

        assert_eq!(left.add(&right), Err(MatrixError::DifferentFields));
        assert_eq!(left.multiply(&right), Err(MatrixError::DifferentFields));
    }
}
