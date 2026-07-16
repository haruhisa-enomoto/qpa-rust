//! Finite-dimensional algebras given by structure constants.
//!
//! QPA constructs these algebras through GAP's `AlgebraByStructureConstants`.
//! The ordered standard basis and row-vector coefficient conventions are
//! observable in QPA's higher preprojective-algebra implementation.

use std::sync::Arc;

use numerica::domains::Field;
use thiserror::Error;

use crate::{DenseMatrix, MatrixError};

#[derive(Debug)]
struct StructureConstantsAlgebraData<F: Field> {
    field: F,
    dimension: usize,
    products: Vec<F::Element>,
}

/// A finite-dimensional associative algebra in its ordered standard basis.
///
/// The product at `(left, right)` is stored as a full coordinate row in the
/// same format produced by GAP's `AlgebraByStructureConstants`.
#[derive(Clone, Debug)]
pub struct StructureConstantsAlgebra<F: Field>(Arc<StructureConstantsAlgebraData<F>>);

impl<F: Field> StructureConstantsAlgebra<F> {
    /// Construct an algebra from a flattened `dimension × dimension ×
    /// dimension` multiplication tensor.
    ///
    /// # Errors
    ///
    /// Rejects an overflowing tensor size or a coefficient count different
    /// from the requested cube.
    pub fn new(
        field: F,
        dimension: usize,
        products: Vec<F::Element>,
    ) -> Result<Self, StructureConstantsError> {
        let expected = dimension
            .checked_mul(dimension)
            .and_then(|square| square.checked_mul(dimension))
            .ok_or(StructureConstantsError::TensorShapeOverflow { dimension })?;
        if products.len() != expected {
            return Err(StructureConstantsError::TensorLength {
                actual: products.len(),
                expected,
                dimension,
            });
        }
        Ok(Self(Arc::new(StructureConstantsAlgebraData {
            field,
            dimension,
            products,
        })))
    }

    /// Return the exact coefficient field.
    #[must_use]
    pub fn field(&self) -> &F {
        &self.0.field
    }

    /// Return the dimension of the ordered standard basis.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.0.dimension
    }

    /// Return the full coordinate row of a standard-basis product.
    ///
    /// # Panics
    ///
    /// Panics if either basis index is outside the algebra dimension.
    #[must_use]
    pub fn basis_product(&self, left: usize, right: usize) -> &[F::Element] {
        assert!(left < self.dimension(), "left basis index out of bounds");
        assert!(right < self.dimension(), "right basis index out of bounds");
        let start = (left * self.dimension() + right) * self.dimension();
        &self.0.products[start..start + self.dimension()]
    }

    /// Return the additive zero element.
    #[must_use]
    pub fn zero(&self) -> StructureConstantsAlgebraElement<F> {
        StructureConstantsAlgebraElement {
            algebra: self.clone(),
            coordinates: vec![self.field().zero(); self.dimension()],
        }
    }

    /// Return `BasisVectors(Basis(A))` in GAP's standard-basis order.
    #[must_use]
    pub fn basis_elements(&self) -> Vec<StructureConstantsAlgebraElement<F>> {
        (0..self.dimension())
            .map(|index| {
                let mut element = self.zero();
                element.coordinates[index] = self.field().one();
                element
            })
            .collect()
    }

    /// Construct an element from standard-basis coordinates.
    ///
    /// # Errors
    ///
    /// Rejects a coordinate row whose length differs from the dimension.
    pub fn element(
        &self,
        coordinates: Vec<F::Element>,
    ) -> Result<StructureConstantsAlgebraElement<F>, StructureConstantsError> {
        if coordinates.len() != self.dimension() {
            return Err(StructureConstantsError::ElementLength {
                actual: coordinates.len(),
                expected: self.dimension(),
            });
        }
        Ok(StructureConstantsAlgebraElement {
            algebra: self.clone(),
            coordinates,
        })
    }

    /// Translate QPA's `FindMultiplicativeIdentity(A)` method.
    ///
    /// QPA concatenates the matrices of right multiplication by every basis
    /// vector and solves one row system. Consequently this returns the same
    /// first GAP-compatible left identity, or `None` for GAP's `fail`.
    ///
    /// # Errors
    ///
    /// Propagates checked matrix construction and solution failures.
    pub fn find_multiplicative_identity(
        &self,
    ) -> Result<Option<StructureConstantsAlgebraElement<F>>, StructureConstantsError> {
        let dimension = self.dimension();
        if dimension == 0 {
            return Err(StructureConstantsError::EmptyIdentityBlockMatrix);
        }
        let columns = dimension
            .checked_mul(dimension)
            .ok_or(StructureConstantsError::IdentityMatrixShapeOverflow { dimension })?;
        let mut matrix = Vec::with_capacity(
            dimension
                .checked_mul(columns)
                .ok_or(StructureConstantsError::IdentityMatrixShapeOverflow { dimension })?,
        );
        for unknown in 0..dimension {
            for right in 0..dimension {
                matrix.extend_from_slice(self.basis_product(unknown, right));
            }
        }
        let total_matrix =
            DenseMatrix::from_flat(dimension, columns, matrix, self.field().clone())?;
        let mut target = Vec::with_capacity(columns);
        for basis_index in 0..dimension {
            for coordinate in 0..dimension {
                target.push(if basis_index == coordinate {
                    self.field().one()
                } else {
                    self.field().zero()
                });
            }
        }
        total_matrix
            .solution_row(&target)?
            .map(|coordinates| self.element(coordinates))
            .transpose()
    }
}

impl<F: Field> PartialEq for StructureConstantsAlgebra<F> {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl<F: Field> Eq for StructureConstantsAlgebra<F> {}

/// An element of a [`StructureConstantsAlgebra`].
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StructureConstantsAlgebraElement<F: Field> {
    algebra: StructureConstantsAlgebra<F>,
    coordinates: Vec<F::Element>,
}

impl<F: Field> StructureConstantsAlgebraElement<F> {
    /// Return the parent algebra.
    #[must_use]
    pub const fn algebra(&self) -> &StructureConstantsAlgebra<F> {
        &self.algebra
    }

    /// Return standard-basis coordinates.
    #[must_use]
    pub fn coordinates(&self) -> &[F::Element] {
        &self.coordinates
    }

    /// Multiply elements using the stored structure constants.
    ///
    /// # Errors
    ///
    /// Rejects elements from distinct parent objects.
    pub fn multiply(&self, right: &Self) -> Result<Self, StructureConstantsError> {
        if self.algebra != right.algebra {
            return Err(StructureConstantsError::DifferentParents);
        }
        let field = self.algebra.field();
        let mut coordinates = vec![field.zero(); self.algebra.dimension()];
        for (left_index, left_coefficient) in self.coordinates.iter().enumerate() {
            if field.is_zero(left_coefficient) {
                continue;
            }
            for (right_index, right_coefficient) in right.coordinates.iter().enumerate() {
                if field.is_zero(right_coefficient) {
                    continue;
                }
                let scalar = field.mul(left_coefficient, right_coefficient);
                for (coordinate, product_coefficient) in coordinates
                    .iter_mut()
                    .zip(self.algebra.basis_product(left_index, right_index))
                {
                    field.add_mul_assign(coordinate, &scalar, product_coefficient);
                }
            }
        }
        self.algebra.element(coordinates)
    }
}

/// Failure to construct or operate on a structure-constant algebra.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum StructureConstantsError {
    /// The requested multiplication tensor cube overflowed `usize`.
    #[error("the structure-constant tensor shape for dimension {dimension} overflows")]
    TensorShapeOverflow {
        /// Requested algebra dimension.
        dimension: usize,
    },
    /// The supplied flattened multiplication tensor has the wrong length.
    #[error(
        "the structure-constant tensor has {actual} entries, expected {expected} for dimension {dimension}"
    )]
    TensorLength {
        /// Supplied coefficient count.
        actual: usize,
        /// Required coefficient count.
        expected: usize,
        /// Requested algebra dimension.
        dimension: usize,
    },
    /// An element coordinate row has the wrong length.
    #[error("the algebra element has {actual} coordinates, expected {expected}")]
    ElementLength {
        /// Supplied coordinate count.
        actual: usize,
        /// Required coordinate count.
        expected: usize,
    },
    /// The block matrix used by QPA's identity method overflowed `usize`.
    #[error("the multiplicative-identity matrix shape for dimension {dimension} overflows")]
    IdentityMatrixShapeOverflow {
        /// Algebra dimension.
        dimension: usize,
    },
    /// GAP's `BlockMatrix` call has no block-shape information in dimension zero.
    #[error("need <rpb>, <cpb>, <zero> if <blocks> is empty")]
    EmptyIdentityBlockMatrix,
    /// The entered elements have different parent objects.
    #[error("the algebra elements belong to different structure-constant algebras")]
    DifferentParents,
    /// Exact matrix construction or elimination failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Ring,
        finite_field::{FiniteFieldCore, Z2},
        rational::Q,
    };
    use serde_json::{Value, json};

    use super::StructureConstantsAlgebra;

    #[test]
    fn finds_gap_left_identity_and_rejects_a_nonunital_algebra() {
        // e1^2=e1, e1e2=e2, e2e1=e2, e2^2=0.
        let q = |value: i64| Q.nth(value.into());
        let dual_numbers = StructureConstantsAlgebra::new(
            Q,
            2,
            [
                1, 0, // e1 e1
                0, 1, // e1 e2
                0, 1, // e2 e1
                0, 0, // e2 e2
            ]
            .into_iter()
            .map(q)
            .collect(),
        )
        .unwrap();
        assert_eq!(
            dual_numbers
                .find_multiplicative_identity()
                .unwrap()
                .unwrap()
                .coordinates(),
            &[q(1), q(0)]
        );

        let zero = StructureConstantsAlgebra::new(Z2, 1, vec![Z2.zero()]).unwrap();
        assert!(zero.find_multiplicative_identity().unwrap().is_none());

        let empty = StructureConstantsAlgebra::new(Q, 0, Vec::new()).unwrap();
        assert!(matches!(
            empty.find_multiplicative_identity(),
            Err(super::StructureConstantsError::EmptyIdentityBlockMatrix)
        ));
    }

    #[test]
    fn pinned_qpa_identity_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/algebra_structure_constants.schema-v1.json"
        ))
        .unwrap();

        let rational = StructureConstantsAlgebra::new(
            Q,
            2,
            [1, 0, 0, 1, 0, 0, 0, 0]
                .into_iter()
                .map(|value| Q.nth(value.into()))
                .collect(),
        )
        .unwrap();
        let rational_identity = rational.find_multiplicative_identity().unwrap().unwrap();
        assert_eq!(
            json!(
                rational_identity
                    .coordinates()
                    .iter()
                    .map(|coefficient| coefficient.to_string().parse::<i64>().unwrap())
                    .collect::<Vec<_>>()
            ),
            expected["rational"]["identity"]
        );

        let binary = StructureConstantsAlgebra::new(
            Z2,
            2,
            vec![
                Z2.one(),
                Z2.zero(),
                Z2.zero(),
                Z2.one(),
                Z2.zero(),
                Z2.zero(),
                Z2.zero(),
                Z2.zero(),
            ],
        )
        .unwrap();
        let binary_identity = binary.find_multiplicative_identity().unwrap().unwrap();
        assert_eq!(
            json!(
                binary_identity
                    .coordinates()
                    .iter()
                    .map(|coefficient| Z2
                        .to_integer(coefficient)
                        .to_string()
                        .parse::<i64>()
                        .unwrap())
                    .collect::<Vec<_>>()
            ),
            expected["f2"]["identity"]
        );
    }

    #[test]
    fn multiplies_in_standard_basis_coordinates() {
        let algebra = StructureConstantsAlgebra::new(
            Z2,
            2,
            vec![
                Z2.one(),
                Z2.zero(),
                Z2.zero(),
                Z2.one(),
                Z2.zero(),
                Z2.one(),
                Z2.zero(),
                Z2.zero(),
            ],
        )
        .unwrap();
        let basis = algebra.basis_elements();
        assert_eq!(basis[1].multiply(&basis[0]).unwrap(), basis[1]);
        assert_eq!(basis[1].multiply(&basis[1]).unwrap(), algebra.zero());
    }
}
