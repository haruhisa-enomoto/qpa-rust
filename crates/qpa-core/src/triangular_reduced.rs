//! QPA's finite-dimensional triangular-reduced property.

use numerica::domains::Field;
use thiserror::Error;

use crate::{PathAlgebra, PathAlgebraError, QuotientError, QuotientPathAlgebra};

impl<F: Field> PathAlgebra<F> {
    /// Translate QPA's `IsTriangularReduced` method for a path algebra.
    ///
    /// QPA enumerates every nonempty proper subset of vertices by cardinality,
    /// forms its idempotent sum `e`, and tests every canonical basis vector in
    /// `e A (1-e)`. The first all-zero product list returns `false`.
    ///
    /// # Errors
    ///
    /// QPA raises an error for an infinite-dimensional algebra. Sparse path
    /// arithmetic failures are propagated.
    pub fn is_triangular_reduced(&self) -> Result<bool, TriangularReducedError> {
        let basis = self
            .canonical_basis()
            .ok_or(TriangularReducedError::InfiniteDimensionalAlgebra)?;
        let one = self.one();
        let generators = self
            .quiver()
            .vertices()
            .map(|vertex| Ok(one.multiply(&self.element(vertex.path())?)?))
            .collect::<Result<Vec<_>, TriangularReducedError>>()?;
        for size in 1..self.quiver().number_of_vertices() {
            for combination in combinations(self.quiver().number_of_vertices(), size) {
                let mut idempotent = self.zero();
                for position in combination {
                    idempotent = idempotent.add(&generators[position])?;
                }
                let complement = one.sub(&idempotent)?;
                let products = basis
                    .iter()
                    .map(|vector| Ok(idempotent.multiply(vector)?.multiply(&complement)?))
                    .collect::<Result<Vec<_>, TriangularReducedError>>()?;
                if products.iter().all(crate::PathAlgebraElement::is_zero) {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }
}

impl<F: Field> QuotientPathAlgebra<F> {
    /// Translate QPA's `IsTriangularReduced` method for a quotient path
    /// algebra using its actual quotient basis and products.
    ///
    /// # Errors
    ///
    /// QPA raises an error for an infinite-dimensional algebra. The canonical
    /// basis must be available in the current cache state; quotient arithmetic
    /// failures are propagated.
    pub fn is_triangular_reduced(&self) -> Result<bool, TriangularReducedError> {
        if !self.is_finite_dimensional()? {
            return Err(TriangularReducedError::InfiniteDimensionalAlgebra);
        }
        let basis = self
            .canonical_basis()?
            .ok_or(TriangularReducedError::CanonicalBasisUnavailable)?;
        let one = self.one()?;
        let generators = self
            .original_path_algebra()
            .quiver()
            .vertices()
            .map(|vertex| Ok(one.multiply(&self.path(vertex.path())?)?))
            .collect::<Result<Vec<_>, TriangularReducedError>>()?;
        for size in 1..generators.len() {
            for combination in combinations(generators.len(), size) {
                let mut idempotent = self.zero();
                for position in combination {
                    idempotent = idempotent.add(&generators[position])?;
                }
                let complement = one.sub(&idempotent)?;
                let products = basis
                    .iter()
                    .map(|vector| Ok(idempotent.multiply(vector)?.multiply(&complement)?))
                    .collect::<Result<Vec<_>, TriangularReducedError>>()?;
                if products
                    .iter()
                    .all(crate::QuotientPathAlgebraElement::is_zero)
                {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }
}

fn combinations(length: usize, size: usize) -> Vec<Vec<usize>> {
    fn append(
        length: usize,
        size: usize,
        next: usize,
        current: &mut Vec<usize>,
        result: &mut Vec<Vec<usize>>,
    ) {
        if current.len() == size {
            result.push(current.clone());
            return;
        }
        let needed = size - current.len();
        for position in next..=length - needed {
            current.push(position);
            append(length, size, position + 1, current, result);
            current.pop();
        }
    }

    let mut result = Vec::new();
    append(length, size, 0, &mut Vec::with_capacity(size), &mut result);
    result
}

/// Failure in QPA's triangular-reduced property method.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum TriangularReducedError {
    /// QPA explicitly rejects infinite-dimensional input.
    #[error("the entered algebra is not finite dimensional")]
    InfiniteDimensionalAlgebra,
    /// QPA's cached canonical basis object has no vectors.
    #[error("the quotient canonical basis has no available vectors")]
    CanonicalBasisUnavailable,
    /// Parent path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Quotient basis construction or arithmetic failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
}

#[cfg(test)]
mod tests {
    use super::combinations;

    #[test]
    fn combinations_follow_gap_lexicographic_order() {
        assert_eq!(
            combinations(4, 2),
            vec![
                vec![0, 1],
                vec![0, 2],
                vec![0, 3],
                vec![1, 2],
                vec![1, 3],
                vec![2, 3]
            ]
        );
    }
}
