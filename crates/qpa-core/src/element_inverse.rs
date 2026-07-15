//! QPA inversion and signed powers for finite-dimensional quotient elements.

use numerica::domains::{Field, integer::Integer};
use thiserror::Error;

use crate::{
    DenseMatrix, MatrixError, QuotientError, QuotientNiceBasis, QuotientPathAlgebraElement,
};

impl<F: Field> QuotientPathAlgebraElement<F> {
    /// Translate QPA's `Inverse` method for a finite-dimensional quiver-algebra
    /// element.
    ///
    /// QPA takes coordinates in a basis whose first entries are replaced by
    /// the vertex idempotents. It returns `fail` when any corresponding vertex
    /// coefficient is zero. Otherwise it constructs the diagonal inverse and
    /// evaluates QPA's binomial polynomial after finding the nilpotence index
    /// of `1 - elm * diagonal_inverse`.
    ///
    /// `Ok(None)` is QPA's `fail` result.
    ///
    /// # Errors
    ///
    /// The quotient must carry a normal-form function and a complete finite
    /// canonical basis. Coordinate, exact-matrix, and quotient arithmetic
    /// failures are propagated.
    pub fn inverse(&self) -> Result<Option<Self>, ElementInverseError> {
        let algebra = self.algebra();
        if !algebra.has_normal_form_function() {
            return Err(ElementInverseError::MissingNormalFormFunction);
        }
        if !algebra.is_finite_dimensional()? {
            return Err(ElementInverseError::InfiniteDimensionalAlgebra);
        }
        let mut basis = algebra
            .canonical_basis()?
            .ok_or(ElementInverseError::CanonicalBasisUnavailable)?;
        let number_of_vertices = algebra
            .original_path_algebra()
            .quiver()
            .number_of_vertices();
        let generators = algebra.generators()?;
        if basis.len() < number_of_vertices || generators.len() < number_of_vertices {
            return Err(ElementInverseError::InvalidBasis);
        }
        basis[..number_of_vertices].clone_from_slice(&generators[..number_of_vertices]);

        let nice_basis = QuotientNiceBasis::new(algebra.clone(), &basis)?;
        let rows = basis
            .iter()
            .map(|element| {
                nice_basis
                    .nice_vector(element)?
                    .ok_or(ElementInverseError::InvalidBasis)
            })
            .collect::<Result<Vec<_>, ElementInverseError>>()?;
        let target = nice_basis
            .nice_vector(self)?
            .ok_or(ElementInverseError::ElementOutsideBasis)?;
        let matrix = DenseMatrix::from_rows(rows, algebra.original_path_algebra().field().clone())?;
        let coefficients = matrix
            .solution_row(&target)?
            .ok_or(ElementInverseError::ElementOutsideBasis)?;
        let vertex_coefficients = &coefficients[..number_of_vertices];
        let field = algebra.original_path_algebra().field();
        if vertex_coefficients
            .iter()
            .any(|coefficient| field.is_zero(coefficient))
        {
            return Ok(None);
        }

        let mut diagonal_inverse = algebra.zero();
        for (vertex, coefficient) in generators
            .iter()
            .take(number_of_vertices)
            .zip(vertex_coefficients)
        {
            diagonal_inverse =
                diagonal_inverse.add(&vertex.scale_left(&field.inv(coefficient))?)?;
        }
        let difference = algebra.one()?.sub(&self.multiply(&diagonal_inverse)?)?;
        let mut power = difference.clone();
        let mut nilpotence_index = 1_usize;
        while !power.is_zero() {
            nilpotence_index = nilpotence_index
                .checked_add(1)
                .ok_or(ElementInverseError::NilpotenceIndexOverflow)?;
            power = power.multiply(&difference)?;
        }

        let mut result = algebra.zero();
        let leading_binomial = field.nth(binomial(nilpotence_index, nilpotence_index - 1));
        result = result.add(&diagonal_inverse.scale_left(&leading_binomial)?)?;
        let normalized = self.multiply(&diagonal_inverse)?;
        for index in 2..=nilpotence_index {
            let mut scalar = field.nth(binomial(nilpotence_index, nilpotence_index - index));
            if (index + 1) % 2 == 1 {
                scalar = field.neg(&scalar);
            }
            let term = diagonal_inverse
                .multiply(&normalized.pow((index - 1) as u64)?)?
                .scale_left(&scalar)?;
            result = result.add(&term)?;
        }
        Ok(Some(result))
    }

    /// Translate QPA's negative- and zero-exponent methods together with its
    /// existing nonnegative quotient power.
    ///
    /// `Ok(None)` is QPA's `fail` for a negative power of a nonunit. Exponent
    /// zero returns the algebra identity even when the base is not invertible.
    ///
    /// # Errors
    ///
    /// Propagates source-translated inversion and quotient arithmetic failures.
    pub fn pow_integer(&self, exponent: i64) -> Result<Option<Self>, ElementInverseError> {
        if exponent >= 0 {
            return Ok(Some(self.pow(exponent.unsigned_abs())?));
        }
        let Some(inverse) = self.inverse()? else {
            return Ok(None);
        };
        Ok(Some(inverse.pow(exponent.unsigned_abs())?))
    }
}

fn binomial(n: usize, k: usize) -> Integer {
    let k = k.min(n - k);
    let mut result = Integer::from(1);
    for index in 1..=k {
        result *= Integer::from(n - k + index);
        result /= Integer::from(index);
    }
    result
}

/// Failure in QPA's quotient-element inversion method.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ElementInverseError {
    /// QPA quotient multiplication cannot reduce representatives in this family.
    #[error("the quotient family has no normal-form function")]
    MissingNormalFormFunction,
    /// QPA explicitly rejects inverse computation in an infinite algebra.
    #[error("the element is not in a finite-dimensional algebra")]
    InfiniteDimensionalAlgebra,
    /// QPA's `Basis(A)` has no available vectors in the current cache state.
    #[error("the quotient canonical basis has no available vectors")]
    CanonicalBasisUnavailable,
    /// The supplied canonical vectors did not form QPA's expected basis.
    #[error("the quotient canonical vectors do not form the expected basis")]
    InvalidBasis,
    /// Coordinate extraction unexpectedly failed for an algebra element.
    #[error("the element does not lie in the selected algebra basis")]
    ElementOutsideBasis,
    /// The source loop's nilpotence counter exceeded Rust's addressable size.
    #[error("the nilpotence index exceeds Rust's addressable size")]
    NilpotenceIndexOverflow,
    /// Quotient basis construction or arithmetic failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// GAP-compatible coordinate matrix construction or solving failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
}

#[cfg(test)]
mod tests {
    use super::binomial;
    use numerica::domains::integer::Integer;

    #[test]
    fn arbitrary_precision_binomial_matches_qpa_coefficients() {
        assert_eq!(binomial(1, 0), Integer::from(1));
        assert_eq!(binomial(5, 2), Integer::from(10));
        assert_eq!(
            binomial(100, 50).to_string(),
            "100891344545564193334812497256"
        );
    }
}
