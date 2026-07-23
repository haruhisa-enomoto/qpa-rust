//! QPA's finite-dimensional idempotent algorithms.
//!
//! This module translates `lib/idempotent.gi`.  Its polynomial factorization
//! follows the finite-field `Factors` dependency in GAP's `lib/polyfinf.gi`:
//! square-free reduction, distinct-degree splitting, and randomized
//! equal-degree splitting.

use std::cmp::Ordering;

use numerica::domains::{Field, integer::Integer};
use rand::{Rng, SeedableRng, rngs::StdRng};
use thiserror::Error;

use crate::{
    CoefficientDomainError, DenseMatrix, MatrixError, QpaCoefficientField, QpaFieldCharacteristic,
    StructureConstantsAlgebra, StructureConstantsAlgebraElement, StructureConstantsError,
    StructureConstantsSubalgebra,
    structure_constants::{
        add_rows, evaluate_algebra_polynomial, is_commutative, minimum_polynomial, span_rows,
        subtract_rows,
    },
};

/// Reproducible seed used by the no-argument HOPF and primitive-idempotent
/// entry points.
///
/// QPA consumes GAP's mutable global random source. Rust additionally exposes
/// `_with_rng` variants so callers can reproduce or vary the probabilistic
/// source algorithm explicitly.
pub const HOPF_DEFAULT_SEED: u64 = 0x5150_4148_4f50_4601;

type PolynomialDivision<E> = (Vec<E>, Vec<E>);

impl<F: Field> StructureConstantsAlgebra<F> {
    /// Translate QPA's `HOPF_MinPoly(A, a)` power-span algorithm.
    ///
    /// Coefficients are returned from constant term to leading coefficient,
    /// exactly as passed to GAP's `LaurentPolynomialByCoefficients`.
    ///
    /// # Errors
    ///
    /// Rejects an element from another algebra and propagates the source
    /// identity and exact row-basis computations.
    pub fn hopf_minimum_polynomial(
        &self,
        element: &StructureConstantsAlgebraElement<F>,
    ) -> Result<Vec<F::Element>, HopfError> {
        if element.algebra() != self {
            return Err(StructureConstantsError::DifferentParents.into());
        }
        let identity = self
            .find_multiplicative_identity()?
            .ok_or(StructureConstantsError::MissingMultiplicativeIdentity)?;
        Ok(minimum_polynomial(
            self,
            identity.coordinates(),
            element.coordinates(),
        )?)
    }

    /// Translate QPA's `HOPF_LeftIdentity(A)` row-system solve.
    ///
    /// `None` is GAP's `fail` when the entered finite-dimensional algebra has
    /// no left identity.
    ///
    /// # Errors
    ///
    /// Propagates the exact shaped-matrix construction and solution failures.
    pub fn hopf_left_identity(
        &self,
    ) -> Result<Option<StructureConstantsAlgebraElement<F>>, HopfError> {
        Ok(self.find_multiplicative_identity()?)
    }
}

impl<F: QpaCoefficientField> StructureConstantsAlgebra<F> {
    /// Translate QPA's randomized `HOPF_ZeroDivisor(A)` using an explicit RNG.
    ///
    /// `None` is QPA's empty list for a commutative algebra. Otherwise the two
    /// returned nonzero polynomial values multiply to zero.
    ///
    /// # Errors
    ///
    /// Requires the finite prime-field tier for a noncommutative algebra and
    /// propagates exact minimum-polynomial and finite-field factorization
    /// failures.
    pub fn hopf_zero_divisor_with_rng<R: Rng + ?Sized>(
        &self,
        rng: &mut R,
    ) -> Result<Option<[StructureConstantsAlgebraElement<F>; 2]>, HopfError> {
        if is_commutative(self) {
            return Ok(None);
        }
        let elements = self
            .field()
            .qpa_elements()?
            .ok_or(HopfError::FiniteAlgebraRequired)?;
        let identity = self
            .find_multiplicative_identity()?
            .ok_or(StructureConstantsError::MissingMultiplicativeIdentity)?;

        loop {
            let coordinates = (0..self.dimension())
                .map(|_| elements[rng.random_range(0..elements.len())].clone())
                .collect::<Vec<_>>();
            let polynomial = minimum_polynomial(self, identity.coordinates(), &coordinates)?;
            let factors =
                factor_finite_field_polynomial(self.field(), &elements, &polynomial, rng)?;
            if factors.len() <= 1 {
                continue;
            }
            let first = evaluate_algebra_polynomial(
                self,
                &factors[0],
                &coordinates,
                identity.coordinates(),
            );
            let remaining = factors[1..]
                .iter()
                .fold(vec![self.field().one()], |product, factor| {
                    polynomial_product(self.field(), &product, factor)
                });
            let second =
                evaluate_algebra_polynomial(self, &remaining, &coordinates, identity.coordinates());
            return Ok(Some([self.element(first)?, self.element(second)?]));
        }
    }

    /// Run QPA's HOPF zero-divisor algorithm with a reproducible default RNG.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::hopf_zero_divisor_with_rng`] failures.
    pub fn hopf_zero_divisor(
        &self,
    ) -> Result<Option<[StructureConstantsAlgebraElement<F>; 2]>, HopfError> {
        let mut rng = StdRng::seed_from_u64(HOPF_DEFAULT_SEED);
        self.hopf_zero_divisor_with_rng(&mut rng)
    }

    /// Translate QPA's recursive `HOPF_SingularIdempotent(A)` corner
    /// construction using an explicit RNG.
    ///
    /// # Errors
    ///
    /// Propagates zero-divisor, right-ideal, left-identity, corner-algebra, and
    /// exact coordinate-lifting failures.
    pub fn hopf_singular_idempotent_with_rng<R: Rng + ?Sized>(
        &self,
        rng: &mut R,
    ) -> Result<StructureConstantsAlgebraElement<F>, HopfError> {
        let mut current = self.clone();
        let mut current_to_original = DenseMatrix::identity(self.dimension(), self.field().clone());

        while !is_commutative(&current) {
            let Some([zero_divisor, _]) = current.hopf_zero_divisor_with_rng(rng)? else {
                return Err(HopfError::NoncommutativeZeroDivisorMissing);
            };
            let ideal = right_principal_subalgebra(&current, zero_divisor.coordinates())?;
            let ideal_identity = ideal
                .algebra()
                .hopf_left_identity()?
                .ok_or(HopfError::MissingLeftIdentity)?;
            let corner_identity = ideal.lift_coordinates(ideal_identity.coordinates())?;
            let corner = sandwich_subalgebra(&current, &corner_identity)?;
            current_to_original = corner.basis().multiply(&current_to_original)?;
            current = corner.algebra().clone();
        }

        let identity = current
            .find_multiplicative_identity()?
            .ok_or(StructureConstantsError::MissingMultiplicativeIdentity)?;
        self.element(multiply_row_matrix(
            identity.coordinates(),
            &current_to_original,
        ))
        .map_err(Into::into)
    }

    /// Run QPA's singular-idempotent algorithm with a reproducible default
    /// RNG.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::hopf_singular_idempotent_with_rng`] failures.
    pub fn hopf_singular_idempotent(
        &self,
    ) -> Result<StructureConstantsAlgebraElement<F>, HopfError> {
        let mut rng = StdRng::seed_from_u64(HOPF_DEFAULT_SEED);
        self.hopf_singular_idempotent_with_rng(&mut rng)
    }

    /// Translate QPA's cached `PrimitiveIdempotents(A)` attribute using an
    /// explicit RNG for the first evaluation.
    ///
    /// The method retains QPA's check order on the supported finite-field
    /// domain, recursive singular idempotent, right-ideal dimension, ordered
    /// standard-basis search, and left-identity construction.
    ///
    /// # Errors
    ///
    /// Rejects non-semisimple, nonsimple, or non-finite algebras and
    /// propagates all exact HOPF and structure-constant failures.
    pub fn primitive_idempotents_with_rng<R: Rng + ?Sized>(
        &self,
        rng: &mut R,
    ) -> Result<Vec<StructureConstantsAlgebraElement<F>>, HopfError> {
        if let Some(cached) = self.cached_primitive_idempotents() {
            return cached
                .iter()
                .cloned()
                .map(|coordinates| self.element(coordinates).map_err(Into::into))
                .collect();
        }
        if !self.is_semisimple_algebra()? {
            return Err(HopfError::NotSemisimple);
        }
        if self.field().qpa_elements()?.is_none() {
            return Err(HopfError::FiniteAlgebraRequired);
        }
        if self.central_idempotents_split_prime()?.len() > 1 {
            return Err(HopfError::NotSimple);
        }

        let singular = self.hopf_singular_idempotent_with_rng(rng)?;
        let first_ideal = right_principal_subalgebra(self, singular.coordinates())?;
        let ideal_dimension = first_ideal.algebra().dimension();
        if ideal_dimension == 0 || !self.dimension().is_multiple_of(ideal_dimension) {
            return Err(HopfError::InvalidSimpleAlgebraDimensions {
                algebra: self.dimension(),
                minimal_right_ideal: ideal_dimension,
            });
        }
        let count = self.dimension() / ideal_dimension;
        let identity = self
            .find_multiplicative_identity()?
            .ok_or(StructureConstantsError::MissingMultiplicativeIdentity)?;
        let basis = self.basis_elements();
        let mut idempotents = Vec::with_capacity(count);
        idempotents.push(singular.coordinates().to_vec());
        let mut sum = vec![self.field().zero(); self.dimension()];

        for position in 1..count {
            sum = add_rows(self.field(), &sum, &idempotents[position - 1]);
            let complement = subtract_rows(self.field(), identity.coordinates(), &sum);
            let candidate = basis
                .iter()
                .map(|basis_element| {
                    let left = self.coordinate_product(&complement, basis_element.coordinates());
                    self.coordinate_product(&left, singular.coordinates())
                })
                .find(|candidate| !is_zero_polynomial(self.field(), candidate))
                .ok_or(HopfError::PrimitiveBasisSearchExhausted {
                    position: position + 1,
                })?;
            let ideal = right_principal_subalgebra(self, &candidate)?;
            let ideal_identity = ideal
                .algebra()
                .hopf_left_identity()?
                .ok_or(HopfError::MissingLeftIdentity)?;
            let lifted_identity = ideal.lift_coordinates(ideal_identity.coordinates())?;
            idempotents.push(self.coordinate_product(&lifted_identity, &complement));
        }

        let cached = self.cache_primitive_idempotents(idempotents);
        cached
            .iter()
            .cloned()
            .map(|coordinates| self.element(coordinates).map_err(Into::into))
            .collect()
    }

    /// Run QPA's cached primitive-idempotent algorithm with a reproducible
    /// default RNG.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::primitive_idempotents_with_rng`] failures.
    pub fn primitive_idempotents(
        &self,
    ) -> Result<Vec<StructureConstantsAlgebraElement<F>>, HopfError> {
        let mut rng = StdRng::seed_from_u64(HOPF_DEFAULT_SEED);
        self.primitive_idempotents_with_rng(&mut rng)
    }
}

fn right_principal_subalgebra<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    generator: &[F::Element],
) -> Result<StructureConstantsSubalgebra<F>, HopfError> {
    let rows = algebra
        .basis_elements()
        .iter()
        .map(|basis| algebra.coordinate_product(generator, basis.coordinates()))
        .collect::<Vec<_>>();
    let basis = span_rows(algebra.field(), algebra.dimension(), rows)?;
    Ok(StructureConstantsSubalgebra::from_basis(
        algebra.clone(),
        basis,
    )?)
}

fn sandwich_subalgebra<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    idempotent: &[F::Element],
) -> Result<StructureConstantsSubalgebra<F>, HopfError> {
    let rows = algebra
        .basis_elements()
        .iter()
        .map(|basis| {
            let left = algebra.coordinate_product(idempotent, basis.coordinates());
            algebra.coordinate_product(&left, idempotent)
        })
        .collect::<Vec<_>>();
    let basis = span_rows(algebra.field(), algebra.dimension(), rows)?;
    Ok(StructureConstantsSubalgebra::from_basis(
        algebra.clone(),
        basis,
    )?)
}

fn multiply_row_matrix<F: Field>(row: &[F::Element], matrix: &DenseMatrix<F>) -> Vec<F::Element> {
    debug_assert_eq!(row.len(), matrix.nrows());
    let mut result = vec![matrix.field().zero(); matrix.ncols()];
    for (inner, coefficient) in row.iter().enumerate() {
        if matrix.field().is_zero(coefficient) {
            continue;
        }
        for column in 0..matrix.ncols() {
            matrix.field().add_mul_assign(
                &mut result[column],
                coefficient,
                &matrix[(inner, column)],
            );
        }
    }
    result
}

pub(crate) fn factor_finite_field_polynomial<F: QpaCoefficientField, R: Rng + ?Sized>(
    field: &F,
    elements: &[F::Element],
    polynomial: &[F::Element],
    rng: &mut R,
) -> Result<Vec<Vec<F::Element>>, HopfError> {
    let polynomial = monic_polynomial(field, polynomial)?;
    let Some(degree) = polynomial_degree(field, &polynomial) else {
        return Err(HopfError::ZeroPolynomialFactorization);
    };
    if degree < 2 {
        return Ok(vec![polynomial]);
    }

    let derivative = polynomial_derivative(field, &polynomial);
    let mut factors = if is_zero_polynomial(field, &derivative) {
        let prime = prime_characteristic(field)?;
        let stride = usize::try_from(prime).map_err(|_| HopfError::PolynomialDegreeOverflow)?;
        let root = polynomial
            .iter()
            .step_by(stride)
            .map(|coefficient| {
                field
                    .qpa_frobenius(coefficient, -1)
                    .expect("finite fields have inverse Frobenius")
            })
            .collect::<Vec<_>>();
        let root_factors = factor_finite_field_polynomial(field, elements, &root, rng)?;
        let repetitions = stride;
        (0..repetitions)
            .flat_map(|_| root_factors.iter().cloned())
            .collect::<Vec<_>>()
    } else {
        let mut repeated = polynomial_gcd(field, &polynomial, &derivative)?;
        let square_free = exact_polynomial_quotient(field, &polynomial, &repeated)?;
        let mut factors = factor_square_free(field, elements, &square_free, rng)?;
        for factor in factors.clone() {
            loop {
                let (quotient, remainder) = polynomial_div_rem(field, &repeated, &factor)?;
                if !is_zero_polynomial(field, &remainder) {
                    break;
                }
                factors.push(factor.clone());
                repeated = quotient;
            }
        }
        if polynomial_degree(field, &repeated).is_some_and(|degree| degree > 0) {
            factors.extend(factor_finite_field_polynomial(
                field, elements, &repeated, rng,
            )?);
        }
        factors
    };

    factors.sort_by(|left, right| compare_polynomials(field, elements, left, right));
    Ok(factors)
}

fn factor_square_free<F: QpaCoefficientField, R: Rng + ?Sized>(
    field: &F,
    elements: &[F::Element],
    polynomial: &[F::Element],
    rng: &mut R,
) -> Result<Vec<Vec<F::Element>>, HopfError> {
    let mut remaining = polynomial.to_vec();
    let mut factors = Vec::new();
    let mut power = vec![field.zero(), field.one()];
    let x = power.clone();
    let field_order = field.qpa_order().ok_or(HopfError::FiniteAlgebraRequired)?;
    let mut degree = 0usize;

    while polynomial_degree(field, &remaining)
        .is_some_and(|remaining_degree| 2 * (degree + 1) <= remaining_degree)
    {
        power = polynomial_power_mod_integer(field, &power, field_order.clone(), &remaining)?;
        degree += 1;
        let cyclotomic = polynomial_subtract(field, &power, &x);
        let common = polynomial_gcd(field, &remaining, &cyclotomic)?;
        if polynomial_degree(field, &common).is_some_and(|common_degree| common_degree >= degree) {
            factors.extend(factor_common_degree(field, elements, &common, degree, rng)?);
            remaining = exact_polynomial_quotient(field, &remaining, &common)?;
            if polynomial_degree(field, &remaining).is_some_and(|value| value > 0) {
                power = polynomial_remainder(field, &power, &remaining)?;
            }
        }
    }
    if polynomial_degree(field, &remaining).is_some_and(|degree| degree > 0) {
        factors.push(monic_polynomial(field, &remaining)?);
    }
    Ok(factors)
}

fn factor_common_degree<F: QpaCoefficientField, R: Rng + ?Sized>(
    field: &F,
    elements: &[F::Element],
    polynomial: &[F::Element],
    factor_degree: usize,
    rng: &mut R,
) -> Result<Vec<Vec<F::Element>>, HopfError> {
    let degree = polynomial_degree(field, polynomial).unwrap_or(0);
    if degree < factor_degree {
        return Ok(Vec::new());
    }
    if degree == factor_degree {
        return Ok(vec![monic_polynomial(field, polynomial)?]);
    }

    loop {
        let mut random = (0..factor_degree.saturating_mul(2))
            .map(|_| elements[rng.random_range(0..elements.len())].clone())
            .collect::<Vec<_>>();
        trim_polynomial(field, &mut random);
        if random.is_empty() {
            continue;
        }

        let prime = prime_characteristic(field)?;
        let field_order = field.qpa_order().ok_or(HopfError::FiniteAlgebraRequired)?;
        let separator = if prime == 2 {
            let mut term = polynomial_remainder(field, &random, polynomial)?;
            let mut trace = term.clone();
            for _ in 1..factor_degree {
                term = polynomial_power_mod_integer(field, &term, field_order.clone(), polynomial)?;
                trace = polynomial_add(field, &trace, &term);
                trace = polynomial_remainder(field, &trace, polynomial)?;
            }
            trace
        } else {
            let factor_degree =
                u64::try_from(factor_degree).map_err(|_| HopfError::PolynomialDegreeOverflow)?;
            let mut exponent = field_order.pow(factor_degree);
            exponent -= Integer::from(1);
            exponent /= Integer::from(2);
            let mut result = polynomial_power_mod_integer(field, &random, exponent, polynomial)?;
            if result.is_empty() {
                result.push(field.neg(&field.one()));
            } else {
                field.sub_assign(&mut result[0], &field.one());
            }
            trim_polynomial(field, &mut result);
            result
        };
        let common = polynomial_gcd(field, polynomial, &separator)?;
        let common_degree = polynomial_degree(field, &common).unwrap_or(0);
        if common_degree == 0 || common_degree == degree {
            continue;
        }
        let quotient = exact_polynomial_quotient(field, polynomial, &common)?;
        let mut factors = factor_common_degree(field, elements, &quotient, factor_degree, rng)?;
        factors.extend(factor_common_degree(
            field,
            elements,
            &common,
            factor_degree,
            rng,
        )?);
        return Ok(factors);
    }
}

fn prime_characteristic<F: QpaCoefficientField>(field: &F) -> Result<u64, HopfError> {
    match field.qpa_characteristic() {
        QpaFieldCharacteristic::Prime(prime) => Ok(prime),
        QpaFieldCharacteristic::Zero => Err(HopfError::FiniteAlgebraRequired),
    }
}

fn polynomial_degree<F: Field>(field: &F, polynomial: &[F::Element]) -> Option<usize> {
    polynomial
        .iter()
        .rposition(|coefficient| !field.is_zero(coefficient))
}

fn trim_polynomial<F: Field>(field: &F, polynomial: &mut Vec<F::Element>) {
    while polynomial.last().is_some_and(|entry| field.is_zero(entry)) {
        polynomial.pop();
    }
}

fn is_zero_polynomial<F: Field>(field: &F, polynomial: &[F::Element]) -> bool {
    polynomial.iter().all(|entry| field.is_zero(entry))
}

fn monic_polynomial<F: Field>(
    field: &F,
    polynomial: &[F::Element],
) -> Result<Vec<F::Element>, HopfError> {
    let Some(degree) = polynomial_degree(field, polynomial) else {
        return Err(HopfError::ZeroPolynomialFactorization);
    };
    let inverse = field.inv(&polynomial[degree]);
    Ok(polynomial[..=degree]
        .iter()
        .map(|coefficient| field.mul(coefficient, &inverse))
        .collect())
}

fn polynomial_derivative<F: Field>(field: &F, polynomial: &[F::Element]) -> Vec<F::Element> {
    let mut derivative = polynomial
        .iter()
        .enumerate()
        .skip(1)
        .map(|(degree, coefficient)| field.mul(coefficient, &field.nth(Integer::from(degree))))
        .collect::<Vec<_>>();
    trim_polynomial(field, &mut derivative);
    derivative
}

fn polynomial_add<F: Field>(
    field: &F,
    left: &[F::Element],
    right: &[F::Element],
) -> Vec<F::Element> {
    let mut result = vec![field.zero(); left.len().max(right.len())];
    for (position, entry) in result.iter_mut().enumerate() {
        if let Some(value) = left.get(position) {
            field.add_assign(entry, value);
        }
        if let Some(value) = right.get(position) {
            field.add_assign(entry, value);
        }
    }
    trim_polynomial(field, &mut result);
    result
}

fn polynomial_subtract<F: Field>(
    field: &F,
    left: &[F::Element],
    right: &[F::Element],
) -> Vec<F::Element> {
    let mut result = vec![field.zero(); left.len().max(right.len())];
    for (position, entry) in result.iter_mut().enumerate() {
        if let Some(value) = left.get(position) {
            field.add_assign(entry, value);
        }
        if let Some(value) = right.get(position) {
            field.sub_assign(entry, value);
        }
    }
    trim_polynomial(field, &mut result);
    result
}

pub(crate) fn polynomial_product<F: Field>(
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
            field.add_mul_assign(
                &mut result[left_degree + right_degree],
                left_coefficient,
                right_coefficient,
            );
        }
    }
    trim_polynomial(field, &mut result);
    result
}

fn polynomial_div_rem<F: Field>(
    field: &F,
    dividend: &[F::Element],
    divisor: &[F::Element],
) -> Result<PolynomialDivision<F::Element>, HopfError> {
    let Some(divisor_degree) = polynomial_degree(field, divisor) else {
        return Err(HopfError::PolynomialDivisionByZero);
    };
    let mut remainder = dividend.to_vec();
    trim_polynomial(field, &mut remainder);
    let Some(dividend_degree) = polynomial_degree(field, &remainder) else {
        return Ok((Vec::new(), Vec::new()));
    };
    if dividend_degree < divisor_degree {
        return Ok((Vec::new(), remainder));
    }
    let mut quotient = vec![field.zero(); dividend_degree - divisor_degree + 1];
    let leading_inverse = field.inv(&divisor[divisor_degree]);
    while let Some(remainder_degree) = polynomial_degree(field, &remainder) {
        if remainder_degree < divisor_degree {
            break;
        }
        let shift = remainder_degree - divisor_degree;
        let coefficient = field.mul(&remainder[remainder_degree], &leading_inverse);
        quotient[shift] = coefficient.clone();
        for (degree, divisor_coefficient) in divisor[..=divisor_degree].iter().enumerate() {
            field.sub_mul_assign(
                &mut remainder[shift + degree],
                &coefficient,
                divisor_coefficient,
            );
        }
        trim_polynomial(field, &mut remainder);
    }
    trim_polynomial(field, &mut quotient);
    Ok((quotient, remainder))
}

fn exact_polynomial_quotient<F: Field>(
    field: &F,
    dividend: &[F::Element],
    divisor: &[F::Element],
) -> Result<Vec<F::Element>, HopfError> {
    let (quotient, remainder) = polynomial_div_rem(field, dividend, divisor)?;
    if !is_zero_polynomial(field, &remainder) {
        return Err(HopfError::InexactPolynomialDivision);
    }
    Ok(quotient)
}

fn polynomial_remainder<F: Field>(
    field: &F,
    polynomial: &[F::Element],
    modulus: &[F::Element],
) -> Result<Vec<F::Element>, HopfError> {
    polynomial_div_rem(field, polynomial, modulus).map(|(_, remainder)| remainder)
}

fn polynomial_gcd<F: Field>(
    field: &F,
    left: &[F::Element],
    right: &[F::Element],
) -> Result<Vec<F::Element>, HopfError> {
    let mut left = left.to_vec();
    let mut right = right.to_vec();
    trim_polynomial(field, &mut left);
    trim_polynomial(field, &mut right);
    while !right.is_empty() {
        let remainder = polynomial_remainder(field, &left, &right)?;
        left = right;
        right = remainder;
    }
    monic_polynomial(field, &left)
}

fn polynomial_power_mod_integer<F: Field>(
    field: &F,
    base: &[F::Element],
    mut exponent: Integer,
    modulus: &[F::Element],
) -> Result<Vec<F::Element>, HopfError> {
    let zero = Integer::from(0);
    let one = Integer::from(1);
    let mut result = vec![field.one()];
    let mut base = polynomial_remainder(field, base, modulus)?;
    while exponent != zero {
        if (&exponent & &one) == one {
            result = polynomial_product(field, &result, &base);
            result = polynomial_remainder(field, &result, modulus)?;
        }
        exponent >>= 1_u32;
        if exponent != zero {
            base = polynomial_product(field, &base, &base);
            base = polynomial_remainder(field, &base, modulus)?;
        }
    }
    Ok(result)
}

fn compare_polynomials<F: QpaCoefficientField>(
    field: &F,
    elements: &[F::Element],
    left: &[F::Element],
    right: &[F::Element],
) -> Ordering {
    let left_degree = polynomial_degree(field, left);
    let right_degree = polynomial_degree(field, right);
    match left_degree.cmp(&right_degree) {
        Ordering::Equal => {}
        ordering => return ordering,
    }
    for (left, right) in left.iter().rev().zip(right.iter().rev()) {
        let left_position = elements.iter().position(|element| element == left);
        let right_position = elements.iter().position(|element| element == right);
        match left_position.cmp(&right_position) {
            Ordering::Equal => {}
            ordering => return ordering,
        }
    }
    left.len().cmp(&right.len())
}

/// Failure in QPA's HOPF and primitive-idempotent algorithms.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum HopfError {
    /// `PrimitiveIdempotents` requires a semisimple algebra.
    #[error("the entered algebra is not a semisimple algebra")]
    NotSemisimple,
    /// `PrimitiveIdempotents` requires a single simple block.
    #[error("the entered algebra is not a simple algebra")]
    NotSimple,
    /// A noncommutative HOPF computation requires a finite field.
    #[error("the entered algebra is not over a supported finite field")]
    FiniteAlgebraRequired,
    /// QPA assumes that the right ideal constructed from its zero divisor has
    /// a left identity.
    #[error("the HOPF right ideal has no left identity")]
    MissingLeftIdentity,
    /// An internal commutativity branch contradicted QPA's zero-divisor
    /// result.
    #[error("a noncommutative algebra returned QPA's empty zero-divisor result")]
    NoncommutativeZeroDivisorMissing,
    /// The simple-algebra matrix and minimal-right-ideal dimensions did not
    /// have QPA's required integral ratio.
    #[error(
        "algebra dimension {algebra} is not divisible by minimal right-ideal dimension {minimal_right_ideal}"
    )]
    InvalidSimpleAlgebraDimensions {
        /// Dimension of the entered algebra.
        algebra: usize,
        /// Dimension of the computed minimal right ideal.
        minimal_right_ideal: usize,
    },
    /// QPA's ordered basis scan ran out before finding the next primitive
    /// component.
    #[error("no nonzero basis candidate exists for primitive idempotent {position}")]
    PrimitiveBasisSearchExhausted {
        /// One-based primitive-idempotent position.
        position: usize,
    },
    /// GAP's finite-field factorization was asked to factor zero.
    #[error("cannot factor the zero polynomial")]
    ZeroPolynomialFactorization,
    /// An internal polynomial division had zero divisor.
    #[error("cannot divide a polynomial by zero")]
    PolynomialDivisionByZero,
    /// A division that GAP's factor algorithm requires to be exact had a
    /// nonzero remainder.
    #[error("finite-field polynomial division was unexpectedly inexact")]
    InexactPolynomialDivision,
    /// A finite-field characteristic or degree was not addressable by Rust.
    #[error("finite-field polynomial degree is not addressable")]
    PolynomialDegreeOverflow,
    /// Exact structure-constant algebra construction failed.
    #[error(transparent)]
    Structure(#[from] StructureConstantsError),
    /// Exact shaped-matrix arithmetic failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// The Numerica field adapter could not enumerate the finite field.
    #[error(transparent)]
    Coefficient(#[from] CoefficientDomainError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring, RingOps,
        finite_field::{FiniteFieldCore, Z2, Zp},
        rational::Q,
    };
    use rand::{SeedableRng, rngs::StdRng};
    use serde_json::{Value, json};

    use crate::{QpaCoefficientField, StructureConstantsAlgebra};

    use super::{
        HOPF_DEFAULT_SEED, HopfError, factor_finite_field_polynomial, polynomial_product,
        right_principal_subalgebra,
    };

    fn full_matrix_algebra<F: Field>(field: F, size: usize) -> StructureConstantsAlgebra<F> {
        let dimension = size * size;
        let mut products = vec![field.zero(); dimension * dimension * dimension];
        for left_row in 0..size {
            for inner in 0..size {
                let left = left_row * size + inner;
                for right_column in 0..size {
                    let right = inner * size + right_column;
                    let product = left_row * size + right_column;
                    products[(left * dimension + right) * dimension + product] = field.one();
                }
            }
        }
        StructureConstantsAlgebra::new(field, dimension, products).unwrap()
    }

    fn product_of_fields<F: Field>(field: F, count: usize) -> StructureConstantsAlgebra<F> {
        let mut products = vec![field.zero(); count * count * count];
        for position in 0..count {
            products[(position * count + position) * count + position] = field.one();
        }
        StructureConstantsAlgebra::new(field, count, products).unwrap()
    }

    fn dual_numbers<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        let zero = field.zero();
        let one = field.one();
        StructureConstantsAlgebra::new(
            field,
            2,
            vec![
                one.clone(),
                zero.clone(),
                zero.clone(),
                one.clone(),
                zero.clone(),
                one,
                zero.clone(),
                zero,
            ],
        )
        .unwrap()
    }

    fn zero_algebra<F: Field>(field: &F) -> StructureConstantsAlgebra<F> {
        StructureConstantsAlgebra::new(field.clone(), 1, vec![field.zero()]).unwrap()
    }

    fn assert_primitive_system<F: QpaCoefficientField>(
        algebra: &StructureConstantsAlgebra<F>,
        expected_count: usize,
        expected_right_ideal_dimension: usize,
    ) {
        assert!(algebra.cached_primitive_idempotents().is_none());
        let idempotents = algebra.primitive_idempotents().unwrap();
        assert_eq!(idempotents.len(), expected_count);
        assert!(algebra.cached_primitive_idempotents().is_some());
        assert_eq!(algebra.primitive_idempotents().unwrap(), idempotents);

        let mut sum = vec![algebra.field().zero(); algebra.dimension()];
        for (left_position, left) in idempotents.iter().enumerate() {
            assert_eq!(left.multiply(left).unwrap(), *left);
            assert_eq!(
                right_principal_subalgebra(algebra, left.coordinates())
                    .unwrap()
                    .algebra()
                    .dimension(),
                expected_right_ideal_dimension
            );
            for (right_position, right) in idempotents.iter().enumerate() {
                if left_position != right_position {
                    assert!(
                        left.multiply(right)
                            .unwrap()
                            .coordinates()
                            .iter()
                            .all(|entry| algebra.field().is_zero(entry))
                    );
                }
            }
            for (entry, coefficient) in sum.iter_mut().zip(left.coordinates()) {
                algebra.field().add_assign(entry, coefficient);
            }
        }
        assert_eq!(
            sum,
            algebra
                .find_multiplicative_identity()
                .unwrap()
                .unwrap()
                .coordinates()
        );
    }

    fn finite_matrix_case_json<F: QpaCoefficientField>(
        field: &F,
        size: usize,
        coefficient: impl Fn(&F::Element) -> i64,
    ) -> Value {
        let algebra = full_matrix_algebra(field.clone(), size);
        let e11_polynomial = algebra
            .hopf_minimum_polynomial(&algebra.basis_elements()[0])
            .unwrap();
        let irreducible_polynomial =
            (field.qpa_characteristic() == crate::QpaFieldCharacteristic::Prime(2) && size == 2)
                .then(|| {
                    let mut coordinates = vec![field.zero(); 4];
                    coordinates[1] = field.one();
                    coordinates[2] = field.one();
                    coordinates[3] = field.one();
                    algebra
                        .hopf_minimum_polynomial(&algebra.element(coordinates).unwrap())
                        .unwrap()
                });
        let left_identity = algebra.hopf_left_identity().unwrap().unwrap();
        let zero_divisor = algebra.hopf_zero_divisor().unwrap().unwrap();
        let singular = algebra.hopf_singular_idempotent().unwrap();
        let has_before = algebra.cached_primitive_idempotents().is_some();
        let idempotents = algebra.primitive_idempotents().unwrap();
        let second = algebra.primitive_idempotents().unwrap();
        let identity = algebra.find_multiplicative_identity().unwrap().unwrap();
        let mut sum = vec![field.zero(); algebra.dimension()];
        for idempotent in &idempotents {
            for (entry, coefficient) in sum.iter_mut().zip(idempotent.coordinates()) {
                field.add_assign(entry, coefficient);
            }
        }
        let pairwise_orthogonal = idempotents.iter().enumerate().all(|(left_position, left)| {
            idempotents
                .iter()
                .enumerate()
                .filter(|(right_position, _)| left_position != *right_position)
                .all(|(_, right)| {
                    left.multiply(right)
                        .unwrap()
                        .coordinates()
                        .iter()
                        .all(|entry| field.is_zero(entry))
                })
        });
        let mut result = json!({
            "dimension": algebra.dimension(),
            "e11_minpoly": e11_polynomial.iter().map(&coefficient).collect::<Vec<_>>(),
            "left_identity": left_identity == identity,
            "zero_divisor": {
                "nonzero": zero_divisor.iter().all(|element| {
                    element.coordinates().iter().any(|entry| !field.is_zero(entry))
                }),
                "product_zero": zero_divisor[0]
                    .multiply(&zero_divisor[1])
                    .unwrap()
                    .coordinates()
                    .iter()
                    .all(|entry| field.is_zero(entry))
            },
            "singular": {
                "idempotent": singular.multiply(&singular).unwrap() == singular,
                "right_dimension": right_principal_subalgebra(&algebra, singular.coordinates())
                    .unwrap()
                    .algebra()
                    .dimension()
            },
            "primitive": {
                "has_before": has_before,
                "has_after": algebra.cached_primitive_idempotents().is_some(),
                "cache_identical": idempotents == second,
                "count": idempotents.len(),
                "complete": sum == identity.coordinates(),
                "idempotent": idempotents.iter().all(|element| {
                    element.multiply(element).unwrap() == *element
                }),
                "orthogonal": pairwise_orthogonal,
                "right_dimensions": idempotents.iter().map(|element| {
                    right_principal_subalgebra(&algebra, element.coordinates())
                        .unwrap()
                        .algebra()
                        .dimension()
                }).collect::<Vec<_>>()
            }
        });
        if let Some(polynomial) = irreducible_polynomial {
            result["irreducible_minpoly"] =
                json!(polynomial.iter().map(&coefficient).collect::<Vec<_>>());
        }
        result
    }

    #[test]
    fn minimum_polynomial_and_left_identity_follow_qpa_power_span_and_row_solve() {
        let field3 = Zp::new(3);
        let matrices = full_matrix_algebra(field3.clone(), 2);
        let e11 = matrices.basis_elements()[0].clone();
        assert_eq!(
            matrices.hopf_minimum_polynomial(&e11).unwrap(),
            vec![field3.zero(), field3.neg(&field3.one()), field3.one()]
        );
        assert_eq!(
            matrices
                .hopf_left_identity()
                .unwrap()
                .unwrap()
                .coordinates(),
            &[field3.one(), field3.zero(), field3.zero(), field3.one()]
        );

        let binary_matrices = full_matrix_algebra(Z2, 2);
        let irreducible = binary_matrices
            .element(vec![Z2.zero(), Z2.one(), Z2.one(), Z2.one()])
            .unwrap();
        assert_eq!(
            binary_matrices
                .hopf_minimum_polynomial(&irreducible)
                .unwrap(),
            vec![Z2.one(), Z2.one(), Z2.one()]
        );
        assert!(
            zero_algebra(&field3)
                .hopf_left_identity()
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn finite_field_factors_follow_gap_square_free_and_common_degree_stages() {
        let elements2 = Z2.qpa_elements().unwrap().unwrap();
        let mut rng = StdRng::seed_from_u64(1);
        let repeated = vec![Z2.zero(), Z2.zero(), Z2.one()];
        let factors = factor_finite_field_polynomial(&Z2, &elements2, &repeated, &mut rng).unwrap();
        assert_eq!(factors, vec![vec![Z2.zero(), Z2.one()]; 2]);

        let irreducible = vec![Z2.one(), Z2.one(), Z2.one()];
        let reducible = polynomial_product(&Z2, &[Z2.one(), Z2.one()], &irreducible);
        let factors =
            factor_finite_field_polynomial(&Z2, &elements2, &reducible, &mut rng).unwrap();
        assert_eq!(factors.len(), 2);
        assert_eq!(
            factors
                .iter()
                .fold(vec![Z2.one()], |product, factor| polynomial_product(
                    &Z2, &product, factor
                )),
            reducible
        );

        let field5 = Zp::new(5);
        let elements5 = field5.qpa_elements().unwrap().unwrap();
        let two_quadratics = polynomial_product(
            &field5,
            &[field5.nth(2.into()), field5.zero(), field5.one()],
            &[field5.nth(3.into()), field5.zero(), field5.one()],
        );
        let factors =
            factor_finite_field_polynomial(&field5, &elements5, &two_quadratics, &mut rng).unwrap();
        assert_eq!(factors.len(), 2);
        assert_eq!(
            factors
                .iter()
                .fold(vec![field5.one()], |product, factor| polynomial_product(
                    &field5, &product, factor
                )),
            two_quadratics
        );
    }

    #[test]
    fn zero_divisor_and_recursive_corner_reproduce_qpa_hopf_invariants() {
        for algebra in [full_matrix_algebra(Z2, 2), full_matrix_algebra(Z2, 3)] {
            let pair = algebra.hopf_zero_divisor().unwrap().unwrap();
            assert!(
                pair.iter()
                    .all(|element| element.coordinates().iter().any(|entry| !Z2.is_zero(entry)))
            );
            assert!(
                pair[0]
                    .multiply(&pair[1])
                    .unwrap()
                    .coordinates()
                    .iter()
                    .all(|entry| Z2.is_zero(entry))
            );
            let singular = algebra.hopf_singular_idempotent().unwrap();
            assert_eq!(singular.multiply(&singular).unwrap(), singular);
            assert_eq!(
                right_principal_subalgebra(&algebra, singular.coordinates())
                    .unwrap()
                    .algebra()
                    .dimension(),
                algebra.dimension().isqrt()
            );
        }

        let commutative = product_of_fields(Z2, 2);
        assert!(commutative.hopf_zero_divisor().unwrap().is_none());
        assert_eq!(
            commutative.hopf_singular_idempotent().unwrap(),
            commutative.find_multiplicative_identity().unwrap().unwrap()
        );
    }

    #[test]
    fn primitive_idempotents_are_complete_orthogonal_cached_and_prime_field_generic() {
        assert_primitive_system(&full_matrix_algebra(Z2, 2), 2, 2);
        assert_primitive_system(&full_matrix_algebra(Z2, 3), 3, 3);
        assert_primitive_system(&full_matrix_algebra(Zp::new(3), 2), 2, 2);
        assert_primitive_system(&full_matrix_algebra(Zp::new(5), 2), 2, 2);
        assert_primitive_system(&product_of_fields(Z2, 1), 1, 1);
    }

    #[test]
    fn qpa_preconditions_and_commutative_early_return_are_explicit() {
        assert_eq!(
            dual_numbers(Z2).primitive_idempotents().unwrap_err(),
            HopfError::NotSemisimple
        );
        assert_eq!(
            product_of_fields(Z2, 2)
                .primitive_idempotents()
                .unwrap_err(),
            HopfError::NotSimple
        );
        assert_eq!(
            full_matrix_algebra(Q, 2)
                .primitive_idempotents()
                .unwrap_err(),
            HopfError::FiniteAlgebraRequired
        );
        assert!(
            product_of_fields(Q, 2)
                .hopf_zero_divisor()
                .unwrap()
                .is_none()
        );

        let mut rng = StdRng::seed_from_u64(HOPF_DEFAULT_SEED);
        assert!(
            full_matrix_algebra(Z2, 2)
                .hopf_zero_divisor_with_rng(&mut rng)
                .unwrap()
                .is_some()
        );
    }

    #[test]
    fn pinned_qpa_idempotent_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/idempotent.schema-v1.json"
        ))
        .unwrap();
        let binary = |coefficient: &<Z2 as numerica::domains::Set>::Element| {
            Z2.to_integer(coefficient)
                .to_string()
                .parse::<i64>()
                .unwrap()
        };
        let field3 = Zp::new(3);
        let ternary = |coefficient: &<Zp as numerica::domains::Set>::Element| {
            field3
                .to_integer(coefficient)
                .to_string()
                .parse::<i64>()
                .unwrap()
        };
        let commutative_q = product_of_fields(Q, 2);
        let commutative_f2 = product_of_fields(Z2, 2);
        let simple_f2 = product_of_fields(Z2, 1);
        let actual = json!({
            "schema_version": 1,
            "kind": "idempotent",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "f2_m2": finite_matrix_case_json(&Z2, 2, binary),
            "f3_m2": finite_matrix_case_json(&field3, 2, ternary),
            "f2_m3": finite_matrix_case_json(&Z2, 3, binary),
            "commutative": {
                "rational_zero_divisor_empty": commutative_q
                    .hopf_zero_divisor().unwrap().is_none(),
                "binary_zero_divisor_empty": commutative_f2
                    .hopf_zero_divisor().unwrap().is_none(),
                "rational_singular_is_one": commutative_q
                    .hopf_singular_idempotent().unwrap()
                    == commutative_q.find_multiplicative_identity().unwrap().unwrap(),
                "simple_binary_primitive_count": simple_f2
                    .primitive_idempotents().unwrap().len(),
                "zero_left_identity_is_fail": zero_algebra(&Zp::new(3))
                    .hopf_left_identity().unwrap().is_none()
            },
            "errors": {
                "nonsemisimple_success": dual_numbers(Z2).primitive_idempotents().is_ok(),
                "nonsimple_success": product_of_fields(Z2, 2)
                    .primitive_idempotents().is_ok(),
                "nonfinite_success": full_matrix_algebra(Q, 2)
                    .primitive_idempotents().is_ok(),
                "nonfinite_zero_divisor_success": full_matrix_algebra(Q, 2)
                    .hopf_zero_divisor().is_ok()
            }
        });
        assert_eq!(actual, expected);
    }
}
