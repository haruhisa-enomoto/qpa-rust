//! Exact finite extension fields compatible with Numerica's field contract.
//!
//! The reference representation is the polynomial basis of
//! `F_p[x] / (f(x))`. GAP's designated generator is the residue class of `x`
//! when `f` is the Conway polynomial. This module intentionally starts with a
//! correctness-first representation; optimized Zech logarithms can be placed
//! behind the same field contract later.

use std::{
    cmp::Ordering,
    collections::HashMap,
    fmt::{Debug, Display, Error as FormatError, Formatter, Write},
    hash::Hash,
    str::FromStr,
    sync::{Mutex, OnceLock},
};

use numerica::{
    domains::{
        EuclideanDomain, Field, InternalOrdering, Ring, RingOps, Set,
        finite_field::{FiniteFieldWorkspace, PrimeIteratorU64},
        integer::Integer,
    },
    printer::{PrintOptions, PrintState},
};
use rand::RngCore;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// How the defining polynomial of a finite extension field was selected.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FiniteFieldConstruction {
    /// The polynomial and generator use GAP's Conway convention.
    GapConway,
    /// The caller supplied a primitive irreducible polynomial.
    ExplicitPrimitivePolynomial,
}

/// Failure to construct an exact finite extension field.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum GaloisFieldError {
    /// The characteristic must be prime.
    #[error("{characteristic} is not a prime characteristic")]
    CharacteristicNotPrime {
        /// Rejected characteristic.
        characteristic: u32,
    },
    /// A local extension-field context represents a nontrivial extension.
    #[error("the extension degree must be at least two")]
    DegreeTooSmall,
    /// Numerica's integer-power interface accepts at most a `u32` exponent.
    #[error("the extension degree is too large for exact field construction")]
    DegreeTooLarge,
    /// A modulus for degree `n` has exactly `n + 1` coefficients.
    #[error("a degree-{degree} modulus needs {expected} coefficients, not {actual}")]
    ModulusLength {
        /// Requested extension degree.
        degree: usize,
        /// Required coefficient count.
        expected: usize,
        /// Supplied coefficient count.
        actual: usize,
    },
    /// Every coefficient is represented canonically in `[0, p)`.
    #[error("modulus coefficient {coefficient} is not in [0, {characteristic})")]
    CoefficientOutOfRange {
        /// Rejected coefficient.
        coefficient: u32,
        /// Field characteristic.
        characteristic: u32,
    },
    /// Polynomial-basis reduction currently requires a monic modulus.
    #[error("the defining polynomial must be monic")]
    ModulusNotMonic,
    /// A reducible quotient is not a field.
    #[error("the defining polynomial is reducible over GF({characteristic})")]
    ModulusReducible {
        /// Field characteristic.
        characteristic: u32,
    },
    /// GAP's Conway root must generate the nonzero multiplicative group.
    #[error("the residue class of x is not primitive in GF({characteristic}^{degree})")]
    GeneratorNotPrimitive {
        /// Field characteristic.
        characteristic: u32,
        /// Extension degree.
        degree: usize,
    },
    /// An element must have one coordinate for every polynomial-basis vector.
    #[error("GF({characteristic}^{degree}) needs {degree} coordinates, not {actual}")]
    ElementLength {
        /// Field characteristic.
        characteristic: u32,
        /// Extension degree.
        degree: usize,
        /// Supplied coordinate count.
        actual: usize,
    },
    /// Canonical finite-field embeddings exist exactly along the Conway tower.
    #[error(
        "GF({source_characteristic}^{source_degree}) is not a canonical subfield of GF({target_characteristic}^{target_degree})"
    )]
    IncompatibleSubfield {
        /// Source characteristic.
        source_characteristic: u32,
        /// Source extension degree.
        source_degree: usize,
        /// Target characteristic.
        target_characteristic: u32,
        /// Target extension degree.
        target_degree: usize,
    },
}

/// An element of a [`GaloisField`] in its canonical polynomial basis.
///
/// The containing field performs arithmetic, following Numerica's separation
/// between ring contexts and context-free element values.
#[derive(Clone, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct GaloisFieldElement(Box<[u32]>);

impl GaloisFieldElement {
    /// Return coefficients from the constant coefficient upwards.
    #[must_use]
    pub fn coefficients(&self) -> &[u32] {
        &self.0
    }
}

impl Debug for GaloisFieldElement {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("GaloisFieldElement").field(&self.0).finish()
    }
}

impl Display for GaloisFieldElement {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "[")?;
        for (position, coefficient) in self.0.iter().enumerate() {
            if position != 0 {
                write!(f, ",")?;
            }
            write!(f, "{coefficient}")?;
        }
        write!(f, "]")
    }
}

impl InternalOrdering for GaloisFieldElement {
    fn internal_cmp(&self, other: &Self) -> Ordering {
        self.0.as_ref().cmp(other.0.as_ref())
    }
}

/// The exact field `GF(p^n)` represented as `F_p[x] / (f(x))`.
#[derive(Clone, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct GaloisField {
    characteristic: u32,
    degree: usize,
    modulus: Box<[u32]>,
    construction: FiniteFieldConstruction,
}

impl GaloisField {
    /// Construct GAP's canonical model of `GF(p^n)`.
    ///
    /// This uses the complete precomputed Conway database shipped by the
    /// pinned GAP revision. If that database has no entry, it follows GAP's
    /// deterministic candidate ordering and compatible-subfield test to
    /// construct the Conway polynomial on demand.
    ///
    /// # Errors
    ///
    /// Rejects nonprime characteristics or an invalid extension degree.
    pub fn conway(characteristic: u32, degree: usize) -> Result<Self, GaloisFieldError> {
        validate_characteristic_degree(characteristic, degree)?;
        let cache = gap_conway_field_cache();
        if let Some(field) = cache
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .get(&(characteristic, degree))
            .cloned()
        {
            return Ok(field);
        }

        let modulus = if let Some(encoded) = gap_conway_data_encoding(characteristic, degree) {
            decode_gap_conway_integer(characteristic, degree, encoded)
        } else {
            let small_order = u32::try_from(degree)
                .ok()
                .and_then(|degree| u64::from(characteristic).checked_pow(degree));
            if let Some(encoded) = small_order.and_then(gap_small_conway_encoding) {
                decode_gap_conway(characteristic, degree, encoded)
            } else {
                construct_gap_conway_polynomial(characteristic, degree)?
            }
        };
        // Both the imported entries and the fallback algorithm are the
        // defining data used by GAP. Avoid refactoring and factoring p^n-1 a
        // second time merely to revalidate that trusted result.
        let field = Self {
            characteristic,
            degree,
            modulus: modulus.into_boxed_slice(),
            construction: FiniteFieldConstruction::GapConway,
        };
        cache
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .insert((characteristic, degree), field.clone());
        Ok(field)
    }

    /// Construct `GF(p^n)` from a monic primitive irreducible polynomial.
    ///
    /// Coefficients are ordered from the constant coefficient upwards. The
    /// residue class of `x` is validated as a multiplicative generator.
    ///
    /// # Errors
    ///
    /// Rejects invalid characteristics, malformed or reducible polynomials,
    /// and polynomials whose residue class of `x` is not primitive.
    pub fn from_primitive_polynomial(
        characteristic: u32,
        modulus: Vec<u32>,
    ) -> Result<Self, GaloisFieldError> {
        Self::from_primitive_polynomial_with_construction(
            characteristic,
            modulus,
            FiniteFieldConstruction::ExplicitPrimitivePolynomial,
        )
    }

    fn from_primitive_polynomial_with_construction(
        characteristic: u32,
        modulus: Vec<u32>,
        construction: FiniteFieldConstruction,
    ) -> Result<Self, GaloisFieldError> {
        let degree = modulus.len().saturating_sub(1);
        validate_characteristic_degree(characteristic, degree)?;
        let expected = degree + 1;
        if modulus.len() != expected {
            return Err(GaloisFieldError::ModulusLength {
                degree,
                expected,
                actual: modulus.len(),
            });
        }
        if let Some(&coefficient) = modulus
            .iter()
            .find(|&&coefficient| coefficient >= characteristic)
        {
            return Err(GaloisFieldError::CoefficientOutOfRange {
                coefficient,
                characteristic,
            });
        }
        if modulus.last() != Some(&1) {
            return Err(GaloisFieldError::ModulusNotMonic);
        }
        if !polynomial_is_irreducible(characteristic, &modulus) {
            return Err(GaloisFieldError::ModulusReducible { characteristic });
        }
        let field = Self {
            characteristic,
            degree,
            modulus: modulus.into_boxed_slice(),
            construction,
        };
        if !field.indeterminate_is_primitive() {
            return Err(GaloisFieldError::GeneratorNotPrimitive {
                characteristic,
                degree,
            });
        }
        Ok(field)
    }

    /// Return the prime characteristic `p`.
    #[must_use]
    pub const fn prime(&self) -> u32 {
        self.characteristic
    }

    /// Return the extension degree `n` in `GF(p^n)`.
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.degree
    }

    /// Return `p^n` as an arbitrary-precision integer.
    #[must_use]
    pub fn order(&self) -> Integer {
        let exponent = u64::try_from(self.degree).unwrap_or(u64::MAX);
        Integer::from(self.characteristic).pow(exponent)
    }

    /// Return the monic defining polynomial, constant coefficient first.
    #[must_use]
    pub fn modulus(&self) -> &[u32] {
        &self.modulus
    }

    /// Return how the defining polynomial was selected.
    #[must_use]
    pub const fn construction(&self) -> FiniteFieldConstruction {
        self.construction
    }

    /// Return the designated primitive generator, the residue class of `x`.
    #[must_use]
    pub fn primitive_element(&self) -> GaloisFieldElement {
        let mut polynomial = vec![0, 1];
        reduce_polynomial(self.characteristic, &self.modulus, &mut polynomial);
        polynomial.resize(self.degree, 0);
        GaloisFieldElement(polynomial.into_boxed_slice())
    }

    /// Construct an element from polynomial-basis coefficients.
    ///
    /// # Errors
    ///
    /// Rejects the wrong number of coefficients or a noncanonical coefficient.
    pub fn element(&self, coefficients: Vec<u32>) -> Result<GaloisFieldElement, GaloisFieldError> {
        if coefficients.len() != self.degree {
            return Err(GaloisFieldError::ElementLength {
                characteristic: self.characteristic,
                degree: self.degree,
                actual: coefficients.len(),
            });
        }
        if let Some(&coefficient) = coefficients
            .iter()
            .find(|&&coefficient| coefficient >= self.characteristic)
        {
            return Err(GaloisFieldError::CoefficientOutOfRange {
                coefficient,
                characteristic: self.characteristic,
            });
        }
        Ok(GaloisFieldElement(coefficients.into_boxed_slice()))
    }

    /// Apply the `p^power`-Frobenius automorphism.
    #[must_use]
    pub fn frobenius(&self, element: &GaloisFieldElement, power: usize) -> GaloisFieldElement {
        let reduced_power = power % self.degree;
        let mut result = element.clone();
        for _ in 0..reduced_power {
            result = self.pow(&result, u64::from(self.characteristic));
        }
        result
    }

    /// Embed an element of a canonical Conway subfield into this field.
    ///
    /// GAP sends the designated generator of `GF(p^d)` to
    /// `Z(p^n)^((p^n-1)/(p^d-1))` when `d` divides `n`.
    ///
    /// # Errors
    ///
    /// Rejects different characteristics, a nondividing source degree, or an
    /// explicitly supplied (non-Conway) field model.
    pub fn embed_subfield_element(
        &self,
        source: &Self,
        element: &GaloisFieldElement,
    ) -> Result<GaloisFieldElement, GaloisFieldError> {
        if self.construction != FiniteFieldConstruction::GapConway
            || source.construction != FiniteFieldConstruction::GapConway
            || self.characteristic != source.characteristic
            || !self.degree.is_multiple_of(source.degree)
        {
            return Err(GaloisFieldError::IncompatibleSubfield {
                source_characteristic: source.characteristic,
                source_degree: source.degree,
                target_characteristic: self.characteristic,
                target_degree: self.degree,
            });
        }
        let exponent = (self.order() - 1) / (source.order() - 1);
        let embedded_generator = self.pow_integer(&self.primitive_element(), &exponent);
        let mut result = self.zero();
        for &coefficient in element.coefficients().iter().rev() {
            result = self.mul(&result, &embedded_generator);
            if coefficient != 0 {
                result = self.add(&result, &self.nth(Integer::from(coefficient)));
            }
        }
        Ok(result)
    }

    /// Enumerate zero followed by powers of GAP's designated generator.
    ///
    /// Returns `None` when the cardinality does not fit in this process's
    /// address space.
    #[must_use]
    pub fn elements_gap_order(&self) -> Option<Vec<GaloisFieldElement>> {
        let order = usize::try_from(self.order().to_u64()?).ok()?;
        let mut elements = Vec::with_capacity(order);
        elements.push(self.zero());
        let generator = self.primitive_element();
        let mut power = self.one();
        for _ in 1..order {
            elements.push(power.clone());
            power = self.mul(&power, &generator);
        }
        Some(elements)
    }

    fn add_elements(
        &self,
        left: &GaloisFieldElement,
        right: &GaloisFieldElement,
    ) -> GaloisFieldElement {
        debug_assert_eq!(left.0.len(), self.degree);
        debug_assert_eq!(right.0.len(), self.degree);
        GaloisFieldElement(
            left.0
                .iter()
                .zip(right.0.iter())
                .map(|(&left, &right)| add_mod(left, right, self.characteristic))
                .collect::<Vec<_>>()
                .into_boxed_slice(),
        )
    }

    fn sub_elements(
        &self,
        left: &GaloisFieldElement,
        right: &GaloisFieldElement,
    ) -> GaloisFieldElement {
        GaloisFieldElement(
            left.0
                .iter()
                .zip(right.0.iter())
                .map(|(&left, &right)| sub_mod(left, right, self.characteristic))
                .collect::<Vec<_>>()
                .into_boxed_slice(),
        )
    }

    fn mul_elements(
        &self,
        left: &GaloisFieldElement,
        right: &GaloisFieldElement,
    ) -> GaloisFieldElement {
        let mut product = vec![0; self.degree.saturating_mul(2).saturating_sub(1)];
        for (left_degree, &left_coefficient) in left.0.iter().enumerate() {
            for (right_degree, &right_coefficient) in right.0.iter().enumerate() {
                let position = left_degree + right_degree;
                product[position] = add_mod(
                    product[position],
                    mul_mod(left_coefficient, right_coefficient, self.characteristic),
                    self.characteristic,
                );
            }
        }
        reduce_polynomial(self.characteristic, &self.modulus, &mut product);
        product.resize(self.degree, 0);
        GaloisFieldElement(product.into_boxed_slice())
    }

    fn pow_integer(&self, base: &GaloisFieldElement, exponent: &Integer) -> GaloisFieldElement {
        debug_assert!(exponent >= &Integer::zero());
        let mut exponent = exponent.clone();
        let mut base = base.clone();
        let mut result = self.one();
        while exponent != 0 {
            if &exponent % 2 == 1 {
                result = self.mul(&result, &base);
            }
            exponent /= 2;
            if exponent != 0 {
                base = self.mul(&base, &base);
            }
        }
        result
    }

    fn indeterminate_is_primitive(&self) -> bool {
        let order: Integer = self.order() - Integer::one();
        let generator = self.primitive_element();
        for (prime, _) in order.factor() {
            if self.is_one(&self.pow_integer(&generator, &(&order / prime))) {
                return false;
            }
        }
        true
    }

    fn format_element<W: Write>(
        &self,
        element: &GaloisFieldElement,
        output: &mut W,
    ) -> Result<(), FormatError> {
        if self.is_zero(element) {
            return output.write_char('0');
        }
        let mut first = true;
        for (degree, &coefficient) in element.0.iter().enumerate() {
            if coefficient == 0 {
                continue;
            }
            if !first {
                output.write_char('+')?;
            }
            if coefficient != 1 || degree == 0 {
                write!(output, "{coefficient}")?;
                if degree > 0 {
                    output.write_char('*')?;
                }
            }
            if degree == 1 {
                output.write_char('z')?;
            } else if degree > 1 {
                write!(output, "z^{degree}")?;
            }
            first = false;
        }
        Ok(())
    }
}

impl Display for GaloisField {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "GF({}^{})", self.characteristic, self.degree)
    }
}

impl Set for GaloisField {
    type Element = GaloisFieldElement;

    fn size(&self) -> Option<Integer> {
        Some(self.order())
    }
}

macro_rules! impl_ring_ops {
    ($element:ty) => {
        impl RingOps<$element> for GaloisField {
            fn add(&self, left: $element, right: $element) -> Self::Element {
                self.add_elements(&left, &right)
            }

            fn sub(&self, left: $element, right: $element) -> Self::Element {
                self.sub_elements(&left, &right)
            }

            fn mul(&self, left: $element, right: $element) -> Self::Element {
                self.mul_elements(&left, &right)
            }

            fn neg(&self, element: $element) -> Self::Element {
                self.sub_elements(&self.zero(), &element)
            }

            fn add_assign(&self, left: &mut Self::Element, right: $element) {
                *left = self.add_elements(left, &right);
            }

            fn sub_assign(&self, left: &mut Self::Element, right: $element) {
                *left = self.sub_elements(left, &right);
            }

            fn mul_assign(&self, left: &mut Self::Element, right: $element) {
                *left = self.mul_elements(left, &right);
            }

            fn add_mul_assign(
                &self,
                accumulator: &mut Self::Element,
                left: $element,
                right: $element,
            ) {
                let product = self.mul_elements(&left, &right);
                *accumulator = self.add_elements(accumulator, &product);
            }

            fn sub_mul_assign(
                &self,
                accumulator: &mut Self::Element,
                left: $element,
                right: $element,
            ) {
                let product = self.mul_elements(&left, &right);
                *accumulator = self.sub_elements(accumulator, &product);
            }
        }
    };
}

impl_ring_ops!(GaloisFieldElement);
impl_ring_ops!(&GaloisFieldElement);

impl Ring for GaloisField {
    fn zero(&self) -> Self::Element {
        GaloisFieldElement(vec![0; self.degree].into_boxed_slice())
    }

    fn one(&self) -> Self::Element {
        let mut coefficients = vec![0; self.degree];
        coefficients[0] = 1;
        GaloisFieldElement(coefficients.into_boxed_slice())
    }

    fn nth(&self, integer: Integer) -> Self::Element {
        let characteristic = Integer::from(self.characteristic);
        let mut coefficient = integer % &characteristic;
        if coefficient < 0 {
            coefficient += characteristic;
        }
        let mut coefficients = vec![0; self.degree];
        coefficients[0] = u32::try_from(coefficient.to_u64().expect("residue fits u32"))
            .expect("residue fits u32");
        GaloisFieldElement(coefficients.into_boxed_slice())
    }

    fn pow(&self, base: &Self::Element, mut exponent: u64) -> Self::Element {
        let mut base = base.clone();
        let mut result = self.one();
        while exponent != 0 {
            if exponent % 2 == 1 {
                result = self.mul(&result, &base);
            }
            exponent /= 2;
            if exponent != 0 {
                base = self.mul(&base, &base);
            }
        }
        result
    }

    fn is_zero(&self, element: &Self::Element) -> bool {
        element.0.iter().all(|&coefficient| coefficient == 0)
    }

    fn is_one(&self, element: &Self::Element) -> bool {
        element.0.first() == Some(&1)
            && element
                .0
                .iter()
                .skip(1)
                .all(|&coefficient| coefficient == 0)
    }

    fn one_is_gcd_unit() -> bool {
        true
    }

    fn characteristic(&self) -> Integer {
        self.characteristic.into()
    }

    fn try_inv(&self, element: &Self::Element) -> Option<Self::Element> {
        if self.is_zero(element) {
            return None;
        }
        Some(self.pow_integer(element, &(self.order() - 2)))
    }

    fn try_div(&self, left: &Self::Element, right: &Self::Element) -> Option<Self::Element> {
        self.try_inv(right).map(|inverse| self.mul(left, &inverse))
    }

    fn sample(&self, rng: &mut impl RngCore, _range: (i64, i64)) -> Self::Element {
        GaloisFieldElement(
            (0..self.degree)
                .map(|_| {
                    u32::try_from(rng.next_u64() % u64::from(self.characteristic))
                        .expect("a residue modulo a u32 characteristic fits u32")
                })
                .collect::<Vec<_>>()
                .into_boxed_slice(),
        )
    }

    fn format<W: Write>(
        &self,
        element: &Self::Element,
        _options: &PrintOptions,
        _state: PrintState,
        output: &mut W,
    ) -> Result<bool, FormatError> {
        self.format_element(element, output)?;
        Ok(false)
    }

    fn format_ring<W: Write>(
        &self,
        _options: &PrintOptions,
        _state: PrintState,
        output: &mut W,
    ) -> Result<bool, FormatError> {
        write!(output, "{self}")?;
        Ok(false)
    }
}

impl EuclideanDomain for GaloisField {
    fn rem(&self, _left: &Self::Element, _right: &Self::Element) -> Self::Element {
        self.zero()
    }

    fn quot_rem(
        &self,
        left: &Self::Element,
        right: &Self::Element,
    ) -> (Self::Element, Self::Element) {
        (self.div(left, right), self.zero())
    }

    fn gcd(&self, _left: &Self::Element, _right: &Self::Element) -> Self::Element {
        self.one()
    }
}

impl Field for GaloisField {
    fn div(&self, left: &Self::Element, right: &Self::Element) -> Self::Element {
        self.mul(left, &self.inv(right))
    }

    fn div_assign(&self, left: &mut Self::Element, right: &Self::Element) {
        *left = self.div(left, right);
    }

    fn inv(&self, element: &Self::Element) -> Self::Element {
        self.try_inv(element)
            .unwrap_or_else(|| panic!("zero is not invertible in {self}"))
    }
}

fn validate_characteristic_degree(
    characteristic: u32,
    degree: usize,
) -> Result<(), GaloisFieldError> {
    if PrimeIteratorU64::new(u64::from(characteristic).saturating_sub(1)).next()
        != Some(u64::from(characteristic))
    {
        return Err(GaloisFieldError::CharacteristicNotPrime { characteristic });
    }
    if degree < 2 {
        return Err(GaloisFieldError::DegreeTooSmall);
    }
    if u32::try_from(degree).is_err() {
        return Err(GaloisFieldError::DegreeTooLarge);
    }
    Ok(())
}

fn add_mod(left: u32, right: u32, modulus: u32) -> u32 {
    u32::try_from((u64::from(left) + u64::from(right)) % u64::from(modulus))
        .expect("a residue modulo a u32 modulus fits u32")
}

fn sub_mod(left: u32, right: u32, modulus: u32) -> u32 {
    u32::try_from((u64::from(left) + u64::from(modulus) - u64::from(right)) % u64::from(modulus))
        .expect("a residue modulo a u32 modulus fits u32")
}

fn mul_mod(left: u32, right: u32, modulus: u32) -> u32 {
    u32::try_from((u64::from(left) * u64::from(right)) % u64::from(modulus))
        .expect("a residue modulo a u32 modulus fits u32")
}

fn reduce_polynomial(characteristic: u32, modulus: &[u32], polynomial: &mut Vec<u32>) {
    let degree = modulus.len() - 1;
    while polynomial.len() > degree {
        let leading_degree = polynomial.len() - 1;
        let leading = polynomial[leading_degree];
        if leading != 0 {
            let shift = leading_degree - degree;
            for (position, &coefficient) in modulus[..degree].iter().enumerate() {
                let reduction = mul_mod(leading, coefficient, characteristic);
                polynomial[shift + position] =
                    sub_mod(polynomial[shift + position], reduction, characteristic);
            }
        }
        polynomial.pop();
    }
    trim_polynomial(polynomial);
}

fn trim_polynomial(polynomial: &mut Vec<u32>) {
    while polynomial.last() == Some(&0) {
        polynomial.pop();
    }
}

fn polynomial_sub(characteristic: u32, left: &[u32], right: &[u32]) -> Vec<u32> {
    let mut result = vec![0; left.len().max(right.len())];
    for (position, coefficient) in result.iter_mut().enumerate() {
        *coefficient = sub_mod(
            left.get(position).copied().unwrap_or(0),
            right.get(position).copied().unwrap_or(0),
            characteristic,
        );
    }
    trim_polynomial(&mut result);
    result
}

fn polynomial_multiply_mod(
    characteristic: u32,
    modulus: &[u32],
    left: &[u32],
    right: &[u32],
) -> Vec<u32> {
    if left.is_empty() || right.is_empty() {
        return Vec::new();
    }
    let mut result = vec![0; left.len() + right.len() - 1];
    for (left_degree, &left_coefficient) in left.iter().enumerate() {
        for (right_degree, &right_coefficient) in right.iter().enumerate() {
            let position = left_degree + right_degree;
            result[position] = add_mod(
                result[position],
                mul_mod(left_coefficient, right_coefficient, characteristic),
                characteristic,
            );
        }
    }
    reduce_polynomial(characteristic, modulus, &mut result);
    result
}

fn polynomial_power_mod(
    characteristic: u32,
    modulus: &[u32],
    base: &[u32],
    mut exponent: u64,
) -> Vec<u32> {
    let mut base = base.to_vec();
    let mut result = vec![1];
    while exponent != 0 {
        if exponent % 2 == 1 {
            result = polynomial_multiply_mod(characteristic, modulus, &result, &base);
        }
        exponent /= 2;
        if exponent != 0 {
            base = polynomial_multiply_mod(characteristic, modulus, &base, &base);
        }
    }
    result
}

fn polynomial_power_mod_integer(
    characteristic: u32,
    modulus: &[u32],
    base: &[u32],
    exponent: &Integer,
) -> Vec<u32> {
    let mut exponent = exponent.clone();
    let mut base = base.to_vec();
    let mut result = vec![1];
    while exponent != 0 {
        if &exponent % 2 == 1 {
            result = polynomial_multiply_mod(characteristic, modulus, &result, &base);
        }
        exponent /= 2;
        if exponent != 0 {
            base = polynomial_multiply_mod(characteristic, modulus, &base, &base);
        }
    }
    result
}

fn polynomial_evaluate_mod(
    characteristic: u32,
    modulus: &[u32],
    polynomial: &[u32],
    value: &[u32],
) -> Vec<u32> {
    let mut result = Vec::new();
    for &coefficient in polynomial.iter().rev() {
        result = polynomial_multiply_mod(characteristic, modulus, &result, value);
        if coefficient != 0 {
            if result.is_empty() {
                result.push(coefficient);
            } else {
                result[0] = add_mod(result[0], coefficient, characteristic);
                trim_polynomial(&mut result);
            }
        }
    }
    result
}

fn inverse_mod_prime(value: u32, prime: u32) -> u32 {
    let mut exponent = u64::from(prime - 2);
    let mut base = value;
    let mut result = 1;
    while exponent != 0 {
        if exponent % 2 == 1 {
            result = mul_mod(result, base, prime);
        }
        exponent /= 2;
        if exponent != 0 {
            base = mul_mod(base, base, prime);
        }
    }
    result
}

fn polynomial_remainder(characteristic: u32, dividend: &[u32], divisor: &[u32]) -> Vec<u32> {
    let mut remainder = dividend.to_vec();
    trim_polynomial(&mut remainder);
    let mut divisor = divisor.to_vec();
    trim_polynomial(&mut divisor);
    assert!(!divisor.is_empty(), "polynomial division by zero");
    let divisor_degree = divisor.len() - 1;
    let leading_inverse = inverse_mod_prime(divisor[divisor_degree], characteristic);
    while remainder.len() >= divisor.len() {
        let shift = remainder.len() - divisor.len();
        let factor = mul_mod(
            *remainder.last().expect("nonempty remainder"),
            leading_inverse,
            characteristic,
        );
        for (position, &coefficient) in divisor.iter().enumerate() {
            let index = shift + position;
            remainder[index] = sub_mod(
                remainder[index],
                mul_mod(factor, coefficient, characteristic),
                characteristic,
            );
        }
        trim_polynomial(&mut remainder);
    }
    remainder
}

fn polynomial_gcd(characteristic: u32, left: &[u32], right: &[u32]) -> Vec<u32> {
    let mut left = left.to_vec();
    let mut right = right.to_vec();
    trim_polynomial(&mut left);
    trim_polynomial(&mut right);
    while !right.is_empty() {
        let remainder = polynomial_remainder(characteristic, &left, &right);
        left = right;
        right = remainder;
    }
    if let Some(&leading) = left.last() {
        let inverse = inverse_mod_prime(leading, characteristic);
        for coefficient in &mut left {
            *coefficient = mul_mod(*coefficient, inverse, characteristic);
        }
    }
    left
}

fn prime_divisors(mut value: usize) -> Vec<usize> {
    let mut result = Vec::new();
    let mut divisor = 2;
    while divisor <= value / divisor {
        if value.is_multiple_of(divisor) {
            result.push(divisor);
            while value.is_multiple_of(divisor) {
                value /= divisor;
            }
        }
        divisor += usize::from(divisor == 2) + 1;
    }
    if value > 1 {
        result.push(value);
    }
    result
}

fn polynomial_is_irreducible(characteristic: u32, modulus: &[u32]) -> bool {
    let degree = modulus.len() - 1;
    let indeterminate = vec![0, 1];
    for divisor in prime_divisors(degree) {
        let mut power = indeterminate.clone();
        for _ in 0..degree / divisor {
            power =
                polynomial_power_mod(characteristic, modulus, &power, u64::from(characteristic));
        }
        let difference = polynomial_sub(characteristic, &power, &indeterminate);
        if polynomial_gcd(characteristic, modulus, &difference).len() != 1 {
            return false;
        }
    }
    let mut power = indeterminate.clone();
    for _ in 0..degree {
        power = polynomial_power_mod(characteristic, modulus, &power, u64::from(characteristic));
    }
    polynomial_sub(characteristic, &power, &indeterminate).is_empty()
}

fn gap_conway_field_cache() -> &'static Mutex<HashMap<(u32, usize), GaloisField>> {
    static CACHE: OnceLock<Mutex<HashMap<(u32, usize), GaloisField>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

fn gap_conway_data_encoding(characteristic: u32, degree: usize) -> Option<&'static str> {
    static DATA: OnceLock<HashMap<(u32, usize), &'static str>> = OnceLock::new();
    DATA.get_or_init(|| {
        include_str!("conway_table.tsv")
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                let mut columns = line.split('\t');
                let prime = columns
                    .next()
                    .expect("generated Conway row has a characteristic")
                    .parse()
                    .expect("generated Conway characteristic is an integer");
                let degree = columns
                    .next()
                    .expect("generated Conway row has a degree")
                    .parse()
                    .expect("generated Conway degree is an integer");
                let encoded = columns
                    .next()
                    .expect("generated Conway row has encoded coefficients");
                assert!(
                    columns.next().is_none(),
                    "generated Conway row has three columns"
                );
                ((prime, degree), encoded)
            })
            .collect()
    })
    .get(&(characteristic, degree))
    .copied()
}

fn decode_gap_conway_integer(characteristic: u32, degree: usize, encoded: &str) -> Vec<u32> {
    let mut encoded = Integer::from_str(encoded).expect("generated Conway encoding is an integer");
    let prime = Integer::from(characteristic);
    let mut coefficients = Vec::with_capacity(degree + 1);
    for _ in 0..degree {
        coefficients.push(
            u32::try_from(
                (&encoded % &prime)
                    .to_u64()
                    .expect("Conway coefficient fits u64"),
            )
            .expect("Conway coefficient fits u32"),
        );
        encoded /= &prime;
    }
    debug_assert_eq!(encoded, 0);
    coefficients.push(1);
    coefficients
}

fn construct_gap_conway_polynomial(
    characteristic: u32,
    degree: usize,
) -> Result<Vec<u32>, GaloisFieldError> {
    let degree_prime_divisors = prime_divisors(degree);
    let subfield_degrees = degree_prime_divisors
        .iter()
        .map(|divisor| degree / divisor)
        .collect::<Vec<_>>();
    let subfield_polynomials = if subfield_degrees == [1] {
        Vec::new()
    } else {
        subfield_degrees
            .iter()
            .map(|&subdegree| {
                GaloisField::conway(characteristic, subdegree).map(|field| field.modulus)
            })
            .collect::<Result<Vec<_>, _>>()?
    };
    let order_minus_one: Integer =
        Integer::from(characteristic).pow(u64::try_from(degree).expect("validated degree")) - 1;
    let multiplicative_prime_divisors = order_minus_one
        .clone()
        .factor()
        .into_iter()
        .map(|(prime, _)| prime)
        .collect::<Vec<_>>();
    let compatibility_exponents = subfield_degrees
        .iter()
        .map(|&subdegree| {
            &order_minus_one
                / (Integer::from(characteristic)
                    .pow(u64::try_from(subdegree).expect("subdegree fits u64"))
                    - 1)
        })
        .collect::<Vec<_>>();

    let epsilon = if degree.is_multiple_of(2) {
        1
    } else {
        characteristic - 1
    };
    let primitive_root = smallest_primitive_root_mod_prime(characteristic);
    let mut candidate = vec![0; degree + 1];
    candidate[0] = mul_mod(epsilon, primitive_root, characteristic);
    candidate[degree] = 1;
    let indeterminate = [0, 1];

    loop {
        let primitive = polynomial_power_mod_integer(
            characteristic,
            &candidate,
            &indeterminate,
            &order_minus_one,
        ) == [1]
            && multiplicative_prime_divisors.iter().all(|prime| {
                polynomial_power_mod_integer(
                    characteristic,
                    &candidate,
                    &indeterminate,
                    &(&order_minus_one / prime),
                ) != [1]
            });
        let compatible = primitive
            && subfield_polynomials
                .iter()
                .zip(compatibility_exponents.iter())
                .all(|(subfield_polynomial, exponent)| {
                    let subfield_generator = polynomial_power_mod_integer(
                        characteristic,
                        &candidate,
                        &indeterminate,
                        exponent,
                    );
                    polynomial_evaluate_mod(
                        characteristic,
                        &candidate,
                        subfield_polynomial,
                        &subfield_generator,
                    )
                    .is_empty()
                });
        if compatible {
            return Ok(candidate);
        }
        next_gap_conway_candidate(characteristic, epsilon, &mut candidate);
    }
}

fn smallest_primitive_root_mod_prime(prime: u32) -> u32 {
    let order = Integer::from(prime - 1);
    let prime_divisors = order
        .clone()
        .factor()
        .into_iter()
        .map(|(factor, _)| factor)
        .collect::<Vec<_>>();
    (1..prime)
        .find(|&candidate| {
            prime_divisors.iter().all(|factor| {
                let exponent = (&order / factor)
                    .to_u64()
                    .expect("a prime-field exponent fits u64");
                modular_power(candidate, exponent, prime) != 1
            })
        })
        .expect("every finite prime field has a primitive root")
}

fn modular_power(mut base: u32, mut exponent: u64, modulus: u32) -> u32 {
    let mut result = 1;
    while exponent != 0 {
        if exponent % 2 == 1 {
            result = mul_mod(result, base, modulus);
        }
        exponent /= 2;
        if exponent != 0 {
            base = mul_mod(base, base, modulus);
        }
    }
    result
}

fn next_gap_conway_candidate(characteristic: u32, epsilon: u32, candidate: &mut [u32]) {
    let mut coefficient = 1;
    let mut increment = epsilon;
    loop {
        increment = sub_mod(0, increment, characteristic);
        candidate[coefficient] = add_mod(candidate[coefficient], increment, characteristic);
        if candidate[coefficient] != 0 {
            return;
        }
        coefficient += 1;
        assert!(
            coefficient + 1 < candidate.len(),
            "GAP candidate ordering must find a Conway polynomial before changing its leading term"
        );
    }
}

fn decode_gap_conway(characteristic: u32, degree: usize, mut encoded: u64) -> Vec<u32> {
    let mut coefficients = Vec::with_capacity(degree + 1);
    for _ in 0..degree {
        coefficients.push(
            u32::try_from(encoded % u64::from(characteristic))
                .expect("a residue modulo a u32 characteristic fits u32"),
        );
        encoded /= u64::from(characteristic);
    }
    coefficients.push(1);
    coefficients
}

// Source: pinned GAP `src/finfield_conway.h`, `PolsFF`.
// Keep the source expressions separate even when some evaluate to equal u64s.
#[allow(clippy::match_same_arms)]
fn gap_small_conway_encoding(order: u64) -> Option<u64> {
    Some(match order {
        4 | 8 | 16 => 1 + 2,
        32 => 1 + 4,
        64 => 1 + 2 + 8 + 16,
        128 => 1 + 2,
        256 => 1 + 4 + 8 + 16,
        512 => 1 + 16,
        1024 => 1 + 2 + 4 + 8 + 32 + 64,
        2048 => 1 + 4,
        4096 => 1 + 2 + 8 + 32 + 64 + 128,
        8192 => 1 + 2 + 8 + 16,
        16384 => 1 + 8 + 32 + 128,
        32768 => 1 + 4 + 16 + 32,
        65536 => 1 + 4 + 8 + 32,
        9 => 2 + 2 * 3,
        27 => 1 + 2 * 3,
        81 => 2 + 2 * 27,
        243 => 1 + 2 * 3,
        729 => 2 + 2 * 3 + 9 + 2 * 81,
        2187 => 1 + 2 * 9,
        6561 => 2 + 2 * 3 + 2 * 9 + 81 + 2 * 243,
        19683 => 1 + 3 + 2 * 9 + 2 * 27,
        59049 => 2 + 3 + 2 * 81 + 2 * 243 + 2 * 729,
        25 => 2 + 4 * 5,
        125 => 3 + 3 * 5,
        625 => 2 + 4 * 5 + 4 * 25,
        3125 => 3 + 4 * 5,
        15625 => 2 + 25 + 4 * 125 + 625,
        49 => 3 + 6 * 7,
        343 => 4 + 6 * 49,
        2401 => 3 + 4 * 7 + 5 * 49,
        16807 => 4 + 7,
        121 => 2 + 7 * 11,
        1331 => 9 + 2 * 11,
        14641 => 2 + 10 * 11 + 8 * 121,
        169 => 2 + 12 * 13,
        2197 => 11 + 2 * 13,
        28561 => 2 + 12 * 13 + 3 * 169,
        289 => 3 + 16 * 17,
        4913 => 14 + 17,
        361 => 2 + 18 * 19,
        6859 => 17 + 4 * 19,
        529 => 5 + 21 * 23,
        12167 => 18 + 2 * 23,
        841 => 2 + 24 * 29,
        24389 => 27 + 2 * 29,
        961 => 3 + 29 * 31,
        29791 => 28 + 31,
        1369 => 2 + 33 * 37,
        50653 => 35 + 6 * 37,
        1681 => 6 + 38 * 41,
        1849 => 3 + 42 * 43,
        2209 => 5 + 45 * 47,
        2809 => 2 + 49 * 53,
        3481 => 2 + 58 * 59,
        3721 => 2 + 60 * 61,
        4489 => 2 + 63 * 67,
        5041 => 7 + 69 * 71,
        5329 => 5 + 70 * 73,
        6241 => 3 + 78 * 79,
        6889 => 2 + 82 * 83,
        7921 => 3 + 82 * 89,
        9409 => 5 + 96 * 97,
        10201 => 2 + 97 * 101,
        10609 => 5 + 102 * 103,
        11449 => 2 + 103 * 107,
        11881 => 6 + 108 * 109,
        12769 => 3 + 101 * 113,
        16129 => 3 + 126 * 127,
        17161 => 2 + 127 * 131,
        18769 => 3 + 131 * 137,
        19321 => 2 + 138 * 139,
        22201 => 2 + 145 * 149,
        22801 => 6 + 149 * 151,
        24649 => 5 + 152 * 157,
        26569 => 2 + 159 * 163,
        27889 => 5 + 166 * 167,
        29929 => 2 + 169 * 173,
        32041 => 2 + 172 * 179,
        32761 => 2 + 177 * 181,
        36481 => 19 + 190 * 191,
        37249 => 5 + 192 * 193,
        38809 => 2 + 192 * 197,
        39601 => 3 + 193 * 199,
        44521 => 2 + 207 * 211,
        49729 => 3 + 221 * 223,
        51529 => 2 + 220 * 227,
        52441 => 6 + 228 * 229,
        54289 => 3 + 232 * 233,
        57121 => 7 + 237 * 239,
        58081 => 7 + 238 * 241,
        63001 => 6 + 242 * 251,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use numerica::domains::{Field, Ring, RingOps, Set, finite_field::FiniteFieldWorkspace as _};
    use serde_json::Value;

    use super::{FiniteFieldConstruction, GaloisField, GaloisFieldError};

    #[test]
    fn gap_conway_polynomials_decode_in_three_characteristics() {
        let f4 = GaloisField::conway(2, 2).unwrap();
        assert_eq!(f4.modulus(), [1, 1, 1]);
        assert_eq!(f4.order(), 4);
        assert_eq!(f4.construction(), FiniteFieldConstruction::GapConway);

        let f9 = GaloisField::conway(3, 2).unwrap();
        assert_eq!(f9.modulus(), [2, 2, 1]);
        assert_eq!(f9.order(), 9);

        let f25 = GaloisField::conway(5, 2).unwrap();
        assert_eq!(f25.modulus(), [2, 4, 1]);
        assert_eq!(f25.order(), 25);

        // GF(2^17) lies beyond GAP's immediate-field kernel table and proves
        // that the imported library database is active.
        assert_eq!(
            GaloisField::conway(2, 17).unwrap().modulus(),
            [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        );
    }

    #[test]
    fn missing_database_entry_uses_gap_candidate_ordering() {
        // The pinned conwdat2.g has no stored degree-two entry for p=307.
        // Pinned GAP constructs x^2-x+5 on demand.
        assert_eq!(GaloisField::conway(307, 2).unwrap().modulus(), [5, 306, 1]);
    }

    #[test]
    fn gf4_arithmetic_uses_polynomial_reduction() {
        let field = GaloisField::conway(2, 2).unwrap();
        let z = field.primitive_element();
        let z_squared = field.mul(&z, &z);
        assert_eq!(z.coefficients(), [0, 1]);
        assert_eq!(z_squared.coefficients(), [1, 1]);
        assert_eq!(field.add(&z_squared, &z), field.one());
        assert_eq!(field.inv(&z), z_squared);
        assert_eq!(field.div(&field.one(), &z), z_squared);
    }

    #[test]
    fn every_nonzero_small_field_element_has_an_inverse() {
        for (prime, degree) in [(2, 3), (3, 2), (5, 2)] {
            let field = GaloisField::conway(prime, degree).unwrap();
            let elements = field.elements_gap_order().unwrap();
            assert_eq!(
                elements.len(),
                usize::try_from(field.size().unwrap().to_u64().unwrap()).unwrap()
            );
            assert!(field.is_zero(&elements[0]));
            for element in elements.iter().skip(1) {
                assert_eq!(field.mul(element, &field.inv(element)), field.one());
            }
        }
    }

    #[test]
    fn exhaustive_small_field_axioms_hold() {
        for (prime, degree) in [(2, 2), (2, 3), (3, 2)] {
            let field = GaloisField::conway(prime, degree).unwrap();
            let elements = field.elements_gap_order().unwrap();
            for left in &elements {
                for middle in &elements {
                    assert_eq!(field.add(left, middle), field.add(middle, left));
                    assert_eq!(field.mul(left, middle), field.mul(middle, left));
                    for right in &elements {
                        assert_eq!(
                            field.add(&field.add(left, middle), right),
                            field.add(left, &field.add(middle, right))
                        );
                        assert_eq!(
                            field.mul(&field.mul(left, middle), right),
                            field.mul(left, &field.mul(middle, right))
                        );
                        assert_eq!(
                            field.mul(left, &field.add(middle, right)),
                            field.add(&field.mul(left, middle), &field.mul(left, right))
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn frobenius_has_the_expected_order() {
        for (prime, degree) in [(2, 3), (3, 2), (5, 2)] {
            let field = GaloisField::conway(prime, degree).unwrap();
            for element in field.elements_gap_order().unwrap() {
                assert_eq!(field.frobenius(&element, degree), element);
                assert_eq!(
                    field.frobenius(&element, 1),
                    field.pow(&element, u64::from(prime))
                );
            }
        }
    }

    #[test]
    fn conway_subfield_embedding_uses_gap_generator_relation() {
        let source = GaloisField::conway(2, 2).unwrap();
        let target = GaloisField::conway(2, 4).unwrap();
        let embedded_generator = target
            .embed_subfield_element(&source, &source.primitive_element())
            .unwrap();
        assert_eq!(embedded_generator.coefficients(), [0, 1, 1, 0]);

        let source_elements = source.elements_gap_order().unwrap();
        for left in &source_elements {
            for right in &source_elements {
                assert_eq!(
                    target
                        .embed_subfield_element(&source, &source.add(left, right))
                        .unwrap(),
                    target.add(
                        &target.embed_subfield_element(&source, left).unwrap(),
                        &target.embed_subfield_element(&source, right).unwrap()
                    )
                );
                assert_eq!(
                    target
                        .embed_subfield_element(&source, &source.mul(left, right))
                        .unwrap(),
                    target.mul(
                        &target.embed_subfield_element(&source, left).unwrap(),
                        &target.embed_subfield_element(&source, right).unwrap()
                    )
                );
            }
        }
        assert!(matches!(
            GaloisField::conway(3, 2)
                .unwrap()
                .embed_subfield_element(&source, &source.one()),
            Err(GaloisFieldError::IncompatibleSubfield { .. })
        ));
    }

    #[test]
    fn scalar_arithmetic_matches_pinned_gap_conway_fixture() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/finite_extension_field_scalars.schema-v1.json"
        ))
        .unwrap();
        for case in fixture["cases"].as_array().unwrap() {
            let prime = u32::try_from(case["characteristic"].as_u64().unwrap()).unwrap();
            let degree = usize::try_from(case["degree"].as_u64().unwrap()).unwrap();
            let field = GaloisField::conway(prime, degree).unwrap();
            let generator = field.primitive_element();
            let expected_modulus = case["conway"]
                .as_array()
                .unwrap()
                .iter()
                .map(|value| u32::try_from(value.as_u64().unwrap()).unwrap())
                .collect::<Vec<_>>();
            assert_eq!(field.modulus(), expected_modulus);

            let mut power = field.one();
            let powers = (0..case["generator_order"].as_u64().unwrap())
                .map(|_| {
                    let coefficients = power.coefficients().to_vec();
                    power = field.mul(&power, &generator);
                    coefficients
                })
                .collect::<Vec<_>>();
            assert_eq!(serde_json::to_value(powers).unwrap(), case["powers"]);
            assert_eq!(power, field.one());
            assert_eq!(
                serde_json::to_value(field.inv(&generator).coefficients()).unwrap(),
                case["inverse_generator"]
            );
            assert_eq!(
                serde_json::to_value(field.frobenius(&generator, 1).coefficients()).unwrap(),
                case["frobenius_generator"]
            );
            for subfield in case["proper_extension_subfields"].as_array().unwrap() {
                let subdegree = usize::try_from(subfield["degree"].as_u64().unwrap()).unwrap();
                let source = GaloisField::conway(prime, subdegree).unwrap();
                let embedded = field
                    .embed_subfield_element(&source, &source.primitive_element())
                    .unwrap();
                assert_eq!(
                    serde_json::to_value(embedded.coefficients()).unwrap(),
                    subfield["generator"]
                );
            }
        }
    }

    #[test]
    fn malformed_or_nonprimitive_moduli_are_rejected() {
        assert_eq!(
            GaloisField::from_primitive_polynomial(4, vec![1, 1, 1]),
            Err(GaloisFieldError::CharacteristicNotPrime { characteristic: 4 })
        );
        assert_eq!(
            GaloisField::from_primitive_polynomial(2, vec![0, 0, 1]),
            Err(GaloisFieldError::ModulusReducible { characteristic: 2 })
        );
        assert_eq!(
            GaloisField::from_primitive_polynomial(3, vec![1, 0, 1]),
            Err(GaloisFieldError::GeneratorNotPrimitive {
                characteristic: 3,
                degree: 2,
            })
        );
    }

    #[test]
    fn exact_integer_embedding_reduces_modulo_the_characteristic() {
        let field = GaloisField::conway(3, 2).unwrap();
        assert_eq!(field.nth((-1).into()).coefficients(), [2, 0]);
        assert_eq!(field.nth(7.into()).coefficients(), [1, 0]);
    }
}
