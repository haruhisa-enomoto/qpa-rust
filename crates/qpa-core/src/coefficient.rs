//! Runtime descriptions of the initially supported exact coefficient fields.

use std::cmp::Ordering;

use numerica::domains::{
    Field, Ring, RingOps,
    finite_field::{FiniteFieldCore, PrimitiveRootIterator, Z2, Zp},
    integer::Integer,
    rational::RationalField,
};
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::finite_field::{FiniteFieldConstruction, GaloisField};

/// An exact coefficient field admitted by the Numerica-backed port.
///
/// This public descriptor is independent of Numerica's concrete Rust types so
/// that the QPA API can retain stable field identity and serialization.
#[derive(Clone, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum CoefficientDomain {
    /// The rational field.
    Rational,
    /// The prime field with the given characteristic.
    PrimeField {
        /// The prime characteristic, stored without a machine-word bound.
        characteristic: Box<str>,
    },
    /// A finite extension field represented by a primitive polynomial.
    FiniteExtensionField {
        /// The prime characteristic, stored without a descriptor word bound.
        characteristic: Box<str>,
        /// Degree over the prime field.
        degree: usize,
        /// Defining polynomial coefficients, constant coefficient first.
        modulus: Vec<Box<str>>,
        /// Whether the defining data uses GAP's Conway convention.
        gap_conway: bool,
    },
}

impl CoefficientDomain {
    /// Construct the rational field descriptor.
    #[must_use]
    pub const fn rational() -> Self {
        Self::Rational
    }

    /// Return an unchecked prime-field descriptor.
    ///
    /// Validation is deliberately separate because Numerica selects different
    /// concrete backends according to the size and value of the prime. The
    /// checked constructor will be added with that adapter.
    pub fn prime_field_unchecked(characteristic: impl Into<Box<str>>) -> Self {
        Self::PrimeField {
            characteristic: characteristic.into(),
        }
    }

    /// Describe a validated local finite extension field.
    #[must_use]
    pub fn finite_extension_field(field: &GaloisField) -> Self {
        Self::FiniteExtensionField {
            characteristic: field.prime().to_string().into(),
            degree: field.degree(),
            modulus: field
                .modulus()
                .iter()
                .map(|coefficient| coefficient.to_string().into())
                .collect(),
            gap_conway: field.construction() == FiniteFieldConstruction::GapConway,
        }
    }
}

/// Failure to construct or use an exact coefficient domain.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum CoefficientDomainError {
    /// The requested characteristic is not a prime admitted by the backend.
    #[error("{characteristic} is not a supported prime characteristic")]
    UnsupportedPrime {
        /// The rejected textual characteristic.
        characteristic: Box<str>,
    },
    /// GAP's algorithms enumerate all field elements, but the cardinality is
    /// not addressable by this Rust process.
    #[error("the finite field of size {size} is too large to enumerate")]
    FieldTooLarge {
        /// Rejected field size.
        size: u64,
    },
    /// A prime field must have a primitive generator of its multiplicative group.
    #[error("no primitive root was found for the prime field of size {size}")]
    MissingPrimitiveRoot {
        /// Prime field size.
        size: u64,
    },
    /// GAP's algorithms enumerate all field elements, but an extension-field
    /// cardinality is not addressable by this Rust process.
    #[error("the finite field of size {size} is too large to enumerate")]
    ExtensionFieldTooLarge {
        /// Rejected exact field size.
        size: Box<str>,
    },
}

/// Characteristic information used by the pinned GAP algebra algorithms.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum QpaFieldCharacteristic {
    /// Characteristic zero.
    Zero,
    /// Positive prime characteristic.
    Prime(u64),
}

/// Exact coefficient fields admitted by QPA's transitive GAP algorithms.
///
/// `qpa_elements` follows GAP's canonical `Elements(GF(q))` order: zero,
/// followed by successive powers of the least primitive root. It returns
/// `None` for infinite fields.
pub trait QpaCoefficientField: Field {
    /// Return the field characteristic used by GAP's radical algorithm.
    fn qpa_characteristic(&self) -> QpaFieldCharacteristic;

    /// Return the field order, or `None` for an infinite coefficient field.
    fn qpa_order(&self) -> Option<Integer> {
        self.size()
    }

    /// Return the degree over the prime field, or `None` in characteristic zero.
    fn qpa_degree_over_prime_field(&self) -> Option<usize> {
        match self.qpa_characteristic() {
            QpaFieldCharacteristic::Zero => None,
            QpaFieldCharacteristic::Prime(_) => Some(1),
        }
    }

    /// Apply a possibly negative power of the prime-field Frobenius.
    ///
    /// Returns `None` in characteristic zero. The default is the identity on a
    /// prime field; nontrivial extensions override it.
    fn qpa_frobenius(&self, element: &Self::Element, _power: i64) -> Option<Self::Element> {
        match self.qpa_characteristic() {
            QpaFieldCharacteristic::Zero => None,
            QpaFieldCharacteristic::Prime(_) => Some(element.clone()),
        }
    }

    /// Return every element in GAP order, or `None` for an infinite field.
    ///
    /// # Errors
    ///
    /// Rejects a nonprime Numerica modular domain, an unaddressable field, or
    /// a missing multiplicative generator.
    fn qpa_elements(&self) -> Result<Option<Vec<Self::Element>>, CoefficientDomainError>;

    /// Compare two field elements using GAP's object order.
    ///
    /// # Errors
    ///
    /// Propagates finite-field enumeration failures needed to recover GAP's
    /// canonical finite-field element order.
    fn qpa_compare(
        &self,
        left: &Self::Element,
        right: &Self::Element,
    ) -> Result<Ordering, CoefficientDomainError>;
}

impl QpaCoefficientField for RationalField {
    fn qpa_characteristic(&self) -> QpaFieldCharacteristic {
        QpaFieldCharacteristic::Zero
    }

    fn qpa_elements(&self) -> Result<Option<Vec<Self::Element>>, CoefficientDomainError> {
        Ok(None)
    }

    fn qpa_compare(
        &self,
        left: &Self::Element,
        right: &Self::Element,
    ) -> Result<Ordering, CoefficientDomainError> {
        let left_scaled = left.numerator_ref().clone() * right.denominator_ref();
        let right_scaled = right.numerator_ref().clone() * left.denominator_ref();
        Ok(left_scaled.cmp(&right_scaled))
    }
}

impl QpaCoefficientField for Z2 {
    fn qpa_characteristic(&self) -> QpaFieldCharacteristic {
        QpaFieldCharacteristic::Prime(2)
    }

    fn qpa_elements(&self) -> Result<Option<Vec<Self::Element>>, CoefficientDomainError> {
        Ok(Some(vec![self.zero(), self.one()]))
    }

    fn qpa_compare(
        &self,
        left: &Self::Element,
        right: &Self::Element,
    ) -> Result<Ordering, CoefficientDomainError> {
        Ok(match (self.is_zero(left), self.is_zero(right)) {
            (true, true) | (false, false) => Ordering::Equal,
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
        })
    }
}

impl QpaCoefficientField for Zp {
    fn qpa_characteristic(&self) -> QpaFieldCharacteristic {
        QpaFieldCharacteristic::Prime(u64::from(self.get_prime()))
    }

    fn qpa_elements(&self) -> Result<Option<Vec<Self::Element>>, CoefficientDomainError> {
        let prime = self.get_prime();
        if prime < 3 || prime.is_multiple_of(2) || !self.is_prime_field(0) {
            return Err(CoefficientDomainError::UnsupportedPrime {
                characteristic: prime.to_string().into(),
            });
        }
        let size = u64::from(prime);
        let capacity =
            usize::try_from(size).map_err(|_| CoefficientDomainError::FieldTooLarge { size })?;
        let primitive_root = PrimitiveRootIterator::new(size)
            .next()
            .ok_or(CoefficientDomainError::MissingPrimitiveRoot { size })?;
        let primitive_root = u32::try_from(primitive_root)
            .map_err(|_| CoefficientDomainError::MissingPrimitiveRoot { size })?;
        let generator = self.to_element(primitive_root);
        let mut elements = Vec::with_capacity(capacity);
        elements.push(self.zero());
        let mut power = self.one();
        for _ in 0..size - 1 {
            elements.push(power);
            power = self.mul(&power, &generator);
        }
        Ok(Some(elements))
    }

    fn qpa_compare(
        &self,
        left: &Self::Element,
        right: &Self::Element,
    ) -> Result<Ordering, CoefficientDomainError> {
        let elements = self.qpa_elements()?.expect("Zp is finite");
        let left_position = elements.iter().position(|element| element == left);
        let right_position = elements.iter().position(|element| element == right);
        Ok(left_position.cmp(&right_position))
    }
}

impl QpaCoefficientField for GaloisField {
    fn qpa_characteristic(&self) -> QpaFieldCharacteristic {
        QpaFieldCharacteristic::Prime(u64::from(self.prime()))
    }

    fn qpa_degree_over_prime_field(&self) -> Option<usize> {
        Some(self.degree())
    }

    fn qpa_frobenius(&self, element: &Self::Element, power: i64) -> Option<Self::Element> {
        let degree = i64::try_from(self.degree()).ok()?;
        Some(self.frobenius(element, usize::try_from(power.rem_euclid(degree)).ok()?))
    }

    fn qpa_elements(&self) -> Result<Option<Vec<Self::Element>>, CoefficientDomainError> {
        self.elements_gap_order().map(Some).ok_or_else(|| {
            CoefficientDomainError::ExtensionFieldTooLarge {
                size: self.order().to_string().into(),
            }
        })
    }

    fn qpa_compare(
        &self,
        left: &Self::Element,
        right: &Self::Element,
    ) -> Result<Ordering, CoefficientDomainError> {
        let elements = self.qpa_elements()?.expect("GaloisField is finite");
        let left_position = elements.iter().position(|element| element == left);
        let right_position = elements.iter().position(|element| element == right);
        Ok(left_position.cmp(&right_position))
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use numerica::domains::{
        Field, Ring,
        finite_field::{FiniteFieldCore, Z2, Zp},
        rational::Q,
    };

    use crate::finite_field::GaloisField;

    use super::{CoefficientDomain, QpaCoefficientField};

    #[test]
    fn coefficient_descriptors_have_stable_json() {
        assert_eq!(
            serde_json::to_string(&CoefficientDomain::rational()).unwrap(),
            r#"{"kind":"rational"}"#
        );
        assert_eq!(
            serde_json::to_string(&CoefficientDomain::prime_field_unchecked("2")).unwrap(),
            r#"{"kind":"prime_field","characteristic":"2"}"#
        );
        let f4 = GaloisField::conway(2, 2).unwrap();
        assert_eq!(
            serde_json::to_string(&CoefficientDomain::finite_extension_field(&f4)).unwrap(),
            r#"{"kind":"finite_extension_field","characteristic":"2","degree":2,"modulus":["1","1","1"],"gap_conway":true}"#
        );
    }

    #[test]
    fn prime_field_elements_follow_gap_primitive_root_order() {
        assert_eq!(Z2.qpa_elements().unwrap().unwrap().len(), 2);
        assert!(Q.qpa_elements().unwrap().is_none());

        let field5 = Zp::new(5);
        let elements5 = field5.qpa_elements().unwrap().unwrap();
        assert_eq!(
            elements5
                .iter()
                .map(|element| field5.from_element(element))
                .collect::<Vec<_>>(),
            [0, 1, 2, 4, 3]
        );
        assert_eq!(
            field5
                .qpa_compare(&field5.to_element(4), &field5.to_element(3))
                .unwrap(),
            Ordering::Less
        );
        let half = Q.div(&Q.one(), &Q.nth(2.into()));
        let two_thirds = Q.div(&Q.nth(2.into()), &Q.nth(3.into()));
        assert_eq!(Q.qpa_compare(&half, &two_thirds).unwrap(), Ordering::Less);
        let field7 = Zp::new(7);
        assert_eq!(
            field7
                .qpa_elements()
                .unwrap()
                .unwrap()
                .iter()
                .map(|element| field7.from_element(element))
                .collect::<Vec<_>>(),
            [0, 1, 3, 2, 6, 4, 5]
        );
    }

    #[test]
    fn extension_metadata_and_frobenius_are_not_specialized_to_prime_fields() {
        let field = GaloisField::conway(2, 3).unwrap();
        let generator = field.primitive_element();
        assert_eq!(field.qpa_order(), Some(8.into()));
        assert_eq!(field.qpa_degree_over_prime_field(), Some(3));
        assert_eq!(
            field.qpa_frobenius(&generator, -1).unwrap(),
            field.frobenius(&generator, 2)
        );
        assert_eq!(field.qpa_elements().unwrap().unwrap().len(), 8);

        let unaddressable = GaloisField::conway(2, 120).unwrap();
        assert!(matches!(
            unaddressable.qpa_elements(),
            Err(super::CoefficientDomainError::ExtensionFieldTooLarge { .. })
        ));
    }
}
