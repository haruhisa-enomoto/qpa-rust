//! Runtime descriptions of the initially supported exact coefficient fields.

use numerica::domains::{
    Field, Ring, RingOps,
    finite_field::{FiniteFieldCore, PrimitiveRootIterator, Z2, Zp},
    rational::RationalField,
};
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// An exact coefficient field admitted by the initial Numerica-backed port.
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
/// `qpa_elements` follows GAP's canonical `Elements(GF(p))` order: zero,
/// followed by successive powers of the least primitive root. It returns
/// `None` for infinite fields.
pub trait QpaCoefficientField: Field {
    /// Return the field characteristic used by GAP's radical algorithm.
    fn qpa_characteristic(&self) -> QpaFieldCharacteristic;

    /// Return every element in GAP order, or `None` for an infinite field.
    ///
    /// # Errors
    ///
    /// Rejects a nonprime Numerica modular domain, an unaddressable field, or
    /// a missing multiplicative generator.
    fn qpa_elements(&self) -> Result<Option<Vec<Self::Element>>, CoefficientDomainError>;
}

impl QpaCoefficientField for RationalField {
    fn qpa_characteristic(&self) -> QpaFieldCharacteristic {
        QpaFieldCharacteristic::Zero
    }

    fn qpa_elements(&self) -> Result<Option<Vec<Self::Element>>, CoefficientDomainError> {
        Ok(None)
    }
}

impl QpaCoefficientField for Z2 {
    fn qpa_characteristic(&self) -> QpaFieldCharacteristic {
        QpaFieldCharacteristic::Prime(2)
    }

    fn qpa_elements(&self) -> Result<Option<Vec<Self::Element>>, CoefficientDomainError> {
        Ok(Some(vec![self.zero(), self.one()]))
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
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        finite_field::{FiniteFieldCore, Z2, Zp},
        rational::Q,
    };

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
    }

    #[test]
    fn prime_field_elements_follow_gap_primitive_root_order() {
        assert_eq!(Z2.qpa_elements().unwrap().unwrap().len(), 2);
        assert!(Q.qpa_elements().unwrap().is_none());

        let field5 = Zp::new(5);
        assert_eq!(
            field5
                .qpa_elements()
                .unwrap()
                .unwrap()
                .iter()
                .map(|element| field5.from_element(element))
                .collect::<Vec<_>>(),
            [0, 1, 2, 4, 3]
        );
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
}
