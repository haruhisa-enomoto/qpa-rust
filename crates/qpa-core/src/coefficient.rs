//! Runtime descriptions of the initially supported exact coefficient fields.

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
}

#[cfg(test)]
mod tests {
    use super::CoefficientDomain;

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
}
