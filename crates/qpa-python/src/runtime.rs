//! Shared runtime dispatch and exact coefficient conversion for Python wrappers.

use numerica::domains::{
    Ring as _,
    finite_field::{FiniteFieldCore, ToFiniteField, Z2, Zp},
    integer::Integer,
    rational::RationalField,
};
use pyo3::PyResult;
use qpa_core::GaloisField;
use qpa_core::QpaCoefficientField;

use crate::{FieldKind, PyScalar, ScalarKind, invalid};

/// Extract one monomorphized Rust value from a runtime wrapper enum.
pub(crate) trait Extract<T> {
    /// Return the matching variant or a Python family error.
    fn extract(self) -> PyResult<T>;
}

/// Extract a homogeneous Python sequence into one monomorphized Rust sequence.
pub(crate) fn extract_all<S, T>(values: Vec<S>) -> PyResult<Vec<T>>
where
    S: Extract<T>,
{
    values.into_iter().map(Extract::extract).collect()
}

/// Wrap one monomorphized Rust value in its runtime-dispatched enum.
pub(crate) trait Wrap<T> {
    /// Preserve the concrete coefficient field in the runtime wrapper.
    fn wrap(value: T) -> Self;
}

/// Wrap a homogeneous Rust sequence while preserving its coefficient field.
pub(crate) fn wrap_all<S, T>(values: Vec<T>) -> Vec<S>
where
    S: Wrap<T>,
{
    values.into_iter().map(S::wrap).collect()
}

/// Exact field behavior required by every runtime-dispatched Python wrapper.
pub(crate) trait PythonField: QpaCoefficientField {
    /// Return the corresponding runtime descriptor.
    fn python_kind(&self) -> FieldKind;

    /// Convert one Python scalar after enforcing coefficient-family identity.
    fn scalar_value(&self, scalar: &PyScalar) -> PyResult<Self::Element>;

    /// Convert an arbitrary-precision integer into this exact field.
    fn integer_value(&self, value: &Integer) -> Self::Element;

    /// Convert an exact field element back to its Python representation.
    fn wrap_scalar(&self, value: &Self::Element) -> PyScalar;

    /// Format a field element using the stable Python representative.
    fn scalar_text(&self, value: &Self::Element) -> String;
}

impl PythonField for RationalField {
    fn python_kind(&self) -> FieldKind {
        FieldKind::Rational
    }

    fn scalar_value(&self, scalar: &PyScalar) -> PyResult<Self::Element> {
        match &scalar.value {
            ScalarKind::Rational(value) => Ok(value.clone()),
            ScalarKind::Prime { .. } => Err(invalid("coefficient belongs to a prime field")),
            ScalarKind::Extension { .. } => {
                Err(invalid("coefficient belongs to a finite extension field"))
            }
        }
    }

    fn integer_value(&self, value: &Integer) -> Self::Element {
        self.nth(value.clone())
    }

    fn wrap_scalar(&self, value: &Self::Element) -> PyScalar {
        PyScalar {
            value: ScalarKind::Rational(value.clone()),
        }
    }

    fn scalar_text(&self, value: &Self::Element) -> String {
        value.to_string()
    }
}

impl PythonField for Z2 {
    fn python_kind(&self) -> FieldKind {
        FieldKind::Binary
    }

    fn scalar_value(&self, scalar: &PyScalar) -> PyResult<Self::Element> {
        match scalar.value {
            ScalarKind::Prime {
                characteristic: 2,
                representative,
            } => Ok(representative.to_finite_field(self)),
            _ => Err(invalid("coefficient does not belong to GF(2)")),
        }
    }

    fn integer_value(&self, value: &Integer) -> Self::Element {
        value.to_finite_field(self)
    }

    fn wrap_scalar(&self, value: &Self::Element) -> PyScalar {
        PyScalar::prime(2, u32::from(!self.is_zero(value)))
    }

    fn scalar_text(&self, value: &Self::Element) -> String {
        u32::from(!self.is_zero(value)).to_string()
    }
}

impl PythonField for Zp {
    fn python_kind(&self) -> FieldKind {
        FieldKind::Prime(self.get_prime())
    }

    fn scalar_value(&self, scalar: &PyScalar) -> PyResult<Self::Element> {
        match scalar.value {
            ScalarKind::Prime {
                characteristic,
                representative,
            } if characteristic == self.get_prime() => Ok(self.to_element(representative)),
            _ => Err(invalid(format!(
                "coefficient does not belong to GF({})",
                self.get_prime()
            ))),
        }
    }

    fn integer_value(&self, value: &Integer) -> Self::Element {
        value.to_finite_field(self)
    }

    fn wrap_scalar(&self, value: &Self::Element) -> PyScalar {
        PyScalar::prime(self.get_prime(), self.from_element(value))
    }

    fn scalar_text(&self, value: &Self::Element) -> String {
        self.from_element(value).to_string()
    }
}

impl PythonField for GaloisField {
    fn python_kind(&self) -> FieldKind {
        FieldKind::Extension {
            characteristic: self.prime(),
            degree: self.degree(),
        }
    }

    fn scalar_value(&self, scalar: &PyScalar) -> PyResult<Self::Element> {
        match &scalar.value {
            ScalarKind::Extension {
                characteristic,
                degree,
                coefficients,
            } if *characteristic == self.prime() && *degree == self.degree() => self
                .element(coefficients.clone())
                .map_err(|error| invalid(error.to_string())),
            _ => Err(invalid(format!(
                "coefficient does not belong to GF({}^{})",
                self.prime(),
                self.degree()
            ))),
        }
    }

    fn integer_value(&self, value: &Integer) -> Self::Element {
        self.nth(value.clone())
    }

    fn wrap_scalar(&self, value: &Self::Element) -> PyScalar {
        PyScalar::extension(self.prime(), self.degree(), value.coefficients().to_vec())
    }

    fn scalar_text(&self, value: &Self::Element) -> String {
        self.printer(value).to_string()
    }
}

/// Dispatch the same borrowed operation over the three supported field variants.
macro_rules! dispatch_ref {
    ($value:expr, $kind:ident, |$inner:ident| $body:expr) => {
        match $value {
            $kind::Rational($inner) => $body,
            $kind::Binary($inner) => $body,
            $kind::Prime($inner) => $body,
            $kind::Extension($inner) => $body,
        }
    };
}

/// Map one runtime wrapper enum into another while preserving its field variant.
macro_rules! map_ref {
    ($value:expr, $source:ident => $target:ident, |$inner:ident| $body:expr) => {
        match $value {
            $source::Rational($inner) => $target::Rational($body),
            $source::Binary($inner) => $target::Binary($body),
            $source::Prime($inner) => $target::Prime($body),
            $source::Extension($inner) => $target::Extension($body),
        }
    };
}

/// Dispatch two wrappers only when their runtime field variants agree.
macro_rules! dispatch_pair_ref {
    ($left:expr, $right:expr, $kind:ident, |$lhs:ident, $rhs:ident| $body:expr, $mismatch:expr) => {
        match ($left, $right) {
            ($kind::Rational($lhs), $kind::Rational($rhs)) => $body,
            ($kind::Binary($lhs), $kind::Binary($rhs)) => $body,
            ($kind::Prime($lhs), $kind::Prime($rhs)) => $body,
            ($kind::Extension($lhs), $kind::Extension($rhs)) => $body,
            _ => $mismatch,
        }
    };
}

/// Map two same-field wrappers into another runtime wrapper or return a family error.
macro_rules! map_pair_ref {
    (
        $left:expr,
        $right:expr,
        $source:ident => $target:ident,
        |$lhs:ident, $rhs:ident| $body:expr,
        $mismatch:expr
    ) => {
        match ($left, $right) {
            ($source::Rational($lhs), $source::Rational($rhs)) => Ok($target::Rational($body)),
            ($source::Binary($lhs), $source::Binary($rhs)) => Ok($target::Binary($body)),
            ($source::Prime($lhs), $source::Prime($rhs)) => Ok($target::Prime($body)),
            ($source::Extension($lhs), $source::Extension($rhs)) => Ok($target::Extension($body)),
            _ => Err($mismatch),
        }
    };
}

pub(crate) use {dispatch_pair_ref, dispatch_ref, map_pair_ref, map_ref};

/// Construct a scalar from an exact numerator and denominator in one field.
pub(crate) fn scalar_fraction<F: PythonField>(
    field: &F,
    numerator: &Integer,
    denominator: &Integer,
) -> PyResult<PyScalar> {
    let numerator = field.integer_value(numerator);
    let denominator = field.integer_value(denominator);
    if field.is_zero(&denominator) {
        return Err(invalid(format!(
            "coefficient denominator is zero in {}",
            field.python_kind().display()
        )));
    }
    Ok(field.wrap_scalar(&field.div(&numerator, &denominator)))
}
