//! Python power syntax backed by the existing exact element operations.

use pyo3::exceptions::PyTypeError;

use super::tensor::QuiverAlgebraElementKind;
#[allow(clippy::wildcard_imports)]
// Element runtime variants belong to the parent binding domain.
use super::*;

fn nonnegative_exponent(
    exponent: &Bound<'_, PyAny>,
    modulo: Option<&Bound<'_, PyAny>>,
) -> PyResult<u64> {
    if modulo.is_some() {
        return Err(PyTypeError::new_err(
            "modular powers are not defined for algebra elements",
        ));
    }
    // Match the existing pow(u64) conversion, including its rejection of
    // negative integers, overflow, and nonintegral exponents.
    exponent.extract::<u64>()
}

#[pymethods]
impl PyAlgebraElement {
    /// Raise to a nonnegative integer power; modular powers are unsupported.
    fn __pow__(
        &self,
        exponent: &Bound<'_, PyAny>,
        modulo: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        self.pow(nonnegative_exponent(exponent, modulo)?)
    }
}

#[pymethods]
impl PyQuotientElement {
    /// Raise to a nonnegative integer power; use `pow_integer` for inverses.
    fn __pow__(
        &self,
        exponent: &Bound<'_, PyAny>,
        modulo: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let exponent = nonnegative_exponent(exponent, modulo)?;
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| element.pow(exponent).map_err(invalid)?
        );
        Ok(Self { kind })
    }
}

#[pymethods]
impl PyQuiverAlgebraElement {
    /// Raise to a nonnegative integer power; modular powers are unsupported.
    fn __pow__(
        &self,
        exponent: &Bound<'_, PyAny>,
        modulo: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let exponent = nonnegative_exponent(exponent, modulo)?;
        let kind = map_ref!(
            &self.kind,
            QuiverAlgebraElementKind => QuiverAlgebraElementKind,
            |element| match element {
                QuiverAlgebraElement::Path(value) =>
                    QuiverAlgebraElement::Path(value.pow(exponent).map_err(invalid)?),
                QuiverAlgebraElement::Quotient(value) =>
                    QuiverAlgebraElement::Quotient(value.pow(exponent).map_err(invalid)?),
            }
        );
        Ok(Self { kind })
    }

    fn __neg__(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuiverAlgebraElementKind => QuiverAlgebraElementKind,
            |element| match element {
                QuiverAlgebraElement::Path(value) =>
                    QuiverAlgebraElement::Path(value.neg().map_err(invalid)?),
                QuiverAlgebraElement::Quotient(value) =>
                    QuiverAlgebraElement::Quotient(value.neg().map_err(invalid)?),
            }
        );
        Ok(Self { kind })
    }
}
