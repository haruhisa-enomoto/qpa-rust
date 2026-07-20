//! Python descriptors and exact conversions for Numerica coefficient fields.

use std::str::FromStr;

use numerica::domains::{
    finite_field::{PrimeIteratorU64, Z2, Zp},
    integer::Integer,
    rational::{Q, Rational},
};
use pyo3::{
    Bound, PyResult,
    exceptions::PyTypeError,
    prelude::*,
    types::{PyAny, PyInt},
};

use crate::{InvalidInputError, runtime};

pub(crate) fn parse_python_integer(value: &Bound<'_, PyAny>) -> PyResult<Integer> {
    if !value.is_instance_of::<PyInt>() {
        return Err(PyTypeError::new_err("an exact Python integer is required"));
    }
    Integer::from_str(value.str()?.to_str()?)
        .map_err(|error| PyTypeError::new_err(format!("invalid Python integer: {error}")))
}

/// Convert a Numerica arbitrary-precision integer to Python's exact `int`.
pub(crate) fn wrap_python_integer(py: Python<'_>, value: &Integer) -> PyResult<Py<PyAny>> {
    Ok(py
        .import("builtins")?
        .getattr("int")?
        .call1((value.to_string(),))?
        .unbind())
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum FieldKind {
    Rational,
    Binary,
    Prime(u32),
}

impl FieldKind {
    pub(crate) fn display(self) -> String {
        match self {
            Self::Rational => "Rationals".into(),
            Self::Binary => "GF(2)".into(),
            Self::Prime(characteristic) => format!("GF({characteristic})"),
        }
    }

    fn characteristic(self) -> Option<u32> {
        match self {
            Self::Rational => None,
            Self::Binary => Some(2),
            Self::Prime(characteristic) => Some(characteristic),
        }
    }
}

/// Runtime exact field descriptor used by the Python binding.
#[pyclass(name = "Field", frozen, skip_from_py_object)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct PyField {
    pub(crate) kind: FieldKind,
}

#[pymethods]
impl PyField {
    #[staticmethod]
    fn rational() -> Self {
        Self {
            kind: FieldKind::Rational,
        }
    }

    #[staticmethod]
    fn prime(characteristic: u32) -> PyResult<Self> {
        let kind = match characteristic {
            2 => FieldKind::Binary,
            3.. => {
                if PrimeIteratorU64::new(u64::from(characteristic) - 1).next()
                    != Some(u64::from(characteristic))
                {
                    return Err(InvalidInputError::new_err(format!(
                        "{characteristic} is not prime"
                    )));
                }
                FieldKind::Prime(characteristic)
            }
            _ => {
                return Err(InvalidInputError::new_err(format!(
                    "{characteristic} is not prime"
                )));
            }
        };
        Ok(Self { kind })
    }

    #[getter]
    fn characteristic(&self) -> Option<u32> {
        self.kind.characteristic()
    }

    #[getter]
    fn is_finite(&self) -> bool {
        !matches!(self.kind, FieldKind::Rational)
    }

    #[pyo3(signature = (numerator, denominator=None))]
    fn element(
        &self,
        numerator: &Bound<'_, PyAny>,
        denominator: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyScalar> {
        let numerator = parse_python_integer(numerator)?;
        let denominator = denominator
            .map(parse_python_integer)
            .transpose()?
            .unwrap_or_else(Integer::one);
        if denominator.is_zero() {
            return Err(InvalidInputError::new_err(
                "coefficient denominator is zero",
            ));
        }
        match self.kind {
            FieldKind::Rational => runtime::scalar_fraction(&Q, &numerator, &denominator),
            FieldKind::Binary => runtime::scalar_fraction(&Z2, &numerator, &denominator),
            FieldKind::Prime(characteristic) => {
                runtime::scalar_fraction(&Zp::new(characteristic), &numerator, &denominator)
            }
        }
    }

    fn __repr__(&self) -> String {
        self.kind.display()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum ScalarKind {
    Rational(Rational),
    Prime {
        characteristic: u32,
        representative: u32,
    },
}

/// Exact scalar transferred between Python and a Numerica-backed Rust field.
#[pyclass(name = "Scalar", frozen, skip_from_py_object)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct PyScalar {
    pub(crate) value: ScalarKind,
}

impl PyScalar {
    pub(crate) fn prime(characteristic: u32, representative: u32) -> Self {
        Self {
            value: ScalarKind::Prime {
                characteristic,
                representative,
            },
        }
    }

    pub(crate) fn value_for<F: runtime::PythonField>(&self, field: &F) -> PyResult<F::Element> {
        field.scalar_value(self)
    }
}

#[pymethods]
impl PyScalar {
    #[getter]
    fn field(&self) -> PyField {
        let kind = match self.value {
            ScalarKind::Rational(_) => FieldKind::Rational,
            ScalarKind::Prime {
                characteristic: 2, ..
            } => FieldKind::Binary,
            ScalarKind::Prime { characteristic, .. } => FieldKind::Prime(characteristic),
        };
        PyField { kind }
    }

    #[getter]
    fn numerator(&self) -> Option<String> {
        match &self.value {
            ScalarKind::Rational(value) => Some(value.numerator_ref().to_string()),
            ScalarKind::Prime { .. } => None,
        }
    }

    #[getter]
    fn denominator(&self) -> Option<String> {
        match &self.value {
            ScalarKind::Rational(value) => Some(value.denominator_ref().to_string()),
            ScalarKind::Prime { .. } => None,
        }
    }

    #[getter]
    fn representative(&self) -> Option<u32> {
        match self.value {
            ScalarKind::Rational(_) => None,
            ScalarKind::Prime { representative, .. } => Some(representative),
        }
    }

    fn __repr__(&self) -> String {
        match &self.value {
            ScalarKind::Rational(value) => value.to_string(),
            ScalarKind::Prime {
                characteristic,
                representative,
            } => format!("{representative} mod {characteristic}"),
        }
    }

    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }
}
