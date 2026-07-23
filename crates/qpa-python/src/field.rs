//! Python descriptors and exact conversions for Numerica coefficient fields.

use std::{
    collections::hash_map::DefaultHasher,
    hash::{Hash, Hasher},
    num::NonZeroUsize,
    str::FromStr,
};

use numerica::domains::{
    Ring,
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
use qpa_core::GaloisField;

use crate::{InvalidInputError, runtime, runtime::PythonField as _};

fn finish_python_hash(hasher: DefaultHasher) -> isize {
    let maximum = u64::try_from(isize::MAX).expect("isize::MAX always fits u64");
    isize::try_from(hasher.finish() % maximum)
        .expect("a nonnegative residue below isize::MAX fits isize")
}

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

pub(crate) fn parse_python_integer_vector(
    py: Python<'_>,
    values: Vec<Py<PyAny>>,
) -> PyResult<Vec<Integer>> {
    values
        .into_iter()
        .map(|value| parse_python_integer(value.bind(py)))
        .collect()
}

pub(crate) fn parse_python_integer_matrix(
    py: Python<'_>,
    rows: Vec<Vec<Py<PyAny>>>,
) -> PyResult<Vec<Vec<Integer>>> {
    rows.into_iter()
        .map(|row| parse_python_integer_vector(py, row))
        .collect()
}

pub(crate) fn wrap_python_integer_vector(
    py: Python<'_>,
    values: &[Integer],
) -> PyResult<Vec<Py<PyAny>>> {
    values
        .iter()
        .map(|value| wrap_python_integer(py, value))
        .collect()
}

pub(crate) fn wrap_python_integer_matrix(
    py: Python<'_>,
    rows: impl Iterator<Item = Vec<Integer>>,
) -> PyResult<Vec<Vec<Py<PyAny>>>> {
    rows.map(|row| wrap_python_integer_vector(py, &row))
        .collect()
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum FieldKind {
    Rational,
    Binary,
    Prime(u32),
    Extension { characteristic: u32, degree: usize },
}

impl FieldKind {
    pub(crate) fn display(self) -> String {
        match self {
            Self::Rational => "Rationals".into(),
            Self::Binary => "GF(2)".into(),
            Self::Prime(characteristic) => format!("GF({characteristic})"),
            Self::Extension {
                characteristic,
                degree,
            } => format!("GF({characteristic}^{degree})"),
        }
    }

    fn characteristic(self) -> Option<u32> {
        match self {
            Self::Rational => None,
            Self::Binary => Some(2),
            Self::Prime(characteristic) | Self::Extension { characteristic, .. } => {
                Some(characteristic)
            }
        }
    }

    pub(crate) fn galois_field(self) -> Option<GaloisField> {
        match self {
            Self::Extension {
                characteristic,
                degree,
            } => Some(
                GaloisField::conway(characteristic, degree)
                    .expect("validated Python extension-field descriptor"),
            ),
            _ => None,
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
    /// Return the exact rational field.
    #[staticmethod]
    fn rational() -> Self {
        Self {
            kind: FieldKind::Rational,
        }
    }

    /// Return the prime field `F_p` after verifying that `p` is prime.
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

    /// Return GAP's Conway model of the finite field `GF(p^n)`.
    #[staticmethod]
    fn extension(characteristic: u32, degree: usize) -> PyResult<Self> {
        GaloisField::conway(characteristic, degree)
            .map_err(|error| InvalidInputError::new_err(error.to_string()))?;
        Ok(Self {
            kind: FieldKind::Extension {
                characteristic,
                degree,
            },
        })
    }

    #[getter]
    fn characteristic(&self) -> Option<u32> {
        self.kind.characteristic()
    }

    /// Return the degree over the prime field, or `None` for the rationals.
    #[getter]
    fn degree(&self) -> Option<usize> {
        match self.kind {
            FieldKind::Rational => None,
            FieldKind::Binary | FieldKind::Prime(_) => Some(1),
            FieldKind::Extension { degree, .. } => Some(degree),
        }
    }

    /// Return the Conway polynomial coefficients, constant coefficient first.
    #[getter]
    fn conway_polynomial(&self) -> Option<Vec<u32>> {
        self.kind
            .galois_field()
            .map(|field| field.modulus().to_vec())
    }

    #[getter]
    fn is_finite(&self) -> bool {
        !matches!(self.kind, FieldKind::Rational)
    }

    /// Construct an exact field element from an integer numerator and denominator.
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
            FieldKind::Extension { .. } => runtime::scalar_fraction(
                &self
                    .kind
                    .galois_field()
                    .expect("extension variant has a field"),
                &numerator,
                &denominator,
            ),
        }
    }

    /// Construct an extension-field element from polynomial-basis coefficients.
    fn element_from_coefficients(&self, coefficients: Vec<u32>) -> PyResult<PyScalar> {
        let field = self.kind.galois_field().ok_or_else(|| {
            InvalidInputError::new_err(
                "polynomial-basis coefficients require a finite extension field",
            )
        })?;
        let element = field
            .element(coefficients)
            .map_err(|error| InvalidInputError::new_err(error.to_string()))?;
        Ok(field.wrap_scalar(&element))
    }

    /// Return GAP's designated primitive generator of an extension field.
    fn primitive_element(&self) -> PyResult<PyScalar> {
        let field = self.kind.galois_field().ok_or_else(|| {
            InvalidInputError::new_err("a primitive extension generator requires GF(p^n), n > 1")
        })?;
        Ok(field.wrap_scalar(&field.primitive_element()))
    }

    /// Apply a possibly negative power of the prime-field Frobenius.
    fn frobenius(&self, scalar: &PyScalar, power: i64) -> PyResult<PyScalar> {
        let field = self.kind.galois_field().ok_or_else(|| {
            InvalidInputError::new_err("Frobenius inspection requires a finite extension field")
        })?;
        let value = field.scalar_value(scalar)?;
        let degree = i64::try_from(field.degree()).expect("validated extension degree fits i64");
        Ok(field.wrap_scalar(&field.frobenius(
            &value,
            usize::try_from(power.rem_euclid(degree)).expect("nonnegative"),
        )))
    }

    /// Embed a scalar from a canonical Conway subfield into this field.
    fn embed(&self, scalar: &PyScalar) -> PyResult<PyScalar> {
        let target = self.kind.galois_field().ok_or_else(|| {
            InvalidInputError::new_err("a canonical embedding target must be GF(p^n), n > 1")
        })?;
        match &scalar.value {
            ScalarKind::Prime {
                characteristic,
                representative,
            } if *characteristic == target.prime() => {
                Ok(target.wrap_scalar(&target.nth(Integer::from(*representative))))
            }
            ScalarKind::Extension {
                characteristic,
                degree,
                ..
            } => {
                let source = GaloisField::conway(*characteristic, *degree)
                    .map_err(|error| InvalidInputError::new_err(error.to_string()))?;
                let value = source.scalar_value(scalar)?;
                let embedded = target
                    .embed_subfield_element(&source, &value)
                    .map_err(|error| InvalidInputError::new_err(error.to_string()))?;
                Ok(target.wrap_scalar(&embedded))
            }
            _ => Err(InvalidInputError::new_err(
                "the scalar field is not a canonical subfield of the target",
            )),
        }
    }

    /// Translate QPA's `NullList(length, field)` operation.
    fn null_list(&self, length: usize) -> PyResult<Vec<PyScalar>> {
        let length = NonZeroUsize::new(length)
            .ok_or_else(|| InvalidInputError::new_err("NullList length must be positive"))?;
        Ok(match self.kind {
            FieldKind::Rational => qpa_core::null_list(length, &Q)
                .iter()
                .map(|value| Q.wrap_scalar(value))
                .collect(),
            FieldKind::Binary => qpa_core::null_list(length, &Z2)
                .iter()
                .map(|value| Z2.wrap_scalar(value))
                .collect(),
            FieldKind::Prime(characteristic) => {
                let field = Zp::new(characteristic);
                qpa_core::null_list(length, &field)
                    .iter()
                    .map(|value| field.wrap_scalar(value))
                    .collect()
            }
            FieldKind::Extension { .. } => {
                let field = self
                    .kind
                    .galois_field()
                    .expect("extension variant has a field");
                qpa_core::null_list(length, &field)
                    .iter()
                    .map(|value| field.wrap_scalar(value))
                    .collect()
            }
        })
    }

    fn __repr__(&self) -> String {
        self.kind.display()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }

    fn __hash__(&self) -> isize {
        let mut hasher = DefaultHasher::new();
        self.kind.hash(&mut hasher);
        finish_python_hash(hasher)
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum ScalarKind {
    Rational(Rational),
    Prime {
        characteristic: u32,
        representative: u32,
    },
    Extension {
        characteristic: u32,
        degree: usize,
        coefficients: Vec<u32>,
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

    pub(crate) fn extension(characteristic: u32, degree: usize, coefficients: Vec<u32>) -> Self {
        Self {
            value: ScalarKind::Extension {
                characteristic,
                degree,
                coefficients,
            },
        }
    }

    pub(crate) fn value_for<F: runtime::PythonField>(&self, field: &F) -> PyResult<F::Element> {
        field.scalar_value(self)
    }

    pub(crate) fn is_zero_value(&self) -> bool {
        match &self.value {
            ScalarKind::Rational(value) => value.is_zero(),
            ScalarKind::Prime { representative, .. } => *representative == 0,
            ScalarKind::Extension { coefficients, .. } => {
                coefficients.iter().all(|&coefficient| coefficient == 0)
            }
        }
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
            ScalarKind::Extension {
                characteristic,
                degree,
                ..
            } => FieldKind::Extension {
                characteristic,
                degree,
            },
        };
        PyField { kind }
    }

    #[getter]
    fn numerator(&self) -> Option<String> {
        match &self.value {
            ScalarKind::Rational(value) => Some(value.numerator_ref().to_string()),
            ScalarKind::Prime { .. } | ScalarKind::Extension { .. } => None,
        }
    }

    #[getter]
    fn denominator(&self) -> Option<String> {
        match &self.value {
            ScalarKind::Rational(value) => Some(value.denominator_ref().to_string()),
            ScalarKind::Prime { .. } | ScalarKind::Extension { .. } => None,
        }
    }

    #[getter]
    fn representative(&self) -> Option<u32> {
        match self.value {
            ScalarKind::Prime { representative, .. } => Some(representative),
            ScalarKind::Rational(_) | ScalarKind::Extension { .. } => None,
        }
    }

    /// Return polynomial-basis coefficients for an extension element.
    #[getter]
    fn coefficients(&self) -> Option<Vec<u32>> {
        match &self.value {
            ScalarKind::Extension { coefficients, .. } => Some(coefficients.clone()),
            _ => None,
        }
    }

    #[getter]
    fn is_zero(&self) -> bool {
        self.is_zero_value()
    }

    fn __repr__(&self) -> String {
        match &self.value {
            ScalarKind::Rational(value) => value.to_string(),
            ScalarKind::Prime {
                characteristic,
                representative,
            } => format!("{representative} mod {characteristic}"),
            ScalarKind::Extension {
                characteristic,
                degree,
                coefficients,
            } => GaloisField::conway(*characteristic, *degree)
                .and_then(|field| {
                    field
                        .element(coefficients.clone())
                        .map(|value| (field, value))
                })
                .map_or_else(
                    |_| format!("{coefficients:?} in GF({characteristic}^{degree})"),
                    |(field, value)| field.scalar_text(&value),
                ),
        }
    }

    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }

    fn __hash__(&self) -> isize {
        let mut hasher = DefaultHasher::new();
        match &self.value {
            ScalarKind::Rational(value) => {
                value.numerator_ref().to_string().hash(&mut hasher);
                value.denominator_ref().to_string().hash(&mut hasher);
            }
            ScalarKind::Prime {
                characteristic,
                representative,
            } => {
                characteristic.hash(&mut hasher);
                representative.hash(&mut hasher);
            }
            ScalarKind::Extension {
                characteristic,
                degree,
                coefficients,
            } => {
                characteristic.hash(&mut hasher);
                degree.hash(&mut hasher);
                coefficients.hash(&mut hasher);
            }
        }
        finish_python_hash(hasher)
    }
}
