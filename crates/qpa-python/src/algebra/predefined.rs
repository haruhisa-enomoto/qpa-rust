//! Python dispatch for QPA's basic predefined algebra constructors.

use pyo3::{PyAny, pyfunction, pymethods};
use qpa_core::{PathAlgebra, QuiverAlgebra};
use qpa_gbnp::{NakayamaAlgebraResult, PredefinedAlgebra};

#[allow(clippy::wildcard_imports)]
// Runtime fields and the quiver-algebra wrapper belong to the parent module.
use super::*;
use crate::field::{
    parse_python_integer, parse_python_integer_vector, wrap_python_integer_matrix,
    wrap_python_integer_vector,
};

/// A QPA predefined algebra together with its constructor-installed filters.
#[pyclass(name = "PredefinedAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPredefinedAlgebra {
    algebra: PyQuiverAlgebra,
    kind: &'static str,
    is_nakayama: Option<bool>,
    is_canonical: Option<bool>,
    is_kronecker: Option<bool>,
    is_finite_type: Option<bool>,
    is_finite_global_dimension: Option<bool>,
    is_hereditary: Option<bool>,
}

#[pymethods]
impl PyPredefinedAlgebra {
    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        self.algebra.clone()
    }

    #[getter]
    fn kind(&self) -> &'static str {
        self.kind
    }

    #[getter]
    fn is_nakayama_algebra(&self) -> Option<bool> {
        self.is_nakayama
    }

    #[getter]
    fn is_canonical_algebra(&self) -> Option<bool> {
        self.is_canonical
    }

    #[getter]
    fn is_kronecker_algebra(&self) -> Option<bool> {
        self.is_kronecker
    }

    #[getter]
    fn is_finite_type_algebra(&self) -> Option<bool> {
        self.is_finite_type
    }

    #[getter]
    fn is_finite_global_dimension_algebra(&self) -> Option<bool> {
        self.is_finite_global_dimension
    }

    #[getter]
    fn is_hereditary_algebra(&self) -> Option<bool> {
        self.is_hereditary
    }

    fn __repr__(&self) -> String {
        format!("<QPA {} algebra>", self.kind)
    }
}

/// QPA's algebra-or-original-sequence result from `NakayamaAlgebra`.
#[pyclass(name = "NakayamaAlgebraResult", frozen, skip_from_py_object)]
pub(crate) struct PyNakayamaAlgebraResult {
    algebra: Option<PyPredefinedAlgebra>,
    invalid_sequence: Option<Vec<Integer>>,
    diagnostic: Option<&'static str>,
}

#[pymethods]
impl PyNakayamaAlgebraResult {
    #[getter]
    fn algebra(&self) -> Option<PyPredefinedAlgebra> {
        self.algebra.clone()
    }

    #[getter]
    fn invalid_sequence(&self, py: Python<'_>) -> PyResult<Option<Vec<Py<PyAny>>>> {
        self.invalid_sequence
            .as_deref()
            .map(|sequence| wrap_python_integer_vector(py, sequence))
            .transpose()
    }

    #[getter]
    fn diagnostic(&self) -> Option<&'static str> {
        self.diagnostic
    }

    #[getter]
    fn is_valid(&self) -> bool {
        self.algebra.is_some()
    }

    fn __repr__(&self) -> &'static str {
        if self.is_valid() {
            "<valid QPA Nakayama algebra result>"
        } else {
            "<invalid QPA Nakayama sequence result>"
        }
    }
}

fn wrap_predefined<F: runtime::PythonField>(algebra: PredefinedAlgebra<F>) -> PyPredefinedAlgebra
where
    super::conversion::QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    let kind = match algebra.kind() {
        qpa_gbnp::PredefinedAlgebraKind::Nakayama => "Nakayama",
        qpa_gbnp::PredefinedAlgebraKind::Canonical => "Canonical",
        qpa_gbnp::PredefinedAlgebraKind::Kronecker => "Kronecker",
    };
    PyPredefinedAlgebra {
        is_nakayama: algebra.is_nakayama_algebra(),
        is_canonical: algebra.is_canonical_algebra(),
        is_kronecker: algebra.is_kronecker_algebra(),
        is_finite_type: algebra.is_finite_type_algebra(),
        is_finite_global_dimension: algebra.is_finite_global_dimension_algebra(),
        is_hereditary: algebra.is_hereditary_algebra(),
        algebra: super::conversion::wrap_quiver_algebra(algebra.into_algebra()),
        kind,
    }
}

fn nakayama_for<F: runtime::PythonField>(
    field: F,
    sequence: &[Integer],
) -> PyResult<PyNakayamaAlgebraResult>
where
    super::conversion::QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    let result = qpa_gbnp::nakayama_algebra(field, sequence).map_err(invalid)?;
    let diagnostic = result.diagnostic();
    Ok(match result {
        NakayamaAlgebraResult::Algebra(algebra) => PyNakayamaAlgebraResult {
            algebra: Some(wrap_predefined(algebra)),
            invalid_sequence: None,
            diagnostic: None,
        },
        NakayamaAlgebraResult::InvalidSequence(sequence) => PyNakayamaAlgebraResult {
            algebra: None,
            invalid_sequence: Some(sequence),
            diagnostic,
        },
    })
}

/// Construct QPA's Nakayama algebra or preserve its invalid-sequence result.
#[pyfunction]
#[pyo3(name = "nakayama_algebra")]
pub(crate) fn py_nakayama_algebra(
    py: Python<'_>,
    field: &PyField,
    admissible_sequence: Vec<Py<PyAny>>,
) -> PyResult<PyNakayamaAlgebraResult> {
    let sequence = parse_python_integer_vector(py, admissible_sequence)?;
    match field.kind {
        FieldKind::Rational => nakayama_for(Q, &sequence),
        FieldKind::Binary => nakayama_for(Z2, &sequence),
        FieldKind::Prime(characteristic) => nakayama_for(Zp::new(characteristic), &sequence),
    }
}

fn canonical_for<F: runtime::PythonField>(
    py: Python<'_>,
    field: F,
    weights: &[Integer],
    relation_coefficients: Option<Vec<Py<PyScalar>>>,
) -> PyResult<PyPredefinedAlgebra>
where
    super::conversion::QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    let algebra = if let Some(coefficients) = relation_coefficients {
        let coefficients = coefficients
            .into_iter()
            .map(|coefficient| field.scalar_value(&coefficient.borrow(py)))
            .collect::<PyResult<Vec<_>>>()?;
        qpa_gbnp::canonical_algebra(field, weights, &coefficients)
    } else {
        qpa_gbnp::canonical_algebra_two_weights(field, weights)
    }
    .map_err(invalid)?;
    Ok(wrap_predefined(algebra))
}

/// Construct QPA's canonical algebra, using its two-weight overload when omitted.
#[pyfunction]
#[pyo3(name = "canonical_algebra", signature = (field, weights, relation_coefficients=None))]
pub(crate) fn py_canonical_algebra(
    py: Python<'_>,
    field: &PyField,
    weights: Vec<Py<PyAny>>,
    relation_coefficients: Option<Vec<Py<PyScalar>>>,
) -> PyResult<PyPredefinedAlgebra> {
    let weights = parse_python_integer_vector(py, weights)?;
    match field.kind {
        FieldKind::Rational => canonical_for(py, Q, &weights, relation_coefficients),
        FieldKind::Binary => canonical_for(py, Z2, &weights, relation_coefficients),
        FieldKind::Prime(characteristic) => {
            canonical_for(py, Zp::new(characteristic), &weights, relation_coefficients)
        }
    }
}

fn kronecker_for<F: runtime::PythonField>(
    field: F,
    number_of_arrows: &Integer,
) -> PyResult<PyPredefinedAlgebra>
where
    super::conversion::QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    Ok(wrap_predefined(
        qpa_gbnp::kronecker_algebra(field, number_of_arrows).map_err(invalid)?,
    ))
}

/// Construct QPA's generalized Kronecker algebra.
#[pyfunction]
#[pyo3(name = "kronecker_algebra")]
pub(crate) fn py_kronecker_algebra(
    field: &PyField,
    number_of_arrows: &Bound<'_, PyAny>,
) -> PyResult<PyPredefinedAlgebra> {
    let number_of_arrows = parse_python_integer(number_of_arrows)?;
    match field.kind {
        FieldKind::Rational => kronecker_for(Q, &number_of_arrows),
        FieldKind::Binary => kronecker_for(Z2, &number_of_arrows),
        FieldKind::Prime(characteristic) => {
            kronecker_for(Zp::new(characteristic), &number_of_arrows)
        }
    }
}

pub(super) fn preprojective_for<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
) -> PyResult<PyQuiverAlgebra>
where
    super::conversion::QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    Ok(super::conversion::wrap_quiver_algebra(
        qpa_gbnp::preprojective_algebra(algebra).map_err(invalid)?,
    ))
}

fn finite_global_for<F: runtime::PythonField>(
    field: F,
    number: &Integer,
    constructor: impl FnOnce(F, &Integer) -> Result<QuiverAlgebra<F>, qpa_gbnp::PredefinedAlgebraError>,
) -> PyResult<PyQuiverAlgebra>
where
    super::conversion::QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    Ok(super::conversion::wrap_quiver_algebra(
        constructor(field, number).map_err(invalid)?,
    ))
}

macro_rules! finite_global_constructor {
    ($python:ident, $python_name:literal, $qpa:path) => {
        #[pyfunction]
        #[pyo3(name = $python_name)]
        pub(crate) fn $python(
            field: &PyField,
            number: &Bound<'_, PyAny>,
        ) -> PyResult<PyQuiverAlgebra> {
            let number = parse_python_integer(number)?;
            match field.kind {
                FieldKind::Rational => finite_global_for(Q, &number, $qpa),
                FieldKind::Binary => finite_global_for(Z2, &number, $qpa),
                FieldKind::Prime(characteristic) => {
                    finite_global_for(Zp::new(characteristic), &number, $qpa)
                }
            }
        }
    };
}

finite_global_constructor!(
    py_finite_global_dim_green,
    "finite_global_dim_green",
    qpa_gbnp::finite_global_dim_green
);
finite_global_constructor!(
    py_finite_global_dim_kirk_kuz,
    "finite_global_dim_kirk_kuz",
    qpa_gbnp::finite_global_dim_kirk_kuz
);

/// Enumerate QPA's admissible sequences with arbitrary-precision parameters.
#[pyfunction]
#[pyo3(name = "admissible_sequence_generator")]
pub(crate) fn py_admissible_sequence_generator(
    py: Python<'_>,
    number_of_vertices: &Bound<'_, PyAny>,
    height: &Bound<'_, PyAny>,
) -> PyResult<Vec<Vec<Py<PyAny>>>> {
    let sequences = qpa_gbnp::admissible_sequence_generator(
        &parse_python_integer(number_of_vertices)?,
        &parse_python_integer(height)?,
    )
    .map_err(invalid)?;
    wrap_python_integer_matrix(py, sequences.into_iter())
}
