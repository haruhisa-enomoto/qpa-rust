//! Python dispatch for QPA string validity and local AR moves.

use pyo3::{
    PyAny, pyfunction, pymethods,
    types::{PyList, PyString},
};
use qpa_gbnp::{QpaStringExtension, StringSigmaEpsilon};

#[allow(clippy::wildcard_imports)]
// The quiver-algebra runtime union belongs to the parent binding module.
use super::*;
use crate::quiver::PyQuiverArrow;

/// QPA's concrete sigma and epsilon sign vectors in quiver-arrow order.
#[pyclass(name = "StringSigmaEpsilon", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyStringSigmaEpsilon {
    inner: StringSigmaEpsilon,
}

#[pymethods]
impl PyStringSigmaEpsilon {
    #[getter]
    fn sigma(&self) -> Vec<i8> {
        self.inner.sigma().to_vec()
    }

    #[getter]
    fn epsilon(&self) -> Vec<i8> {
        self.inner.epsilon().to_vec()
    }

    fn __repr__(&self) -> String {
        format!(
            "StringSigmaEpsilon(sigma={:?}, epsilon={:?})",
            self.inner.sigma(),
            self.inner.epsilon()
        )
    }
}

fn extension_to_python(py: Python<'_>, extension: QpaStringExtension) -> PyResult<Py<PyAny>> {
    match extension {
        QpaStringExtension::CannotPerform => Ok(PyString::new(py, "Cannot Perform The Operation")
            .into_any()
            .unbind()),
        QpaStringExtension::Word(word) => Ok(PyString::new(py, &word).into_any().unbind()),
        QpaStringExtension::ArrowList(arrows) => {
            let arrows = arrows
                .into_iter()
                .map(|inner| Py::new(py, PyQuiverArrow { inner }))
                .collect::<PyResult<Vec<_>>>()?;
            Ok(PyList::new(py, arrows)?.into_any().unbind())
        }
    }
}

pub(super) fn sigma_epsilon_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<PyStringSigmaEpsilon> {
    let inner = match algebra {
        QuiverAlgebra::Path(algebra) => {
            qpa_gbnp::qpa_string_sigma_epsilon(algebra).map_err(invalid)?
        }
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::qpa_string_sigma_epsilon(algebra).map_err(invalid)?
        }
    };
    Ok(PyStringSigmaEpsilon { inner })
}

pub(super) fn is_valid_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
    input: &str,
) -> PyResult<bool> {
    match algebra {
        QuiverAlgebra::Path(algebra) => qpa_gbnp::is_valid_string(algebra, input).map_err(invalid),
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::is_valid_string(algebra, input).map_err(invalid)
        }
    }
}

macro_rules! extension_dispatch {
    ($name:ident, $qpa:path) => {
        pub(super) fn $name<F: runtime::PythonField>(
            py: Python<'_>,
            algebra: &QuiverAlgebra<F>,
            input: &str,
            signs: &PyStringSigmaEpsilon,
        ) -> PyResult<Py<PyAny>> {
            let extension = match algebra {
                QuiverAlgebra::Path(algebra) => $qpa(algebra, input, &signs.inner),
                QuiverAlgebra::Quotient(algebra) => $qpa(algebra, input, &signs.inner),
            }
            .map_err(invalid)?;
            extension_to_python(py, extension)
        }
    };
}

extension_dispatch!(direct_left_for, qpa_gbnp::qpa_string_direct_left);
extension_dispatch!(direct_right_for, qpa_gbnp::qpa_string_direct_right);
extension_dispatch!(inverse_left_for, qpa_gbnp::qpa_string_inverse_left);
extension_dispatch!(inverse_right_for, qpa_gbnp::qpa_string_inverse_right);

pub(super) fn local_ar_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
    input: &str,
) -> PyResult<Vec<String>> {
    match algebra {
        QuiverAlgebra::Path(algebra) => qpa_gbnp::local_ar_quiver(algebra, input).map_err(invalid),
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::local_ar_quiver(algebra, input).map_err(invalid)
        }
    }
}

/// Translate QPA's relation-word-to-quiver-path helper.
#[pyfunction]
#[pyo3(name = "string_quiver_rho")]
pub(crate) fn py_string_quiver_rho(
    quiver: &PyQuiver,
    relations: Vec<String>,
) -> PyResult<Vec<PyPath>> {
    Ok(qpa_gbnp::qpa_string_quiver_rho(&quiver.inner, &relations)
        .map_err(invalid)?
        .into_iter()
        .map(|inner| PyPath { inner })
        .collect())
}
