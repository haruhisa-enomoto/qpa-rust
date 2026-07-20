//! Python adapters for QPA's generic list utilities.

use pyo3::{
    Py, PyAny, PyResult, Python, pyfunction,
    types::{PyAnyMethods, PyList, PyListMethods, PyTuple, PyTupleMethods},
};

use crate::PyScalar;

fn python_is_zero(value: &pyo3::Bound<'_, PyAny>) -> PyResult<bool> {
    if let Ok(scalar) = value.extract::<pyo3::PyRef<'_, PyScalar>>() {
        return Ok(scalar.is_zero_value());
    }
    if let Ok(list) = value.cast::<PyList>() {
        for item in list.iter() {
            if !python_is_zero(&item)? {
                return Ok(false);
            }
        }
        return Ok(true);
    }
    if let Ok(tuple) = value.cast::<PyTuple>() {
        for item in tuple.iter() {
            if !python_is_zero(&item)? {
                return Ok(false);
            }
        }
        return Ok(true);
    }
    if let Ok(attribute) = value.getattr("is_zero") {
        let result = if attribute.is_callable() {
            attribute.call0()?
        } else {
            attribute
        };
        if let Ok(is_zero) = result.extract::<bool>() {
            return Ok(is_zero);
        }
    }
    value.eq(0)
}

/// Translate QPA's generic `PositionsNonzero` method.
///
/// Positions remain one-based because they are QPA positions, not Python
/// sequence indices.
#[pyfunction(name = "positions_nonzero")]
pub(crate) fn py_positions_nonzero(py: Python<'_>, values: Vec<Py<PyAny>>) -> PyResult<Vec<usize>> {
    let zeroes = values
        .iter()
        .map(|value| python_is_zero(value.bind(py)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(qpa_core::positions_nonzero(&zeroes, |is_zero| *is_zero))
}
