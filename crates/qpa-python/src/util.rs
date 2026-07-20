//! Python adapters for QPA's generic list utilities.

use pyo3::{
    Py, PyAny, PyResult, Python, pyfunction,
    types::{
        PyAnyMethods, PyList, PyListMethods, PyString, PyStringMethods, PyTuple, PyTupleMethods,
    },
};

use crate::PyScalar;

pub(crate) fn sequence_symbols(
    py: Python<'_>,
    sequence: &pyo3::Bound<'_, PyAny>,
) -> PyResult<Vec<Py<PyAny>>> {
    if let Ok(string) = sequence.cast::<PyString>() {
        return string
            .to_str()?
            .chars()
            .map(|character| {
                Ok(PyString::new(py, &character.to_string())
                    .into_any()
                    .unbind())
            })
            .collect();
    }
    sequence
        .try_iter()?
        .map(|item| item.map(pyo3::Bound::unbind))
        .collect()
}

pub(crate) fn intern_symbol(
    py: Python<'_>,
    representatives: &mut Vec<Py<PyAny>>,
    symbol: &Py<PyAny>,
) -> PyResult<usize> {
    if let Some(index) = symbol_index(py, representatives, symbol)? {
        return Ok(index);
    }
    let index = representatives.len();
    representatives.push(symbol.clone_ref(py));
    Ok(index)
}

pub(crate) fn symbol_index(
    py: Python<'_>,
    representatives: &[Py<PyAny>],
    symbol: &Py<PyAny>,
) -> PyResult<Option<usize>> {
    for (index, representative) in representatives.iter().enumerate() {
        if representative.bind(py).eq(symbol.bind(py))? {
            return Ok(Some(index));
        }
    }
    Ok(None)
}

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
