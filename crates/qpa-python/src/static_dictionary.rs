//! Python wrappers for QPA's generic and quiver static dictionaries.

use pyo3::{
    Bound, Py, PyAny, PyResult, Python,
    prelude::{PyAnyMethods, PyStringMethods, pyclass, pyfunction, pymethods},
    types::PyString,
};
use qpa_core::{DictionaryMatch, QuiverStaticDictionary, StaticDictionary};

use crate::{
    field::wrap_python_integer,
    invalid,
    quiver::{PyPath, PyQuiver},
};

fn sequence_symbols(py: Python<'_>, sequence: &Bound<'_, PyAny>) -> PyResult<Vec<Py<PyAny>>> {
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
        .map(|item| item.map(Bound::unbind))
        .collect()
}

fn intern_symbol(
    py: Python<'_>,
    representatives: &mut Vec<Py<PyAny>>,
    symbol: &Py<PyAny>,
) -> PyResult<usize> {
    for (index, representative) in representatives.iter().enumerate() {
        if representative.bind(py).eq(symbol.bind(py))? {
            return Ok(index);
        }
    }
    let index = representatives.len();
    representatives.push(symbol.clone_ref(py));
    Ok(index)
}

fn matches_to_python(matches: Vec<DictionaryMatch>) -> Vec<(usize, Vec<usize>)> {
    matches
        .into_iter()
        .map(|matched| (matched.end_position, matched.pattern_indices))
        .collect()
}

/// QPA's generic dense-sequence `StaticDictionary`.
#[pyclass(name = "StaticDictionary", frozen, skip_from_py_object)]
pub(crate) struct PyStaticDictionary {
    patterns: Vec<Py<PyAny>>,
    representatives: Vec<Py<PyAny>>,
    inner: StaticDictionary<usize>,
}

#[pymethods]
impl PyStaticDictionary {
    #[new]
    fn new(py: Python<'_>, patterns: Vec<Py<PyAny>>) -> PyResult<Self> {
        let mut representatives = Vec::new();
        let mut encoded = Vec::with_capacity(patterns.len());
        for pattern in &patterns {
            encoded.push(
                sequence_symbols(py, pattern.bind(py))?
                    .iter()
                    .map(|symbol| intern_symbol(py, &mut representatives, symbol))
                    .collect::<PyResult<Vec<_>>>()?,
            );
        }
        Ok(Self {
            patterns,
            representatives,
            inner: StaticDictionary::new(&encoded),
        })
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_static_dictionary(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn patterns(&self, py: Python<'_>) -> Vec<Py<PyAny>> {
        self.patterns
            .iter()
            .map(|pattern| pattern.clone_ref(py))
            .collect()
    }

    fn search(
        &self,
        py: Python<'_>,
        sequence: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<(usize, Vec<usize>)>> {
        let symbols = sequence_symbols(py, sequence)?;
        let unknown = self.representatives.len();
        let mut encoded = Vec::with_capacity(symbols.len());
        for symbol in symbols {
            let mut token = unknown;
            for (index, representative) in self.representatives.iter().enumerate() {
                if representative.bind(py).eq(symbol.bind(py))? {
                    token = index;
                    break;
                }
            }
            encoded.push(token);
        }
        Ok(matches_to_python(
            self.inner.search(&encoded).map_err(invalid)?,
        ))
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }
}

/// QPA's path-specialized finite-automaton dictionary.
#[pyclass(name = "QuiverStaticDictionary", frozen, skip_from_py_object)]
pub(crate) struct PyQuiverStaticDictionary {
    inner: QuiverStaticDictionary,
}

#[pymethods]
impl PyQuiverStaticDictionary {
    #[new]
    fn new(py: Python<'_>, quiver: &PyQuiver, patterns: Vec<Py<PyPath>>) -> PyResult<Self> {
        let patterns = patterns
            .into_iter()
            .map(|pattern| pattern.borrow(py).inner.clone())
            .collect::<Vec<_>>();
        Ok(Self {
            inner: QuiverStaticDictionary::new(&quiver.inner, &patterns).map_err(invalid)?,
        })
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_quiver_static_dictionary(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn patterns(&self) -> Vec<PyPath> {
        self.inner
            .patterns()
            .iter()
            .cloned()
            .map(|inner| PyPath { inner })
            .collect()
    }

    fn search(&self, path: &PyPath) -> PyResult<Vec<(usize, Vec<usize>)>> {
        Ok(matches_to_python(
            self.inner.search(&path.inner).map_err(invalid)?,
        ))
    }

    fn prefix_search(&self, path: &PyPath) -> PyResult<Vec<(usize, Vec<usize>)>> {
        Ok(matches_to_python(
            self.inner.prefix_search(&path.inner).map_err(invalid)?,
        ))
    }

    #[getter]
    fn is_finite_difference(&self) -> bool {
        self.inner.is_finite_difference()
    }

    #[getter]
    fn difference_size(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match self.inner.difference_size() {
            Some(size) => wrap_python_integer(py, &size),
            None => Ok(PyString::new(py, "infinity").into_any().unbind()),
        }
    }

    #[getter]
    fn difference_words(&self) -> Option<Vec<PyPath>> {
        self.inner
            .difference_words()
            .map(|paths| paths.into_iter().map(|inner| PyPath { inner }).collect())
    }

    fn __repr__(&self) -> String {
        self.inner.view_string()
    }
}

/// Exact level-one messages emitted through QPA's `InfoStaticDictionary`.
#[pyfunction(name = "static_dictionary_info_messages")]
pub(crate) fn py_static_dictionary_info_messages() -> (&'static str, &'static str) {
    (
        qpa_core::STATIC_DICTIONARY_CREATING_AUTOMATA_INFO,
        qpa_core::STATIC_DICTIONARY_CHECKING_CYCLE_INFO,
    )
}
