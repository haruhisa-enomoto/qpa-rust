//! Python representation of QPA's finite posets.

use std::fmt;

use numerica::domains::integer::Integer;
use pyo3::{
    prelude::*,
    types::{PyInt, PyString},
};
use qpa_core::Poset;

use crate::{
    field::{parse_python_integer, wrap_python_integer},
    invalid,
};

/// Exact Python point types supported by the QPA poset boundary.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) enum PyPosetElement {
    Integer(Integer),
    String(String),
}

impl fmt::Display for PyPosetElement {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Integer(value) => value.fmt(formatter),
            Self::String(value) => value.fmt(formatter),
        }
    }
}

fn parse_element(value: &Bound<'_, PyAny>) -> PyResult<PyPosetElement> {
    if value.is_instance_of::<PyInt>() {
        return Ok(PyPosetElement::Integer(parse_python_integer(value)?));
    }
    if value.is_instance_of::<PyString>() {
        return Ok(PyPosetElement::String(value.extract()?));
    }
    Err(invalid("a poset point must be an exact integer or string"))
}

fn parse_elements(py: Python<'_>, values: Vec<Py<PyAny>>) -> PyResult<Vec<PyPosetElement>> {
    values
        .into_iter()
        .map(|value| parse_element(value.bind(py)))
        .collect()
}

fn wrap_element(py: Python<'_>, value: &PyPosetElement) -> PyResult<Py<PyAny>> {
    match value {
        PyPosetElement::Integer(value) => wrap_python_integer(py, value),
        PyPosetElement::String(value) => Ok(PyString::new(py, value).into_any().unbind()),
    }
}

fn wrap_elements(py: Python<'_>, values: &[PyPosetElement]) -> PyResult<Vec<Py<PyAny>>> {
    values.iter().map(|value| wrap_element(py, value)).collect()
}

/// QPA's callable `PartialOrderOfPoset` attribute.
#[pyclass(name = "PosetOrder", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPosetOrder {
    poset: PyPoset,
}

#[pymethods]
impl PyPosetOrder {
    fn __call__(&self, left: &Bound<'_, PyAny>, right: &Bound<'_, PyAny>) -> PyResult<bool> {
        self.poset
            .inner
            .less_than_or_equal(&parse_element(left)?, &parse_element(right)?)
            .map_err(invalid)
    }

    #[allow(clippy::unused_self)]
    fn __repr__(&self) -> &'static str {
        "<partial order of a QPA poset>"
    }
}

/// QPA's finite `IsPoset` object over exact Python integers and strings.
#[pyclass(name = "Poset", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPoset {
    pub(crate) inner: Poset<PyPosetElement>,
}

#[pymethods]
impl PyPoset {
    #[new]
    fn new(
        py: Python<'_>,
        points: Vec<Py<PyAny>>,
        relations: Vec<Vec<Py<PyAny>>>,
    ) -> PyResult<Self> {
        let points = parse_elements(py, points)?;
        let relations = relations
            .into_iter()
            .map(|row| parse_elements(py, row))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self {
            inner: Poset::new(points, &relations).map_err(invalid)?,
        })
    }

    #[classattr]
    #[allow(non_upper_case_globals)]
    const family_name: &'static str = "TheFamilyOfPosets";

    #[classattr]
    #[allow(non_upper_case_globals)]
    const type_name: &'static str = "TheTypePoset";

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_poset(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    fn __len__(&self) -> usize {
        self.inner.size()
    }

    #[getter]
    fn underlying_set(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        wrap_elements(py, self.inner.underlying_set())
    }

    #[getter]
    fn minimal_relations(&self, py: Python<'_>) -> PyResult<Vec<(Py<PyAny>, Py<PyAny>)>> {
        self.inner
            .minimal_relations()
            .iter()
            .map(|(left, right)| Ok((wrap_element(py, left)?, wrap_element(py, right)?)))
            .collect()
    }

    #[getter]
    fn partial_order(&self) -> PyPosetOrder {
        PyPosetOrder {
            poset: self.clone(),
        }
    }

    fn less_than_or_equal(
        &self,
        left: &Bound<'_, PyAny>,
        right: &Bound<'_, PyAny>,
    ) -> PyResult<bool> {
        self.inner
            .less_than_or_equal(&parse_element(left)?, &parse_element(right)?)
            .map_err(invalid)
    }

    fn is_identical(&self, other: &Self) -> bool {
        self.inner.is_identical(&other.inner)
    }

    fn display(&self) -> String {
        let underlying = self
            .inner
            .underlying_set()
            .iter()
            .map(|element| match element {
                PyPosetElement::Integer(value) => value.to_string(),
                PyPosetElement::String(value) => format!("{value:?}"),
            })
            .collect::<Vec<_>>()
            .join(", ");
        let minimal = self
            .inner
            .minimal_relations()
            .iter()
            .map(|(left, right)| format!("DirectProductElement( [ {left}, {right} ] )"))
            .collect::<Vec<_>>()
            .join(", ");
        format!(
            "A poset with underlying set:\n[ {underlying} ]\nand partial order:\n[ {minimal} ]\n"
        )
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    fn view(&self) -> String {
        self.__repr__()
    }
}
