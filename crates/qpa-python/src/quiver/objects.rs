//! Concrete Python objects for QPA quiver generators and path iteration.

use pyo3::{
    Py, PyAny, PyResult, Python,
    exceptions::PyIndexError,
    prelude::{PyRef, PyRefMut, pyclass, pymethods},
};
use qpa_core::{Arrow, Generator, PathEndpoint, PathIterator, Quiver, Vertex};

use super::{PyPath, PyQuiver};
use crate::invalid;

/// A concrete QPA quiver vertex rather than only its display name.
#[pyclass(name = "QuiverVertex", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuiverVertex {
    pub(crate) inner: Vertex,
}

impl PyQuiverVertex {
    pub(super) const fn from_inner(inner: Vertex) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyQuiverVertex {
    #[getter]
    fn quiver(&self) -> PyQuiver {
        PyQuiver::from_inner(self.inner.quiver().clone())
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.name().to_owned()
    }

    #[getter]
    fn index(&self) -> usize {
        self.inner.index().get()
    }

    #[getter]
    fn generator_index(&self) -> usize {
        self.inner.index().get()
    }

    #[getter]
    fn in_degree(&self) -> usize {
        self.inner.in_degree()
    }

    #[getter]
    fn out_degree(&self) -> usize {
        self.inner.out_degree()
    }

    #[getter]
    fn incoming_arrows(&self) -> PyResult<Vec<PyQuiverArrow>> {
        Ok(self
            .inner
            .quiver()
            .incoming_arrows(&self.inner)
            .map_err(invalid)?
            .into_iter()
            .map(PyQuiverArrow::from_inner)
            .collect())
    }

    #[getter]
    fn outgoing_arrows(&self) -> PyResult<Vec<PyQuiverArrow>> {
        Ok(self
            .inner
            .quiver()
            .outgoing_arrows(&self.inner)
            .map_err(invalid)?
            .into_iter()
            .map(PyQuiverArrow::from_inner)
            .collect())
    }

    #[getter]
    fn neighbors(&self) -> PyResult<Vec<Self>> {
        Ok(self
            .inner
            .quiver()
            .neighbors(&self.inner)
            .map_err(invalid)?
            .into_iter()
            .map(Self::from_inner)
            .collect())
    }

    fn path(&self) -> PyPath {
        PyPath {
            inner: self.inner.path(),
        }
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

/// A concrete QPA quiver arrow rather than only its display name and endpoints.
#[pyclass(name = "QuiverArrow", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuiverArrow {
    pub(crate) inner: Arrow,
}

impl PyQuiverArrow {
    pub(super) const fn from_inner(inner: Arrow) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyQuiverArrow {
    #[getter]
    fn quiver(&self) -> PyQuiver {
        PyQuiver::from_inner(self.inner.quiver().clone())
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.name().to_owned()
    }

    #[getter]
    fn index(&self) -> usize {
        self.inner.index().get()
    }

    #[getter]
    fn generator_index(&self) -> usize {
        self.inner.generator_index()
    }

    #[getter]
    fn source(&self) -> PyQuiverVertex {
        PyQuiverVertex::from_inner(self.inner.source())
    }

    #[getter]
    fn target(&self) -> PyQuiverVertex {
        PyQuiverVertex::from_inner(self.inner.target())
    }

    fn path(&self) -> PyPath {
        PyPath {
            inner: self.inner.path(),
        }
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

pub(super) fn generator_to_python(py: Python<'_>, generator: Generator) -> PyResult<Py<PyAny>> {
    match generator {
        Generator::Vertex(vertex) => {
            Ok(Py::new(py, PyQuiverVertex::from_inner(vertex))?.into_any())
        }
        Generator::Arrow(arrow) => Ok(Py::new(py, PyQuiverArrow::from_inner(arrow))?.into_any()),
    }
}

pub(super) fn endpoint_to_python(endpoint: PathEndpoint) -> Option<PyQuiverVertex> {
    match endpoint {
        PathEndpoint::Vertex(vertex) => Some(PyQuiverVertex::from_inner(vertex)),
        PathEndpoint::Zero(_) => None,
    }
}

/// Mutable QPA breadth-first path iterator.
#[pyclass(name = "QuiverPathIterator", skip_from_py_object)]
pub(crate) struct PyQuiverPathIterator {
    inner: PathIterator,
}

impl PyQuiverPathIterator {
    pub(super) const fn from_inner(inner: PathIterator) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyQuiverPathIterator {
    #[getter]
    fn is_done(&self) -> bool {
        self.inner.is_done()
    }

    fn next_path(&mut self) -> Option<PyPath> {
        self.inner.next().map(|inner| PyPath { inner })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyPath> {
        slf.next_path()
    }
}

/// Stateless QPA quiver enumerator with one-based `at` and Python indexing.
#[pyclass(name = "QuiverPathEnumerator", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuiverPathEnumerator {
    quiver: Quiver,
}

impl PyQuiverPathEnumerator {
    pub(super) const fn new(quiver: Quiver) -> Self {
        Self { quiver }
    }
}

#[pymethods]
impl PyQuiverPathEnumerator {
    #[getter]
    fn quiver(&self) -> PyQuiver {
        PyQuiver::from_inner(self.quiver.clone())
    }

    fn at(&self, position: usize) -> Option<PyPath> {
        self.quiver.path_at(position).map(|inner| PyPath { inner })
    }

    fn __getitem__(&self, index: isize) -> PyResult<PyPath> {
        let position = usize::try_from(index)
            .ok()
            .and_then(|index| index.checked_add(1))
            .ok_or_else(|| PyIndexError::new_err("negative path indices are not supported"))?;
        self.at(position)
            .ok_or_else(|| PyIndexError::new_err("path index is outside this quiver enumerator"))
    }

    fn __iter__(&self) -> PyQuiverPathIterator {
        PyQuiverPathIterator::from_inner(self.quiver.paths())
    }

    fn __repr__(&self) -> String {
        format!("<path enumerator of {}>", self.quiver)
    }
}
