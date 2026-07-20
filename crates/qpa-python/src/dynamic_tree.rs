//! Python wrappers for QPA's mutable dynamic trees.

use std::sync::Arc;

use pyo3::{Py, PyAny, PyResult, Python, prelude::*};
use qpa_core::{DynamicForest, DynamicTree};

use crate::invalid;

type PythonDynamicForest = DynamicForest<Arc<Py<PyAny>>>;
type PythonDynamicTree = DynamicTree<Arc<Py<PyAny>>>;

/// Ownership arena for QPA dynamic-tree nodes.
#[pyclass(name = "DynamicForest", frozen, unsendable, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyDynamicForest {
    inner: PythonDynamicForest,
}

/// One identity-bearing QPA dynamic-tree node.
#[pyclass(name = "DynamicTree", frozen, unsendable, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyDynamicTree {
    inner: PythonDynamicTree,
}

impl PyDynamicTree {
    fn wrap_related(value: Option<PythonDynamicTree>) -> Option<Self> {
        value.map(|inner| Self { inner })
    }
}

#[pymethods]
impl PyDynamicForest {
    #[new]
    fn new() -> Self {
        Self {
            inner: DynamicForest::new(),
        }
    }

    #[classattr]
    #[pyo3(name = "family_name")]
    const FAMILY_NAME: &'static str = "DynamicTreeFamily";

    #[pyo3(signature = (item=None))]
    fn create(&self, item: Option<Py<PyAny>>) -> PyDynamicTree {
        PyDynamicTree {
            inner: self.inner.create(item.map(Arc::new)),
        }
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    #[getter]
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

#[pymethods]
impl PyDynamicTree {
    #[getter]
    #[allow(clippy::unused_self)]
    fn is_dynamic_tree(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn has_item(&self) -> bool {
        self.inner.has_item()
    }

    #[getter]
    fn item(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.inner.item().map(|item| item.as_ref().clone_ref(py))
    }

    #[setter]
    fn set_item(&self, item: Option<Py<PyAny>>) {
        self.inner.replace_item(item.map(Arc::new));
    }

    #[getter]
    fn parent(&self) -> Option<Self> {
        Self::wrap_related(self.inner.parent())
    }

    #[setter]
    fn set_parent(&self, parent: Option<&Self>) -> PyResult<()> {
        self.inner
            .set_parent(parent.map(|parent| &parent.inner))
            .map_err(invalid)
    }

    #[getter]
    fn left(&self) -> Option<Self> {
        Self::wrap_related(self.inner.left())
    }

    #[setter]
    fn set_left(&self, left: Option<&Self>) -> PyResult<()> {
        self.inner
            .set_left(left.map(|left| &left.inner))
            .map_err(invalid)
    }

    #[getter]
    fn right(&self) -> Option<Self> {
        Self::wrap_related(self.inner.right())
    }

    #[setter]
    fn set_right(&self, right: Option<&Self>) -> PyResult<()> {
        self.inner
            .set_right(right.map(|right| &right.inner))
            .map_err(invalid)
    }

    fn right_rotate(&self) -> PyResult<()> {
        self.inner.right_rotate().map_err(invalid)
    }

    fn left_rotate(&self) -> PyResult<()> {
        self.inner.left_rotate().map_err(invalid)
    }

    fn splice(&self) -> PyResult<()> {
        self.inner.splice().map_err(invalid)
    }

    fn splay(&self) -> PyResult<()> {
        self.inner.splay().map_err(invalid)
    }

    fn root(&self) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.root().map_err(invalid)?,
        })
    }

    fn link(&self, parent: &Self) -> PyResult<()> {
        self.inner.link(&parent.inner).map_err(invalid)
    }

    fn cut(&self) -> PyResult<()> {
        self.inner.cut().map_err(invalid)
    }

    fn is_identical(&self, other: &Self) -> bool {
        self.inner.is_identical(&other.inner)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}
