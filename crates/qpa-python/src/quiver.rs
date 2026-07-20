//! Python wrappers for QPA quivers and paths.

mod objects;
mod ordering;

use std::cmp::Ordering;

use pyo3::{Py, PyAny, prelude::*};
use qpa_core::{
    Arrow, ArrowIndex, ArrowSpec, DynkinOrientation, Generator, NamedArrowSpec, Path, PathEndpoint,
    Quiver, QuiverCardinality, QuiverOrdering, QuiverProduct, QuiverProductDecomposition, Vertex,
    VertexIndex, include_in_product_quiver, project_from_product_quiver, quiver_product,
    reverse_path, unique_quiver_name, walk_of_path_or_vertex,
};

pub(crate) use objects::{
    PyQuiverArrow, PyQuiverPathEnumerator, PyQuiverPathIterator, PyQuiverVertex,
};
use objects::{endpoint_to_python, generator_to_python};
pub(crate) use ordering::PyQuiverOrdering;
use ordering::generators_from_positions;

use crate::invalid;

/// Allocate the same process-global name produced by QPA's `UniqueQuiverName`.
#[pyfunction]
#[pyo3(name = "unique_quiver_name")]
pub(crate) fn py_unique_quiver_name() -> String {
    unique_quiver_name()
}

/// Python-owned handle to a Rust QPA quiver identity.
#[pyclass(name = "Quiver", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuiver {
    pub(crate) inner: Quiver,
    ordering: Option<QuiverOrdering>,
}

/// A quiver constructed by QPA's Cartesian `QuiverProduct` method.
#[pyclass(name = "QuiverProduct", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuiverProduct {
    inner: QuiverProduct,
}

impl PyQuiverProduct {
    pub(crate) fn from_inner(inner: QuiverProduct) -> Self {
        Self { inner }
    }
}

/// QPA's two-factor decomposition attached to a product quiver.
#[pyclass(name = "QuiverProductDecomposition", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuiverProductDecomposition {
    inner: QuiverProductDecomposition,
}

#[pymethods]
impl PyQuiver {
    #[new]
    fn new(vertices: Vec<String>, arrows: Vec<(String, String, String)>) -> PyResult<Self> {
        let arrows = arrows
            .into_iter()
            .map(|(source, target, name)| NamedArrowSpec::named(source, target, name))
            .collect();
        let inner = Quiver::from_named(
            vertices.into_iter().map(String::into_boxed_str).collect(),
            arrows,
        )
        .map_err(invalid)?;
        Ok(Self::from_inner(inner))
    }

    fn product(&self, right: &Self) -> PyResult<PyQuiverProduct> {
        Ok(PyQuiverProduct {
            inner: quiver_product(&self.inner, &right.inner).map_err(invalid)?,
        })
    }

    #[staticmethod]
    fn named(
        vertices: Vec<String>,
        arrows: Vec<(String, String, Option<String>)>,
    ) -> PyResult<Self> {
        let arrows = arrows
            .into_iter()
            .map(|(source, target, name)| match name {
                Some(name) => NamedArrowSpec::named(source, target, name),
                None => NamedArrowSpec::unnamed(source, target),
            })
            .collect();
        let inner = Quiver::from_named(
            vertices.into_iter().map(String::into_boxed_str).collect(),
            arrows,
        )
        .map_err(invalid)?;
        Ok(Self::from_inner(inner))
    }

    #[staticmethod]
    fn indexed(
        number_of_vertices: usize,
        arrows: Vec<(usize, usize, Option<String>)>,
    ) -> PyResult<Self> {
        let arrows = arrows
            .into_iter()
            .map(|(source, target, name)| {
                let source = VertexIndex::new(source)
                    .ok_or_else(|| invalid("QPA vertex indices are one-based"))?;
                let target = VertexIndex::new(target)
                    .ok_or_else(|| invalid("QPA vertex indices are one-based"))?;
                Ok(name.map_or_else(
                    || ArrowSpec::unnamed(source, target),
                    |name| ArrowSpec::named(source, target, name),
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self::from_inner(
            Quiver::from_indexed(number_of_vertices, arrows).map_err(invalid)?,
        ))
    }

    #[staticmethod]
    fn from_adjacency(adjacency: Vec<Vec<usize>>) -> PyResult<Self> {
        Ok(Self::from_inner(
            Quiver::from_adjacency(adjacency).map_err(invalid)?,
        ))
    }

    #[staticmethod]
    fn dynkin(family: &str, rank: usize, orientation: Vec<String>) -> PyResult<Option<Self>> {
        let orientation = orientation
            .iter()
            .map(|symbol| parse_dynkin_orientation(symbol))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Quiver::dynkin_from_name(family, rank, &orientation)
            .map_err(invalid)?
            .map(Self::from_inner))
    }

    #[staticmethod]
    fn dynkin_a(rank: usize, orientation: Vec<String>) -> PyResult<Self> {
        Ok(Self::from_inner(
            Quiver::dynkin_a(rank, &parse_dynkin_orientations(&orientation)?).map_err(invalid)?,
        ))
    }

    #[staticmethod]
    fn dynkin_d(rank: usize, orientation: Vec<String>) -> PyResult<Self> {
        Ok(Self::from_inner(
            Quiver::dynkin_d(rank, &parse_dynkin_orientations(&orientation)?).map_err(invalid)?,
        ))
    }

    #[staticmethod]
    fn dynkin_e(rank: usize, orientation: Vec<String>) -> PyResult<Self> {
        Ok(Self::from_inner(
            Quiver::dynkin_e(rank, &parse_dynkin_orientations(&orientation)?).map_err(invalid)?,
        ))
    }

    #[getter]
    fn vertices(&self) -> Vec<String> {
        self.inner
            .vertices()
            .map(|vertex| vertex.name().to_owned())
            .collect()
    }

    #[getter]
    fn arrows(&self) -> Vec<(String, String, String)> {
        self.inner
            .arrows()
            .map(|arrow| {
                (
                    arrow.source().name().to_owned(),
                    arrow.target().name().to_owned(),
                    arrow.name().to_owned(),
                )
            })
            .collect()
    }

    /// Concrete `VerticesOfQuiver` objects in QPA generator order.
    #[getter]
    fn vertex_objects(&self) -> Vec<PyQuiverVertex> {
        self.inner
            .vertices()
            .map(PyQuiverVertex::from_inner)
            .collect()
    }

    /// Concrete `ArrowsOfQuiver` objects in QPA generator order.
    #[getter]
    fn arrow_objects(&self) -> Vec<PyQuiverArrow> {
        self.inner.arrows().map(PyQuiverArrow::from_inner).collect()
    }

    /// Vertices followed by arrows, exactly as in `GeneratorsOfQuiver`.
    #[getter]
    fn generators(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        self.inner
            .generators()
            .map(|generator| generator_to_python(py, generator))
            .collect()
    }

    fn vertex(&self, name: &str) -> PyResult<PyQuiverVertex> {
        Ok(PyQuiverVertex::from_inner(self.vertex_by_name(name)?))
    }

    fn vertex_at(&self, index: usize) -> PyResult<Option<PyQuiverVertex>> {
        let Some(index) = VertexIndex::new(index) else {
            return Err(invalid("QPA vertex indices are one-based"));
        };
        Ok(self.inner.vertex(index).map(PyQuiverVertex::from_inner))
    }

    fn arrow(&self, name: &str) -> PyResult<PyQuiverArrow> {
        Ok(PyQuiverArrow::from_inner(self.arrow_by_name(name)?))
    }

    fn arrow_at(&self, index: usize) -> PyResult<Option<PyQuiverArrow>> {
        let Some(index) = ArrowIndex::new(index) else {
            return Err(invalid("QPA arrow indices are one-based"));
        };
        Ok(self.inner.arrow(index).map(PyQuiverArrow::from_inner))
    }

    #[getter]
    fn number_of_vertices(&self) -> usize {
        self.inner.number_of_vertices()
    }

    #[getter]
    fn number_of_arrows(&self) -> usize {
        self.inner.number_of_arrows()
    }

    #[getter]
    fn adjacency_matrix(&self) -> Vec<Vec<usize>> {
        self.inner
            .adjacency_matrix()
            .iter()
            .map(|row| row.to_vec())
            .collect()
    }

    #[getter]
    fn is_acyclic(&self) -> bool {
        self.inner.is_acyclic()
    }

    #[getter]
    fn is_finite(&self) -> bool {
        self.inner.is_finite()
    }

    #[getter]
    fn is_connected(&self) -> bool {
        self.inner.is_connected()
    }

    #[getter]
    fn is_unoriented_acyclic(&self) -> bool {
        self.inner.is_unoriented_acyclic()
    }

    #[getter]
    fn is_tree(&self) -> bool {
        self.inner.is_tree()
    }

    #[getter]
    fn has_special_arrow_names(&self) -> bool {
        self.inner.has_special_single_letter_arrow_names()
    }

    #[getter]
    fn cardinality(&self) -> Option<String> {
        match self.inner.cardinality() {
            QuiverCardinality::Finite(value) => Some(value.to_string()),
            QuiverCardinality::Infinite => None,
        }
    }

    #[getter]
    fn topological_order(&self) -> Option<Vec<String>> {
        self.inner.topological_order().map(|vertices| {
            vertices
                .into_iter()
                .map(|vertex| vertex.name().to_owned())
                .collect()
        })
    }

    #[getter]
    fn dynkin_type(&self) -> Option<String> {
        self.inner.dynkin_type().map(|value| value.to_string())
    }

    fn incoming_arrows(&self, vertex: &str) -> PyResult<Vec<String>> {
        Ok(self
            .inner
            .incoming_arrows(&self.vertex_by_name(vertex)?)
            .map_err(invalid)?
            .into_iter()
            .map(|arrow| arrow.name().to_owned())
            .collect())
    }

    fn outgoing_arrows(&self, vertex: &str) -> PyResult<Vec<String>> {
        Ok(self
            .inner
            .outgoing_arrows(&self.vertex_by_name(vertex)?)
            .map_err(invalid)?
            .into_iter()
            .map(|arrow| arrow.name().to_owned())
            .collect())
    }

    fn neighbors(&self, vertex: &str) -> PyResult<Vec<String>> {
        Ok(self
            .inner
            .neighbors(&self.vertex_by_name(vertex)?)
            .map_err(invalid)?
            .into_iter()
            .map(|vertex| vertex.name().to_owned())
            .collect())
    }

    fn compare_paths(&self, left: &PyPath, right: &PyPath) -> PyResult<i8> {
        Ok(
            match self
                .inner
                .compare_paths(&left.inner, &right.inner)
                .map_err(invalid)?
            {
                Ordering::Less => -1,
                Ordering::Equal => 0,
                Ordering::Greater => 1,
            },
        )
    }

    #[getter]
    fn ordering(&self) -> PyResult<PyQuiverOrdering> {
        Ok(PyQuiverOrdering::from_inner(match &self.ordering {
            Some(ordering) => ordering.clone(),
            None => self.inner.ordering().map_err(invalid)?,
        }))
    }

    fn ordered_by(&self, ordering: &PyQuiverOrdering) -> PyResult<Self> {
        let ordered = self
            .inner
            .ordered_by(ordering.inner.clone())
            .map_err(invalid)?;
        let (inner, ordering) = ordered.into_parts().map_err(invalid)?;
        Ok(Self {
            inner,
            ordering: Some(ordering),
        })
    }

    fn vertex_path(&self, name: &str) -> PyResult<PyPath> {
        match self.inner.generator_by_name(name) {
            Some(Generator::Vertex(vertex)) => Ok(PyPath {
                inner: vertex.path(),
            }),
            _ => Err(invalid(format!("unknown vertex {name:?}"))),
        }
    }

    fn arrow_path(&self, name: &str) -> PyResult<PyPath> {
        match self.inner.generator_by_name(name) {
            Some(Generator::Arrow(arrow)) => Ok(PyPath {
                inner: arrow.path(),
            }),
            _ => Err(invalid(format!("unknown arrow {name:?}"))),
        }
    }

    fn path(&self, arrows: Vec<String>) -> PyResult<PyPath> {
        let mut arrows = arrows.into_iter();
        let first = arrows
            .next()
            .ok_or_else(|| invalid("an empty arrow list does not determine a vertex path"))?;
        let mut path = self.arrow_path(&first)?.inner;
        for name in arrows {
            path = path
                .multiply(&self.arrow_path(&name)?.inner)
                .map_err(invalid)?;
        }
        Ok(PyPath { inner: path })
    }

    fn path_from_generators(&self, positions: Vec<usize>) -> PyResult<PyPath> {
        Ok(PyPath {
            inner: Path::from_generator_list(
                &self.inner,
                &generators_from_positions(self, &positions)?,
            )
            .map_err(invalid)?,
        })
    }

    fn path_at(&self, position: usize) -> Option<PyPath> {
        self.inner.path_at(position).map(|inner| PyPath { inner })
    }

    fn iterator(&self) -> PyQuiverPathIterator {
        PyQuiverPathIterator::from_inner(self.inner.paths())
    }

    fn enumerator(&self) -> PyQuiverPathEnumerator {
        PyQuiverPathEnumerator::new(self.inner.clone())
    }

    fn path_from_external_representation(
        &self,
        representation: Vec<usize>,
    ) -> PyResult<Option<PyPath>> {
        Ok(
            Path::from_external_representation(&self.inner, &representation)
                .map_err(invalid)?
                .map(|inner| PyPath { inner }),
        )
    }

    fn full_subquiver(&self, vertices: Vec<String>) -> PyResult<Self> {
        let vertices = vertices
            .iter()
            .map(|name| self.vertex_by_name(name))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self::from_inner(
            self.inner.full_subquiver(&vertices).map_err(invalid)?,
        ))
    }

    fn connected_components(&self) -> PyResult<Vec<Self>> {
        Ok(self
            .inner
            .connected_components()
            .map_err(invalid)?
            .into_iter()
            .map(Self::from_inner)
            .collect())
    }

    fn opposite(&self) -> PyResult<Self> {
        Ok(Self::from_inner(self.inner.opposite().map_err(invalid)?))
    }

    fn opposite_name(&self, name: &str) -> PyResult<String> {
        Ok(self
            .inner
            .opposite_quiver_name(name)
            .map_err(invalid)?
            .into())
    }

    fn separated(&self) -> PyResult<Self> {
        Ok(Self::from_inner(self.inner.separated().map_err(invalid)?))
    }

    fn double(&self) -> PyResult<Self> {
        Ok(Self::from_inner(self.inner.double().map_err(invalid)?))
    }

    #[getter]
    fn is_special_biserial(&self) -> bool {
        qpa_gbnp::is_special_biserial_quiver(&self.inner)
    }

    fn zero_path(&self) -> PyPath {
        PyPath {
            inner: self.inner.zero(),
        }
    }

    fn is_identical(&self, other: &Self) -> bool {
        self.inner.is_identical(&other.inner)
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    fn view(&self) -> String {
        self.inner.view_string()
    }

    fn __iter__(&self) -> PyQuiverPathIterator {
        self.iterator()
    }
}

#[pymethods]
impl PyQuiverProduct {
    #[getter]
    fn quiver(&self) -> PyQuiver {
        PyQuiver::from_inner(self.inner.quiver().clone())
    }

    #[getter]
    fn decomposition(&self) -> PyQuiverProductDecomposition {
        PyQuiverProductDecomposition {
            inner: self.inner.decomposition().clone(),
        }
    }

    fn __repr__(&self) -> String {
        self.inner.quiver().to_string()
    }
}

#[pymethods]
impl PyQuiverProductDecomposition {
    #[getter]
    #[allow(clippy::unused_self)]
    fn is_quiver_product_decomposition(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    fn is_bound(&self, one_based_index: usize) -> bool {
        self.inner.is_bound(one_based_index)
    }

    fn factor(&self, one_based_index: usize) -> PyResult<PyQuiver> {
        Ok(PyQuiver::from_inner(
            self.inner.factor(one_based_index).map_err(invalid)?.clone(),
        ))
    }

    #[allow(clippy::unused_self)]
    fn __len__(&self) -> usize {
        2
    }

    #[allow(clippy::unused_self)]
    fn __repr__(&self) -> &'static str {
        "<quiver product decomposition>"
    }
}

impl PyQuiver {
    pub(crate) const fn from_inner(inner: Quiver) -> Self {
        Self {
            inner,
            ordering: None,
        }
    }

    pub(crate) const fn from_inner_with_ordering(inner: Quiver, ordering: QuiverOrdering) -> Self {
        Self {
            inner,
            ordering: Some(ordering),
        }
    }

    pub(crate) const fn rust_quiver(&self) -> &Quiver {
        &self.inner
    }

    pub(crate) const fn rust_ordering(&self) -> Option<&QuiverOrdering> {
        self.ordering.as_ref()
    }

    fn vertex_by_name(&self, name: &str) -> PyResult<Vertex> {
        match self.inner.generator_by_name(name) {
            Some(Generator::Vertex(vertex)) => Ok(vertex),
            _ => Err(invalid(format!("unknown vertex {name:?}"))),
        }
    }

    fn arrow_by_name(&self, name: &str) -> PyResult<Arrow> {
        match self.inner.generator_by_name(name) {
            Some(Generator::Arrow(arrow)) => Ok(arrow),
            _ => Err(invalid(format!("unknown arrow {name:?}"))),
        }
    }
}

fn parse_dynkin_orientations(symbols: &[String]) -> PyResult<Vec<DynkinOrientation>> {
    symbols
        .iter()
        .map(|symbol| parse_dynkin_orientation(symbol))
        .collect()
}

fn parse_dynkin_orientation(symbol: &str) -> PyResult<DynkinOrientation> {
    match symbol {
        "l" => Ok(DynkinOrientation::Left),
        "r" => Ok(DynkinOrientation::Right),
        "u" => Ok(DynkinOrientation::Up),
        "d" => Ok(DynkinOrientation::Down),
        _ => Err(invalid(format!(
            "invalid Dynkin orientation {symbol:?}; expected l, r, u, or d"
        ))),
    }
}

pub(crate) const fn ordering_value(ordering: Ordering) -> i8 {
    match ordering {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

/// A path using QPA's left-to-right multiplication convention.
#[pyclass(name = "Path", frozen, skip_from_py_object)]
pub(crate) struct PyPath {
    pub(crate) inner: Path,
}

fn endpoint_name(endpoint: PathEndpoint) -> Option<String> {
    match endpoint {
        PathEndpoint::Vertex(vertex) => Some(vertex.name().to_owned()),
        PathEndpoint::Zero(_) => None,
    }
}

#[pymethods]
impl PyPath {
    #[getter]
    fn quiver(&self) -> PyQuiver {
        PyQuiver::from_inner(self.inner.quiver().clone())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }

    #[getter]
    fn is_vertex(&self) -> bool {
        self.inner.is_vertex()
    }

    #[getter]
    fn is_arrow(&self) -> bool {
        self.inner.is_arrow()
    }

    #[getter]
    fn length(&self) -> Option<usize> {
        self.inner.length()
    }

    #[getter]
    fn source(&self) -> Option<String> {
        endpoint_name(self.inner.source())
    }

    #[getter]
    fn target(&self) -> Option<String> {
        endpoint_name(self.inner.target())
    }

    #[getter]
    fn external_representation(&self) -> Vec<usize> {
        self.inner.external_representation()
    }

    /// Concrete QPA generators in `WalkOfPath` order.
    #[getter]
    fn walk(&self, py: Python<'_>) -> PyResult<Option<Vec<Py<PyAny>>>> {
        self.inner
            .generator_walk()
            .map(|walk| {
                walk.into_iter()
                    .map(|generator| generator_to_python(py, generator))
                    .collect()
            })
            .transpose()
    }

    #[getter]
    fn source_object(&self) -> Option<PyQuiverVertex> {
        endpoint_to_python(self.inner.source())
    }

    #[getter]
    fn target_object(&self) -> Option<PyQuiverVertex> {
        endpoint_to_python(self.inner.target())
    }

    #[getter]
    fn generator_index(&self) -> PyResult<usize> {
        if self.inner.is_vertex() || self.inner.is_arrow() {
            self.inner
                .external_representation()
                .first()
                .copied()
                .ok_or_else(|| invalid("this path has no generator index"))
        } else {
            Err(invalid(
                "GeneratorIndex is defined only for vertices and arrows",
            ))
        }
    }

    fn multiply(&self, right: &Self) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.multiply(&right.inner).map_err(invalid)?,
        })
    }

    fn opposite(&self) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .quiver()
                .opposite_path(&self.inner)
                .map_err(invalid)?,
        })
    }

    fn walk_or_vertex(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        walk_of_path_or_vertex(&self.inner)
            .map_err(invalid)?
            .into_iter()
            .map(|generator| generator_to_python(py, generator))
            .collect()
    }

    fn reverse_path(&self) -> PyResult<Self> {
        Ok(Self {
            inner: reverse_path(&self.inner).map_err(invalid)?,
        })
    }

    fn include_in_product_quiver(&self, right: &Self, product: &PyQuiverProduct) -> PyResult<Self> {
        Ok(Self {
            inner: include_in_product_quiver([&self.inner, &right.inner], &product.inner)
                .map_err(invalid)?,
        })
    }

    fn project_from_product_quiver(
        &self,
        one_based_factor: usize,
        product: &PyQuiverProduct,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: project_from_product_quiver(one_based_factor, &self.inner, &product.inner)
                .map_err(invalid)?,
        })
    }

    fn __mul__(&self, right: &Self) -> PyResult<Self> {
        self.multiply(right)
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}
