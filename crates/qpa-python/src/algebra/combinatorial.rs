//! Python wrappers for QPA combinatorial maps of gentle algebras.

use pyo3::pymethods;
use qpa_core::{CombinatorialMap, Permutation};

use super::conversion::QuiverAlgebraKind;
#[allow(clippy::wildcard_imports)]
// The gentle-algebra runtime union belongs to the parent binding module.
use super::*;
use crate::ScalarKind;

/// QPA's combinatorial map on a finite half-edge set.
#[pyclass(name = "CombinatorialMap", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyCombinatorialMap {
    inner: CombinatorialMap,
}

type PyCutNonSeparatingResult = (PyCombinatorialMap, Vec<usize>, Vec<usize>, Vec<usize>);

impl PyCombinatorialMap {
    fn from_inner(inner: CombinatorialMap) -> Self {
        Self { inner }
    }
}

pub(super) fn maximal_paths_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Vec<PyPath>> {
    Ok(qpa_gbnp::maximal_paths_of_gentle_algebra(algebra)
        .map_err(invalid)?
        .iter()
        .cloned()
        .map(|inner| PyPath { inner })
        .collect())
}

pub(super) fn combinatorial_map_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<PyCombinatorialMap> {
    Ok(PyCombinatorialMap::from_inner(
        qpa_gbnp::combinatorial_map_of_gentle_algebra(algebra)
            .map_err(invalid)?
            .clone(),
    ))
}

pub(super) fn derived_equivalent(
    left: &QuiverAlgebraKind,
    right: &QuiverAlgebraKind,
) -> PyResult<bool> {
    match (left, right) {
        (QuiverAlgebraKind::Rational(left), QuiverAlgebraKind::Rational(right)) => {
            qpa_gbnp::are_derived_equivalent(left, right).map_err(invalid)
        }
        (QuiverAlgebraKind::Rational(left), QuiverAlgebraKind::Binary(right)) => {
            qpa_gbnp::are_derived_equivalent(left, right).map_err(invalid)
        }
        (QuiverAlgebraKind::Rational(left), QuiverAlgebraKind::Prime(right)) => {
            qpa_gbnp::are_derived_equivalent(left, right).map_err(invalid)
        }
        (QuiverAlgebraKind::Binary(left), QuiverAlgebraKind::Rational(right)) => {
            qpa_gbnp::are_derived_equivalent(left, right).map_err(invalid)
        }
        (QuiverAlgebraKind::Binary(left), QuiverAlgebraKind::Binary(right)) => {
            qpa_gbnp::are_derived_equivalent(left, right).map_err(invalid)
        }
        (QuiverAlgebraKind::Binary(left), QuiverAlgebraKind::Prime(right)) => {
            qpa_gbnp::are_derived_equivalent(left, right).map_err(invalid)
        }
        (QuiverAlgebraKind::Prime(left), QuiverAlgebraKind::Rational(right)) => {
            qpa_gbnp::are_derived_equivalent(left, right).map_err(invalid)
        }
        (QuiverAlgebraKind::Prime(left), QuiverAlgebraKind::Binary(right)) => {
            qpa_gbnp::are_derived_equivalent(left, right).map_err(invalid)
        }
        (QuiverAlgebraKind::Prime(left), QuiverAlgebraKind::Prime(right)) => {
            qpa_gbnp::are_derived_equivalent(left, right).map_err(invalid)
        }
    }
}

fn permutation_images(permutation: &Permutation, size: usize) -> PyResult<Vec<usize>> {
    permutation.images_on(size).map_err(invalid)
}

#[pymethods]
impl PyCombinatorialMap {
    #[new]
    fn new(
        size: usize,
        ordering: Vec<usize>,
        pairing: Vec<usize>,
        marked_edges: Vec<usize>,
    ) -> PyResult<Self> {
        Ok(Self::from_inner(
            CombinatorialMap::new(
                size,
                Permutation::from_images(ordering).map_err(invalid)?,
                Permutation::from_images(pairing).map_err(invalid)?,
                marked_edges,
            )
            .map_err(invalid)?,
        ))
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_combinatorial_map(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn size(&self) -> usize {
        self.inner.size()
    }

    #[getter]
    fn ordering(&self) -> PyResult<Vec<usize>> {
        permutation_images(self.inner.ordering(), self.inner.size())
    }

    #[getter]
    fn pairing(&self) -> PyResult<Vec<usize>> {
        permutation_images(self.inner.pairing(), self.inner.size())
    }

    #[getter]
    fn marked_edges(&self) -> Vec<usize> {
        self.inner.marked_edges().to_vec()
    }

    #[getter]
    fn faces(&self) -> Vec<Vec<usize>> {
        self.inner.faces().to_vec()
    }

    #[getter]
    fn genus(&self) -> PyScalar {
        PyScalar {
            value: ScalarKind::Rational(self.inner.genus().clone()),
        }
    }

    fn dual(&self) -> Self {
        Self::from_inner(self.inner.dual().clone())
    }

    fn is_identical(&self, other: &Self) -> bool {
        self.inner.is_identical(&other.inner)
    }

    fn remove_edge(&self, edge: Vec<usize>) -> PyResult<Self> {
        Ok(Self::from_inner(
            self.inner.remove_edge(&edge).map_err(invalid)?,
        ))
    }

    #[getter]
    fn marked_boundaries(&self) -> Vec<(Vec<usize>, usize)> {
        self.inner.marked_boundaries().to_vec()
    }

    #[getter]
    fn boundary_curves(&self) -> PyResult<Vec<(Vec<usize>, usize)>> {
        Ok(self.inner.boundary_curves().map_err(invalid)?.to_vec())
    }

    fn winding_number(&self, curve: Vec<usize>) -> PyResult<i64> {
        self.inner.winding_number(&curve).map_err(invalid)
    }

    fn depth_search(&self, half_edge: usize) -> PyResult<Vec<Vec<usize>>> {
        self.inner.depth_search(half_edge).map_err(invalid)
    }

    fn width_search(&self, half_edge: usize) -> PyResult<Vec<Vec<usize>>> {
        self.inner.width_search(half_edge).map_err(invalid)
    }

    #[getter]
    fn non_separating_curve(&self) -> PyResult<Vec<usize>> {
        Ok(self.inner.non_separating_curve().map_err(invalid)?.to_vec())
    }

    fn cut_non_separating(
        &self,
        curve: Vec<usize>,
        index: Vec<usize>,
    ) -> PyResult<PyCutNonSeparatingResult> {
        let (map, first, second, index) = self
            .inner
            .cut_non_separating(&curve, &index)
            .map_err(invalid)?
            .into_parts();
        Ok((Self::from_inner(map), first, second, index))
    }

    fn join_curve(
        &self,
        first_boundary: Vec<usize>,
        second_boundary: Vec<usize>,
        index: Vec<usize>,
    ) -> PyResult<Vec<usize>> {
        self.inner
            .join_curve(&first_boundary, &second_boundary, &index)
            .map_err(invalid)
    }

    fn cut_join_curve(
        &self,
        first_boundary: Vec<usize>,
        second_boundary: Vec<usize>,
        curve: Vec<usize>,
        index: Vec<usize>,
    ) -> PyResult<(Self, Vec<usize>)> {
        let (map, index) = self
            .inner
            .cut_join_curve(&first_boundary, &second_boundary, &curve, &index)
            .map_err(invalid)?
            .into_parts();
        Ok((Self::from_inner(map), index))
    }

    #[getter]
    fn homology_basis(&self) -> PyResult<Vec<(Vec<usize>, Vec<usize>)>> {
        Ok(self.inner.homology_basis().map_err(invalid)?.to_vec())
    }

    fn __repr__(&self) -> String {
        format!("<combinatorial map on {} half-edges>", self.inner.size())
    }
}
