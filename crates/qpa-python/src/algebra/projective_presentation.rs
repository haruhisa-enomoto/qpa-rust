//! Python wrapper for QPA's legacy `present2` representation.

use numerica::domains::{
    finite_field::{Z2, Zp},
    rational::RationalField,
};
use pyo3::{Py, PyResult, Python, prelude::*};
use qpa_core::{PathAlgebraElement, ProjectivePresentation};

use super::{AlgebraElementKind, AlgebraKind, PyAlgebraElement, PyPathAlgebra};
use crate::{invalid, quiver::PyQuiverVertex};

#[derive(Clone)]
enum ProjectivePresentationKind {
    Rational(ProjectivePresentation<RationalField>),
    Binary(ProjectivePresentation<Z2>),
    Prime(ProjectivePresentation<Zp>),
}

/// QPA's declared projective-presentation attribute object.
#[pyclass(name = "ProjectivePresentation", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyProjectivePresentation {
    kind: ProjectivePresentationKind,
}

fn rational_matrix(
    py: Python<'_>,
    matrix: &[Vec<Py<PyAlgebraElement>>],
) -> PyResult<Vec<Vec<PathAlgebraElement<RationalField>>>> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|element| match &element.borrow(py).kind {
                    AlgebraElementKind::Rational(value) => Ok(value.clone()),
                    _ => Err(invalid(
                        "matrix element belongs to another coefficient field",
                    )),
                })
                .collect()
        })
        .collect()
}

fn binary_matrix(
    py: Python<'_>,
    matrix: &[Vec<Py<PyAlgebraElement>>],
) -> PyResult<Vec<Vec<PathAlgebraElement<Z2>>>> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|element| match &element.borrow(py).kind {
                    AlgebraElementKind::Binary(value) => Ok(value.clone()),
                    _ => Err(invalid(
                        "matrix element belongs to another coefficient field",
                    )),
                })
                .collect()
        })
        .collect()
}

fn prime_matrix(
    py: Python<'_>,
    matrix: &[Vec<Py<PyAlgebraElement>>],
) -> PyResult<Vec<Vec<PathAlgebraElement<Zp>>>> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|element| match &element.borrow(py).kind {
                    AlgebraElementKind::Prime(value) => Ok(value.clone()),
                    _ => Err(invalid(
                        "matrix element belongs to another coefficient field",
                    )),
                })
                .collect()
        })
        .collect()
}

fn wrap_rational_matrix(
    matrix: &[Vec<PathAlgebraElement<RationalField>>],
) -> Vec<Vec<PyAlgebraElement>> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .cloned()
                .map(|value| PyAlgebraElement {
                    kind: AlgebraElementKind::Rational(value),
                })
                .collect()
        })
        .collect()
}

fn wrap_binary_matrix(matrix: &[Vec<PathAlgebraElement<Z2>>]) -> Vec<Vec<PyAlgebraElement>> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .cloned()
                .map(|value| PyAlgebraElement {
                    kind: AlgebraElementKind::Binary(value),
                })
                .collect()
        })
        .collect()
}

fn wrap_prime_matrix(matrix: &[Vec<PathAlgebraElement<Zp>>]) -> Vec<Vec<PyAlgebraElement>> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .cloned()
                .map(|value| PyAlgebraElement {
                    kind: AlgebraElementKind::Prime(value),
                })
                .collect()
        })
        .collect()
}

#[pymethods]
impl PyProjectivePresentation {
    #[new]
    fn new(
        py: Python<'_>,
        parent_algebra: &PyPathAlgebra,
        matrix: Vec<Vec<Py<PyAlgebraElement>>>,
    ) -> PyResult<Self> {
        let kind = match &parent_algebra.kind {
            AlgebraKind::Rational(algebra) => ProjectivePresentationKind::Rational(
                ProjectivePresentation::from_matrix(algebra, &rational_matrix(py, &matrix)?)
                    .map_err(invalid)?,
            ),
            AlgebraKind::Binary(algebra) => ProjectivePresentationKind::Binary(
                ProjectivePresentation::from_matrix(algebra, &binary_matrix(py, &matrix)?)
                    .map_err(invalid)?,
            ),
            AlgebraKind::Prime(algebra) => ProjectivePresentationKind::Prime(
                ProjectivePresentation::from_matrix(algebra, &prime_matrix(py, &matrix)?)
                    .map_err(invalid)?,
            ),
        };
        Ok(Self { kind })
    }

    #[staticmethod]
    #[pyo3(signature = (parent_algebra, matrix, matrix_num_rows, matrix_num_cols, p0_vertex_list, p1_vertex_list=None))]
    #[allow(clippy::too_many_arguments)]
    fn from_attributes(
        py: Python<'_>,
        parent_algebra: &PyPathAlgebra,
        matrix: Vec<Vec<Py<PyAlgebraElement>>>,
        matrix_num_rows: usize,
        matrix_num_cols: usize,
        p0_vertex_list: Vec<Py<PyQuiverVertex>>,
        p1_vertex_list: Option<Vec<Py<PyQuiverVertex>>>,
    ) -> PyResult<Self> {
        let p0 = p0_vertex_list
            .iter()
            .map(|vertex| vertex.borrow(py).inner.clone())
            .collect::<Vec<_>>();
        let p1 = p1_vertex_list.map(|vertices| {
            vertices
                .iter()
                .map(|vertex| vertex.borrow(py).inner.clone())
                .collect::<Vec<_>>()
        });
        let kind = match &parent_algebra.kind {
            AlgebraKind::Rational(algebra) => {
                ProjectivePresentationKind::Rational(ProjectivePresentation::from_attributes(
                    algebra.clone(),
                    rational_matrix(py, &matrix)?,
                    matrix_num_rows,
                    matrix_num_cols,
                    p0,
                    p1,
                ))
            }
            AlgebraKind::Binary(algebra) => {
                ProjectivePresentationKind::Binary(ProjectivePresentation::from_attributes(
                    algebra.clone(),
                    binary_matrix(py, &matrix)?,
                    matrix_num_rows,
                    matrix_num_cols,
                    p0,
                    p1,
                ))
            }
            AlgebraKind::Prime(algebra) => {
                ProjectivePresentationKind::Prime(ProjectivePresentation::from_attributes(
                    algebra.clone(),
                    prime_matrix(py, &matrix)?,
                    matrix_num_rows,
                    matrix_num_cols,
                    p0,
                    p1,
                ))
            }
        };
        Ok(Self { kind })
    }

    #[classattr]
    #[pyo3(name = "family_name")]
    const FAMILY_NAME: &'static str = "ProjectivePresentationFamily";

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_projective_presentation(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn parent_algebra(&self) -> PyPathAlgebra {
        let kind = match &self.kind {
            ProjectivePresentationKind::Rational(value) => {
                AlgebraKind::Rational(value.parent_algebra().clone())
            }
            ProjectivePresentationKind::Binary(value) => {
                AlgebraKind::Binary(value.parent_algebra().clone())
            }
            ProjectivePresentationKind::Prime(value) => {
                AlgebraKind::Prime(value.parent_algebra().clone())
            }
        };
        PyPathAlgebra { kind }
    }

    #[getter]
    fn matrix_representation(&self) -> Vec<Vec<PyAlgebraElement>> {
        match &self.kind {
            ProjectivePresentationKind::Rational(value) => {
                wrap_rational_matrix(value.matrix_representation())
            }
            ProjectivePresentationKind::Binary(value) => {
                wrap_binary_matrix(value.matrix_representation())
            }
            ProjectivePresentationKind::Prime(value) => {
                wrap_prime_matrix(value.matrix_representation())
            }
        }
    }

    #[getter]
    fn matrix_num_rows(&self) -> usize {
        match &self.kind {
            ProjectivePresentationKind::Rational(value) => value.matrix_num_rows(),
            ProjectivePresentationKind::Binary(value) => value.matrix_num_rows(),
            ProjectivePresentationKind::Prime(value) => value.matrix_num_rows(),
        }
    }

    #[getter]
    fn matrix_num_cols(&self) -> usize {
        match &self.kind {
            ProjectivePresentationKind::Rational(value) => value.matrix_num_cols(),
            ProjectivePresentationKind::Binary(value) => value.matrix_num_cols(),
            ProjectivePresentationKind::Prime(value) => value.matrix_num_cols(),
        }
    }

    #[getter]
    fn p0_vertex_list(&self) -> Vec<PyQuiverVertex> {
        match &self.kind {
            ProjectivePresentationKind::Rational(value) => value.p0_vertex_list(),
            ProjectivePresentationKind::Binary(value) => value.p0_vertex_list(),
            ProjectivePresentationKind::Prime(value) => value.p0_vertex_list(),
        }
        .iter()
        .cloned()
        .map(|inner| PyQuiverVertex { inner })
        .collect()
    }

    #[getter]
    fn p1_vertex_list(&self) -> Option<Vec<PyQuiverVertex>> {
        match &self.kind {
            ProjectivePresentationKind::Rational(value) => value.p1_vertex_list(),
            ProjectivePresentationKind::Binary(value) => value.p1_vertex_list(),
            ProjectivePresentationKind::Prime(value) => value.p1_vertex_list(),
        }
        .map(|vertices| {
            vertices
                .iter()
                .cloned()
                .map(|inner| PyQuiverVertex { inner })
                .collect()
        })
    }

    #[getter]
    fn has_p1_vertex_list(&self) -> bool {
        match &self.kind {
            ProjectivePresentationKind::Rational(value) => value.p1_vertex_list().is_some(),
            ProjectivePresentationKind::Binary(value) => value.p1_vertex_list().is_some(),
            ProjectivePresentationKind::Prime(value) => value.p1_vertex_list().is_some(),
        }
    }
}
