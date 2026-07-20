//! Python dispatch for QPA's Brauer-configuration algebra constructor.

use pyo3::{PyAny, pyfunction, pymethods};
use qpa_core::QuiverAlgebra;
use qpa_gbnp::{
    BrauerConfiguration, BrauerConfigurationAlgebraResult, BrauerConfigurationPolygon,
    BrauerConfigurationValidationError, BrauerConfigurationVertex,
};

#[allow(clippy::wildcard_imports)]
// Runtime fields and the quiver-algebra wrapper belong to the parent module.
use super::*;
use crate::field::{parse_python_integer, wrap_python_integer};

type PyBrauerVertex = (String, Py<PyAny>);
type PyBrauerPolygon = (String, Vec<String>);

/// Typed Python form of QPA's vertices, polygons, and orientation rows.
#[pyclass(name = "BrauerConfiguration", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyBrauerConfiguration {
    inner: BrauerConfiguration,
}

#[pymethods]
impl PyBrauerConfiguration {
    #[new]
    fn new(
        py: Python<'_>,
        vertices: Vec<PyBrauerVertex>,
        polygons: Vec<PyBrauerPolygon>,
        orientations: Vec<Vec<String>>,
    ) -> PyResult<Self> {
        let vertices = vertices
            .into_iter()
            .map(|(name, multiplicity)| {
                Ok(BrauerConfigurationVertex {
                    name,
                    multiplicity: parse_python_integer(multiplicity.bind(py))?,
                })
            })
            .collect::<PyResult<Vec<_>>>()?;
        let polygons = polygons
            .into_iter()
            .map(|(name, vertices)| BrauerConfigurationPolygon { name, vertices })
            .collect();
        Ok(Self {
            inner: BrauerConfiguration {
                vertices,
                polygons,
                orientations,
            },
        })
    }

    #[getter]
    fn vertices(&self, py: Python<'_>) -> PyResult<Vec<PyBrauerVertex>> {
        self.inner
            .vertices
            .iter()
            .map(|vertex| {
                Ok((
                    vertex.name.clone(),
                    wrap_python_integer(py, &vertex.multiplicity)?,
                ))
            })
            .collect()
    }

    #[getter]
    fn polygons(&self) -> Vec<PyBrauerPolygon> {
        self.inner
            .polygons
            .iter()
            .map(|polygon| (polygon.name.clone(), polygon.vertices.clone()))
            .collect()
    }

    #[getter]
    fn orientations(&self) -> Vec<Vec<String>> {
        self.inner.orientations.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "<Brauer configuration with {} vertices and {} polygons>",
            self.inner.vertices.len(),
            self.inner.polygons.len()
        )
    }
}

/// A valid QPA Brauer-configuration algebra with its constructor tag.
#[pyclass(name = "BrauerConfigurationAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyBrauerConfigurationAlgebra {
    algebra: PyQuiverAlgebra,
}

#[pymethods]
impl PyBrauerConfigurationAlgebra {
    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        self.algebra.clone()
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_symmetric_algebra(&self) -> bool {
        true
    }

    #[allow(clippy::unused_self)]
    fn __repr__(&self) -> &'static str {
        "<QPA Brauer-configuration algebra>"
    }
}

/// QPA's algebra-or-`fail` Brauer-configuration result.
#[pyclass(name = "BrauerConfigurationAlgebraResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyBrauerConfigurationAlgebraResult {
    algebra: Option<PyBrauerConfigurationAlgebra>,
    diagnostic: Option<&'static str>,
    validation_error: Option<&'static str>,
}

#[pymethods]
impl PyBrauerConfigurationAlgebraResult {
    #[getter]
    fn algebra(&self) -> Option<PyBrauerConfigurationAlgebra> {
        self.algebra.clone()
    }

    #[getter]
    fn diagnostic(&self) -> Option<&'static str> {
        self.diagnostic
    }

    #[getter]
    fn validation_error(&self) -> Option<&'static str> {
        self.validation_error
    }

    #[getter]
    fn is_valid(&self) -> bool {
        self.algebra.is_some()
    }

    fn __repr__(&self) -> &'static str {
        if self.is_valid() {
            "<valid QPA Brauer-configuration result>"
        } else {
            "<invalid QPA Brauer-configuration result>"
        }
    }
}

fn validation_error_name(error: &BrauerConfigurationValidationError) -> &'static str {
    match error {
        BrauerConfigurationValidationError::NoVertices => "no_vertices",
        BrauerConfigurationValidationError::NonPositiveMultiplicity { .. } => {
            "non_positive_multiplicity"
        }
        BrauerConfigurationValidationError::DuplicateVertexName { .. } => "duplicate_vertex_name",
        BrauerConfigurationValidationError::PolygonWithoutVertex { .. } => "polygon_without_vertex",
        BrauerConfigurationValidationError::PolygonNameMatchesVertex { .. } => {
            "polygon_name_matches_vertex"
        }
        BrauerConfigurationValidationError::PolygonContainsUnknownVertex { .. } => {
            "polygon_contains_unknown_vertex"
        }
        BrauerConfigurationValidationError::DuplicatePolygonName { .. } => "duplicate_polygon_name",
        BrauerConfigurationValidationError::OrientationCountMismatch { .. } => {
            "orientation_count_mismatch"
        }
        BrauerConfigurationValidationError::EmptyOrientation { .. } => "empty_orientation",
        BrauerConfigurationValidationError::OrientationContainsUnknownPolygon { .. } => {
            "orientation_contains_unknown_polygon"
        }
        BrauerConfigurationValidationError::OrientationContainmentMismatch { .. } => {
            "orientation_containment_mismatch"
        }
    }
}

fn brauer_for<F: runtime::PythonField>(
    field: F,
    configuration: &mut BrauerConfiguration,
) -> PyResult<PyBrauerConfigurationAlgebraResult>
where
    super::conversion::QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    let result = qpa_gbnp::brauer_configuration_algebra(field, configuration).map_err(invalid)?;
    Ok(match result {
        BrauerConfigurationAlgebraResult::Algebra(algebra) => PyBrauerConfigurationAlgebraResult {
            algebra: Some(PyBrauerConfigurationAlgebra {
                algebra: super::conversion::wrap_quiver_algebra(algebra.into_algebra()),
            }),
            diagnostic: None,
            validation_error: None,
        },
        BrauerConfigurationAlgebraResult::InvalidConfiguration(error) => {
            PyBrauerConfigurationAlgebraResult {
                algebra: None,
                diagnostic: Some(error.diagnostic()),
                validation_error: Some(validation_error_name(&error)),
            }
        }
    })
}

/// Run QPA's Brauer-configuration constructor, including its `fail` branches.
#[pyfunction]
#[pyo3(name = "brauer_configuration_algebra")]
pub(crate) fn py_brauer_configuration_algebra(
    field: &PyField,
    mut configuration: PyRefMut<'_, PyBrauerConfiguration>,
) -> PyResult<PyBrauerConfigurationAlgebraResult> {
    match field.kind {
        FieldKind::Rational => brauer_for(Q, &mut configuration.inner),
        FieldKind::Binary => brauer_for(Z2, &mut configuration.inner),
        FieldKind::Prime(characteristic) => {
            brauer_for(Zp::new(characteristic), &mut configuration.inner)
        }
    }
}
