//! Python wrappers for QPA's finite-type numerical degeneration order.

use std::sync::Arc;

use pyo3::{PyAny, pyfunction, pymethods, types::PyInt};
use qpa_gbnp::{
    ArMesh, ArOrbit, ArOrbitPosition, ArQuiverDescription, ArQuiverNumerical,
    CodimensionConjectureReport, DimensionVectorInput, NumericalModule, PredefinedArQuiver,
    PredefinedArQuiverData,
};

#[allow(clippy::wildcard_imports)]
// The Python error and integer-conversion conventions belong to the parent module.
use super::*;
use crate::field::{
    parse_python_integer_vector, wrap_python_integer, wrap_python_integer_matrix,
    wrap_python_integer_vector,
};

type PyOrbitDescription = Vec<(usize, Vec<(usize, usize)>)>;
type PyIntegerMatrix = Vec<Vec<Py<PyAny>>>;
type PyIntegerMatrixPair = (PyIntegerMatrix, PyIntegerMatrix);

fn positive_index(value: &Bound<'_, PyAny>) -> PyResult<usize> {
    let index = value
        .extract::<usize>()
        .map_err(|_| invalid("an indecomposable module number must be a positive integer"))?;
    index
        .checked_sub(1)
        .ok_or_else(|| invalid("an indecomposable module number must be positive"))
}

fn numerical_module(py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<NumericalModule> {
    if value.is_instance_of::<PyInt>() {
        return Ok(NumericalModule::Indecomposable(positive_index(value)?));
    }
    let values = value.extract::<Vec<Py<PyAny>>>().map_err(|_| {
        invalid("a module must be a positive indecomposable number or multiplicity vector")
    })?;
    Ok(NumericalModule::Multiplicities(
        parse_python_integer_vector(py, values)?,
    ))
}

fn dimension_vector_input(
    py: Python<'_>,
    value: &Bound<'_, PyAny>,
) -> PyResult<DimensionVectorInput> {
    if value.is_instance_of::<PyInt>() {
        return Ok(DimensionVectorInput::Indecomposable(positive_index(value)?));
    }
    let values = value.extract::<Vec<Py<PyAny>>>().map_err(|_| {
        invalid("the input must be a positive indecomposable number or dimension vector")
    })?;
    Ok(DimensionVectorInput::DimensionVector(
        parse_python_integer_vector(py, values)?,
    ))
}

fn mesh_description(description: &ArQuiverDescription) -> Option<Vec<Vec<usize>>> {
    let ArQuiverDescription::Meshes(meshes) = description else {
        return None;
    };
    Some(
        meshes
            .iter()
            .map(|mesh| {
                mesh.direct_predecessors
                    .iter()
                    .map(|position| position + 1)
                    .chain(std::iter::once(mesh.tau.map_or(0, |position| position + 1)))
                    .collect()
            })
            .collect(),
    )
}

fn orbit_description(description: &ArQuiverDescription) -> Option<PyOrbitDescription> {
    let ArQuiverDescription::Orbits(orbits) = description else {
        return None;
    };
    Some(
        orbits
            .iter()
            .map(|orbit| {
                (
                    orbit.length_minus_one,
                    orbit
                        .direct_predecessors
                        .iter()
                        .map(|position| (position.orbit + 1, position.shift))
                        .collect(),
                )
            })
            .collect(),
    )
}

fn selector(
    name: &str,
    number: Option<usize>,
    parameter: Option<usize>,
) -> PyResult<PredefinedArQuiver> {
    let fixed = match name {
        "D4 subspace" => Some(PredefinedArQuiver::D4Subspace),
        "E6 subspace" => Some(PredefinedArQuiver::E6Subspace),
        "A3 zero" => Some(PredefinedArQuiver::A3Zero),
        "A1,2 zero" => Some(PredefinedArQuiver::A12Zero),
        "A2,2 comm" => Some(PredefinedArQuiver::A22Commutative),
        "R nilp" => Some(PredefinedArQuiver::RNilpotent),
        _ => None,
    };
    if let Some(selector) = fixed {
        if number.is_none() && parameter.is_none() {
            return Ok(selector);
        }
        return Err(invalid("wrong number or type of arguments"));
    }
    match (name, number, parameter) {
        ("BG", Some(number), parameter) => {
            Ok(PredefinedArQuiver::BongartzGabriel { number, parameter })
        }
        ("BGo", Some(number), None) => Ok(PredefinedArQuiver::BongartzGabrielOrbits { number }),
        _ => Err(invalid("wrong number or type of arguments")),
    }
}

fn meshes_from_python(rows: Vec<Vec<usize>>) -> PyResult<Vec<ArMesh>> {
    rows.into_iter()
        .map(|mut row| {
            let tau = row
                .pop()
                .ok_or_else(|| invalid("an AR-quiver mesh row must be nonempty"))?;
            let direct_predecessors = row
                .into_iter()
                .map(|position| {
                    position
                        .checked_sub(1)
                        .ok_or_else(|| invalid("mesh predecessor numbers must be positive"))
                })
                .collect::<PyResult<Vec<_>>>()?;
            let tau = if tau == 0 { None } else { Some(tau - 1) };
            Ok(ArMesh::new(direct_predecessors, tau))
        })
        .collect()
}

fn orbits_from_python(orbits: PyOrbitDescription) -> PyResult<Vec<ArOrbit>> {
    orbits
        .into_iter()
        .map(|(length_minus_one, predecessors)| {
            let predecessors = predecessors
                .into_iter()
                .map(|(orbit, shift)| {
                    Ok(ArOrbitPosition::new(
                        orbit
                            .checked_sub(1)
                            .ok_or_else(|| invalid("orbit predecessor numbers must be positive"))?,
                        shift,
                    ))
                })
                .collect::<PyResult<Vec<_>>>()?;
            Ok(ArOrbit::new(length_minus_one, predecessors))
        })
        .collect()
}

/// QPA's numerical finite-type Auslander--Reiten quiver.
#[pyclass(name = "ArQuiverNumerical", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyArQuiverNumerical {
    inner: Arc<ArQuiverNumerical>,
}

impl PyArQuiverNumerical {
    fn from_inner(inner: ArQuiverNumerical) -> Self {
        Self {
            inner: Arc::new(inner),
        }
    }
}

#[pymethods]
impl PyArQuiverNumerical {
    #[new]
    fn new(
        number_of_indecomposables: usize,
        number_of_projectives: usize,
        mesh_description: Vec<Vec<usize>>,
    ) -> PyResult<Self> {
        Ok(Self::from_inner(
            ArQuiverNumerical::new(
                number_of_indecomposables,
                number_of_projectives,
                ArQuiverDescription::Meshes(meshes_from_python(mesh_description)?),
            )
            .map_err(invalid)?,
        ))
    }

    #[staticmethod]
    fn from_orbits(
        number_of_indecomposables: usize,
        number_of_projectives: usize,
        orbit_description: PyOrbitDescription,
    ) -> PyResult<Self> {
        Ok(Self::from_inner(
            ArQuiverNumerical::new(
                number_of_indecomposables,
                number_of_projectives,
                ArQuiverDescription::Orbits(orbits_from_python(orbit_description)?),
            )
            .map_err(invalid)?,
        ))
    }

    #[staticmethod]
    #[pyo3(signature = (name, number=None, parameter=None))]
    fn predefined(name: &str, number: Option<usize>, parameter: Option<usize>) -> PyResult<Self> {
        Ok(Self::from_inner(
            qpa_gbnp::predefined_ar_quiver(selector(name, number, parameter)?).map_err(invalid)?,
        ))
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_ar_quiver_numerical(&self) -> bool {
        true
    }

    #[getter]
    fn number_of_indecomposables(&self) -> usize {
        self.inner.number_of_indecomposables()
    }

    #[getter]
    fn number_of_projectives(&self) -> usize {
        self.inner.number_of_projectives()
    }

    #[getter]
    fn mesh_description(&self) -> Option<Vec<Vec<usize>>> {
        mesh_description(self.inner.description())
    }

    #[getter]
    fn orbit_description(&self) -> Option<PyOrbitDescription> {
        orbit_description(self.inner.description())
    }

    #[getter]
    fn dim_hom_matrix(&self, py: Python<'_>) -> PyResult<Vec<Vec<Py<PyAny>>>> {
        wrap_python_integer_matrix(py, self.inner.dim_hom_matrix().iter().cloned())
    }

    #[getter]
    fn simples(&self) -> Vec<usize> {
        self.inner
            .simples()
            .iter()
            .map(|position| position + 1)
            .collect()
    }

    fn dimension_vector(
        &self,
        py: Python<'_>,
        module: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Py<PyAny>>> {
        let module = numerical_module(py, module)?;
        wrap_python_integer_vector(py, &self.inner.dimension_vector(&module).map_err(invalid)?)
    }

    fn modules_of_dimension_vector(
        &self,
        py: Python<'_>,
        which: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Vec<Py<PyAny>>>> {
        let which = dimension_vector_input(py, which)?;
        wrap_python_integer_matrix(
            py,
            self.inner
                .modules_of_dimension_vector(&which)
                .map_err(invalid)?
                .into_iter(),
        )
    }

    fn dim_hom(
        &self,
        py: Python<'_>,
        source: &Bound<'_, PyAny>,
        target: &Bound<'_, PyAny>,
    ) -> PyResult<Py<PyAny>> {
        let source = numerical_module(py, source)?;
        let target = numerical_module(py, target)?;
        wrap_python_integer(py, &self.inner.dim_hom(&source, &target).map_err(invalid)?)
    }

    fn dim_end(&self, py: Python<'_>, module: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let module = numerical_module(py, module)?;
        wrap_python_integer(py, &self.inner.dim_end(&module).map_err(invalid)?)
    }

    fn orbit_dimension(&self, py: Python<'_>, module: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let module = numerical_module(py, module)?;
        wrap_python_integer(py, &self.inner.orbit_dimension(&module).map_err(invalid)?)
    }

    fn orbit_codimension(
        &self,
        py: Python<'_>,
        source: &Bound<'_, PyAny>,
        target: &Bound<'_, PyAny>,
    ) -> PyResult<Py<PyAny>> {
        let source = numerical_module(py, source)?;
        let target = numerical_module(py, target)?;
        wrap_python_integer(
            py,
            &self
                .inner
                .orbit_codimension(&source, &target)
                .map_err(invalid)?,
        )
    }

    fn degeneration_order_leq(
        &self,
        py: Python<'_>,
        source: &Bound<'_, PyAny>,
        target: &Bound<'_, PyAny>,
    ) -> PyResult<bool> {
        self.inner
            .degeneration_order_leq(
                &numerical_module(py, source)?,
                &numerical_module(py, target)?,
            )
            .map_err(invalid)
    }

    fn degeneration_order_leq_unchecked(
        &self,
        py: Python<'_>,
        source: &Bound<'_, PyAny>,
        target: &Bound<'_, PyAny>,
    ) -> PyResult<bool> {
        self.inner
            .degeneration_order_leq_unchecked(
                &numerical_module(py, source)?,
                &numerical_module(py, target)?,
            )
            .map_err(invalid)
    }

    fn degeneration_predecessors(
        &self,
        py: Python<'_>,
        module: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Vec<Py<PyAny>>>> {
        wrap_python_integer_matrix(
            py,
            self.inner
                .degeneration_predecessors(&numerical_module(py, module)?)
                .map_err(invalid)?
                .into_iter(),
        )
    }

    fn direct_degeneration_predecessors(
        &self,
        py: Python<'_>,
        module: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Vec<Py<PyAny>>>> {
        wrap_python_integer_matrix(
            py,
            self.inner
                .direct_degeneration_predecessors(&numerical_module(py, module)?)
                .map_err(invalid)?
                .into_iter(),
        )
    }

    fn degeneration_predecessors_with_direct(
        &self,
        py: Python<'_>,
        module: &Bound<'_, PyAny>,
    ) -> PyResult<PyIntegerMatrixPair> {
        let (all, direct) = self
            .inner
            .degeneration_predecessors_with_direct(&numerical_module(py, module)?)
            .map_err(invalid)?;
        Ok((
            wrap_python_integer_matrix(py, all.into_iter())?,
            wrap_python_integer_matrix(py, direct.into_iter())?,
        ))
    }

    fn degeneration_successors(
        &self,
        py: Python<'_>,
        module: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Vec<Py<PyAny>>>> {
        wrap_python_integer_matrix(
            py,
            self.inner
                .degeneration_successors(&numerical_module(py, module)?)
                .map_err(invalid)?
                .into_iter(),
        )
    }

    fn direct_degeneration_successors(
        &self,
        py: Python<'_>,
        module: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Vec<Py<PyAny>>>> {
        wrap_python_integer_matrix(
            py,
            self.inner
                .direct_degeneration_successors(&numerical_module(py, module)?)
                .map_err(invalid)?
                .into_iter(),
        )
    }

    fn degeneration_successors_with_direct(
        &self,
        py: Python<'_>,
        module: &Bound<'_, PyAny>,
    ) -> PyResult<PyIntegerMatrixPair> {
        let (all, direct) = self
            .inner
            .degeneration_successors_with_direct(&numerical_module(py, module)?)
            .map_err(invalid)?;
        Ok((
            wrap_python_integer_matrix(py, all.into_iter())?,
            wrap_python_integer_matrix(py, direct.into_iter())?,
        ))
    }

    fn test_codimension_conjecture(
        &self,
        print_progress: bool,
        break_on_first: bool,
    ) -> PyResult<PyCodimensionConjectureReport> {
        Ok(PyCodimensionConjectureReport {
            inner: self
                .inner
                .test_codimension_conjecture(print_progress, break_on_first)
                .map_err(invalid)?,
        })
    }

    fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.inner, &other.inner)
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }
}

/// Raw QPA predefined AR-quiver data before numerical construction.
#[pyclass(name = "PredefinedArQuiverData", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPredefinedArQuiverData {
    inner: PredefinedArQuiverData,
}

#[pymethods]
impl PyPredefinedArQuiverData {
    #[getter]
    fn number_of_indecomposables(&self) -> usize {
        self.inner.number_of_indecomposables
    }

    #[getter]
    fn number_of_projectives(&self) -> usize {
        self.inner.number_of_projectives
    }

    #[getter]
    fn mesh_description(&self) -> Option<Vec<Vec<usize>>> {
        mesh_description(&self.inner.description)
    }

    #[getter]
    fn orbit_description(&self) -> Option<PyOrbitDescription> {
        orbit_description(&self.inner.description)
    }

    fn ar_quiver(&self) -> PyResult<PyArQuiverNumerical> {
        Ok(PyArQuiverNumerical::from_inner(
            self.inner.clone().into_ar_quiver().map_err(invalid)?,
        ))
    }
}

/// Structured result of QPA's codimension-conjecture diagnostic loop.
#[pyclass(name = "CodimensionConjectureReport", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyCodimensionConjectureReport {
    inner: CodimensionConjectureReport,
}

#[pymethods]
impl PyCodimensionConjectureReport {
    #[getter]
    fn counterexamples(&self, py: Python<'_>) -> PyResult<Vec<(Vec<Py<PyAny>>, usize)>> {
        self.inner
            .counterexamples
            .iter()
            .map(|(vector, indecomposable)| {
                Ok((wrap_python_integer_vector(py, vector)?, indecomposable + 1))
            })
            .collect()
    }

    #[getter]
    fn predecessor_count(&self) -> usize {
        self.inner.predecessor_count
    }

    #[getter]
    fn direct_predecessor_count(&self) -> usize {
        self.inner.direct_predecessor_count
    }

    #[getter]
    fn stopped_early(&self) -> bool {
        self.inner.stopped_early
    }
}

/// QPA's raw `PredefARQuivers` dispatch.
#[pyfunction]
#[pyo3(name = "predefined_ar_quiver_data", signature = (name, number=None, parameter=None))]
pub(crate) fn py_predefined_ar_quiver_data(
    name: &str,
    number: Option<usize>,
    parameter: Option<usize>,
) -> PyResult<PyPredefinedArQuiverData> {
    Ok(PyPredefinedArQuiverData {
        inner: qpa_gbnp::predefined_ar_quiver_data(selector(name, number, parameter)?)
            .map_err(invalid)?,
    })
}

/// QPA's stable `PredefARQuivers("what")` catalog output.
#[pyfunction]
#[pyo3(name = "predefined_ar_quiver_catalog")]
pub(crate) fn py_predefined_ar_quiver_catalog() -> &'static str {
    qpa_gbnp::PREDEFINED_AR_QUIVER_CATALOG
}

/// Return QPA's `PrintMultiplicityVector` rendering instead of writing stdout.
#[pyfunction]
#[pyo3(name = "format_multiplicity_vector")]
pub(crate) fn py_format_multiplicity_vector(
    py: Python<'_>,
    vector: Vec<Py<PyAny>>,
) -> PyResult<String> {
    Ok(qpa_gbnp::format_multiplicity_vector(
        &parse_python_integer_vector(py, vector)?,
    ))
}

/// Return QPA's `PrintMultiplicityVectors` rendering instead of writing stdout.
#[pyfunction]
#[pyo3(name = "format_multiplicity_vectors")]
pub(crate) fn py_format_multiplicity_vectors(
    py: Python<'_>,
    vectors: Vec<Vec<Py<PyAny>>>,
) -> PyResult<String> {
    let vectors = vectors
        .into_iter()
        .map(|vector| parse_python_integer_vector(py, vector))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(qpa_gbnp::format_multiplicity_vectors(&vectors))
}
