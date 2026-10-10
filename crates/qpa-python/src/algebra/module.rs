//! Matrix modules, module elements, and module homomorphisms.

mod generic_hom;
mod homomorphism;
mod remaining;
mod results;
mod submodule;
mod values;

pub(crate) use generic_hom::{
    PyGenericModuleEndomorphismAlgebra, PyGenericModuleHomSpace, PyGenericModuleHomomorphism,
};
pub(super) use generic_hom::{generic_end, generic_hom};
pub(crate) use remaining::{
    PyBasedRightAlgebraModule, PyBongartzTestDecision, PyFiniteTypeDecision, PyModuleTensorProduct,
    py_complement_in_full_row_space,
};

use pyo3::pymethods;
use qpa_core::{ModuleDimensionSeries, PathAlgebraModuleMatrix, PathAlgebraVector, VertexIndex};

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

fn batch_ext_dimensions<F: runtime::PythonField>(
    module: &qpa_core::PathAlgebraModule<F>,
    targets: Vec<MatrixModuleKind>,
) -> PyResult<Vec<usize>>
where
    MatrixModuleKind: runtime::Extract<qpa_core::PathAlgebraModule<F>>,
{
    let targets: Vec<qpa_core::PathAlgebraModule<F>> = runtime::extract_all(targets)?;
    if targets.is_empty() {
        return Ok(Vec::new());
    }
    if targets
        .iter()
        .any(|target| target.right_acting_algebra() != module.right_acting_algebra())
    {
        return Err(invalid("Ext modules belong to different algebras"));
    }
    let prepared = qpa_gbnp::ExtDimensionSource::new(module).map_err(invalid)?;
    targets
        .iter()
        .map(|target| prepared.dimension(target).map_err(invalid))
        .collect()
}

#[derive(Clone)]
enum ProjectivePathPresentationKind {
    Rational(qpa_gbnp::ProjectivePathAlgebraPresentation<RationalField>),
    Binary(qpa_gbnp::ProjectivePathAlgebraPresentation<Z2>),
    Prime(qpa_gbnp::ProjectivePathAlgebraPresentation<Zp>),
    Extension(qpa_gbnp::ProjectivePathAlgebraPresentation<GaloisField>),
}

/// QPA's list-valued `ProjectivePathAlgebraPresentation` result.
#[pyclass(
    name = "ProjectivePathAlgebraPresentation",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PyProjectivePathAlgebraPresentation {
    kind: ProjectivePathPresentationKind,
}

#[derive(Clone)]
enum GszProjectiveResolutionKind {
    Rational(qpa_core::ProjectiveResolutionFpPathAlgebraModule<RationalField>),
    Binary(qpa_core::ProjectiveResolutionFpPathAlgebraModule<Z2>),
    Prime(qpa_core::ProjectiveResolutionFpPathAlgebraModule<Zp>),
    Extension(qpa_core::ProjectiveResolutionFpPathAlgebraModule<GaloisField>),
}

/// QPA's GSZ `ProjectiveResolutionFpPathAlgebraModule` attribute object.
#[pyclass(name = "GszProjectiveResolution", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyGszProjectiveResolution {
    kind: GszProjectiveResolutionKind,
}

/// Exact vertexwise homology of a projective resolution modulo the arrow
/// ideal.
#[pyclass(name = "SemisimpleResolutionHomology", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PySemisimpleResolutionHomology {
    projective_vertex_multiplicities: Vec<Vec<usize>>,
    differential_vertex_ranks: Vec<Vec<usize>>,
    homology_vertex_dimensions: Vec<Vec<usize>>,
}

macro_rules! impl_projective_presentation_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<qpa_gbnp::ProjectivePathAlgebraPresentation<$field>>
            for ProjectivePathPresentationKind
        {
            fn wrap(value: qpa_gbnp::ProjectivePathAlgebraPresentation<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::ProjectiveResolutionFpPathAlgebraModule<$field>>
            for GszProjectiveResolutionKind
        {
            fn wrap(value: qpa_core::ProjectiveResolutionFpPathAlgebraModule<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_projective_presentation_runtime!(RationalField, Rational);
impl_projective_presentation_runtime!(Z2, Binary);
impl_projective_presentation_runtime!(Zp, Prime);
impl_projective_presentation_runtime!(GaloisField, Extension);

fn clone_quotient_element_matrix(
    py: Python<'_>,
    rows: Vec<Vec<Py<PyQuotientElement>>>,
) -> Vec<Vec<QuotientElementKind>> {
    rows.into_iter()
        .map(|row| {
            row.into_iter()
                .map(|element| element.borrow(py).kind.clone())
                .collect()
        })
        .collect()
}

fn extract_quotient_element_matrix<F: runtime::PythonField>(
    rows: Vec<Vec<QuotientElementKind>>,
) -> PyResult<Vec<Vec<QuotientPathAlgebraElement<F>>>>
where
    QuotientElementKind: runtime::Extract<QuotientPathAlgebraElement<F>>,
{
    rows.into_iter().map(runtime::extract_all).collect()
}

fn clone_projective_vector_lists(
    py: Python<'_>,
    rows: Vec<Vec<Py<PyPathAlgebraVector>>>,
) -> Vec<Vec<projective::ProjectiveVectorKind>> {
    rows.into_iter()
        .map(|row| {
            row.into_iter()
                .map(|vector| vector.borrow(py).kind.clone())
                .collect()
        })
        .collect()
}

fn extract_projective_vector_lists<F: runtime::PythonField>(
    rows: Vec<Vec<projective::ProjectiveVectorKind>>,
) -> PyResult<Vec<Vec<PathAlgebraVector<F>>>>
where
    projective::ProjectiveVectorKind: runtime::Extract<PathAlgebraVector<F>>,
{
    rows.into_iter().map(runtime::extract_all).collect()
}

pub(super) fn explicit_gsz_resolution(
    algebra: &PyQuotientAlgebra,
    py: Python<'_>,
    ideal: &PyPathAlgebraIdeal,
    presentation_map: Vec<Vec<Py<PyQuotientElement>>>,
    term: Option<usize>,
) -> PyResult<PyGszProjectiveResolution> {
    let presentation_map = clone_quotient_element_matrix(py, presentation_map);
    let kind = match (&algebra.kind, &ideal.kind) {
        (QuotientKind::Rational(algebra), IdealKind::Rational(ideal)) => {
            let map = extract_quotient_element_matrix(presentation_map)?;
            GszProjectiveResolutionKind::Rational(match term {
                Some(term) => qpa_core::projective_resolution_fp_path_algebra_module(
                    algebra, ideal, &map, term,
                )
                .map_err(invalid)?,
                None => qpa_core::projective_resolution_fp_path_algebra_module_without_term(
                    algebra, ideal, &map,
                )
                .map_err(invalid)?,
            })
        }
        (QuotientKind::Binary(algebra), IdealKind::Binary(ideal)) => {
            let map = extract_quotient_element_matrix(presentation_map)?;
            GszProjectiveResolutionKind::Binary(match term {
                Some(term) => qpa_core::projective_resolution_fp_path_algebra_module(
                    algebra, ideal, &map, term,
                )
                .map_err(invalid)?,
                None => qpa_core::projective_resolution_fp_path_algebra_module_without_term(
                    algebra, ideal, &map,
                )
                .map_err(invalid)?,
            })
        }
        (QuotientKind::Prime(algebra), IdealKind::Prime(ideal)) => {
            let map = extract_quotient_element_matrix(presentation_map)?;
            GszProjectiveResolutionKind::Prime(match term {
                Some(term) => qpa_core::projective_resolution_fp_path_algebra_module(
                    algebra, ideal, &map, term,
                )
                .map_err(invalid)?,
                None => qpa_core::projective_resolution_fp_path_algebra_module_without_term(
                    algebra, ideal, &map,
                )
                .map_err(invalid)?,
            })
        }
        (QuotientKind::Extension(algebra), IdealKind::Extension(ideal)) => {
            let map = extract_quotient_element_matrix(presentation_map)?;
            GszProjectiveResolutionKind::Extension(match term {
                Some(term) => qpa_core::projective_resolution_fp_path_algebra_module(
                    algebra, ideal, &map, term,
                )
                .map_err(invalid)?,
                None => qpa_core::projective_resolution_fp_path_algebra_module_without_term(
                    algebra, ideal, &map,
                )
                .map_err(invalid)?,
            })
        }
        _ => {
            return Err(invalid(
                "the quotient algebra and resolution ideal use different coefficient fields",
            ));
        }
    };
    Ok(PyGszProjectiveResolution { kind })
}

type InducedMapOperation<F> =
    fn(
        &PathAlgebraModuleHomomorphism<F>,
        &PathAlgebraModuleHomomorphism<F>,
        &PathAlgebraModuleHomomorphism<F>,
        &PathAlgebraModuleHomomorphism<F>,
    ) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, qpa_core::PathAlgebraModuleError>;

pub(super) fn wrap_module<F: runtime::PythonField>(module: PathAlgebraModule<F>) -> PyMatrixModule
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    PyMatrixModule {
        kind: runtime::Wrap::wrap(module),
    }
}

fn wrap_module_element<F: runtime::PythonField>(
    element: PathAlgebraModuleElement<F>,
) -> PyModuleElement
where
    ModuleElementKind: runtime::Wrap<PathAlgebraModuleElement<F>>,
{
    PyModuleElement {
        kind: runtime::Wrap::wrap(element),
    }
}

fn wrap_module_elements<F: runtime::PythonField>(
    elements: Vec<PathAlgebraModuleElement<F>>,
) -> Vec<PyModuleElement>
where
    ModuleElementKind: runtime::Wrap<PathAlgebraModuleElement<F>>,
{
    elements.into_iter().map(wrap_module_element).collect()
}

fn wrap_module_homomorphism<F: runtime::PythonField>(
    homomorphism: PathAlgebraModuleHomomorphism<F>,
) -> PyModuleHomomorphism
where
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    PyModuleHomomorphism {
        kind: runtime::Wrap::wrap(homomorphism),
    }
}

fn wrap_module_homomorphisms<F: runtime::PythonField>(
    homomorphisms: Vec<PathAlgebraModuleHomomorphism<F>>,
) -> Vec<PyModuleHomomorphism>
where
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    homomorphisms
        .into_iter()
        .map(wrap_module_homomorphism)
        .collect()
}

fn wrap_module_homomorphism_pair<F: runtime::PythonField>(
    pair: Option<(
        PathAlgebraModuleHomomorphism<F>,
        PathAlgebraModuleHomomorphism<F>,
    )>,
) -> Option<(PyModuleHomomorphism, PyModuleHomomorphism)>
where
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    pair.map(|(first, second)| {
        (
            wrap_module_homomorphism(first),
            wrap_module_homomorphism(second),
        )
    })
}

fn scalar_rows<F: runtime::PythonField>(field: &F, rows: &[Vec<F::Element>]) -> Vec<Vec<PyScalar>> {
    rows.iter()
        .map(|row| row.iter().map(|value| field.wrap_scalar(value)).collect())
        .collect()
}

fn matrix_rows<F: runtime::PythonField>(matrix: &DenseMatrix<F>) -> Vec<Vec<PyScalar>> {
    scalar_rows(
        matrix.field(),
        &(0..matrix.nrows())
            .map(|row| matrix.row(row).to_vec())
            .collect::<Vec<_>>(),
    )
}

fn extract_scalar_rows<F: runtime::PythonField>(
    py: Python<'_>,
    field: &F,
    rows: Vec<Vec<Py<PyScalar>>>,
) -> PyResult<Vec<Vec<F::Element>>> {
    rows.into_iter()
        .map(|row| {
            row.into_iter()
                .map(|value| field.scalar_value(&value.borrow(py)))
                .collect()
        })
        .collect()
}

fn extract_matrix<F: runtime::PythonField>(
    py: Python<'_>,
    field: &F,
    rows: Vec<Vec<Py<PyScalar>>>,
) -> PyResult<DenseMatrix<F>> {
    DenseMatrix::from_rows(extract_scalar_rows(py, field, rows)?, field.clone()).map_err(invalid)
}

fn is_in_full_matrix_ring_for<F: runtime::PythonField>(
    py: Python<'_>,
    field: &F,
    rows: Vec<Vec<Py<PyScalar>>>,
) -> PyResult<bool> {
    let Ok(rows) = extract_scalar_rows(py, field, rows) else {
        return Ok(false);
    };
    let matrix = DenseMatrix::from_rows(rows, field.clone()).map_err(invalid)?;
    Ok(qpa_core::is_in_full_matrix_ring(&matrix, field))
}

#[pyfunction(name = "is_in_full_matrix_ring")]
pub(crate) fn py_is_in_full_matrix_ring(
    py: Python<'_>,
    rows: Vec<Vec<Py<PyScalar>>>,
    field: PyRef<'_, PyField>,
) -> PyResult<bool> {
    match field.kind {
        FieldKind::Rational => is_in_full_matrix_ring_for(py, &Q, rows),
        FieldKind::Binary => is_in_full_matrix_ring_for(py, &Z2, rows),
        FieldKind::Prime(characteristic) => {
            is_in_full_matrix_ring_for(py, &Zp::new(characteristic), rows)
        }
        FieldKind::Extension { .. } => is_in_full_matrix_ring_for(
            py,
            &field
                .kind
                .galois_field()
                .ok_or_else(|| invalid("expected an extension field"))?,
            rows,
        ),
    }
}

fn extract_scalars<F: runtime::PythonField>(
    py: Python<'_>,
    field: &F,
    values: Vec<Py<PyScalar>>,
) -> PyResult<Vec<F::Element>> {
    values
        .into_iter()
        .map(|value| field.scalar_value(&value.borrow(py)))
        .collect()
}

fn wrap_scalars<F: runtime::PythonField>(field: &F, values: &[F::Element]) -> Vec<PyScalar> {
    values
        .iter()
        .map(|value| field.wrap_scalar(value))
        .collect()
}

fn build_module_element<F: runtime::PythonField>(
    py: Python<'_>,
    module: &PathAlgebraModule<F>,
    coordinates: Vec<Vec<Py<PyScalar>>>,
) -> PyResult<PathAlgebraModuleElement<F>> {
    let field = module.right_acting_algebra().field();
    module
        .element_from_vertex_coordinates(extract_scalar_rows(py, field, coordinates)?)
        .map_err(invalid)
}

fn build_module_homomorphism<F: runtime::PythonField>(
    py: Python<'_>,
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
    matrices: Vec<Vec<Vec<Py<PyAny>>>>,
) -> PyResult<PathAlgebraModuleHomomorphism<F>> {
    let field = source.right_acting_algebra().field();
    let matrices = matrices
        .into_iter()
        .map(|rows| {
            let rows = rows
                .into_iter()
                .map(|row| {
                    row.into_iter()
                        .map(|value| crate::field::exact_scalar_for(field, value.bind(py)))
                        .collect::<PyResult<Vec<_>>>()
                })
                .collect::<PyResult<Vec<_>>>()?;
            DenseMatrix::from_rows(rows, field.clone()).map_err(invalid)
        })
        .collect::<PyResult<Vec<_>>>()?;
    qpa_core::right_module_hom_over_algebra(source.clone(), range.clone(), matrices)
        .map_err(invalid)
}

fn direct_sum_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>> + runtime::Wrap<PathAlgebraModule<F>>,
{
    Ok(
        qpa_core::direct_sum_of_qpa_modules(&runtime::extract_all(modules)?)
            .map_err(invalid)?
            .map(wrap_module),
    )
}

fn right_approximation_by_add_list_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
    target: &PathAlgebraModule<F>,
) -> PyResult<PyModuleHomomorphism>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>>,
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    Ok(wrap_module_homomorphism(
        qpa_gbnp::right_approximation_by_add_m_list(&runtime::extract_all(modules)?, target)
            .map_err(invalid)?,
    ))
}

fn make_basic_list_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
) -> PyResult<Vec<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>> + runtime::Wrap<PathAlgebraModule<F>>,
{
    Ok(
        qpa_core::make_basic_list_of_modules(&runtime::extract_all(modules)?)
            .map_err(invalid)?
            .into_iter()
            .map(wrap_module)
            .collect(),
    )
}

fn end_of_basic_module_as_quiver_algebra_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
    bound: Option<usize>,
) -> PyResult<PyEndOfBasicModuleAsQuiverAlgebraResult>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>>,
    EndOfBasicModuleAsQuiverAlgebraResultKind:
        runtime::Wrap<qpa_gbnp::EndOfBasicModuleAsQuiverAlgebraResult<F>>,
{
    let modules = runtime::extract_all(modules)?;
    let result = match bound {
        Some(bound) => qpa_gbnp::end_of_basic_module_as_quiver_algebra_with_bound(&modules, bound),
        None => qpa_gbnp::end_of_basic_module_as_quiver_algebra(&modules),
    }
    .map_err(invalid)?;
    Ok(PyEndOfBasicModuleAsQuiverAlgebraResult {
        kind: runtime::Wrap::wrap(result),
    })
}

fn radical_right_approximation_by_add_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
    index: &Integer,
) -> PyResult<PyModuleHomomorphism>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>>,
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    Ok(wrap_module_homomorphism(
        qpa_gbnp::radical_right_approximation_by_add_m(&runtime::extract_all(modules)?, index)
            .map_err(invalid)?,
    ))
}

fn simple_endo_projective_resolution_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
    index: &Integer,
    length: &Integer,
) -> PyResult<PySimpleEndoProjectiveResolutionResult>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>>,
    SimpleEndoProjectiveResolutionKind:
        runtime::Wrap<qpa_gbnp::ProjectiveResolutionOfSimpleModuleOverEndoResult<F>>,
{
    Ok(PySimpleEndoProjectiveResolutionResult {
        kind: runtime::Wrap::wrap(
            qpa_gbnp::projective_resolution_of_simple_module_over_endo(
                &runtime::extract_all(modules)?,
                index,
                length,
            )
            .map_err(invalid)?,
        ),
    })
}

fn faithful_dimension_result_text(result: qpa_gbnp::FaithfulDimensionResult) -> String {
    match result {
        qpa_gbnp::FaithfulDimensionResult::Finite(value) => value.to_string(),
        qpa_gbnp::FaithfulDimensionResult::Infinite => "infinity".into(),
    }
}

fn number_of_complements_result_text(result: qpa_gbnp::NumberOfComplementsResult) -> String {
    match result {
        qpa_gbnp::NumberOfComplementsResult::Finite(value) => value.to_string(),
        qpa_gbnp::NumberOfComplementsResult::Infinite => "infinity".into(),
    }
}

fn clone_module_elements(
    py: Python<'_>,
    elements: Vec<Py<PyModuleElement>>,
) -> Vec<ModuleElementKind> {
    elements
        .into_iter()
        .map(|element| element.borrow(py).kind.clone())
        .collect()
}

fn dimension_series_to_python(
    py: Python<'_>,
    series: ModuleDimensionSeries,
) -> PyResult<Py<PyAny>> {
    match series {
        ModuleDimensionSeries::ZeroDimensionVector(vector) => {
            Ok(PyList::new(py, vector)?.into_any().unbind())
        }
        ModuleDimensionSeries::Layers(layers) => Ok(PyList::new(py, layers)?.into_any().unbind()),
    }
}

fn matrices_of_module_for<F: runtime::PythonField>(
    py: Python<'_>,
    module: &PathAlgebraModule<F>,
) -> PyResult<Vec<Py<PyAny>>> {
    module
        .matrices_of_path_algebra_module()
        .into_iter()
        .map(|matrix| match matrix {
            PathAlgebraModuleMatrix::ZeroScalar => {
                crate::field::wrap_python_integer(py, &Integer::from(0))
            }
            PathAlgebraModuleMatrix::Matrix(matrix) => {
                Ok(PyList::new(py, matrix_rows(&matrix))?.into_any().unbind())
            }
        })
        .collect()
}

fn induced_map_for<F: runtime::PythonField>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: ModuleHomomorphismKind,
    beta: ModuleHomomorphismKind,
    alpha: ModuleHomomorphismKind,
    operation: InducedMapOperation<F>,
) -> PyResult<Option<PyModuleHomomorphism>>
where
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>
        + runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    let second = runtime::Extract::extract(second)?;
    let beta = runtime::Extract::extract(beta)?;
    let alpha = runtime::Extract::extract(alpha)?;
    Ok(operation(first, &second, &beta, &alpha)
        .map_err(invalid)?
        .map(wrap_module_homomorphism))
}

fn multiply_lists_of_maps_for<F: runtime::PythonField>(
    projections: Vec<ModuleHomomorphismKind>,
    matrix: Vec<Vec<ModuleHomomorphismKind>>,
    inclusions: Vec<ModuleHomomorphismKind>,
) -> PyResult<PyModuleHomomorphism>
where
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>
        + runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    let projections = runtime::extract_all(projections)?;
    let matrix = matrix
        .into_iter()
        .map(runtime::extract_all)
        .collect::<PyResult<Vec<_>>>()?;
    let inclusions = runtime::extract_all(inclusions)?;
    Ok(wrap_module_homomorphism(
        qpa_core::multiply_lists_of_maps(&projections, &matrix, &inclusions).map_err(invalid)?,
    ))
}

fn hom_basis_coefficients_for<F: runtime::PythonField>(
    basis: &qpa_core::PathAlgebraModuleHomomorphismBasis<F>,
    homomorphism: ModuleHomomorphismKind,
) -> PyResult<Option<Vec<PyScalar>>>
where
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>,
{
    let homomorphism = runtime::Extract::extract(homomorphism)?;
    Ok(basis
        .coefficients(&homomorphism)
        .map_err(invalid)?
        .map(|coefficients| wrap_scalars(basis.field(), &coefficients)))
}

fn restrict_module_for<F: runtime::PythonField>(
    module: &PathAlgebraModule<F>,
    homomorphism: AlgebraHomomorphismKind,
) -> PyResult<PyMatrixModule>
where
    AlgebraHomomorphismKind: runtime::Extract<qpa_core::QuiverAlgebraHomomorphism<F>>,
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    let homomorphism = runtime::Extract::extract(homomorphism)?;
    Ok(wrap_module(
        qpa_core::restriction_via_algebra_homomorphism(&homomorphism, module).map_err(invalid)?,
    ))
}

fn restrict_module_homomorphism_for<F: runtime::PythonField>(
    module_homomorphism: &PathAlgebraModuleHomomorphism<F>,
    algebra_homomorphism: AlgebraHomomorphismKind,
) -> PyResult<PyModuleHomomorphism>
where
    AlgebraHomomorphismKind: runtime::Extract<qpa_core::QuiverAlgebraHomomorphism<F>>,
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    let algebra_homomorphism = runtime::Extract::extract(algebra_homomorphism)?;
    Ok(wrap_module_homomorphism(
        qpa_core::restriction_via_algebra_homomorphism_map(
            &algebra_homomorphism,
            module_homomorphism,
        )
        .map_err(invalid)?,
    ))
}

fn wrap_module_submodule<F: runtime::PythonField>(
    submodule: PathAlgebraModuleSubmodule<F>,
) -> PyModuleSubmodule
where
    ModuleSubmoduleKind: runtime::Wrap<PathAlgebraModuleSubmodule<F>>,
{
    PyModuleSubmodule {
        kind: runtime::Wrap::wrap(submodule),
    }
}

#[pymethods]
impl PyProjectivePathAlgebraPresentation {
    #[getter]
    fn projective(&self) -> PyRightProjectiveModule {
        dispatch_ref!(&self.kind, ProjectivePathPresentationKind, |value| {
            projective::wrap_projective_module(value.projective().clone())
        })
    }

    #[getter]
    fn submodule(&self) -> Option<PyRightProjectiveSubmodule> {
        dispatch_ref!(&self.kind, ProjectivePathPresentationKind, |value| value
            .submodule()
            .cloned()
            .map(projective::wrap_projective_submodule))
    }

    #[getter]
    fn projective_generators(&self) -> Vec<PyPathAlgebraVector> {
        dispatch_ref!(&self.kind, ProjectivePathPresentationKind, |value| {
            projective::wrap_vectors(value.projective_generators().to_vec())
        })
    }

    #[getter]
    fn syzygy_generators(&self) -> Vec<PyPathAlgebraVector> {
        dispatch_ref!(&self.kind, ProjectivePathPresentationKind, |value| {
            projective::wrap_vectors(value.syzygy_generators().to_vec())
        })
    }

    #[getter]
    fn matrix(&self) -> Vec<Vec<PyAlgebraElement>> {
        dispatch_ref!(&self.kind, ProjectivePathPresentationKind, |value| value
            .matrix()
            .iter()
            .map(|row| wrap_algebra_elements(row.clone()))
            .collect())
    }

    #[getter]
    fn reduced_groebner_basis(&self) -> Vec<PyPathAlgebraVector> {
        dispatch_ref!(&self.kind, ProjectivePathPresentationKind, |value| {
            projective::wrap_vectors(value.reduced_groebner_basis().to_vec())
        })
    }

    #[getter]
    fn relation_ideal(&self) -> PyPathAlgebraIdeal {
        let kind = map_ref!(
            &self.kind,
            ProjectivePathPresentationKind => IdealKind,
            |value| value.relation_ideal().clone()
        );
        PyPathAlgebraIdeal { kind }
    }

    fn __len__(&self) -> usize {
        dispatch_ref!(&self.kind, ProjectivePathPresentationKind, |value| value
            .result_length())
    }
}

#[pymethods]
impl PyGszProjectiveResolution {
    #[getter]
    #[allow(clippy::unused_self)]
    fn is_projective_resolution_fp_path_algebra_module(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn parent_algebra(&self) -> PyQuotientAlgebra {
        let kind = map_ref!(
            &self.kind,
            GszProjectiveResolutionKind => QuotientKind,
            |value| value.parent_algebra().clone()
        );
        PyQuotientAlgebra { kind }
    }

    #[getter]
    fn ring_ideal(&self) -> PyPathAlgebraIdeal {
        let kind = map_ref!(
            &self.kind,
            GszProjectiveResolutionKind => IdealKind,
            |value| value.ring_ideal().clone()
        );
        PyPathAlgebraIdeal { kind }
    }

    #[getter]
    fn projectives(&self) -> Vec<PyRightProjectiveQuotientModule> {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| value
            .projectives()
            .iter()
            .cloned()
            .map(quotient_projective::wrap_quotient_projective)
            .collect())
    }

    #[getter]
    fn projectives_attribute_is_bound(&self) -> bool {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| value
            .projectives_attribute_is_bound())
    }

    #[getter]
    fn r_projectives(&self) -> Vec<PyRightProjectiveModule> {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| value
            .r_projectives()
            .iter()
            .cloned()
            .map(projective::wrap_projective_module)
            .collect())
    }

    #[getter]
    fn maps(&self) -> Vec<Vec<Vec<PyAlgebraElement>>> {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| value
            .maps()
            .iter()
            .map(|matrix| matrix
                .iter()
                .map(|row| wrap_algebra_elements(row.clone()))
                .collect())
            .collect())
    }

    #[getter]
    fn f_lists(&self) -> Vec<Vec<PyPathAlgebraVector>> {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| value
            .f_lists()
            .iter()
            .map(|vectors| projective::wrap_vectors(vectors.clone()))
            .collect())
    }

    #[getter]
    fn f_prime_lists(&self) -> Vec<Vec<PyPathAlgebraVector>> {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| value
            .f_prime_lists()
            .iter()
            .map(|vectors| projective::wrap_vectors(vectors.clone()))
            .collect())
    }

    #[getter]
    fn has_r_projectives_vertex_list(&self) -> bool {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| value
            .r_projectives_vertex_list()
            .is_bound())
    }

    #[getter]
    fn r_projectives_vertex_list(&self) -> Option<Vec<Vec<usize>>> {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| {
            match value.r_projectives_vertex_list() {
                qpa_core::CachedAttribute::Value(rows) => Some(
                    rows.iter()
                        .map(|row| row.iter().map(|vertex| vertex.get()).collect())
                        .collect(),
                ),
                _ => None,
            }
        })
    }

    fn set_r_projectives_vertex_list(&mut self, rows: Vec<Vec<usize>>) -> PyResult<()> {
        let rows = rows
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|vertex| {
                        VertexIndex::new(vertex)
                            .ok_or_else(|| invalid("vertex indices must be one-based and positive"))
                    })
                    .collect()
            })
            .collect::<PyResult<Vec<_>>>()?;
        match &mut self.kind {
            GszProjectiveResolutionKind::Rational(value) => {
                value.set_r_projectives_vertex_list(rows);
            }
            GszProjectiveResolutionKind::Binary(value) => {
                value.set_r_projectives_vertex_list(rows);
            }
            GszProjectiveResolutionKind::Prime(value) => {
                value.set_r_projectives_vertex_list(rows);
            }
            GszProjectiveResolutionKind::Extension(value) => {
                value.set_r_projectives_vertex_list(rows);
            }
        }
        Ok(())
    }

    #[getter]
    fn has_projectives_f_prime_list(&self) -> bool {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| value
            .projectives_f_prime_list()
            .is_bound())
    }

    #[getter]
    fn projectives_f_prime_list(&self) -> Option<Vec<Vec<PyPathAlgebraVector>>> {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| {
            match value.projectives_f_prime_list() {
                qpa_core::CachedAttribute::Value(rows) => Some(
                    rows.iter()
                        .map(|row| projective::wrap_vectors(row.clone()))
                        .collect(),
                ),
                _ => None,
            }
        })
    }

    fn set_projectives_f_prime_list(
        &mut self,
        py: Python<'_>,
        rows: Vec<Vec<Py<PyPathAlgebraVector>>>,
    ) -> PyResult<()> {
        let rows = clone_projective_vector_lists(py, rows);
        match &mut self.kind {
            GszProjectiveResolutionKind::Rational(value) => {
                value.set_projectives_f_prime_list(extract_projective_vector_lists(rows)?);
            }
            GszProjectiveResolutionKind::Binary(value) => {
                value.set_projectives_f_prime_list(extract_projective_vector_lists(rows)?);
            }
            GszProjectiveResolutionKind::Prime(value) => {
                value.set_projectives_f_prime_list(extract_projective_vector_lists(rows)?);
            }
            GszProjectiveResolutionKind::Extension(value) => {
                value.set_projectives_f_prime_list(extract_projective_vector_lists(rows)?);
            }
        }
        Ok(())
    }

    fn find_next_r_projective(&self, ring: &PyPathAlgebra, term: usize) -> PyResult<()> {
        match (&self.kind, &ring.kind) {
            (GszProjectiveResolutionKind::Rational(resolution), AlgebraKind::Rational(ring)) => {
                qpa_core::find_next_r_projective(resolution, ring, term).map_err(invalid)
            }
            (GszProjectiveResolutionKind::Binary(resolution), AlgebraKind::Binary(ring)) => {
                qpa_core::find_next_r_projective(resolution, ring, term).map_err(invalid)
            }
            (GszProjectiveResolutionKind::Prime(resolution), AlgebraKind::Prime(ring)) => {
                qpa_core::find_next_r_projective(resolution, ring, term).map_err(invalid)
            }
            (GszProjectiveResolutionKind::Extension(resolution), AlgebraKind::Extension(ring)) => {
                qpa_core::find_next_r_projective(resolution, ring, term).map_err(invalid)
            }
            _ => Err(invalid(
                "QPA declares FindNextRProjective but installs no method",
            )),
        }
    }

    fn find_next_syzygy(&self, term: usize) -> PyResult<()> {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |resolution| {
            qpa_core::find_next_syzygy(resolution, term).map_err(invalid)
        })
    }

    #[getter]
    fn finite_at_projective(&self) -> Option<usize> {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| value
            .finite_at_projective())
    }

    fn semisimple_homology(&self) -> PyResult<PySemisimpleResolutionHomology> {
        dispatch_ref!(&self.kind, GszProjectiveResolutionKind, |value| {
            let result = value.semisimple_homology().map_err(invalid)?;
            Ok(PySemisimpleResolutionHomology {
                projective_vertex_multiplicities: (result
                    .projective_vertex_multiplicities()
                    .to_vec()),
                differential_vertex_ranks: (result.differential_vertex_ranks().to_vec()),
                homology_vertex_dimensions: (result.homology_vertex_dimensions().to_vec()),
            })
        })
    }
}

#[pymethods]
impl PySemisimpleResolutionHomology {
    #[getter]
    fn projective_vertex_multiplicities(&self) -> Vec<Vec<usize>> {
        self.projective_vertex_multiplicities.clone()
    }

    #[getter]
    fn differential_vertex_ranks(&self) -> Vec<Vec<usize>> {
        self.differential_vertex_ranks.clone()
    }

    #[getter]
    fn homology_vertex_dimensions(&self) -> Vec<Vec<usize>> {
        self.homology_vertex_dimensions.clone()
    }

    #[getter]
    fn minimal_projective_betti_multiplicities(&self) -> Vec<Vec<usize>> {
        self.homology_vertex_dimensions.clone()
    }
}

#[pymethods]
impl PyMatrixModule {
    fn projective_path_algebra_presentation(
        &self,
    ) -> PyResult<PyProjectivePathAlgebraPresentation> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ProjectivePathPresentationKind,
            |module| qpa_gbnp::projective_path_algebra_presentation(module).map_err(invalid)?
        );
        Ok(PyProjectivePathAlgebraPresentation { kind })
    }

    fn projective_resolution_of_path_algebra_module(
        &self,
        term: usize,
    ) -> PyResult<PyGszProjectiveResolution> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => GszProjectiveResolutionKind,
            |module| qpa_gbnp::projective_resolution_of_path_algebra_module(module, term)
                .map_err(invalid)?
        );
        Ok(PyGszProjectiveResolution { kind })
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_algebra_module(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_path_algebra_mat_module(&self) -> bool {
        true
    }

    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .right_acting_algebra()
            .field()
            .python_kind());
        PyField { kind }
    }

    #[getter]
    fn right_acting_algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            conversion::wrap_quiver_algebra(module.right_acting_algebra().clone())
        })
    }

    #[getter]
    fn dimension_vector(&self) -> Vec<usize> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .dimension_vector()
            .to_vec())
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module.dimension())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module.is_zero())
    }

    #[getter]
    fn arrow_matrix_shapes(&self) -> Vec<(usize, usize)> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| matrix_shapes(
            module.arrow_matrices()
        ))
    }

    #[getter]
    fn arrow_matrices(&self) -> Vec<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .arrow_matrices()
            .iter()
            .map(matrix_rows)
            .collect())
    }

    fn zero_element(&self) -> PyModuleElement {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleElementKind,
            |module| module.zero_element()
        );
        PyModuleElement { kind }
    }

    fn element(
        &self,
        py: Python<'_>,
        coordinates: Vec<Vec<Py<PyScalar>>>,
    ) -> PyResult<PyModuleElement> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleElementKind,
            |module| build_module_element(py, module, coordinates)?
        );
        Ok(PyModuleElement { kind })
    }

    fn basis_elements(&self) -> Vec<PyModuleElement> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| wrap_module_elements(
            module.basis_elements()
        ))
    }

    fn element_space(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyModuleElement>>,
    ) -> PyResult<PyModuleElementSpace> {
        let generators = clone_module_elements(py, generators);
        let kind = match &self.kind {
            MatrixModuleKind::Rational(module) => ModuleElementSpaceKind::Rational(
                PathAlgebraModuleElementSpace::by_generators(
                    module.clone(),
                    runtime::extract_all(generators)?,
                )
                .map_err(invalid)?,
            ),
            MatrixModuleKind::Binary(module) => ModuleElementSpaceKind::Binary(
                PathAlgebraModuleElementSpace::by_generators(
                    module.clone(),
                    runtime::extract_all(generators)?,
                )
                .map_err(invalid)?,
            ),
            MatrixModuleKind::Prime(module) => ModuleElementSpaceKind::Prime(
                PathAlgebraModuleElementSpace::by_generators(
                    module.clone(),
                    runtime::extract_all(generators)?,
                )
                .map_err(invalid)?,
            ),
            MatrixModuleKind::Extension(module) => ModuleElementSpaceKind::Extension(
                PathAlgebraModuleElementSpace::by_generators(
                    module.clone(),
                    runtime::extract_all(generators)?,
                )
                .map_err(invalid)?,
            ),
        };
        Ok(PyModuleElementSpace { kind })
    }

    fn canonical_basis(&self) -> PyResult<PyModuleElementBasis> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleElementBasisKind,
            |module| PathAlgebraModuleElementBasis::canonical_for_module(module)
                .map_err(invalid)?
        );
        Ok(PyModuleElementBasis { kind })
    }

    fn element_from_external_representation(
        &self,
        py: Python<'_>,
        coordinates: Vec<Vec<Py<PyScalar>>>,
    ) -> PyResult<PyModuleElement> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleElementKind,
            |module| module
                .element_from_external_representation(extract_scalar_rows(
                    py,
                    module.right_acting_algebra().field(),
                    coordinates,
                )?)
                .map_err(invalid)?
        );
        Ok(PyModuleElement { kind })
    }

    fn matrices_of_path_algebra_module(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            matrices_of_module_for(py, module)
        })
    }

    fn dimension_vector_leq(&self, other: &Self) -> bool {
        // QPA's method filter only requires two matrix modules: neither the
        // coefficient field nor the acting algebra has to agree. GAP list
        // subtraction also pads the shorter vector with zeros.
        let left = self.dimension_vector();
        let right = other.dimension_vector();
        let length = left.len().max(right.len());
        (0..length).all(|index| {
            left.get(index).copied().unwrap_or(0) <= right.get(index).copied().unwrap_or(0)
        })
    }

    fn annihilator(&self) -> PyResult<Vec<PyQuiverAlgebraElement>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            Ok(tensor::wrap_quiver_elements(
                &qpa_gbnp::annihilator_of_module(module).map_err(invalid)?,
            ))
        })
    }

    fn radical_series(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            dimension_series_to_python(py, qpa_core::radical_series(module).map_err(invalid)?)
        })
    }

    fn socle_series(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            dimension_series_to_python(py, qpa_gbnp::socle_series(module).map_err(invalid)?)
        })
    }

    fn lift_idempotents(&self) -> PyResult<()> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_core::lift_idempotents(module).map_err(invalid)
        })
    }

    #[staticmethod]
    fn direct_sum(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
    ) -> PyResult<Option<PyMatrixModule>> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        match modules.first() {
            None => Ok(None),
            Some(MatrixModuleKind::Rational(_)) => direct_sum_for::<RationalField>(modules),
            Some(MatrixModuleKind::Binary(_)) => direct_sum_for::<Z2>(modules),
            Some(MatrixModuleKind::Prime(_)) => direct_sum_for::<Zp>(modules),
            Some(MatrixModuleKind::Extension(_)) => direct_sum_for::<GaloisField>(modules),
        }
    }

    #[getter]
    fn is_direct_sum(&self) -> bool {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .is_direct_sum_of_modules())
    }

    fn direct_sum_summands(&self) -> Option<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .direct_sum_summands()
            .map(|summands| summands.iter().cloned().map(wrap_module).collect()))
    }

    fn direct_sum_projections(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(qpa_core::direct_sum_projections(module).map_err(invalid)?,)
        ))
    }

    fn direct_sum_inclusions(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(qpa_core::direct_sum_inclusions(module).map_err(invalid)?,)
        ))
    }

    fn decompose(&self) -> PyResult<Vec<Self>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_core::decompose_module(module,)
                .map_err(invalid)?
                .into_iter()
                .map(wrap_module)
                .collect()
        ))
    }

    #[getter]
    fn is_indecomposable(&self) -> PyResult<Option<bool>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_core::is_indecomposable_module(module).map_err(invalid)
        })
    }

    fn decompose_with_inclusions(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(
                qpa_core::decompose_module_with_inclusions(module).map_err(invalid)?,
            )
        ))
    }

    fn decompose_with_multiplicities(
        &self,
    ) -> PyResult<Option<PyModuleDecompositionWithMultiplicities>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_core::decompose_module_with_multiplicities(module)
                .map_err(invalid)?
                .map(|result| PyModuleDecompositionWithMultiplicities {
                    kind: runtime::Wrap::wrap(result),
                })
        ))
    }

    fn block_splitting_idempotents(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(
                qpa_core::block_splitting_idempotents(module).map_err(invalid)?,
            )
        ))
    }

    fn block_decomposition(&self) -> PyResult<Vec<Self>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_core::block_decomposition_of_module(module)
                .map_err(invalid)?
                .into_iter()
                .map(wrap_module)
                .collect()
        ))
    }

    fn basic_version(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_core::basic_version_of_module(module).map_err(invalid)?,
        )))
    }

    #[staticmethod]
    fn make_basic(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
    ) -> PyResult<Vec<PyMatrixModule>> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        match modules.first() {
            Some(MatrixModuleKind::Rational(_)) | None => {
                make_basic_list_for::<RationalField>(modules)
            }
            Some(MatrixModuleKind::Binary(_)) => make_basic_list_for::<Z2>(modules),
            Some(MatrixModuleKind::Prime(_)) => make_basic_list_for::<Zp>(modules),
            Some(MatrixModuleKind::Extension(_)) => make_basic_list_for::<GaloisField>(modules),
        }
    }

    fn nonisomorphic_direct_summand_data(&self) -> PyResult<PyNonisomorphicDirectSummandData> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            PyNonisomorphicDirectSummandData {
                value: qpa_core::number_of_nonisomorphic_direct_summands(module)
                    .map_err(invalid)?,
            }
        ))
    }

    fn isomorphic_to(&self, other: &Self) -> PyResult<bool> {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            MatrixModuleKind,
            |source, range| qpa_core::isomorphic_modules(source, range).map_err(invalid),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn isomorphism_to(&self, other: &Self) -> PyResult<Option<PyModuleHomomorphism>> {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            MatrixModuleKind,
            |source, range| Ok(qpa_core::isomorphism_of_modules(source, range)
                .map_err(invalid)?
                .map(wrap_module_homomorphism)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn is_direct_summand_of(&self, module: &Self) -> PyResult<bool> {
        dispatch_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind,
            |candidate, module| qpa_core::is_direct_summand(candidate, module).map_err(invalid),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn is_in_additive_closure_of(&self, generator: &Self) -> PyResult<Option<bool>> {
        dispatch_pair_ref!(
            &self.kind,
            &generator.kind,
            MatrixModuleKind,
            |candidate, generator| qpa_core::is_in_additive_closure(candidate, generator)
                .map_err(invalid),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn common_direct_summand(&self, other: &Self) -> PyResult<PyCommonDirectSummandResult> {
        let kind = map_pair_ref!(
            &self.kind,
            &other.kind,
            MatrixModuleKind => CommonDirectSummandResultKind,
            |source, range| qpa_core::common_direct_summand(source, range).map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyCommonDirectSummandResult { kind })
    }

    fn maximal_common_direct_summand(
        &self,
        other: &Self,
    ) -> PyResult<Option<PyMaximalCommonDirectSummand>> {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            MatrixModuleKind,
            |source, range| Ok(qpa_core::maximal_common_direct_summand(source, range)
                .map_err(invalid)?
                .map(|result| PyMaximalCommonDirectSummand {
                    kind: runtime::Wrap::wrap(result),
                })),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn identity(&self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module| qpa_core::identity_module_mapping(module.clone()).map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    fn zero_map_to(&self, range: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, range| qpa_core::zero_module_mapping(source.clone(), range.clone())
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    /// Return QPA's ordered basis of homomorphisms to `range`, or `None`.
    fn hom_basis_to(&self, range: &Self) -> PyResult<Option<Vec<PyModuleHomomorphism>>> {
        dispatch_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind,
            |source, range| Ok(qpa_core::hom_over_algebra(source, range)
                .map_err(invalid)?
                .map(wrap_module_homomorphisms)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    /// Return the Hom dimension without constructing basis maps, or `None`.
    fn hom_dimension_to(&self, range: &Self) -> PyResult<Option<usize>> {
        dispatch_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind,
            |source, range| qpa_core::hom_dimension_over_algebra(source, range).map_err(invalid),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn hom_basis_with_coefficients_to(
        &self,
        range: &Self,
    ) -> PyResult<Option<PyModuleHomomorphismBasis>> {
        dispatch_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind,
            |source, range| Ok(
                qpa_core::hom_over_algebra_with_basis_function(source, range)
                    .map_err(invalid)?
                    .map(|basis| PyModuleHomomorphismBasis {
                        kind: runtime::Wrap::wrap(basis),
                    })
            ),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn endomorphism_algebra(&self) -> PyResult<PyModuleEndomorphismAlgebra> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleEndomorphismAlgebraKind,
            |module| qpa_core::end_over_algebra(module).map_err(invalid)?
        );
        Ok(PyModuleEndomorphismAlgebra { kind })
    }

    fn endomorphism_as_quiver_algebra(&self) -> PyResult<PyEndOfModuleAsQuiverAlgebraResult> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => EndOfModuleAsQuiverAlgebraResultKind,
            |module| qpa_gbnp::end_of_module_as_quiver_algebra(module).map_err(invalid)?
        );
        Ok(PyEndOfModuleAsQuiverAlgebraResult { kind })
    }

    #[staticmethod]
    #[pyo3(signature = (modules, bound=None))]
    fn endomorphism_as_quiver_algebra_of_basic(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
        bound: Option<usize>,
    ) -> PyResult<PyEndOfBasicModuleAsQuiverAlgebraResult> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        match modules.first() {
            Some(MatrixModuleKind::Rational(_)) | None => {
                end_of_basic_module_as_quiver_algebra_for::<RationalField>(modules, bound)
            }
            Some(MatrixModuleKind::Binary(_)) => {
                end_of_basic_module_as_quiver_algebra_for::<Z2>(modules, bound)
            }
            Some(MatrixModuleKind::Prime(_)) => {
                end_of_basic_module_as_quiver_algebra_for::<Zp>(modules, bound)
            }
            Some(MatrixModuleKind::Extension(_)) => {
                end_of_basic_module_as_quiver_algebra_for::<GaloisField>(modules, bound)
            }
        }
    }

    #[staticmethod]
    fn adr_algebra_of_basic(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
    ) -> PyResult<PyQuiverAlgebra> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        match modules.first() {
            Some(MatrixModuleKind::Rational(_)) | None => {
                adr_algebra_of_basic_module_for::<RationalField>(modules)
            }
            Some(MatrixModuleKind::Binary(_)) => adr_algebra_of_basic_module_for::<Z2>(modules),
            Some(MatrixModuleKind::Prime(_)) => adr_algebra_of_basic_module_for::<Zp>(modules),
            Some(MatrixModuleKind::Extension(_)) => {
                adr_algebra_of_basic_module_for::<GaloisField>(modules)
            }
        }
    }

    fn stable_endomorphism_quotient(&self) -> PyResult<PyModuleEndomorphismQuotientMap> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleEndomorphismQuotientMapKind,
            |module| qpa_gbnp::end_modulo_projectives(module).map_err(invalid)?
        );
        Ok(PyModuleEndomorphismQuotientMap { kind })
    }

    fn endomorphism_from_matrix(
        &self,
        py: Python<'_>,
        matrix: Vec<Vec<Py<PyScalar>>>,
    ) -> PyResult<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            let field = module.right_acting_algebra().field();
            Ok(wrap_module_homomorphism(
                qpa_core::from_end_m_to_hom_mm(module, &extract_matrix(py, field, matrix)?)
                    .map_err(invalid)?,
            ))
        })
    }

    fn hom_factoring_through_projectives_to(
        &self,
        range: &Self,
    ) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind,
            |source, range| Ok(wrap_module_homomorphisms(
                qpa_gbnp::hom_factoring_through_projectives(source, range).map_err(invalid)?,
            )),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn reject_of(&self, rejected_by: &Self) -> PyResult<PyModuleHomomorphism> {
        dispatch_pair_ref!(
            &self.kind,
            &rejected_by.kind,
            MatrixModuleKind,
            |module, rejected_by| Ok(wrap_module_homomorphism(
                qpa_core::reject_of_module(module, rejected_by).map_err(invalid)?,
            )),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn simple_submodule_inclusions(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(
                qpa_core::all_simple_submodules_of_module(module).map_err(invalid)?,
            )
        ))
    }

    fn all_submodules(&self) -> PyResult<PyAllSubmodulesOfModuleResult> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => AllSubmodulesOfModuleResultKind,
            |module| qpa_core::all_submodules_of_module(module).map_err(invalid)?
        );
        Ok(PyAllSubmodulesOfModuleResult { kind })
    }

    fn restrict_scalars(&self, homomorphism: &PyAlgebraHomomorphism) -> PyResult<PyMatrixModule> {
        let homomorphism = homomorphism.kind.clone();
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            restrict_module_for(module, homomorphism)
        })
    }

    fn homomorphism_from_images(
        &self,
        py: Python<'_>,
        range: &Self,
        images: Vec<Py<PyModuleElement>>,
    ) -> PyResult<PyModuleHomomorphism> {
        let images = clone_module_elements(py, images);
        dispatch_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind,
            |source, range| Ok(wrap_module_homomorphism(
                qpa_core::homomorphism_from_images(source, range, &runtime::extract_all(images)?,)
                    .map_err(invalid)?,
            )),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    /// Construct a checked map from exact Scalar, int, or Fraction matrices.
    ///
    /// Matrices follow vertex order and row action. Floats and scalars from
    /// other fields are rejected; every arrow square must commute.
    fn homomorphism_to(
        &self,
        py: Python<'_>,
        range: &Self,
        vertex_matrices: Vec<Vec<Vec<Py<PyAny>>>>,
    ) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, range| build_module_homomorphism(py, source, range, vertex_matrices)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn radical(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_core::radical_of_module(module.clone()).map_err(invalid)?,
        )))
    }

    fn top(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_core::top_of_module(module.clone()).map_err(invalid)?,
        )))
    }

    fn socle(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_core::socle_of_module(module.clone()).map_err(invalid)?,
        )))
    }

    #[getter]
    fn is_simple(&self) -> bool {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_core::is_simple_qpa_module(module)
        })
    }

    #[getter]
    fn is_semisimple(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_core::is_semisimple_module(module).map_err(invalid)
        })
    }

    #[getter]
    fn loewy_length(&self) -> PyResult<usize> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_core::loewy_length_of_module(module.clone()).map_err(invalid)
        })
    }

    fn minimal_generating_set(&self) -> PyResult<Vec<PyModuleElement>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_elements(
                qpa_core::minimal_generating_set_of_module(module.clone()).map_err(invalid)?,
            )
        ))
    }

    fn submodule(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyModuleElement>>,
    ) -> PyResult<PyModuleSubmodule> {
        let generators = clone_module_elements(py, generators);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_submodule(
                qpa_core::subalgebra_module(module.clone(), &runtime::extract_all(generators)?,)
                    .map_err(invalid)?,
            )
        ))
    }

    fn subrepresentation(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyModuleElement>>,
    ) -> PyResult<Self> {
        let generators = clone_module_elements(py, generators);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_core::subrepresentation(module.clone(), &runtime::extract_all(generators)?,)
                .map_err(invalid)?,
        )))
    }

    fn subrepresentation_inclusion(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyModuleElement>>,
    ) -> PyResult<PyModuleHomomorphism> {
        let generators = clone_module_elements(py, generators);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(
                qpa_core::subrepresentation_inclusion(
                    module.clone(),
                    &runtime::extract_all(generators)?,
                )
                .map_err(invalid)?,
            )
        ))
    }

    fn radical_inclusion(&self) -> PyResult<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(
                qpa_core::radical_of_module_inclusion(module.clone()).map_err(invalid)?,
            )
        ))
    }

    fn top_projection(&self) -> PyResult<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(
                qpa_core::top_of_module_projection(module.clone()).map_err(invalid)?,
            )
        ))
    }

    fn socle_inclusion(&self) -> PyResult<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(
                qpa_core::socle_of_module_inclusion(module.clone()).map_err(invalid)?,
            )
        ))
    }

    #[getter]
    fn is_projective(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_projective_module(module).map_err(invalid)
        })
    }

    #[getter]
    fn is_injective(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_injective_module(module).map_err(invalid)
        })
    }

    fn compare_with_indecomposable_projective(&self) -> PyResult<Option<usize>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::compare_with_indec_projective(module).map_err(invalid)
        })
    }

    fn compare_with_indecomposable_injective(&self) -> PyResult<Option<usize>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::compare_with_indec_injective(module).map_err(invalid)
        })
    }

    fn projective_description(&self) -> PyResult<Vec<Option<usize>>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::description_of_projective_module(module).map_err(invalid)
        })
    }

    fn direct_sum_minus_summands(
        &self,
        removed: Vec<usize>,
    ) -> PyResult<(PyModuleHomomorphism, PyModuleHomomorphism)> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            let (inclusion, projection) =
                qpa_gbnp::direct_sum_minus_summands(module, &removed).map_err(invalid)?;
            Ok((
                wrap_module_homomorphism(inclusion),
                wrap_module_homomorphism(projection),
            ))
        })
    }

    /// Return QPA's projective-cover map onto this module.
    fn projective_cover(&self) -> PyResult<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(qpa_gbnp::projective_cover(module).map_err(invalid)?,)
        ))
    }

    fn injective_envelope(&self) -> PyResult<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(qpa_gbnp::injective_envelope(module).map_err(invalid)?,)
        ))
    }

    /// Return the kernel of QPA's projective cover.
    fn first_syzygy(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_gbnp::first_syzygy(module).map_err(invalid)?,
        )))
    }

    fn nth_syzygy(&self, degree: i64) -> PyResult<Self> {
        let degree = Integer::from(degree);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_gbnp::nth_syzygy(module, &degree).map_err(invalid)?,
        )))
    }

    /// Return QPA's lazy projective resolution, with this module in degree `-1`.
    fn projective_resolution(&self) -> PyResult<PyProjectiveResolution> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ProjectiveResolutionKind,
            |module| qpa_gbnp::projective_resolution(module).map_err(invalid)?
        );
        Ok(PyProjectiveResolution { kind })
    }

    fn injective_resolution(&self) -> PyResult<PyInjectiveResolution> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => InjectiveResolutionKind,
            |module| qpa_gbnp::injective_resolution(module).map_err(invalid)?
        );
        Ok(PyInjectiveResolution { kind })
    }

    fn zero_complex(&self) -> PyResult<PyModuleComplex> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleComplexKind,
            |module| qpa_gbnp::zero_complex(
                qpa_gbnp::cat_of_right_algebra_modules(module.right_acting_algebra())
                    .map_err(invalid)?,
            )
            .map_err(invalid)?
        );
        Ok(PyModuleComplex { kind })
    }

    fn stalk_complex(&self, degree: i64) -> PyResult<PyModuleComplex> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleComplexKind,
            |module| qpa_gbnp::stalk_complex(
                qpa_gbnp::cat_of_right_algebra_modules(module.right_acting_algebra())
                    .map_err(invalid)?,
                module,
                degree,
            )
            .map_err(invalid)?
        );
        Ok(PyModuleComplex { kind })
    }

    fn finite_complex(
        &self,
        py: Python<'_>,
        base_degree: i64,
        differentials: Vec<Py<PyModuleHomomorphism>>,
    ) -> PyResult<PyModuleComplex> {
        let differentials = clone_module_homomorphisms(py, differentials);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            let category = qpa_gbnp::cat_of_right_algebra_modules(module.right_acting_algebra())
                .map_err(invalid)?;
            let kind = runtime::Wrap::wrap(
                qpa_gbnp::finite_complex(
                    category,
                    Integer::from(base_degree),
                    runtime::extract_all(differentials)?,
                )
                .map_err(invalid)?,
            );
            Ok(PyModuleComplex { kind })
        })
    }

    fn projective_resolution_complex(&self) -> PyResult<PyModuleComplex> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleComplexKind,
            |module| qpa_gbnp::projective_resolution_complex(module).map_err(invalid)?
        );
        Ok(PyModuleComplex { kind })
    }

    fn injective_resolution_complex(&self) -> PyResult<PyModuleComplex> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleComplexKind,
            |module| qpa_gbnp::injective_resolution_complex(module).map_err(invalid)?
        );
        Ok(PyModuleComplex { kind })
    }

    #[getter]
    fn cached_projective_dimension(&self) -> Option<usize> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::proj_dimension(module)
        })
    }

    /// Compute projective dimension up to `bound`, or return `None` if unresolved.
    fn projective_dimension(&self, bound: i64) -> PyResult<Option<usize>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::proj_dimension_of_module(module, &bound).map_err(invalid)
        })
    }

    #[getter]
    fn cached_injective_dimension(&self) -> Option<usize> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::inj_dimension(module)
        })
    }

    fn injective_dimension(&self, bound: i64) -> PyResult<Option<usize>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::inj_dimension_of_module(module, &bound).map_err(invalid)
        })
    }

    fn dominant_dimension(&self, bound: i64) -> PyResult<Option<String>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            dominant_dimension_text(
                qpa_gbnp::dominant_dimension_of_module(module, &bound).map_err(invalid)?,
            )
        ))
    }

    fn is_nth_syzygy(&self, degree: i64) -> PyResult<bool> {
        let degree = Integer::from(degree);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_nth_syzygy(module, &degree).map_err(invalid)
        })
    }

    fn omega_period(&self, bound: i64) -> PyResult<Option<String>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_gbnp::is_omega_periodic(module, &bound)
                .map_err(invalid)?
                .map(|degree| degree.to_string())
        ))
    }

    fn tau_period(&self, bound: i64) -> PyResult<Option<String>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_gbnp::is_tau_periodic(module, &bound)
                .map_err(invalid)?
                .map(|degree| degree.to_string())
        ))
    }

    /// Compute QPA's degree-one Ext result with coefficients in the argument.
    fn ext1(&self, coefficients_in: &Self) -> PyResult<PyExtOverAlgebraResult> {
        let kind = map_pair_ref!(
            &self.kind,
            &coefficients_in.kind,
            MatrixModuleKind => ExtOverAlgebraKind,
            |module, coefficients_in| qpa_gbnp::ext_over_algebra(module, coefficients_in)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyExtOverAlgebraResult { kind })
    }

    /// Compute the degree-one Ext dimension without quotient representatives.
    fn ext1_dimension(&self, coefficients_in: &Self) -> PyResult<usize> {
        dispatch_pair_ref!(
            &self.kind,
            &coefficients_in.kind,
            MatrixModuleKind,
            |module, target| qpa_gbnp::ext_dimension_over_algebra(module, target).map_err(invalid),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    /// Compute Ext dimensions in input order, reusing this source's cover and syzygy.
    fn ext1_dimensions(&self, coefficients_in: Vec<PyRef<'_, Self>>) -> PyResult<Vec<usize>> {
        let targets = coefficients_in
            .iter()
            .map(|target| target.kind.clone())
            .collect();
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| batch_ext_dimensions(
            module, targets
        ))
    }

    fn ext_algebra_generators(
        &self,
        maximum_degree: usize,
    ) -> PyResult<PyExtAlgebraGeneratorsResult> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ExtAlgebraGeneratorsKind,
            |module| qpa_gbnp::ext_algebra_generators(module, maximum_degree).map_err(invalid)?
        );
        Ok(PyExtAlgebraGeneratorsResult { kind })
    }

    fn partial_iyama_generator(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_gbnp::partial_iyama_generator(module).map_err(invalid)?,
        )))
    }

    fn iyama_generator(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_gbnp::iyama_generator(module).map_err(invalid)?,
        )))
    }

    fn complexity(&self, degree: i64) -> PyResult<usize> {
        let degree = Integer::from(degree);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::complexity_of_module(module, &degree).map_err(invalid)
        })
    }

    #[getter]
    fn is_rigid(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_rigid_module(module).map_err(invalid)
        })
    }

    #[getter]
    fn is_tau_rigid(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_tau_rigid_module(module).map_err(invalid)
        })
    }

    #[getter]
    fn is_exceptional(&self) -> PyResult<Option<bool>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_exceptional_module(module).map_err(invalid)
        })
    }

    fn is_n_rigid(&self, bound: i64) -> PyResult<bool> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::n_rigid_module(module, &bound).map_err(invalid)
        })
    }

    #[getter]
    fn cached_is_tilting(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_tilting_module(module)
        })
    }

    fn finite_coresolution_in_add(
        &self,
        module: &Self,
        bound: i64,
    ) -> PyResult<Option<PyAddMCoresolution>> {
        let bound = Integer::from(bound);
        dispatch_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind,
            |source, module| Ok(qpa_gbnp::have_finite_coresolution_in_add_m(
                source, module, &bound,
            )
            .map_err(invalid)?
            .map(|result| PyAddMCoresolution {
                kind: runtime::Wrap::wrap(result),
            })),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn tilting(&self, bound: i64) -> PyResult<Option<PyTiltingModuleResult>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_gbnp::tilting_module(module, &bound)
                .map_err(invalid)?
                .map(|result| PyTiltingModuleResult {
                    kind: runtime::Wrap::wrap(result),
                })
        ))
    }

    #[getter]
    fn cached_is_cotilting(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_cotilting_module(module)
        })
    }

    fn finite_resolution_in_add(
        &self,
        module: &Self,
        bound: i64,
    ) -> PyResult<Option<PyFiniteAddMResolutionResult>> {
        let bound = Integer::from(bound);
        dispatch_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind,
            |target, module| Ok(
                qpa_gbnp::have_finite_resolution_in_add_m(target, module, &bound,)
                    .map_err(invalid)?
                    .map(|result| PyFiniteAddMResolutionResult {
                        kind: runtime::Wrap::wrap(result),
                    })
            ),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn cotilting(&self, bound: i64) -> PyResult<Option<PyCotiltingModuleResult>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_gbnp::cotilting_module(module, &bound)
                .map_err(invalid)?
                .map(|result| PyCotiltingModuleResult {
                    kind: runtime::Wrap::wrap(result),
                })
        ))
    }

    fn all_complements(&self, complement: &Self) -> PyResult<PyAllComplementsResult> {
        let kind = map_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind => AllComplementsResultKind,
            |module, complement| qpa_gbnp::all_complements_of_almost_complete_tilting_module(
                module, complement,
            )
            .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyAllComplementsResult { kind })
    }

    fn all_cotilting_complements(&self, complement: &Self) -> PyResult<PyAllComplementsResult> {
        let kind = map_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind => AllComplementsResultKind,
            |module, complement| qpa_gbnp::all_complements_of_almost_complete_cotilting_module(
                module, complement,
            )
            .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyAllComplementsResult { kind })
    }

    #[getter]
    fn cached_faithful_dimension(&self) -> Option<String> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .cached_faithful_dimension()
            .map(
                |dimension| dimension.map_or_else(|| "infinity".into(), |value| value.to_string())
            ))
    }

    fn faithful_dimension(&self) -> PyResult<String> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            faithful_dimension_result_text(qpa_gbnp::faithful_dimension(module).map_err(invalid)?)
        ))
    }

    fn number_of_tilting_complements(&self) -> PyResult<String> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            number_of_complements_result_text(
                qpa_gbnp::number_of_complements_of_almost_complete_tilting_module(module)
                    .map_err(invalid)?,
            )
        ))
    }

    fn number_of_cotilting_complements(&self) -> PyResult<String> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            number_of_complements_result_text(
                qpa_gbnp::number_of_complements_of_almost_complete_cotilting_module(module)
                    .map_err(invalid)?,
            )
        ))
    }

    fn left_tilting_complement_mutation(&self, complement: &Self) -> PyResult<Option<Self>> {
        dispatch_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind,
            |module, complement| Ok(qpa_gbnp::left_mutation_of_tilting_module_complement(
                module, complement
            )
            .map_err(invalid)?
            .map(wrap_module)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn right_tilting_complement_mutation(&self, complement: &Self) -> PyResult<Option<Self>> {
        dispatch_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind,
            |module, complement| Ok(qpa_gbnp::right_mutation_of_tilting_module_complement(
                module, complement
            )
            .map_err(invalid)?
            .map(wrap_module)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn left_cotilting_complement_mutation(&self, complement: &Self) -> PyResult<Option<Self>> {
        dispatch_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind,
            |module, complement| Ok(qpa_gbnp::left_mutation_of_cotilting_module_complement(
                module, complement
            )
            .map_err(invalid)?
            .map(wrap_module)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn right_cotilting_complement_mutation(&self, complement: &Self) -> PyResult<Option<Self>> {
        dispatch_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind,
            |module, complement| Ok(qpa_gbnp::right_mutation_of_cotilting_module_complement(
                module, complement
            )
            .map_err(invalid)?
            .map(wrap_module)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn right_approximation_by_add(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &target.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module, target| qpa_gbnp::right_approximation_by_add_m(module, target)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    #[staticmethod]
    fn right_approximation_by_add_list(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
        target: &Self,
    ) -> PyResult<PyModuleHomomorphism> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        match &target.kind {
            MatrixModuleKind::Rational(target) => {
                right_approximation_by_add_list_for::<RationalField>(modules, target)
            }
            MatrixModuleKind::Binary(target) => {
                right_approximation_by_add_list_for::<Z2>(modules, target)
            }
            MatrixModuleKind::Prime(target) => {
                right_approximation_by_add_list_for::<Zp>(modules, target)
            }
            MatrixModuleKind::Extension(target) => {
                right_approximation_by_add_list_for::<GaloisField>(modules, target)
            }
        }
    }

    fn minimal_right_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &target.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module, target| qpa_gbnp::minimal_right_approximation(module, target)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_right_add_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        self.minimal_right_approximation(target)
    }

    fn left_approximation_by_add(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, module| qpa_gbnp::left_approximation_by_add_m(source, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_left_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, module| qpa_gbnp::minimal_left_approximation(source, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_left_add_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        self.minimal_left_approximation(module)
    }

    fn right_fac_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &target.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module, target| qpa_core::right_fac_approximation(module, target).map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn right_fac_m_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        self.right_fac_approximation(target)
    }

    fn minimal_right_fac_m_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        self.right_fac_approximation(target)
    }

    fn left_fac_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, module| qpa_gbnp::left_fac_m_approximation(source, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_left_fac_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, module| qpa_gbnp::minimal_left_fac_m_approximation(source, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn right_sub_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &target.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module, target| qpa_gbnp::right_sub_m_approximation(module, target)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_right_sub_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &target.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module, target| qpa_gbnp::minimal_right_sub_m_approximation(module, target)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn left_sub_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, module| qpa_gbnp::left_sub_m_approximation(source, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_left_sub_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        self.left_sub_approximation(module)
    }

    fn right_approximation_by_perp(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |cotilting, module| qpa_gbnp::right_approximation_by_perp_t(cotilting, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn left_approximation_by_add_t_hat(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |cotilting, module| qpa_gbnp::left_approximation_by_add_t_hat(cotilting, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    #[staticmethod]
    fn radical_right_approximation_by_add(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
        one_based_index: i64,
    ) -> PyResult<PyModuleHomomorphism> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        let index = Integer::from(one_based_index);
        match modules.first() {
            Some(MatrixModuleKind::Rational(_)) => {
                radical_right_approximation_by_add_for::<RationalField>(modules, &index)
            }
            Some(MatrixModuleKind::Binary(_)) => {
                radical_right_approximation_by_add_for::<Z2>(modules, &index)
            }
            Some(MatrixModuleKind::Prime(_)) => {
                radical_right_approximation_by_add_for::<Zp>(modules, &index)
            }
            Some(MatrixModuleKind::Extension(_)) => {
                radical_right_approximation_by_add_for::<GaloisField>(modules, &index)
            }
            // The Rust translation rejects the one-based index before it needs
            // a coefficient field, so any monomorphization preserves QPA's
            // empty-list error instead of inventing a binding-only shortcut.
            None => radical_right_approximation_by_add_for::<RationalField>(modules, &index),
        }
    }

    #[staticmethod]
    fn projective_resolution_of_simple_module_over_endo(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
        one_based_index: i64,
        length: i64,
    ) -> PyResult<PySimpleEndoProjectiveResolutionResult> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        let index = Integer::from(one_based_index);
        let length = Integer::from(length);
        match modules.first() {
            Some(MatrixModuleKind::Rational(_)) | None => {
                simple_endo_projective_resolution_for::<RationalField>(modules, &index, &length)
            }
            Some(MatrixModuleKind::Binary(_)) => {
                simple_endo_projective_resolution_for::<Z2>(modules, &index, &length)
            }
            Some(MatrixModuleKind::Prime(_)) => {
                simple_endo_projective_resolution_for::<Zp>(modules, &index, &length)
            }
            Some(MatrixModuleKind::Extension(_)) => {
                simple_endo_projective_resolution_for::<GaloisField>(modules, &index, &length)
            }
        }
    }

    #[pyo3(signature = (side = "ending"))]
    fn almost_split_sequence(
        &self,
        side: &str,
    ) -> PyResult<Option<(PyModuleHomomorphism, PyModuleHomomorphism)>> {
        let side = match side {
            "ending" | "r" => qpa_gbnp::AlmostSplitSide::Ending,
            "starting" | "l" => qpa_gbnp::AlmostSplitSide::Starting,
            _ => {
                return Err(invalid(
                    "almost-split side must be 'ending'/'r' or 'starting'/'l'",
                ));
            }
        };
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism_pair(
                qpa_gbnp::almost_split_sequence_from_side(module, side).map_err(invalid)?,
            )
        ))
    }

    fn almost_split_sequence_in_perp(
        &self,
        module: &Self,
    ) -> PyResult<Option<(PyModuleHomomorphism, PyModuleHomomorphism)>> {
        dispatch_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind,
            |cotilting, module| Ok(wrap_module_homomorphism_pair(
                qpa_gbnp::almost_split_sequence_in_perp_t(cotilting, module).map_err(invalid)?,
            )),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn irreducible_morphisms_ending_in(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(
                qpa_gbnp::irreducible_morphisms_ending_in(module).map_err(invalid)?,
            )
        ))
    }

    fn irreducible_morphisms_starting_in(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(
                qpa_gbnp::irreducible_morphisms_starting_in(module).map_err(invalid)?,
            )
        ))
    }

    fn predecessors(&self, bound: i64) -> PyResult<PyPredecessorsOfModuleResult> {
        let bound = Integer::from(bound);
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => PredecessorsOfModuleResultKind,
            |module| qpa_gbnp::predecessors_of_module(module, &bound).map_err(invalid)?
        );
        Ok(PyPredecessorsOfModuleResult { kind })
    }

    fn dual(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::dual_of_module(module).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn transpose(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::transpose_of_module(module).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[pyo3(signature = (exponent=1))]
    fn d_tr(&self, exponent: i64) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::d_tr_iterated(module, exponent).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[pyo3(signature = (exponent=1))]
    fn tr_d(&self, exponent: i64) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::tr_d_iterated(module, exponent).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn star(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::star_of_module(module).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn identity_to_double_star(&self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module| qpa_gbnp::from_identity_to_double_star_homomorphism(module)
                .map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    fn nakayama(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::nakayama_functor_of_module(module).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn opposite_nakayama(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::opposite_nakayama_functor_of_module(module).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            MatrixModuleKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn view(&self) -> String {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module.to_string())
    }

    fn print_obj(&self) -> String {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .qpa_print_string())
    }

    fn display(&self) -> PyResult<String> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .qpa_display_string_with(|value| module
                .right_acting_algebra()
                .field()
                .scalar_text(value))
            .map_err(invalid))
    }

    fn __eq__(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            MatrixModuleKind,
            |left, right| left == right,
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .qpa_print_string())
    }
}

fn matrix_shapes<F: numerica::domains::Field>(matrices: &[DenseMatrix<F>]) -> Vec<(usize, usize)> {
    matrices
        .iter()
        .map(|matrix| (matrix.nrows(), matrix.ncols()))
        .collect()
}
