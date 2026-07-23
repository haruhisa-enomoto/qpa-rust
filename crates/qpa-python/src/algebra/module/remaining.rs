//! Final source-faithful QPA module operations exposed to Python.

use std::num::NonZeroUsize;

use numerica::domains::finite_field::FiniteFieldCore as _;
use pyo3::{pymethods, types::PyList};
use qpa_core::{
    Arrow, BasedRightAlgebraModule, IndexedPathAlgebraModuleArrowAction,
    NamedPathAlgebraModuleArrowAction, PathAlgebraModuleArrowAction, PathAlgebraModuleArrowActions,
    ProbabilisticModuleDecomposition,
};
use qpa_gbnp::TensorProductOfModulesResult;
use rand::{SeedableRng, rngs::StdRng};

use super::super::conversion::QuiverAlgebraKind;
#[allow(clippy::wildcard_imports)]
// This binding domain intentionally shares its parent module runtime enums.
use super::*;
use crate::ScalarKind;

#[derive(Clone)]
enum BasedRightAlgebraModuleKind {
    Rational(BasedRightAlgebraModule<RationalField>),
    Binary(BasedRightAlgebraModule<Z2>),
    Prime(BasedRightAlgebraModule<Zp>),
    Extension(BasedRightAlgebraModule<GaloisField>),
}

/// A finite-dimensional based right algebra module accepted by QPA's generic
/// conversion to a path-algebra matrix module.
#[pyclass(name = "BasedRightAlgebraModule", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyBasedRightAlgebraModule {
    kind: BasedRightAlgebraModuleKind,
}

#[derive(Clone)]
enum ModuleTensorProductKind {
    Rational(TensorProductOfModulesResult<RationalField>),
    Binary(TensorProductOfModulesResult<Z2>),
    Prime(TensorProductOfModulesResult<Zp>),
    Extension(TensorProductOfModulesResult<GaloisField>),
}

enum PythonInferredArrowAction {
    Matrix(Vec<Vec<Py<PyScalar>>>),
    ZeroSpaceDimensions { source: usize, target: usize },
}

/// QPA's vector space and elementary-tensor closure returned by
/// `TensorProductOfModules`.
#[pyclass(name = "ModuleTensorProduct", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleTensorProduct {
    kind: ModuleTensorProductKind,
}

/// QPA's `false`/`fail` Bongartz result together with its printed message.
#[pyclass(name = "BongartzTestDecision", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyBongartzTestDecision {
    value: Option<bool>,
    messages: Vec<String>,
}

/// QPA finite-type value, source branch, and printed diagnostics.
#[pyclass(name = "FiniteTypeDecision", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyFiniteTypeDecision {
    value: Option<bool>,
    reason: String,
    messages: Vec<String>,
}

macro_rules! impl_remaining_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<BasedRightAlgebraModule<$field>> for BasedRightAlgebraModuleKind {
            fn wrap(value: BasedRightAlgebraModule<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<TensorProductOfModulesResult<$field>> for ModuleTensorProductKind {
            fn wrap(value: TensorProductOfModulesResult<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_remaining_runtime!(RationalField, Rational);
impl_remaining_runtime!(Z2, Binary);
impl_remaining_runtime!(Zp, Prime);
impl_remaining_runtime!(GaloisField, Extension);

fn wrap_modules<F: runtime::PythonField>(modules: Vec<PathAlgebraModule<F>>) -> Vec<PyMatrixModule>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    modules.into_iter().map(wrap_module).collect()
}

fn square_matrices<F: runtime::PythonField>(
    py: Python<'_>,
    field: &F,
    dimension: usize,
    matrices: Vec<Vec<Vec<Py<PyScalar>>>>,
) -> PyResult<Vec<DenseMatrix<F>>> {
    matrices
        .into_iter()
        .map(|rows| {
            if rows.len() != dimension || rows.iter().any(|row| row.len() != dimension) {
                return Err(invalid(format!(
                    "based-module action has shape {}x{}, expected {dimension}x{dimension}",
                    rows.len(),
                    rows.first().map_or(0, Vec::len)
                )));
            }
            let entries = rows
                .into_iter()
                .flatten()
                .map(|entry| field.scalar_value(&entry.borrow(py)))
                .collect::<PyResult<Vec<_>>>()?;
            DenseMatrix::from_flat(dimension, dimension, entries, field.clone()).map_err(invalid)
        })
        .collect()
}

fn based_module_for<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: QuiverAlgebra<F>,
    dimension: usize,
    vertex_actions: Vec<Vec<Vec<Py<PyScalar>>>>,
    arrow_actions: Vec<Vec<Vec<Py<PyScalar>>>>,
) -> PyResult<PyBasedRightAlgebraModule>
where
    BasedRightAlgebraModuleKind: runtime::Wrap<BasedRightAlgebraModule<F>>,
{
    let field = algebra.field().clone();
    Ok(PyBasedRightAlgebraModule {
        kind: runtime::Wrap::wrap(
            BasedRightAlgebraModule::new(
                algebra,
                square_matrices(py, &field, dimension, vertex_actions)?,
                square_matrices(py, &field, dimension, arrow_actions)?,
            )
            .map_err(invalid)?,
        ),
    })
}

fn all_modules_at_most_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
    bound: &Integer,
) -> PyResult<Vec<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    Ok(wrap_modules(
        qpa_gbnp::all_modules_of_length_at_most(algebra, bound).map_err(invalid)?,
    ))
}

fn all_indec_modules_at_most_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
    bound: &Integer,
) -> PyResult<Vec<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    Ok(wrap_modules(
        qpa_gbnp::all_indec_modules_of_length_at_most(algebra, bound).map_err(invalid)?,
    ))
}

fn modules_plus_one_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
) -> PyResult<Vec<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>> + runtime::Wrap<PathAlgebraModule<F>>,
{
    Ok(wrap_modules(
        qpa_gbnp::all_modules_of_length_plus_one(&runtime::extract_all(modules)?)
            .map_err(invalid)?,
    ))
}

fn basis_of_projectives_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Option<Vec<Vec<Vec<PyQuiverAlgebraElement>>>>>
where
    tensor::QuiverAlgebraElementKind: runtime::Wrap<QuiverAlgebraElement<F>>,
{
    Ok(qpa_gbnp::basis_of_projectives(algebra)
        .map_err(invalid)?
        .map(|projectives| {
            projectives
                .into_iter()
                .map(|components| {
                    components
                        .into_iter()
                        .map(|elements| {
                            elements
                                .into_iter()
                                .map(tensor::wrap_quiver_element)
                                .collect()
                        })
                        .collect()
                })
                .collect()
        }))
}

fn random_module_for<F: runtime::PythonField + qpa_gbnp::QpaRandomScalar>(
    algebra: &QuiverAlgebra<F>,
    source_summands: NonZeroUsize,
    range_summands: NonZeroUsize,
    seed: u64,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    let mut rng = StdRng::seed_from_u64(seed);
    Ok(
        qpa_gbnp::random_module(algebra, source_summands, range_summands, &mut rng)
            .map_err(invalid)?
            .map(wrap_module),
    )
}

fn explicit_module_for<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: QuiverAlgebra<F>,
    dimension_vector: Vec<usize>,
    actions: Vec<(String, Vec<Vec<Py<PyScalar>>>)>,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    let field = algebra.field().clone();
    let actions = actions
        .into_iter()
        .map(|(label, rows)| {
            Ok(NamedArrowMatrix::new(
                label,
                extract_matrix(py, &field, rows)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(
        qpa_core::right_module_over_path_algebra(algebra, dimension_vector, &actions)
            .map_err(invalid)?
            .map(wrap_module),
    )
}

fn inferred_checked_for<F: runtime::PythonField>(
    algebra: QuiverAlgebra<F>,
    actions: PathAlgebraModuleArrowActions<F>,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    Ok(
        qpa_core::right_module_over_path_algebra_infer_dimension(algebra, &actions)
            .map_err(invalid)?
            .map(wrap_module),
    )
}

fn inferred_nc_for<F: runtime::PythonField>(
    algebra: QuiverAlgebra<F>,
    actions: PathAlgebraModuleArrowActions<F>,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    Ok(
        qpa_core::right_module_over_path_algebra_infer_dimension_nc(algebra, &actions)
            .map_err(invalid)?
            .map(wrap_module),
    )
}

fn unchecked_prime_action_field<'a>(
    py: Python<'_>,
    scalars: impl IntoIterator<Item = &'a Py<PyScalar>>,
    default_characteristic: u32,
) -> PyResult<Zp> {
    let mut characteristic = None;
    for scalar in scalars {
        let ScalarKind::Prime {
            characteristic: scalar_characteristic,
            ..
        } = scalar.borrow(py).value
        else {
            return Err(invalid(
                "a prime-field NC module requires prime-field action entries",
            ));
        };
        match characteristic {
            Some(previous) if previous != scalar_characteristic => {
                return Err(invalid(
                    "all NC module action matrices must use one coefficient field",
                ));
            }
            None => characteristic = Some(scalar_characteristic),
            Some(_) => {}
        }
    }
    Ok(Zp::new(characteristic.unwrap_or(default_characteristic)))
}

fn parse_inferred_action(py: Python<'_>, action: Py<PyAny>) -> PyResult<PythonInferredArrowAction> {
    let action = action.bind(py);
    if let Ok(rows) = action.extract::<Vec<Vec<Py<PyScalar>>>>() {
        return Ok(PythonInferredArrowAction::Matrix(rows));
    }
    if let Ok(dimensions) = action.extract::<Vec<usize>>()
        && let [source, target] = dimensions.as_slice()
    {
        return Ok(PythonInferredArrowAction::ZeroSpaceDimensions {
            source: *source,
            target: *target,
        });
    }
    Err(invalid(
        "an inferred arrow action must be a scalar matrix or a two-entry zero-space dimension list",
    ))
}

fn inferred_action_for<F: runtime::PythonField>(
    py: Python<'_>,
    field: &F,
    action: PythonInferredArrowAction,
) -> PyResult<PathAlgebraModuleArrowAction<F>> {
    match action {
        PythonInferredArrowAction::Matrix(rows) => Ok(PathAlgebraModuleArrowAction::matrix(
            extract_matrix(py, field, rows)?,
        )),
        PythonInferredArrowAction::ZeroSpaceDimensions { source, target } => Ok(
            PathAlgebraModuleArrowAction::zero_space_dimensions(source, target),
        ),
    }
}

fn named_nc_for<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: QuiverAlgebra<F>,
    field: &F,
    actions: Vec<(String, PythonInferredArrowAction)>,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    let actions = actions
        .into_iter()
        .map(|(label, action)| {
            Ok(NamedPathAlgebraModuleArrowAction::new(
                label,
                inferred_action_for(py, field, action)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    inferred_nc_for(algebra, PathAlgebraModuleArrowActions::Named(actions))
}

fn named_checked_for<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: QuiverAlgebra<F>,
    field: &F,
    actions: Vec<(String, PythonInferredArrowAction)>,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    let actions = actions
        .into_iter()
        .map(|(label, action)| {
            Ok(NamedPathAlgebraModuleArrowAction::new(
                label,
                inferred_action_for(py, field, action)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    inferred_checked_for(algebra, PathAlgebraModuleArrowActions::Named(actions))
}

fn indexed_nc_for<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: QuiverAlgebra<F>,
    field: &F,
    actions: Vec<(Arrow, PythonInferredArrowAction)>,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    let actions = actions
        .into_iter()
        .map(|(arrow, action)| {
            if arrow.quiver() != algebra.quiver() {
                return Err(invalid(
                    "an indexed NC arrow does not belong to the algebra quiver",
                ));
            }
            Ok(IndexedPathAlgebraModuleArrowAction::new(
                arrow.index(),
                inferred_action_for(py, field, action)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    inferred_nc_for(algebra, PathAlgebraModuleArrowActions::Indexed(actions))
}

fn indexed_checked_for<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: QuiverAlgebra<F>,
    field: &F,
    actions: Vec<(Arrow, PythonInferredArrowAction)>,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    let actions = actions
        .into_iter()
        .map(|(arrow, action)| {
            if arrow.quiver() != algebra.quiver() {
                return Err(invalid(
                    "an indexed inferred arrow does not belong to the algebra quiver",
                ));
            }
            Ok(IndexedPathAlgebraModuleArrowAction::new(
                arrow.index(),
                inferred_action_for(py, field, action)?,
            ))
        })
        .collect::<PyResult<Vec<_>>>()?;
    inferred_checked_for(algebra, PathAlgebraModuleArrowActions::Indexed(actions))
}

fn inferred_action_scalars<'a>(
    actions: impl IntoIterator<Item = &'a PythonInferredArrowAction>,
) -> Vec<&'a Py<PyScalar>> {
    actions
        .into_iter()
        .flat_map(|action| match action {
            PythonInferredArrowAction::Matrix(rows) => rows.iter().flatten().collect::<Vec<_>>(),
            PythonInferredArrowAction::ZeroSpaceDimensions { .. } => Vec::new(),
        })
        .collect()
}

fn decomposition_probabilistic_for<F: runtime::PythonField>(
    py: Python<'_>,
    module: &PathAlgebraModule<F>,
    homomorphisms: Vec<ModuleHomomorphismKind>,
) -> PyResult<Py<PyAny>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>,
{
    match qpa_core::decompose_module_probabilistic(&runtime::extract_all(homomorphisms)?, module)
        .map_err(invalid)?
    {
        ProbabilisticModuleDecomposition::Unsplit(module) => {
            Ok(Py::new(py, wrap_module(module))?.into_any())
        }
        ProbabilisticModuleDecomposition::Split(modules) => Ok(PyList::new(
            py,
            modules.into_iter().map(wrap_module).collect::<Vec<_>>(),
        )?
        .into_any()
        .unbind()),
    }
}

fn tensor_product_modules_for<F: runtime::PythonField>(
    left: &PathAlgebraModule<F>,
    right: &PathAlgebraModule<F>,
) -> PyResult<PyModuleTensorProduct>
where
    ModuleTensorProductKind: runtime::Wrap<TensorProductOfModulesResult<F>>,
{
    Ok(PyModuleTensorProduct {
        kind: runtime::Wrap::wrap(
            qpa_gbnp::tensor_product_of_modules(left, right).map_err(invalid)?,
        ),
    })
}

#[pyfunction(name = "complement_in_full_row_space")]
pub(crate) fn py_complement_in_full_row_space(
    py: Python<'_>,
    rows: Vec<Vec<Py<PyScalar>>>,
    field: PyRef<'_, PyField>,
) -> PyResult<Vec<Vec<PyScalar>>> {
    match field.kind {
        FieldKind::Rational => Ok(matrix_rows(
            &qpa_core::complement_in_full_row_space(&extract_matrix(py, &Q, rows)?)
                .map_err(invalid)?,
        )),
        FieldKind::Binary => Ok(matrix_rows(
            &qpa_core::complement_in_full_row_space(&extract_matrix(py, &Z2, rows)?)
                .map_err(invalid)?,
        )),
        FieldKind::Prime(characteristic) => {
            let field = Zp::new(characteristic);
            Ok(matrix_rows(
                &qpa_core::complement_in_full_row_space(&extract_matrix(py, &field, rows)?)
                    .map_err(invalid)?,
            ))
        }
        FieldKind::Extension { .. } => {
            let extension = field
                .kind
                .galois_field()
                .ok_or_else(|| invalid("expected an extension field"))?;
            Ok(matrix_rows(
                &qpa_core::complement_in_full_row_space(&extract_matrix(py, &extension, rows)?)
                    .map_err(invalid)?,
            ))
        }
    }
}

#[pymethods]
impl PyBongartzTestDecision {
    #[getter]
    fn value(&self) -> Option<bool> {
        self.value
    }

    #[getter]
    fn messages(&self) -> Vec<String> {
        self.messages.clone()
    }
}

#[pymethods]
impl PyFiniteTypeDecision {
    #[getter]
    const fn value(&self) -> Option<bool> {
        self.value
    }

    #[getter]
    fn reason(&self) -> String {
        self.reason.clone()
    }

    #[getter]
    fn messages(&self) -> Vec<String> {
        self.messages.clone()
    }
}

#[pymethods]
impl PyBasedRightAlgebraModule {
    #[staticmethod]
    fn from_actions(
        py: Python<'_>,
        algebra: &PyQuiverAlgebra,
        dimension: usize,
        vertex_actions: Vec<Vec<Vec<Py<PyScalar>>>>,
        arrow_actions: Vec<Vec<Vec<Py<PyScalar>>>>,
    ) -> PyResult<Self> {
        match &algebra.kind {
            conversion::QuiverAlgebraKind::Rational(algebra) => based_module_for(
                py,
                algebra.clone(),
                dimension,
                vertex_actions,
                arrow_actions,
            ),
            conversion::QuiverAlgebraKind::Binary(algebra) => based_module_for(
                py,
                algebra.clone(),
                dimension,
                vertex_actions,
                arrow_actions,
            ),
            conversion::QuiverAlgebraKind::Prime(algebra) => based_module_for(
                py,
                algebra.clone(),
                dimension,
                vertex_actions,
                arrow_actions,
            ),
            conversion::QuiverAlgebraKind::Extension(algebra) => based_module_for(
                py,
                algebra.clone(),
                dimension,
                vertex_actions,
                arrow_actions,
            ),
        }
    }

    #[staticmethod]
    fn from_matrix_module(module: &PyMatrixModule) -> PyResult<Self> {
        let kind = map_ref!(
            &module.kind,
            MatrixModuleKind => BasedRightAlgebraModuleKind,
            |module| BasedRightAlgebraModule::from_path_algebra_module(module).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, BasedRightAlgebraModuleKind, |module| {
            conversion::wrap_quiver_algebra(module.algebra().clone())
        })
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, BasedRightAlgebraModuleKind, |module| module
            .dimension())
    }

    #[getter]
    fn vertex_actions(&self) -> Vec<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, BasedRightAlgebraModuleKind, |module| module
            .vertex_actions()
            .iter()
            .map(matrix_rows)
            .collect())
    }

    #[getter]
    fn arrow_actions(&self) -> Vec<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, BasedRightAlgebraModuleKind, |module| module
            .arrow_actions()
            .iter()
            .map(matrix_rows)
            .collect())
    }

    fn to_path_algebra_matrix_module(&self) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, BasedRightAlgebraModuleKind, |module| Ok(
            wrap_module(
                qpa_core::right_algebra_module_to_path_algebra_mat_module(module)
                    .map_err(invalid)?,
            )
        ))
    }
}

#[pymethods]
impl PyModuleTensorProduct {
    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, ModuleTensorProductKind, |tensor| tensor
            .dimension())
    }

    #[getter]
    fn ambient_dimension(&self) -> usize {
        dispatch_ref!(&self.kind, ModuleTensorProductKind, |tensor| tensor
            .ambient_dimension())
    }

    #[getter]
    fn relation_basis(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, ModuleTensorProductKind, |tensor| matrix_rows(
            tensor.relation_basis()
        ))
    }

    #[getter]
    fn basis_lifts(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, ModuleTensorProductKind, |tensor| matrix_rows(
            tensor.basis_lifts()
        ))
    }

    #[getter]
    fn projection(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, ModuleTensorProductKind, |tensor| matrix_rows(
            tensor.projection()
        ))
    }

    #[getter]
    fn basis_vectors(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, ModuleTensorProductKind, |tensor| scalar_rows(
            tensor.projection().field(),
            &tensor.basis_vectors(),
        ))
    }

    fn elementary_tensor(
        &self,
        left: &PyModuleElement,
        right: &PyModuleElement,
    ) -> PyResult<Vec<PyScalar>> {
        match (&self.kind, &left.kind, &right.kind) {
            (
                ModuleTensorProductKind::Rational(tensor),
                ModuleElementKind::Rational(left),
                ModuleElementKind::Rational(right),
            ) => Ok(tensor
                .elementary_tensor(left, right)
                .map_err(invalid)?
                .iter()
                .map(|value| Q.wrap_scalar(value))
                .collect()),
            (
                ModuleTensorProductKind::Binary(tensor),
                ModuleElementKind::Binary(left),
                ModuleElementKind::Binary(right),
            ) => Ok(tensor
                .elementary_tensor(left, right)
                .map_err(invalid)?
                .iter()
                .map(|value| Z2.wrap_scalar(value))
                .collect()),
            (
                ModuleTensorProductKind::Prime(tensor),
                ModuleElementKind::Prime(left),
                ModuleElementKind::Prime(right),
            ) => {
                let field = tensor.projection().field();
                Ok(tensor
                    .elementary_tensor(left, right)
                    .map_err(invalid)?
                    .iter()
                    .map(|value| field.wrap_scalar(value))
                    .collect())
            }
            (
                ModuleTensorProductKind::Extension(tensor),
                ModuleElementKind::Extension(left),
                ModuleElementKind::Extension(right),
            ) => {
                let field = tensor.projection().field();
                Ok(tensor
                    .elementary_tensor(left, right)
                    .map_err(invalid)?
                    .iter()
                    .map(|value| field.wrap_scalar(value))
                    .collect())
            }
            _ => Err(invalid(
                "tensor product and module elements use different coefficient fields",
            )),
        }
    }
}

#[pymethods]
impl PyQuiverAlgebra {
    #[getter]
    fn loewy_length(&self) -> PyResult<usize> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            qpa_gbnp::loewy_length_of_algebra(algebra).map_err(invalid)
        })
    }

    fn basis_of_projectives(&self) -> PyResult<Option<Vec<Vec<Vec<PyQuiverAlgebraElement>>>>> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            basis_of_projectives_for(algebra)
        })
    }

    fn all_modules_of_length_at_most(
        &self,
        bound: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<PyMatrixModule>> {
        let bound = crate::field::parse_python_integer(bound)?;
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            all_modules_at_most_for(algebra, &bound)
        })
    }

    fn all_indec_modules_of_length_at_most(
        &self,
        bound: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<PyMatrixModule>> {
        let bound = crate::field::parse_python_integer(bound)?;
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            all_indec_modules_at_most_for(algebra, &bound)
        })
    }

    fn bongartz_test(&self, bound: i64) -> PyResult<PyBongartzTestDecision> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            let decision = qpa_gbnp::bongartz_test(algebra, bound).map_err(invalid)?;
            Ok(PyBongartzTestDecision {
                value: decision.value(),
                messages: decision
                    .messages()
                    .iter()
                    .map(ToString::to_string)
                    .collect(),
            })
        })
    }

    #[pyo3(signature = (source_summands, range_summands, seed=1))]
    fn random_module(
        &self,
        source_summands: usize,
        range_summands: usize,
        seed: u64,
    ) -> PyResult<Option<PyMatrixModule>> {
        let source_summands = NonZeroUsize::new(source_summands)
            .ok_or_else(|| invalid("RandomModule source summand count must be positive"))?;
        let range_summands = NonZeroUsize::new(range_summands)
            .ok_or_else(|| invalid("RandomModule range summand count must be positive"))?;
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            random_module_for(algebra, source_summands, range_summands, seed)
        })
    }

    fn right_module_over_path_algebra_with_dimensions(
        &self,
        py: Python<'_>,
        dimension_vector: Vec<usize>,
        actions: Vec<(String, Vec<Vec<Py<PyScalar>>>)>,
    ) -> PyResult<Option<PyMatrixModule>> {
        match &self.kind {
            QuiverAlgebraKind::Rational(algebra) => {
                explicit_module_for(py, algebra.clone(), dimension_vector, actions)
            }
            QuiverAlgebraKind::Binary(algebra) => {
                explicit_module_for(py, algebra.clone(), dimension_vector, actions)
            }
            QuiverAlgebraKind::Prime(algebra) => {
                explicit_module_for(py, algebra.clone(), dimension_vector, actions)
            }
            QuiverAlgebraKind::Extension(algebra) => {
                explicit_module_for(py, algebra.clone(), dimension_vector, actions)
            }
        }
    }

    fn right_module_over_path_algebra(
        &self,
        py: Python<'_>,
        actions: Vec<Vec<Vec<Py<PyScalar>>>>,
    ) -> PyResult<Option<PyMatrixModule>> {
        match &self.kind {
            QuiverAlgebraKind::Rational(algebra) => inferred_checked_for(
                algebra.clone(),
                PathAlgebraModuleArrowActions::Positional(
                    actions
                        .into_iter()
                        .map(|rows| extract_matrix(py, &Q, rows))
                        .collect::<PyResult<_>>()?,
                ),
            ),
            QuiverAlgebraKind::Binary(algebra) => inferred_checked_for(
                algebra.clone(),
                PathAlgebraModuleArrowActions::Positional(
                    actions
                        .into_iter()
                        .map(|rows| extract_matrix(py, &Z2, rows))
                        .collect::<PyResult<_>>()?,
                ),
            ),
            QuiverAlgebraKind::Prime(algebra) => inferred_checked_for(
                algebra.clone(),
                PathAlgebraModuleArrowActions::Positional(
                    actions
                        .into_iter()
                        .map(|rows| extract_matrix(py, algebra.field(), rows))
                        .collect::<PyResult<_>>()?,
                ),
            ),
            QuiverAlgebraKind::Extension(algebra) => inferred_checked_for(
                algebra.clone(),
                PathAlgebraModuleArrowActions::Positional(
                    actions
                        .into_iter()
                        .map(|rows| extract_matrix(py, algebra.field(), rows))
                        .collect::<PyResult<_>>()?,
                ),
            ),
        }
    }

    fn right_module_over_path_algebra_named(
        &self,
        py: Python<'_>,
        actions: Vec<(String, Py<PyAny>)>,
    ) -> PyResult<Option<PyMatrixModule>> {
        let actions = actions
            .into_iter()
            .map(|(label, action)| Ok((label, parse_inferred_action(py, action)?)))
            .collect::<PyResult<Vec<_>>>()?;
        match &self.kind {
            QuiverAlgebraKind::Rational(algebra) => {
                named_checked_for(py, algebra.clone(), &Q, actions)
            }
            QuiverAlgebraKind::Binary(algebra) => {
                named_checked_for(py, algebra.clone(), &Z2, actions)
            }
            QuiverAlgebraKind::Prime(algebra) => {
                named_checked_for(py, algebra.clone(), algebra.field(), actions)
            }
            QuiverAlgebraKind::Extension(algebra) => {
                named_checked_for(py, algebra.clone(), algebra.field(), actions)
            }
        }
    }

    fn right_module_over_path_algebra_indexed(
        &self,
        py: Python<'_>,
        actions: Vec<(Py<PyQuiverArrow>, Py<PyAny>)>,
    ) -> PyResult<Option<PyMatrixModule>> {
        let actions = actions
            .into_iter()
            .map(|(arrow, action)| {
                Ok((
                    arrow.borrow(py).inner.clone(),
                    parse_inferred_action(py, action)?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        match &self.kind {
            QuiverAlgebraKind::Rational(algebra) => {
                indexed_checked_for(py, algebra.clone(), &Q, actions)
            }
            QuiverAlgebraKind::Binary(algebra) => {
                indexed_checked_for(py, algebra.clone(), &Z2, actions)
            }
            QuiverAlgebraKind::Prime(algebra) => {
                indexed_checked_for(py, algebra.clone(), algebra.field(), actions)
            }
            QuiverAlgebraKind::Extension(algebra) => {
                indexed_checked_for(py, algebra.clone(), algebra.field(), actions)
            }
        }
    }

    fn right_module_over_path_algebra_nc(
        &self,
        py: Python<'_>,
        actions: Vec<Vec<Vec<Py<PyScalar>>>>,
    ) -> PyResult<Option<PyMatrixModule>> {
        match &self.kind {
            QuiverAlgebraKind::Rational(algebra) => inferred_nc_for(
                algebra.clone(),
                PathAlgebraModuleArrowActions::Positional(
                    actions
                        .into_iter()
                        .map(|rows| extract_matrix(py, &Q, rows))
                        .collect::<PyResult<_>>()?,
                ),
            ),
            QuiverAlgebraKind::Binary(algebra) => inferred_nc_for(
                algebra.clone(),
                PathAlgebraModuleArrowActions::Positional(
                    actions
                        .into_iter()
                        .map(|rows| extract_matrix(py, &Z2, rows))
                        .collect::<PyResult<_>>()?,
                ),
            ),
            QuiverAlgebraKind::Prime(algebra) => inferred_nc_for(algebra.clone(), {
                // QPA's NC method intentionally skips the coefficient-family
                // check. `Zp` carries its prime at runtime, so retain that
                // observable behavior for foreign odd-prime matrices.
                let action_field = unchecked_prime_action_field(
                    py,
                    actions.iter().flatten().flatten(),
                    algebra.field().get_prime(),
                )?;
                PathAlgebraModuleArrowActions::Positional(
                    actions
                        .into_iter()
                        .map(|rows| extract_matrix(py, &action_field, rows))
                        .collect::<PyResult<_>>()?,
                )
            }),
            QuiverAlgebraKind::Extension(algebra) => inferred_nc_for(
                algebra.clone(),
                PathAlgebraModuleArrowActions::Positional(
                    actions
                        .into_iter()
                        .map(|rows| extract_matrix(py, algebra.field(), rows))
                        .collect::<PyResult<_>>()?,
                ),
            ),
        }
    }

    fn right_module_over_path_algebra_nc_named(
        &self,
        py: Python<'_>,
        actions: Vec<(String, Py<PyAny>)>,
    ) -> PyResult<Option<PyMatrixModule>> {
        let actions = actions
            .into_iter()
            .map(|(label, action)| Ok((label, parse_inferred_action(py, action)?)))
            .collect::<PyResult<Vec<_>>>()?;
        match &self.kind {
            QuiverAlgebraKind::Rational(algebra) => named_nc_for(py, algebra.clone(), &Q, actions),
            QuiverAlgebraKind::Binary(algebra) => named_nc_for(py, algebra.clone(), &Z2, actions),
            QuiverAlgebraKind::Prime(algebra) => {
                let action_field = unchecked_prime_action_field(
                    py,
                    inferred_action_scalars(actions.iter().map(|(_, action)| action)),
                    algebra.field().get_prime(),
                )?;
                named_nc_for(py, algebra.clone(), &action_field, actions)
            }
            QuiverAlgebraKind::Extension(algebra) => {
                named_nc_for(py, algebra.clone(), algebra.field(), actions)
            }
        }
    }

    fn right_module_over_path_algebra_nc_indexed(
        &self,
        py: Python<'_>,
        actions: Vec<(Py<PyQuiverArrow>, Py<PyAny>)>,
    ) -> PyResult<Option<PyMatrixModule>> {
        let actions = actions
            .into_iter()
            .map(|(arrow, action)| {
                Ok((
                    arrow.borrow(py).inner.clone(),
                    parse_inferred_action(py, action)?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        match &self.kind {
            QuiverAlgebraKind::Rational(algebra) => {
                indexed_nc_for(py, algebra.clone(), &Q, actions)
            }
            QuiverAlgebraKind::Binary(algebra) => indexed_nc_for(py, algebra.clone(), &Z2, actions),
            QuiverAlgebraKind::Prime(algebra) => {
                let action_field = unchecked_prime_action_field(
                    py,
                    inferred_action_scalars(actions.iter().map(|(_, action)| action)),
                    algebra.field().get_prime(),
                )?;
                indexed_nc_for(py, algebra.clone(), &action_field, actions)
            }
            QuiverAlgebraKind::Extension(algebra) => {
                indexed_nc_for(py, algebra.clone(), algebra.field(), actions)
            }
        }
    }

    fn element_as_indec_projective(
        &self,
        element: &PyQuiverAlgebraElement,
    ) -> PyResult<PyModuleElement> {
        let kind = match (&self.kind, &element.kind) {
            (
                QuiverAlgebraKind::Rational(algebra),
                tensor::QuiverAlgebraElementKind::Rational(element),
            ) => ModuleElementKind::Rational(
                qpa_gbnp::element_in_v_a_as_element_in_indec_proj(algebra, element)
                    .map_err(invalid)?,
            ),
            (
                QuiverAlgebraKind::Binary(algebra),
                tensor::QuiverAlgebraElementKind::Binary(element),
            ) => ModuleElementKind::Binary(
                qpa_gbnp::element_in_v_a_as_element_in_indec_proj(algebra, element)
                    .map_err(invalid)?,
            ),
            (
                QuiverAlgebraKind::Prime(algebra),
                tensor::QuiverAlgebraElementKind::Prime(element),
            ) => ModuleElementKind::Prime(
                qpa_gbnp::element_in_v_a_as_element_in_indec_proj(algebra, element)
                    .map_err(invalid)?,
            ),
            (
                QuiverAlgebraKind::Extension(algebra),
                tensor::QuiverAlgebraElementKind::Extension(element),
            ) => ModuleElementKind::Extension(
                qpa_gbnp::element_in_v_a_as_element_in_indec_proj(algebra, element)
                    .map_err(invalid)?,
            ),
            _ => {
                return Err(invalid(
                    "algebra and element use different coefficient fields",
                ));
            }
        };
        Ok(PyModuleElement { kind })
    }
}

#[pymethods]
impl PyPathAlgebra {
    fn finite_type_decision(&self) -> PyResult<PyFiniteTypeDecision> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            let decision = qpa_gbnp::path_algebra_finite_type(algebra).map_err(invalid)?;
            Ok(PyFiniteTypeDecision {
                value: Some(decision.is_finite()),
                reason: format!("{:?}", decision.reason()),
                messages: decision
                    .messages()
                    .iter()
                    .map(ToString::to_string)
                    .collect(),
            })
        })
    }

    #[getter]
    fn is_finite_type_algebra(&self) -> PyResult<bool> {
        Ok(self
            .finite_type_decision()?
            .value
            .expect("path-algebra finite-type decisions are Boolean"))
    }
}

#[pymethods]
impl PyQuotientAlgebra {
    fn finite_type_decision(&self) -> PyResult<PyFiniteTypeDecision> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            let decision = qpa_gbnp::quotient_algebra_finite_type(algebra).map_err(invalid)?;
            Ok(PyFiniteTypeDecision {
                value: decision.value(),
                reason: format!("{:?}", decision.reason()),
                messages: decision
                    .messages()
                    .iter()
                    .map(ToString::to_string)
                    .collect(),
            })
        })
    }

    #[getter]
    fn is_finite_type_algebra(&self) -> PyResult<Option<bool>> {
        Ok(self.finite_type_decision()?.value)
    }
}

#[pymethods]
impl PyMatrixModule {
    #[staticmethod]
    fn all_modules_of_length_plus_one(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
    ) -> PyResult<Vec<PyMatrixModule>> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        match modules.first() {
            Some(MatrixModuleKind::Rational(_)) | None => {
                modules_plus_one_for::<RationalField>(modules)
            }
            Some(MatrixModuleKind::Binary(_)) => modules_plus_one_for::<Z2>(modules),
            Some(MatrixModuleKind::Prime(_)) => modules_plus_one_for::<Zp>(modules),
            Some(MatrixModuleKind::Extension(_)) => modules_plus_one_for::<GaloisField>(modules),
        }
    }

    fn as_based_right_algebra_module(&self) -> PyResult<PyBasedRightAlgebraModule> {
        PyBasedRightAlgebraModule::from_matrix_module(self)
    }

    fn decompose_via_top(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_modules(
            qpa_core::decompose_module_via_top(module).map_err(invalid)?,
        )))
    }

    fn decompose_via_char_poly(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_modules(
            qpa_core::decompose_module_via_char_poly(module).map_err(invalid)?,
        )))
    }

    fn decompose_probabilistic(
        &self,
        py: Python<'_>,
        endomorphism_basis: Vec<Py<PyModuleHomomorphism>>,
    ) -> PyResult<Py<PyAny>> {
        let homomorphisms = clone_module_homomorphisms(py, endomorphism_basis);
        match &self.kind {
            MatrixModuleKind::Rational(module) => {
                decomposition_probabilistic_for(py, module, homomorphisms)
            }
            MatrixModuleKind::Binary(module) => {
                decomposition_probabilistic_for(py, module, homomorphisms)
            }
            MatrixModuleKind::Prime(module) => {
                decomposition_probabilistic_for(py, module, homomorphisms)
            }
            MatrixModuleKind::Extension(module) => {
                decomposition_probabilistic_for(py, module, homomorphisms)
            }
        }
    }

    fn tensor_product(&self, right: &Self) -> PyResult<PyModuleTensorProduct> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            MatrixModuleKind,
            |left, right| tensor_product_modules_for(left, right),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }
}

#[pymethods]
impl PyModuleElement {
    fn element_in_indec_projective(
        &self,
        algebra: &PyQuiverAlgebra,
        projective_vertex: usize,
    ) -> PyResult<PyAlgebraElement> {
        let kind = match (&algebra.kind, &self.kind) {
            (QuiverAlgebraKind::Rational(algebra), ModuleElementKind::Rational(element)) => {
                AlgebraElementKind::Rational(
                    qpa_gbnp::element_in_indec_projective(algebra, element, projective_vertex)
                        .map_err(invalid)?,
                )
            }
            (QuiverAlgebraKind::Binary(algebra), ModuleElementKind::Binary(element)) => {
                AlgebraElementKind::Binary(
                    qpa_gbnp::element_in_indec_projective(algebra, element, projective_vertex)
                        .map_err(invalid)?,
                )
            }
            (QuiverAlgebraKind::Prime(algebra), ModuleElementKind::Prime(element)) => {
                AlgebraElementKind::Prime(
                    qpa_gbnp::element_in_indec_projective(algebra, element, projective_vertex)
                        .map_err(invalid)?,
                )
            }
            (QuiverAlgebraKind::Extension(algebra), ModuleElementKind::Extension(element)) => {
                AlgebraElementKind::Extension(
                    qpa_gbnp::element_in_indec_projective(algebra, element, projective_vertex)
                        .map_err(invalid)?,
                )
            }
            _ => {
                return Err(invalid(
                    "algebra and element use different coefficient fields",
                ));
            }
        };
        Ok(PyAlgebraElement { kind })
    }
}
