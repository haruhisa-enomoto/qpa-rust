//! Shared construction conveniences for concrete and wrapped quiver algebras.

use super::conversion::QuiverAlgebraKind;
#[allow(clippy::wildcard_imports)]
// These methods share the algebra binding domain's runtime dispatch helpers.
use super::*;

fn build_exact_matrix_module<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: QuiverAlgebra<F>,
    dimension_vector: Vec<usize>,
    actions: Vec<(String, Vec<Vec<Py<PyAny>>>)>,
) -> PyResult<PathAlgebraModule<F>> {
    let actions = actions
        .into_iter()
        .map(|(label, rows)| {
            let rows = rows
                .into_iter()
                .map(|row| {
                    row.into_iter()
                        .map(|value| {
                            crate::field::exact_scalar_for(algebra.field(), value.bind(py))
                        })
                        .collect::<PyResult<Vec<_>>>()
                })
                .collect::<PyResult<Vec<_>>>()?;
            let matrix = DenseMatrix::from_rows(rows, algebra.field().clone()).map_err(invalid)?;
            Ok(NamedArrowMatrix::new(label, matrix))
        })
        .collect::<PyResult<Vec<_>>>()?;
    right_module_over_path_algebra(algebra, dimension_vector, &actions)
        .map_err(invalid)?
        .ok_or_else(|| {
            RelationNotSatisfiedError::new_err("module actions do not satisfy relations")
        })
}

fn extract_algebra(value: &Bound<'_, PyAny>) -> PyResult<PyQuiverAlgebra> {
    if let Ok(algebra) = value.extract::<PyRef<'_, PyQuiverAlgebra>>() {
        return Ok(algebra.clone());
    }
    if let Ok(algebra) = value.extract::<PyRef<'_, PyPathAlgebra>>() {
        return Ok(algebra.as_quiver_algebra());
    }
    if let Ok(algebra) = value.extract::<PyRef<'_, PyQuotientAlgebra>>() {
        return Ok(algebra.as_quiver_algebra());
    }
    Err(pyo3::exceptions::PyTypeError::new_err(
        "tensor factor must be a PathAlgebra, QuotientAlgebra, or QuiverAlgebra",
    ))
}

pub(super) fn tensor_product(
    left: &PyQuiverAlgebra,
    right: &Bound<'_, PyAny>,
) -> PyResult<PyTensorProductAlgebra> {
    let right = extract_algebra(right)?;
    dispatch_pair_ref!(
        &left.kind,
        &right.kind,
        QuiverAlgebraKind,
        |left, right| tensor::tensor_product_for(left, right),
        Err(invalid(
            "tensor factors belong to different coefficient fields"
        ))
    )
}

#[pymethods]
impl PyPathAlgebra {
    /// Construct a checked right module from exact scalar, integer, or fraction entries.
    fn module_from_matrices(
        &self,
        py: Python<'_>,
        dimension_vector: Vec<usize>,
        actions: Vec<(String, Vec<Vec<Py<PyAny>>>)>,
    ) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            Ok(module::wrap_module(build_exact_matrix_module(
                py,
                QuiverAlgebra::Path(algebra.clone()),
                dimension_vector,
                actions,
            )?))
        })
    }

    /// Wrap this algebra while retaining its identity and caches.
    fn as_quiver_algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            wrap_quiver_algebra(QuiverAlgebra::Path(algebra.clone()))
        })
    }

    /// Form a tensor product with a concrete or wrapped quiver algebra.
    fn tensor_product(&self, right: &Bound<'_, PyAny>) -> PyResult<PyTensorProductAlgebra> {
        tensor_product(&self.as_quiver_algebra(), right)
    }

    /// Return the cached enveloping algebra and its tensor-factor data.
    fn enveloping_algebra(&self) -> PyResult<PyTensorProductAlgebra> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            tensor::enveloping_algebra_for(&QuiverAlgebra::Path(algebra.clone()))
        })
    }
}

#[pymethods]
impl PyQuotientAlgebra {
    /// Construct a right module from exact entries and check all quotient relations.
    fn module_from_matrices(
        &self,
        py: Python<'_>,
        dimension_vector: Vec<usize>,
        actions: Vec<(String, Vec<Vec<Py<PyAny>>>)>,
    ) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            Ok(module::wrap_module(build_exact_matrix_module(
                py,
                QuiverAlgebra::Quotient(algebra.clone()),
                dimension_vector,
                actions,
            )?))
        })
    }

    /// Wrap this quotient while retaining its identity, relations, and caches.
    fn as_quiver_algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            wrap_quiver_algebra(QuiverAlgebra::Quotient(algebra.clone()))
        })
    }

    /// Form a tensor product with a concrete or wrapped quiver algebra.
    fn tensor_product(&self, right: &Bound<'_, PyAny>) -> PyResult<PyTensorProductAlgebra> {
        tensor_product(&self.as_quiver_algebra(), right)
    }

    /// Return the cached enveloping algebra and its tensor-factor data.
    fn enveloping_algebra(&self) -> PyResult<PyTensorProductAlgebra> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            tensor::enveloping_algebra_for(&QuiverAlgebra::Quotient(algebra.clone()))
        })
    }
}

#[pymethods]
impl PyQuiverAlgebra {
    /// Construct a checked right module from exact scalar, integer, or fraction entries.
    fn module_from_matrices(
        &self,
        py: Python<'_>,
        dimension_vector: Vec<usize>,
        actions: Vec<(String, Vec<Vec<Py<PyAny>>>)>,
    ) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            Ok(module::wrap_module(build_exact_matrix_module(
                py,
                algebra.clone(),
                dimension_vector,
                actions,
            )?))
        })
    }

    /// Return standard simple modules in vertex order.
    fn simple_modules(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            Ok(wrap_matrix_modules(
                qpa_gbnp::simple_modules(algebra).map_err(invalid)?,
            ))
        })
    }

    /// Return indecomposable projectives in vertex order.
    fn indecomposable_projectives(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            Ok(wrap_matrix_modules(
                qpa_gbnp::indec_projective_modules(algebra).map_err(invalid)?,
            ))
        })
    }

    /// Return indecomposable injectives in vertex order.
    fn indecomposable_injectives(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            Ok(wrap_matrix_modules(
                qpa_gbnp::indec_injective_modules(algebra).map_err(invalid)?,
            ))
        })
    }

    fn zero_module(&self) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            Ok(module::wrap_module(
                qpa_core::zero_module(algebra).map_err(invalid)?,
            ))
        })
    }

    /// Construct a checked right module, including all quotient relations.
    fn module_from_integer_matrices(
        &self,
        dimension_vector: Vec<usize>,
        actions: Vec<(String, Vec<Vec<i64>>)>,
    ) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            Ok(module::wrap_module(build_matrix_module(
                algebra.clone(),
                dimension_vector,
                actions,
            )?))
        })
    }

    /// Return a named vertex or arrow in this algebra.
    fn generator(&self, name: &str) -> PyResult<PyQuiverAlgebraElement> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            let path = algebra
                .quiver()
                .generator_by_name(name)
                .ok_or_else(|| invalid(format!("unknown quiver generator {name:?}")))?
                .path();
            Ok(tensor::wrap_quiver_element(
                qpa_core::include_in_path_algebra(&path, algebra).map_err(invalid)?,
            ))
        })
    }

    /// Return vertices followed by arrows, matching the concrete Python algebra.
    #[getter]
    fn generators(&self) -> PyResult<Vec<PyQuiverAlgebraElement>> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            Ok(tensor::wrap_quiver_elements(
                &algebra.generators_with_one().map_err(invalid)?,
            ))
        })
    }

    fn zero(&self) -> PyQuiverAlgebraElement {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            tensor::wrap_quiver_element(algebra.zero())
        })
    }

    fn one(&self) -> PyResult<PyQuiverAlgebraElement> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            Ok(tensor::wrap_quiver_element(algebra.one().map_err(invalid)?))
        })
    }
}
