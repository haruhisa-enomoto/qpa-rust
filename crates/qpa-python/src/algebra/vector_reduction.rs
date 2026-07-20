//! Python exposure of QPA's path-vector reduction helpers.

use pyo3::types::PyList;
use qpa_core::{VectorCoefficients, VectorTip};

#[allow(clippy::wildcard_imports)]
// Runtime wrappers and coefficient variants belong to the parent algebra domain.
use super::*;

/// Result of QPA's `CoefficientsOfVectors` path-element overload.
#[pyclass(name = "PathVectorCoefficientsResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPathVectorCoefficientsResult {
    coefficients: Option<Vec<PyScalar>>,
    non_span_basis: Option<Vec<Vec<PyAlgebraElement>>>,
    non_span_vector: Option<Vec<PyAlgebraElement>>,
}

/// Result of QPA's `CoefficientsOfVectors` quotient-element overload.
#[pyclass(name = "QuotientVectorCoefficientsResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuotientVectorCoefficientsResult {
    coefficients: Option<Vec<PyScalar>>,
    non_span_basis: Option<Vec<Vec<PyQuotientElement>>>,
    non_span_vector: Option<Vec<PyQuotientElement>>,
}

fn wrap_tip<F: runtime::PythonField>(
    field: &F,
    tip: VectorTip<F::Element>,
) -> (usize, PyPath, PyScalar) {
    (
        tip.position(),
        PyPath {
            inner: tip.monomial().clone(),
        },
        field.wrap_scalar(tip.coefficient()),
    )
}

pub(super) fn path_tip_for<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    vector: Vec<AlgebraElementKind>,
) -> PyResult<(usize, PyPath, PyScalar)>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let vector = runtime::extract_all(vector)?;
    Ok(wrap_tip(
        algebra.field(),
        algebra
            .tip_monomial_and_coefficient_of_vector(&vector)
            .map_err(invalid)?,
    ))
}

pub(super) fn quotient_tip_for<F: runtime::PythonField>(
    algebra: &QuotientPathAlgebra<F>,
    vector: Vec<QuotientElementKind>,
) -> PyResult<(usize, PyPath, PyScalar)>
where
    QuotientElementKind: runtime::Extract<QuotientPathAlgebraElement<F>>,
{
    let vector = runtime::extract_all(vector)?;
    Ok(wrap_tip(
        algebra.original_path_algebra().field(),
        algebra
            .tip_monomial_and_coefficient_of_vector(&vector)
            .map_err(invalid)?,
    ))
}

pub(super) fn path_rows(vectors: &Bound<'_, PyList>) -> PyResult<Vec<Vec<AlgebraElementKind>>> {
    vectors
        .iter()
        .map(|row| {
            row.cast::<PyList>()?
                .iter()
                .map(|element| {
                    element
                        .extract::<PyRef<'_, PyAlgebraElement>>()
                        .map(|element| element.kind.clone())
                        .map_err(Into::into)
                })
                .collect()
        })
        .collect()
}

pub(super) fn quotient_rows(
    vectors: &Bound<'_, PyList>,
) -> PyResult<Vec<Vec<QuotientElementKind>>> {
    vectors
        .iter()
        .map(|row| {
            row.cast::<PyList>()?
                .iter()
                .map(|element| {
                    element
                        .extract::<PyRef<'_, PyQuotientElement>>()
                        .map(|element| element.kind.clone())
                        .map_err(Into::into)
                })
                .collect()
        })
        .collect()
}

pub(super) fn path_reduce_for<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    vectors: Vec<Vec<AlgebraElementKind>>,
) -> PyResult<Option<Vec<Vec<PyAlgebraElement>>>>
where
    AlgebraElementKind:
        runtime::Extract<PathAlgebraElement<F>> + runtime::Wrap<PathAlgebraElement<F>>,
{
    let mut vectors = vectors
        .into_iter()
        .map(runtime::extract_all)
        .collect::<PyResult<Vec<_>>>()?;
    if !algebra.tip_reduce_vectors(&mut vectors).map_err(invalid)? {
        return Ok(None);
    }
    Ok(Some(
        vectors.into_iter().map(wrap_algebra_elements).collect(),
    ))
}

pub(super) fn quotient_reduce_for<F: runtime::PythonField>(
    algebra: &QuotientPathAlgebra<F>,
    vectors: Vec<Vec<QuotientElementKind>>,
) -> PyResult<Option<Vec<Vec<PyQuotientElement>>>>
where
    QuotientElementKind: runtime::Extract<QuotientPathAlgebraElement<F>>
        + runtime::Wrap<QuotientPathAlgebraElement<F>>,
{
    let mut vectors = vectors
        .into_iter()
        .map(runtime::extract_all)
        .collect::<PyResult<Vec<_>>>()?;
    if !algebra.tip_reduce_vectors(&mut vectors).map_err(invalid)? {
        return Ok(None);
    }
    Ok(Some(
        vectors.into_iter().map(wrap_quotient_elements).collect(),
    ))
}

pub(super) fn replace_path_rows(
    vectors: &Bound<'_, PyList>,
    rows: Vec<Vec<PyAlgebraElement>>,
) -> PyResult<()> {
    for (row_index, row) in rows.into_iter().enumerate() {
        let target = vectors.get_item(row_index)?.cast_into::<PyList>()?;
        for (column, element) in row.into_iter().enumerate() {
            target.set_item(column, Py::new(vectors.py(), element)?)?;
        }
    }
    Ok(())
}

pub(super) fn replace_quotient_rows(
    vectors: &Bound<'_, PyList>,
    rows: Vec<Vec<PyQuotientElement>>,
) -> PyResult<()> {
    for (row_index, row) in rows.into_iter().enumerate() {
        let target = vectors.get_item(row_index)?.cast_into::<PyList>()?;
        for (column, element) in row.into_iter().enumerate() {
            target.set_item(column, Py::new(vectors.py(), element)?)?;
        }
    }
    Ok(())
}

pub(super) fn path_coefficients_for<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    basis: Vec<Vec<AlgebraElementKind>>,
    vector: Vec<AlgebraElementKind>,
) -> PyResult<PyPathVectorCoefficientsResult>
where
    AlgebraElementKind:
        runtime::Extract<PathAlgebraElement<F>> + runtime::Wrap<PathAlgebraElement<F>>,
{
    let basis = basis
        .into_iter()
        .map(runtime::extract_all)
        .collect::<PyResult<Vec<_>>>()?;
    let vector = runtime::extract_all(vector)?;
    Ok(
        match algebra
            .coefficients_of_vectors(&basis, &vector)
            .map_err(invalid)?
        {
            VectorCoefficients::Coefficients(coefficients) => PyPathVectorCoefficientsResult {
                coefficients: Some(
                    coefficients
                        .iter()
                        .map(|coefficient| algebra.field().wrap_scalar(coefficient))
                        .collect(),
                ),
                non_span_basis: None,
                non_span_vector: None,
            },
            VectorCoefficients::NonSpanWitness { basis, vector } => {
                PyPathVectorCoefficientsResult {
                    coefficients: None,
                    non_span_basis: Some(basis.into_iter().map(wrap_algebra_elements).collect()),
                    non_span_vector: Some(wrap_algebra_elements(vector)),
                }
            }
        },
    )
}

pub(super) fn quotient_coefficients_for<F: runtime::PythonField>(
    algebra: &QuotientPathAlgebra<F>,
    basis: Vec<Vec<QuotientElementKind>>,
    vector: Vec<QuotientElementKind>,
) -> PyResult<PyQuotientVectorCoefficientsResult>
where
    QuotientElementKind: runtime::Extract<QuotientPathAlgebraElement<F>>
        + runtime::Wrap<QuotientPathAlgebraElement<F>>,
{
    let basis = basis
        .into_iter()
        .map(runtime::extract_all)
        .collect::<PyResult<Vec<_>>>()?;
    let vector = runtime::extract_all(vector)?;
    let field = algebra.original_path_algebra().field();
    Ok(
        match algebra
            .coefficients_of_vectors(&basis, &vector)
            .map_err(invalid)?
        {
            VectorCoefficients::Coefficients(coefficients) => PyQuotientVectorCoefficientsResult {
                coefficients: Some(
                    coefficients
                        .iter()
                        .map(|coefficient| field.wrap_scalar(coefficient))
                        .collect(),
                ),
                non_span_basis: None,
                non_span_vector: None,
            },
            VectorCoefficients::NonSpanWitness { basis, vector } => {
                PyQuotientVectorCoefficientsResult {
                    coefficients: None,
                    non_span_basis: Some(basis.into_iter().map(wrap_quotient_elements).collect()),
                    non_span_vector: Some(wrap_quotient_elements(vector)),
                }
            }
        },
    )
}

#[pymethods]
impl PyPathVectorCoefficientsResult {
    #[getter]
    fn coefficients(&self) -> Option<Vec<PyScalar>> {
        self.coefficients.clone()
    }

    #[getter]
    fn is_non_span_witness(&self) -> bool {
        self.coefficients.is_none()
    }

    #[getter]
    fn non_span_basis(&self) -> Option<Vec<Vec<PyAlgebraElement>>> {
        self.non_span_basis.clone()
    }

    #[getter]
    fn non_span_vector(&self) -> Option<Vec<PyAlgebraElement>> {
        self.non_span_vector.clone()
    }
}

#[pymethods]
impl PyQuotientVectorCoefficientsResult {
    #[getter]
    fn coefficients(&self) -> Option<Vec<PyScalar>> {
        self.coefficients.clone()
    }

    #[getter]
    fn is_non_span_witness(&self) -> bool {
        self.coefficients.is_none()
    }

    #[getter]
    fn non_span_basis(&self) -> Option<Vec<Vec<PyQuotientElement>>> {
        self.non_span_basis.clone()
    }

    #[getter]
    fn non_span_vector(&self) -> Option<Vec<PyQuotientElement>> {
        self.non_span_vector.clone()
    }
}
