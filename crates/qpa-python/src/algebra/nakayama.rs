//! Python closures for QPA's Nakayama and Frobenius attributes.

use pyo3::types::PyString;
use qpa_core::{DenseMatrix, PathAlgebraModule, QuiverAlgebra, QuiverAlgebraElement};
use qpa_gbnp::{
    FrobeniusForm, FrobeniusLinearFunctional, NakayamaAutomorphism, NakayamaAutomorphismOrder,
    NakayamaPermutation,
};

#[allow(clippy::wildcard_imports)]
// Runtime field variants and wrappers belong to the parent algebra domain.
use super::*;
use super::{
    conversion::wrap_quiver_algebra,
    tensor::{QuiverAlgebraElementKind, wrap_quiver_element},
};
use crate::field::{parse_python_integer, wrap_python_integer};

#[derive(Clone)]
pub(super) enum NakayamaPermutationKind {
    Rational(NakayamaPermutation<RationalField>),
    Binary(NakayamaPermutation<Z2>),
    Prime(NakayamaPermutation<Zp>),
}

#[derive(Clone)]
pub(super) enum NakayamaAutomorphismKind {
    Rational(NakayamaAutomorphism<RationalField>),
    Binary(NakayamaAutomorphism<Z2>),
    Prime(NakayamaAutomorphism<Zp>),
}

#[derive(Clone)]
pub(super) enum FrobeniusFormKind {
    Rational(FrobeniusForm<RationalField>),
    Binary(FrobeniusForm<Z2>),
    Prime(FrobeniusForm<Zp>),
}

#[derive(Clone)]
pub(super) enum FrobeniusLinearFunctionalKind {
    Rational(FrobeniusLinearFunctional<RationalField>),
    Binary(FrobeniusLinearFunctional<Z2>),
    Prime(FrobeniusLinearFunctional<Zp>),
}

/// QPA's two captured Nakayama-permutation closures.
#[pyclass(name = "NakayamaPermutation", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyNakayamaPermutation {
    kind: NakayamaPermutationKind,
}

/// QPA's captured Nakayama-automorphism closure.
#[pyclass(name = "NakayamaAutomorphism", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyNakayamaAutomorphism {
    kind: NakayamaAutomorphismKind,
}

/// QPA's captured Frobenius bilinear-form closure.
#[pyclass(name = "FrobeniusForm", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyFrobeniusForm {
    kind: FrobeniusFormKind,
}

/// QPA's captured Frobenius linear-functional closure.
#[pyclass(name = "FrobeniusLinearFunctional", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyFrobeniusLinearFunctional {
    kind: FrobeniusLinearFunctionalKind,
}

macro_rules! impl_nakayama_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<NakayamaPermutation<$field>> for NakayamaPermutationKind {
            fn wrap(value: NakayamaPermutation<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<NakayamaAutomorphism<$field>> for NakayamaAutomorphismKind {
            fn wrap(value: NakayamaAutomorphism<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<FrobeniusForm<$field>> for FrobeniusFormKind {
            fn wrap(value: FrobeniusForm<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<FrobeniusLinearFunctional<$field>> for FrobeniusLinearFunctionalKind {
            fn wrap(value: FrobeniusLinearFunctional<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_nakayama_runtime!(RationalField, Rational);
impl_nakayama_runtime!(Z2, Binary);
impl_nakayama_runtime!(Zp, Prime);

fn matrix_rows<F: runtime::PythonField>(matrix: &DenseMatrix<F>) -> Vec<Vec<PyScalar>> {
    (0..matrix.nrows())
        .map(|row| {
            matrix
                .row(row)
                .iter()
                .map(|value| matrix.field().wrap_scalar(value))
                .collect()
        })
        .collect()
}

pub(super) fn permutation_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Option<PyNakayamaPermutation>>
where
    NakayamaPermutationKind: runtime::Wrap<NakayamaPermutation<F>>,
{
    Ok(qpa_gbnp::nakayama_permutation(algebra)
        .map_err(invalid)?
        .map(|permutation| PyNakayamaPermutation {
            kind: runtime::Wrap::wrap(permutation),
        }))
}

pub(super) fn automorphism_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Option<PyNakayamaAutomorphism>>
where
    NakayamaAutomorphismKind: runtime::Wrap<NakayamaAutomorphism<F>>,
{
    Ok(qpa_gbnp::nakayama_automorphism(algebra)
        .map_err(invalid)?
        .map(|automorphism| PyNakayamaAutomorphism {
            kind: runtime::Wrap::wrap(automorphism),
        }))
}

pub(super) fn order_for<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Option<Py<PyAny>>> {
    qpa_gbnp::order_of_nakayama_automorphism(algebra)
        .map_err(invalid)?
        .map(|order| match order {
            NakayamaAutomorphismOrder::Finite(value) => wrap_python_integer(py, &value),
            NakayamaAutomorphismOrder::Infinite => {
                Ok(PyString::new(py, "infinity").into_any().unbind())
            }
        })
        .transpose()
}

pub(super) fn form_for<F: runtime::PythonField>(
    quotient: &QuotientPathAlgebra<F>,
) -> PyResult<Option<PyFrobeniusForm>>
where
    FrobeniusFormKind: runtime::Wrap<FrobeniusForm<F>>,
{
    Ok(qpa_gbnp::frobenius_form(quotient)
        .map_err(invalid)?
        .map(|form| PyFrobeniusForm {
            kind: runtime::Wrap::wrap(form),
        }))
}

pub(super) fn functional_for<F: runtime::PythonField>(
    quotient: &QuotientPathAlgebra<F>,
) -> PyResult<Option<PyFrobeniusLinearFunctional>>
where
    FrobeniusLinearFunctionalKind: runtime::Wrap<FrobeniusLinearFunctional<F>>,
{
    Ok(qpa_gbnp::frobenius_linear_functional(quotient)
        .map_err(invalid)?
        .map(|functional| PyFrobeniusLinearFunctional {
            kind: runtime::Wrap::wrap(functional),
        }))
}

fn permutation_module_image<F: runtime::PythonField>(
    permutation: &NakayamaPermutation<F>,
    module: MatrixModuleKind,
) -> PyResult<PyMatrixModule>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>> + runtime::Wrap<PathAlgebraModule<F>>,
{
    let module = runtime::Extract::extract(module)?;
    Ok(module::wrap_module(
        permutation.image_of_module(&module).map_err(invalid)?,
    ))
}

fn automorphism_image<F: runtime::PythonField>(
    automorphism: &NakayamaAutomorphism<F>,
    element: QuiverAlgebraElementKind,
) -> PyResult<PyQuiverAlgebraElement>
where
    QuiverAlgebraElementKind:
        runtime::Extract<QuiverAlgebraElement<F>> + runtime::Wrap<QuiverAlgebraElement<F>>,
{
    let element = runtime::Extract::extract(element)?;
    Ok(wrap_quiver_element(
        automorphism.image(&element).map_err(invalid)?,
    ))
}

fn form_evaluate<F: runtime::PythonField>(
    form: &FrobeniusForm<F>,
    left: QuiverAlgebraElementKind,
    right: QuiverAlgebraElementKind,
) -> PyResult<PyScalar>
where
    QuiverAlgebraElementKind: runtime::Extract<QuiverAlgebraElement<F>>,
{
    let left = runtime::Extract::extract(left)?;
    let right = runtime::Extract::extract(right)?;
    Ok(form
        .algebra()
        .original_path_algebra()
        .field()
        .wrap_scalar(&form.evaluate(&left, &right).map_err(invalid)?))
}

fn functional_evaluate<F: runtime::PythonField>(
    functional: &FrobeniusLinearFunctional<F>,
    element: QuiverAlgebraElementKind,
) -> PyResult<PyScalar>
where
    QuiverAlgebraElementKind: runtime::Extract<QuiverAlgebraElement<F>>,
{
    let element = runtime::Extract::extract(element)?;
    Ok(functional
        .algebra()
        .original_path_algebra()
        .field()
        .wrap_scalar(&functional.evaluate(&element).map_err(invalid)?))
}

#[pymethods]
impl PyNakayamaPermutation {
    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, NakayamaPermutationKind, |permutation| {
            wrap_quiver_algebra(permutation.algebra().clone())
        })
    }

    #[getter]
    fn index_images(&self) -> Vec<usize> {
        dispatch_ref!(&self.kind, NakayamaPermutationKind, |permutation| {
            permutation.index_images().to_vec()
        })
    }

    fn image_of_index(&self, index: &Bound<'_, PyAny>) -> PyResult<usize> {
        let index = parse_python_integer(index)?;
        dispatch_ref!(&self.kind, NakayamaPermutationKind, |permutation| {
            permutation.image_of_index(&index).map_err(invalid)
        })
    }

    fn image_of_module(&self, module: &PyMatrixModule) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, NakayamaPermutationKind, |permutation| {
            permutation_module_image(permutation, module.kind.clone())
        })
    }

    fn __repr__(&self) -> String {
        format!("NakayamaPermutation({:?})", self.index_images())
    }
}

#[pymethods]
impl PyNakayamaAutomorphism {
    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, NakayamaAutomorphismKind, |automorphism| {
            wrap_quiver_algebra(automorphism.algebra().clone())
        })
    }

    #[getter]
    fn matrix(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, NakayamaAutomorphismKind, |automorphism| {
            matrix_rows(automorphism.matrix())
        })
    }

    fn image(&self, element: &PyQuiverAlgebraElement) -> PyResult<PyQuiverAlgebraElement> {
        dispatch_ref!(&self.kind, NakayamaAutomorphismKind, |automorphism| {
            automorphism_image(automorphism, element.kind.clone())
        })
    }

    fn __call__(&self, element: &PyQuiverAlgebraElement) -> PyResult<PyQuiverAlgebraElement> {
        self.image(element)
    }
}

#[pymethods]
impl PyFrobeniusForm {
    #[getter]
    fn algebra(&self) -> PyQuotientAlgebra {
        PyQuotientAlgebra {
            kind: dispatch_ref!(&self.kind, FrobeniusFormKind, |form| {
                runtime::Wrap::wrap(form.algebra().clone())
            }),
        }
    }

    #[getter]
    fn matrix(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, FrobeniusFormKind, |form| {
            matrix_rows(form.matrix())
        })
    }

    fn evaluate(
        &self,
        left: &PyQuiverAlgebraElement,
        right: &PyQuiverAlgebraElement,
    ) -> PyResult<PyScalar> {
        dispatch_ref!(&self.kind, FrobeniusFormKind, |form| {
            form_evaluate(form, left.kind.clone(), right.kind.clone())
        })
    }

    fn __call__(
        &self,
        left: &PyQuiverAlgebraElement,
        right: &PyQuiverAlgebraElement,
    ) -> PyResult<PyScalar> {
        self.evaluate(left, right)
    }
}

#[pymethods]
impl PyFrobeniusLinearFunctional {
    #[getter]
    fn algebra(&self) -> PyQuotientAlgebra {
        PyQuotientAlgebra {
            kind: dispatch_ref!(&self.kind, FrobeniusLinearFunctionalKind, |functional| {
                runtime::Wrap::wrap(functional.algebra().clone())
            }),
        }
    }

    #[getter]
    fn canonical_basis_values(&self) -> Vec<PyScalar> {
        dispatch_ref!(&self.kind, FrobeniusLinearFunctionalKind, |functional| {
            functional
                .canonical_basis_values()
                .iter()
                .map(|value| {
                    functional
                        .algebra()
                        .original_path_algebra()
                        .field()
                        .wrap_scalar(value)
                })
                .collect()
        })
    }

    fn evaluate(&self, element: &PyQuiverAlgebraElement) -> PyResult<PyScalar> {
        dispatch_ref!(&self.kind, FrobeniusLinearFunctionalKind, |functional| {
            functional_evaluate(functional, element.kind.clone())
        })
    }

    fn __call__(&self, element: &PyQuiverAlgebraElement) -> PyResult<PyScalar> {
        self.evaluate(element)
    }
}
