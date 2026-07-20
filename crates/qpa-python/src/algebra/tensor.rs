//! Python wrappers for QPA tensor products and enveloping algebras.

use pyo3::pymethods;
use qpa_core::TensorProductAlgebra;

use super::conversion::QuiverAlgebraKind;
#[allow(clippy::wildcard_imports)]
// This binding domain intentionally shares the parent algebra runtime enums.
use super::*;
use crate::quiver::PyQuiverProduct;

#[derive(Clone)]
pub(super) enum TensorProductAlgebraKind {
    Rational(TensorProductAlgebra<RationalField>),
    Binary(TensorProductAlgebra<Z2>),
    Prime(TensorProductAlgebra<Zp>),
}

/// QPA's tensor-product quotient with its retained factors and product quiver.
#[pyclass(name = "TensorProductAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyTensorProductAlgebra {
    kind: TensorProductAlgebraKind,
}

macro_rules! impl_tensor_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<TensorProductAlgebra<$field>> for TensorProductAlgebraKind {
            fn wrap(value: TensorProductAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_tensor_runtime!(RationalField, Rational);
impl_tensor_runtime!(Z2, Binary);
impl_tensor_runtime!(Zp, Prime);

pub(super) fn tensor_product_for<F: runtime::PythonField>(
    left: &QuiverAlgebra<F>,
    right: &QuiverAlgebra<F>,
) -> PyResult<PyTensorProductAlgebra>
where
    TensorProductAlgebraKind: runtime::Wrap<TensorProductAlgebra<F>>,
{
    Ok(PyTensorProductAlgebra {
        kind: runtime::Wrap::wrap(
            qpa_gbnp::tensor_product_of_algebras(left, right).map_err(invalid)?,
        ),
    })
}

pub(super) fn enveloping_algebra_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<PyTensorProductAlgebra>
where
    TensorProductAlgebraKind: runtime::Wrap<TensorProductAlgebra<F>>,
{
    Ok(PyTensorProductAlgebra {
        kind: runtime::Wrap::wrap(qpa_gbnp::enveloping_algebra(algebra).map_err(invalid)?),
    })
}

#[pymethods]
impl PyTensorProductAlgebra {
    /// Return the tensor product as QPA's path-or-quotient algebra union.
    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        PyQuiverAlgebra {
            kind: map_ref!(
                &self.kind,
                TensorProductAlgebraKind => QuiverAlgebraKind,
                |tensor| tensor.algebra().clone()
            ),
        }
    }

    /// Return the quotient algebra constructed by QPA's tensor method.
    #[getter]
    fn quotient(&self) -> PyQuotientAlgebra {
        PyQuotientAlgebra {
            kind: map_ref!(
                &self.kind,
                TensorProductAlgebraKind => QuotientKind,
                |tensor| tensor.quotient().clone()
            ),
        }
    }

    /// Return QPA's retained product quiver and decomposition metadata.
    #[getter]
    fn product_quiver(&self) -> PyQuiverProduct {
        let product = dispatch_ref!(&self.kind, TensorProductAlgebraKind, |tensor| {
            tensor.product_quiver().clone()
        });
        PyQuiverProduct::from_inner(product)
    }

    /// Return QPA's retained `TensorProductDecomposition` list.
    #[getter]
    fn factors(&self) -> Vec<PyQuiverAlgebra> {
        dispatch_ref!(&self.kind, TensorProductAlgebraKind, |tensor| {
            tensor
                .factors()
                .iter()
                .cloned()
                .map(conversion::wrap_quiver_algebra)
                .collect()
        })
    }

    /// Alias the retained factor list by its QPA attribute role.
    #[getter]
    fn decomposition(&self) -> Vec<PyQuiverAlgebra> {
        self.factors()
    }

    /// Select one retained factor using QPA's one-based indexing.
    fn factor(&self, one_based_index: usize) -> PyResult<PyQuiverAlgebra> {
        dispatch_ref!(&self.kind, TensorProductAlgebraKind, |tensor| {
            Ok(conversion::wrap_quiver_algebra(
                tensor.factor(one_based_index).map_err(invalid)?.clone(),
            ))
        })
    }

    /// Apply QPA's `TensorAlgebraInclusion(T, n)` convention.
    fn inclusion(&self, factor_number: i64) -> PyResult<PyAlgebraHomomorphism> {
        let kind = map_ref!(
            &self.kind,
            TensorProductAlgebraKind => AlgebraHomomorphismKind,
            |tensor| qpa_gbnp::tensor_algebra_inclusion(tensor, factor_number)
                .map_err(invalid)?
        );
        Ok(PyAlgebraHomomorphism { kind })
    }

    /// Return whether two handles are the identical retained tensor object.
    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            TensorProductAlgebraKind,
            |left, right| left.is_identical(right),
            false
        )
    }
}
