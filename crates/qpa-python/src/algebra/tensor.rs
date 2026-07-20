//! Python wrappers for QPA tensor products and enveloping algebras.

use pyo3::pymethods;
use qpa_core::{QuiverAlgebraElement, TensorProductAlgebra};

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

#[derive(Clone)]
pub(super) enum QuiverAlgebraElementKind {
    Rational(QuiverAlgebraElement<RationalField>),
    Binary(QuiverAlgebraElement<Z2>),
    Prime(QuiverAlgebraElement<Zp>),
}

/// QPA's tensor-product quotient with its retained factors and product quiver.
#[pyclass(name = "TensorProductAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyTensorProductAlgebra {
    kind: TensorProductAlgebraKind,
}

/// An element of QPA's path-algebra-or-quotient union.
#[pyclass(name = "QuiverAlgebraElement", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuiverAlgebraElement {
    kind: QuiverAlgebraElementKind,
}

macro_rules! impl_tensor_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<TensorProductAlgebra<$field>> for TensorProductAlgebraKind {
            fn wrap(value: TensorProductAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<QuiverAlgebraElement<$field>> for QuiverAlgebraElementKind {
            fn wrap(value: QuiverAlgebraElement<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Extract<QuiverAlgebraElement<$field>> for QuiverAlgebraElementKind {
            fn extract(self) -> PyResult<QuiverAlgebraElement<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "quiver-algebra element belongs to another coefficient field",
                    )),
                }
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

pub(super) fn include_path_for<F: runtime::PythonField>(
    path: &PyPath,
    algebra: &QuiverAlgebra<F>,
) -> PyResult<PyQuiverAlgebraElement>
where
    QuiverAlgebraElementKind: runtime::Wrap<QuiverAlgebraElement<F>>,
{
    Ok(PyQuiverAlgebraElement {
        kind: runtime::Wrap::wrap(
            qpa_core::include_in_path_algebra(&path.inner, algebra).map_err(invalid)?,
        ),
    })
}

pub(super) fn vertices_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Vec<PyQuiverAlgebraElement>>
where
    QuiverAlgebraElementKind: runtime::Wrap<QuiverAlgebraElement<F>>,
{
    Ok(qpa_core::vertices_of_path_algebra(algebra)
        .map_err(invalid)?
        .into_iter()
        .map(|element| PyQuiverAlgebraElement {
            kind: runtime::Wrap::wrap(element),
        })
        .collect())
}

pub(super) fn regular_enveloping_module_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    Ok(qpa_gbnp::algebra_as_module_over_enveloping_algebra(algebra)
        .map_err(invalid)?
        .map(module::wrap_module))
}

pub(super) fn dual_enveloping_module_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    Ok(
        qpa_gbnp::dual_of_algebra_as_module_over_enveloping_algebra(algebra)
            .map_err(invalid)?
            .map(module::wrap_module),
    )
}

fn simple_tensor_for<F: runtime::PythonField>(
    tensor: &TensorProductAlgebra<F>,
    left: QuiverAlgebraElementKind,
    right: QuiverAlgebraElementKind,
) -> PyResult<PyQuiverAlgebraElement>
where
    QuiverAlgebraElementKind:
        runtime::Extract<QuiverAlgebraElement<F>> + runtime::Wrap<QuiverAlgebraElement<F>>,
{
    let left = runtime::Extract::extract(left)?;
    let right = runtime::Extract::extract(right)?;
    Ok(PyQuiverAlgebraElement {
        kind: runtime::Wrap::wrap(tensor.simple_tensor([&left, &right]).map_err(invalid)?),
    })
}

fn element_algebra<F: runtime::PythonField>(element: &QuiverAlgebraElement<F>) -> QuiverAlgebra<F> {
    match element {
        QuiverAlgebraElement::Path(element) => QuiverAlgebra::Path(element.algebra().clone()),
        QuiverAlgebraElement::Quotient(element) => {
            QuiverAlgebra::Quotient(element.algebra().clone())
        }
    }
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

    /// Apply QPA's bilinear `SimpleTensor` expansion.
    fn simple_tensor(
        &self,
        left: &PyQuiverAlgebraElement,
        right: &PyQuiverAlgebraElement,
    ) -> PyResult<PyQuiverAlgebraElement> {
        dispatch_ref!(&self.kind, TensorProductAlgebraKind, |tensor| {
            simple_tensor_for(tensor, left.kind.clone(), right.kind.clone())
        })
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

#[pymethods]
impl PyQuiverAlgebraElement {
    /// Lift a concrete path-algebra element into QPA's algebra union.
    #[staticmethod]
    fn from_path_element(element: &PyAlgebraElement) -> Self {
        Self {
            kind: map_ref!(
                &element.kind,
                AlgebraElementKind => QuiverAlgebraElementKind,
                |element| QuiverAlgebraElement::Path(element.clone())
            ),
        }
    }

    /// Lift a concrete quotient element into QPA's algebra union.
    #[staticmethod]
    fn from_quotient_element(element: &PyQuotientElement) -> Self {
        Self {
            kind: map_ref!(
                &element.kind,
                QuotientElementKind => QuiverAlgebraElementKind,
                |element| QuiverAlgebraElement::Quotient(element.clone())
            ),
        }
    }

    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        let kind = map_ref!(
            &self.kind,
            QuiverAlgebraElementKind => QuiverAlgebraKind,
            |element| element_algebra(element)
        );
        PyQuiverAlgebra { kind }
    }

    #[getter]
    fn path_element(&self) -> Option<PyAlgebraElement> {
        dispatch_ref!(&self.kind, QuiverAlgebraElementKind, |element| {
            match element {
                QuiverAlgebraElement::Path(element) => Some(PyAlgebraElement {
                    kind: runtime::Wrap::wrap(element.clone()),
                }),
                QuiverAlgebraElement::Quotient(_) => None,
            }
        })
    }

    #[getter]
    fn quotient_element(&self) -> Option<PyQuotientElement> {
        dispatch_ref!(&self.kind, QuiverAlgebraElementKind, |element| {
            match element {
                QuiverAlgebraElement::Path(_) => None,
                QuiverAlgebraElement::Quotient(element) => Some(PyQuotientElement {
                    kind: runtime::Wrap::wrap(element.clone()),
                }),
            }
        })
    }

    #[getter]
    fn terms(&self) -> PyResult<Vec<(PyPath, PyScalar)>> {
        dispatch_ref!(&self.kind, QuiverAlgebraElementKind, |element| {
            Ok(qpa_core::path_algebra_element_terms(element)
                .map_err(invalid)?
                .into_iter()
                .map(|term| {
                    (
                        PyPath {
                            inner: term.monomial().clone(),
                        },
                        element.field().wrap_scalar(term.coefficient()),
                    )
                })
                .collect())
        })
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, QuiverAlgebraElementKind, |element| {
            element.is_zero()
        })
    }

    #[getter]
    fn leading_term(&self) -> Self {
        Self {
            kind: map_ref!(
                &self.kind,
                QuiverAlgebraElementKind => QuiverAlgebraElementKind,
                |element| element.leading_term()
            ),
        }
    }

    #[getter]
    fn leading_coefficient(&self) -> PyScalar {
        dispatch_ref!(&self.kind, QuiverAlgebraElementKind, |element| {
            element.field().wrap_scalar(&element.leading_coefficient())
        })
    }

    #[getter]
    fn leading_monomial(&self) -> PyPath {
        let inner = dispatch_ref!(&self.kind, QuiverAlgebraElementKind, |element| {
            element.leading_monomial()
        });
        PyPath { inner }
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            QuiverAlgebraElementKind => QuiverAlgebraElementKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn subtract(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            QuiverAlgebraElementKind => QuiverAlgebraElementKind,
            |left, right| left.sub(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn multiply(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            QuiverAlgebraElementKind => QuiverAlgebraElementKind,
            |left, right| left.multiply(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn scale_left(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuiverAlgebraElementKind => QuiverAlgebraElementKind,
            |element| element
                .scale_left(&scalar.value_for(element.field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn __add__(&self, right: &Self) -> PyResult<Self> {
        self.add(right)
    }

    fn __sub__(&self, right: &Self) -> PyResult<Self> {
        self.subtract(right)
    }

    fn __mul__(&self, right: &Self) -> PyResult<Self> {
        self.multiply(right)
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, QuiverAlgebraElementKind, |element| {
            match element {
                QuiverAlgebraElement::Path(element) => {
                    element.print_obj_with(|value| element.algebra().field().scalar_text(value))
                }
                QuiverAlgebraElement::Quotient(element) => element.print_obj_with(|value| {
                    element
                        .algebra()
                        .original_path_algebra()
                        .field()
                        .scalar_text(value)
                }),
            }
        })
    }

    fn __eq__(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            QuiverAlgebraElementKind,
            |left, right| left == right,
            false
        )
    }
}
