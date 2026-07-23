//! Python wrapper for QPA's nice-basis handling of quotient-element spaces.

use pyo3::pymethods;

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

#[derive(Clone)]
pub(super) enum QuotientNiceBasisKind {
    Rational(QuotientNiceBasis<RationalField>),
    Binary(QuotientNiceBasis<Z2>),
    Prime(QuotientNiceBasis<Zp>),
    Extension(QuotientNiceBasis<GaloisField>),
}

macro_rules! impl_nice_basis_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<QuotientNiceBasis<$field>> for QuotientNiceBasisKind {
            fn wrap(value: QuotientNiceBasis<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_nice_basis_runtime!(RationalField, Rational);
impl_nice_basis_runtime!(Z2, Binary);
impl_nice_basis_runtime!(Zp, Prime);
impl_nice_basis_runtime!(GaloisField, Extension);

/// QPA's nice-basis conversion data for a space of quotient elements.
#[pyclass(name = "QuotientNiceBasis", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuotientNiceBasis {
    kind: QuotientNiceBasisKind,
}

pub(super) fn nice_basis_for<F: runtime::PythonField>(
    quotient: &QuotientPathAlgebra<F>,
    generators: Vec<QuotientElementKind>,
) -> PyResult<PyQuotientNiceBasis>
where
    QuotientElementKind: runtime::Extract<QuotientPathAlgebraElement<F>>,
    QuotientNiceBasisKind: runtime::Wrap<QuotientNiceBasis<F>>,
{
    Ok(PyQuotientNiceBasis {
        kind: runtime::Wrap::wrap(
            QuotientNiceBasis::new(quotient.clone(), &runtime::extract_all(generators)?)
                .map_err(invalid)?,
        ),
    })
}

#[pymethods]
impl PyQuotientNiceBasis {
    #[getter]
    fn quotient(&self) -> PyQuotientAlgebra {
        let kind = map_ref!(
            &self.kind,
            QuotientNiceBasisKind => QuotientKind,
            |basis| basis.quotient().clone()
        );
        PyQuotientAlgebra { kind }
    }

    #[getter]
    fn monomials(&self) -> Vec<PyPath> {
        dispatch_ref!(&self.kind, QuotientNiceBasisKind, |basis| basis
            .monomials()
            .iter()
            .cloned()
            .map(|inner| PyPath { inner })
            .collect())
    }

    fn nice_vector(&self, element: &PyQuotientElement) -> PyResult<Option<Vec<PyScalar>>> {
        let element = element.kind.clone();
        dispatch_ref!(&self.kind, QuotientNiceBasisKind, |basis| {
            let element = runtime::Extract::extract(element)?;
            Ok(basis
                .nice_vector(&element)
                .map_err(invalid)?
                .map(|coordinates| {
                    coordinates
                        .iter()
                        .map(|coefficient| {
                            basis
                                .quotient()
                                .original_path_algebra()
                                .field()
                                .wrap_scalar(coefficient)
                        })
                        .collect()
                }))
        })
    }

    fn element_from_nice_vector(
        &self,
        py: Python<'_>,
        coordinates: Vec<Py<PyScalar>>,
    ) -> PyResult<Option<PyQuotientElement>> {
        dispatch_ref!(&self.kind, QuotientNiceBasisKind, |basis| {
            let field = basis.quotient().original_path_algebra().field();
            let coordinates = coordinates
                .into_iter()
                .map(|coefficient| field.scalar_value(&coefficient.borrow(py)))
                .collect::<PyResult<Vec<_>>>()?;
            Ok(basis
                .element_from_nice_vector(&coordinates)
                .map_err(invalid)?
                .map(runtime::Wrap::wrap)
                .map(|kind| PyQuotientElement { kind }))
        })
    }

    /// QPA names the inverse nice-basis conversion `UglyVector`.
    fn ugly_vector(
        &self,
        py: Python<'_>,
        coordinates: Vec<Py<PyScalar>>,
    ) -> PyResult<Option<PyQuotientElement>> {
        self.element_from_nice_vector(py, coordinates)
    }
}
