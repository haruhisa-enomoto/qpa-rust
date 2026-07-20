//! Python wrappers for QPA's finite-dimensional algebra conversions.

use numerica::domains::{
    finite_field::{Z2, Zp},
    rational::{Q, RationalField},
};
use pyo3::prelude::*;
use qpa_core::{PathAlgebraModule, QuiverAlgebra, StructureConstantsAlgebra};

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

#[derive(Clone)]
pub(super) enum StructureConstantsAlgebraKind {
    Rational(StructureConstantsAlgebra<RationalField>),
    Binary(StructureConstantsAlgebra<Z2>),
    Prime(StructureConstantsAlgebra<Zp>),
}

/// A finite-dimensional algebra in QPA's ordered structure-constant basis.
#[pyclass(name = "StructureConstantsAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyStructureConstantsAlgebra {
    pub(super) kind: StructureConstantsAlgebraKind,
}

#[derive(Clone)]
pub(super) enum QuiverAlgebraKind {
    Rational(QuiverAlgebra<RationalField>),
    Binary(QuiverAlgebra<Z2>),
    Prime(QuiverAlgebra<Zp>),
}

/// A QPA path algebra or quotient path algebra returned by a conversion.
#[pyclass(name = "QuiverAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuiverAlgebra {
    pub(super) kind: QuiverAlgebraKind,
}

#[derive(Clone)]
pub(super) enum AlgebraAsQuiverAlgebraResultKind {
    Rational(qpa_gbnp::AlgebraAsQuiverAlgebraResult<RationalField>),
    Binary(qpa_gbnp::AlgebraAsQuiverAlgebraResult<Z2>),
    Prime(qpa_gbnp::AlgebraAsQuiverAlgebraResult<Zp>),
}

/// QPA's converted quiver algebra and ordered generator images.
#[pyclass(name = "AlgebraAsQuiverAlgebraResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAlgebraAsQuiverAlgebraResult {
    pub(super) kind: AlgebraAsQuiverAlgebraResultKind,
}

#[derive(Clone)]
pub(super) enum AdrAlgebraOfAlgebraResultKind {
    Rational(qpa_gbnp::AdrAlgebraOfAlgebraResult<RationalField>),
    Binary(qpa_gbnp::AdrAlgebraOfAlgebraResult<Z2>),
    Prime(qpa_gbnp::AdrAlgebraOfAlgebraResult<Zp>),
}

/// QPA's ADR algebra and installed finite-global-dimension fact.
#[pyclass(name = "AdrAlgebraOfAlgebraResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAdrAlgebraOfAlgebraResult {
    pub(super) kind: AdrAlgebraOfAlgebraResultKind,
}

macro_rules! impl_conversion_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<StructureConstantsAlgebra<$field>> for StructureConstantsAlgebraKind {
            fn wrap(value: StructureConstantsAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<QuiverAlgebra<$field>> for QuiverAlgebraKind {
            fn wrap(value: QuiverAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::AlgebraAsQuiverAlgebraResult<$field>>
            for AlgebraAsQuiverAlgebraResultKind
        {
            fn wrap(value: qpa_gbnp::AlgebraAsQuiverAlgebraResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::AdrAlgebraOfAlgebraResult<$field>>
            for AdrAlgebraOfAlgebraResultKind
        {
            fn wrap(value: qpa_gbnp::AdrAlgebraOfAlgebraResult<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_conversion_runtime!(RationalField, Rational);
impl_conversion_runtime!(Z2, Binary);
impl_conversion_runtime!(Zp, Prime);

fn structure_constants_for<F: runtime::PythonField>(
    py: Python<'_>,
    field: F,
    products: Vec<Vec<Vec<Py<PyScalar>>>>,
) -> PyResult<PyStructureConstantsAlgebra>
where
    StructureConstantsAlgebraKind: runtime::Wrap<StructureConstantsAlgebra<F>>,
{
    let dimension = products.len();
    if products.iter().any(|left| {
        left.len() != dimension || left.iter().any(|product| product.len() != dimension)
    }) {
        return Err(invalid(format!(
            "structure constants must form a {dimension} x {dimension} x {dimension} tensor"
        )));
    }
    let products = products
        .into_iter()
        .flatten()
        .flatten()
        .map(|coefficient| field.scalar_value(&coefficient.borrow(py)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyStructureConstantsAlgebra {
        kind: runtime::Wrap::wrap(
            StructureConstantsAlgebra::new(field, dimension, products).map_err(invalid)?,
        ),
    })
}

pub(super) fn wrap_quiver_algebra<F: runtime::PythonField>(
    algebra: QuiverAlgebra<F>,
) -> PyQuiverAlgebra
where
    QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    PyQuiverAlgebra {
        kind: runtime::Wrap::wrap(algebra),
    }
}

pub(super) fn adr_algebra_of_basic_module_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
) -> PyResult<PyQuiverAlgebra>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>>,
    QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    Ok(wrap_quiver_algebra(
        qpa_gbnp::adr_algebra_of_basic_module(&runtime::extract_all(modules)?).map_err(invalid)?,
    ))
}

#[pymethods]
impl PyStructureConstantsAlgebra {
    #[new]
    fn new(
        py: Python<'_>,
        field: PyRef<'_, PyField>,
        products: Vec<Vec<Vec<Py<PyScalar>>>>,
    ) -> PyResult<Self> {
        match field.kind {
            FieldKind::Rational => structure_constants_for(py, Q, products),
            FieldKind::Binary => structure_constants_for(py, Z2, products),
            FieldKind::Prime(characteristic) => {
                structure_constants_for(py, Zp::new(characteristic), products)
            }
        }
    }

    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            algebra.field().python_kind()
        });
        PyField { kind }
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            algebra.dimension()
        })
    }

    #[getter]
    fn multiplication_tensor(&self) -> Vec<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            (0..algebra.dimension())
                .map(|left| {
                    (0..algebra.dimension())
                        .map(|right| {
                            algebra
                                .basis_product(left, right)
                                .iter()
                                .map(|coefficient| algebra.field().wrap_scalar(coefficient))
                                .collect()
                        })
                        .collect()
                })
                .collect()
        })
    }

    fn as_quiver_algebra(&self) -> PyResult<PyAlgebraAsQuiverAlgebraResult> {
        let kind = map_ref!(
            &self.kind,
            StructureConstantsAlgebraKind => AlgebraAsQuiverAlgebraResultKind,
            |algebra| qpa_gbnp::algebra_as_quiver_algebra(algebra).map_err(invalid)?
        );
        Ok(PyAlgebraAsQuiverAlgebraResult { kind })
    }
}

#[pymethods]
impl PyQuiverAlgebra {
    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            algebra.field().python_kind()
        });
        PyField { kind }
    }

    #[getter]
    fn quiver(&self) -> PyQuiver {
        let inner = dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            algebra.quiver().clone()
        });
        PyQuiver { inner }
    }

    #[getter]
    fn is_quotient(&self) -> bool {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            matches!(algebra, QuiverAlgebra::Quotient(_))
        })
    }

    #[getter]
    fn path_algebra(&self) -> Option<PyPathAlgebra> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| match algebra {
            QuiverAlgebra::Path(algebra) => Some(PyPathAlgebra {
                kind: runtime::Wrap::wrap(algebra.clone()),
            }),
            QuiverAlgebra::Quotient(_) => None,
        })
    }

    #[getter]
    fn quotient_algebra(&self) -> Option<PyQuotientAlgebra> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| match algebra {
            QuiverAlgebra::Path(_) => None,
            QuiverAlgebra::Quotient(algebra) => Some(PyQuotientAlgebra {
                kind: runtime::Wrap::wrap(algebra.clone()),
            }),
        })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            QuiverAlgebraKind,
            |left, right| left.is_identical(right),
            false
        )
    }
}

#[pymethods]
impl PyAlgebraAsQuiverAlgebraResult {
    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, AlgebraAsQuiverAlgebraResultKind, |result| {
            wrap_quiver_algebra(result.algebra().clone())
        })
    }

    #[getter]
    fn generator_images(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, AlgebraAsQuiverAlgebraResultKind, |result| {
            result
                .generator_images()
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|coefficient| result.algebra().field().wrap_scalar(coefficient))
                        .collect()
                })
                .collect()
        })
    }

    #[getter]
    fn loewy_length(&self) -> usize {
        dispatch_ref!(&self.kind, AlgebraAsQuiverAlgebraResultKind, |result| {
            result.loewy_length()
        })
    }
}

#[pymethods]
impl PyAdrAlgebraOfAlgebraResult {
    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, AdrAlgebraOfAlgebraResultKind, |result| {
            wrap_quiver_algebra(result.algebra().clone())
        })
    }

    #[getter]
    fn is_finite_global_dimension_algebra(&self) -> bool {
        dispatch_ref!(&self.kind, AdrAlgebraOfAlgebraResultKind, |result| {
            result.is_finite_global_dimension_algebra()
        })
    }
}
