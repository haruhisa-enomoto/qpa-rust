//! Python wrappers for QPA centers and finite unit groups.

use qpa_core::{
    AlgebraDimension, PathAlgebraCentre, QuiverAlgebra, QuiverAlgebraElement,
    QuiverAlgebraUnitGroup, QuotientPathAlgebraCentre,
};

use super::conversion::QuiverAlgebraKind;
#[allow(clippy::wildcard_imports)]
// Runtime wrappers and coefficient variants belong to the parent algebra domain.
use super::*;

#[derive(Clone)]
pub(super) enum PathAlgebraCentreKind {
    Rational(PathAlgebraCentre<RationalField>),
    Binary(PathAlgebraCentre<Z2>),
    Prime(PathAlgebraCentre<Zp>),
}

/// QPA's center subalgebra for an unfactored path algebra.
#[pyclass(name = "PathAlgebraCentre", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPathAlgebraCentre {
    kind: PathAlgebraCentreKind,
}

#[derive(Clone)]
pub(super) enum QuotientPathAlgebraCentreKind {
    Rational(QuotientPathAlgebraCentre<RationalField>),
    Binary(QuotientPathAlgebraCentre<Z2>),
    Prime(QuotientPathAlgebraCentre<Zp>),
}

/// QPA's known-basis center subalgebra for a quotient path algebra.
#[pyclass(name = "QuotientPathAlgebraCentre", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuotientPathAlgebraCentre {
    kind: QuotientPathAlgebraCentreKind,
}

#[derive(Clone)]
pub(super) enum QuiverAlgebraUnitGroupKind {
    Rational(QuiverAlgebraUnitGroup<RationalField>),
    Binary(QuiverAlgebraUnitGroup<Z2>),
    Prime(QuiverAlgebraUnitGroup<Zp>),
}

/// QPA's finite table group returned by `Units(A)`.
#[pyclass(name = "QuiverAlgebraUnitGroup", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuiverAlgebraUnitGroup {
    kind: QuiverAlgebraUnitGroupKind,
}

macro_rules! impl_finite_structure_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<PathAlgebraCentre<$field>> for PathAlgebraCentreKind {
            fn wrap(value: PathAlgebraCentre<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<QuotientPathAlgebraCentre<$field>> for QuotientPathAlgebraCentreKind {
            fn wrap(value: QuotientPathAlgebraCentre<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<QuiverAlgebraUnitGroup<$field>> for QuiverAlgebraUnitGroupKind {
            fn wrap(value: QuiverAlgebraUnitGroup<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_finite_structure_runtime!(RationalField, Rational);
impl_finite_structure_runtime!(Z2, Binary);
impl_finite_structure_runtime!(Zp, Prime);

pub(super) fn path_centre_for<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
) -> PyResult<PyPathAlgebraCentre>
where
    PathAlgebraCentreKind: runtime::Wrap<PathAlgebraCentre<F>>,
{
    Ok(PyPathAlgebraCentre {
        kind: runtime::Wrap::wrap(algebra.centre().map_err(invalid)?),
    })
}

pub(super) fn quotient_centre_for<F: runtime::PythonField>(
    algebra: &QuotientPathAlgebra<F>,
) -> PyResult<PyQuotientPathAlgebraCentre>
where
    QuotientPathAlgebraCentreKind: runtime::Wrap<QuotientPathAlgebraCentre<F>>,
{
    Ok(PyQuotientPathAlgebraCentre {
        kind: runtime::Wrap::wrap(qpa_gbnp::quotient_centre(algebra).map_err(invalid)?),
    })
}

pub(super) fn units_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<PyQuiverAlgebraUnitGroup>
where
    QuiverAlgebraUnitGroupKind: runtime::Wrap<QuiverAlgebraUnitGroup<F>>,
{
    Ok(PyQuiverAlgebraUnitGroup {
        kind: runtime::Wrap::wrap(qpa_gbnp::units(algebra).map_err(invalid)?),
    })
}

pub(super) fn elements_of_units_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Vec<PyQuiverAlgebraElement>>
where
    tensor::QuiverAlgebraElementKind: runtime::Wrap<QuiverAlgebraElement<F>>,
{
    Ok(algebra
        .elements_of_units()
        .map_err(invalid)?
        .into_iter()
        .map(tensor::wrap_quiver_element)
        .collect())
}

fn path_centre_contains_for<F: runtime::PythonField>(
    centre: &PathAlgebraCentre<F>,
    element: AlgebraElementKind,
) -> PyResult<bool>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let element = runtime::Extract::extract(element)?;
    centre.contains(&element).map_err(invalid)
}

fn quotient_centre_contains_for<F: runtime::PythonField>(
    centre: &QuotientPathAlgebraCentre<F>,
    element: QuotientElementKind,
) -> PyResult<bool>
where
    QuotientElementKind: runtime::Extract<QuotientPathAlgebraElement<F>>,
{
    let element = runtime::Extract::extract(element)?;
    centre.contains(&element).map_err(invalid)
}

#[pymethods]
impl PyPathAlgebraCentre {
    #[getter]
    fn ambient(&self) -> PyPathAlgebra {
        let kind = map_ref!(
            &self.kind,
            PathAlgebraCentreKind => AlgebraKind,
            |centre| centre.ambient().clone()
        );
        PyPathAlgebra { kind }
    }

    #[getter]
    fn generators(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, PathAlgebraCentreKind, |centre| {
            wrap_algebra_elements(centre.generators().to_vec())
        })
    }

    #[getter]
    fn dimension(&self) -> Option<String> {
        dispatch_ref!(&self.kind, PathAlgebraCentreKind, |centre| {
            match centre.dimension() {
                AlgebraDimension::Finite(value) => Some(value.to_string()),
                AlgebraDimension::Infinite => None,
            }
        })
    }

    fn contains(&self, element: &PyAlgebraElement) -> PyResult<bool> {
        dispatch_ref!(&self.kind, PathAlgebraCentreKind, |centre| {
            path_centre_contains_for(centre, element.kind.clone())
        })
    }

    fn __contains__(&self, element: &PyAlgebraElement) -> PyResult<bool> {
        self.contains(element)
    }
}

#[pymethods]
impl PyQuotientPathAlgebraCentre {
    #[getter]
    fn ambient(&self) -> PyQuotientAlgebra {
        let kind = map_ref!(
            &self.kind,
            QuotientPathAlgebraCentreKind => QuotientKind,
            |centre| centre.ambient().clone()
        );
        PyQuotientAlgebra { kind }
    }

    #[getter]
    fn basis(&self) -> Vec<PyQuotientElement> {
        dispatch_ref!(&self.kind, QuotientPathAlgebraCentreKind, |centre| {
            wrap_quotient_elements(centre.basis().to_vec())
        })
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, QuotientPathAlgebraCentreKind, |centre| centre
            .dimension())
    }

    fn contains(&self, element: &PyQuotientElement) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientPathAlgebraCentreKind, |centre| {
            quotient_centre_contains_for(centre, element.kind.clone())
        })
    }

    fn __contains__(&self, element: &PyQuotientElement) -> PyResult<bool> {
        self.contains(element)
    }
}

#[pymethods]
impl PyQuiverAlgebraUnitGroup {
    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        let kind = map_ref!(
            &self.kind,
            QuiverAlgebraUnitGroupKind => QuiverAlgebraKind,
            |group| group.algebra().clone()
        );
        PyQuiverAlgebra { kind }
    }

    #[getter]
    fn elements(&self) -> Vec<PyQuiverAlgebraElement> {
        dispatch_ref!(&self.kind, QuiverAlgebraUnitGroupKind, |group| {
            tensor::wrap_quiver_elements(group.elements())
        })
    }

    #[getter]
    fn coordinates(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, QuiverAlgebraUnitGroupKind, |group| group
            .coordinates()
            .iter()
            .map(|row| {
                row.iter()
                    .map(|coefficient| group.algebra().field().wrap_scalar(coefficient))
                    .collect()
            })
            .collect())
    }

    #[getter]
    fn order(&self) -> usize {
        dispatch_ref!(&self.kind, QuiverAlgebraUnitGroupKind, |group| group
            .order())
    }

    #[getter]
    fn multiplication_table(&self) -> Vec<Vec<usize>> {
        dispatch_ref!(&self.kind, QuiverAlgebraUnitGroupKind, |group| group
            .multiplication_table()
            .iter()
            .map(|row| row.to_vec())
            .collect())
    }

    #[getter]
    fn identity_index(&self) -> usize {
        dispatch_ref!(&self.kind, QuiverAlgebraUnitGroupKind, |group| group
            .identity_index())
    }

    fn product_index(&self, left: usize, right: usize) -> Option<usize> {
        dispatch_ref!(&self.kind, QuiverAlgebraUnitGroupKind, |group| group
            .product_index(left, right))
    }

    fn inverse_index(&self, index: usize) -> Option<usize> {
        dispatch_ref!(&self.kind, QuiverAlgebraUnitGroupKind, |group| group
            .inverse_index(index))
    }

    fn __len__(&self) -> usize {
        self.order()
    }
}
