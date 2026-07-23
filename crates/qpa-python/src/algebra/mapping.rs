//! Python wrappers for QPA path-algebra mappings and namespace assignments.

use qpa_core::{NaturalQuotientMap, QuiverAlgebra, QuiverAlgebraElement};

#[allow(clippy::wildcard_imports)]
// Runtime wrappers and field variants belong to the parent algebra domain.
use super::*;

#[derive(Clone)]
pub(super) enum NaturalQuotientMapKind {
    Rational(NaturalQuotientMap<RationalField>),
    Binary(NaturalQuotientMap<Z2>),
    Prime(NaturalQuotientMap<Zp>),
    Extension(NaturalQuotientMap<GaloisField>),
}

/// QPA's natural surjection from a path algebra to its quotient by an ideal.
#[pyclass(name = "NaturalQuotientMap", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyNaturalQuotientMap {
    kind: NaturalQuotientMapKind,
}

macro_rules! impl_natural_map_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<NaturalQuotientMap<$field>> for NaturalQuotientMapKind {
            fn wrap(value: NaturalQuotientMap<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_natural_map_runtime!(RationalField, Rational);
impl_natural_map_runtime!(Z2, Binary);
impl_natural_map_runtime!(Zp, Prime);
impl_natural_map_runtime!(GaloisField, Extension);

pub(super) fn natural_quotient_map_for<F: runtime::PythonField>(
    ideal: &PathAlgebraIdeal<F>,
) -> PyResult<PyNaturalQuotientMap>
where
    NaturalQuotientMapKind: runtime::Wrap<NaturalQuotientMap<F>>,
{
    Ok(PyNaturalQuotientMap {
        kind: runtime::Wrap::wrap(ideal.natural_homomorphism().map_err(invalid)?),
    })
}

pub(super) fn generator_variable_assignments_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
    read_only_names: &[String],
) -> PyResult<Vec<(String, PyQuiverAlgebraElement)>>
where
    tensor::QuiverAlgebraElementKind: runtime::Wrap<QuiverAlgebraElement<F>>,
{
    algebra
        .generator_variable_assignments_with_read_only(|name| {
            read_only_names.iter().any(|candidate| candidate == name)
        })
        .map_err(invalid)?
        .into_iter()
        .map(|assignment| {
            let (name, element) = assignment.into_parts();
            Ok((name.into(), tensor::wrap_quiver_element(element)))
        })
        .collect()
}

fn image_for<F: runtime::PythonField>(
    mapping: &NaturalQuotientMap<F>,
    element: AlgebraElementKind,
) -> PyResult<PyQuotientElement>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
    QuotientElementKind: runtime::Wrap<QuotientPathAlgebraElement<F>>,
{
    let element = runtime::Extract::extract(element)?;
    Ok(PyQuotientElement {
        kind: runtime::Wrap::wrap(mapping.image(&element).map_err(invalid)?),
    })
}

#[pymethods]
impl PyNaturalQuotientMap {
    #[getter]
    fn source(&self) -> PyPathAlgebra {
        let kind = map_ref!(
            &self.kind,
            NaturalQuotientMapKind => AlgebraKind,
            |mapping| mapping.source().clone()
        );
        PyPathAlgebra { kind }
    }

    #[getter]
    fn range(&self) -> PyQuotientAlgebra {
        let kind = map_ref!(
            &self.kind,
            NaturalQuotientMapKind => QuotientKind,
            |mapping| mapping.range().clone()
        );
        PyQuotientAlgebra { kind }
    }

    #[getter]
    fn kernel(&self) -> PyPathAlgebraIdeal {
        let kind = map_ref!(
            &self.kind,
            NaturalQuotientMapKind => IdealKind,
            |mapping| mapping.kernel().clone()
        );
        PyPathAlgebraIdeal { kind }
    }

    #[getter]
    fn is_surjective(&self) -> bool {
        dispatch_ref!(&self.kind, NaturalQuotientMapKind, |mapping| mapping
            .is_surjective())
    }

    #[getter]
    fn generator_images(&self) -> PyResult<Vec<PyQuotientElement>> {
        dispatch_ref!(&self.kind, NaturalQuotientMapKind, |mapping| Ok(
            wrap_quotient_elements(mapping.generator_images().map_err(invalid)?)
        ))
    }

    fn image(&self, element: &PyAlgebraElement) -> PyResult<PyQuotientElement> {
        dispatch_ref!(&self.kind, NaturalQuotientMapKind, |mapping| {
            image_for(mapping, element.kind.clone())
        })
    }

    fn mapped_expression(&self, element: &PyAlgebraElement) -> PyResult<PyQuotientElement> {
        self.image(element)
    }

    fn images_representative(&self, element: &PyAlgebraElement) -> PyResult<PyQuotientElement> {
        self.image(element)
    }
}
