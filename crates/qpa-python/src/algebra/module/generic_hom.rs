//! QPA's generic `Hom(R, M, N)` and `End(R, M)` path-ring overloads.

use pyo3::{exceptions::PyTypeError, pymethods, types::PyAny};
use qpa_core::{
    GenericPathAlgebraModuleEndomorphismAlgebra, GenericPathAlgebraModuleHomSpace,
    GenericPathAlgebraModuleHomomorphism, PathAlgebraMatModuleObject,
};

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

#[derive(Clone)]
enum GenericModuleObjectKind {
    Rational(PathAlgebraMatModuleObject<RationalField>),
    Binary(PathAlgebraMatModuleObject<Z2>),
    Prime(PathAlgebraMatModuleObject<Zp>),
}

#[derive(Clone)]
enum GenericModuleHomomorphismKind {
    Rational(GenericPathAlgebraModuleHomomorphism<RationalField>),
    Binary(GenericPathAlgebraModuleHomomorphism<Z2>),
    Prime(GenericPathAlgebraModuleHomomorphism<Zp>),
}

#[derive(Clone)]
enum GenericModuleHomSpaceKind {
    Rational(GenericPathAlgebraModuleHomSpace<RationalField>),
    Binary(GenericPathAlgebraModuleHomSpace<Z2>),
    Prime(GenericPathAlgebraModuleHomSpace<Zp>),
}

#[derive(Clone)]
enum GenericModuleEndomorphismAlgebraKind {
    Rational(GenericPathAlgebraModuleEndomorphismAlgebra<RationalField>),
    Binary(GenericPathAlgebraModuleEndomorphismAlgebra<Z2>),
    Prime(GenericPathAlgebraModuleEndomorphismAlgebra<Zp>),
}

/// One matrix mapping returned by QPA's generic path-ring `Hom` method.
#[pyclass(name = "GenericModuleHomomorphism", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyGenericModuleHomomorphism {
    kind: GenericModuleHomomorphismKind,
}

/// QPA's vector space returned by generic `Hom(R, M, N)`.
#[pyclass(name = "GenericModuleHomSpace", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyGenericModuleHomSpace {
    kind: GenericModuleHomSpaceKind,
}

/// QPA's algebra-with-one returned by generic `End(R, M)`.
#[pyclass(name = "GenericModuleEndomorphismAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyGenericModuleEndomorphismAlgebra {
    kind: GenericModuleEndomorphismAlgebraKind,
}

macro_rules! impl_generic_hom_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<GenericPathAlgebraModuleHomomorphism<$field>>
            for GenericModuleHomomorphismKind
        {
            fn wrap(value: GenericPathAlgebraModuleHomomorphism<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<GenericPathAlgebraModuleHomSpace<$field>> for GenericModuleHomSpaceKind {
            fn wrap(value: GenericPathAlgebraModuleHomSpace<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<GenericPathAlgebraModuleEndomorphismAlgebra<$field>>
            for GenericModuleEndomorphismAlgebraKind
        {
            fn wrap(value: GenericPathAlgebraModuleEndomorphismAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_generic_hom_runtime!(RationalField, Rational);
impl_generic_hom_runtime!(Z2, Binary);
impl_generic_hom_runtime!(Zp, Prime);

fn object_from_python(value: &Bound<'_, PyAny>) -> PyResult<GenericModuleObjectKind> {
    if let Ok(module) = value.extract::<PyRef<'_, PyMatrixModule>>() {
        return Ok(match &module.kind {
            MatrixModuleKind::Rational(module) => {
                GenericModuleObjectKind::Rational(module.clone().into())
            }
            MatrixModuleKind::Binary(module) => {
                GenericModuleObjectKind::Binary(module.clone().into())
            }
            MatrixModuleKind::Prime(module) => {
                GenericModuleObjectKind::Prime(module.clone().into())
            }
        });
    }
    if let Ok(submodule) = value.extract::<PyRef<'_, PyModuleSubmodule>>() {
        return Ok(match &submodule.kind {
            ModuleSubmoduleKind::Rational(submodule) => {
                GenericModuleObjectKind::Rational(submodule.clone().into())
            }
            ModuleSubmoduleKind::Binary(submodule) => {
                GenericModuleObjectKind::Binary(submodule.clone().into())
            }
            ModuleSubmoduleKind::Prime(submodule) => {
                GenericModuleObjectKind::Prime(submodule.clone().into())
            }
        });
    }
    Err(PyTypeError::new_err(
        "expected a MatrixModule or ModuleSubmodule",
    ))
}

fn object_to_python<F: runtime::PythonField>(
    py: Python<'_>,
    object: &PathAlgebraMatModuleObject<F>,
) -> PyResult<Py<PyAny>>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
    ModuleSubmoduleKind: runtime::Wrap<PathAlgebraModuleSubmodule<F>>,
{
    match object {
        PathAlgebraMatModuleObject::Module(module) => {
            Ok(Py::new(py, wrap_module(module.clone()))?.into_any())
        }
        PathAlgebraMatModuleObject::Submodule(submodule) => {
            Ok(Py::new(py, wrap_module_submodule(submodule.clone()))?.into_any())
        }
    }
}

fn wrap_generic_homomorphism<F: runtime::PythonField>(
    homomorphism: GenericPathAlgebraModuleHomomorphism<F>,
) -> PyGenericModuleHomomorphism
where
    GenericModuleHomomorphismKind: runtime::Wrap<GenericPathAlgebraModuleHomomorphism<F>>,
{
    PyGenericModuleHomomorphism {
        kind: runtime::Wrap::wrap(homomorphism),
    }
}

pub(crate) fn generic_hom(
    algebra: &PyPathAlgebra,
    source: &Bound<'_, PyAny>,
    range: &Bound<'_, PyAny>,
) -> PyResult<PyGenericModuleHomSpace> {
    let source = object_from_python(source)?;
    let range = object_from_python(range)?;
    let kind = match (&algebra.kind, source, range) {
        (
            AlgebraKind::Rational(algebra),
            GenericModuleObjectKind::Rational(source),
            GenericModuleObjectKind::Rational(range),
        ) => GenericModuleHomSpaceKind::Rational(
            qpa_core::generic_hom_over_path_ring(
                &QuiverAlgebra::Path(algebra.clone()),
                source,
                range,
            )
            .map_err(invalid)?,
        ),
        (
            AlgebraKind::Binary(algebra),
            GenericModuleObjectKind::Binary(source),
            GenericModuleObjectKind::Binary(range),
        ) => GenericModuleHomSpaceKind::Binary(
            qpa_core::generic_hom_over_path_ring(
                &QuiverAlgebra::Path(algebra.clone()),
                source,
                range,
            )
            .map_err(invalid)?,
        ),
        (
            AlgebraKind::Prime(algebra),
            GenericModuleObjectKind::Prime(source),
            GenericModuleObjectKind::Prime(range),
        ) => GenericModuleHomSpaceKind::Prime(
            qpa_core::generic_hom_over_path_ring(
                &QuiverAlgebra::Path(algebra.clone()),
                source,
                range,
            )
            .map_err(invalid)?,
        ),
        _ => return Err(invalid("Hom inputs use different coefficient fields")),
    };
    Ok(PyGenericModuleHomSpace { kind })
}

pub(crate) fn generic_end(
    algebra: &PyPathAlgebra,
    module: &Bound<'_, PyAny>,
) -> PyResult<PyGenericModuleEndomorphismAlgebra> {
    let module = object_from_python(module)?;
    let kind = match (&algebra.kind, module) {
        (AlgebraKind::Rational(algebra), GenericModuleObjectKind::Rational(module)) => {
            GenericModuleEndomorphismAlgebraKind::Rational(
                qpa_core::generic_end_over_path_ring(&QuiverAlgebra::Path(algebra.clone()), module)
                    .map_err(invalid)?,
            )
        }
        (AlgebraKind::Binary(algebra), GenericModuleObjectKind::Binary(module)) => {
            GenericModuleEndomorphismAlgebraKind::Binary(
                qpa_core::generic_end_over_path_ring(&QuiverAlgebra::Path(algebra.clone()), module)
                    .map_err(invalid)?,
            )
        }
        (AlgebraKind::Prime(algebra), GenericModuleObjectKind::Prime(module)) => {
            GenericModuleEndomorphismAlgebraKind::Prime(
                qpa_core::generic_end_over_path_ring(&QuiverAlgebra::Path(algebra.clone()), module)
                    .map_err(invalid)?,
            )
        }
        _ => return Err(invalid("End inputs use different coefficient fields")),
    };
    Ok(PyGenericModuleEndomorphismAlgebra { kind })
}

#[pymethods]
impl PyGenericModuleHomomorphism {
    #[getter]
    fn source(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        dispatch_ref!(&self.kind, GenericModuleHomomorphismKind, |mapping| {
            object_to_python(py, mapping.source())
        })
    }

    #[getter]
    fn range(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        dispatch_ref!(&self.kind, GenericModuleHomomorphismKind, |mapping| {
            object_to_python(py, mapping.range())
        })
    }

    #[getter]
    fn matrix(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, GenericModuleHomomorphismKind, |mapping| {
            matrix_rows(mapping.matrix())
        })
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_algebra_module_homomorphism(&self) -> bool {
        true
    }

    fn image(&self, element: &PyModuleElement) -> PyResult<PyModuleElement> {
        match (&self.kind, &element.kind) {
            (
                GenericModuleHomomorphismKind::Rational(mapping),
                ModuleElementKind::Rational(element),
            ) => Ok(wrap_module_element(
                mapping.image(element).map_err(invalid)?,
            )),
            (
                GenericModuleHomomorphismKind::Binary(mapping),
                ModuleElementKind::Binary(element),
            ) => Ok(wrap_module_element(
                mapping.image(element).map_err(invalid)?,
            )),
            (GenericModuleHomomorphismKind::Prime(mapping), ModuleElementKind::Prime(element)) => {
                Ok(wrap_module_element(
                    mapping.image(element).map_err(invalid)?,
                ))
            }
            _ => Err(invalid(
                "mapping and element use different coefficient fields",
            )),
        }
    }

    fn __call__(&self, element: &PyModuleElement) -> PyResult<PyModuleElement> {
        self.image(element)
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            GenericModuleHomomorphismKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __eq__(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            GenericModuleHomomorphismKind,
            |left, right| left == right,
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, GenericModuleHomomorphismKind, |mapping| {
            format!(
                "<generic path-algebra module homomorphism with {}x{} matrix>",
                mapping.matrix().nrows(),
                mapping.matrix().ncols(),
            )
        })
    }
}

#[pymethods]
impl PyGenericModuleHomSpace {
    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(&self.kind, GenericModuleHomSpaceKind, |space| {
            space.field().python_kind()
        });
        PyField { kind }
    }

    #[getter]
    fn source(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        dispatch_ref!(&self.kind, GenericModuleHomSpaceKind, |space| {
            object_to_python(py, space.source())
        })
    }

    #[getter]
    fn range(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        dispatch_ref!(&self.kind, GenericModuleHomSpaceKind, |space| {
            object_to_python(py, space.range())
        })
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, GenericModuleHomSpaceKind, |space| {
            space.dimension()
        })
    }

    #[getter]
    fn generators(&self) -> Vec<PyGenericModuleHomomorphism> {
        dispatch_ref!(&self.kind, GenericModuleHomSpaceKind, |space| {
            space
                .generators()
                .iter()
                .cloned()
                .map(wrap_generic_homomorphism)
                .collect()
        })
    }

    fn coefficients(
        &self,
        py: Python<'_>,
        matrix: Vec<Vec<Py<PyScalar>>>,
    ) -> PyResult<Option<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, GenericModuleHomSpaceKind, |space| {
            let matrix = extract_matrix(py, space.field(), matrix)?;
            Ok(space
                .coefficients(&matrix)
                .map_err(invalid)?
                .map(|values| wrap_scalars(space.field(), &values)))
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "<generic path-algebra module Hom space of dimension {} over {}>",
            self.dimension(),
            self.field().kind.display(),
        )
    }
}

#[pymethods]
impl PyGenericModuleEndomorphismAlgebra {
    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(
            &self.kind,
            GenericModuleEndomorphismAlgebraKind,
            |algebra| algebra.field().python_kind()
        );
        PyField { kind }
    }

    #[getter]
    fn module(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        dispatch_ref!(
            &self.kind,
            GenericModuleEndomorphismAlgebraKind,
            |algebra| object_to_python(py, algebra.module())
        )
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(
            &self.kind,
            GenericModuleEndomorphismAlgebraKind,
            |algebra| algebra.dimension()
        )
    }

    #[getter]
    fn generators(&self) -> Vec<PyGenericModuleHomomorphism> {
        dispatch_ref!(
            &self.kind,
            GenericModuleEndomorphismAlgebraKind,
            |algebra| algebra
                .generators()
                .iter()
                .cloned()
                .map(wrap_generic_homomorphism)
                .collect()
        )
    }

    #[getter]
    fn one(&self) -> PyGenericModuleHomomorphism {
        dispatch_ref!(
            &self.kind,
            GenericModuleEndomorphismAlgebraKind,
            |algebra| wrap_generic_homomorphism(algebra.one().clone())
        )
    }

    #[getter]
    fn zero(&self) -> PyGenericModuleHomomorphism {
        dispatch_ref!(
            &self.kind,
            GenericModuleEndomorphismAlgebraKind,
            |algebra| wrap_generic_homomorphism(algebra.zero())
        )
    }

    fn multiply(
        &self,
        left: &PyGenericModuleHomomorphism,
        right: &PyGenericModuleHomomorphism,
    ) -> PyResult<PyGenericModuleHomomorphism> {
        match (&self.kind, &left.kind, &right.kind) {
            (
                GenericModuleEndomorphismAlgebraKind::Rational(algebra),
                GenericModuleHomomorphismKind::Rational(left),
                GenericModuleHomomorphismKind::Rational(right),
            ) => Ok(wrap_generic_homomorphism(
                algebra.multiply(left, right).map_err(invalid)?,
            )),
            (
                GenericModuleEndomorphismAlgebraKind::Binary(algebra),
                GenericModuleHomomorphismKind::Binary(left),
                GenericModuleHomomorphismKind::Binary(right),
            ) => Ok(wrap_generic_homomorphism(
                algebra.multiply(left, right).map_err(invalid)?,
            )),
            (
                GenericModuleEndomorphismAlgebraKind::Prime(algebra),
                GenericModuleHomomorphismKind::Prime(left),
                GenericModuleHomomorphismKind::Prime(right),
            ) => Ok(wrap_generic_homomorphism(
                algebra.multiply(left, right).map_err(invalid)?,
            )),
            _ => Err(invalid(
                "endomorphism algebra and mappings use different coefficient fields",
            )),
        }
    }

    fn coefficients(
        &self,
        py: Python<'_>,
        matrix: Vec<Vec<Py<PyScalar>>>,
    ) -> PyResult<Option<Vec<PyScalar>>> {
        dispatch_ref!(
            &self.kind,
            GenericModuleEndomorphismAlgebraKind,
            |algebra| {
                let matrix = extract_matrix(py, algebra.field(), matrix)?;
                Ok(algebra
                    .coefficients(&matrix)
                    .map_err(invalid)?
                    .map(|values| wrap_scalars(algebra.field(), &values)))
            }
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "<generic path-algebra module endomorphism algebra of dimension {} over {}>",
            self.dimension(),
            self.field().kind.display(),
        )
    }
}
