//! QPA's source-faithful uniform-element projective module family.

use pyo3::pymethods;
use qpa_core::{
    PathAlgebraElement, RightUniformPathAlgebraModule, RightUniformPathAlgebraVector,
    right_uniform_path_algebra_module,
};

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

#[derive(Clone)]
enum UniformModuleKind {
    Rational(RightUniformPathAlgebraModule<RationalField>),
    Binary(RightUniformPathAlgebraModule<Z2>),
    Prime(RightUniformPathAlgebraModule<Zp>),
    Extension(RightUniformPathAlgebraModule<GaloisField>),
}

/// QPA's projective module constructed from uniform path-algebra elements.
#[pyclass(name = "RightUniformPathAlgebraModule", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyRightUniformPathAlgebraModule {
    kind: UniformModuleKind,
}

#[derive(Clone)]
enum UniformVectorKind {
    Rational(RightUniformPathAlgebraVector<RationalField>),
    Binary(RightUniformPathAlgebraVector<Z2>),
    Prime(RightUniformPathAlgebraVector<Zp>),
    Extension(RightUniformPathAlgebraVector<GaloisField>),
}

/// A positional vector in QPA's uniform-element projective family.
#[pyclass(name = "RightUniformPathAlgebraVector", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyRightUniformPathAlgebraVector {
    kind: UniformVectorKind,
}

macro_rules! impl_uniform_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<RightUniformPathAlgebraModule<$field>> for UniformModuleKind {
            fn wrap(value: RightUniformPathAlgebraModule<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<RightUniformPathAlgebraVector<$field>> for UniformVectorKind {
            fn wrap(value: RightUniformPathAlgebraVector<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Extract<RightUniformPathAlgebraVector<$field>> for UniformVectorKind {
            fn extract(self) -> PyResult<RightUniformPathAlgebraVector<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "uniform path vector uses a different coefficient field",
                    )),
                }
            }
        }
    };
}

impl_uniform_runtime!(RationalField, Rational);
impl_uniform_runtime!(Z2, Binary);
impl_uniform_runtime!(Zp, Prime);
impl_uniform_runtime!(GaloisField, Extension);

fn clone_element_kinds(
    py: Python<'_>,
    elements: Vec<Py<PyAlgebraElement>>,
) -> Vec<AlgebraElementKind> {
    elements
        .iter()
        .map(|element| element.borrow(py).kind.clone())
        .collect()
}

fn wrap_elements<F: runtime::PythonField>(
    elements: impl IntoIterator<Item = PathAlgebraElement<F>>,
) -> Vec<PyAlgebraElement>
where
    AlgebraElementKind: runtime::Wrap<PathAlgebraElement<F>>,
{
    elements
        .into_iter()
        .map(|element| PyAlgebraElement {
            kind: runtime::Wrap::wrap(element),
        })
        .collect()
}

fn wrap_vectors<F: runtime::PythonField>(
    vectors: impl IntoIterator<Item = RightUniformPathAlgebraVector<F>>,
) -> Vec<PyRightUniformPathAlgebraVector>
where
    UniformVectorKind: runtime::Wrap<RightUniformPathAlgebraVector<F>>,
{
    vectors
        .into_iter()
        .map(|vector| PyRightUniformPathAlgebraVector {
            kind: runtime::Wrap::wrap(vector),
        })
        .collect()
}

#[pymethods]
impl PyPathAlgebra {
    fn right_projective_module_from_uniform_elements(
        &self,
        py: Python<'_>,
        elements: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyRightUniformPathAlgebraModule> {
        let elements = clone_element_kinds(py, elements);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => UniformModuleKind,
            |algebra| right_uniform_path_algebra_module(
                algebra,
                &runtime::extract_all(elements)?
            )
            .map_err(invalid)?
        );
        Ok(PyRightUniformPathAlgebraModule { kind })
    }

    fn is_left_uniform_list_at_path(
        &self,
        py: Python<'_>,
        elements: Vec<Py<PyAlgebraElement>>,
        path: &PyPath,
    ) -> PyResult<bool> {
        let elements = clone_element_kinds(py, elements);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .is_left_uniform_list_at_path(&runtime::extract_all(elements)?, &path.inner)
            .map_err(invalid))
    }

    fn is_right_uniform_list_at_path(
        &self,
        py: Python<'_>,
        elements: Vec<Py<PyAlgebraElement>>,
        path: &PyPath,
    ) -> PyResult<bool> {
        let elements = clone_element_kinds(py, elements);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .is_right_uniform_list_at_path(&runtime::extract_all(elements)?, &path.inner)
            .map_err(invalid))
    }
}

#[pymethods]
impl PyRightUniformPathAlgebraModule {
    #[staticmethod]
    fn from_elements(
        algebra: &PyPathAlgebra,
        py: Python<'_>,
        elements: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<Self> {
        algebra.right_projective_module_from_uniform_elements(py, elements)
    }

    #[getter]
    fn algebra(&self) -> PyPathAlgebra {
        let kind = map_ref!(&self.kind, UniformModuleKind => AlgebraKind, |module| module
            .algebra()
            .clone());
        PyPathAlgebra { kind }
    }

    #[getter]
    fn elements(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, UniformModuleKind, |module| wrap_elements(
            module.elements().iter().cloned()
        ))
    }

    fn zero(&self) -> PyRightUniformPathAlgebraVector {
        let kind = map_ref!(
            &self.kind,
            UniformModuleKind => UniformVectorKind,
            |module| module.zero()
        );
        PyRightUniformPathAlgebraVector { kind }
    }

    #[getter]
    fn generators(&self) -> Vec<PyRightUniformPathAlgebraVector> {
        dispatch_ref!(&self.kind, UniformModuleKind, |module| wrap_vectors(
            module.generators()
        ))
    }

    fn vector(
        &self,
        py: Python<'_>,
        components: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyRightUniformPathAlgebraVector> {
        let components = clone_element_kinds(py, components);
        let kind = map_ref!(
            &self.kind,
            UniformModuleKind => UniformVectorKind,
            |module| module
                .vector(runtime::extract_all(components)?)
                .map_err(invalid)?
        );
        Ok(PyRightUniformPathAlgebraVector { kind })
    }

    fn vector_nc(
        &self,
        py: Python<'_>,
        components: Vec<Py<PyAlgebraElement>>,
        leading_position: usize,
        normal: bool,
    ) -> PyResult<PyRightUniformPathAlgebraVector> {
        let components = clone_element_kinds(py, components);
        let kind = map_ref!(
            &self.kind,
            UniformModuleKind => UniformVectorKind,
            |module| module.vector_nc(
                runtime::extract_all(components)?,
                leading_position,
                normal
            )
        );
        Ok(PyRightUniformPathAlgebraVector { kind })
    }

    #[getter]
    fn is_vertex_projective(&self) -> bool {
        dispatch_ref!(&self.kind, UniformModuleKind, |module| module
            .is_vertex_projective_module())
    }

    fn canonical_basis(&self) -> PyResult<Option<Vec<PyRightUniformPathAlgebraVector>>> {
        dispatch_ref!(&self.kind, UniformModuleKind, |module| Ok(module
            .canonical_basis()
            .map_err(invalid)?
            .map(wrap_vectors)))
    }

    fn uniform_generators(&self) -> PyResult<Vec<PyRightUniformPathAlgebraVector>> {
        dispatch_ref!(&self.kind, UniformModuleKind, |module| Ok(wrap_vectors(
            module.uniform_generators().map_err(invalid)?
        )))
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            UniformModuleKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.is_identical(other)
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, UniformModuleKind, |module| format!(
            "<right-module over <{}[<quiver with {} vertices and {} arrows>]>>",
            module.algebra().field().python_kind().display(),
            module.algebra().quiver().number_of_vertices(),
            module.algebra().quiver().number_of_arrows(),
        ))
    }
}

#[pymethods]
impl PyRightUniformPathAlgebraVector {
    #[getter]
    fn module(&self) -> PyRightUniformPathAlgebraModule {
        let kind = map_ref!(
            &self.kind,
            UniformVectorKind => UniformModuleKind,
            |vector| vector.module().clone()
        );
        PyRightUniformPathAlgebraModule { kind }
    }

    #[getter]
    fn components(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, UniformVectorKind, |vector| wrap_elements(
            vector.components().iter().cloned()
        ))
    }

    #[getter]
    fn leading_position(&self) -> usize {
        dispatch_ref!(&self.kind, UniformVectorKind, |vector| vector
            .leading_position())
    }

    #[getter]
    fn is_normal_form(&self) -> bool {
        dispatch_ref!(&self.kind, UniformVectorKind, |vector| vector
            .is_normal_form())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, UniformVectorKind, |vector| vector.is_zero())
    }

    #[getter]
    fn leading_term(&self) -> Self {
        let kind = map_ref!(
            &self.kind,
            UniformVectorKind => UniformVectorKind,
            |vector| vector.leading_term()
        );
        Self { kind }
    }

    #[getter]
    fn leading_coefficient(&self) -> PyResult<PyScalar> {
        dispatch_ref!(&self.kind, UniformVectorKind, |vector| {
            Ok(vector
                .module()
                .algebra()
                .field()
                .wrap_scalar(&vector.leading_coefficient().map_err(invalid)?))
        })
    }

    #[getter]
    fn leading_component(&self) -> PyResult<PyAlgebraElement> {
        let kind = map_ref!(
            &self.kind,
            UniformVectorKind => AlgebraElementKind,
            |vector| vector.leading_component().map_err(invalid)?.clone()
        );
        Ok(PyAlgebraElement { kind })
    }

    #[getter]
    fn target_vertex(&self) -> Option<PyPath> {
        dispatch_ref!(&self.kind, UniformVectorKind, |vector| vector
            .target_vertex()
            .map(|inner| PyPath { inner }))
    }

    #[getter]
    fn source_vertex(&self) -> Option<PyPath> {
        dispatch_ref!(&self.kind, UniformVectorKind, |vector| vector
            .source_vertex()
            .map(|inner| PyPath { inner }))
    }

    #[getter]
    fn is_left_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, UniformVectorKind, |vector| vector
            .is_left_uniform())
    }

    #[getter]
    fn is_right_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, UniformVectorKind, |vector| vector
            .is_right_uniform())
    }

    #[getter]
    fn is_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, UniformVectorKind, |vector| vector.is_uniform())
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            UniformVectorKind => UniformVectorKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("uniform path vectors use different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn subtract(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            UniformVectorKind => UniformVectorKind,
            |left, right| left.sub(right).map_err(invalid)?,
            invalid("uniform path vectors use different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn negate(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            UniformVectorKind => UniformVectorKind,
            |vector| vector.neg().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_left(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            UniformVectorKind => UniformVectorKind,
            |vector| vector
                .scale_left(&scalar.value_for(vector.module().algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_right(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            UniformVectorKind => UniformVectorKind,
            |vector| vector
                .scale_right(&scalar.value_for(vector.module().algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn act(&self, element: &PyAlgebraElement) -> PyResult<Self> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            UniformVectorKind => UniformVectorKind,
            |vector| vector
                .act(&runtime::Extract::extract(element)?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn compare(&self, right: &Self) -> PyResult<i8> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            UniformVectorKind,
            |left, right| Ok(crate::quiver::ordering_value(
                left.compare(right).map_err(invalid)?
            )),
            Err(invalid(
                "uniform path vectors use different coefficient fields"
            ))
        )
    }

    fn is_left_divisible(&self, divisor: &Self) -> PyResult<bool> {
        dispatch_pair_ref!(
            &self.kind,
            &divisor.kind,
            UniformVectorKind,
            |vector, divisor| vector.is_left_divisible(divisor).map_err(invalid),
            Err(invalid(
                "uniform path vectors use different coefficient fields"
            ))
        )
    }

    fn __add__(&self, right: &Self) -> PyResult<Self> {
        self.add(right)
    }

    fn __sub__(&self, right: &Self) -> PyResult<Self> {
        self.subtract(right)
    }

    fn __neg__(&self) -> PyResult<Self> {
        self.negate()
    }

    fn __lt__(&self, right: &Self) -> PyResult<bool> {
        Ok(self.compare(right)? < 0)
    }

    fn __eq__(&self, right: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            UniformVectorKind,
            |left, right| left == right,
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, UniformVectorKind, |vector| vector
            .print_obj_with(|value| vector
                .module()
                .algebra()
                .field()
                .scalar_text(value)))
    }

    fn print_obj(&self) -> String {
        self.__repr__()
    }
}
