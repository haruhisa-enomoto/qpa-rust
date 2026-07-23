//! Python wrappers for QPA's quotient-algebra vertex projectives.

use numerica::domains::{
    Field,
    finite_field::{Z2, Zp},
    rational::RationalField,
};
use pyo3::{Py, PyResult, Python, prelude::*};
use qpa_core::{
    LiftedProjectivePathAlgebraModule, QuotientPathAlgebra, QuotientPathAlgebraVector,
    RightProjectiveQuotientPathAlgebraModule, right_projective_quotient_path_algebra_module,
};

#[allow(clippy::wildcard_imports)]
use super::*;

#[derive(Clone)]
pub(super) enum QuotientProjectiveKind {
    Rational(RightProjectiveQuotientPathAlgebraModule<RationalField>),
    Binary(RightProjectiveQuotientPathAlgebraModule<Z2>),
    Prime(RightProjectiveQuotientPathAlgebraModule<Zp>),
    Extension(RightProjectiveQuotientPathAlgebraModule<GaloisField>),
}

/// QPA's quotient-path-algebra vertex projective.
#[pyclass(name = "RightProjectiveQuotientModule", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyRightProjectiveQuotientModule {
    kind: QuotientProjectiveKind,
}

#[derive(Clone)]
enum QuotientProjectiveVectorKind {
    Rational(QuotientPathAlgebraVector<RationalField>),
    Binary(QuotientPathAlgebraVector<Z2>),
    Prime(QuotientPathAlgebraVector<Zp>),
    Extension(QuotientPathAlgebraVector<GaloisField>),
}

/// QPA's positional vector for a quotient vertex projective.
#[pyclass(name = "QuotientPathAlgebraVector", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuotientPathAlgebraVector {
    kind: QuotientProjectiveVectorKind,
}

#[derive(Clone)]
enum LiftedProjectiveKind {
    Rational(LiftedProjectivePathAlgebraModule<RationalField>),
    Binary(LiftedProjectivePathAlgebraModule<Z2>),
    Prime(LiftedProjectivePathAlgebraModule<Zp>),
    Extension(LiftedProjectivePathAlgebraModule<GaloisField>),
}

/// QPA's path-algebra lift of a quotient vertex projective.
#[pyclass(
    name = "LiftedProjectivePathAlgebraModule",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PyLiftedProjectivePathAlgebraModule {
    kind: LiftedProjectiveKind,
}

macro_rules! impl_quotient_projective_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<RightProjectiveQuotientPathAlgebraModule<$field>>
            for QuotientProjectiveKind
        {
            fn wrap(value: RightProjectiveQuotientPathAlgebraModule<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<QuotientPathAlgebraVector<$field>> for QuotientProjectiveVectorKind {
            fn wrap(value: QuotientPathAlgebraVector<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Extract<QuotientPathAlgebraVector<$field>> for QuotientProjectiveVectorKind {
            fn extract(self) -> PyResult<QuotientPathAlgebraVector<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "quotient projective vector belongs to another coefficient field",
                    )),
                }
            }
        }

        impl runtime::Wrap<LiftedProjectivePathAlgebraModule<$field>> for LiftedProjectiveKind {
            fn wrap(value: LiftedProjectivePathAlgebraModule<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_quotient_projective_runtime!(RationalField, Rational);
impl_quotient_projective_runtime!(Z2, Binary);
impl_quotient_projective_runtime!(Zp, Prime);
impl_quotient_projective_runtime!(GaloisField, Extension);

fn quotient_projective<F: runtime::PythonField>(
    algebra: &QuotientPathAlgebra<F>,
    vertices: &[usize],
) -> PyResult<RightProjectiveQuotientPathAlgebraModule<F>> {
    let number_of_vertices = algebra
        .original_path_algebra()
        .quiver()
        .number_of_vertices();
    let generators = algebra.generators().map_err(invalid)?;
    let vertices = vertices
        .iter()
        .map(|&position| {
            if position == 0 || position > number_of_vertices {
                return Err(invalid(format!(
                    "vertex index {position} is outside 1..={number_of_vertices}"
                )));
            }
            Ok(generators[position - 1].clone())
        })
        .collect::<PyResult<Vec<_>>>()?;
    right_projective_quotient_path_algebra_module(algebra, &vertices).map_err(invalid)
}

fn clone_quotient_element_kinds(
    py: Python<'_>,
    elements: Vec<Py<PyQuotientElement>>,
) -> Vec<QuotientElementKind> {
    elements
        .into_iter()
        .map(|element| element.borrow(py).kind.clone())
        .collect()
}

fn wrap_vectors<F: Field>(
    values: Vec<QuotientPathAlgebraVector<F>>,
) -> Vec<PyQuotientPathAlgebraVector>
where
    QuotientProjectiveVectorKind: runtime::Wrap<QuotientPathAlgebraVector<F>>,
{
    runtime::wrap_all(values)
        .into_iter()
        .map(|kind| PyQuotientPathAlgebraVector { kind })
        .collect()
}

pub(super) fn right_projective_module(
    algebra: &PyQuotientAlgebra,
    vertices: Vec<usize>,
) -> PyResult<PyRightProjectiveQuotientModule> {
    let kind = map_ref!(
        &algebra.kind,
        QuotientKind => QuotientProjectiveKind,
        |algebra| quotient_projective(algebra, &vertices)?
    );
    Ok(PyRightProjectiveQuotientModule { kind })
}

pub(super) fn wrap_quotient_projective<F>(
    value: RightProjectiveQuotientPathAlgebraModule<F>,
) -> PyRightProjectiveQuotientModule
where
    F: Field,
    QuotientProjectiveKind: runtime::Wrap<RightProjectiveQuotientPathAlgebraModule<F>>,
{
    PyRightProjectiveQuotientModule {
        kind: runtime::Wrap::wrap(value),
    }
}

#[pymethods]
impl PyRightProjectiveQuotientModule {
    #[staticmethod]
    fn from_vertices(algebra: &PyQuotientAlgebra, vertices: Vec<usize>) -> PyResult<Self> {
        right_projective_module(algebra, vertices)
    }

    #[getter]
    fn algebra(&self) -> PyQuotientAlgebra {
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveKind => QuotientKind,
            |module| module.algebra().clone()
        );
        PyQuotientAlgebra { kind }
    }

    #[getter]
    fn vertex_indices(&self) -> Vec<usize> {
        dispatch_ref!(&self.kind, QuotientProjectiveKind, |module| module
            .vertex_indices()
            .iter()
            .map(|vertex| vertex.get())
            .collect())
    }

    fn zero(&self) -> PyQuotientPathAlgebraVector {
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveKind => QuotientProjectiveVectorKind,
            |module| module.zero()
        );
        PyQuotientPathAlgebraVector { kind }
    }

    #[getter]
    fn generators(&self) -> PyResult<Vec<PyQuotientPathAlgebraVector>> {
        dispatch_ref!(&self.kind, QuotientProjectiveKind, |module| Ok(
            wrap_vectors(module.generators().map_err(invalid)?)
        ))
    }

    fn vector(
        &self,
        py: Python<'_>,
        components: Vec<Py<PyQuotientElement>>,
    ) -> PyResult<PyQuotientPathAlgebraVector> {
        let components = clone_quotient_element_kinds(py, components);
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveKind => QuotientProjectiveVectorKind,
            |module| module
                .vector(runtime::extract_all(components)?)
                .map_err(invalid)?
        );
        Ok(PyQuotientPathAlgebraVector { kind })
    }

    fn vector_nc(
        &self,
        py: Python<'_>,
        components: Vec<Py<PyQuotientElement>>,
        leading_position: usize,
        normal: bool,
    ) -> PyResult<PyQuotientPathAlgebraVector> {
        let components = clone_quotient_element_kinds(py, components);
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveKind => QuotientProjectiveVectorKind,
            |module| module.vector_nc(
                runtime::extract_all(components)?,
                leading_position,
                normal,
            )
        );
        Ok(PyQuotientPathAlgebraVector { kind })
    }

    fn canonical_basis(&self) -> PyResult<Option<Vec<PyQuotientPathAlgebraVector>>> {
        dispatch_ref!(&self.kind, QuotientProjectiveKind, |module| Ok(module
            .canonical_basis()
            .map_err(invalid)?
            .map(wrap_vectors)))
    }

    #[getter]
    fn dimension(&self) -> PyResult<Option<usize>> {
        dispatch_ref!(&self.kind, QuotientProjectiveKind, |module| module
            .dimension()
            .map_err(invalid))
    }

    #[getter]
    fn is_vertex_projective(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientProjectiveKind, |module| module
            .is_vertex_projective_module())
    }

    fn lift(&self) -> PyResult<PyLiftedProjectivePathAlgebraModule> {
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveKind => LiftedProjectiveKind,
            |module| module.lift_path_algebra_module().map_err(invalid)?
        );
        Ok(PyLiftedProjectivePathAlgebraModule { kind })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            QuotientProjectiveKind,
            |left, right| left.is_identical(right),
            false
        )
    }
}

#[pymethods]
impl PyQuotientPathAlgebraVector {
    #[getter]
    fn module(&self) -> PyRightProjectiveQuotientModule {
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveVectorKind => QuotientProjectiveKind,
            |vector| vector.module().clone()
        );
        PyRightProjectiveQuotientModule { kind }
    }

    #[getter]
    fn components(&self) -> Vec<PyQuotientElement> {
        dispatch_ref!(&self.kind, QuotientProjectiveVectorKind, |vector| {
            vector
                .components()
                .iter()
                .cloned()
                .map(|value| PyQuotientElement {
                    kind: runtime::Wrap::wrap(value),
                })
                .collect()
        })
    }

    #[getter]
    fn leading_position(&self) -> usize {
        dispatch_ref!(&self.kind, QuotientProjectiveVectorKind, |vector| vector
            .leading_position())
    }

    #[getter]
    fn is_normal_form(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientProjectiveVectorKind, |vector| vector
            .is_normal_form())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientProjectiveVectorKind, |vector| vector
            .is_zero())
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            QuotientProjectiveVectorKind => QuotientProjectiveVectorKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("quotient projective vectors use different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn subtract(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            QuotientProjectiveVectorKind => QuotientProjectiveVectorKind,
            |left, right| left.sub(right).map_err(invalid)?,
            invalid("quotient projective vectors use different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn negate(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveVectorKind => QuotientProjectiveVectorKind,
            |vector| vector.neg().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_left(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveVectorKind => QuotientProjectiveVectorKind,
            |vector| vector
                .scale_left(
                    &scalar.value_for(vector.module().algebra().original_path_algebra().field())?,
                )
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_right(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveVectorKind => QuotientProjectiveVectorKind,
            |vector| vector
                .scale_right(
                    &scalar.value_for(vector.module().algebra().original_path_algebra().field())?,
                )
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn act(&self, element: &PyQuotientElement) -> PyResult<Self> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveVectorKind => QuotientProjectiveVectorKind,
            |vector| vector
                .act(&runtime::Extract::extract(element)?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[getter]
    fn leading_term(&self) -> Self {
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveVectorKind => QuotientProjectiveVectorKind,
            |vector| vector.leading_term()
        );
        Self { kind }
    }

    #[getter]
    fn leading_coefficient(&self) -> PyResult<PyScalar> {
        dispatch_ref!(&self.kind, QuotientProjectiveVectorKind, |vector| {
            let field = vector.module().algebra().original_path_algebra().field();
            Ok(field.wrap_scalar(&vector.leading_coefficient().map_err(invalid)?))
        })
    }

    #[getter]
    fn tip(&self) -> Self {
        self.leading_term()
    }

    #[getter]
    fn tip_coefficient(&self) -> PyResult<PyScalar> {
        self.leading_coefficient()
    }

    #[getter]
    fn leading_component(&self) -> PyResult<PyQuotientElement> {
        let kind = map_ref!(
            &self.kind,
            QuotientProjectiveVectorKind => QuotientElementKind,
            |vector| vector.leading_component().map_err(invalid)?.clone()
        );
        Ok(PyQuotientElement { kind })
    }

    #[getter]
    fn target_vertex(&self) -> Option<PyPath> {
        dispatch_ref!(&self.kind, QuotientProjectiveVectorKind, |vector| vector
            .target_vertex()
            .map(|inner| PyPath { inner }))
    }

    #[getter]
    fn source_vertex(&self) -> Option<PyPath> {
        dispatch_ref!(&self.kind, QuotientProjectiveVectorKind, |vector| vector
            .source_vertex()
            .map(|inner| PyPath { inner }))
    }

    #[getter]
    fn is_left_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientProjectiveVectorKind, |vector| vector
            .is_left_uniform())
    }

    #[getter]
    fn is_right_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientProjectiveVectorKind, |vector| vector
            .is_right_uniform())
    }

    #[getter]
    fn is_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientProjectiveVectorKind, |vector| vector
            .is_uniform())
    }

    fn compare(&self, right: &Self) -> PyResult<i8> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            QuotientProjectiveVectorKind,
            |left, right| Ok(ordering_value(left.compare(right).map_err(invalid)?)),
            Err(invalid(
                "quotient projective vectors use different coefficient fields"
            ))
        )
    }

    fn is_left_divisible(&self, divisor: &Self) -> PyResult<bool> {
        dispatch_pair_ref!(
            &self.kind,
            &divisor.kind,
            QuotientProjectiveVectorKind,
            |vector, divisor| vector.is_left_divisible(divisor).map_err(invalid),
            Err(invalid(
                "quotient projective vectors use different coefficient fields"
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
            QuotientProjectiveVectorKind,
            |left, right| left == right,
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, QuotientProjectiveVectorKind, |vector| vector
            .print_obj_with(|value| vector
                .module()
                .algebra()
                .original_path_algebra()
                .field()
                .scalar_text(value)))
    }

    fn print_obj(&self) -> String {
        self.__repr__()
    }
}

#[pymethods]
impl PyLiftedProjectivePathAlgebraModule {
    #[getter]
    fn parent_module(&self) -> PyRightProjectiveModule {
        dispatch_ref!(&self.kind, LiftedProjectiveKind, |lifted| {
            projective::wrap_projective_module(lifted.parent_module().clone())
        })
    }

    #[getter]
    fn submodule(&self) -> Option<PyRightProjectiveSubmodule> {
        dispatch_ref!(&self.kind, LiftedProjectiveKind, |lifted| lifted
            .submodule()
            .cloned()
            .map(projective::wrap_projective_submodule))
    }

    #[getter]
    fn is_submodule(&self) -> bool {
        dispatch_ref!(&self.kind, LiftedProjectiveKind, |lifted| lifted
            .is_submodule())
    }
}
