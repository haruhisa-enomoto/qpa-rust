//! Python wrappers for QPA's vertex-projective path modules and vectors.

use numerica::domains::{
    Field,
    finite_field::{Z2, Zp},
    rational::RationalField,
};
use pyo3::{Py, PyResult, Python, prelude::*};
use qpa_core::{
    PathAlgebra, PathAlgebraVector, RightPathAlgebraModuleGroebnerBasis,
    RightProjectivePathAlgebraModule, completely_reduce_path_algebra_module_groebner_basis,
    right_projective_path_algebra_module,
};

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

#[derive(Clone)]
enum ProjectiveModuleKind {
    Rational(RightProjectivePathAlgebraModule<RationalField>),
    Binary(RightProjectivePathAlgebraModule<Z2>),
    Prime(RightProjectivePathAlgebraModule<Zp>),
}

/// QPA's vertex-projective polynomial right module.
#[pyclass(name = "RightProjectiveModule", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyRightProjectiveModule {
    kind: ProjectiveModuleKind,
}

#[derive(Clone)]
enum ProjectiveVectorKind {
    Rational(PathAlgebraVector<RationalField>),
    Binary(PathAlgebraVector<Z2>),
    Prime(PathAlgebraVector<Zp>),
}

/// QPA's sparse polynomial `PathAlgebraVector`.
#[pyclass(name = "PathAlgebraVector", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPathAlgebraVector {
    kind: ProjectiveVectorKind,
}

#[derive(Clone)]
enum ProjectiveGroebnerKind {
    Rational(RightPathAlgebraModuleGroebnerBasis<RationalField>),
    Binary(RightPathAlgebraModuleGroebnerBasis<Z2>),
    Prime(RightPathAlgebraModuleGroebnerBasis<Zp>),
}

/// QPA's right Gröbner basis for a vertex-projective path module.
#[pyclass(name = "PathAlgebraModuleGroebnerBasis", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPathAlgebraModuleGroebnerBasis {
    kind: ProjectiveGroebnerKind,
}

macro_rules! impl_projective_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<RightProjectivePathAlgebraModule<$field>> for ProjectiveModuleKind {
            fn wrap(value: RightProjectivePathAlgebraModule<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<PathAlgebraVector<$field>> for ProjectiveVectorKind {
            fn wrap(value: PathAlgebraVector<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Extract<PathAlgebraVector<$field>> for ProjectiveVectorKind {
            fn extract(self) -> PyResult<PathAlgebraVector<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "path-algebra vector belongs to another coefficient field",
                    )),
                }
            }
        }

        impl runtime::Wrap<RightPathAlgebraModuleGroebnerBasis<$field>> for ProjectiveGroebnerKind {
            fn wrap(value: RightPathAlgebraModuleGroebnerBasis<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_projective_runtime!(RationalField, Rational);
impl_projective_runtime!(Z2, Binary);
impl_projective_runtime!(Zp, Prime);

fn wrap_vectors<F: Field>(values: Vec<PathAlgebraVector<F>>) -> Vec<PyPathAlgebraVector>
where
    ProjectiveVectorKind: runtime::Wrap<PathAlgebraVector<F>>,
{
    runtime::wrap_all(values)
        .into_iter()
        .map(|kind| PyPathAlgebraVector { kind })
        .collect()
}

fn projective_module<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    vertices: &[usize],
) -> PyResult<RightProjectivePathAlgebraModule<F>> {
    let number_of_vertices = algebra.quiver().number_of_vertices();
    let generators = algebra.generators();
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
    right_projective_path_algebra_module(algebra, &vertices).map_err(invalid)
}

fn clone_element_kinds(
    py: Python<'_>,
    elements: Vec<Py<PyAlgebraElement>>,
) -> Vec<AlgebraElementKind> {
    elements
        .into_iter()
        .map(|element| element.borrow(py).kind.clone())
        .collect()
}

pub(super) fn right_projective_module(
    algebra: &PyPathAlgebra,
    vertices: Vec<usize>,
) -> PyResult<PyRightProjectiveModule> {
    let kind = map_ref!(
        &algebra.kind,
        AlgebraKind => ProjectiveModuleKind,
        |algebra| projective_module(algebra, &vertices)?
    );
    Ok(PyRightProjectiveModule { kind })
}

#[pymethods]
impl PyRightProjectiveModule {
    #[staticmethod]
    fn from_vertices(algebra: &PyPathAlgebra, vertices: Vec<usize>) -> PyResult<Self> {
        right_projective_module(algebra, vertices)
    }

    #[getter]
    fn algebra(&self) -> PyPathAlgebra {
        let kind = map_ref!(
            &self.kind,
            ProjectiveModuleKind => AlgebraKind,
            |module| module.algebra().clone()
        );
        PyPathAlgebra { kind }
    }

    #[getter]
    fn vertex_indices(&self) -> Vec<usize> {
        dispatch_ref!(&self.kind, ProjectiveModuleKind, |module| module
            .vertex_indices()
            .iter()
            .map(|index| index.get())
            .collect())
    }

    fn zero(&self) -> PyPathAlgebraVector {
        let kind = map_ref!(
            &self.kind,
            ProjectiveModuleKind => ProjectiveVectorKind,
            |module| module.zero()
        );
        PyPathAlgebraVector { kind }
    }

    #[getter]
    fn generators(&self) -> PyResult<Vec<PyPathAlgebraVector>> {
        dispatch_ref!(&self.kind, ProjectiveModuleKind, |module| Ok(wrap_vectors(
            module.generators().map_err(invalid)?
        )))
    }

    fn canonical_basis(&self) -> PyResult<Option<Vec<PyPathAlgebraVector>>> {
        dispatch_ref!(&self.kind, ProjectiveModuleKind, |module| Ok(module
            .canonical_basis()
            .map_err(invalid)?
            .map(wrap_vectors)))
    }

    fn vector(
        &self,
        py: Python<'_>,
        components: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyPathAlgebraVector> {
        let components = clone_element_kinds(py, components);
        let kind = map_ref!(
            &self.kind,
            ProjectiveModuleKind => ProjectiveVectorKind,
            |module| module
                .vector(runtime::extract_all(components)?)
                .map_err(invalid)?
        );
        Ok(PyPathAlgebraVector { kind })
    }

    fn vector_nc(
        &self,
        py: Python<'_>,
        components: Vec<Py<PyAlgebraElement>>,
        leading_position: usize,
        normal: bool,
    ) -> PyResult<PyPathAlgebraVector> {
        let components = clone_element_kinds(py, components);
        let kind = map_ref!(
            &self.kind,
            ProjectiveModuleKind => ProjectiveVectorKind,
            |module| module.vector_nc(
                runtime::extract_all(components)?,
                leading_position,
                normal,
            )
        );
        Ok(PyPathAlgebraVector { kind })
    }

    #[getter]
    fn dimension(&self) -> Option<usize> {
        dispatch_ref!(&self.kind, ProjectiveModuleKind, |module| module
            .dimension())
    }

    #[getter]
    fn is_vertex_projective(&self) -> bool {
        dispatch_ref!(&self.kind, ProjectiveModuleKind, |module| module
            .is_vertex_projective_module())
    }

    fn lift(&self) -> Self {
        let kind = map_ref!(
            &self.kind,
            ProjectiveModuleKind => ProjectiveModuleKind,
            |module| module.lift_path_algebra_module()
        );
        Self { kind }
    }

    fn uniform_generators(&self) -> PyResult<Vec<PyPathAlgebraVector>> {
        dispatch_ref!(&self.kind, ProjectiveModuleKind, |module| Ok(wrap_vectors(
            module.uniform_generators().map_err(invalid)?
        )))
    }

    fn right_groebner_basis(&self) -> PyResult<PyPathAlgebraModuleGroebnerBasis> {
        let kind = map_ref!(
            &self.kind,
            ProjectiveModuleKind => ProjectiveGroebnerKind,
            |module| module.right_groebner_basis().map_err(invalid)?
        );
        Ok(PyPathAlgebraModuleGroebnerBasis { kind })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            ProjectiveModuleKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, ProjectiveModuleKind, |module| format!(
            "<right-module over <{}[<quiver with {} vertices and {} arrows>]>>",
            module.algebra().field().python_kind().display(),
            module.algebra().quiver().number_of_vertices(),
            module.algebra().quiver().number_of_arrows(),
        ))
    }
}

#[pymethods]
impl PyPathAlgebraVector {
    #[getter]
    fn module(&self) -> PyRightProjectiveModule {
        let kind = map_ref!(
            &self.kind,
            ProjectiveVectorKind => ProjectiveModuleKind,
            |vector| vector.module().clone()
        );
        PyRightProjectiveModule { kind }
    }

    #[getter]
    fn components(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, ProjectiveVectorKind, |vector| {
            wrap_algebra_elements(vector.components().to_vec())
        })
    }

    #[getter]
    fn leading_position(&self) -> usize {
        dispatch_ref!(&self.kind, ProjectiveVectorKind, |vector| vector
            .leading_position())
    }

    #[getter]
    fn is_normal_form(&self) -> bool {
        dispatch_ref!(&self.kind, ProjectiveVectorKind, |vector| vector
            .is_normal_form())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, ProjectiveVectorKind, |vector| vector.is_zero())
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            ProjectiveVectorKind => ProjectiveVectorKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("path-algebra vectors use different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn subtract(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            ProjectiveVectorKind => ProjectiveVectorKind,
            |left, right| left.sub(right).map_err(invalid)?,
            invalid("path-algebra vectors use different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn negate(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ProjectiveVectorKind => ProjectiveVectorKind,
            |vector| vector.neg().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_left(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ProjectiveVectorKind => ProjectiveVectorKind,
            |vector| vector
                .scale_left(&scalar.value_for(vector.module().algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn act(&self, element: &PyAlgebraElement) -> PyResult<Self> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            ProjectiveVectorKind => ProjectiveVectorKind,
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
            ProjectiveVectorKind => ProjectiveVectorKind,
            |vector| vector.leading_term()
        );
        Self { kind }
    }

    #[getter]
    fn leading_coefficient(&self) -> PyResult<PyScalar> {
        dispatch_ref!(&self.kind, ProjectiveVectorKind, |vector| Ok(vector
            .module()
            .algebra()
            .field()
            .wrap_scalar(&vector.leading_coefficient().map_err(invalid)?)))
    }

    #[getter]
    fn leading_component(&self) -> PyResult<PyAlgebraElement> {
        let kind = map_ref!(
            &self.kind,
            ProjectiveVectorKind => AlgebraElementKind,
            |vector| vector.leading_component().map_err(invalid)?.clone()
        );
        Ok(PyAlgebraElement { kind })
    }

    #[getter]
    fn target_vertex(&self) -> Option<PyPath> {
        dispatch_ref!(&self.kind, ProjectiveVectorKind, |vector| vector
            .target_vertex()
            .map(|inner| PyPath { inner }))
    }

    #[getter]
    fn source_vertex(&self) -> Option<PyPath> {
        dispatch_ref!(&self.kind, ProjectiveVectorKind, |vector| vector
            .source_vertex()
            .map(|inner| PyPath { inner }))
    }

    #[getter]
    fn is_left_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, ProjectiveVectorKind, |vector| vector
            .is_left_uniform())
    }

    #[getter]
    fn is_right_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, ProjectiveVectorKind, |vector| vector
            .is_right_uniform())
    }

    #[getter]
    fn is_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, ProjectiveVectorKind, |vector| vector
            .is_uniform())
    }

    fn compare(&self, right: &Self) -> PyResult<i8> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            ProjectiveVectorKind,
            |left, right| Ok(ordering_value(left.compare(right).map_err(invalid)?)),
            Err(invalid(
                "path-algebra vectors use different coefficient fields"
            ))
        )
    }

    fn is_left_divisible(&self, divisor: &Self) -> PyResult<bool> {
        dispatch_pair_ref!(
            &self.kind,
            &divisor.kind,
            ProjectiveVectorKind,
            |vector, divisor| vector.is_left_divisible(divisor).map_err(invalid),
            Err(invalid(
                "path-algebra vectors use different coefficient fields"
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
            ProjectiveVectorKind,
            |left, right| left == right,
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, ProjectiveVectorKind, |vector| vector
            .print_obj_with(|value| vector
                .module()
                .algebra()
                .field()
                .scalar_text(value)))
    }
}

#[pymethods]
impl PyPathAlgebraModuleGroebnerBasis {
    #[getter]
    fn parent_module(&self) -> PyRightProjectiveModule {
        let kind = map_ref!(
            &self.kind,
            ProjectiveGroebnerKind => ProjectiveModuleKind,
            |basis| basis.parent_module().clone()
        );
        PyRightProjectiveModule { kind }
    }

    #[getter]
    fn is_for_whole_module(&self) -> bool {
        dispatch_ref!(&self.kind, ProjectiveGroebnerKind, |basis| basis
            .is_for_whole_module())
    }

    #[getter]
    fn basis_vectors(&self) -> Vec<PyPathAlgebraVector> {
        dispatch_ref!(&self.kind, ProjectiveGroebnerKind, |basis| {
            wrap_vectors(basis.basis_vectors().to_vec())
        })
    }

    #[getter]
    fn is_right(&self) -> bool {
        dispatch_ref!(&self.kind, ProjectiveGroebnerKind, |basis| basis
            .side()
            .is_right())
    }

    #[getter]
    fn is_left(&self) -> bool {
        dispatch_ref!(&self.kind, ProjectiveGroebnerKind, |basis| basis
            .side()
            .is_left())
    }

    fn completely_reduce(&self, vector: &PyPathAlgebraVector) -> PyResult<PyPathAlgebraVector> {
        let vector = vector.kind.clone();
        let kind = map_ref!(
            &self.kind,
            ProjectiveGroebnerKind => ProjectiveVectorKind,
            |basis| basis
                .completely_reduce(&runtime::Extract::extract(vector)?)
                .map_err(invalid)?
        );
        Ok(PyPathAlgebraVector { kind })
    }

    fn completely_reduced_basis_vectors(&self) -> PyResult<Vec<PyPathAlgebraVector>> {
        dispatch_ref!(&self.kind, ProjectiveGroebnerKind, |basis| Ok(
            wrap_vectors(
                completely_reduce_path_algebra_module_groebner_basis(basis).map_err(invalid)?,
            )
        ))
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, ProjectiveGroebnerKind, |basis| basis
            .view_obj_with_basis(
                &basis
                    .parent_module()
                    .algebra()
                    .field()
                    .python_kind()
                    .display(),
                |value| basis.parent_module().algebra().field().scalar_text(value),
            ))
    }
}
