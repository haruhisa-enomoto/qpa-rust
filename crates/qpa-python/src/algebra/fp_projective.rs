//! Python wrappers for QPA's finitely presented path modules `P/N`.

use numerica::domains::{
    Field,
    finite_field::{Z2, Zp},
    rational::RationalField,
};
use pyo3::{Py, PyResult, Python, prelude::*};
use qpa_core::{
    FpPathAlgebraByImagesHomomorphism, FpPathAlgebraMatrixHomomorphism, FpPathAlgebraVector,
    NaturalSubmoduleHomomorphism, PathAlgebra, PathAlgebraElement, RightFpPathAlgebraModule,
    RightFpPathAlgebraSubmodule, end_of_fp_path_algebra_module, hom_of_fp_path_algebra_modules,
    natural_homomorphism_by_projective_submodule, vertex_projective_presentation,
    vertex_projective_presentation_with_vertices,
};

use super::projective::{ProjectiveModuleKind, ProjectiveSubmoduleKind, ProjectiveVectorKind};
#[allow(clippy::wildcard_imports)]
use super::*;

#[derive(Clone)]
enum FpModuleKind {
    Rational(RightFpPathAlgebraModule<RationalField>),
    Binary(RightFpPathAlgebraModule<Z2>),
    Prime(RightFpPathAlgebraModule<Zp>),
}

#[pyclass(name = "RightFpPathAlgebraModule", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyRightFpPathAlgebraModule {
    kind: FpModuleKind,
}

#[derive(Clone)]
enum FpVectorKind {
    Rational(FpPathAlgebraVector<RationalField>),
    Binary(FpPathAlgebraVector<Z2>),
    Prime(FpPathAlgebraVector<Zp>),
}

#[pyclass(name = "FpPathAlgebraVector", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyFpPathAlgebraVector {
    kind: FpVectorKind,
}

#[derive(Clone)]
enum FpSubmoduleKind {
    Rational(RightFpPathAlgebraSubmodule<RationalField>),
    Binary(RightFpPathAlgebraSubmodule<Z2>),
    Prime(RightFpPathAlgebraSubmodule<Zp>),
}

#[pyclass(name = "RightFpPathAlgebraSubmodule", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyRightFpPathAlgebraSubmodule {
    kind: FpSubmoduleKind,
}

#[derive(Clone)]
enum NaturalHomomorphismKind {
    Rational(NaturalSubmoduleHomomorphism<RationalField>),
    Binary(NaturalSubmoduleHomomorphism<Z2>),
    Prime(NaturalSubmoduleHomomorphism<Zp>),
}

#[pyclass(name = "NaturalSubmoduleHomomorphism", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyNaturalSubmoduleHomomorphism {
    kind: NaturalHomomorphismKind,
}

#[derive(Clone)]
enum FpByImagesHomomorphismKind {
    Rational(FpPathAlgebraByImagesHomomorphism<RationalField>),
    Binary(FpPathAlgebraByImagesHomomorphism<Z2>),
    Prime(FpPathAlgebraByImagesHomomorphism<Zp>),
}

/// QPA's by-images representation used by its first `ImagesSource` overload.
#[pyclass(
    name = "FpPathAlgebraByImagesHomomorphism",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PyFpPathAlgebraByImagesHomomorphism {
    kind: FpByImagesHomomorphismKind,
}

#[derive(Clone)]
enum FpMatrixHomomorphismKind {
    Rational(FpPathAlgebraMatrixHomomorphism<RationalField>),
    Binary(FpPathAlgebraMatrixHomomorphism<Z2>),
    Prime(FpPathAlgebraMatrixHomomorphism<Zp>),
}

/// QPA's matrix representation used by its second `ImagesSource` overload.
#[pyclass(name = "FpPathAlgebraMatrixHomomorphism", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyFpPathAlgebraMatrixHomomorphism {
    kind: FpMatrixHomomorphismKind,
}

macro_rules! impl_fp_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<RightFpPathAlgebraModule<$field>> for FpModuleKind {
            fn wrap(value: RightFpPathAlgebraModule<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<FpPathAlgebraVector<$field>> for FpVectorKind {
            fn wrap(value: FpPathAlgebraVector<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Extract<FpPathAlgebraVector<$field>> for FpVectorKind {
            fn extract(self) -> PyResult<FpPathAlgebraVector<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "finitely presented vector belongs to another coefficient field",
                    )),
                }
            }
        }

        impl runtime::Wrap<RightFpPathAlgebraSubmodule<$field>> for FpSubmoduleKind {
            fn wrap(value: RightFpPathAlgebraSubmodule<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<NaturalSubmoduleHomomorphism<$field>> for NaturalHomomorphismKind {
            fn wrap(value: NaturalSubmoduleHomomorphism<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<FpPathAlgebraByImagesHomomorphism<$field>>
            for FpByImagesHomomorphismKind
        {
            fn wrap(value: FpPathAlgebraByImagesHomomorphism<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<FpPathAlgebraMatrixHomomorphism<$field>> for FpMatrixHomomorphismKind {
            fn wrap(value: FpPathAlgebraMatrixHomomorphism<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_fp_runtime!(RationalField, Rational);
impl_fp_runtime!(Z2, Binary);
impl_fp_runtime!(Zp, Prime);

fn wrap_vectors<F: Field>(values: Vec<FpPathAlgebraVector<F>>) -> Vec<PyFpPathAlgebraVector>
where
    FpVectorKind: runtime::Wrap<FpPathAlgebraVector<F>>,
{
    runtime::wrap_all(values)
        .into_iter()
        .map(|kind| PyFpPathAlgebraVector { kind })
        .collect()
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

fn clone_element_matrix(
    py: Python<'_>,
    rows: Vec<Vec<Py<PyAlgebraElement>>>,
) -> Vec<Vec<AlgebraElementKind>> {
    rows.into_iter()
        .map(|row| clone_element_kinds(py, row))
        .collect()
}

fn extract_element_matrix<F: Field>(
    rows: Vec<Vec<AlgebraElementKind>>,
) -> PyResult<Vec<Vec<PathAlgebraElement<F>>>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    rows.into_iter().map(runtime::extract_all).collect()
}

fn embedded_vertices<F: Field>(
    algebra: &PathAlgebra<F>,
    vertices: &[usize],
) -> PyResult<Vec<PathAlgebraElement<F>>> {
    let count = algebra.quiver().number_of_vertices();
    let generators = algebra.generators();
    vertices
        .iter()
        .map(|&position| {
            if position == 0 || position > count {
                Err(invalid(format!(
                    "vertex index {position} is outside 1..={count}"
                )))
            } else {
                Ok(generators[position - 1].clone())
            }
        })
        .collect()
}

fn clone_vector_kinds(
    py: Python<'_>,
    vectors: Vec<Py<PyFpPathAlgebraVector>>,
) -> Vec<FpVectorKind> {
    vectors
        .into_iter()
        .map(|vector| vector.borrow(py).kind.clone())
        .collect()
}

fn scalar_matrix<F: runtime::PythonField>(
    py: Python<'_>,
    field: &F,
    rows: Vec<Vec<Py<PyScalar>>>,
) -> PyResult<DenseMatrix<F>> {
    let rows = rows
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|value| field.scalar_value(&value.borrow(py)))
                .collect()
        })
        .collect::<PyResult<Vec<_>>>()?;
    DenseMatrix::from_rows(rows, field.clone()).map_err(invalid)
}

fn wrap_scalar_matrix<F: runtime::PythonField>(matrix: &DenseMatrix<F>) -> Vec<Vec<PyScalar>> {
    (0..matrix.nrows())
        .map(|row| {
            matrix
                .row(row)
                .iter()
                .map(|value| matrix.field().wrap_scalar(value))
                .collect()
        })
        .collect()
}

pub(super) fn hom(
    algebra: &PyPathAlgebra,
    source: &PyRightFpPathAlgebraModule,
    range: &PyRightFpPathAlgebraModule,
) -> PyResult<()> {
    match (&algebra.kind, &source.kind, &range.kind) {
        (
            AlgebraKind::Rational(algebra),
            FpModuleKind::Rational(source),
            FpModuleKind::Rational(range),
        ) => {
            hom_of_fp_path_algebra_modules(algebra, source, range).map_err(invalid)?;
        }
        (
            AlgebraKind::Binary(algebra),
            FpModuleKind::Binary(source),
            FpModuleKind::Binary(range),
        ) => {
            hom_of_fp_path_algebra_modules(algebra, source, range).map_err(invalid)?;
        }
        (AlgebraKind::Prime(algebra), FpModuleKind::Prime(source), FpModuleKind::Prime(range)) => {
            hom_of_fp_path_algebra_modules(algebra, source, range).map_err(invalid)?;
        }
        _ => return Err(invalid("Hom inputs use different coefficient fields")),
    }
    Ok(())
}

pub(super) fn end(algebra: &PyPathAlgebra, module: &PyRightFpPathAlgebraModule) -> PyResult<()> {
    match (&algebra.kind, &module.kind) {
        (AlgebraKind::Rational(algebra), FpModuleKind::Rational(module)) => {
            end_of_fp_path_algebra_module(algebra, module).map_err(invalid)?;
        }
        (AlgebraKind::Binary(algebra), FpModuleKind::Binary(module)) => {
            end_of_fp_path_algebra_module(algebra, module).map_err(invalid)?;
        }
        (AlgebraKind::Prime(algebra), FpModuleKind::Prime(module)) => {
            end_of_fp_path_algebra_module(algebra, module).map_err(invalid)?;
        }
        _ => return Err(invalid("End inputs use different coefficient fields")),
    }
    Ok(())
}

pub(super) fn natural_quotient(
    parent: &PyRightProjectiveModule,
    submodule: &PyRightProjectiveSubmodule,
) -> PyResult<PyNaturalSubmoduleHomomorphism> {
    let kind = match (&parent.kind, &submodule.kind) {
        (ProjectiveModuleKind::Rational(parent), ProjectiveSubmoduleKind::Rational(submodule)) => {
            NaturalHomomorphismKind::Rational(
                natural_homomorphism_by_projective_submodule(parent, submodule).map_err(invalid)?,
            )
        }
        (ProjectiveModuleKind::Binary(parent), ProjectiveSubmoduleKind::Binary(submodule)) => {
            NaturalHomomorphismKind::Binary(
                natural_homomorphism_by_projective_submodule(parent, submodule).map_err(invalid)?,
            )
        }
        (ProjectiveModuleKind::Prime(parent), ProjectiveSubmoduleKind::Prime(submodule)) => {
            NaturalHomomorphismKind::Prime(
                natural_homomorphism_by_projective_submodule(parent, submodule).map_err(invalid)?,
            )
        }
        _ => {
            return Err(invalid(
                "the projective module and submodule use different coefficient fields",
            ));
        }
    };
    Ok(PyNaturalSubmoduleHomomorphism { kind })
}

pub(super) fn vertex_presentation(
    algebra: &PyPathAlgebra,
    py: Python<'_>,
    rows: Vec<Vec<Py<PyAlgebraElement>>>,
) -> PyResult<PyRightFpPathAlgebraModule> {
    let rows = clone_element_matrix(py, rows);
    let kind = map_ref!(
        &algebra.kind,
        AlgebraKind => FpModuleKind,
        |algebra| vertex_projective_presentation(algebra, &extract_element_matrix(rows)?)
            .map_err(invalid)?
    );
    Ok(PyRightFpPathAlgebraModule { kind })
}

pub(super) fn vertex_presentation_with_vertices(
    algebra: &PyPathAlgebra,
    py: Python<'_>,
    rows: Vec<Vec<Py<PyAlgebraElement>>>,
    vertices: Vec<usize>,
) -> PyResult<PyRightFpPathAlgebraModule> {
    let rows = clone_element_matrix(py, rows);
    let kind = map_ref!(
        &algebra.kind,
        AlgebraKind => FpModuleKind,
        |algebra| vertex_projective_presentation_with_vertices(
            algebra,
            &extract_element_matrix(rows)?,
            &embedded_vertices(algebra, &vertices)?,
        )
        .map_err(invalid)?
    );
    Ok(PyRightFpPathAlgebraModule { kind })
}

#[pymethods]
impl PyNaturalSubmoduleHomomorphism {
    #[getter]
    fn source(&self) -> PyRightProjectiveModule {
        dispatch_ref!(&self.kind, NaturalHomomorphismKind, |mapping| {
            projective::wrap_projective_module(mapping.source().clone())
        })
    }

    #[getter]
    fn range(&self) -> PyRightFpPathAlgebraModule {
        let kind = map_ref!(
            &self.kind,
            NaturalHomomorphismKind => FpModuleKind,
            |mapping| mapping.range().clone()
        );
        PyRightFpPathAlgebraModule { kind }
    }

    #[getter]
    fn kernel(&self) -> PyRightProjectiveSubmodule {
        dispatch_ref!(&self.kind, NaturalHomomorphismKind, |mapping| {
            projective::wrap_projective_submodule(mapping.kernel().clone())
        })
    }

    #[getter]
    fn is_surjective(&self) -> bool {
        dispatch_ref!(&self.kind, NaturalHomomorphismKind, |mapping| mapping
            .is_surjective())
    }

    fn image(&self, vector: &PyPathAlgebraVector) -> PyResult<PyFpPathAlgebraVector> {
        let kind = match (&self.kind, &vector.kind) {
            (
                NaturalHomomorphismKind::Rational(mapping),
                ProjectiveVectorKind::Rational(vector),
            ) => FpVectorKind::Rational(mapping.image(vector).map_err(invalid)?),
            (NaturalHomomorphismKind::Binary(mapping), ProjectiveVectorKind::Binary(vector)) => {
                FpVectorKind::Binary(mapping.image(vector).map_err(invalid)?)
            }
            (NaturalHomomorphismKind::Prime(mapping), ProjectiveVectorKind::Prime(vector)) => {
                FpVectorKind::Prime(mapping.image(vector).map_err(invalid)?)
            }
            _ => {
                return Err(invalid(
                    "the mapping and vector use different coefficient fields",
                ));
            }
        };
        Ok(PyFpPathAlgebraVector { kind })
    }
}

#[pymethods]
impl PyFpPathAlgebraByImagesHomomorphism {
    #[staticmethod]
    fn from_images(
        py: Python<'_>,
        range: &PyRightFpPathAlgebraModule,
        images: Vec<Py<PyFpPathAlgebraVector>>,
    ) -> PyResult<Self> {
        let images = clone_vector_kinds(py, images);
        let kind = map_ref!(
            &range.kind,
            FpModuleKind => FpByImagesHomomorphismKind,
            |range| FpPathAlgebraByImagesHomomorphism::new(
                range,
                runtime::extract_all(images)?,
            )
            .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[getter]
    fn range(&self) -> PyRightFpPathAlgebraModule {
        let kind = map_ref!(
            &self.kind,
            FpByImagesHomomorphismKind => FpModuleKind,
            |mapping| mapping.range().clone()
        );
        PyRightFpPathAlgebraModule { kind }
    }

    #[getter]
    fn images(&self) -> Vec<PyFpPathAlgebraVector> {
        dispatch_ref!(&self.kind, FpByImagesHomomorphismKind, |mapping| {
            wrap_vectors(mapping.images().to_vec())
        })
    }

    #[getter]
    fn has_basis_image(&self) -> bool {
        dispatch_ref!(&self.kind, FpByImagesHomomorphismKind, |mapping| mapping
            .has_basis_image())
    }

    fn images_source(&self) -> PyResult<PyRightFpPathAlgebraSubmodule> {
        let kind = map_ref!(
            &self.kind,
            FpByImagesHomomorphismKind => FpSubmoduleKind,
            |mapping| mapping.images_source().map_err(invalid)?
        );
        Ok(PyRightFpPathAlgebraSubmodule { kind })
    }
}

#[pymethods]
impl PyFpPathAlgebraMatrixHomomorphism {
    #[staticmethod]
    fn from_matrix(
        py: Python<'_>,
        range: &PyRightFpPathAlgebraModule,
        matrix: Vec<Vec<Py<PyScalar>>>,
        basis_range: Vec<Py<PyFpPathAlgebraVector>>,
    ) -> PyResult<Self> {
        let basis_range = clone_vector_kinds(py, basis_range);
        let kind = map_ref!(
            &range.kind,
            FpModuleKind => FpMatrixHomomorphismKind,
            |range| FpPathAlgebraMatrixHomomorphism::new(
                range,
                scalar_matrix(py, range.algebra().field(), matrix)?,
                runtime::extract_all(basis_range)?,
            )
            .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[getter]
    fn range(&self) -> PyRightFpPathAlgebraModule {
        let kind = map_ref!(
            &self.kind,
            FpMatrixHomomorphismKind => FpModuleKind,
            |mapping| mapping.range().clone()
        );
        PyRightFpPathAlgebraModule { kind }
    }

    #[getter]
    fn matrix(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, FpMatrixHomomorphismKind, |mapping| {
            wrap_scalar_matrix(mapping.matrix())
        })
    }

    #[getter]
    fn basis_range(&self) -> Vec<PyFpPathAlgebraVector> {
        dispatch_ref!(&self.kind, FpMatrixHomomorphismKind, |mapping| {
            wrap_vectors(mapping.basis_range().to_vec())
        })
    }

    #[getter]
    fn has_basis_image(&self) -> bool {
        dispatch_ref!(&self.kind, FpMatrixHomomorphismKind, |mapping| mapping
            .has_basis_image())
    }

    fn images_source(&self) -> PyResult<PyRightFpPathAlgebraSubmodule> {
        let kind = map_ref!(
            &self.kind,
            FpMatrixHomomorphismKind => FpSubmoduleKind,
            |mapping| mapping.images_source().map_err(invalid)?
        );
        Ok(PyRightFpPathAlgebraSubmodule { kind })
    }
}

#[pymethods]
impl PyRightFpPathAlgebraModule {
    #[getter]
    fn algebra(&self) -> PyPathAlgebra {
        let kind = map_ref!(&self.kind, FpModuleKind => AlgebraKind, |module| module
            .algebra()
            .clone());
        PyPathAlgebra { kind }
    }

    #[getter]
    fn parent_projective(&self) -> PyRightProjectiveModule {
        dispatch_ref!(&self.kind, FpModuleKind, |module| {
            projective::wrap_projective_module(module.parent_projective().clone())
        })
    }

    #[getter]
    fn bottom_submodule(&self) -> PyRightProjectiveSubmodule {
        dispatch_ref!(&self.kind, FpModuleKind, |module| {
            projective::wrap_projective_submodule(module.bottom_submodule().clone())
        })
    }

    fn zero(&self) -> PyFpPathAlgebraVector {
        let kind = map_ref!(&self.kind, FpModuleKind => FpVectorKind, |module| module
            .zero());
        PyFpPathAlgebraVector { kind }
    }

    #[getter]
    fn generators(&self) -> PyResult<Vec<PyFpPathAlgebraVector>> {
        dispatch_ref!(&self.kind, FpModuleKind, |module| Ok(wrap_vectors(
            module.generators().map_err(invalid)?
        )))
    }

    #[getter]
    fn basis(&self) -> Vec<PyFpPathAlgebraVector> {
        dispatch_ref!(&self.kind, FpModuleKind, |module| wrap_vectors(
            module.basis()
        ))
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, FpModuleKind, |module| module.dimension())
    }

    #[getter]
    fn is_fp_path_algebra_module(&self) -> bool {
        dispatch_ref!(&self.kind, FpModuleKind, |module| module
            .is_fp_path_algebra_module())
    }

    #[getter]
    fn is_path_algebra_module(&self) -> bool {
        dispatch_ref!(&self.kind, FpModuleKind, |module| module
            .is_path_algebra_module())
    }

    #[getter]
    fn is_whole_family(&self) -> bool {
        dispatch_ref!(&self.kind, FpModuleKind, |module| module.is_whole_family())
    }

    fn vector(
        &self,
        py: Python<'_>,
        components: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyFpPathAlgebraVector> {
        let components = clone_element_kinds(py, components);
        let kind = map_ref!(
            &self.kind,
            FpModuleKind => FpVectorKind,
            |module| module
                .vector(runtime::extract_all(components)?)
                .map_err(invalid)?
        );
        Ok(PyFpPathAlgebraVector { kind })
    }

    fn vector_nc(
        &self,
        py: Python<'_>,
        components: Vec<Py<PyAlgebraElement>>,
        leading_position: usize,
        normal: bool,
    ) -> PyResult<PyFpPathAlgebraVector> {
        let components = clone_element_kinds(py, components);
        let kind = map_ref!(
            &self.kind,
            FpModuleKind => FpVectorKind,
            |module| module
                .vector_nc(
                    runtime::extract_all(components)?,
                    leading_position,
                    normal,
                )
                .map_err(invalid)?
        );
        Ok(PyFpPathAlgebraVector { kind })
    }

    fn submodule(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyFpPathAlgebraVector>>,
    ) -> PyResult<PyRightFpPathAlgebraSubmodule> {
        let generators = clone_vector_kinds(py, generators);
        let kind = map_ref!(
            &self.kind,
            FpModuleKind => FpSubmoduleKind,
            |module| module
                .submodule(runtime::extract_all(generators)?)
                .map_err(invalid)?
        );
        Ok(PyRightFpPathAlgebraSubmodule { kind })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            FpModuleKind,
            |left, right| left.is_identical(right),
            false
        )
    }
}

#[pymethods]
impl PyRightFpPathAlgebraSubmodule {
    #[getter]
    fn parent(&self) -> PyRightFpPathAlgebraModule {
        let kind = map_ref!(&self.kind, FpSubmoduleKind => FpModuleKind, |module| module
            .parent()
            .clone());
        PyRightFpPathAlgebraModule { kind }
    }

    #[getter]
    fn generators(&self) -> Vec<PyFpPathAlgebraVector> {
        dispatch_ref!(&self.kind, FpSubmoduleKind, |module| wrap_vectors(
            module.generators().to_vec()
        ))
    }

    #[getter]
    fn dimension(&self) -> PyResult<usize> {
        dispatch_ref!(&self.kind, FpSubmoduleKind, |module| module
            .dimension()
            .map_err(invalid))
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            FpSubmoduleKind,
            |left, right| left.is_identical(right),
            false
        )
    }
}

#[pymethods]
impl PyFpPathAlgebraVector {
    #[getter]
    fn module(&self) -> PyRightFpPathAlgebraModule {
        let kind = map_ref!(&self.kind, FpVectorKind => FpModuleKind, |vector| vector
            .module()
            .clone());
        PyRightFpPathAlgebraModule { kind }
    }

    #[getter]
    fn components(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, FpVectorKind, |vector| wrap_algebra_elements(
            vector.components().to_vec()
        ))
    }

    #[getter]
    fn leading_position(&self) -> usize {
        dispatch_ref!(&self.kind, FpVectorKind, |vector| vector.leading_position())
    }

    #[getter]
    fn is_normal_form(&self) -> bool {
        dispatch_ref!(&self.kind, FpVectorKind, |vector| vector.is_normal_form())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, FpVectorKind, |vector| vector.is_zero())
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            FpVectorKind => FpVectorKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("finitely presented vectors use different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn subtract(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            FpVectorKind => FpVectorKind,
            |left, right| left.sub(right).map_err(invalid)?,
            invalid("finitely presented vectors use different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn negate(&self) -> PyResult<Self> {
        let kind = map_ref!(&self.kind, FpVectorKind => FpVectorKind, |vector| vector
            .neg()
            .map_err(invalid)?);
        Ok(Self { kind })
    }

    fn scale_left(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(&self.kind, FpVectorKind => FpVectorKind, |vector| vector
            .scale_left(&scalar.value_for(vector.module().algebra().field())?)
            .map_err(invalid)?);
        Ok(Self { kind })
    }

    fn scale_right(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(&self.kind, FpVectorKind => FpVectorKind, |vector| vector
            .scale_right(&scalar.value_for(vector.module().algebra().field())?)
            .map_err(invalid)?);
        Ok(Self { kind })
    }

    fn act(&self, element: &PyAlgebraElement) -> PyResult<Self> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            FpVectorKind => FpVectorKind,
            |vector| vector
                .act(&runtime::Extract::extract(element)?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[getter]
    fn leading_term(&self) -> PyResult<Self> {
        let kind = map_ref!(&self.kind, FpVectorKind => FpVectorKind, |vector| vector
            .leading_term()
            .map_err(invalid)?);
        Ok(Self { kind })
    }

    #[getter]
    fn leading_coefficient(&self) -> PyResult<PyScalar> {
        dispatch_ref!(&self.kind, FpVectorKind, |vector| Ok(vector
            .module()
            .algebra()
            .field()
            .wrap_scalar(&vector.leading_coefficient().map_err(invalid)?)))
    }

    #[getter]
    fn leading_component(&self) -> PyResult<PyAlgebraElement> {
        let kind = map_ref!(
            &self.kind,
            FpVectorKind => AlgebraElementKind,
            |vector| vector.leading_component().map_err(invalid)?.clone()
        );
        Ok(PyAlgebraElement { kind })
    }

    #[getter]
    fn target_vertex(&self) -> Option<PyPath> {
        dispatch_ref!(&self.kind, FpVectorKind, |vector| vector
            .target_vertex()
            .map(|inner| PyPath { inner }))
    }

    #[getter]
    fn source_vertex(&self) -> Option<PyPath> {
        dispatch_ref!(&self.kind, FpVectorKind, |vector| vector
            .source_vertex()
            .map(|inner| PyPath { inner }))
    }

    #[getter]
    fn is_left_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, FpVectorKind, |vector| vector.is_left_uniform())
    }

    #[getter]
    fn is_right_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, FpVectorKind, |vector| vector.is_right_uniform())
    }

    #[getter]
    fn is_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, FpVectorKind, |vector| vector.is_uniform())
    }

    fn compare(&self, right: &Self) -> PyResult<i8> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            FpVectorKind,
            |left, right| Ok(ordering_value(left.compare(right).map_err(invalid)?)),
            Err(invalid(
                "finitely presented vectors use different coefficient fields"
            ))
        )
    }

    fn is_left_divisible(&self, divisor: &Self) -> PyResult<bool> {
        dispatch_pair_ref!(
            &self.kind,
            &divisor.kind,
            FpVectorKind,
            |vector, divisor| vector.is_left_divisible(divisor).map_err(invalid),
            Err(invalid(
                "finitely presented vectors use different coefficient fields"
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
            FpVectorKind,
            |left, right| left == right,
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, FpVectorKind, |vector| vector.print_obj_with(
            |value| vector.module().algebra().field().scalar_text(value)
        ))
    }

    fn print_obj(&self) -> String {
        self.__repr__()
    }
}
