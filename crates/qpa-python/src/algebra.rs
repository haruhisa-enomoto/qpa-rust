//! Python wrappers for path algebras, ideals, Gröbner bases, and quotients.

mod groebner;
mod ideal;
mod module;
mod quotient;

use numerica::domains::{
    finite_field::{Z2, Zp},
    integer::Integer,
    rational::{Q, RationalField},
};
use pyo3::{Py, PyRef, PyResult, Python, prelude::*};
use qpa_core::{
    AlgebraDimension, DenseMatrix, GroebnerBasis, GroebnerBasisKind, NamedArrowMatrix, PathAlgebra,
    PathAlgebraElement, PathAlgebraIdeal, PathAlgebraModule, PathAlgebraModuleElement,
    PathAlgebraModuleHomomorphism, PathAlgebraModuleSubmodule, QuiverAlgebra, QuiverAlgebraElement,
    QuotientPathAlgebra, QuotientPathAlgebraElement, right_module_over_path_algebra,
};

use crate::{
    FieldKind, PyField, PyScalar, RelationNotSatisfiedError, invalid,
    quiver::{PyPath, PyQuiver, ordering_value},
    runtime::{self, PythonField as _, dispatch_pair_ref, dispatch_ref, map_pair_ref, map_ref},
};
#[derive(Clone)]
enum AlgebraKind {
    Rational(PathAlgebra<RationalField>),
    Binary(PathAlgebra<Z2>),
    Prime(PathAlgebra<Zp>),
}

/// Runtime-dispatched path algebra backed by the existing generic Rust type.
#[pyclass(name = "PathAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPathAlgebra {
    kind: AlgebraKind,
}

#[derive(Clone)]
enum AlgebraElementKind {
    Rational(PathAlgebraElement<RationalField>),
    Binary(PathAlgebraElement<Z2>),
    Prime(PathAlgebraElement<Zp>),
}

macro_rules! impl_algebra_element_extract {
    ($field:ty, $variant:ident) => {
        impl runtime::Extract<PathAlgebraElement<$field>> for AlgebraElementKind {
            fn extract(self) -> PyResult<PathAlgebraElement<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid("relation belongs to another coefficient field")),
                }
            }
        }

        impl runtime::Wrap<PathAlgebraElement<$field>> for AlgebraElementKind {
            fn wrap(value: PathAlgebraElement<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_algebra_element_extract!(RationalField, Rational);
impl_algebra_element_extract!(Z2, Binary);
impl_algebra_element_extract!(Zp, Prime);

/// Sparse element of a runtime-dispatched Python path algebra.
#[pyclass(name = "AlgebraElement", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAlgebraElement {
    kind: AlgebraElementKind,
}

#[derive(Clone)]
enum IdealKind {
    Rational(PathAlgebraIdeal<RationalField>),
    Binary(PathAlgebraIdeal<Z2>),
    Prime(PathAlgebraIdeal<Zp>),
}

/// Two-sided path-algebra ideal using QPA's lazy Gröbner-basis behavior.
#[pyclass(name = "PathAlgebraIdeal", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPathAlgebraIdeal {
    kind: IdealKind,
}

#[derive(Clone)]
enum GroebnerKind {
    Rational(GroebnerBasis<RationalField>),
    Binary(GroebnerBasis<Z2>),
    Prime(GroebnerBasis<Zp>),
}

/// QPA's two-sided or right Gröbner-basis object.
#[pyclass(name = "GroebnerBasis", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyGroebnerBasis {
    kind: GroebnerKind,
}

macro_rules! impl_runtime_wrap {
    ($wrapper:ty, $value:ty, $variant:ident) => {
        impl runtime::Wrap<$value> for $wrapper {
            fn wrap(value: $value) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_runtime_wrap!(IdealKind, PathAlgebraIdeal<RationalField>, Rational);
impl_runtime_wrap!(IdealKind, PathAlgebraIdeal<Z2>, Binary);
impl_runtime_wrap!(IdealKind, PathAlgebraIdeal<Zp>, Prime);
impl_runtime_wrap!(GroebnerKind, GroebnerBasis<RationalField>, Rational);
impl_runtime_wrap!(GroebnerKind, GroebnerBasis<Z2>, Binary);
impl_runtime_wrap!(GroebnerKind, GroebnerBasis<Zp>, Prime);

#[derive(Clone)]
enum QuotientKind {
    Rational(QuotientPathAlgebra<RationalField>),
    Binary(QuotientPathAlgebra<Z2>),
    Prime(QuotientPathAlgebra<Zp>),
}

/// Quotient path algebra constructed by QPA's high-level Groebner route.
#[pyclass(name = "QuotientAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuotientAlgebra {
    kind: QuotientKind,
}

#[derive(Clone)]
enum QuotientElementKind {
    Rational(QuotientPathAlgebraElement<RationalField>),
    Binary(QuotientPathAlgebraElement<Z2>),
    Prime(QuotientPathAlgebraElement<Zp>),
}

macro_rules! impl_quotient_element_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Extract<QuotientPathAlgebraElement<$field>> for QuotientElementKind {
            fn extract(self) -> PyResult<QuotientPathAlgebraElement<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "quotient element belongs to another coefficient field",
                    )),
                }
            }
        }

        impl runtime::Wrap<QuotientPathAlgebraElement<$field>> for QuotientElementKind {
            fn wrap(value: QuotientPathAlgebraElement<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_quotient_element_runtime!(RationalField, Rational);
impl_quotient_element_runtime!(Z2, Binary);
impl_quotient_element_runtime!(Zp, Prime);

/// Element of a runtime-dispatched quotient path algebra.
#[pyclass(name = "QuotientElement", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuotientElement {
    kind: QuotientElementKind,
}

#[derive(Clone)]
enum MatrixModuleKind {
    Rational(PathAlgebraModule<RationalField>),
    Binary(PathAlgebraModule<Z2>),
    Prime(PathAlgebraModule<Zp>),
}

macro_rules! impl_module_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Extract<PathAlgebraModule<$field>> for MatrixModuleKind {
            fn extract(self) -> PyResult<PathAlgebraModule<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid("module belongs to another coefficient field")),
                }
            }
        }

        impl runtime::Wrap<PathAlgebraModule<$field>> for MatrixModuleKind {
            fn wrap(value: PathAlgebraModule<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_module_runtime!(RationalField, Rational);
impl_module_runtime!(Z2, Binary);
impl_module_runtime!(Zp, Prime);

/// Checked QPA right matrix module over a path algebra or quotient.
#[pyclass(name = "MatrixModule", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyMatrixModule {
    kind: MatrixModuleKind,
}

#[derive(Clone)]
enum ModuleElementKind {
    Rational(PathAlgebraModuleElement<RationalField>),
    Binary(PathAlgebraModuleElement<Z2>),
    Prime(PathAlgebraModuleElement<Zp>),
}

/// Element of a QPA matrix module over a path algebra or quotient.
#[pyclass(name = "ModuleElement", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleElement {
    kind: ModuleElementKind,
}

#[derive(Clone)]
enum ModuleHomomorphismKind {
    Rational(PathAlgebraModuleHomomorphism<RationalField>),
    Binary(PathAlgebraModuleHomomorphism<Z2>),
    Prime(PathAlgebraModuleHomomorphism<Zp>),
}

/// QPA homomorphism between right matrix modules.
#[pyclass(name = "ModuleHomomorphism", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleHomomorphism {
    kind: ModuleHomomorphismKind,
}

#[derive(Clone)]
enum ModuleSubmoduleKind {
    Rational(PathAlgebraModuleSubmodule<RationalField>),
    Binary(PathAlgebraModuleSubmodule<Z2>),
    Prime(PathAlgebraModuleSubmodule<Zp>),
}

/// QPA's generator-defined concrete submodule of a path-algebra matrix module.
#[pyclass(name = "ModuleSubmodule", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleSubmodule {
    kind: ModuleSubmoduleKind,
}

macro_rules! impl_module_value_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Extract<PathAlgebraModuleElement<$field>> for ModuleElementKind {
            fn extract(self) -> PyResult<PathAlgebraModuleElement<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "module element belongs to another coefficient field",
                    )),
                }
            }
        }

        impl runtime::Wrap<PathAlgebraModuleElement<$field>> for ModuleElementKind {
            fn wrap(value: PathAlgebraModuleElement<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<PathAlgebraModuleHomomorphism<$field>> for ModuleHomomorphismKind {
            fn wrap(value: PathAlgebraModuleHomomorphism<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<PathAlgebraModuleSubmodule<$field>> for ModuleSubmoduleKind {
            fn wrap(value: PathAlgebraModuleSubmodule<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_module_value_runtime!(RationalField, Rational);
impl_module_value_runtime!(Z2, Binary);
impl_module_value_runtime!(Zp, Prime);

fn clone_relations(
    py: Python<'_>,
    relations: Vec<Py<PyAlgebraElement>>,
) -> Vec<AlgebraElementKind> {
    relations
        .into_iter()
        .map(|relation| relation.borrow(py).kind.clone())
        .collect()
}

fn wrap_algebra_elements<F: runtime::PythonField>(
    elements: Vec<PathAlgebraElement<F>>,
) -> Vec<PyAlgebraElement>
where
    AlgebraElementKind: runtime::Wrap<PathAlgebraElement<F>>,
{
    runtime::wrap_all::<AlgebraElementKind, _>(elements)
        .into_iter()
        .map(|kind| PyAlgebraElement { kind })
        .collect()
}

fn wrap_groebner_basis<F: runtime::PythonField>(basis: GroebnerBasis<F>) -> PyGroebnerBasis
where
    GroebnerKind: runtime::Wrap<GroebnerBasis<F>>,
{
    PyGroebnerBasis {
        kind: runtime::Wrap::wrap(basis),
    }
}

fn wrap_quotient_elements<F: runtime::PythonField>(
    elements: Vec<QuotientPathAlgebraElement<F>>,
) -> Vec<PyQuotientElement>
where
    QuotientElementKind: runtime::Wrap<QuotientPathAlgebraElement<F>>,
{
    runtime::wrap_all::<QuotientElementKind, _>(elements)
        .into_iter()
        .map(|kind| PyQuotientElement { kind })
        .collect()
}

fn wrap_matrix_modules<F: runtime::PythonField>(
    modules: Vec<PathAlgebraModule<F>>,
) -> Vec<PyMatrixModule>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    runtime::wrap_all::<MatrixModuleKind, _>(modules)
        .into_iter()
        .map(|kind| PyMatrixModule { kind })
        .collect()
}

fn build_ideal<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    generators: Vec<AlgebraElementKind>,
) -> PyResult<PathAlgebraIdeal<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    PathAlgebraIdeal::new(algebra.clone(), runtime::extract_all(generators)?).map_err(invalid)
}

fn build_groebner_basis<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<GroebnerBasis<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    GroebnerBasis::high_level(algebra, &runtime::extract_all(relations)?).map_err(invalid)
}

#[pymethods]
impl PyPathAlgebra {
    #[new]
    fn new(field: PyRef<'_, PyField>, quiver: PyRef<'_, PyQuiver>) -> PyResult<Self> {
        let kind = match field.kind {
            FieldKind::Rational => AlgebraKind::Rational(
                PathAlgebra::new(Q, quiver.rust_quiver().clone()).map_err(invalid)?,
            ),
            FieldKind::Binary => AlgebraKind::Binary(
                PathAlgebra::new(Z2, quiver.rust_quiver().clone()).map_err(invalid)?,
            ),
            FieldKind::Prime(characteristic) => AlgebraKind::Prime(
                PathAlgebra::new(Zp::new(characteristic), quiver.rust_quiver().clone())
                    .map_err(invalid)?,
            ),
        };
        Ok(Self { kind })
    }

    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .field()
            .python_kind());
        PyField { kind }
    }

    #[getter]
    fn quiver(&self) -> PyQuiver {
        let inner = dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra.quiver().clone());
        PyQuiver { inner }
    }

    #[getter]
    fn global_dimension(&self) -> usize {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .qpa_global_dimension())
    }

    #[getter]
    fn has_finite_canonical_basis(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .has_finite_canonical_basis())
    }

    #[getter]
    fn is_monomial(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_monomial(algebra)
        })
    }

    #[getter]
    fn is_special_biserial(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_special_biserial(algebra)
        })
    }

    #[getter]
    fn is_string(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_string(algebra)
        })
    }

    #[getter]
    fn is_gentle(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_gentle(algebra)
        })
    }

    #[getter]
    fn generators(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| wrap_algebra_elements(
            algebra.generators()
        ))
    }

    #[getter]
    fn relations(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| wrap_algebra_elements(
            algebra.relations()
        ))
    }

    fn canonical_basis(&self) -> Option<Vec<PyAlgebraElement>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .canonical_basis()
            .map(wrap_algebra_elements))
    }

    fn opposite(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraKind,
            |algebra| algebra.opposite().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn generator(&self, name: &str) -> PyResult<PyAlgebraElement> {
        let path = dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .quiver()
            .generator_by_name(name)
            .ok_or_else(|| invalid(format!("unknown quiver generator {name:?}")))?
            .path());
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra.element(path).map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn term(&self, path: &PyPath, coefficient: &PyScalar) -> PyResult<PyAlgebraElement> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra
                .from_terms([(
                    path.inner.clone(),
                    coefficient.value_for(algebra.field())?,
                )])
                .map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn zero(&self) -> PyAlgebraElement {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra.zero()
        );
        PyAlgebraElement { kind }
    }

    fn one(&self) -> PyAlgebraElement {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra.one()
        );
        PyAlgebraElement { kind }
    }

    fn quotient(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyQuotientAlgebra> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => QuotientKind,
            |algebra| build_quotient(algebra, relations)?
        );
        Ok(PyQuotientAlgebra { kind })
    }

    fn quotient_high_level(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyQuotientAlgebra> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => QuotientKind,
            |algebra| build_high_level_quotient(algebra, relations)?
        );
        Ok(PyQuotientAlgebra { kind })
    }

    fn ideal(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyPathAlgebraIdeal> {
        let generators = clone_relations(py, generators);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => IdealKind,
            |algebra| build_ideal(algebra, generators)?
        );
        Ok(PyPathAlgebraIdeal { kind })
    }

    fn groebner_basis(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyGroebnerBasis> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => GroebnerKind,
            |algebra| build_groebner_basis(algebra, relations)?
        );
        Ok(PyGroebnerBasis { kind })
    }

    fn gbnp_groebner_basis(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyGroebnerBasis> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => GroebnerKind,
            |algebra| build_gbnp_groebner_basis(algebra, relations)?
        );
        Ok(PyGroebnerBasis { kind })
    }

    fn make_uniform_on_right(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<Vec<PyAlgebraElement>> {
        let generators = clone_relations(py, generators);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            let generators = runtime::extract_all(generators)?;
            Ok(wrap_algebra_elements(
                algebra
                    .make_uniform_on_right(&generators)
                    .map_err(invalid)?,
            ))
        })
    }

    fn generators_times_arrows_on_right(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<Vec<PyAlgebraElement>> {
        let generators = clone_relations(py, generators);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            let generators = runtime::extract_all(generators)?;
            Ok(wrap_algebra_elements(
                algebra
                    .generators_times_arrows_on_right(&generators)
                    .map_err(invalid)?,
            ))
        })
    }

    fn nth_power_of_arrow_ideal(&self, exponent: i64) -> PyResult<Vec<PyAlgebraElement>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| Ok(
            wrap_algebra_elements(
                algebra
                    .nth_power_of_arrow_ideal(exponent)
                    .map_err(invalid)?,
            )
        ))
    }

    fn module_from_integer_matrices(
        &self,
        dimension_vector: Vec<usize>,
        actions: Vec<(String, Vec<Vec<i64>>)>,
    ) -> PyResult<PyMatrixModule> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => MatrixModuleKind,
            |algebra| build_path_module(algebra, dimension_vector, actions)?
        );
        Ok(PyMatrixModule { kind })
    }

    fn zero_module(&self) -> PyResult<PyMatrixModule> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => MatrixModuleKind,
            |algebra| qpa_core::zero_module(&QuiverAlgebra::Path(algebra.clone()))
                .map_err(invalid)?
        );
        Ok(PyMatrixModule { kind })
    }

    fn simple_modules(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| Ok(wrap_matrix_modules(
            qpa_gbnp::simple_modules(&QuiverAlgebra::Path(algebra.clone())).map_err(invalid)?,
        )))
    }

    fn indecomposable_projectives(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| Ok(wrap_matrix_modules(
            qpa_gbnp::indec_projective_modules(&QuiverAlgebra::Path(algebra.clone()))
                .map_err(invalid)?,
        )))
    }

    fn indecomposable_injectives(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| Ok(wrap_matrix_modules(
            qpa_gbnp::indec_injective_modules(&QuiverAlgebra::Path(algebra.clone()))
                .map_err(invalid)?,
        )))
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            AlgebraKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .view_obj_with(&algebra.field().python_kind().display()))
    }
}

fn build_quotient<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<QuotientPathAlgebra<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let relations = runtime::extract_all(relations)?;
    qpa_gbnp::quotient_by_relations(algebra.clone(), &relations).map_err(invalid)
}

fn build_high_level_quotient<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<QuotientPathAlgebra<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let relations = runtime::extract_all(relations)?;
    QuotientPathAlgebra::by_relations_high_level(algebra.clone(), &relations).map_err(invalid)
}

fn build_gbnp_groebner_basis<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<GroebnerBasis<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let relations = runtime::extract_all(relations)?;
    let complete = qpa_gbnp::gbnp_groebner_basis(&relations, algebra).map_err(invalid)?;
    GroebnerBasis::from_complete_relations(algebra, &complete).map_err(invalid)
}

fn build_ideal_quotient<F: runtime::PythonField>(
    ideal: &PathAlgebraIdeal<F>,
) -> PyResult<QuotientPathAlgebra<F>> {
    if ideal.groebner_basis().is_none() {
        let complete = qpa_gbnp::gbnp_groebner_basis_of_ideal(ideal).map_err(invalid)?;
        ideal
            .with_complete_groebner_basis(&complete)
            .map_err(invalid)?;
    }
    QuotientPathAlgebra::by_ideal(ideal).map_err(invalid)
}

fn build_path_module<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    dimension_vector: Vec<usize>,
    actions: Vec<(String, Vec<Vec<i64>>)>,
) -> PyResult<PathAlgebraModule<F>> {
    build_matrix_module(
        QuiverAlgebra::Path(algebra.clone()),
        dimension_vector,
        actions,
    )
}

fn build_matrix_module<F: runtime::PythonField>(
    algebra: QuiverAlgebra<F>,
    dimension_vector: Vec<usize>,
    actions: Vec<(String, Vec<Vec<i64>>)>,
) -> PyResult<PathAlgebraModule<F>> {
    let actions = actions
        .into_iter()
        .map(|(label, rows)| {
            let rows = rows
                .into_iter()
                .map(|row| {
                    row.into_iter()
                        .map(|value| algebra.field().integer_value(&Integer::from(value)))
                        .collect()
                })
                .collect();
            let matrix = DenseMatrix::from_rows(rows, algebra.field().clone()).map_err(invalid)?;
            Ok(NamedArrowMatrix::new(label, matrix))
        })
        .collect::<PyResult<Vec<_>>>()?;
    right_module_over_path_algebra(algebra, dimension_vector, &actions)
        .map_err(invalid)?
        .ok_or_else(|| {
            RelationNotSatisfiedError::new_err("module actions do not satisfy relations")
        })
}

#[pymethods]
impl PyAlgebraElement {
    #[getter]
    fn algebra(&self) -> PyPathAlgebra {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraKind,
            |element| element.algebra().clone()
        );
        PyPathAlgebra { kind }
    }

    #[getter]
    fn terms(&self) -> Vec<(PyPath, PyScalar)> {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .terms()
            .iter()
            .map(|term| (
                PyPath {
                    inner: term.path().clone(),
                },
                element.algebra().field().wrap_scalar(term.coefficient()),
            ))
            .collect())
    }

    #[getter]
    fn external_representation(&self) -> Vec<(Vec<usize>, PyScalar)> {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .external_representation()
            .terms()
            .iter()
            .map(|(path, coefficient)| (
                path.clone(),
                element.algebra().field().wrap_scalar(coefficient),
            ))
            .collect())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element.is_zero())
    }

    #[getter]
    fn is_left_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .is_left_uniform())
    }

    #[getter]
    fn is_right_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .is_right_uniform())
    }

    #[getter]
    fn is_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .is_uniform())
    }

    #[getter]
    fn leading_term(&self) -> Self {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.leading_term()
        );
        Self { kind }
    }

    #[getter]
    fn leading_coefficient(&self) -> PyScalar {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .algebra()
            .field()
            .wrap_scalar(&element.leading_coefficient()))
    }

    #[getter]
    fn leading_monomial(&self) -> PyPath {
        let inner = dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .leading_monomial());
        PyPath { inner }
    }

    #[getter]
    fn vertex_position(&self) -> PyResult<usize> {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .vertex_position()
            .map_err(invalid))
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind => AlgebraElementKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn negate(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.neg().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn subtract(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind => AlgebraElementKind,
            |left, right| left.sub(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn multiply(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind => AlgebraElementKind,
            |left, right| left.multiply(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn pow(&self, exponent: u64) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.pow(exponent).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_left(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element
                .scale_left(&scalar.value_for(element.algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_right(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element
                .scale_right(&scalar.value_for(element.algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn divide_scalar(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element
                .divide_scalar(&scalar.value_for(element.algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn compare(&self, right: &Self) -> PyResult<i8> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind,
            |left, right| Ok(ordering_value(left.compare(right).map_err(invalid)?)),
            Err(invalid("elements belong to different coefficient fields"))
        )
    }

    fn opposite(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.opposite().map_err(invalid)?
        );
        Ok(Self { kind })
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

    fn __mul__(&self, right: &Self) -> PyResult<Self> {
        self.multiply(right)
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .print_obj_with(|value| element.algebra().field().scalar_text(value)))
    }

    fn __eq__(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            AlgebraElementKind,
            |left, right| left == right,
            false
        )
    }
}
