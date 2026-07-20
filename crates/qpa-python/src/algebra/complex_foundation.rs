//! Python wrappers for QPA's module category and lazy morphism lists.

use numerica::domains::{
    finite_field::{Z2, Zp},
    integer::Integer,
    rational::RationalField,
};
use pyo3::{Py, PyResult, Python, prelude::*};
use qpa_core::{InfiniteList, PathAlgebraModuleHomomorphism, QuiverAlgebra};
use qpa_gbnp::{ComplexError, QpaCategory, RightModuleCategory};

#[allow(clippy::wildcard_imports)]
use super::*;

fn clone_homomorphism_kinds(
    py: Python<'_>,
    morphisms: Vec<Py<PyModuleHomomorphism>>,
) -> Vec<ModuleHomomorphismKind> {
    morphisms
        .into_iter()
        .map(|morphism| morphism.borrow(py).kind.clone())
        .collect()
}

fn complex_half<F: runtime::PythonField + Send + Sync + 'static>(
    values: Option<Vec<ModuleHomomorphismKind>>,
) -> PyResult<qpa_gbnp::ComplexHalf<RightModuleCategory<F>>>
where
    F::Element: Send + Sync,
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>,
{
    Ok(match values {
        None => qpa_gbnp::ComplexHalf::Zero,
        Some(values) => qpa_gbnp::ComplexHalf::Repeat(runtime::extract_all(values)?),
    })
}

fn chain_map_half<F: runtime::PythonField + Send + Sync + 'static>(
    values: Option<Vec<ModuleHomomorphismKind>>,
) -> PyResult<qpa_gbnp::ChainMapHalf<RightModuleCategory<F>>>
where
    F::Element: Send + Sync,
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>,
{
    Ok(match values {
        None => qpa_gbnp::ChainMapHalf::Zero,
        Some(values) => qpa_gbnp::ChainMapHalf::Repeat(runtime::extract_all(values)?),
    })
}

#[derive(Clone)]
pub(super) enum ModuleCategoryKind {
    Rational(RightModuleCategory<RationalField>),
    Binary(RightModuleCategory<Z2>),
    Prime(RightModuleCategory<Zp>),
}

/// QPA's `CatOfRightAlgebraModules` record.
#[pyclass(name = "ModuleCategory", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleCategory {
    pub(super) kind: ModuleCategoryKind,
}

#[derive(Clone)]
pub(super) enum ModuleMorphismListKind {
    Rational(InfiniteList<PathAlgebraModuleHomomorphism<RationalField>, ComplexError>),
    Binary(InfiniteList<PathAlgebraModuleHomomorphism<Z2>, ComplexError>),
    Prime(InfiniteList<PathAlgebraModuleHomomorphism<Zp>, ComplexError>),
}

/// QPA's lazy `DifferentialsOfComplex`/`MorphismsOfChainMap` infinite list.
#[pyclass(name = "ModuleMorphismList", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleMorphismList {
    pub(super) kind: ModuleMorphismListKind,
}

macro_rules! impl_complex_foundation_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<RightModuleCategory<$field>> for ModuleCategoryKind {
            fn wrap(value: RightModuleCategory<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<InfiniteList<PathAlgebraModuleHomomorphism<$field>, ComplexError>>
            for ModuleMorphismListKind
        {
            fn wrap(
                value: InfiniteList<PathAlgebraModuleHomomorphism<$field>, ComplexError>,
            ) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_complex_foundation_runtime!(RationalField, Rational);
impl_complex_foundation_runtime!(Z2, Binary);
impl_complex_foundation_runtime!(Zp, Prime);

fn extended_integer_text(value: qpa_core::ExtendedInteger) -> String {
    match value {
        qpa_core::ExtendedInteger::NegativeInfinity => "-infinity".into(),
        qpa_core::ExtendedInteger::Finite(value) => value.to_string(),
        qpa_core::ExtendedInteger::PositiveInfinity => "infinity".into(),
    }
}

pub(super) fn category_of_complex(complex: &PyModuleComplex) -> PyModuleCategory {
    let kind = map_ref!(
        &complex.kind,
        ModuleComplexKind => ModuleCategoryKind,
        |complex| complex.category().clone()
    );
    PyModuleCategory { kind }
}

pub(super) fn differentials_of_complex(complex: &PyModuleComplex) -> PyModuleMorphismList {
    let kind = map_ref!(
        &complex.kind,
        ModuleComplexKind => ModuleMorphismListKind,
        |complex| complex.differentials()
    );
    PyModuleMorphismList { kind }
}

pub(super) fn complex_by_differential_list(
    category: &PyModuleCategory,
    differentials: &PyModuleMorphismList,
) -> PyResult<PyModuleComplex> {
    let kind = match (&category.kind, &differentials.kind) {
        (ModuleCategoryKind::Rational(category), ModuleMorphismListKind::Rational(list)) => {
            ModuleComplexKind::Rational(qpa_gbnp::QpaComplex::by_differential_list(
                category.clone(),
                list.clone(),
            ))
        }
        (ModuleCategoryKind::Binary(category), ModuleMorphismListKind::Binary(list)) => {
            ModuleComplexKind::Binary(qpa_gbnp::QpaComplex::by_differential_list(
                category.clone(),
                list.clone(),
            ))
        }
        (ModuleCategoryKind::Prime(category), ModuleMorphismListKind::Prime(list)) => {
            ModuleComplexKind::Prime(qpa_gbnp::QpaComplex::by_differential_list(
                category.clone(),
                list.clone(),
            ))
        }
        _ => {
            return Err(invalid(
                "the category and differential list use different coefficient fields",
            ));
        }
    };
    Ok(PyModuleComplex { kind })
}

fn complex_from_parts_for<F>(
    category: &RightModuleCategory<F>,
    base_degree: i64,
    middle: Vec<ModuleHomomorphismKind>,
    positive: Option<Vec<ModuleHomomorphismKind>>,
    negative: Option<Vec<ModuleHomomorphismKind>>,
) -> PyResult<PyModuleComplex>
where
    F: runtime::PythonField + Send + Sync + 'static,
    F::Element: Send + Sync,
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>,
    ModuleComplexKind: runtime::Wrap<qpa_gbnp::ModuleComplex<F>>,
{
    Ok(PyModuleComplex {
        kind: runtime::Wrap::wrap(
            qpa_gbnp::complex(
                category.clone(),
                Integer::from(base_degree),
                runtime::extract_all(middle)?,
                complex_half::<F>(positive)?,
                complex_half::<F>(negative)?,
            )
            .map_err(invalid)?,
        ),
    })
}

pub(super) fn complex_from_parts(
    category: &PyModuleCategory,
    py: Python<'_>,
    base_degree: i64,
    middle: Vec<Py<PyModuleHomomorphism>>,
    positive: Option<Vec<Py<PyModuleHomomorphism>>>,
    negative: Option<Vec<Py<PyModuleHomomorphism>>>,
) -> PyResult<PyModuleComplex> {
    let middle = clone_homomorphism_kinds(py, middle);
    let positive = positive.map(|values| clone_homomorphism_kinds(py, values));
    let negative = negative.map(|values| clone_homomorphism_kinds(py, values));
    match &category.kind {
        ModuleCategoryKind::Rational(category) => {
            complex_from_parts_for(category, base_degree, middle, positive, negative)
        }
        ModuleCategoryKind::Binary(category) => {
            complex_from_parts_for(category, base_degree, middle, positive, negative)
        }
        ModuleCategoryKind::Prime(category) => {
            complex_from_parts_for(category, base_degree, middle, positive, negative)
        }
    }
}

fn chain_map_from_parts_for<F>(
    source: &qpa_gbnp::ModuleComplex<F>,
    range: &qpa_gbnp::ModuleComplex<F>,
    base_degree: i64,
    middle: Vec<ModuleHomomorphismKind>,
    positive: Option<Vec<ModuleHomomorphismKind>>,
    negative: Option<Vec<ModuleHomomorphismKind>>,
) -> PyResult<PyModuleChainMap>
where
    F: runtime::PythonField + Send + Sync + 'static,
    F::Element: Send + Sync,
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>,
    ModuleChainMapKind: runtime::Wrap<qpa_gbnp::ModuleChainMap<F>>,
{
    Ok(PyModuleChainMap {
        kind: runtime::Wrap::wrap(
            qpa_gbnp::chain_map(
                source.clone(),
                range.clone(),
                Integer::from(base_degree),
                runtime::extract_all(middle)?,
                chain_map_half::<F>(positive)?,
                chain_map_half::<F>(negative)?,
            )
            .map_err(invalid)?,
        ),
    })
}

pub(super) fn chain_map_from_parts(
    source: &PyModuleComplex,
    range: &PyModuleComplex,
    py: Python<'_>,
    base_degree: i64,
    middle: Vec<Py<PyModuleHomomorphism>>,
    positive: Option<Vec<Py<PyModuleHomomorphism>>>,
    negative: Option<Vec<Py<PyModuleHomomorphism>>>,
) -> PyResult<PyModuleChainMap> {
    let middle = clone_homomorphism_kinds(py, middle);
    let positive = positive.map(|values| clone_homomorphism_kinds(py, values));
    let negative = negative.map(|values| clone_homomorphism_kinds(py, values));
    match (&source.kind, &range.kind) {
        (ModuleComplexKind::Rational(source), ModuleComplexKind::Rational(range)) => {
            chain_map_from_parts_for(source, range, base_degree, middle, positive, negative)
        }
        (ModuleComplexKind::Binary(source), ModuleComplexKind::Binary(range)) => {
            chain_map_from_parts_for(source, range, base_degree, middle, positive, negative)
        }
        (ModuleComplexKind::Prime(source), ModuleComplexKind::Prime(range)) => {
            chain_map_from_parts_for(source, range, base_degree, middle, positive, negative)
        }
        _ => Err(invalid(
            "source and range complexes use different coefficient fields",
        )),
    }
}

pub(super) fn morphisms_of_chain_map(map: &PyModuleChainMap) -> PyModuleMorphismList {
    let kind = map_ref!(
        &map.kind,
        ModuleChainMapKind => ModuleMorphismListKind,
        |map| map.morphisms()
    );
    PyModuleMorphismList { kind }
}

#[pymethods]
impl PyModuleCategory {
    #[staticmethod]
    fn from_path_algebra(algebra: &PyPathAlgebra) -> PyResult<Self> {
        let kind = map_ref!(
            &algebra.kind,
            AlgebraKind => ModuleCategoryKind,
            |algebra| qpa_gbnp::cat_of_right_algebra_modules(&QuiverAlgebra::Path(
                algebra.clone(),
            ))
            .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[staticmethod]
    fn from_quotient_algebra(algebra: &PyQuotientAlgebra) -> PyResult<Self> {
        let kind = map_ref!(
            &algebra.kind,
            QuotientKind => ModuleCategoryKind,
            |algebra| qpa_gbnp::cat_of_right_algebra_modules(&QuiverAlgebra::Quotient(
                algebra.clone(),
            ))
            .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_cat(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn name(&self) -> String {
        dispatch_ref!(&self.kind, ModuleCategoryKind, |category| category
            .name()
            .to_owned())
    }

    fn zero_object(&self) -> PyMatrixModule {
        let kind = map_ref!(
            &self.kind,
            ModuleCategoryKind => MatrixModuleKind,
            |category| category.zero_object()
        );
        PyMatrixModule { kind }
    }

    fn is_zero_object(&self, object: &PyMatrixModule) -> PyResult<bool> {
        match (&self.kind, &object.kind) {
            (ModuleCategoryKind::Rational(category), MatrixModuleKind::Rational(object)) => {
                Ok(category.is_zero_object(object))
            }
            (ModuleCategoryKind::Binary(category), MatrixModuleKind::Binary(object)) => {
                Ok(category.is_zero_object(object))
            }
            (ModuleCategoryKind::Prime(category), MatrixModuleKind::Prime(object)) => {
                Ok(category.is_zero_object(object))
            }
            _ => Err(invalid(
                "the category and module use different coefficient fields",
            )),
        }
    }

    fn is_identical(&self, right: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleCategoryKind,
            |left, right| left == right,
            false
        )
    }

    fn __repr__(&self) -> String {
        self.name()
    }
}

#[pymethods]
impl PyModuleMorphismList {
    fn value(&self, degree: i64) -> PyResult<PyModuleHomomorphism> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleMorphismListKind => ModuleHomomorphismKind,
            |list| list.value(&degree).map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    #[getter]
    fn highest_known_position(&self) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, ModuleMorphismListKind, |list| Ok(list
            .highest_known_position()
            .map_err(invalid)?
            .map(extended_integer_text)))
    }

    #[getter]
    fn lowest_known_position(&self) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, ModuleMorphismListKind, |list| Ok(list
            .lowest_known_position()
            .map_err(invalid)?
            .map(extended_integer_text)))
    }
}
