//! Python wrappers for QPA's module category and lazy morphism lists.

use std::sync::Arc;

use numerica::domains::{
    finite_field::{Z2, Zp},
    integer::Integer,
    rational::RationalField,
};
use pyo3::{Py, PyAny, PyResult, Python, prelude::*, types::PyTuple};
use qpa_core::{InfiniteList, PathAlgebraModuleHomomorphism, QuiverAlgebra};
use qpa_gbnp::{ComplexError, QpaCategory, RightModuleCategory};

#[allow(clippy::wildcard_imports)]
use super::*;
use crate::field::wrap_python_integer;

enum PythonMorphismHalf {
    Zero,
    Repeat(Vec<ModuleHomomorphismKind>),
    Next {
        callback: Py<PyAny>,
        initial: ModuleHomomorphismKind,
        repeatify: bool,
    },
    Position {
        callback: Py<PyAny>,
        store_values: bool,
    },
}

enum PythonCombinedHalf {
    Zero,
    Repeat(Vec<Vec<ModuleHomomorphismKind>>),
    Next {
        callback: Py<PyAny>,
        initial: Vec<ModuleHomomorphismKind>,
        repeatify: bool,
    },
    Position {
        callback: Py<PyAny>,
        store_values: bool,
    },
}

fn specification_items<'py>(value: &Bound<'py, PyAny>) -> PyResult<Vec<Bound<'py, PyAny>>> {
    value.try_iter()?.collect()
}

fn clone_one_homomorphism(value: &Bound<'_, PyAny>) -> PyResult<ModuleHomomorphismKind> {
    Ok(value
        .extract::<PyRef<'_, PyModuleHomomorphism>>()?
        .kind
        .clone())
}

fn parse_morphism_half(py: Python<'_>, value: Option<Py<PyAny>>) -> PyResult<PythonMorphismHalf> {
    let Some(value) = value else {
        return Ok(PythonMorphismHalf::Zero);
    };
    let value = value.bind(py);
    if value.extract::<String>().is_ok_and(|value| value == "zero") {
        return Ok(PythonMorphismHalf::Zero);
    }
    if let Ok(values) = value.extract::<Vec<Py<PyModuleHomomorphism>>>() {
        return Ok(PythonMorphismHalf::Repeat(clone_homomorphism_kinds(
            py, values,
        )));
    }
    let items = specification_items(value)?;
    let Some(kind) = items.first() else {
        return Err(invalid("a complex half specification cannot be empty"));
    };
    let kind = kind.extract::<String>()?;
    match kind.as_str() {
        "repeat" if items.len() == 2 => {
            let values = items[1].extract::<Vec<Py<PyModuleHomomorphism>>>()?;
            Ok(PythonMorphismHalf::Repeat(clone_homomorphism_kinds(
                py, values,
            )))
        }
        "next" | "next/repeat" if items.len() == 3 => {
            if !items[1].is_callable() {
                return Err(invalid("the next-morphism function must be callable"));
            }
            Ok(PythonMorphismHalf::Next {
                callback: items[1].clone().unbind(),
                initial: clone_one_homomorphism(&items[2])?,
                repeatify: kind == "next/repeat",
            })
        }
        "pos" if (2..=3).contains(&items.len()) => {
            if !items[1].is_callable() {
                return Err(invalid("the position-morphism function must be callable"));
            }
            Ok(PythonMorphismHalf::Position {
                callback: items[1].clone().unbind(),
                store_values: if items.len() == 3 {
                    items[2].extract::<bool>()?
                } else {
                    true
                },
            })
        }
        _ => Err(invalid("invalid QPA complex half specification")),
    }
}

fn parse_combined_half(py: Python<'_>, value: Option<Py<PyAny>>) -> PyResult<PythonCombinedHalf> {
    let Some(value) = value else {
        return Ok(PythonCombinedHalf::Zero);
    };
    let value = value.bind(py);
    if value.extract::<String>().is_ok_and(|value| value == "zero") {
        return Ok(PythonCombinedHalf::Zero);
    }
    if let Ok(rows) = value.extract::<Vec<Vec<Py<PyModuleHomomorphism>>>>() {
        return Ok(PythonCombinedHalf::Repeat(clone_homomorphism_rows(
            py, rows,
        )));
    }
    let items = specification_items(value)?;
    let Some(kind) = items.first() else {
        return Err(invalid(
            "a combined complex half specification cannot be empty",
        ));
    };
    let kind = kind.extract::<String>()?;
    match kind.as_str() {
        "repeat" if items.len() == 2 => {
            let rows = items[1].extract::<Vec<Vec<Py<PyModuleHomomorphism>>>>()?;
            Ok(PythonCombinedHalf::Repeat(clone_homomorphism_rows(
                py, rows,
            )))
        }
        "next" | "next/repeat" if items.len() == 3 => {
            if !items[1].is_callable() {
                return Err(invalid("the next-tuple function must be callable"));
            }
            let initial = items[2].extract::<Vec<Py<PyModuleHomomorphism>>>()?;
            Ok(PythonCombinedHalf::Next {
                callback: items[1].clone().unbind(),
                initial: clone_homomorphism_kinds(py, initial),
                repeatify: kind == "next/repeat",
            })
        }
        "pos" if (2..=3).contains(&items.len()) => {
            if !items[1].is_callable() {
                return Err(invalid("the position-tuple function must be callable"));
            }
            Ok(PythonCombinedHalf::Position {
                callback: items[1].clone().unbind(),
                store_values: if items.len() == 3 {
                    items[2].extract::<bool>()?
                } else {
                    true
                },
            })
        }
        _ => Err(invalid("invalid QPA combined complex half specification")),
    }
}

fn callback_error(error: PyErr) -> ComplexError {
    ComplexError::CategoryCallback(error.to_string())
}

fn extract_callback_homomorphism<F: runtime::PythonField>(
    value: &Bound<'_, PyAny>,
) -> PyResult<PathAlgebraModuleHomomorphism<F>>
where
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>,
{
    runtime::Extract::extract(clone_one_homomorphism(value)?)
}

fn extract_callback_row<F: runtime::PythonField>(
    py: Python<'_>,
    value: &Bound<'_, PyAny>,
) -> PyResult<Vec<PathAlgebraModuleHomomorphism<F>>>
where
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>,
{
    runtime::extract_all(clone_homomorphism_kinds(
        py,
        value.extract::<Vec<Py<PyModuleHomomorphism>>>()?,
    ))
}

fn callback_homomorphism<F: runtime::PythonField>(
    py: Python<'_>,
    value: PathAlgebraModuleHomomorphism<F>,
) -> PyResult<Py<PyModuleHomomorphism>>
where
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    Py::new(
        py,
        PyModuleHomomorphism {
            kind: runtime::Wrap::wrap(value),
        },
    )
}

fn callback_homomorphism_row<F: runtime::PythonField>(
    py: Python<'_>,
    values: &[PathAlgebraModuleHomomorphism<F>],
) -> PyResult<Vec<Py<PyModuleHomomorphism>>>
where
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    values
        .iter()
        .cloned()
        .map(|value| callback_homomorphism(py, value))
        .collect()
}

fn clone_homomorphism_kinds(
    py: Python<'_>,
    morphisms: Vec<Py<PyModuleHomomorphism>>,
) -> Vec<ModuleHomomorphismKind> {
    morphisms
        .into_iter()
        .map(|morphism| morphism.borrow(py).kind.clone())
        .collect()
}

fn clone_complex_kinds(
    py: Python<'_>,
    complexes: Vec<Py<PyModuleComplex>>,
) -> Vec<ModuleComplexKind> {
    complexes
        .into_iter()
        .map(|complex| complex.borrow(py).kind.clone())
        .collect()
}

fn clone_homomorphism_rows(
    py: Python<'_>,
    rows: Vec<Vec<Py<PyModuleHomomorphism>>>,
) -> Vec<Vec<ModuleHomomorphismKind>> {
    rows.into_iter()
        .map(|row| clone_homomorphism_kinds(py, row))
        .collect()
}

fn complex_half<F: runtime::PythonField + Send + Sync + 'static>(
    specification: PythonMorphismHalf,
) -> PyResult<qpa_gbnp::ComplexHalf<RightModuleCategory<F>>>
where
    F::Element: Send + Sync,
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>
        + runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
    ModuleComplexKind: runtime::Wrap<qpa_gbnp::ModuleComplex<F>>,
{
    Ok(match specification {
        PythonMorphismHalf::Zero => qpa_gbnp::ComplexHalf::Zero,
        PythonMorphismHalf::Repeat(values) => {
            qpa_gbnp::ComplexHalf::Repeat(runtime::extract_all(values)?)
        }
        PythonMorphismHalf::Next {
            callback,
            initial,
            repeatify,
        } => qpa_gbnp::ComplexHalf::Next {
            function: Arc::new(move |previous| {
                Python::attach(|py| {
                    let previous = callback_homomorphism(py, previous.clone())?;
                    let value = callback.call1(py, (previous,))?;
                    extract_callback_homomorphism(value.bind(py))
                })
                .map_err(callback_error)
            }),
            initial_value: runtime::Extract::extract(initial)?,
            repeatify,
        },
        PythonMorphismHalf::Position {
            callback,
            store_values,
        } => qpa_gbnp::ComplexHalf::Position {
            function: Arc::new(move |complex, degree| {
                Python::attach(|py| {
                    let complex = Py::new(
                        py,
                        PyModuleComplex {
                            kind: runtime::Wrap::wrap(complex.clone()),
                        },
                    )?;
                    let degree = wrap_python_integer(py, degree)?;
                    let value = callback.call1(py, (complex, degree))?;
                    extract_callback_homomorphism(value.bind(py))
                })
                .map_err(callback_error)
            }),
            store_values,
        },
    })
}

fn chain_map_half<F: runtime::PythonField + Send + Sync + 'static>(
    specification: PythonMorphismHalf,
) -> PyResult<qpa_gbnp::ChainMapHalf<RightModuleCategory<F>>>
where
    F::Element: Send + Sync,
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>
        + runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
    ModuleChainMapKind: runtime::Wrap<qpa_gbnp::ModuleChainMap<F>>,
{
    Ok(match specification {
        PythonMorphismHalf::Zero => qpa_gbnp::ChainMapHalf::Zero,
        PythonMorphismHalf::Repeat(values) => {
            qpa_gbnp::ChainMapHalf::Repeat(runtime::extract_all(values)?)
        }
        PythonMorphismHalf::Next {
            callback,
            initial,
            repeatify,
        } => qpa_gbnp::ChainMapHalf::Next {
            function: Arc::new(move |previous| {
                Python::attach(|py| {
                    let previous = callback_homomorphism(py, previous.clone())?;
                    let value = callback.call1(py, (previous,))?;
                    extract_callback_homomorphism(value.bind(py))
                })
                .map_err(callback_error)
            }),
            initial_value: runtime::Extract::extract(initial)?,
            repeatify,
        },
        PythonMorphismHalf::Position {
            callback,
            store_values,
        } => qpa_gbnp::ChainMapHalf::Position {
            function: Arc::new(move |map, degree| {
                Python::attach(|py| {
                    let map = Py::new(
                        py,
                        PyModuleChainMap {
                            kind: runtime::Wrap::wrap(map.clone()),
                        },
                    )?;
                    let degree = wrap_python_integer(py, degree)?;
                    let value = callback.call1(py, (map, degree))?;
                    extract_callback_homomorphism(value.bind(py))
                })
                .map_err(callback_error)
            }),
            store_values,
        },
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
    positive: PythonMorphismHalf,
    negative: PythonMorphismHalf,
) -> PyResult<PyModuleComplex>
where
    F: runtime::PythonField + Send + Sync + 'static,
    F::Element: Send + Sync,
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>
        + runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
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
    positive: Option<Py<PyAny>>,
    negative: Option<Py<PyAny>>,
) -> PyResult<PyModuleComplex> {
    let middle = clone_homomorphism_kinds(py, middle);
    let positive = parse_morphism_half(py, positive)?;
    let negative = parse_morphism_half(py, negative)?;
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
    positive: PythonMorphismHalf,
    negative: PythonMorphismHalf,
) -> PyResult<PyModuleChainMap>
where
    F: runtime::PythonField + Send + Sync + 'static,
    F::Element: Send + Sync,
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>
        + runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
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
    positive: Option<Py<PyAny>>,
    negative: Option<Py<PyAny>>,
) -> PyResult<PyModuleChainMap> {
    let middle = clone_homomorphism_kinds(py, middle);
    let positive = parse_morphism_half(py, positive)?;
    let negative = parse_morphism_half(py, negative)?;
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

fn combined_half<F>(
    specification: PythonCombinedHalf,
) -> PyResult<qpa_gbnp::ComplexAndChainMapsHalf<RightModuleCategory<F>>>
where
    F: runtime::PythonField + Send + Sync + 'static,
    F::Element: Send + Sync,
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>
        + runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
    ModuleComplexKind: runtime::Wrap<qpa_gbnp::ModuleComplex<F>>,
    ModuleChainMapKind: runtime::Wrap<qpa_gbnp::ModuleChainMap<F>>,
{
    Ok(match specification {
        PythonCombinedHalf::Zero => qpa_gbnp::ComplexAndChainMapsHalf::Zero,
        PythonCombinedHalf::Repeat(values) => qpa_gbnp::ComplexAndChainMapsHalf::Repeat(
            values
                .into_iter()
                .map(runtime::extract_all)
                .collect::<PyResult<Vec<_>>>()?,
        ),
        PythonCombinedHalf::Next {
            callback,
            initial,
            repeatify,
        } => qpa_gbnp::ComplexAndChainMapsHalf::Next {
            function: Arc::new(move |previous| {
                Python::attach(|py| {
                    let previous = callback_homomorphism_row(py, previous)?;
                    let value = callback.call1(py, (previous,))?;
                    extract_callback_row(py, value.bind(py))
                })
                .map_err(callback_error)
            }),
            initial_value: runtime::extract_all(initial)?,
            repeatify,
        },
        PythonCombinedHalf::Position {
            callback,
            store_values,
        } => qpa_gbnp::ComplexAndChainMapsHalf::Position {
            function: Arc::new(move |complex, incoming, outgoing, degree| {
                Python::attach(|py| {
                    let mut arguments = Vec::with_capacity(2 + incoming.len() + outgoing.len());
                    arguments.push(
                        Py::new(
                            py,
                            PyModuleComplex {
                                kind: runtime::Wrap::wrap(complex.clone()),
                            },
                        )?
                        .into_any(),
                    );
                    for map in incoming.iter().chain(outgoing) {
                        arguments.push(
                            Py::new(
                                py,
                                PyModuleChainMap {
                                    kind: runtime::Wrap::wrap(map.clone()),
                                },
                            )?
                            .into_any(),
                        );
                    }
                    arguments.push(wrap_python_integer(py, degree)?);
                    let arguments = PyTuple::new(py, arguments)?;
                    let value = callback.bind(py).call1(arguments)?;
                    extract_callback_row(py, &value)
                })
                .map_err(callback_error)
            }),
            store_values,
        },
    })
}

fn complex_and_chain_maps_for<F>(
    source_complexes: &[ModuleComplexKind],
    range_complexes: &[ModuleComplexKind],
    base_degree: i64,
    middle: &[Vec<ModuleHomomorphismKind>],
    positive: PythonCombinedHalf,
    negative: PythonCombinedHalf,
) -> PyResult<(
    PyModuleComplex,
    Vec<PyModuleChainMap>,
    Vec<PyModuleChainMap>,
)>
where
    F: runtime::PythonField + Send + Sync + 'static,
    F::Element: Send + Sync,
    ModuleComplexKind:
        runtime::Extract<qpa_gbnp::ModuleComplex<F>> + runtime::Wrap<qpa_gbnp::ModuleComplex<F>>,
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>
        + runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
    ModuleChainMapKind: runtime::Wrap<qpa_gbnp::ModuleChainMap<F>>,
{
    let source_complexes = source_complexes
        .iter()
        .cloned()
        .map(runtime::Extract::extract)
        .collect::<PyResult<Vec<qpa_gbnp::ModuleComplex<F>>>>()?;
    let range_complexes = range_complexes
        .iter()
        .cloned()
        .map(runtime::Extract::extract)
        .collect::<PyResult<Vec<qpa_gbnp::ModuleComplex<F>>>>()?;
    let middle = middle
        .iter()
        .cloned()
        .map(runtime::extract_all)
        .collect::<PyResult<Vec<Vec<PathAlgebraModuleHomomorphism<F>>>>>()?;
    let result = qpa_gbnp::complex_and_chain_maps(
        &source_complexes,
        &range_complexes,
        Integer::from(base_degree),
        middle,
        combined_half::<F>(positive)?,
        combined_half::<F>(negative)?,
    )
    .map_err(invalid)?;
    Ok((
        PyModuleComplex {
            kind: runtime::Wrap::wrap(result.complex),
        },
        runtime::wrap_all(result.incoming)
            .into_iter()
            .map(|kind| PyModuleChainMap { kind })
            .collect(),
        runtime::wrap_all(result.outgoing)
            .into_iter()
            .map(|kind| PyModuleChainMap { kind })
            .collect(),
    ))
}

pub(super) fn complex_and_chain_maps(
    py: Python<'_>,
    source_complexes: Vec<Py<PyModuleComplex>>,
    range_complexes: Vec<Py<PyModuleComplex>>,
    base_degree: i64,
    middle: Vec<Vec<Py<PyModuleHomomorphism>>>,
    positive: Option<Py<PyAny>>,
    negative: Option<Py<PyAny>>,
) -> PyResult<(
    PyModuleComplex,
    Vec<PyModuleChainMap>,
    Vec<PyModuleChainMap>,
)> {
    let source_complexes = clone_complex_kinds(py, source_complexes);
    let range_complexes = clone_complex_kinds(py, range_complexes);
    let middle = clone_homomorphism_rows(py, middle);
    let positive = parse_combined_half(py, positive)?;
    let negative = parse_combined_half(py, negative)?;
    let first = source_complexes
        .first()
        .or_else(|| range_complexes.first())
        .ok_or_else(|| invalid(qpa_gbnp::ComplexError::CombinedComplexListEmpty))?;
    match first {
        ModuleComplexKind::Rational(_) => complex_and_chain_maps_for::<RationalField>(
            &source_complexes,
            &range_complexes,
            base_degree,
            &middle,
            positive,
            negative,
        ),
        ModuleComplexKind::Binary(_) => complex_and_chain_maps_for::<Z2>(
            &source_complexes,
            &range_complexes,
            base_degree,
            &middle,
            positive,
            negative,
        ),
        ModuleComplexKind::Prime(_) => complex_and_chain_maps_for::<Zp>(
            &source_complexes,
            &range_complexes,
            base_degree,
            &middle,
            positive,
            negative,
        ),
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
        format!("<cat: {}>", self.name())
    }

    fn print_obj(&self) -> String {
        self.__repr__()
    }

    fn view(&self) -> String {
        self.__repr__()
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
