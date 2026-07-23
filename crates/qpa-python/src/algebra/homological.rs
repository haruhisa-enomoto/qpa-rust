//! Lazy projective and injective resolutions.

use pyo3::{PyAny, pyfunction, pymethods, types::PyString};

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

fn simple_endo_projective_dimension(
    py: Python<'_>,
    status: &qpa_gbnp::ProjectiveResolutionOfSimpleModuleOverEndoStatus,
) -> PyResult<Py<PyAny>> {
    match status {
        qpa_gbnp::ProjectiveResolutionOfSimpleModuleOverEndoStatus::Finite(dimension) => Ok(py
            .import("builtins")?
            .getattr("int")?
            .call1((dimension.to_string(),))?
            .unbind()),
        qpa_gbnp::ProjectiveResolutionOfSimpleModuleOverEndoStatus::GreaterThanModuleIndex(
            index,
        ) => Ok(PyString::new(py, &format!("projdim > {index}"))
            .into_any()
            .unbind()),
    }
}

#[pyfunction(name = "complex_and_chain_maps")]
#[pyo3(signature = (source_complexes, range_complexes, base_degree, middle, positive=None, negative=None))]
pub(crate) fn py_complex_and_chain_maps(
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
    complex_foundation::complex_and_chain_maps(
        py,
        source_complexes,
        range_complexes,
        base_degree,
        middle,
        positive,
        negative,
    )
}

#[pymethods]
impl PyProjectiveResolution {
    /// Compute or retrieve the differential at the requested degree.
    fn differential(&self, degree: i64) -> PyResult<PyModuleHomomorphism> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ProjectiveResolutionKind => ModuleHomomorphismKind,
            |resolution| resolution.differential(&degree).map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    /// Compute or retrieve the resolution object at the requested degree.
    fn object(&self, degree: i64) -> PyResult<PyMatrixModule> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ProjectiveResolutionKind => MatrixModuleKind,
            |resolution| resolution.object(&degree).map_err(invalid)?
        );
        Ok(PyMatrixModule { kind })
    }
}

#[pymethods]
impl PyInjectiveResolution {
    fn differential(&self, degree: i64) -> PyResult<PyModuleHomomorphism> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            InjectiveResolutionKind => ModuleHomomorphismKind,
            |resolution| resolution.differential(&degree).map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    fn object(&self, degree: i64) -> PyResult<PyMatrixModule> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            InjectiveResolutionKind => MatrixModuleKind,
            |resolution| resolution.object(&degree).map_err(invalid)?
        );
        Ok(PyMatrixModule { kind })
    }
}

#[pymethods]
impl PySimpleEndoProjectiveResolutionResult {
    #[getter]
    fn projective_dimension(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        dispatch_ref!(&self.kind, SimpleEndoProjectiveResolutionKind, |result| {
            simple_endo_projective_dimension(py, result.projective_dimension())
        })
    }

    #[getter]
    fn syzygies(&self) -> Vec<PyMatrixModule> {
        dispatch_ref!(&self.kind, SimpleEndoProjectiveResolutionKind, |result| {
            result
                .syzygies()
                .iter()
                .cloned()
                .map(module::wrap_module)
                .collect()
        })
    }
}

fn extended_integer_text(value: qpa_core::ExtendedInteger) -> String {
    match value {
        qpa_core::ExtendedInteger::NegativeInfinity => "-infinity".into(),
        qpa_core::ExtendedInteger::Finite(value) => value.to_string(),
        qpa_core::ExtendedInteger::PositiveInfinity => "infinity".into(),
    }
}

fn integer_text(values: Vec<Integer>) -> Vec<String> {
    values.into_iter().map(|value| value.to_string()).collect()
}

fn wrap_mapping_cone<F>(
    result: qpa_gbnp::ComplexAndChainMaps<qpa_gbnp::RightModuleCategory<F>>,
) -> PyResult<(PyModuleComplex, PyModuleChainMap, PyModuleChainMap)>
where
    F: runtime::PythonField + Send + Sync + 'static,
    F::Element: Send + Sync,
    ModuleComplexKind: runtime::Wrap<qpa_gbnp::ModuleComplex<F>>,
    ModuleChainMapKind: runtime::Wrap<qpa_gbnp::ModuleChainMap<F>>,
{
    let inclusion = result
        .incoming
        .into_iter()
        .next()
        .ok_or_else(|| invalid("QPA mapping cone did not return its inclusion chain map"))?;
    let projection = result
        .outgoing
        .into_iter()
        .next()
        .ok_or_else(|| invalid("QPA mapping cone did not return its projection chain map"))?;
    Ok((
        PyModuleComplex {
            kind: runtime::Wrap::wrap(result.complex),
        },
        PyModuleChainMap {
            kind: runtime::Wrap::wrap(inclusion),
        },
        PyModuleChainMap {
            kind: runtime::Wrap::wrap(projection),
        },
    ))
}

fn for_every_degree_for<F>(
    complex: &qpa_gbnp::ModuleComplex<F>,
    function: &Bound<'_, PyAny>,
) -> PyResult<Option<bool>>
where
    F: runtime::PythonField + Send + Sync + 'static,
    F::Element: Send + Sync,
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    complex
        .for_every_degree(|left, right| {
            function
                .call1((
                    PyModuleHomomorphism {
                        kind: runtime::Wrap::wrap(left.clone()),
                    },
                    PyModuleHomomorphism {
                        kind: runtime::Wrap::wrap(right.clone()),
                    },
                ))
                .and_then(|result| result.extract::<bool>())
                .map_err(|error| qpa_gbnp::ComplexError::CategoryCallback(error.to_string()))
        })
        .map_err(invalid)
}

fn complex_print_obj_for<F>(complex: &qpa_gbnp::ModuleComplex<F>) -> PyResult<String>
where
    F: runtime::PythonField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    if qpa_gbnp::is_projective_complex(complex).map_err(invalid)? {
        qpa_gbnp::qpa_projective_complex_string(complex).map_err(invalid)
    } else if qpa_gbnp::is_injective_complex(complex).map_err(invalid)? {
        qpa_gbnp::qpa_injective_complex_string(complex).map_err(invalid)
    } else {
        complex.qpa_string().map_err(invalid)
    }
}

#[pymethods]
impl PyModuleComplex {
    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn category(&self) -> PyModuleCategory {
        complex_foundation::category_of_complex(self)
    }

    #[getter]
    fn differentials(&self) -> PyModuleMorphismList {
        complex_foundation::differentials_of_complex(self)
    }

    #[staticmethod]
    fn from_differential_list(
        category: &PyModuleCategory,
        differentials: &PyModuleMorphismList,
    ) -> PyResult<Self> {
        complex_foundation::complex_by_differential_list(category, differentials)
    }

    #[staticmethod]
    #[pyo3(signature = (category, base_degree, middle, positive=None, negative=None))]
    fn from_parts(
        category: &PyModuleCategory,
        py: Python<'_>,
        base_degree: i64,
        middle: Vec<Py<PyModuleHomomorphism>>,
        positive: Option<Py<PyAny>>,
        negative: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        complex_foundation::complex_from_parts(
            category,
            py,
            base_degree,
            middle,
            positive,
            negative,
        )
    }

    fn for_every_degree(&self, function: &Bound<'_, PyAny>) -> PyResult<Option<bool>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| {
            for_every_degree_for(complex, function)
        })
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| complex
            .is_zero_complex())
    }

    fn differential(&self, degree: i64) -> PyResult<PyModuleHomomorphism> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleHomomorphismKind,
            |complex| complex.differential(&degree).map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    fn object(&self, degree: i64) -> PyResult<PyMatrixModule> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => MatrixModuleKind,
            |complex| complex.object(&degree).map_err(invalid)?
        );
        Ok(PyMatrixModule { kind })
    }

    fn cycles(&self, degree: i64) -> PyResult<PyMatrixModule> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => MatrixModuleKind,
            |complex| complex.cycles(&degree).map_err(invalid)?
        );
        Ok(PyMatrixModule { kind })
    }

    fn boundaries(&self, degree: i64) -> PyResult<PyMatrixModule> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => MatrixModuleKind,
            |complex| complex.boundaries(&degree).map_err(invalid)?
        );
        Ok(PyMatrixModule { kind })
    }

    fn homology(&self, degree: i64) -> PyResult<PyMatrixModule> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => MatrixModuleKind,
            |complex| complex.homology(&degree).map_err(invalid)?
        );
        Ok(PyMatrixModule { kind })
    }

    #[getter]
    fn highest_known_degree(&self) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| Ok(complex
            .highest_known_degree()
            .map_err(invalid)?
            .map(extended_integer_text)))
    }

    #[getter]
    fn lowest_known_degree(&self) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| Ok(complex
            .lowest_known_degree()
            .map_err(invalid)?
            .map(extended_integer_text)))
    }

    #[getter]
    fn is_positive_repeating(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| complex
            .is_positive_repeating()
            .map_err(invalid))
    }

    #[getter]
    fn is_negative_repeating(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| complex
            .is_negative_repeating()
            .map_err(invalid))
    }

    fn positive_repeat_degrees(&self) -> PyResult<Option<Vec<String>>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| Ok(complex
            .positive_repeat_degrees()
            .map_err(invalid)?
            .map(integer_text)))
    }

    fn negative_repeat_degrees(&self) -> PyResult<Option<Vec<String>>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| Ok(complex
            .negative_repeat_degrees()
            .map_err(invalid)?
            .map(integer_text)))
    }

    #[getter]
    fn upper_bound(&self) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| Ok(complex
            .upper_bound()
            .map_err(invalid)?
            .map(extended_integer_text)))
    }

    #[getter]
    fn lower_bound(&self) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| Ok(complex
            .lower_bound()
            .map_err(invalid)?
            .map(extended_integer_text)))
    }

    #[getter]
    fn is_finite(&self) -> PyResult<Option<bool>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| complex
            .is_finite_complex()
            .map_err(invalid))
    }

    #[getter]
    fn length(&self) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| Ok(complex
            .length()
            .map_err(invalid)?
            .map(extended_integer_text)))
    }

    fn is_exact_in_degree(&self, degree: i64) -> PyResult<bool> {
        let degree = Integer::from(degree);
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| complex
            .is_exact_in_degree(&degree)
            .map_err(invalid))
    }

    #[getter]
    fn is_exact_sequence(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| complex
            .is_exact_sequence()
            .map_err(invalid))
    }

    #[getter]
    fn is_short_exact_sequence(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| complex
            .is_short_exact_sequence()
            .map_err(invalid))
    }

    #[getter]
    fn is_projective_complex(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| {
            qpa_gbnp::is_projective_complex(complex).map_err(invalid)
        })
    }

    #[getter]
    fn is_injective_complex(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| {
            qpa_gbnp::is_injective_complex(complex).map_err(invalid)
        })
    }

    fn description_in_degree(&self, degree: i64, projective: bool) -> PyResult<Vec<usize>> {
        let degree = Integer::from(degree);
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| {
            qpa_gbnp::description_of_proj_or_inj_complex_in_degree(complex, &degree, projective)
                .map_err(invalid)
        })
    }

    fn projective_description_in_degree(&self, degree: i64) -> PyResult<Vec<usize>> {
        let degree = Integer::from(degree);
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| {
            qpa_gbnp::description_of_proj_complex_in_degree(complex, &degree).map_err(invalid)
        })
    }

    fn injective_description_in_degree(&self, degree: i64) -> PyResult<Vec<usize>> {
        let degree = Integer::from(degree);
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| {
            qpa_gbnp::description_of_inj_complex_in_degree(complex, &degree).map_err(invalid)
        })
    }

    fn finite_description(&self, projective: bool) -> PyResult<Vec<Vec<usize>>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| {
            qpa_gbnp::description_of_finite_proj_or_inj_complex(complex, projective)
                .map_err(invalid)
        })
    }

    fn finite_projective_description(&self) -> PyResult<Vec<Vec<usize>>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| {
            qpa_gbnp::description_of_finite_proj_complex(complex).map_err(invalid)
        })
    }

    fn finite_injective_description(&self) -> PyResult<Vec<Vec<usize>>> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| {
            qpa_gbnp::description_of_finite_inj_complex(complex).map_err(invalid)
        })
    }

    fn shift(&self, amount: i64) -> PyResult<Self> {
        let amount = Integer::from(amount);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.shift(&amount).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn shift_unsigned(&self, amount: i64) -> PyResult<Self> {
        let amount = Integer::from(amount);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.shift_unsigned(&amount).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn good_truncation_below(&self, degree: i64) -> PyResult<Self> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.good_truncation_below(&degree).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn good_truncation_above(&self, degree: i64) -> PyResult<Self> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.good_truncation_above(&degree).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn brutal_truncation_below(&self, degree: i64) -> PyResult<Self> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.brutal_truncation_below(&degree).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn brutal_truncation_above(&self, degree: i64) -> PyResult<Self> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.brutal_truncation_above(&degree).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn syzygy_truncation(&self, degree: i64) -> PyResult<Self> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.syzygy_truncation(&degree).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn cosyzygy_truncation(&self, degree: i64) -> PyResult<Self> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.cosyzygy_truncation(&degree).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn good_truncation(&self, upper: i64, lower: i64) -> PyResult<Self> {
        let upper = Integer::from(upper);
        let lower = Integer::from(lower);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.good_truncation(&upper, &lower).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn brutal_truncation(&self, upper: i64, lower: i64) -> PyResult<Self> {
        let upper = Integer::from(upper);
        let lower = Integer::from(lower);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.brutal_truncation(&upper, &lower).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn syzygy_cosyzygy_truncation(&self, upper: i64, lower: i64) -> PyResult<Self> {
        let upper = Integer::from(upper);
        let lower = Integer::from(lower);
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex
                .syzygy_cosyzygy_truncation(&upper, &lower)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn cut_above(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.cut_complex_above().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn cut_below(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| complex.cut_complex_below().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn yoneda_product(&self, second: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &second.kind,
            ModuleComplexKind => ModuleComplexKind,
            |first, second| qpa_gbnp::yoneda_product(first, second).map_err(invalid)?,
            invalid("complexes belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn projective_resolution(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| qpa_gbnp::projective_resolution_of_complex(complex).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn projective_to_injective_finite(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| qpa_gbnp::projective_to_injective_finite_complex(complex)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn projective_to_injective(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| qpa_gbnp::projective_to_injective_complex(complex).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn tau(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleComplexKind => ModuleComplexKind,
            |complex| qpa_gbnp::tau_of_complex(complex).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn finite_chain_map_to(
        &self,
        py: Python<'_>,
        range: &Self,
        base_degree: i64,
        morphisms: Vec<Py<PyModuleHomomorphism>>,
    ) -> PyResult<PyModuleChainMap> {
        let morphisms = clone_module_homomorphisms(py, morphisms);
        let kind = match (&self.kind, &range.kind) {
            (ModuleComplexKind::Rational(source), ModuleComplexKind::Rational(range)) => {
                ModuleChainMapKind::Rational(
                    qpa_gbnp::finite_chain_map(
                        source.clone(),
                        range.clone(),
                        Integer::from(base_degree),
                        runtime::extract_all(morphisms)?,
                    )
                    .map_err(invalid)?,
                )
            }
            (ModuleComplexKind::Binary(source), ModuleComplexKind::Binary(range)) => {
                ModuleChainMapKind::Binary(
                    qpa_gbnp::finite_chain_map(
                        source.clone(),
                        range.clone(),
                        Integer::from(base_degree),
                        runtime::extract_all(morphisms)?,
                    )
                    .map_err(invalid)?,
                )
            }
            (ModuleComplexKind::Prime(source), ModuleComplexKind::Prime(range)) => {
                ModuleChainMapKind::Prime(
                    qpa_gbnp::finite_chain_map(
                        source.clone(),
                        range.clone(),
                        Integer::from(base_degree),
                        runtime::extract_all(morphisms)?,
                    )
                    .map_err(invalid)?,
                )
            }
            (ModuleComplexKind::Extension(source), ModuleComplexKind::Extension(range)) => {
                ModuleChainMapKind::Extension(
                    qpa_gbnp::finite_chain_map(
                        source.clone(),
                        range.clone(),
                        Integer::from(base_degree),
                        runtime::extract_all(morphisms)?,
                    )
                    .map_err(invalid)?,
                )
            }
            _ => {
                return Err(invalid(
                    "source and range complexes belong to different coefficient fields",
                ));
            }
        };
        Ok(PyModuleChainMap { kind })
    }

    fn zero_chain_map_to(&self, range: &Self) -> PyResult<PyModuleChainMap> {
        let kind = map_pair_ref!(
            &self.kind,
            &range.kind,
            ModuleComplexKind => ModuleChainMapKind,
            |source, range| qpa_gbnp::zero_chain_map(source.clone(), range.clone())
                .map_err(invalid)?,
            invalid("source and range complexes belong to different coefficient fields")
        )?;
        Ok(PyModuleChainMap { kind })
    }

    fn is_identical(&self, right: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleComplexKind,
            |left, right| left == right,
            false
        )
    }

    fn __eq__(&self, right: &Self) -> bool {
        self.is_identical(right)
    }

    fn print_obj(&self) -> PyResult<String> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| {
            complex_print_obj_for(complex)
        })
    }

    fn __repr__(&self) -> PyResult<String> {
        self.print_obj()
    }
}

#[pymethods]
impl PyModuleChainMap {
    #[staticmethod]
    #[pyo3(signature = (source, range, base_degree, middle, positive=None, negative=None))]
    fn from_parts(
        source: &PyModuleComplex,
        range: &PyModuleComplex,
        py: Python<'_>,
        base_degree: i64,
        middle: Vec<Py<PyModuleHomomorphism>>,
        positive: Option<Py<PyAny>>,
        negative: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        complex_foundation::chain_map_from_parts(
            source,
            range,
            py,
            base_degree,
            middle,
            positive,
            negative,
        )
    }

    #[getter]
    fn morphisms(&self) -> PyModuleMorphismList {
        complex_foundation::morphisms_of_chain_map(self)
    }

    #[getter]
    fn source(&self) -> PyModuleComplex {
        let kind = map_ref!(
            &self.kind,
            ModuleChainMapKind => ModuleComplexKind,
            |map| map.source().clone()
        );
        PyModuleComplex { kind }
    }

    #[getter]
    fn range(&self) -> PyModuleComplex {
        let kind = map_ref!(
            &self.kind,
            ModuleChainMapKind => ModuleComplexKind,
            |map| map.range().clone()
        );
        PyModuleComplex { kind }
    }

    fn morphism(&self, degree: i64) -> PyResult<PyModuleHomomorphism> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ModuleChainMapKind => ModuleHomomorphismKind,
            |map| map.morphism(&degree).map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    #[getter]
    fn highest_known_degree(&self) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, ModuleChainMapKind, |map| Ok(map
            .highest_known_degree()
            .map_err(invalid)?
            .map(extended_integer_text)))
    }

    #[getter]
    fn lowest_known_degree(&self) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, ModuleChainMapKind, |map| Ok(map
            .lowest_known_degree()
            .map_err(invalid)?
            .map(extended_integer_text)))
    }

    fn mapping_cone(&self) -> PyResult<(PyModuleComplex, PyModuleChainMap, PyModuleChainMap)> {
        dispatch_ref!(&self.kind, ModuleChainMapKind, |map| wrap_mapping_cone(
            qpa_gbnp::mapping_cone(map).map_err(invalid)?
        ))
    }

    fn is_identical(&self, right: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleChainMapKind,
            |left, right| left == right,
            false
        )
    }

    fn __eq__(&self, right: &Self) -> bool {
        self.is_identical(right)
    }

    fn print_obj(&self) -> &'static str {
        dispatch_ref!(&self.kind, ModuleChainMapKind, |map| map.qpa_string())
    }

    fn __repr__(&self) -> &'static str {
        self.print_obj()
    }
}

#[pymethods]
impl PyProjectiveIsomorphismComponent {
    #[getter]
    fn map(&self) -> PyModuleHomomorphism {
        let kind = map_ref!(
            &self.kind,
            ProjectiveIsomorphismComponentKind => ModuleHomomorphismKind,
            |component| component.map().clone()
        );
        PyModuleHomomorphism { kind }
    }

    #[getter]
    fn vertex(&self) -> usize {
        dispatch_ref!(
            &self.kind,
            ProjectiveIsomorphismComponentKind,
            |component| component.vertex()
        )
    }

    #[getter]
    fn source_position(&self) -> usize {
        dispatch_ref!(
            &self.kind,
            ProjectiveIsomorphismComponentKind,
            |component| component.source_position()
        )
    }

    #[getter]
    fn range_position(&self) -> usize {
        dispatch_ref!(
            &self.kind,
            ProjectiveIsomorphismComponentKind,
            |component| component.range_position()
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "<projective isomorphism at vertex {}, source summand {}, range summand {}>",
            self.vertex(),
            self.source_position(),
            self.range_position(),
        )
    }
}

#[pymethods]
impl PyExtOverAlgebraResult {
    #[getter]
    fn syzygy_inclusion(&self) -> PyModuleHomomorphism {
        let kind = map_ref!(
            &self.kind,
            ExtOverAlgebraKind => ModuleHomomorphismKind,
            |result| result.syzygy_inclusion().clone()
        );
        PyModuleHomomorphism { kind }
    }

    #[getter]
    fn basis(&self) -> Vec<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, ExtOverAlgebraKind, |result| result
            .basis()
            .iter()
            .cloned()
            .map(|homomorphism| PyModuleHomomorphism {
                kind: runtime::Wrap::wrap(homomorphism),
            })
            .collect())
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, ExtOverAlgebraKind, |result| result.dimension())
    }

    fn coefficients(&self, homomorphism: &PyModuleHomomorphism) -> PyResult<Option<Vec<PyScalar>>> {
        match (&self.kind, &homomorphism.kind) {
            (ExtOverAlgebraKind::Rational(result), ModuleHomomorphismKind::Rational(map)) => {
                let field = result
                    .syzygy_inclusion()
                    .source()
                    .right_acting_algebra()
                    .field();
                Ok(result.coefficients(map).map_err(invalid)?.map(|values| {
                    values
                        .iter()
                        .map(|value| field.wrap_scalar(value))
                        .collect()
                }))
            }
            (ExtOverAlgebraKind::Binary(result), ModuleHomomorphismKind::Binary(map)) => {
                let field = result
                    .syzygy_inclusion()
                    .source()
                    .right_acting_algebra()
                    .field();
                Ok(result.coefficients(map).map_err(invalid)?.map(|values| {
                    values
                        .iter()
                        .map(|value| field.wrap_scalar(value))
                        .collect()
                }))
            }
            (ExtOverAlgebraKind::Prime(result), ModuleHomomorphismKind::Prime(map)) => {
                let field = result
                    .syzygy_inclusion()
                    .source()
                    .right_acting_algebra()
                    .field();
                Ok(result.coefficients(map).map_err(invalid)?.map(|values| {
                    values
                        .iter()
                        .map(|value| field.wrap_scalar(value))
                        .collect()
                }))
            }
            (ExtOverAlgebraKind::Extension(result), ModuleHomomorphismKind::Extension(map)) => {
                let field = result
                    .syzygy_inclusion()
                    .source()
                    .right_acting_algebra()
                    .field();
                Ok(result.coefficients(map).map_err(invalid)?.map(|values| {
                    values
                        .iter()
                        .map(|value| field.wrap_scalar(value))
                        .collect()
                }))
            }
            _ => Err(invalid(
                "Ext result and homomorphism belong to different coefficient fields",
            )),
        }
    }
}

#[pymethods]
impl PyExtAlgebraGeneratorsResult {
    #[getter]
    fn dimensions(&self) -> Vec<usize> {
        dispatch_ref!(&self.kind, ExtAlgebraGeneratorsKind, |result| result
            .dimensions()
            .to_vec())
    }

    #[getter]
    fn generator_counts(&self) -> Vec<usize> {
        dispatch_ref!(&self.kind, ExtAlgebraGeneratorsKind, |result| result
            .generator_counts()
            .to_vec())
    }

    #[getter]
    fn generators(&self) -> Vec<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, ExtAlgebraGeneratorsKind, |result| result
            .generators()
            .iter()
            .map(|degree| degree
                .iter()
                .cloned()
                .map(|homomorphism| PyModuleHomomorphism {
                    kind: runtime::Wrap::wrap(homomorphism),
                })
                .collect())
            .collect())
    }
}

#[pymethods]
impl PyAddMCoresolution {
    #[getter]
    fn base_degree(&self) -> String {
        dispatch_ref!(&self.kind, AddMCoresolutionKind, |result| result
            .base_position()
            .to_string())
    }

    #[getter]
    fn differentials(&self) -> Vec<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, AddMCoresolutionKind, |result| result
            .differentials()
            .iter()
            .cloned()
            .map(|homomorphism| PyModuleHomomorphism {
                kind: runtime::Wrap::wrap(homomorphism),
            })
            .collect())
    }

    fn differential(&self, degree: i64) -> PyResult<PyModuleHomomorphism> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            AddMCoresolutionKind => ModuleHomomorphismKind,
            |result| result.differential(&degree).map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    fn object(&self, degree: i64) -> PyResult<PyMatrixModule> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            AddMCoresolutionKind => MatrixModuleKind,
            |result| result.object(&degree).map_err(invalid)?
        );
        Ok(PyMatrixModule { kind })
    }
}

#[pymethods]
impl PyFiniteAddMResolutionResult {
    #[getter]
    fn is_in_add(&self) -> bool {
        dispatch_ref!(&self.kind, FiniteAddMResolutionKind, |result| matches!(
            result,
            qpa_gbnp::FiniteAddMResolutionResult::InAddM
        ))
    }

    #[getter]
    fn differentials(&self) -> Vec<PyModuleHomomorphism> {
        dispatch_ref!(
            &self.kind,
            FiniteAddMResolutionKind,
            |result| match result {
                qpa_gbnp::FiniteAddMResolutionResult::InAddM => Vec::new(),
                qpa_gbnp::FiniteAddMResolutionResult::Resolution(resolution) => resolution
                    .differentials()
                    .iter()
                    .cloned()
                    .map(|homomorphism| PyModuleHomomorphism {
                        kind: runtime::Wrap::wrap(homomorphism),
                    })
                    .collect(),
            }
        )
    }

    fn differential(&self, degree: i64) -> PyResult<Option<PyModuleHomomorphism>> {
        let degree = Integer::from(degree);
        dispatch_ref!(
            &self.kind,
            FiniteAddMResolutionKind,
            |result| match result {
                qpa_gbnp::FiniteAddMResolutionResult::InAddM => Ok(None),
                qpa_gbnp::FiniteAddMResolutionResult::Resolution(resolution) => {
                    Ok(Some(PyModuleHomomorphism {
                        kind: runtime::Wrap::wrap(
                            resolution.differential(&degree).map_err(invalid)?,
                        ),
                    }))
                }
            }
        )
    }

    fn object(&self, degree: i64) -> PyResult<Option<PyMatrixModule>> {
        let degree = Integer::from(degree);
        dispatch_ref!(
            &self.kind,
            FiniteAddMResolutionKind,
            |result| match result {
                qpa_gbnp::FiniteAddMResolutionResult::InAddM => Ok(None),
                qpa_gbnp::FiniteAddMResolutionResult::Resolution(resolution) => {
                    Ok(Some(PyMatrixModule {
                        kind: runtime::Wrap::wrap(resolution.object(&degree).map_err(invalid)?),
                    }))
                }
            }
        )
    }
}

#[pymethods]
impl PyTiltingModuleResult {
    #[getter]
    fn projective_dimension(&self) -> usize {
        dispatch_ref!(&self.kind, TiltingModuleResultKind, |result| result
            .projective_dimension())
    }

    #[getter]
    fn coresolutions(&self) -> Vec<PyAddMCoresolution> {
        dispatch_ref!(&self.kind, TiltingModuleResultKind, |result| result
            .coresolutions()
            .iter()
            .cloned()
            .map(|coresolution| PyAddMCoresolution {
                kind: runtime::Wrap::wrap(coresolution),
            })
            .collect())
    }
}

#[pymethods]
impl PyCotiltingModuleResult {
    #[getter]
    fn injective_dimension(&self) -> usize {
        dispatch_ref!(&self.kind, CotiltingModuleResultKind, |result| result
            .injective_dimension())
    }

    #[getter]
    fn resolutions(&self) -> Vec<PyFiniteAddMResolutionResult> {
        dispatch_ref!(&self.kind, CotiltingModuleResultKind, |result| result
            .resolutions()
            .iter()
            .cloned()
            .map(|resolution| PyFiniteAddMResolutionResult {
                kind: runtime::Wrap::wrap(resolution),
            })
            .collect())
    }
}

#[pymethods]
impl PyAddMResolution {
    #[getter]
    fn differentials(&self) -> Vec<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, AddMResolutionKind, |result| result
            .differentials()
            .iter()
            .cloned()
            .map(|homomorphism| PyModuleHomomorphism {
                kind: runtime::Wrap::wrap(homomorphism),
            })
            .collect())
    }

    fn differential(&self, degree: i64) -> PyResult<PyModuleHomomorphism> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            AddMResolutionKind => ModuleHomomorphismKind,
            |result| result.differential(&degree).map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    fn object(&self, degree: i64) -> PyResult<PyMatrixModule> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            AddMResolutionKind => MatrixModuleKind,
            |result| result.object(&degree).map_err(invalid)?
        );
        Ok(PyMatrixModule { kind })
    }
}

#[pymethods]
impl PyAllComplementsResult {
    #[getter]
    fn resolution(&self) -> Option<PyAddMResolution> {
        dispatch_ref!(&self.kind, AllComplementsResultKind, |result| result
            .resolution()
            .cloned()
            .map(|resolution| PyAddMResolution {
                kind: runtime::Wrap::wrap(resolution),
            }))
    }

    #[getter]
    fn coresolution(&self) -> Option<PyAddMCoresolution> {
        dispatch_ref!(&self.kind, AllComplementsResultKind, |result| result
            .coresolution()
            .cloned()
            .map(|coresolution| PyAddMCoresolution {
                kind: runtime::Wrap::wrap(coresolution),
            }))
    }
}
