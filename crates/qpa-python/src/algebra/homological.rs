//! Lazy projective and injective resolutions.

use pyo3::pymethods;

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

#[pymethods]
impl PyProjectiveResolution {
    fn differential(&self, degree: i64) -> PyResult<PyModuleHomomorphism> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            ProjectiveResolutionKind => ModuleHomomorphismKind,
            |resolution| resolution.differential(&degree).map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

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

#[pymethods]
impl PyModuleComplex {
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

    fn __repr__(&self) -> PyResult<String> {
        dispatch_ref!(&self.kind, ModuleComplexKind, |complex| complex
            .qpa_string()
            .map_err(invalid))
    }
}
