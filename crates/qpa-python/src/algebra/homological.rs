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
