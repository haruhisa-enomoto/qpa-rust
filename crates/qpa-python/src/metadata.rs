//! Python exposure for QPA's information classes and documentation plan.

use std::path::PathBuf;

use pyo3::prelude::{pyclass, pyfunction, pymethods};

/// Fixed inputs used by QPA's `MakeQPADocumentation` helper.
#[pyclass(name = "QpaDocumentationPlan", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQpaDocumentationPlan {
    inner: qpa_core::QpaDocumentationPlan,
}

#[pymethods]
impl PyQpaDocumentationPlan {
    #[getter]
    fn directory(&self) -> String {
        self.inner.directory().to_string_lossy().into_owned()
    }

    #[getter]
    fn main_file(&self) -> &'static str {
        self.inner.main_file()
    }

    #[getter]
    fn additional_files(&self) -> Vec<String> {
        self.inner
            .additional_files()
            .iter()
            .map(|path| path.to_string_lossy().into_owned())
            .collect()
    }

    #[getter]
    fn book_name(&self) -> &'static str {
        self.inner.book_name()
    }

    #[getter]
    fn renderer(&self) -> &'static str {
        self.inner.renderer()
    }

    #[getter]
    fn copy_html_style_files(&self) -> bool {
        self.inner.copy_html_style_files()
    }

    fn __repr__(&self) -> String {
        format!(
            "<QPA documentation plan: {}/{}>",
            self.directory(),
            self.main_file()
        )
    }
}

/// Return the source-faithful build inputs of QPA's `MakeQPADocumentation`.
#[pyfunction(name = "make_qpa_documentation")]
pub(crate) fn py_make_qpa_documentation(package_root: PathBuf) -> PyQpaDocumentationPlan {
    PyQpaDocumentationPlan {
        inner: qpa_core::make_qpa_documentation(package_root),
    }
}

/// Name of QPA's `InfoGroebnerBasis` information class.
#[pyfunction(name = "info_groebner_basis_name")]
pub(crate) const fn py_info_groebner_basis_name() -> &'static str {
    qpa_core::INFO_GROEBNER_BASIS
}

/// Name of QPA's `InfoPathAlgebraModule` information class.
#[pyfunction(name = "info_path_algebra_module_name")]
pub(crate) const fn py_info_path_algebra_module_name() -> &'static str {
    qpa_core::INFO_PATH_ALGEBRA_MODULE
}

/// Initial `InfoLevel` installed for `InfoPathAlgebraModule`.
#[pyfunction(name = "info_path_algebra_module_default_level")]
pub(crate) const fn py_info_path_algebra_module_default_level() -> usize {
    qpa_core::INFO_PATH_ALGEBRA_MODULE_DEFAULT_LEVEL
}

/// Name of QPA's projective-resolution information class.
#[pyfunction(name = "info_projective_resolution_fp_path_algebra_module_name")]
pub(crate) const fn py_info_projective_resolution_name() -> &'static str {
    qpa_core::INFO_PROJECTIVE_RESOLUTION_FP_PATH_ALGEBRA_MODULE
}

/// Initial `InfoLevel` installed for the projective-resolution information class.
#[pyfunction(name = "info_projective_resolution_fp_path_algebra_module_default_level")]
pub(crate) const fn py_info_projective_resolution_default_level() -> usize {
    qpa_core::INFO_PROJECTIVE_RESOLUTION_FP_PATH_ALGEBRA_MODULE_DEFAULT_LEVEL
}

/// Name of QPA's quiver information class.
#[pyfunction(name = "info_quiver_name")]
pub(crate) const fn py_info_quiver_name() -> &'static str {
    qpa_core::INFO_QUIVER
}
