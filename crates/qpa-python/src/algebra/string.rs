//! Python dispatch for QPA's string-algebra enumeration algorithms.

use pyo3::pyfunction;

#[allow(clippy::wildcard_imports)]
// The path/quotient runtime union is owned by the parent binding module.
use super::*;

/// Translate QPA's `computeLPSarray` global helper.
#[pyfunction]
#[pyo3(name = "compute_lps_array")]
pub(crate) fn py_compute_lps_array(pattern: &str) -> Vec<usize> {
    qpa_gbnp::compute_lps_array(pattern)
}

/// Translate QPA's `CyclicPermutationOfABand` global helper.
#[pyfunction]
#[pyo3(name = "cyclic_permutations_of_band")]
pub(crate) fn py_cyclic_permutations_of_band(band: &str) -> Vec<String> {
    qpa_gbnp::cyclic_permutations_of_band(band)
}

pub(super) fn maximum_direct_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<usize> {
    match algebra {
        QuiverAlgebra::Path(algebra) => {
            qpa_gbnp::maximum_length_of_direct_string(algebra).map_err(invalid)
        }
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::maximum_length_of_direct_string(algebra).map_err(invalid)
        }
    }
}

pub(super) fn number_of_joints_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<usize> {
    match algebra {
        QuiverAlgebra::Path(algebra) => qpa_gbnp::number_of_joints(algebra).map_err(invalid),
        QuiverAlgebra::Quotient(algebra) => qpa_gbnp::number_of_joints(algebra).map_err(invalid),
    }
}

pub(super) fn is_a_band_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
    word: &str,
) -> PyResult<bool> {
    match algebra {
        QuiverAlgebra::Path(algebra) => qpa_gbnp::is_a_band(algebra, word).map_err(invalid),
        QuiverAlgebra::Quotient(algebra) => qpa_gbnp::is_a_band(algebra, word).map_err(invalid),
    }
}

pub(super) fn strings_less_than_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
    level: usize,
) -> PyResult<Vec<String>> {
    match algebra {
        QuiverAlgebra::Path(algebra) => {
            qpa_gbnp::strings_less_than(algebra, level).map_err(invalid)
        }
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::strings_less_than(algebra, level).map_err(invalid)
        }
    }
}

pub(super) fn bands_less_than_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
    level: usize,
) -> PyResult<Vec<String>> {
    match algebra {
        QuiverAlgebra::Path(algebra) => qpa_gbnp::bands_less_than(algebra, level).map_err(invalid),
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::bands_less_than(algebra, level).map_err(invalid)
        }
    }
}

pub(super) fn band_representatives_less_than_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
    level: usize,
) -> PyResult<Vec<String>> {
    match algebra {
        QuiverAlgebra::Path(algebra) => {
            qpa_gbnp::band_representatives_less_than(algebra, level).map_err(invalid)
        }
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::band_representatives_less_than(algebra, level).map_err(invalid)
        }
    }
}

pub(super) fn domestic_for<F: runtime::PythonField>(algebra: &QuiverAlgebra<F>) -> PyResult<bool> {
    match algebra {
        QuiverAlgebra::Path(algebra) => {
            qpa_gbnp::is_domestic_string_algebra(algebra).map_err(invalid)
        }
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::is_domestic_string_algebra(algebra).map_err(invalid)
        }
    }
}

pub(super) fn bridge_quiver_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<PyQuiver> {
    let quiver = match algebra {
        QuiverAlgebra::Path(algebra) => qpa_gbnp::bridge_quiver(algebra).map_err(invalid)?,
        QuiverAlgebra::Quotient(algebra) => qpa_gbnp::bridge_quiver(algebra).map_err(invalid)?,
    };
    Ok(PyQuiver::from_inner(quiver))
}
