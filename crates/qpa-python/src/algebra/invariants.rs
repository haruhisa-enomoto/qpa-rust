//! Python conversion for QPA's algebra-invariant matrices and polynomials.

use numerica::domains::rational::Q;
use qpa_core::QuiverAlgebra;

#[allow(clippy::wildcard_imports)]
// Runtime fields and scalar wrappers belong to the parent module.
use super::*;

pub(super) fn cartan_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Option<Vec<Vec<usize>>>> {
    qpa_gbnp::cartan_matrix(algebra).map_err(invalid)
}

pub(super) fn coxeter_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Option<Vec<Vec<PyScalar>>>> {
    Ok(qpa_gbnp::coxeter_matrix(algebra)
        .map_err(invalid)?
        .map(|matrix| {
            (0..matrix.nrows())
                .map(|row| {
                    matrix
                        .row(row)
                        .iter()
                        .map(|coefficient| Q.wrap_scalar(coefficient))
                        .collect()
                })
                .collect()
        }))
}

pub(super) fn coxeter_polynomial_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Option<Vec<PyScalar>>> {
    Ok(qpa_gbnp::coxeter_polynomial(algebra)
        .map_err(invalid)?
        .map(|coefficients| {
            coefficients
                .iter()
                .map(|coefficient| Q.wrap_scalar(coefficient))
                .collect()
        }))
}

pub(super) fn selfinjective_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Option<bool>> {
    qpa_gbnp::is_selfinjective_algebra(algebra).map_err(invalid)
}

macro_rules! optional_property {
    ($name:ident, $qpa:path) => {
        pub(super) fn $name<F: runtime::PythonField>(
            algebra: &QuiverAlgebra<F>,
        ) -> PyResult<Option<bool>> {
            $qpa(algebra).map_err(invalid)
        }
    };
}

optional_property!(symmetric_for, qpa_gbnp::is_symmetric_algebra);
optional_property!(weakly_symmetric_for, qpa_gbnp::is_weakly_symmetric_algebra);
optional_property!(schurian_for, qpa_gbnp::is_schurian_algebra);
optional_property!(semicommutative_for, qpa_gbnp::is_semicommutative_algebra);

pub(super) fn distributive_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<bool> {
    match algebra {
        QuiverAlgebra::Path(algebra) => Ok(qpa_gbnp::path_algebra_is_distributive(algebra)),
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::quotient_is_distributive(algebra).map_err(invalid)
        }
    }
}
