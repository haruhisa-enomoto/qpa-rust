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

// Concrete and common algebra handles expose the same invariant queries.
// Only the wrapper changes; all calculations retain the QPA dispatch below.
macro_rules! concrete_invariants {
    ($wrapper:ty, $kind:ident, $variant:ident) => {
        #[pymethods]
        impl $wrapper {
            /// Return QPA's Cartan matrix, or `None` outside its finite domain.
            fn cartan_matrix(&self) -> PyResult<Option<Vec<Vec<usize>>>> {
                dispatch_ref!(&self.kind, $kind, |algebra| {
                    cartan_for(&QuiverAlgebra::$variant(algebra.clone()))
                })
            }

            /// Return QPA's rational Coxeter matrix, or `None` when unavailable.
            fn coxeter_matrix(&self) -> PyResult<Option<Vec<Vec<PyScalar>>>> {
                dispatch_ref!(&self.kind, $kind, |algebra| {
                    coxeter_for(&QuiverAlgebra::$variant(algebra.clone()))
                })
            }

            /// Return Coxeter polynomial coefficients in increasing degree order.
            fn coxeter_polynomial(&self) -> PyResult<Option<Vec<PyScalar>>> {
                dispatch_ref!(&self.kind, $kind, |algebra| {
                    coxeter_polynomial_for(&QuiverAlgebra::$variant(algebra.clone()))
                })
            }

            fn is_selfinjective_algebra(&self) -> PyResult<Option<bool>> {
                dispatch_ref!(&self.kind, $kind, |algebra| {
                    selfinjective_for(&QuiverAlgebra::$variant(algebra.clone()))
                })
            }

            fn is_symmetric_algebra(&self) -> PyResult<Option<bool>> {
                dispatch_ref!(&self.kind, $kind, |algebra| {
                    symmetric_for(&QuiverAlgebra::$variant(algebra.clone()))
                })
            }

            fn is_weakly_symmetric_algebra(&self) -> PyResult<Option<bool>> {
                dispatch_ref!(&self.kind, $kind, |algebra| {
                    weakly_symmetric_for(&QuiverAlgebra::$variant(algebra.clone()))
                })
            }

            fn is_schurian_algebra(&self) -> PyResult<Option<bool>> {
                dispatch_ref!(&self.kind, $kind, |algebra| {
                    schurian_for(&QuiverAlgebra::$variant(algebra.clone()))
                })
            }

            fn is_semicommutative_algebra(&self) -> PyResult<Option<bool>> {
                dispatch_ref!(&self.kind, $kind, |algebra| {
                    semicommutative_for(&QuiverAlgebra::$variant(algebra.clone()))
                })
            }

            fn is_distributive_algebra(&self) -> PyResult<bool> {
                dispatch_ref!(&self.kind, $kind, |algebra| {
                    distributive_for(&QuiverAlgebra::$variant(algebra.clone()))
                })
            }
        }
    };
}

concrete_invariants!(PyPathAlgebra, AlgebraKind, Path);
concrete_invariants!(PyQuotientAlgebra, QuotientKind, Quotient);
