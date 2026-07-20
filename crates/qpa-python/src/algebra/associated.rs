//! Python dispatch for QPA's associated monomial and graded algebras.

use qpa_core::QuiverAlgebra;

use super::conversion::QuiverAlgebraKind;
#[allow(clippy::wildcard_imports)]
// Runtime fields and the quiver-algebra wrapper belong to the parent domain.
use super::*;

pub(super) fn monomial_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<PyQuiverAlgebra>
where
    QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    let associated = match algebra {
        QuiverAlgebra::Path(algebra) => qpa_gbnp::associated_monomial_path_algebra(algebra),
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::associated_monomial_quotient(algebra).map_err(invalid)?
        }
    };
    Ok(conversion::wrap_quiver_algebra(associated))
}

pub(super) fn graded_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<PyQuiverAlgebra>
where
    QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    let associated = match algebra {
        QuiverAlgebra::Path(algebra) => {
            qpa_gbnp::associated_graded_path_algebra(algebra).map_err(invalid)?
        }
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::associated_graded_quotient(algebra).map_err(invalid)?
        }
    };
    Ok(conversion::wrap_quiver_algebra(associated))
}
