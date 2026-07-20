//! Python dispatch for QPA's `eAe` and `A/AeA` quiver presentations.

use qpa_core::{QuiverAlgebra, QuiverAlgebraElement};

#[allow(clippy::wildcard_imports)]
// Runtime field variants and wrappers belong to the parent algebra domain.
use super::*;
use super::{
    conversion::{QuiverAlgebraKind, wrap_quiver_algebra},
    tensor::QuiverAlgebraElementKind,
};
use crate::field::parse_python_integer;

pub(super) fn vertex_positions(py: Python<'_>, positions: Vec<Py<PyAny>>) -> PyResult<Vec<usize>> {
    Ok(positions
        .into_iter()
        .map(|position| parse_python_integer(position.bind(py)))
        .collect::<PyResult<Vec<_>>>()?
        .into_iter()
        .filter_map(|position| position.to_string().parse::<usize>().ok())
        .collect())
}

pub(super) fn quotient_by_vertices_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
    positions: &[usize],
) -> PyResult<PyQuiverAlgebra>
where
    QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    let result = match algebra {
        QuiverAlgebra::Path(algebra) => {
            qpa_gbnp::quiver_algebra_of_path_algebra_mod_aea(algebra, positions).map_err(invalid)?
        }
        QuiverAlgebra::Quotient(algebra) => {
            qpa_gbnp::quiver_algebra_of_quotient_mod_aea(algebra, positions).map_err(invalid)?
        }
    };
    Ok(wrap_quiver_algebra(result))
}

pub(super) fn corner_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
    idempotent: QuiverAlgebraElementKind,
) -> PyResult<PyQuiverAlgebra>
where
    QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
    QuiverAlgebraElementKind: runtime::Extract<QuiverAlgebraElement<F>>,
{
    let idempotent = runtime::Extract::extract(idempotent)?;
    Ok(wrap_quiver_algebra(
        qpa_gbnp::quiver_algebra_of_eae(algebra, &idempotent).map_err(invalid)?,
    ))
}
