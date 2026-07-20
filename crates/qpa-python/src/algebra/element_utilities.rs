//! Python dispatch for QPA's foundational path-algebra element utilities.

#[allow(clippy::wildcard_imports)]
// Element and algebra runtime variants belong to the parent binding domain.
use super::*;

#[derive(Clone, Copy)]
pub(super) enum SupportSide {
    Left,
    Right,
    Both,
}

pub(super) fn path_support_for<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    element: AlgebraElementKind,
    side: SupportSide,
) -> PyResult<Vec<Vec<usize>>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let element = runtime::Extract::extract(element)?;
    match side {
        SupportSide::Left => Ok(vec![
            algebra
                .left_support_of_quiver_algebra_element(&element)
                .map_err(invalid)?,
        ]),
        SupportSide::Right => Ok(vec![
            algebra
                .right_support_of_quiver_algebra_element(&element)
                .map_err(invalid)?,
        ]),
        SupportSide::Both => Ok(algebra
            .support_of_quiver_algebra_element(&element)
            .map_err(invalid)?
            .into()),
    }
}

pub(super) fn quotient_support_for<F: runtime::PythonField>(
    algebra: &QuotientPathAlgebra<F>,
    element: QuotientElementKind,
    side: SupportSide,
) -> PyResult<Vec<Vec<usize>>>
where
    QuotientElementKind: runtime::Extract<QuotientPathAlgebraElement<F>>,
{
    let element = runtime::Extract::extract(element)?;
    match side {
        SupportSide::Left => Ok(vec![
            algebra
                .left_support_of_quiver_algebra_element(&element)
                .map_err(invalid)?,
        ]),
        SupportSide::Right => Ok(vec![
            algebra
                .right_support_of_quiver_algebra_element(&element)
                .map_err(invalid)?,
        ]),
        SupportSide::Both => Ok(algebra
            .support_of_quiver_algebra_element(&element)
            .map_err(invalid)?
            .into()),
    }
}

pub(super) fn path_removal_for<F: runtime::PythonField>(
    element: &PathAlgebraElement<F>,
    positions: &[usize],
) -> PyResult<PyAlgebraElement>
where
    AlgebraElementKind: runtime::Wrap<PathAlgebraElement<F>>,
{
    Ok(PyAlgebraElement {
        kind: runtime::Wrap::wrap(element.path_removal(positions).map_err(invalid)?),
    })
}

pub(super) fn add_power_for<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
    exponent: i64,
) -> PyResult<Vec<PyAlgebraElement>>
where
    AlgebraElementKind:
        runtime::Extract<PathAlgebraElement<F>> + runtime::Wrap<PathAlgebraElement<F>>,
{
    let mut relations = runtime::extract_all(relations)?;
    let existing = relations.len();
    algebra
        .add_nth_power_to_relations(&mut relations, exponent)
        .map_err(invalid)?;
    Ok(relations
        .into_iter()
        .skip(existing)
        .map(|element| PyAlgebraElement {
            kind: runtime::Wrap::wrap(element),
        })
        .collect())
}

fn opposite_relations_for<F: runtime::PythonField>(
    relations: Vec<AlgebraElementKind>,
) -> PyResult<Vec<PyAlgebraElement>>
where
    AlgebraElementKind:
        runtime::Extract<PathAlgebraElement<F>> + runtime::Wrap<PathAlgebraElement<F>>,
{
    Ok(wrap_algebra_elements(
        qpa_core::opposite_relations(&runtime::extract_all(relations)?).map_err(invalid)?,
    ))
}

/// Apply QPA's elementwise `OppositeRelations` list operation.
#[pyfunction(name = "opposite_relations")]
pub(crate) fn py_opposite_relations(
    py: Python<'_>,
    relations: Vec<Py<PyAlgebraElement>>,
) -> PyResult<Vec<PyAlgebraElement>> {
    let relations = clone_relations(py, relations);
    match relations.first() {
        None => Ok(Vec::new()),
        Some(AlgebraElementKind::Rational(_)) => opposite_relations_for::<RationalField>(relations),
        Some(AlgebraElementKind::Binary(_)) => opposite_relations_for::<Z2>(relations),
        Some(AlgebraElementKind::Prime(_)) => opposite_relations_for::<Zp>(relations),
    }
}
