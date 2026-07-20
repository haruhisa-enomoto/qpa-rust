//! Python surface for QPA's GSZ path-vector helper operations.

use numerica::domains::{
    Field,
    finite_field::{Z2, Zp},
    rational::RationalField,
};
use pyo3::{Py, PyResult, Python, prelude::*};
use qpa_core::{
    PathAlgebraElement, PathAlgebraIdeal, PathAlgebraVector, PathAlgebraVectorXSet, first_part,
    left_division, tip_reduce_path_algebra_vector, tip_reduce_path_algebra_vectors,
    x_set_of_path_algebra_vector,
};

use super::projective::ProjectiveVectorKind;
#[allow(clippy::wildcard_imports)]
use super::*;

#[derive(Clone)]
enum PathAlgebraVectorXSetKind {
    Rational(PathAlgebraVectorXSet<RationalField>),
    Binary(PathAlgebraVectorXSet<Z2>),
    Prime(PathAlgebraVectorXSet<Zp>),
}

/// QPA's positional O-set/N-set result from `XSetOfPathAlgebraVector`.
#[pyclass(name = "PathAlgebraVectorXSet", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPathAlgebraVectorXSet {
    kind: PathAlgebraVectorXSetKind,
}

macro_rules! impl_x_set_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<PathAlgebraVectorXSet<$field>> for PathAlgebraVectorXSetKind {
            fn wrap(value: PathAlgebraVectorXSet<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_x_set_runtime!(RationalField, Rational);
impl_x_set_runtime!(Z2, Binary);
impl_x_set_runtime!(Zp, Prime);

fn clone_vector_kinds(
    py: Python<'_>,
    vectors: Vec<Py<PyPathAlgebraVector>>,
) -> Vec<ProjectiveVectorKind> {
    vectors
        .into_iter()
        .map(|vector| vector.borrow(py).kind.clone())
        .collect()
}

fn left_division_for<F: Field>(
    dividend: &PathAlgebraVector<F>,
    divisor: &PathAlgebraVector<F>,
) -> PyResult<Option<PyAlgebraElement>>
where
    AlgebraElementKind: runtime::Wrap<PathAlgebraElement<F>>,
{
    Ok(left_division(dividend, divisor)
        .map_err(invalid)?
        .map(|value| PyAlgebraElement {
            kind: runtime::Wrap::wrap(value),
        }))
}

fn tip_reduce_list_for<F: Field>(
    vectors: Vec<ProjectiveVectorKind>,
) -> PyResult<Vec<PyPathAlgebraVector>>
where
    ProjectiveVectorKind:
        runtime::Extract<PathAlgebraVector<F>> + runtime::Wrap<PathAlgebraVector<F>>,
{
    Ok(projective::wrap_vectors(
        tip_reduce_path_algebra_vectors(&runtime::extract_all(vectors)?).map_err(invalid)?,
    ))
}

fn tip_reduce_vector_for<F: Field>(
    vector: &PathAlgebraVector<F>,
    reducers: Vec<ProjectiveVectorKind>,
) -> PyResult<PyPathAlgebraVector>
where
    ProjectiveVectorKind:
        runtime::Extract<PathAlgebraVector<F>> + runtime::Wrap<PathAlgebraVector<F>>,
{
    Ok(PyPathAlgebraVector {
        kind: runtime::Wrap::wrap(
            tip_reduce_path_algebra_vector(&runtime::extract_all(reducers)?, vector)
                .map_err(invalid)?,
        ),
    })
}

fn first_part_for<F: runtime::PythonField>(
    vector: &PathAlgebraVector<F>,
    f: Vec<ProjectiveVectorKind>,
    f_prime: Vec<ProjectiveVectorKind>,
) -> PyResult<Vec<PyAlgebraElement>>
where
    ProjectiveVectorKind: runtime::Extract<PathAlgebraVector<F>>,
    AlgebraElementKind: runtime::Wrap<PathAlgebraElement<F>>,
{
    Ok(wrap_algebra_elements(
        first_part(
            &runtime::extract_all(f)?,
            &runtime::extract_all(f_prime)?,
            vector,
        )
        .map_err(invalid)?,
    ))
}

fn x_set_for<F: Field>(
    vector: &PathAlgebraVector<F>,
    ideal: &PathAlgebraIdeal<F>,
    f: Vec<ProjectiveVectorKind>,
    f_prime: Vec<ProjectiveVectorKind>,
) -> PyResult<PyPathAlgebraVectorXSet>
where
    ProjectiveVectorKind: runtime::Extract<PathAlgebraVector<F>>,
    PathAlgebraVectorXSetKind: runtime::Wrap<PathAlgebraVectorXSet<F>>,
{
    Ok(PyPathAlgebraVectorXSet {
        kind: runtime::Wrap::wrap(
            x_set_of_path_algebra_vector(
                &runtime::extract_all(f)?,
                &runtime::extract_all(f_prime)?,
                ideal,
                vector,
            )
            .map_err(invalid)?,
        ),
    })
}

pub(super) fn vector_left_division(
    vector: &PyPathAlgebraVector,
    divisor: &PyPathAlgebraVector,
) -> PyResult<Option<PyAlgebraElement>> {
    match (&vector.kind, &divisor.kind) {
        (ProjectiveVectorKind::Rational(left), ProjectiveVectorKind::Rational(right)) => {
            left_division_for(left, right)
        }
        (ProjectiveVectorKind::Binary(left), ProjectiveVectorKind::Binary(right)) => {
            left_division_for(left, right)
        }
        (ProjectiveVectorKind::Prime(left), ProjectiveVectorKind::Prime(right)) => {
            left_division_for(left, right)
        }
        _ => Err(invalid(
            "path-algebra vectors use different coefficient fields",
        )),
    }
}

pub(super) fn tip_reduce(
    py: Python<'_>,
    vectors: Vec<Py<PyPathAlgebraVector>>,
) -> PyResult<Vec<PyPathAlgebraVector>> {
    let vectors = clone_vector_kinds(py, vectors);
    let Some(first) = vectors.first() else {
        return Err(invalid(
            "QPA has no TipReduce method for an empty path-vector list",
        ));
    };
    match first {
        ProjectiveVectorKind::Rational(_) => tip_reduce_list_for::<RationalField>(vectors),
        ProjectiveVectorKind::Binary(_) => tip_reduce_list_for::<Z2>(vectors),
        ProjectiveVectorKind::Prime(_) => tip_reduce_list_for::<Zp>(vectors),
    }
}

pub(super) fn tip_reduce_by(
    vector: &PyPathAlgebraVector,
    py: Python<'_>,
    reducers: Vec<Py<PyPathAlgebraVector>>,
) -> PyResult<PyPathAlgebraVector> {
    let reducers = clone_vector_kinds(py, reducers);
    match &vector.kind {
        ProjectiveVectorKind::Rational(vector) => tip_reduce_vector_for(vector, reducers),
        ProjectiveVectorKind::Binary(vector) => tip_reduce_vector_for(vector, reducers),
        ProjectiveVectorKind::Prime(vector) => tip_reduce_vector_for(vector, reducers),
    }
}

pub(super) fn first_part_of(
    vector: &PyPathAlgebraVector,
    py: Python<'_>,
    f: Vec<Py<PyPathAlgebraVector>>,
    f_prime: Vec<Py<PyPathAlgebraVector>>,
) -> PyResult<Vec<PyAlgebraElement>> {
    let f = clone_vector_kinds(py, f);
    let f_prime = clone_vector_kinds(py, f_prime);
    match &vector.kind {
        ProjectiveVectorKind::Rational(vector) => first_part_for(vector, f, f_prime),
        ProjectiveVectorKind::Binary(vector) => first_part_for(vector, f, f_prime),
        ProjectiveVectorKind::Prime(vector) => first_part_for(vector, f, f_prime),
    }
}

pub(super) fn x_set(
    vector: &PyPathAlgebraVector,
    py: Python<'_>,
    f: Vec<Py<PyPathAlgebraVector>>,
    f_prime: Vec<Py<PyPathAlgebraVector>>,
    ideal: &PyPathAlgebraIdeal,
) -> PyResult<PyPathAlgebraVectorXSet> {
    let f = clone_vector_kinds(py, f);
    let f_prime = clone_vector_kinds(py, f_prime);
    match (&vector.kind, &ideal.kind) {
        (ProjectiveVectorKind::Rational(vector), IdealKind::Rational(ideal)) => {
            x_set_for(vector, ideal, f, f_prime)
        }
        (ProjectiveVectorKind::Binary(vector), IdealKind::Binary(ideal)) => {
            x_set_for(vector, ideal, f, f_prime)
        }
        (ProjectiveVectorKind::Prime(vector), IdealKind::Prime(ideal)) => {
            x_set_for(vector, ideal, f, f_prime)
        }
        _ => Err(invalid(
            "the path-algebra vector and ideal use different coefficient fields",
        )),
    }
}

#[pymethods]
impl PyPathAlgebraVectorXSet {
    #[getter]
    fn o_set(&self) -> Vec<(PyPath, usize, PyAlgebraElement, PyAlgebraElement)> {
        dispatch_ref!(&self.kind, PathAlgebraVectorXSetKind, |value| value
            .o_set
            .iter()
            .map(|entry| (
                PyPath {
                    inner: entry.right_tip.clone(),
                },
                entry.generator_position,
                PyAlgebraElement {
                    kind: runtime::Wrap::wrap(entry.q_prime.clone()),
                },
                PyAlgebraElement {
                    kind: runtime::Wrap::wrap(entry.relation.clone()),
                },
            ))
            .collect())
    }

    #[getter]
    fn n_set(&self) -> Vec<(PyPath, PyAlgebraElement, PyAlgebraElement)> {
        dispatch_ref!(&self.kind, PathAlgebraVectorXSetKind, |value| value
            .n_set
            .iter()
            .map(|entry| (
                PyPath {
                    inner: entry.right_tip.clone(),
                },
                PyAlgebraElement {
                    kind: runtime::Wrap::wrap(entry.z.clone()),
                },
                PyAlgebraElement {
                    kind: runtime::Wrap::wrap(entry.relation.clone()),
                },
            ))
            .collect())
    }
}
