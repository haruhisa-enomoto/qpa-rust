//! Python wrapper for QPA's source-translated quiver orderings.

use pyo3::{Py, PyAny, PyResult, Python, prelude::*};
use qpa_core::{
    ArrowIndex, Generator, QuiverOrdering, QuiverOrderingKind, VertexIndex,
    check_blocks_for_cycles, check_quiver_subset_for_cycles,
};

use super::{PyPath, PyQuiver, objects::generator_to_python, ordering_value};
use crate::{PyScalar, ScalarKind, invalid};

/// A QPA quiver ordering, including its comparison closure and object identity.
#[pyclass(name = "QuiverOrdering", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuiverOrdering {
    pub(crate) inner: QuiverOrdering,
}

impl PyQuiverOrdering {
    pub(crate) const fn from_inner(inner: QuiverOrdering) -> Self {
        Self { inner }
    }
}

fn optional_ordering(
    py: Python<'_>,
    ordering: Option<Py<PyQuiverOrdering>>,
) -> Option<QuiverOrdering> {
    ordering.map(|ordering| ordering.borrow(py).inner.clone())
}

pub(super) fn generators_from_positions(
    quiver: &PyQuiver,
    positions: &[usize],
) -> PyResult<Vec<Generator>> {
    let number_of_vertices = quiver.inner.number_of_vertices();
    let number_of_generators = number_of_vertices + quiver.inner.number_of_arrows();
    positions
        .iter()
        .map(|&position| {
            if position == 0 || position > number_of_generators {
                return Err(invalid(format!(
                    "generator index {position} is outside 1..={number_of_generators}"
                )));
            }
            if position <= number_of_vertices {
                let index = VertexIndex::new(position).expect("positive generator position");
                Ok(Generator::Vertex(
                    quiver
                        .inner
                        .vertex(index)
                        .expect("validated vertex generator position"),
                ))
            } else {
                let index = ArrowIndex::new(position - number_of_vertices)
                    .expect("positive arrow position");
                Ok(Generator::Arrow(
                    quiver
                        .inner
                        .arrow(index)
                        .expect("validated arrow generator position"),
                ))
            }
        })
        .collect()
}

fn ordering_kind_name(kind: QuiverOrderingKind) -> &'static str {
    match kind {
        QuiverOrderingKind::Length => "length",
        QuiverOrderingKind::LeftLexicographic => "left_lexicographic",
        QuiverOrderingKind::RightLexicographic => "right_lexicographic",
        QuiverOrderingKind::Reverse => "reverse",
        QuiverOrderingKind::LeftVector => "left_vector",
        QuiverOrderingKind::RightVector => "right_vector",
        QuiverOrderingKind::Weight => "weight",
        QuiverOrderingKind::Block => "block",
        QuiverOrderingKind::Wreath => "wreath",
    }
}

#[pymethods]
impl PyQuiverOrdering {
    #[staticmethod]
    fn default(quiver: &PyQuiver) -> PyResult<Self> {
        Ok(Self::from_inner(
            QuiverOrdering::default_for(&quiver.inner).map_err(invalid)?,
        ))
    }

    #[staticmethod]
    #[pyo3(signature = (quiver, next=None))]
    fn length(quiver: &PyQuiver, next: Option<Py<Self>>, py: Python<'_>) -> PyResult<Self> {
        Ok(Self::from_inner(
            QuiverOrdering::length(&quiver.inner, optional_ordering(py, next)).map_err(invalid)?,
        ))
    }

    #[staticmethod]
    #[pyo3(signature = (quiver, preferred=Vec::new()))]
    fn left_lexicographic(quiver: &PyQuiver, preferred: Vec<usize>) -> PyResult<Self> {
        Ok(Self::from_inner(
            QuiverOrdering::left_lexicographic(
                &quiver.inner,
                &generators_from_positions(quiver, &preferred)?,
            )
            .map_err(invalid)?,
        ))
    }

    #[staticmethod]
    #[pyo3(signature = (quiver, preferred=Vec::new()))]
    fn right_lexicographic(quiver: &PyQuiver, preferred: Vec<usize>) -> PyResult<Self> {
        Ok(Self::from_inner(
            QuiverOrdering::right_lexicographic(
                &quiver.inner,
                &generators_from_positions(quiver, &preferred)?,
            )
            .map_err(invalid)?,
        ))
    }

    #[staticmethod]
    fn reverse(quiver: &PyQuiver, next: &Self) -> PyResult<Self> {
        Ok(Self::from_inner(
            QuiverOrdering::reverse(&quiver.inner, next.inner.clone()).map_err(invalid)?,
        ))
    }

    #[staticmethod]
    #[pyo3(signature = (quiver, preferred=Vec::new(), next=None))]
    fn left_vector(
        quiver: &PyQuiver,
        preferred: Vec<usize>,
        next: Option<Py<Self>>,
        py: Python<'_>,
    ) -> PyResult<Self> {
        Ok(Self::from_inner(
            QuiverOrdering::left_vector(
                &quiver.inner,
                &generators_from_positions(quiver, &preferred)?,
                optional_ordering(py, next),
            )
            .map_err(invalid)?,
        ))
    }

    #[staticmethod]
    #[pyo3(signature = (quiver, preferred=Vec::new(), next=None))]
    fn right_vector(
        quiver: &PyQuiver,
        preferred: Vec<usize>,
        next: Option<Py<Self>>,
        py: Python<'_>,
    ) -> PyResult<Self> {
        Ok(Self::from_inner(
            QuiverOrdering::right_vector(
                &quiver.inner,
                &generators_from_positions(quiver, &preferred)?,
                optional_ordering(py, next),
            )
            .map_err(invalid)?,
        ))
    }

    #[staticmethod]
    #[pyo3(signature = (quiver, weights, next=None))]
    fn weight(
        quiver: &PyQuiver,
        weights: Vec<Py<PyScalar>>,
        next: Option<Py<Self>>,
        py: Python<'_>,
    ) -> PyResult<Self> {
        let weights = weights
            .into_iter()
            .map(|weight| match &weight.borrow(py).value {
                ScalarKind::Rational(value) => Ok(value.clone()),
                ScalarKind::Prime { .. } | ScalarKind::Extension { .. } => Err(invalid(
                    "WeightOrdering requires exact rational Scalar values",
                )),
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self::from_inner(
            QuiverOrdering::weight(&quiver.inner, &weights, optional_ordering(py, next))
                .map_err(invalid)?,
        ))
    }

    #[staticmethod]
    #[pyo3(signature = (quiver, blocks, next=None))]
    fn block(
        quiver: &PyQuiver,
        blocks: Vec<(Py<Self>, Vec<usize>)>,
        next: Option<Py<Self>>,
        py: Python<'_>,
    ) -> PyResult<Self> {
        let blocks = blocks
            .into_iter()
            .map(|(ordering, positions)| {
                Ok((
                    ordering.borrow(py).inner.clone(),
                    generators_from_positions(quiver, &positions)?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self::from_inner(
            QuiverOrdering::block(&quiver.inner, blocks, optional_ordering(py, next))
                .map_err(invalid)?,
        ))
    }

    #[staticmethod]
    #[pyo3(signature = (quiver, blocks, next=None))]
    fn wreath(
        quiver: &PyQuiver,
        blocks: Vec<(Py<Self>, Vec<usize>)>,
        next: Option<Py<Self>>,
        py: Python<'_>,
    ) -> PyResult<Self> {
        let blocks = blocks
            .into_iter()
            .map(|(ordering, positions)| {
                Ok((
                    ordering.borrow(py).inner.clone(),
                    generators_from_positions(quiver, &positions)?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self::from_inner(
            QuiverOrdering::wreath(&quiver.inner, blocks, optional_ordering(py, next))
                .map_err(invalid)?,
        ))
    }

    #[getter]
    fn quiver(&self) -> PyQuiver {
        PyQuiver::from_inner(self.inner.quiver().clone())
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.name().to_owned()
    }

    #[getter]
    fn kind(&self) -> &'static str {
        ordering_kind_name(self.inner.kind())
    }

    #[getter]
    fn is_length(&self) -> bool {
        self.inner.kind() == QuiverOrderingKind::Length
    }

    #[getter]
    fn is_left_lexicographic(&self) -> bool {
        self.inner.kind() == QuiverOrderingKind::LeftLexicographic
    }

    #[getter]
    fn is_right_lexicographic(&self) -> bool {
        self.inner.kind() == QuiverOrderingKind::RightLexicographic
    }

    #[getter]
    fn is_lexicographic(&self) -> bool {
        self.inner.kind().is_lexicographic()
    }

    #[getter]
    fn is_reverse(&self) -> bool {
        self.inner.kind() == QuiverOrderingKind::Reverse
    }

    #[getter]
    fn is_left_vector(&self) -> bool {
        self.inner.kind() == QuiverOrderingKind::LeftVector
    }

    #[getter]
    fn is_right_vector(&self) -> bool {
        self.inner.kind() == QuiverOrderingKind::RightVector
    }

    #[getter]
    fn is_vector(&self) -> bool {
        self.inner.kind().is_vector()
    }

    #[getter]
    fn is_weight(&self) -> bool {
        self.inner.kind() == QuiverOrderingKind::Weight
    }

    #[getter]
    fn is_block(&self) -> bool {
        self.inner.kind() == QuiverOrderingKind::Block
    }

    #[getter]
    fn is_wreath(&self) -> bool {
        self.inner.kind() == QuiverOrderingKind::Wreath
    }

    #[getter]
    fn next_ordering(&self) -> Option<Self> {
        self.inner.next_ordering().cloned().map(Self::from_inner)
    }

    #[getter]
    fn lexicographic_table(&self, py: Python<'_>) -> PyResult<Option<Vec<Py<PyAny>>>> {
        self.inner
            .lexicographic_table()
            .map(|table| {
                table
                    .into_iter()
                    .map(|generator| generator_to_python(py, generator))
                    .collect()
            })
            .transpose()
    }

    fn shares_lexicographic_table_with(&self, other: &Self) -> bool {
        self.inner.shares_lexicographic_table_with(&other.inner)
    }

    fn compare(&self, left: &PyPath, right: &PyPath) -> PyResult<i8> {
        Ok(ordering_value(
            self.inner
                .compare(&left.inner, &right.inner)
                .map_err(invalid)?,
        ))
    }

    fn comparison_value(&self, left: &PyPath, right: &PyPath) -> PyResult<PyScalar> {
        Ok(PyScalar {
            value: ScalarKind::Rational(
                self.inner
                    .comparison_value(&left.inner, &right.inner)
                    .map_err(invalid)?,
            ),
        })
    }

    fn less_than(&self, left: &PyPath, right: &PyPath) -> PyResult<bool> {
        self.inner
            .less_than(&left.inner, &right.inner)
            .map_err(invalid)
    }

    fn compare_arrow_lists(&self, left: Vec<usize>, right: Vec<usize>) -> PyResult<i8> {
        let quiver = PyQuiver::from_inner(self.inner.quiver().clone());
        Ok(ordering_value(
            self.inner
                .compare_arrow_lists(
                    &generators_from_positions(&quiver, &left)?,
                    &generators_from_positions(&quiver, &right)?,
                )
                .map_err(invalid)?,
        ))
    }

    #[getter]
    fn is_well(&self) -> PyResult<bool> {
        self.inner.is_well_ordering().map_err(invalid)
    }

    #[getter]
    fn is_well_reversed(&self) -> PyResult<bool> {
        self.inner.is_well_reversed_ordering().map_err(invalid)
    }

    #[getter]
    fn is_total(&self) -> PyResult<bool> {
        self.inner.is_total_ordering().map_err(invalid)
    }

    #[getter]
    fn is_admissible(&self) -> PyResult<bool> {
        self.inner.is_admissible_ordering().map_err(invalid)
    }

    fn is_well_subset(&self, positions: Vec<usize>) -> PyResult<bool> {
        let quiver = PyQuiver::from_inner(self.inner.quiver().clone());
        self.inner
            .is_well_subset(&generators_from_positions(&quiver, &positions)?)
            .map_err(invalid)
    }

    fn is_well_reversed_subset(&self, positions: Vec<usize>) -> PyResult<bool> {
        let quiver = PyQuiver::from_inner(self.inner.quiver().clone());
        self.inner
            .is_well_reversed_subset(&generators_from_positions(&quiver, &positions)?)
            .map_err(invalid)
    }

    fn is_total_subset(&self, positions: Vec<usize>) -> PyResult<bool> {
        let quiver = PyQuiver::from_inner(self.inner.quiver().clone());
        self.inner
            .is_total_subset(&generators_from_positions(&quiver, &positions)?)
            .map_err(invalid)
    }

    fn check_subset_for_cycles(&self, positions: Vec<usize>) -> PyResult<bool> {
        let quiver = PyQuiver::from_inner(self.inner.quiver().clone());
        check_quiver_subset_for_cycles(
            self.inner.quiver(),
            &generators_from_positions(&quiver, &positions)?,
        )
        .map_err(invalid)
    }

    fn check_blocks_for_cycles(&self, blocks: Vec<usize>) -> PyResult<bool> {
        check_blocks_for_cycles(&self.inner, &blocks).map_err(invalid)
    }

    fn is_identical(&self, other: &Self) -> bool {
        self.inner.is_identical(&other.inner)
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    fn print_obj(&self) -> String {
        self.__repr__()
    }

    fn view(&self) -> String {
        self.__repr__()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}
