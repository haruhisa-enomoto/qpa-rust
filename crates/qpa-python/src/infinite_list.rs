//! Python wrappers for QPA's lazy one- and two-sided infinite lists.

use std::{fmt, sync::Arc};

use numerica::domains::integer::Integer;
use pyo3::{Bound, Py, PyAny, PyErr, PyResult, Python, prelude::*};
use qpa_core::{
    ExtendedInteger, HalfInfiniteList, HalfInfiniteListElementFunction,
    HalfInfiniteListSpecification, HalfInfiniteListType, InfiniteList, InfiniteListDirection,
    InfiniteListError, InfiniteListNewValueCallback, InfiniteListNextFunction,
    InfiniteListPositionFunction, InfiniteListRepeatifyCallback, shifted_splice,
};

use crate::{
    field::{parse_python_integer, wrap_python_integer},
    invalid,
};

#[derive(Clone)]
struct PythonValue(Arc<Py<PyAny>>);

impl PythonValue {
    fn new(value: Py<PyAny>) -> Self {
        Self(Arc::new(value))
    }

    fn to_python(&self, py: Python<'_>) -> Py<PyAny> {
        self.0.as_ref().clone_ref(py)
    }
}

impl PartialEq for PythonValue {
    fn eq(&self, other: &Self) -> bool {
        Python::attach(|py| self.0.bind(py).eq(other.0.bind(py)).unwrap_or(false))
    }
}

#[derive(Debug)]
enum PythonInfiniteListError {
    Core(InfiniteListError),
    Python(PyErr),
}

impl From<InfiniteListError> for PythonInfiniteListError {
    fn from(value: InfiniteListError) -> Self {
        Self::Core(value)
    }
}

impl fmt::Display for PythonInfiniteListError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Core(error) => error.fmt(formatter),
            Self::Python(error) => error.fmt(formatter),
        }
    }
}

type PythonHalfList = HalfInfiniteList<PythonValue, PythonInfiniteListError>;
type PythonFullList = InfiniteList<PythonValue, PythonInfiniteListError>;

fn python_error(error: PythonInfiniteListError) -> PyErr {
    match error {
        PythonInfiniteListError::Core(error) => invalid(error),
        PythonInfiniteListError::Python(error) => error,
    }
}

fn parse_direction(direction: i8) -> PyResult<InfiniteListDirection> {
    match direction {
        1 => Ok(InfiniteListDirection::Positive),
        -1 => Ok(InfiniteListDirection::Negative),
        _ => Err(invalid("infinite-list direction must be 1 or -1")),
    }
}

fn direction_value(direction: InfiniteListDirection) -> i8 {
    match direction {
        InfiniteListDirection::Positive => 1,
        InfiniteListDirection::Negative => -1,
    }
}

fn type_name(list_type: HalfInfiniteListType) -> &'static str {
    match list_type {
        HalfInfiniteListType::Repeat => "repeat",
        HalfInfiniteListType::Next => "next",
        HalfInfiniteListType::Position => "pos",
    }
}

fn parse_type_name(value: &str) -> PyResult<HalfInfiniteListType> {
    match value {
        "repeat" => Ok(HalfInfiniteListType::Repeat),
        "next" | "next/repeat" => Ok(HalfInfiniteListType::Next),
        "pos" => Ok(HalfInfiniteListType::Position),
        _ => Err(invalid("unknown infinite-list type")),
    }
}

fn values_from_python(values: Vec<Py<PyAny>>) -> Vec<PythonValue> {
    values.into_iter().map(PythonValue::new).collect()
}

fn values_to_python(py: Python<'_>, values: Vec<PythonValue>) -> Vec<Py<PyAny>> {
    values.iter().map(|value| value.to_python(py)).collect()
}

fn optional_position_to_python(
    py: Python<'_>,
    value: Option<ExtendedInteger>,
) -> PyResult<Py<PyAny>> {
    match value {
        None => Ok(pyo3::types::PyString::new(py, "none").into_any().unbind()),
        Some(value) => position_to_python(py, value),
    }
}

fn position_to_python(py: Python<'_>, value: ExtendedInteger) -> PyResult<Py<PyAny>> {
    match value {
        ExtendedInteger::NegativeInfinity => Ok(pyo3::types::PyString::new(py, "-infinity")
            .into_any()
            .unbind()),
        ExtendedInteger::Finite(value) => wrap_python_integer(py, &value),
        ExtendedInteger::PositiveInfinity => Ok(pyo3::types::PyString::new(py, "infinity")
            .into_any()
            .unbind()),
    }
}

fn next_function(
    callback: Py<PyAny>,
) -> InfiniteListNextFunction<PythonValue, PythonInfiniteListError> {
    Arc::new(move |value| {
        Python::attach(|py| {
            callback
                .call1(py, (value.to_python(py),))
                .map(PythonValue::new)
                .map_err(PythonInfiniteListError::Python)
        })
    })
}

fn position_function(
    callback: Py<PyAny>,
) -> InfiniteListPositionFunction<PythonValue, PythonInfiniteListError> {
    Arc::new(move |position| {
        Python::attach(|py| {
            let position =
                wrap_python_integer(py, position).map_err(PythonInfiniteListError::Python)?;
            callback
                .call1(py, (position,))
                .map(PythonValue::new)
                .map_err(PythonInfiniteListError::Python)
        })
    })
}

fn new_value_callback(
    callback: Option<Py<PyAny>>,
) -> Option<InfiniteListNewValueCallback<PythonInfiniteListError>> {
    callback.map(|callback| {
        Arc::new(move |position: &Integer, direction, list_type| {
            Python::attach(|py| {
                let position =
                    wrap_python_integer(py, position).map_err(PythonInfiniteListError::Python)?;
                callback
                    .call1(
                        py,
                        (position, direction_value(direction), type_name(list_type)),
                    )
                    .map(|_| ())
                    .map_err(PythonInfiniteListError::Python)
            })
        }) as InfiniteListNewValueCallback<PythonInfiniteListError>
    })
}

fn make_repeatify_callback(
    callback: Option<Py<PyAny>>,
) -> Option<InfiniteListRepeatifyCallback<PythonValue, PythonInfiniteListError>> {
    callback.map(|callback| {
        Arc::new(move |tail: Vec<PythonValue>, repeating: PythonHalfList| {
            Python::attach(|py| {
                let repeating = Py::new(py, PyHalfInfiniteList { inner: repeating })
                    .map_err(PythonInfiniteListError::Python)?;
                callback
                    .call1(py, (values_to_python(py, tail), repeating))
                    .map(|_| ())
                    .map_err(PythonInfiniteListError::Python)
            })
        }) as InfiniteListRepeatifyCallback<PythonValue, PythonInfiniteListError>
    })
}

fn specification_items<'py>(specification: &Bound<'py, PyAny>) -> PyResult<Vec<Bound<'py, PyAny>>> {
    specification.try_iter()?.collect()
}

fn parse_specification(
    specification: &Bound<'_, PyAny>,
) -> PyResult<HalfInfiniteListSpecification<PythonValue, PythonInfiniteListError>> {
    let items = specification_items(specification)?;
    let Some(kind) = items.first() else {
        return Err(invalid("an infinite-list specification cannot be empty"));
    };
    let kind = kind.extract::<String>()?;
    match kind.as_str() {
        "repeat" if items.len() == 2 => Ok(HalfInfiniteListSpecification::Repeat(
            values_from_python(items[1].extract::<Vec<Py<PyAny>>>()?),
        )),
        "next" | "next/repeat" if items.len() == 3 => {
            if !items[1].is_callable() {
                return Err(invalid("the next-list element function must be callable"));
            }
            Ok(HalfInfiniteListSpecification::Next {
                function: next_function(items[1].clone().unbind()),
                initial_value: PythonValue::new(items[2].clone().unbind()),
                repeatify: kind == "next/repeat",
            })
        }
        "pos" if (2..=3).contains(&items.len()) => {
            if !items[1].is_callable() {
                return Err(invalid(
                    "the position-list element function must be callable",
                ));
            }
            Ok(HalfInfiniteListSpecification::Position {
                function: position_function(items[1].clone().unbind()),
                store_values: if items.len() == 3 {
                    items[2].extract::<bool>()?
                } else {
                    true
                },
            })
        }
        _ => Err(invalid("invalid QPA infinite-list specification")),
    }
}

/// One callable returned by QPA's `ElementFunction` accessor.
#[pyclass(name = "InfiniteListElementFunction", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyInfiniteListElementFunction {
    inner: HalfInfiniteListElementFunction<PythonValue, PythonInfiniteListError>,
}

/// One callable returned by QPA's `NewValueCallback` accessor.
#[pyclass(name = "InfiniteListNewValueCallback", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyInfiniteListNewValueCallback {
    inner: InfiniteListNewValueCallback<PythonInfiniteListError>,
}

/// QPA's lazy half-infinite list.
#[pyclass(name = "HalfInfiniteList", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyHalfInfiniteList {
    inner: PythonHalfList,
}

/// QPA's two-sided infinite-list object.
#[pyclass(name = "InfiniteList", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyInfiniteList {
    inner: PythonFullList,
}

#[pymethods]
impl PyInfiniteListElementFunction {
    fn __call__(&self, py: Python<'_>, argument: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let value = match &self.inner {
            HalfInfiniteListElementFunction::Next(function) => {
                function(&PythonValue::new(argument.clone().unbind()))
            }
            HalfInfiniteListElementFunction::Position(function) => {
                function(&parse_python_integer(argument)?)
            }
        }
        .map_err(python_error)?;
        Ok(value.to_python(py))
    }
}

#[pymethods]
impl PyInfiniteListNewValueCallback {
    fn __call__(
        &self,
        position: &Bound<'_, PyAny>,
        direction: i8,
        list_type: &str,
    ) -> PyResult<()> {
        (self.inner)(
            &parse_python_integer(position)?,
            parse_direction(direction)?,
            parse_type_name(list_type)?,
        )
        .map_err(python_error)
    }
}

#[pymethods]
impl PyHalfInfiniteList {
    #[new]
    #[pyo3(signature = (start, direction, specification, callback=None, repeatify_callback=None))]
    fn new(
        start: &Bound<'_, PyAny>,
        direction: i8,
        specification: &Bound<'_, PyAny>,
        callback: Option<Py<PyAny>>,
        repeatify_callback: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: HalfInfiniteList::new(
                parse_python_integer(start)?,
                parse_direction(direction)?,
                parse_specification(specification)?,
                new_value_callback(callback),
                make_repeatify_callback(repeatify_callback),
            ),
        })
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_half_infinite_list(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn start_position(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        wrap_python_integer(py, self.inner.start_position())
    }

    #[getter]
    fn direction(&self) -> i8 {
        direction_value(self.inner.direction())
    }

    #[getter]
    fn list_type(&self) -> PyResult<&'static str> {
        self.inner.list_type().map(type_name).map_err(python_error)
    }

    #[getter]
    fn repeating_list(&self, py: Python<'_>) -> PyResult<Option<Vec<Py<PyAny>>>> {
        Ok(self
            .inner
            .repeating_list()
            .map_err(python_error)?
            .map(|values| values_to_python(py, values)))
    }

    #[getter]
    fn element_function(&self) -> PyResult<Option<PyInfiniteListElementFunction>> {
        Ok(self
            .inner
            .element_function()
            .map_err(python_error)?
            .map(|inner| PyInfiniteListElementFunction { inner }))
    }

    #[getter]
    fn is_storing_values(&self) -> PyResult<bool> {
        self.inner.is_storing_values().map_err(python_error)
    }

    #[getter]
    fn new_value_callback(&self) -> PyResult<Option<PyInfiniteListNewValueCallback>> {
        Ok(self
            .inner
            .new_value_callback()
            .map_err(python_error)?
            .map(|inner| PyInfiniteListNewValueCallback { inner }))
    }

    #[getter]
    fn is_repeating(&self) -> PyResult<bool> {
        self.inner.is_repeating().map_err(python_error)
    }

    #[getter]
    fn initial_value(&self, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        Ok(self
            .inner
            .initial_value()
            .map_err(python_error)?
            .map(|value| value.to_python(py)))
    }

    #[getter]
    fn stored_values(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        Ok(values_to_python(
            py,
            self.inner.stored_values().map_err(python_error)?,
        ))
    }

    fn value(&self, py: Python<'_>, position: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        Ok(self
            .inner
            .value(&parse_python_integer(position)?)
            .map_err(python_error)?
            .to_python(py))
    }

    fn __getitem__(&self, py: Python<'_>, position: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        self.value(py, position)
    }

    fn make_repeating_list(
        &self,
        py: Python<'_>,
        first: usize,
        second: usize,
    ) -> PyResult<(Vec<Py<PyAny>>, Self)> {
        let (tail, inner) = self
            .inner
            .make_repeating_list(first, second)
            .map_err(python_error)?;
        Ok((values_to_python(py, tail), Self { inner }))
    }

    #[getter]
    fn lowest_known_position(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        optional_position_to_python(
            py,
            self.inner.lowest_known_position().map_err(python_error)?,
        )
    }

    #[getter]
    fn highest_known_position(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        optional_position_to_python(
            py,
            self.inner.highest_known_position().map_err(python_error)?,
        )
    }

    fn shift(&self, shift: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .shift(&parse_python_integer(shift)?)
                .map_err(python_error)?,
        })
    }

    fn cut(&self, position: &Bound<'_, PyAny>) -> PyResult<PyInfiniteList> {
        Ok(PyInfiniteList {
            inner: self
                .inner
                .cut(&parse_python_integer(position)?)
                .map_err(python_error)?,
        })
    }

    fn map(&self, py: Python<'_>, function: Py<PyAny>) -> PyResult<Self> {
        if !function.bind(py).is_callable() {
            return Err(invalid("infinite-list mapping function must be callable"));
        }
        Ok(Self {
            inner: self
                .inner
                .map(next_function(function))
                .map_err(python_error)?,
        })
    }

    fn is_identical(&self, other: &Self) -> bool {
        self.inner.is_identical(&other.inner)
    }
}

#[pymethods]
impl PyInfiniteList {
    #[new]
    #[pyo3(signature = (base_position, middle, positive, negative, callback=None))]
    fn new(
        base_position: &Bound<'_, PyAny>,
        middle: Vec<Py<PyAny>>,
        positive: &Bound<'_, PyAny>,
        negative: &Bound<'_, PyAny>,
        callback: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: InfiniteList::new(
                parse_python_integer(base_position)?,
                values_from_python(middle),
                parse_specification(positive)?,
                parse_specification(negative)?,
                new_value_callback(callback),
            ),
        })
    }

    #[staticmethod]
    fn from_half_lists(
        base_position: &Bound<'_, PyAny>,
        middle: Vec<Py<PyAny>>,
        positive: Option<&PyHalfInfiniteList>,
        negative: Option<&PyHalfInfiniteList>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: InfiniteList::from_parts(
                parse_python_integer(base_position)?,
                values_from_python(middle),
                positive.map(|value| value.inner.clone()),
                negative.map(|value| value.inner.clone()),
            )
            .map_err(python_error)?,
        })
    }

    #[staticmethod]
    fn function(function: Py<PyAny>) -> PyResult<Self> {
        if !Python::attach(|py| function.bind(py).is_callable()) {
            return Err(invalid("infinite-list element function must be callable"));
        }
        let function = position_function(function);
        Ok(Self {
            inner: InfiniteList::function(move |position| function(position)),
        })
    }

    #[staticmethod]
    fn constant(value: Py<PyAny>) -> Self {
        Self {
            inner: InfiniteList::constant(PythonValue::new(value)),
        }
    }

    #[staticmethod]
    fn finite(base_position: &Bound<'_, PyAny>, values: Vec<Py<PyAny>>) -> PyResult<Self> {
        Ok(Self {
            inner: InfiniteList::finite(
                parse_python_integer(base_position)?,
                values_from_python(values),
            ),
        })
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_infinite_list(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn middle_start(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        wrap_python_integer(py, self.inner.middle_start())
    }

    #[getter]
    fn middle_end(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        wrap_python_integer(py, &self.inner.middle_end().map_err(python_error)?)
    }

    #[getter]
    fn middle_part(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        Ok(values_to_python(
            py,
            self.inner.middle_part().map_err(python_error)?,
        ))
    }

    #[getter]
    fn positive_part(&self) -> PyResult<Option<PyHalfInfiniteList>> {
        Ok(self
            .inner
            .positive_part()
            .map_err(python_error)?
            .map(|inner| PyHalfInfiniteList { inner }))
    }

    #[getter]
    fn negative_part(&self) -> PyResult<Option<PyHalfInfiniteList>> {
        Ok(self
            .inner
            .negative_part()
            .map_err(python_error)?
            .map(|inner| PyHalfInfiniteList { inner }))
    }

    #[getter]
    fn lowest_known_position(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        optional_position_to_python(
            py,
            self.inner.lowest_known_position().map_err(python_error)?,
        )
    }

    #[getter]
    fn highest_known_position(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        optional_position_to_python(
            py,
            self.inner.highest_known_position().map_err(python_error)?,
        )
    }

    #[getter]
    fn upper_bound(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        position_to_python(py, self.inner.upper_bound().map_err(python_error)?)
    }

    #[getter]
    fn lower_bound(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        position_to_python(py, self.inner.lower_bound().map_err(python_error)?)
    }

    fn value(&self, py: Python<'_>, position: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        Ok(self
            .inner
            .value(&parse_python_integer(position)?)
            .map_err(python_error)?
            .to_python(py))
    }

    fn __getitem__(&self, py: Python<'_>, position: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        self.value(py, position)
    }

    fn shift(&self, shift: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .shift(&parse_python_integer(shift)?)
                .map_err(python_error)?,
        })
    }

    fn finite_part_as_list(
        &self,
        py: Python<'_>,
        start: &Bound<'_, PyAny>,
        end: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Py<PyAny>>> {
        Ok(values_to_python(
            py,
            self.inner
                .finite_part_as_list(&parse_python_integer(start)?, &parse_python_integer(end)?)
                .map_err(python_error)?,
        ))
    }

    fn positive_part_from(&self, position: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .positive_part_from(&parse_python_integer(position)?)
                .map_err(python_error)?,
        })
    }

    fn negative_part_from(&self, position: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .negative_part_from(&parse_python_integer(position)?)
                .map_err(python_error)?,
        })
    }

    #[staticmethod]
    fn splice(
        positive: &Self,
        negative: &Self,
        join_position: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: InfiniteList::splice(
                &positive.inner,
                &negative.inner,
                &parse_python_integer(join_position)?,
            )
            .map_err(python_error)?,
        })
    }

    #[staticmethod]
    fn shifted_splice(
        left: &Self,
        left_shift: &Bound<'_, PyAny>,
        right: &Self,
        right_shift: &Bound<'_, PyAny>,
        join_position: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: shifted_splice(
                &left.inner,
                &parse_python_integer(left_shift)?,
                &right.inner,
                &parse_python_integer(right_shift)?,
                &parse_python_integer(join_position)?,
            )
            .map_err(python_error)?,
        })
    }

    #[staticmethod]
    fn concatenate(lists: Vec<Py<Self>>, py: Python<'_>) -> PyResult<Self> {
        let lists = lists
            .iter()
            .map(|list| list.borrow(py).inner.clone())
            .collect::<Vec<_>>();
        Ok(Self {
            inner: InfiniteList::concatenate(&lists).map_err(python_error)?,
        })
    }

    fn map(&self, py: Python<'_>, function: Py<PyAny>) -> PyResult<Self> {
        if !function.bind(py).is_callable() {
            return Err(invalid("infinite-list mapping function must be callable"));
        }
        let function = next_function(function);
        Ok(Self {
            inner: self
                .inner
                .map(move |value| function(value))
                .map_err(python_error)?,
        })
    }

    fn is_identical(&self, other: &Self) -> bool {
        self.inner.is_identical(&other.inner)
    }
}

/// Return QPA's process-global integer infinite list.
#[pyfunction(name = "integers_list")]
pub(crate) fn py_integers_list() -> PyInfiniteList {
    PyInfiniteList {
        inner: InfiniteList::function(|position| {
            Python::attach(|py| {
                wrap_python_integer(py, position)
                    .map(PythonValue::new)
                    .map_err(PythonInfiniteListError::Python)
            })
        }),
    }
}
