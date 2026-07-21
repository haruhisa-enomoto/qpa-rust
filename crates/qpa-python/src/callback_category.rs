//! Python adapter for QPA's generic `Cat(identity, properties)` constructor.

use std::{fmt, sync::Arc};

use pyo3::{
    Bound, Py, PyAny, PyResult, Python,
    prelude::{pyclass, pyfunction, pymethods},
    types::{PyAnyMethods, PyDict, PyDictMethods},
};
use qpa_gbnp::{CallbackCategory, CallbackCategoryProperties, ComplexError, QpaCategory};

use crate::invalid;

#[derive(Clone)]
struct PythonCategoryValue(Arc<Py<PyAny>>);

impl PythonCategoryValue {
    fn new(value: Py<PyAny>) -> Self {
        Self(Arc::new(value))
    }

    fn to_python(&self, py: Python<'_>) -> Py<PyAny> {
        self.0.as_ref().clone_ref(py)
    }
}

impl fmt::Debug for PythonCategoryValue {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let text = Python::attach(|py| {
            self.0
                .bind(py)
                .repr()
                .and_then(|value| value.extract::<String>())
                .unwrap_or_else(|_| "<python object>".to_owned())
        });
        formatter.write_str(&text)
    }
}

impl PartialEq for PythonCategoryValue {
    fn eq(&self, other: &Self) -> bool {
        Python::attach(|py| self.0.bind(py).eq(other.0.bind(py)).unwrap_or(false))
    }
}

impl Eq for PythonCategoryValue {}

type PythonCallbackCategory =
    CallbackCategory<PythonCategoryValue, PythonCategoryValue, PythonCategoryValue>;

fn callback_error(error: pyo3::PyErr) -> ComplexError {
    ComplexError::CategoryCallback(error.to_string())
}

fn call_one(
    callback: &Py<PyAny>,
    value: &PythonCategoryValue,
) -> Result<PythonCategoryValue, ComplexError> {
    Python::attach(|py| {
        callback
            .call1(py, (value.to_python(py),))
            .map(PythonCategoryValue::new)
            .map_err(callback_error)
    })
}

fn call_two(
    callback: &Py<PyAny>,
    left: &PythonCategoryValue,
    right: &PythonCategoryValue,
) -> Result<PythonCategoryValue, ComplexError> {
    Python::attach(|py| {
        callback
            .call1(py, (left.to_python(py), right.to_python(py)))
            .map(PythonCategoryValue::new)
            .map_err(callback_error)
    })
}

fn call_bool_one(callback: &Py<PyAny>, value: &PythonCategoryValue) -> bool {
    Python::attach(|py| {
        callback
            .call1(py, (value.to_python(py),))
            .and_then(|result| result.extract::<bool>(py))
            .unwrap_or(false)
    })
}

fn call_bool_two(
    callback: &Py<PyAny>,
    left: &PythonCategoryValue,
    right: &PythonCategoryValue,
) -> Result<bool, ComplexError> {
    Python::attach(|py| {
        callback
            .call1(py, (left.to_python(py), right.to_python(py)))
            .and_then(|result| result.extract::<bool>(py))
            .map_err(callback_error)
    })
}

fn required(properties: &Bound<'_, PyDict>, name: &str) -> PyResult<Py<PyAny>> {
    let value = properties
        .get_item(name)?
        .ok_or_else(|| invalid(format!("Cat properties record has no {name} entry")))?;
    Ok(value.unbind())
}

fn required_callable(properties: &Bound<'_, PyDict>, name: &str) -> PyResult<Py<PyAny>> {
    let value = required(properties, name)?;
    if value.bind(properties.py()).is_callable() {
        Ok(value)
    } else {
        Err(invalid(format!("Cat property {name} must be callable")))
    }
}

fn optional_callable(properties: &Bound<'_, PyDict>, name: &str) -> PyResult<Option<Py<PyAny>>> {
    let Some(value) = properties.get_item(name)? else {
        return Ok(None);
    };
    if value.is_callable() {
        Ok(Some(value.unbind()))
    } else {
        Err(invalid(format!("Cat property {name} must be callable")))
    }
}

fn object_member(
    value: &PythonCategoryValue,
    callback: Option<&Py<PyAny>>,
    member: &str,
) -> PythonCategoryValue {
    Python::attach(|py| {
        let result: PyResult<Py<PyAny>> = (|| {
            if let Some(callback) = callback {
                callback.call1(py, (value.to_python(py),))
            } else {
                let member = value.0.bind(py).getattr(member)?;
                if member.is_callable() {
                    member.call0().map(Bound::unbind)
                } else {
                    Ok(member.unbind())
                }
            }
        })();
        result.map_or_else(|_| value.clone(), PythonCategoryValue::new)
    })
}

fn negated(value: &PythonCategoryValue, callback: Option<&Py<PyAny>>) -> PythonCategoryValue {
    Python::attach(|py| {
        let result = if let Some(callback) = callback {
            callback.call1(py, (value.to_python(py),))
        } else {
            value.0.bind(py).call_method0("__neg__").map(Bound::unbind)
        };
        result.map_or_else(|_| value.clone(), PythonCategoryValue::new)
    })
}

/// QPA's generic category record backed by Python callbacks.
#[pyclass(name = "CallbackCategory", frozen, skip_from_py_object)]
pub(crate) struct PyCallbackCategory {
    inner: PythonCallbackCategory,
    properties: Arc<Py<PyDict>>,
}

#[pymethods]
impl PyCallbackCategory {
    #[getter]
    #[allow(clippy::unused_self)]
    fn is_cat(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn identity(&self, py: Python<'_>) -> Py<PyAny> {
        self.inner.identity().to_python(py)
    }

    #[getter]
    fn name(&self) -> String {
        self.inner.name().to_owned()
    }

    fn property(&self, py: Python<'_>, name: &str) -> PyResult<Py<PyAny>> {
        self.properties
            .bind(py)
            .get_item(name)?
            .map(Bound::unbind)
            .ok_or_else(|| invalid(format!("Cat properties record has no {name} entry")))
    }

    fn zero_object(&self, py: Python<'_>) -> Py<PyAny> {
        self.inner.zero_object().to_python(py)
    }

    fn is_zero_object(&self, value: Py<PyAny>) -> bool {
        self.inner.is_zero_object(&PythonCategoryValue::new(value))
    }

    fn zero_map(&self, py: Python<'_>, source: Py<PyAny>, range: Py<PyAny>) -> PyResult<Py<PyAny>> {
        self.inner
            .zero_map(
                &PythonCategoryValue::new(source),
                &PythonCategoryValue::new(range),
            )
            .map(|value| value.to_python(py))
            .map_err(invalid)
    }

    fn is_zero_mapping(&self, value: Py<PyAny>) -> bool {
        self.inner.is_zero_mapping(&PythonCategoryValue::new(value))
    }

    fn compose_maps(
        &self,
        py: Python<'_>,
        left: Py<PyAny>,
        right: Py<PyAny>,
    ) -> PyResult<Py<PyAny>> {
        self.inner
            .compose_maps(
                &PythonCategoryValue::new(left),
                &PythonCategoryValue::new(right),
            )
            .map(|value| value.to_python(py))
            .map_err(invalid)
    }

    fn kernel(&self, py: Python<'_>, value: Py<PyAny>) -> PyResult<Py<PyAny>> {
        self.inner
            .kernel(&PythonCategoryValue::new(value))
            .map(|value| value.to_python(py))
            .map_err(invalid)
    }

    fn image(&self, py: Python<'_>, value: Py<PyAny>) -> PyResult<Py<PyAny>> {
        self.inner
            .image(&PythonCategoryValue::new(value))
            .map(|value| value.to_python(py))
            .map_err(invalid)
    }

    fn is_exact(&self, left: Py<PyAny>, right: Py<PyAny>) -> PyResult<bool> {
        self.inner
            .is_exact(
                &PythonCategoryValue::new(left),
                &PythonCategoryValue::new(right),
            )
            .map_err(invalid)
    }

    fn object_string(&self, value: Py<PyAny>) -> String {
        self.inner.object_string(&PythonCategoryValue::new(value))
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __repr__(&self) -> String {
        format!("<cat: {}>", self.inner.name())
    }

    fn print_obj(&self) -> String {
        self.__repr__()
    }

    fn view(&self) -> String {
        self.__repr__()
    }
}

/// Construct QPA's generic `Cat(identity, properties)` callback record.
#[pyfunction(name = "cat")]
pub(crate) fn py_cat(
    identity: Py<PyAny>,
    properties: &Bound<'_, PyDict>,
) -> PyResult<PyCallbackCategory> {
    let name = required(properties, "name")?
        .extract::<String>(properties.py())?
        .into();
    let zero_object = PythonCategoryValue::new(required(properties, "zeroObj")?);
    let is_zero_object = Arc::new(required_callable(properties, "isZeroObj")?);
    let zero_map = Arc::new(required_callable(properties, "zeroMap")?);
    let is_zero_mapping = Arc::new(required_callable(properties, "isZeroMapping")?);
    let compose_maps = Arc::new(required_callable(properties, "composeMaps")?);
    let kernel = Arc::new(required_callable(properties, "ker")?);
    let image = Arc::new(required_callable(properties, "im")?);
    let is_exact = Arc::new(required_callable(properties, "isExact")?);
    let object_string = Arc::new(required_callable(properties, "objStr")?);
    let source = optional_callable(properties, "source")?.map(Arc::new);
    let range = optional_callable(properties, "range")?.map(Arc::new);
    let negate = optional_callable(properties, "negate")?.map(Arc::new);

    let inner = CallbackCategory::new(
        PythonCategoryValue::new(identity),
        CallbackCategoryProperties {
            name,
            zero_object,
            is_zero_object: Arc::new(move |value| call_bool_one(&is_zero_object, value)),
            zero_map: Arc::new(move |left, right| call_two(&zero_map, left, right)),
            is_zero_mapping: Arc::new(move |value| call_bool_one(&is_zero_mapping, value)),
            compose_maps: Arc::new(move |left, right| call_two(&compose_maps, left, right)),
            source: Arc::new(move |value| object_member(value, source.as_deref(), "source")),
            range: Arc::new(move |value| object_member(value, range.as_deref(), "range")),
            negate: Arc::new(move |value| Ok(negated(value, negate.as_deref()))),
            kernel: Arc::new(move |value| call_one(&kernel, value)),
            image: Arc::new(move |value| call_one(&image, value)),
            is_exact: Arc::new(move |left, right| call_bool_two(&is_exact, left, right)),
            object_string: Arc::new(move |value| {
                Python::attach(|py| {
                    object_string
                        .call1(py, (value.to_python(py),))
                        .and_then(|result| result.extract::<String>(py))
                        .unwrap_or_else(|_| format!("{value:?}"))
                })
            }),
        },
    );
    Ok(PyCallbackCategory {
        inner,
        properties: Arc::new(properties.clone().unbind()),
    })
}
