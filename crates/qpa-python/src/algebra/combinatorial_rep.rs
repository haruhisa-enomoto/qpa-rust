//! Python wrappers for QPA's integral unit and Euler forms.

use std::sync::Arc;

use pyo3::{PyAny, pymethods};
use qpa_core::UnitForm;
use qpa_gbnp::EulerBilinearForm;

#[allow(clippy::wildcard_imports)]
// The algebra and coefficient runtime unions belong to the parent module.
use super::*;
use crate::field::{
    parse_python_integer_matrix, parse_python_integer_vector, wrap_python_integer,
    wrap_python_integer_matrix, wrap_python_integer_vector,
};

fn scalar_values<F: runtime::PythonField>(
    py: Python<'_>,
    field: &F,
    values: Vec<Py<PyScalar>>,
) -> PyResult<Vec<F::Element>> {
    values
        .into_iter()
        .map(|value| field.scalar_value(&value.borrow(py)))
        .collect()
}

fn unit_bilinear_over<F: runtime::PythonField>(
    py: Python<'_>,
    form: &UnitForm,
    field: &F,
    left: Vec<Py<PyScalar>>,
    right: Vec<Py<PyScalar>>,
) -> PyResult<PyScalar> {
    let left = scalar_values(py, field, left)?;
    let right = scalar_values(py, field, right)?;
    Ok(field.wrap_scalar(&form.bilinear_form(field, &left, &right).map_err(invalid)?))
}

fn unit_quadratic_over<F: runtime::PythonField>(
    py: Python<'_>,
    form: &UnitForm,
    field: &F,
    vector: Vec<Py<PyScalar>>,
) -> PyResult<PyScalar> {
    let vector = scalar_values(py, field, vector)?;
    Ok(field.wrap_scalar(&form.quadratic_form(field, &vector).map_err(invalid)?))
}

fn unit_reflection_over<F: runtime::PythonField>(
    py: Python<'_>,
    form: &UnitForm,
    field: &F,
    coordinate: usize,
    vector: Vec<Py<PyScalar>>,
) -> PyResult<Vec<PyScalar>> {
    let vector = scalar_values(py, field, vector)?;
    Ok(form
        .reflection_by_bilinear_form(field, coordinate, &vector)
        .map_err(invalid)?
        .iter()
        .map(|value| field.wrap_scalar(value))
        .collect())
}

fn euler_evaluate_over<F: runtime::PythonField>(
    py: Python<'_>,
    form: &EulerBilinearForm,
    field: &F,
    left: Vec<Py<PyScalar>>,
    right: Vec<Py<PyScalar>>,
) -> PyResult<PyScalar> {
    let left = scalar_values(py, field, left)?;
    let right = scalar_values(py, field, right)?;
    Ok(field.wrap_scalar(&form.evaluate(field, &left, &right).map_err(invalid)?))
}

/// QPA's integral symmetric unit form and its cached positivity attributes.
#[pyclass(name = "UnitForm", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyUnitForm {
    inner: Arc<UnitForm>,
}

impl PyUnitForm {
    fn from_inner(inner: UnitForm) -> Self {
        Self {
            inner: Arc::new(inner),
        }
    }
}

#[pymethods]
impl PyUnitForm {
    #[new]
    fn new(py: Python<'_>, rows: Vec<Vec<Py<PyAny>>>) -> PyResult<Self> {
        Ok(Self::from_inner(
            UnitForm::new(parse_python_integer_matrix(py, rows)?).map_err(invalid)?,
        ))
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_unit_form(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn symmetric_matrix(&self, py: Python<'_>) -> PyResult<Vec<Vec<Py<PyAny>>>> {
        wrap_python_integer_matrix(
            py,
            self.inner.symmetric_matrix().row_iter().map(<[_]>::to_vec),
        )
    }

    #[getter]
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }

    fn bilinear_form(
        &self,
        py: Python<'_>,
        left: Vec<Py<PyAny>>,
        right: Vec<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let left = parse_python_integer_vector(py, left)?;
        let right = parse_python_integer_vector(py, right)?;
        let value = self
            .inner
            .bilinear_form(&numerica::domains::integer::Z, &left, &right)
            .map_err(invalid)?;
        wrap_python_integer(py, &value)
    }

    fn bilinear_form_over(
        &self,
        py: Python<'_>,
        field: &PyField,
        left: Vec<Py<PyScalar>>,
        right: Vec<Py<PyScalar>>,
    ) -> PyResult<PyScalar> {
        match field.kind {
            FieldKind::Rational => unit_bilinear_over(py, &self.inner, &Q, left, right),
            FieldKind::Binary => unit_bilinear_over(py, &self.inner, &Z2, left, right),
            FieldKind::Prime(characteristic) => {
                unit_bilinear_over(py, &self.inner, &Zp::new(characteristic), left, right)
            }
        }
    }

    fn quadratic_form(&self, py: Python<'_>, vector: Vec<Py<PyAny>>) -> PyResult<Py<PyAny>> {
        let vector = parse_python_integer_vector(py, vector)?;
        let value = self
            .inner
            .quadratic_form(&numerica::domains::integer::Z, &vector)
            .map_err(invalid)?;
        wrap_python_integer(py, &value)
    }

    fn quadratic_form_over(
        &self,
        py: Python<'_>,
        field: &PyField,
        vector: Vec<Py<PyScalar>>,
    ) -> PyResult<PyScalar> {
        match field.kind {
            FieldKind::Rational => unit_quadratic_over(py, &self.inner, &Q, vector),
            FieldKind::Binary => unit_quadratic_over(py, &self.inner, &Z2, vector),
            FieldKind::Prime(characteristic) => {
                unit_quadratic_over(py, &self.inner, &Zp::new(characteristic), vector)
            }
        }
    }

    fn reflection_by_bilinear_form(
        &self,
        py: Python<'_>,
        coordinate: usize,
        vector: Vec<Py<PyAny>>,
    ) -> PyResult<Vec<Py<PyAny>>> {
        let vector = parse_python_integer_vector(py, vector)?;
        let reflection = self
            .inner
            .reflection_by_bilinear_form(&numerica::domains::integer::Z, coordinate, &vector)
            .map_err(invalid)?;
        wrap_python_integer_vector(py, &reflection)
    }

    fn reflection_by_bilinear_form_over(
        &self,
        py: Python<'_>,
        field: &PyField,
        coordinate: usize,
        vector: Vec<Py<PyScalar>>,
    ) -> PyResult<Vec<PyScalar>> {
        match field.kind {
            FieldKind::Rational => unit_reflection_over(py, &self.inner, &Q, coordinate, vector),
            FieldKind::Binary => unit_reflection_over(py, &self.inner, &Z2, coordinate, vector),
            FieldKind::Prime(characteristic) => unit_reflection_over(
                py,
                &self.inner,
                &Zp::new(characteristic),
                coordinate,
                vector,
            ),
        }
    }

    #[getter]
    fn is_weakly_nonnegative(&self) -> bool {
        self.inner.is_weakly_nonnegative()
    }

    #[getter]
    fn is_weakly_positive(&self) -> bool {
        self.inner.is_weakly_positive()
    }

    #[getter]
    fn positive_roots(&self, py: Python<'_>) -> PyResult<Option<Vec<Vec<Py<PyAny>>>>> {
        self.inner
            .positive_roots()
            .map(|roots| wrap_python_integer_matrix(py, roots.iter().cloned()))
            .transpose()
    }

    fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.inner, &other.inner)
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    fn view(&self) -> String {
        self.__repr__()
    }
}

/// QPA's captured Euler bilinear closure with its inverse Cartan matrix.
#[pyclass(name = "EulerBilinearForm", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyEulerBilinearForm {
    inner: EulerBilinearForm,
}

#[pymethods]
impl PyEulerBilinearForm {
    #[getter]
    fn inverse_cartan_matrix(&self, py: Python<'_>) -> PyResult<Vec<Vec<Py<PyAny>>>> {
        wrap_python_integer_matrix(
            py,
            self.inner
                .inverse_cartan_matrix()
                .row_iter()
                .map(<[_]>::to_vec),
        )
    }

    #[getter]
    fn dimension(&self) -> usize {
        self.inner.dimension()
    }

    fn evaluate(
        &self,
        py: Python<'_>,
        left: Vec<Py<PyAny>>,
        right: Vec<Py<PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        let left = parse_python_integer_vector(py, left)?;
        let right = parse_python_integer_vector(py, right)?;
        let value = self
            .inner
            .evaluate(&numerica::domains::integer::Z, &left, &right)
            .map_err(invalid)?;
        wrap_python_integer(py, &value)
    }

    fn evaluate_over(
        &self,
        py: Python<'_>,
        field: &PyField,
        left: Vec<Py<PyScalar>>,
        right: Vec<Py<PyScalar>>,
    ) -> PyResult<PyScalar> {
        match field.kind {
            FieldKind::Rational => euler_evaluate_over(py, &self.inner, &Q, left, right),
            FieldKind::Binary => euler_evaluate_over(py, &self.inner, &Z2, left, right),
            FieldKind::Prime(characteristic) => {
                euler_evaluate_over(py, &self.inner, &Zp::new(characteristic), left, right)
            }
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "<Euler bilinear form of dimension {}>",
            self.inner.dimension()
        )
    }
}

pub(super) fn tits_unit_form_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<PyUnitForm> {
    Ok(PyUnitForm::from_inner(
        qpa_gbnp::tits_unit_form_of_algebra(algebra).map_err(invalid)?,
    ))
}

pub(super) fn euler_bilinear_form_for<F: runtime::PythonField>(
    algebra: &QuiverAlgebra<F>,
) -> PyResult<Option<PyEulerBilinearForm>> {
    Ok(qpa_gbnp::euler_bilinear_form_of_algebra(algebra)
        .map_err(invalid)?
        .map(|inner| PyEulerBilinearForm { inner }))
}
