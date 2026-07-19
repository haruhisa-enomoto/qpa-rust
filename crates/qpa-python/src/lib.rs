//! Small `PyO3` architecture spike over the already ported QPA Rust surface.

#![allow(clippy::needless_pass_by_value)]
#![allow(missing_docs)] // PyO3 exception structs cannot receive rustdoc attributes.

use std::str::FromStr;

use numerica::domains::{
    Field as _, Ring as _, Set,
    finite_field::{FiniteFieldCore, PrimeIteratorU64, ToFiniteField, Z2, Zp},
    integer::Integer,
    rational::{Q, Rational, RationalField},
};
use pyo3::{
    Bound, Py, PyErr, PyRef, PyResult, Python, create_exception,
    exceptions::{PyException, PyTypeError},
    prelude::*,
    types::{PyAny, PyInt, PyModule},
};
use qpa_core::{
    AlgebraDimension, DenseMatrix, Generator, NamedArrowMatrix, NamedArrowSpec, Path, PathAlgebra,
    PathAlgebraElement, PathAlgebraModule, PathEndpoint, Quiver, QuiverAlgebra,
    QuotientPathAlgebra, QuotientPathAlgebraElement, right_module_over_path_algebra,
};

// Base exception for failures reported by the Rust QPA binding.
create_exception!(_qpa, QpaError, PyException);
// Invalid mathematical input, field ownership, or object-family exception.
create_exception!(_qpa, InvalidInputError, QpaError);
// QPA's explicit `fail` result for actions that do not satisfy relations.
create_exception!(_qpa, RelationNotSatisfiedError, QpaError);

type BinaryElement = <Z2 as Set>::Element;
type PrimeElement = <Zp as Set>::Element;

fn invalid(error: impl std::fmt::Display) -> PyErr {
    InvalidInputError::new_err(error.to_string())
}

fn parse_python_integer(value: &Bound<'_, PyAny>) -> PyResult<Integer> {
    if !value.is_instance_of::<PyInt>() {
        return Err(PyTypeError::new_err("an exact Python integer is required"));
    }
    Integer::from_str(value.str()?.to_str()?)
        .map_err(|error| PyTypeError::new_err(format!("invalid Python integer: {error}")))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FieldKind {
    Rational,
    Binary,
    Prime(u32),
}

impl FieldKind {
    fn display(self) -> String {
        match self {
            Self::Rational => "Rationals".into(),
            Self::Binary => "GF(2)".into(),
            Self::Prime(characteristic) => format!("GF({characteristic})"),
        }
    }

    fn characteristic(self) -> Option<u32> {
        match self {
            Self::Rational => None,
            Self::Binary => Some(2),
            Self::Prime(characteristic) => Some(characteristic),
        }
    }
}

/// Runtime exact field descriptor used by the Python architecture spike.
#[pyclass(name = "Field", frozen, skip_from_py_object)]
#[derive(Clone, Debug, Eq, PartialEq)]
struct PyField {
    kind: FieldKind,
}

#[pymethods]
impl PyField {
    #[staticmethod]
    fn rational() -> Self {
        Self {
            kind: FieldKind::Rational,
        }
    }

    #[staticmethod]
    fn prime(characteristic: u32) -> PyResult<Self> {
        let kind = match characteristic {
            2 => FieldKind::Binary,
            3.. => {
                if PrimeIteratorU64::new(u64::from(characteristic) - 1).next()
                    != Some(u64::from(characteristic))
                {
                    return Err(InvalidInputError::new_err(format!(
                        "{characteristic} is not prime"
                    )));
                }
                FieldKind::Prime(characteristic)
            }
            _ => {
                return Err(InvalidInputError::new_err(format!(
                    "{characteristic} is not prime"
                )));
            }
        };
        Ok(Self { kind })
    }

    #[getter]
    fn characteristic(&self) -> Option<u32> {
        self.kind.characteristic()
    }

    #[getter]
    fn is_finite(&self) -> bool {
        !matches!(self.kind, FieldKind::Rational)
    }

    #[pyo3(signature = (numerator, denominator=None))]
    fn element(
        &self,
        numerator: &Bound<'_, PyAny>,
        denominator: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<PyScalar> {
        let numerator = parse_python_integer(numerator)?;
        let denominator = denominator
            .map(parse_python_integer)
            .transpose()?
            .unwrap_or_else(Integer::one);
        if denominator.is_zero() {
            return Err(InvalidInputError::new_err(
                "coefficient denominator is zero",
            ));
        }
        match self.kind {
            FieldKind::Rational => Ok(PyScalar {
                value: ScalarKind::Rational(Q.to_element(numerator, denominator, true)),
            }),
            FieldKind::Binary => {
                let numerator = numerator.to_finite_field(&Z2);
                let denominator = denominator.to_finite_field(&Z2);
                if Z2.is_zero(&denominator) {
                    return Err(InvalidInputError::new_err(
                        "coefficient denominator is zero in GF(2)",
                    ));
                }
                let value = Z2.div(&numerator, &denominator);
                Ok(PyScalar::prime(2, u32::from(!Z2.is_zero(&value))))
            }
            FieldKind::Prime(characteristic) => {
                let field = Zp::new(characteristic);
                let numerator = numerator.to_finite_field(&field);
                let denominator = denominator.to_finite_field(&field);
                if field.is_zero(&denominator) {
                    return Err(InvalidInputError::new_err(format!(
                        "coefficient denominator is zero in GF({characteristic})"
                    )));
                }
                let value = field.div(&numerator, &denominator);
                Ok(PyScalar::prime(characteristic, field.from_element(&value)))
            }
        }
    }

    fn __repr__(&self) -> String {
        self.kind.display()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum ScalarKind {
    Rational(Rational),
    Prime {
        characteristic: u32,
        representative: u32,
    },
}

/// Exact scalar transferred between Python and a Numerica-backed Rust field.
#[pyclass(name = "Scalar", frozen, skip_from_py_object)]
#[derive(Clone, Debug, Eq, PartialEq)]
struct PyScalar {
    value: ScalarKind,
}

impl PyScalar {
    fn prime(characteristic: u32, representative: u32) -> Self {
        Self {
            value: ScalarKind::Prime {
                characteristic,
                representative,
            },
        }
    }

    fn rational_value(&self) -> PyResult<Rational> {
        match &self.value {
            ScalarKind::Rational(value) => Ok(value.clone()),
            ScalarKind::Prime { .. } => Err(invalid("coefficient belongs to a prime field")),
        }
    }

    fn binary_value(&self) -> PyResult<BinaryElement> {
        match self.value {
            ScalarKind::Prime {
                characteristic: 2,
                representative,
            } => Ok(representative.to_finite_field(&Z2)),
            _ => Err(invalid("coefficient does not belong to GF(2)")),
        }
    }

    fn prime_value(&self, field: &Zp) -> PyResult<PrimeElement> {
        let expected = field.get_prime();
        match self.value {
            ScalarKind::Prime {
                characteristic,
                representative,
            } if characteristic == expected => Ok(field.to_element(representative)),
            _ => Err(invalid(format!(
                "coefficient does not belong to GF({expected})"
            ))),
        }
    }
}

#[pymethods]
impl PyScalar {
    #[getter]
    fn field(&self) -> PyField {
        let kind = match self.value {
            ScalarKind::Rational(_) => FieldKind::Rational,
            ScalarKind::Prime {
                characteristic: 2, ..
            } => FieldKind::Binary,
            ScalarKind::Prime { characteristic, .. } => FieldKind::Prime(characteristic),
        };
        PyField { kind }
    }

    #[getter]
    fn numerator(&self) -> Option<String> {
        match &self.value {
            ScalarKind::Rational(value) => Some(value.numerator_ref().to_string()),
            ScalarKind::Prime { .. } => None,
        }
    }

    #[getter]
    fn denominator(&self) -> Option<String> {
        match &self.value {
            ScalarKind::Rational(value) => Some(value.denominator_ref().to_string()),
            ScalarKind::Prime { .. } => None,
        }
    }

    #[getter]
    fn representative(&self) -> Option<u32> {
        match self.value {
            ScalarKind::Rational(_) => None,
            ScalarKind::Prime { representative, .. } => Some(representative),
        }
    }

    fn __repr__(&self) -> String {
        match &self.value {
            ScalarKind::Rational(value) => value.to_string(),
            ScalarKind::Prime {
                characteristic,
                representative,
            } => format!("{representative} mod {characteristic}"),
        }
    }

    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }
}

/// Python-owned handle to a Rust QPA quiver identity.
#[pyclass(name = "Quiver", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyQuiver {
    inner: Quiver,
}

#[pymethods]
impl PyQuiver {
    #[new]
    fn new(vertices: Vec<String>, arrows: Vec<(String, String, String)>) -> PyResult<Self> {
        let arrows = arrows
            .into_iter()
            .map(|(source, target, name)| NamedArrowSpec::named(source, target, name))
            .collect();
        let inner = Quiver::from_named(
            vertices.into_iter().map(String::into_boxed_str).collect(),
            arrows,
        )
        .map_err(invalid)?;
        Ok(Self { inner })
    }

    #[getter]
    fn vertices(&self) -> Vec<String> {
        self.inner
            .vertices()
            .map(|vertex| vertex.name().to_owned())
            .collect()
    }

    #[getter]
    fn arrows(&self) -> Vec<(String, String, String)> {
        self.inner
            .arrows()
            .map(|arrow| {
                (
                    arrow.source().name().to_owned(),
                    arrow.target().name().to_owned(),
                    arrow.name().to_owned(),
                )
            })
            .collect()
    }

    fn vertex_path(&self, name: &str) -> PyResult<PyPath> {
        match self.inner.generator_by_name(name) {
            Some(Generator::Vertex(vertex)) => Ok(PyPath {
                inner: vertex.path(),
            }),
            _ => Err(invalid(format!("unknown vertex {name:?}"))),
        }
    }

    fn arrow_path(&self, name: &str) -> PyResult<PyPath> {
        match self.inner.generator_by_name(name) {
            Some(Generator::Arrow(arrow)) => Ok(PyPath {
                inner: arrow.path(),
            }),
            _ => Err(invalid(format!("unknown arrow {name:?}"))),
        }
    }

    fn path(&self, arrows: Vec<String>) -> PyResult<PyPath> {
        let mut arrows = arrows.into_iter();
        let first = arrows
            .next()
            .ok_or_else(|| invalid("an empty arrow list does not determine a vertex path"))?;
        let mut path = self.arrow_path(&first)?.inner;
        for name in arrows {
            path = path
                .multiply(&self.arrow_path(&name)?.inner)
                .map_err(invalid)?;
        }
        Ok(PyPath { inner: path })
    }

    fn zero_path(&self) -> PyPath {
        PyPath {
            inner: self.inner.zero(),
        }
    }

    fn is_identical(&self, other: &Self) -> bool {
        self.inner.is_identical(&other.inner)
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }
}

/// A path using QPA's left-to-right multiplication convention.
#[pyclass(name = "Path", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyPath {
    inner: Path,
}

fn endpoint_name(endpoint: PathEndpoint) -> Option<String> {
    match endpoint {
        PathEndpoint::Vertex(vertex) => Some(vertex.name().to_owned()),
        PathEndpoint::Zero(_) => None,
    }
}

#[pymethods]
impl PyPath {
    #[getter]
    fn length(&self) -> Option<usize> {
        self.inner.length()
    }

    #[getter]
    fn source(&self) -> Option<String> {
        endpoint_name(self.inner.source())
    }

    #[getter]
    fn target(&self) -> Option<String> {
        endpoint_name(self.inner.target())
    }

    #[getter]
    fn external_representation(&self) -> Vec<usize> {
        self.inner.external_representation()
    }

    fn multiply(&self, right: &Self) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.multiply(&right.inner).map_err(invalid)?,
        })
    }

    fn __mul__(&self, right: &Self) -> PyResult<Self> {
        self.multiply(right)
    }

    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

#[derive(Clone)]
enum AlgebraKind {
    Rational(PathAlgebra<RationalField>),
    Binary(PathAlgebra<Z2>),
    Prime(PathAlgebra<Zp>),
}

/// Runtime-dispatched path algebra backed by the existing generic Rust type.
#[pyclass(name = "PathAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyPathAlgebra {
    kind: AlgebraKind,
}

#[derive(Clone)]
enum AlgebraElementKind {
    Rational(PathAlgebraElement<RationalField>),
    Binary(PathAlgebraElement<Z2>),
    Prime(PathAlgebraElement<Zp>),
}

/// Sparse element of a runtime-dispatched Python path algebra.
#[pyclass(name = "AlgebraElement", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyAlgebraElement {
    kind: AlgebraElementKind,
}

#[derive(Clone)]
enum QuotientKind {
    Rational(QuotientPathAlgebra<RationalField>),
    Binary(QuotientPathAlgebra<Z2>),
    Prime(QuotientPathAlgebra<Zp>),
}

/// Quotient path algebra constructed by QPA's high-level Groebner route.
#[pyclass(name = "QuotientAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyQuotientAlgebra {
    kind: QuotientKind,
}

#[derive(Clone)]
enum QuotientElementKind {
    Rational(QuotientPathAlgebraElement<RationalField>),
    Binary(QuotientPathAlgebraElement<Z2>),
    Prime(QuotientPathAlgebraElement<Zp>),
}

/// Element of a runtime-dispatched quotient path algebra.
#[pyclass(name = "QuotientElement", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyQuotientElement {
    kind: QuotientElementKind,
}

#[derive(Clone)]
enum MatrixModuleKind {
    Rational(PathAlgebraModule<RationalField>),
    Binary(PathAlgebraModule<Z2>),
    Prime(PathAlgebraModule<Zp>),
}

/// Checked QPA right matrix module over a path algebra or quotient.
#[pyclass(name = "MatrixModule", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyMatrixModule {
    kind: MatrixModuleKind,
}

fn clone_relations(
    py: Python<'_>,
    relations: Vec<Py<PyAlgebraElement>>,
) -> Vec<AlgebraElementKind> {
    relations
        .into_iter()
        .map(|relation| relation.borrow(py).kind.clone())
        .collect()
}

#[pymethods]
impl PyPathAlgebra {
    #[new]
    fn new(field: PyRef<'_, PyField>, quiver: PyRef<'_, PyQuiver>) -> PyResult<Self> {
        let kind = match field.kind {
            FieldKind::Rational => {
                AlgebraKind::Rational(PathAlgebra::new(Q, quiver.inner.clone()).map_err(invalid)?)
            }
            FieldKind::Binary => {
                AlgebraKind::Binary(PathAlgebra::new(Z2, quiver.inner.clone()).map_err(invalid)?)
            }
            FieldKind::Prime(characteristic) => AlgebraKind::Prime(
                PathAlgebra::new(Zp::new(characteristic), quiver.inner.clone()).map_err(invalid)?,
            ),
        };
        Ok(Self { kind })
    }

    #[getter]
    fn field(&self) -> PyField {
        let kind = match &self.kind {
            AlgebraKind::Rational(_) => FieldKind::Rational,
            AlgebraKind::Binary(_) => FieldKind::Binary,
            AlgebraKind::Prime(algebra) => FieldKind::Prime(algebra.field().get_prime()),
        };
        PyField { kind }
    }

    #[getter]
    fn quiver(&self) -> PyQuiver {
        let inner = match &self.kind {
            AlgebraKind::Rational(algebra) => algebra.quiver(),
            AlgebraKind::Binary(algebra) => algebra.quiver(),
            AlgebraKind::Prime(algebra) => algebra.quiver(),
        }
        .clone();
        PyQuiver { inner }
    }

    fn generator(&self, name: &str) -> PyResult<PyAlgebraElement> {
        let path = match &self.kind {
            AlgebraKind::Rational(algebra) => algebra.quiver(),
            AlgebraKind::Binary(algebra) => algebra.quiver(),
            AlgebraKind::Prime(algebra) => algebra.quiver(),
        }
        .generator_by_name(name)
        .ok_or_else(|| invalid(format!("unknown quiver generator {name:?}")))?
        .path();
        let kind = match &self.kind {
            AlgebraKind::Rational(algebra) => {
                AlgebraElementKind::Rational(algebra.element(path).map_err(invalid)?)
            }
            AlgebraKind::Binary(algebra) => {
                AlgebraElementKind::Binary(algebra.element(path).map_err(invalid)?)
            }
            AlgebraKind::Prime(algebra) => {
                AlgebraElementKind::Prime(algebra.element(path).map_err(invalid)?)
            }
        };
        Ok(PyAlgebraElement { kind })
    }

    fn term(&self, path: &PyPath, coefficient: &PyScalar) -> PyResult<PyAlgebraElement> {
        let kind = match &self.kind {
            AlgebraKind::Rational(algebra) => AlgebraElementKind::Rational(
                algebra
                    .from_terms([(path.inner.clone(), coefficient.rational_value()?)])
                    .map_err(invalid)?,
            ),
            AlgebraKind::Binary(algebra) => AlgebraElementKind::Binary(
                algebra
                    .from_terms([(path.inner.clone(), coefficient.binary_value()?)])
                    .map_err(invalid)?,
            ),
            AlgebraKind::Prime(algebra) => AlgebraElementKind::Prime(
                algebra
                    .from_terms([(
                        path.inner.clone(),
                        coefficient.prime_value(algebra.field())?,
                    )])
                    .map_err(invalid)?,
            ),
        };
        Ok(PyAlgebraElement { kind })
    }

    fn zero(&self) -> PyAlgebraElement {
        let kind = match &self.kind {
            AlgebraKind::Rational(algebra) => AlgebraElementKind::Rational(algebra.zero()),
            AlgebraKind::Binary(algebra) => AlgebraElementKind::Binary(algebra.zero()),
            AlgebraKind::Prime(algebra) => AlgebraElementKind::Prime(algebra.zero()),
        };
        PyAlgebraElement { kind }
    }

    fn one(&self) -> PyAlgebraElement {
        let kind = match &self.kind {
            AlgebraKind::Rational(algebra) => AlgebraElementKind::Rational(algebra.one()),
            AlgebraKind::Binary(algebra) => AlgebraElementKind::Binary(algebra.one()),
            AlgebraKind::Prime(algebra) => AlgebraElementKind::Prime(algebra.one()),
        };
        PyAlgebraElement { kind }
    }

    fn quotient(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyQuotientAlgebra> {
        let relations = clone_relations(py, relations);
        let kind = match &self.kind {
            AlgebraKind::Rational(algebra) => {
                let relations = relations
                    .into_iter()
                    .map(|relation| match relation {
                        AlgebraElementKind::Rational(relation) => Ok(relation),
                        _ => Err(invalid("relation belongs to another coefficient field")),
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                QuotientKind::Rational(
                    QuotientPathAlgebra::by_relations_high_level(algebra.clone(), &relations)
                        .map_err(invalid)?,
                )
            }
            AlgebraKind::Binary(algebra) => {
                let relations = relations
                    .into_iter()
                    .map(|relation| match relation {
                        AlgebraElementKind::Binary(relation) => Ok(relation),
                        _ => Err(invalid("relation belongs to another coefficient field")),
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                QuotientKind::Binary(
                    QuotientPathAlgebra::by_relations_high_level(algebra.clone(), &relations)
                        .map_err(invalid)?,
                )
            }
            AlgebraKind::Prime(algebra) => {
                let relations = relations
                    .into_iter()
                    .map(|relation| match relation {
                        AlgebraElementKind::Prime(relation) => Ok(relation),
                        _ => Err(invalid("relation belongs to another coefficient field")),
                    })
                    .collect::<PyResult<Vec<_>>>()?;
                QuotientKind::Prime(
                    QuotientPathAlgebra::by_relations_high_level(algebra.clone(), &relations)
                        .map_err(invalid)?,
                )
            }
        };
        Ok(PyQuotientAlgebra { kind })
    }

    fn module_from_integer_matrices(
        &self,
        dimension_vector: Vec<usize>,
        actions: Vec<(String, Vec<Vec<i64>>)>,
    ) -> PyResult<PyMatrixModule> {
        let kind = match &self.kind {
            AlgebraKind::Rational(algebra) => MatrixModuleKind::Rational(build_path_module(
                algebra,
                dimension_vector,
                actions,
                |value| Q.nth(value.clone()),
            )?),
            AlgebraKind::Binary(algebra) => MatrixModuleKind::Binary(build_path_module(
                algebra,
                dimension_vector,
                actions,
                |value| value.to_finite_field(&Z2),
            )?),
            AlgebraKind::Prime(algebra) => {
                let field = algebra.field().clone();
                MatrixModuleKind::Prime(build_path_module(
                    algebra,
                    dimension_vector,
                    actions,
                    |value| value.to_finite_field(&field),
                )?)
            }
        };
        Ok(PyMatrixModule { kind })
    }

    fn is_identical(&self, other: &Self) -> bool {
        match (&self.kind, &other.kind) {
            (AlgebraKind::Rational(left), AlgebraKind::Rational(right)) => left.is_identical(right),
            (AlgebraKind::Binary(left), AlgebraKind::Binary(right)) => left.is_identical(right),
            (AlgebraKind::Prime(left), AlgebraKind::Prime(right)) => left.is_identical(right),
            _ => false,
        }
    }

    fn __repr__(&self) -> String {
        match &self.kind {
            AlgebraKind::Rational(algebra) => algebra.view_obj_with("Rationals"),
            AlgebraKind::Binary(algebra) => algebra.view_obj_with("GF(2)"),
            AlgebraKind::Prime(algebra) => {
                algebra.view_obj_with(&format!("GF({})", algebra.field().get_prime()))
            }
        }
    }
}

fn build_path_module<F: numerica::domains::Field>(
    algebra: &PathAlgebra<F>,
    dimension_vector: Vec<usize>,
    actions: Vec<(String, Vec<Vec<i64>>)>,
    convert: impl Fn(&Integer) -> F::Element,
) -> PyResult<PathAlgebraModule<F>> {
    let actions = actions
        .into_iter()
        .map(|(label, rows)| {
            let rows = rows
                .into_iter()
                .map(|row| {
                    row.into_iter()
                        .map(|value| convert(&Integer::from(value)))
                        .collect()
                })
                .collect();
            let matrix = DenseMatrix::from_rows(rows, algebra.field().clone()).map_err(invalid)?;
            Ok(NamedArrowMatrix::new(label, matrix))
        })
        .collect::<PyResult<Vec<_>>>()?;
    right_module_over_path_algebra(
        QuiverAlgebra::Path(algebra.clone()),
        dimension_vector,
        &actions,
    )
    .map_err(invalid)?
    .ok_or_else(|| RelationNotSatisfiedError::new_err("module actions do not satisfy relations"))
}

#[pymethods]
impl PyAlgebraElement {
    #[getter]
    fn is_zero(&self) -> bool {
        match &self.kind {
            AlgebraElementKind::Rational(element) => element.is_zero(),
            AlgebraElementKind::Binary(element) => element.is_zero(),
            AlgebraElementKind::Prime(element) => element.is_zero(),
        }
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = match (&self.kind, &right.kind) {
            (AlgebraElementKind::Rational(left), AlgebraElementKind::Rational(right)) => {
                AlgebraElementKind::Rational(left.add(right).map_err(invalid)?)
            }
            (AlgebraElementKind::Binary(left), AlgebraElementKind::Binary(right)) => {
                AlgebraElementKind::Binary(left.add(right).map_err(invalid)?)
            }
            (AlgebraElementKind::Prime(left), AlgebraElementKind::Prime(right)) => {
                AlgebraElementKind::Prime(left.add(right).map_err(invalid)?)
            }
            _ => return Err(invalid("elements belong to different coefficient fields")),
        };
        Ok(Self { kind })
    }

    fn multiply(&self, right: &Self) -> PyResult<Self> {
        let kind = match (&self.kind, &right.kind) {
            (AlgebraElementKind::Rational(left), AlgebraElementKind::Rational(right)) => {
                AlgebraElementKind::Rational(left.multiply(right).map_err(invalid)?)
            }
            (AlgebraElementKind::Binary(left), AlgebraElementKind::Binary(right)) => {
                AlgebraElementKind::Binary(left.multiply(right).map_err(invalid)?)
            }
            (AlgebraElementKind::Prime(left), AlgebraElementKind::Prime(right)) => {
                AlgebraElementKind::Prime(left.multiply(right).map_err(invalid)?)
            }
            _ => return Err(invalid("elements belong to different coefficient fields")),
        };
        Ok(Self { kind })
    }

    fn pow(&self, exponent: u64) -> PyResult<Self> {
        let kind = match &self.kind {
            AlgebraElementKind::Rational(element) => {
                AlgebraElementKind::Rational(element.pow(exponent).map_err(invalid)?)
            }
            AlgebraElementKind::Binary(element) => {
                AlgebraElementKind::Binary(element.pow(exponent).map_err(invalid)?)
            }
            AlgebraElementKind::Prime(element) => {
                AlgebraElementKind::Prime(element.pow(exponent).map_err(invalid)?)
            }
        };
        Ok(Self { kind })
    }

    fn __add__(&self, right: &Self) -> PyResult<Self> {
        self.add(right)
    }

    fn __mul__(&self, right: &Self) -> PyResult<Self> {
        self.multiply(right)
    }

    fn __repr__(&self) -> String {
        match &self.kind {
            AlgebraElementKind::Rational(element) => element.print_obj_with(ToString::to_string),
            AlgebraElementKind::Binary(element) => element.print_obj_with(|value| {
                u32::from(!element.algebra().field().is_zero(value)).to_string()
            }),
            AlgebraElementKind::Prime(element) => element
                .print_obj_with(|value| element.algebra().field().from_element(value).to_string()),
        }
    }
}

#[pymethods]
impl PyQuotientAlgebra {
    #[getter]
    fn field(&self) -> PyField {
        let kind = match &self.kind {
            QuotientKind::Rational(_) => FieldKind::Rational,
            QuotientKind::Binary(_) => FieldKind::Binary,
            QuotientKind::Prime(algebra) => {
                FieldKind::Prime(algebra.original_path_algebra().field().get_prime())
            }
        };
        PyField { kind }
    }

    fn generator(&self, name: &str) -> PyResult<PyQuotientElement> {
        let path = match &self.kind {
            QuotientKind::Rational(algebra) => algebra.quiver(),
            QuotientKind::Binary(algebra) => algebra.quiver(),
            QuotientKind::Prime(algebra) => algebra.quiver(),
        }
        .generator_by_name(name)
        .ok_or_else(|| invalid(format!("unknown quiver generator {name:?}")))?
        .path();
        let kind = match &self.kind {
            QuotientKind::Rational(algebra) => {
                QuotientElementKind::Rational(algebra.path(path).map_err(invalid)?)
            }
            QuotientKind::Binary(algebra) => {
                QuotientElementKind::Binary(algebra.path(path).map_err(invalid)?)
            }
            QuotientKind::Prime(algebra) => {
                QuotientElementKind::Prime(algebra.path(path).map_err(invalid)?)
            }
        };
        Ok(PyQuotientElement { kind })
    }

    fn zero(&self) -> PyQuotientElement {
        let kind = match &self.kind {
            QuotientKind::Rational(algebra) => QuotientElementKind::Rational(algebra.zero()),
            QuotientKind::Binary(algebra) => QuotientElementKind::Binary(algebra.zero()),
            QuotientKind::Prime(algebra) => QuotientElementKind::Prime(algebra.zero()),
        };
        PyQuotientElement { kind }
    }

    fn one(&self) -> PyResult<PyQuotientElement> {
        let kind = match &self.kind {
            QuotientKind::Rational(algebra) => {
                QuotientElementKind::Rational(algebra.one().map_err(invalid)?)
            }
            QuotientKind::Binary(algebra) => {
                QuotientElementKind::Binary(algebra.one().map_err(invalid)?)
            }
            QuotientKind::Prime(algebra) => {
                QuotientElementKind::Prime(algebra.one().map_err(invalid)?)
            }
        };
        Ok(PyQuotientElement { kind })
    }

    fn dimension(&self) -> PyResult<Option<String>> {
        let dimension = match &self.kind {
            QuotientKind::Rational(algebra) => algebra.dimension().map_err(invalid)?,
            QuotientKind::Binary(algebra) => algebra.dimension().map_err(invalid)?,
            QuotientKind::Prime(algebra) => algebra.dimension().map_err(invalid)?,
        };
        Ok(match dimension {
            AlgebraDimension::Finite(value) => Some(value.to_string()),
            AlgebraDimension::Infinite => None,
        })
    }

    fn is_identical(&self, other: &Self) -> bool {
        match (&self.kind, &other.kind) {
            (QuotientKind::Rational(left), QuotientKind::Rational(right)) => {
                left.is_identical(right)
            }
            (QuotientKind::Binary(left), QuotientKind::Binary(right)) => left.is_identical(right),
            (QuotientKind::Prime(left), QuotientKind::Prime(right)) => left.is_identical(right),
            _ => false,
        }
    }

    fn __repr__(&self) -> String {
        match &self.kind {
            QuotientKind::Rational(algebra) => algebra.view_obj_with("Rationals"),
            QuotientKind::Binary(algebra) => algebra.view_obj_with("GF(2)"),
            QuotientKind::Prime(algebra) => algebra.view_obj_with(&format!(
                "GF({})",
                algebra.original_path_algebra().field().get_prime()
            )),
        }
    }
}

#[pymethods]
impl PyQuotientElement {
    #[getter]
    fn is_zero(&self) -> bool {
        match &self.kind {
            QuotientElementKind::Rational(element) => element.is_zero(),
            QuotientElementKind::Binary(element) => element.is_zero(),
            QuotientElementKind::Prime(element) => element.is_zero(),
        }
    }

    #[getter]
    fn is_normal_form(&self) -> bool {
        match &self.kind {
            QuotientElementKind::Rational(element) => element.is_normal_form(),
            QuotientElementKind::Binary(element) => element.is_normal_form(),
            QuotientElementKind::Prime(element) => element.is_normal_form(),
        }
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = match (&self.kind, &right.kind) {
            (QuotientElementKind::Rational(left), QuotientElementKind::Rational(right)) => {
                QuotientElementKind::Rational(left.add(right).map_err(invalid)?)
            }
            (QuotientElementKind::Binary(left), QuotientElementKind::Binary(right)) => {
                QuotientElementKind::Binary(left.add(right).map_err(invalid)?)
            }
            (QuotientElementKind::Prime(left), QuotientElementKind::Prime(right)) => {
                QuotientElementKind::Prime(left.add(right).map_err(invalid)?)
            }
            _ => return Err(invalid("quotient elements belong to different fields")),
        };
        Ok(Self { kind })
    }

    fn multiply(&self, right: &Self) -> PyResult<Self> {
        let kind = match (&self.kind, &right.kind) {
            (QuotientElementKind::Rational(left), QuotientElementKind::Rational(right)) => {
                QuotientElementKind::Rational(left.multiply(right).map_err(invalid)?)
            }
            (QuotientElementKind::Binary(left), QuotientElementKind::Binary(right)) => {
                QuotientElementKind::Binary(left.multiply(right).map_err(invalid)?)
            }
            (QuotientElementKind::Prime(left), QuotientElementKind::Prime(right)) => {
                QuotientElementKind::Prime(left.multiply(right).map_err(invalid)?)
            }
            _ => return Err(invalid("quotient elements belong to different fields")),
        };
        Ok(Self { kind })
    }

    fn __add__(&self, right: &Self) -> PyResult<Self> {
        self.add(right)
    }

    fn __mul__(&self, right: &Self) -> PyResult<Self> {
        self.multiply(right)
    }

    fn __repr__(&self) -> String {
        match &self.kind {
            QuotientElementKind::Rational(element) => element.print_obj_with(ToString::to_string),
            QuotientElementKind::Binary(element) => element.print_obj_with(|value| {
                u32::from(
                    !element
                        .algebra()
                        .original_path_algebra()
                        .field()
                        .is_zero(value),
                )
                .to_string()
            }),
            QuotientElementKind::Prime(element) => element.print_obj_with(|value| {
                element
                    .algebra()
                    .original_path_algebra()
                    .field()
                    .from_element(value)
                    .to_string()
            }),
        }
    }
}

#[pymethods]
impl PyMatrixModule {
    #[getter]
    fn dimension_vector(&self) -> Vec<usize> {
        match &self.kind {
            MatrixModuleKind::Rational(module) => module.dimension_vector(),
            MatrixModuleKind::Binary(module) => module.dimension_vector(),
            MatrixModuleKind::Prime(module) => module.dimension_vector(),
        }
        .to_vec()
    }

    #[getter]
    fn dimension(&self) -> usize {
        match &self.kind {
            MatrixModuleKind::Rational(module) => module.dimension(),
            MatrixModuleKind::Binary(module) => module.dimension(),
            MatrixModuleKind::Prime(module) => module.dimension(),
        }
    }

    #[getter]
    fn arrow_matrix_shapes(&self) -> Vec<(usize, usize)> {
        match &self.kind {
            MatrixModuleKind::Rational(module) => matrix_shapes(module.arrow_matrices()),
            MatrixModuleKind::Binary(module) => matrix_shapes(module.arrow_matrices()),
            MatrixModuleKind::Prime(module) => matrix_shapes(module.arrow_matrices()),
        }
    }

    fn is_identical(&self, other: &Self) -> bool {
        match (&self.kind, &other.kind) {
            (MatrixModuleKind::Rational(left), MatrixModuleKind::Rational(right)) => {
                left.is_identical(right)
            }
            (MatrixModuleKind::Binary(left), MatrixModuleKind::Binary(right)) => {
                left.is_identical(right)
            }
            (MatrixModuleKind::Prime(left), MatrixModuleKind::Prime(right)) => {
                left.is_identical(right)
            }
            _ => false,
        }
    }

    fn __repr__(&self) -> String {
        match &self.kind {
            MatrixModuleKind::Rational(module) => module.qpa_print_string(),
            MatrixModuleKind::Binary(module) => module.qpa_print_string(),
            MatrixModuleKind::Prime(module) => module.qpa_print_string(),
        }
    }
}

fn matrix_shapes<F: numerica::domains::Field>(matrices: &[DenseMatrix<F>]) -> Vec<(usize, usize)> {
    matrices
        .iter()
        .map(|matrix| (matrix.nrows(), matrix.ncols()))
        .collect()
}

#[pymodule]
fn _qpa(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("QpaError", module.py().get_type::<QpaError>())?;
    module.add(
        "InvalidInputError",
        module.py().get_type::<InvalidInputError>(),
    )?;
    module.add(
        "RelationNotSatisfiedError",
        module.py().get_type::<RelationNotSatisfiedError>(),
    )?;
    module.add_class::<PyField>()?;
    module.add_class::<PyScalar>()?;
    module.add_class::<PyQuiver>()?;
    module.add_class::<PyPath>()?;
    module.add_class::<PyPathAlgebra>()?;
    module.add_class::<PyAlgebraElement>()?;
    module.add_class::<PyQuotientAlgebra>()?;
    module.add_class::<PyQuotientElement>()?;
    module.add_class::<PyMatrixModule>()?;
    Ok(())
}
