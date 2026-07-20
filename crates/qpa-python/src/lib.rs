//! `PyO3` bindings over the source-faithful QPA Rust surface.

#![allow(clippy::needless_pass_by_value)]
#![allow(missing_docs)] // PyO3 exception structs cannot receive rustdoc attributes.

use std::{cmp::Ordering, str::FromStr};

mod runtime;

use numerica::domains::{
    finite_field::{PrimeIteratorU64, Z2, Zp},
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
    AlgebraDimension, ArrowSpec, DenseMatrix, DynkinOrientation, Generator, GroebnerBasis,
    GroebnerBasisKind, NamedArrowMatrix, NamedArrowSpec, Path, PathAlgebra, PathAlgebraElement,
    PathAlgebraIdeal, PathAlgebraModule, PathEndpoint, Quiver, QuiverAlgebra, QuiverCardinality,
    QuotientPathAlgebra, QuotientPathAlgebraElement, Vertex, VertexIndex,
    right_module_over_path_algebra,
};

use runtime::{PythonField as _, dispatch_pair_ref, dispatch_ref, map_pair_ref, map_ref};

// Base exception for failures reported by the Rust QPA binding.
create_exception!(_qpa, QpaError, PyException);
// Invalid mathematical input, field ownership, or object-family exception.
create_exception!(_qpa, InvalidInputError, QpaError);
// QPA's explicit `fail` result for actions that do not satisfy relations.
create_exception!(_qpa, RelationNotSatisfiedError, QpaError);

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
            FieldKind::Rational => runtime::scalar_fraction(&Q, &numerator, &denominator),
            FieldKind::Binary => runtime::scalar_fraction(&Z2, &numerator, &denominator),
            FieldKind::Prime(characteristic) => {
                runtime::scalar_fraction(&Zp::new(characteristic), &numerator, &denominator)
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

    fn value_for<F: runtime::PythonField>(&self, field: &F) -> PyResult<F::Element> {
        field.scalar_value(self)
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

    #[staticmethod]
    fn named(
        vertices: Vec<String>,
        arrows: Vec<(String, String, Option<String>)>,
    ) -> PyResult<Self> {
        let arrows = arrows
            .into_iter()
            .map(|(source, target, name)| match name {
                Some(name) => NamedArrowSpec::named(source, target, name),
                None => NamedArrowSpec::unnamed(source, target),
            })
            .collect();
        let inner = Quiver::from_named(
            vertices.into_iter().map(String::into_boxed_str).collect(),
            arrows,
        )
        .map_err(invalid)?;
        Ok(Self { inner })
    }

    #[staticmethod]
    fn indexed(
        number_of_vertices: usize,
        arrows: Vec<(usize, usize, Option<String>)>,
    ) -> PyResult<Self> {
        let arrows = arrows
            .into_iter()
            .map(|(source, target, name)| {
                let source = VertexIndex::new(source)
                    .ok_or_else(|| invalid("QPA vertex indices are one-based"))?;
                let target = VertexIndex::new(target)
                    .ok_or_else(|| invalid("QPA vertex indices are one-based"))?;
                Ok(name.map_or_else(
                    || ArrowSpec::unnamed(source, target),
                    |name| ArrowSpec::named(source, target, name),
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self {
            inner: Quiver::from_indexed(number_of_vertices, arrows).map_err(invalid)?,
        })
    }

    #[staticmethod]
    fn from_adjacency(adjacency: Vec<Vec<usize>>) -> PyResult<Self> {
        Ok(Self {
            inner: Quiver::from_adjacency(adjacency).map_err(invalid)?,
        })
    }

    #[staticmethod]
    fn dynkin(family: &str, rank: usize, orientation: Vec<String>) -> PyResult<Option<Self>> {
        let orientation = orientation
            .iter()
            .map(|symbol| parse_dynkin_orientation(symbol))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Quiver::dynkin_from_name(family, rank, &orientation)
            .map_err(invalid)?
            .map(|inner| Self { inner }))
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

    #[getter]
    fn number_of_vertices(&self) -> usize {
        self.inner.number_of_vertices()
    }

    #[getter]
    fn number_of_arrows(&self) -> usize {
        self.inner.number_of_arrows()
    }

    #[getter]
    fn adjacency_matrix(&self) -> Vec<Vec<usize>> {
        self.inner
            .adjacency_matrix()
            .iter()
            .map(|row| row.to_vec())
            .collect()
    }

    #[getter]
    fn is_acyclic(&self) -> bool {
        self.inner.is_acyclic()
    }

    #[getter]
    fn is_finite(&self) -> bool {
        self.inner.is_finite()
    }

    #[getter]
    fn is_connected(&self) -> bool {
        self.inner.is_connected()
    }

    #[getter]
    fn is_unoriented_acyclic(&self) -> bool {
        self.inner.is_unoriented_acyclic()
    }

    #[getter]
    fn is_tree(&self) -> bool {
        self.inner.is_tree()
    }

    #[getter]
    fn has_special_arrow_names(&self) -> bool {
        self.inner.has_special_single_letter_arrow_names()
    }

    #[getter]
    fn cardinality(&self) -> Option<String> {
        match self.inner.cardinality() {
            QuiverCardinality::Finite(value) => Some(value.to_string()),
            QuiverCardinality::Infinite => None,
        }
    }

    #[getter]
    fn topological_order(&self) -> Option<Vec<String>> {
        self.inner.topological_order().map(|vertices| {
            vertices
                .into_iter()
                .map(|vertex| vertex.name().to_owned())
                .collect()
        })
    }

    #[getter]
    fn dynkin_type(&self) -> Option<String> {
        self.inner.dynkin_type().map(|value| value.to_string())
    }

    fn incoming_arrows(&self, vertex: &str) -> PyResult<Vec<String>> {
        Ok(self
            .inner
            .incoming_arrows(&self.vertex_by_name(vertex)?)
            .map_err(invalid)?
            .into_iter()
            .map(|arrow| arrow.name().to_owned())
            .collect())
    }

    fn outgoing_arrows(&self, vertex: &str) -> PyResult<Vec<String>> {
        Ok(self
            .inner
            .outgoing_arrows(&self.vertex_by_name(vertex)?)
            .map_err(invalid)?
            .into_iter()
            .map(|arrow| arrow.name().to_owned())
            .collect())
    }

    fn neighbors(&self, vertex: &str) -> PyResult<Vec<String>> {
        Ok(self
            .inner
            .neighbors(&self.vertex_by_name(vertex)?)
            .map_err(invalid)?
            .into_iter()
            .map(|vertex| vertex.name().to_owned())
            .collect())
    }

    fn compare_paths(&self, left: &PyPath, right: &PyPath) -> PyResult<i8> {
        Ok(
            match self
                .inner
                .compare_paths(&left.inner, &right.inner)
                .map_err(invalid)?
            {
                Ordering::Less => -1,
                Ordering::Equal => 0,
                Ordering::Greater => 1,
            },
        )
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

    fn path_at(&self, position: usize) -> Option<PyPath> {
        self.inner.path_at(position).map(|inner| PyPath { inner })
    }

    fn path_from_external_representation(
        &self,
        representation: Vec<usize>,
    ) -> PyResult<Option<PyPath>> {
        Ok(
            Path::from_external_representation(&self.inner, &representation)
                .map_err(invalid)?
                .map(|inner| PyPath { inner }),
        )
    }

    fn full_subquiver(&self, vertices: Vec<String>) -> PyResult<Self> {
        let vertices = vertices
            .iter()
            .map(|name| self.vertex_by_name(name))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self {
            inner: self.inner.full_subquiver(&vertices).map_err(invalid)?,
        })
    }

    fn connected_components(&self) -> PyResult<Vec<Self>> {
        Ok(self
            .inner
            .connected_components()
            .map_err(invalid)?
            .into_iter()
            .map(|inner| Self { inner })
            .collect())
    }

    fn opposite(&self) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.opposite().map_err(invalid)?,
        })
    }

    fn separated(&self) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.separated().map_err(invalid)?,
        })
    }

    fn double(&self) -> PyResult<Self> {
        Ok(Self {
            inner: self.inner.double().map_err(invalid)?,
        })
    }

    #[getter]
    fn is_special_biserial(&self) -> bool {
        qpa_gbnp::is_special_biserial_quiver(&self.inner)
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

impl PyQuiver {
    fn vertex_by_name(&self, name: &str) -> PyResult<Vertex> {
        match self.inner.generator_by_name(name) {
            Some(Generator::Vertex(vertex)) => Ok(vertex),
            _ => Err(invalid(format!("unknown vertex {name:?}"))),
        }
    }
}

fn parse_dynkin_orientation(symbol: &str) -> PyResult<DynkinOrientation> {
    match symbol {
        "l" => Ok(DynkinOrientation::Left),
        "r" => Ok(DynkinOrientation::Right),
        "u" => Ok(DynkinOrientation::Up),
        "d" => Ok(DynkinOrientation::Down),
        _ => Err(invalid(format!(
            "invalid Dynkin orientation {symbol:?}; expected l, r, u, or d"
        ))),
    }
}

const fn ordering_value(ordering: Ordering) -> i8 {
    match ordering {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
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
    fn quiver(&self) -> PyQuiver {
        PyQuiver {
            inner: self.inner.quiver().clone(),
        }
    }

    #[getter]
    fn is_zero(&self) -> bool {
        self.inner.is_zero()
    }

    #[getter]
    fn is_vertex(&self) -> bool {
        self.inner.is_vertex()
    }

    #[getter]
    fn is_arrow(&self) -> bool {
        self.inner.is_arrow()
    }

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

    fn opposite(&self) -> PyResult<Self> {
        Ok(Self {
            inner: self
                .inner
                .quiver()
                .opposite_path(&self.inner)
                .map_err(invalid)?,
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

macro_rules! impl_algebra_element_extract {
    ($field:ty, $variant:ident) => {
        impl runtime::Extract<PathAlgebraElement<$field>> for AlgebraElementKind {
            fn extract(self) -> PyResult<PathAlgebraElement<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid("relation belongs to another coefficient field")),
                }
            }
        }

        impl runtime::Wrap<PathAlgebraElement<$field>> for AlgebraElementKind {
            fn wrap(value: PathAlgebraElement<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_algebra_element_extract!(RationalField, Rational);
impl_algebra_element_extract!(Z2, Binary);
impl_algebra_element_extract!(Zp, Prime);

/// Sparse element of a runtime-dispatched Python path algebra.
#[pyclass(name = "AlgebraElement", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyAlgebraElement {
    kind: AlgebraElementKind,
}

#[derive(Clone)]
enum IdealKind {
    Rational(PathAlgebraIdeal<RationalField>),
    Binary(PathAlgebraIdeal<Z2>),
    Prime(PathAlgebraIdeal<Zp>),
}

/// Two-sided path-algebra ideal using QPA's lazy Gröbner-basis behavior.
#[pyclass(name = "PathAlgebraIdeal", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyPathAlgebraIdeal {
    kind: IdealKind,
}

#[derive(Clone)]
enum GroebnerKind {
    Rational(GroebnerBasis<RationalField>),
    Binary(GroebnerBasis<Z2>),
    Prime(GroebnerBasis<Zp>),
}

/// QPA's two-sided or right Gröbner-basis object.
#[pyclass(name = "GroebnerBasis", frozen, skip_from_py_object)]
#[derive(Clone)]
struct PyGroebnerBasis {
    kind: GroebnerKind,
}

macro_rules! impl_runtime_wrap {
    ($wrapper:ty, $value:ty, $variant:ident) => {
        impl runtime::Wrap<$value> for $wrapper {
            fn wrap(value: $value) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_runtime_wrap!(IdealKind, PathAlgebraIdeal<RationalField>, Rational);
impl_runtime_wrap!(IdealKind, PathAlgebraIdeal<Z2>, Binary);
impl_runtime_wrap!(IdealKind, PathAlgebraIdeal<Zp>, Prime);
impl_runtime_wrap!(GroebnerKind, GroebnerBasis<RationalField>, Rational);
impl_runtime_wrap!(GroebnerKind, GroebnerBasis<Z2>, Binary);
impl_runtime_wrap!(GroebnerKind, GroebnerBasis<Zp>, Prime);

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

macro_rules! impl_quotient_element_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Extract<QuotientPathAlgebraElement<$field>> for QuotientElementKind {
            fn extract(self) -> PyResult<QuotientPathAlgebraElement<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "quotient element belongs to another coefficient field",
                    )),
                }
            }
        }

        impl runtime::Wrap<QuotientPathAlgebraElement<$field>> for QuotientElementKind {
            fn wrap(value: QuotientPathAlgebraElement<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_quotient_element_runtime!(RationalField, Rational);
impl_quotient_element_runtime!(Z2, Binary);
impl_quotient_element_runtime!(Zp, Prime);

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

fn wrap_algebra_elements<F: runtime::PythonField>(
    elements: Vec<PathAlgebraElement<F>>,
) -> Vec<PyAlgebraElement>
where
    AlgebraElementKind: runtime::Wrap<PathAlgebraElement<F>>,
{
    runtime::wrap_all::<AlgebraElementKind, _>(elements)
        .into_iter()
        .map(|kind| PyAlgebraElement { kind })
        .collect()
}

fn wrap_groebner_basis<F: runtime::PythonField>(basis: GroebnerBasis<F>) -> PyGroebnerBasis
where
    GroebnerKind: runtime::Wrap<GroebnerBasis<F>>,
{
    PyGroebnerBasis {
        kind: runtime::Wrap::wrap(basis),
    }
}

fn wrap_quotient_elements<F: runtime::PythonField>(
    elements: Vec<QuotientPathAlgebraElement<F>>,
) -> Vec<PyQuotientElement>
where
    QuotientElementKind: runtime::Wrap<QuotientPathAlgebraElement<F>>,
{
    runtime::wrap_all::<QuotientElementKind, _>(elements)
        .into_iter()
        .map(|kind| PyQuotientElement { kind })
        .collect()
}

fn build_ideal<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    generators: Vec<AlgebraElementKind>,
) -> PyResult<PathAlgebraIdeal<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    PathAlgebraIdeal::new(algebra.clone(), runtime::extract_all(generators)?).map_err(invalid)
}

fn build_groebner_basis<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<GroebnerBasis<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    GroebnerBasis::high_level(algebra, &runtime::extract_all(relations)?).map_err(invalid)
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
        let kind = dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .field()
            .python_kind());
        PyField { kind }
    }

    #[getter]
    fn quiver(&self) -> PyQuiver {
        let inner = dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra.quiver().clone());
        PyQuiver { inner }
    }

    #[getter]
    fn global_dimension(&self) -> usize {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .qpa_global_dimension())
    }

    #[getter]
    fn has_finite_canonical_basis(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .has_finite_canonical_basis())
    }

    #[getter]
    fn is_monomial(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_monomial(algebra)
        })
    }

    #[getter]
    fn is_special_biserial(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_special_biserial(algebra)
        })
    }

    #[getter]
    fn is_string(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_string(algebra)
        })
    }

    #[getter]
    fn is_gentle(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_gentle(algebra)
        })
    }

    #[getter]
    fn generators(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| wrap_algebra_elements(
            algebra.generators()
        ))
    }

    #[getter]
    fn relations(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| wrap_algebra_elements(
            algebra.relations()
        ))
    }

    fn canonical_basis(&self) -> Option<Vec<PyAlgebraElement>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .canonical_basis()
            .map(wrap_algebra_elements))
    }

    fn opposite(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraKind,
            |algebra| algebra.opposite().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn generator(&self, name: &str) -> PyResult<PyAlgebraElement> {
        let path = dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .quiver()
            .generator_by_name(name)
            .ok_or_else(|| invalid(format!("unknown quiver generator {name:?}")))?
            .path());
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra.element(path).map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn term(&self, path: &PyPath, coefficient: &PyScalar) -> PyResult<PyAlgebraElement> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra
                .from_terms([(
                    path.inner.clone(),
                    coefficient.value_for(algebra.field())?,
                )])
                .map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn zero(&self) -> PyAlgebraElement {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra.zero()
        );
        PyAlgebraElement { kind }
    }

    fn one(&self) -> PyAlgebraElement {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra.one()
        );
        PyAlgebraElement { kind }
    }

    fn quotient(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyQuotientAlgebra> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => QuotientKind,
            |algebra| build_quotient(algebra, relations)?
        );
        Ok(PyQuotientAlgebra { kind })
    }

    fn quotient_high_level(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyQuotientAlgebra> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => QuotientKind,
            |algebra| build_high_level_quotient(algebra, relations)?
        );
        Ok(PyQuotientAlgebra { kind })
    }

    fn ideal(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyPathAlgebraIdeal> {
        let generators = clone_relations(py, generators);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => IdealKind,
            |algebra| build_ideal(algebra, generators)?
        );
        Ok(PyPathAlgebraIdeal { kind })
    }

    fn groebner_basis(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyGroebnerBasis> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => GroebnerKind,
            |algebra| build_groebner_basis(algebra, relations)?
        );
        Ok(PyGroebnerBasis { kind })
    }

    fn gbnp_groebner_basis(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyGroebnerBasis> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => GroebnerKind,
            |algebra| build_gbnp_groebner_basis(algebra, relations)?
        );
        Ok(PyGroebnerBasis { kind })
    }

    fn make_uniform_on_right(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<Vec<PyAlgebraElement>> {
        let generators = clone_relations(py, generators);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            let generators = runtime::extract_all(generators)?;
            Ok(wrap_algebra_elements(
                algebra
                    .make_uniform_on_right(&generators)
                    .map_err(invalid)?,
            ))
        })
    }

    fn generators_times_arrows_on_right(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<Vec<PyAlgebraElement>> {
        let generators = clone_relations(py, generators);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            let generators = runtime::extract_all(generators)?;
            Ok(wrap_algebra_elements(
                algebra
                    .generators_times_arrows_on_right(&generators)
                    .map_err(invalid)?,
            ))
        })
    }

    fn nth_power_of_arrow_ideal(&self, exponent: i64) -> PyResult<Vec<PyAlgebraElement>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| Ok(
            wrap_algebra_elements(
                algebra
                    .nth_power_of_arrow_ideal(exponent)
                    .map_err(invalid)?,
            )
        ))
    }

    fn module_from_integer_matrices(
        &self,
        dimension_vector: Vec<usize>,
        actions: Vec<(String, Vec<Vec<i64>>)>,
    ) -> PyResult<PyMatrixModule> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => MatrixModuleKind,
            |algebra| build_path_module(algebra, dimension_vector, actions)?
        );
        Ok(PyMatrixModule { kind })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            AlgebraKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .view_obj_with(&algebra.field().python_kind().display()))
    }
}

fn build_quotient<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<QuotientPathAlgebra<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let relations = runtime::extract_all(relations)?;
    qpa_gbnp::quotient_by_relations(algebra.clone(), &relations).map_err(invalid)
}

fn build_high_level_quotient<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<QuotientPathAlgebra<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let relations = runtime::extract_all(relations)?;
    QuotientPathAlgebra::by_relations_high_level(algebra.clone(), &relations).map_err(invalid)
}

fn build_gbnp_groebner_basis<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<GroebnerBasis<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let relations = runtime::extract_all(relations)?;
    let complete = qpa_gbnp::gbnp_groebner_basis(&relations, algebra).map_err(invalid)?;
    GroebnerBasis::from_complete_relations(algebra, &complete).map_err(invalid)
}

fn build_ideal_quotient<F: runtime::PythonField>(
    ideal: &PathAlgebraIdeal<F>,
) -> PyResult<QuotientPathAlgebra<F>> {
    if ideal.groebner_basis().is_none() {
        let complete = qpa_gbnp::gbnp_groebner_basis_of_ideal(ideal).map_err(invalid)?;
        ideal
            .with_complete_groebner_basis(&complete)
            .map_err(invalid)?;
    }
    QuotientPathAlgebra::by_ideal(ideal).map_err(invalid)
}

fn build_path_module<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    dimension_vector: Vec<usize>,
    actions: Vec<(String, Vec<Vec<i64>>)>,
) -> PyResult<PathAlgebraModule<F>> {
    let actions = actions
        .into_iter()
        .map(|(label, rows)| {
            let rows = rows
                .into_iter()
                .map(|row| {
                    row.into_iter()
                        .map(|value| algebra.field().integer_value(&Integer::from(value)))
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
    fn algebra(&self) -> PyPathAlgebra {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraKind,
            |element| element.algebra().clone()
        );
        PyPathAlgebra { kind }
    }

    #[getter]
    fn terms(&self) -> Vec<(PyPath, PyScalar)> {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .terms()
            .iter()
            .map(|term| (
                PyPath {
                    inner: term.path().clone(),
                },
                element.algebra().field().wrap_scalar(term.coefficient()),
            ))
            .collect())
    }

    #[getter]
    fn external_representation(&self) -> Vec<(Vec<usize>, PyScalar)> {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .external_representation()
            .terms()
            .iter()
            .map(|(path, coefficient)| (
                path.clone(),
                element.algebra().field().wrap_scalar(coefficient),
            ))
            .collect())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element.is_zero())
    }

    #[getter]
    fn is_left_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .is_left_uniform())
    }

    #[getter]
    fn is_right_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .is_right_uniform())
    }

    #[getter]
    fn is_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .is_uniform())
    }

    #[getter]
    fn leading_term(&self) -> Self {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.leading_term()
        );
        Self { kind }
    }

    #[getter]
    fn leading_coefficient(&self) -> PyScalar {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .algebra()
            .field()
            .wrap_scalar(&element.leading_coefficient()))
    }

    #[getter]
    fn leading_monomial(&self) -> PyPath {
        let inner = dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .leading_monomial());
        PyPath { inner }
    }

    #[getter]
    fn vertex_position(&self) -> PyResult<usize> {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .vertex_position()
            .map_err(invalid))
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind => AlgebraElementKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn negate(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.neg().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn subtract(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind => AlgebraElementKind,
            |left, right| left.sub(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn multiply(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind => AlgebraElementKind,
            |left, right| left.multiply(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn pow(&self, exponent: u64) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.pow(exponent).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_left(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element
                .scale_left(&scalar.value_for(element.algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_right(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element
                .scale_right(&scalar.value_for(element.algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn divide_scalar(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element
                .divide_scalar(&scalar.value_for(element.algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn compare(&self, right: &Self) -> PyResult<i8> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind,
            |left, right| Ok(ordering_value(left.compare(right).map_err(invalid)?)),
            Err(invalid("elements belong to different coefficient fields"))
        )
    }

    fn opposite(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.opposite().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn __add__(&self, right: &Self) -> PyResult<Self> {
        self.add(right)
    }

    fn __sub__(&self, right: &Self) -> PyResult<Self> {
        self.subtract(right)
    }

    fn __neg__(&self) -> PyResult<Self> {
        self.negate()
    }

    fn __mul__(&self, right: &Self) -> PyResult<Self> {
        self.multiply(right)
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .print_obj_with(|value| element.algebra().field().scalar_text(value)))
    }

    fn __eq__(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            AlgebraElementKind,
            |left, right| left == right,
            false
        )
    }
}

#[pymethods]
impl PyPathAlgebraIdeal {
    #[getter]
    fn algebra(&self) -> PyPathAlgebra {
        let kind = map_ref!(
            &self.kind,
            IdealKind => AlgebraKind,
            |ideal| ideal.algebra().clone()
        );
        PyPathAlgebra { kind }
    }

    #[getter]
    fn generators(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| wrap_algebra_elements(
            ideal.generators().to_vec()
        ))
    }

    #[getter]
    fn is_trivial(&self) -> bool {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal.is_trivial())
    }

    #[getter]
    fn generators_are_monomial(&self) -> bool {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .generators_are_monomial())
    }

    #[getter]
    fn attached_groebner_basis(&self) -> Option<PyGroebnerBasis> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .groebner_basis()
            .map(wrap_groebner_basis))
    }

    fn with_high_level_groebner_basis(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            IdealKind => IdealKind,
            |ideal| ideal.with_high_level_groebner_basis().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn with_complete_groebner_basis(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<Self> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(&self.kind, IdealKind => IdealKind, |ideal| {
            let relations = runtime::extract_all(relations)?;
            ideal
                .with_complete_groebner_basis(&relations)
                .map_err(invalid)?
        });
        Ok(Self { kind })
    }

    fn right_groebner_basis(&self) -> PyResult<PyGroebnerBasis> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| Ok(wrap_groebner_basis(
            ideal.right_groebner_basis().map_err(invalid)?,
        )))
    }

    fn contains(&self, element: &PyAlgebraElement) -> PyResult<bool> {
        let element = element.kind.clone();
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .contains(runtime::Extract::extract(element)?)
            .map_err(invalid))
    }

    fn is_monomial(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .is_monomial()
            .map_err(invalid))
    }

    fn is_admissible(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .is_admissible()
            .map_err(invalid))
    }

    fn product(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            IdealKind => IdealKind,
            |left, right| left.product(right).map_err(invalid)?,
            invalid("ideals belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn minimal_generating_set(&self) -> PyResult<Vec<PyAlgebraElement>> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| Ok(wrap_algebra_elements(
            qpa_gbnp::ideal_minimal_generating_set(ideal).map_err(invalid)?,
        )))
    }

    fn quotient(&self) -> PyResult<PyQuotientAlgebra> {
        let kind = map_ref!(
            &self.kind,
            IdealKind => QuotientKind,
            |ideal| build_ideal_quotient(ideal)?
        );
        Ok(PyQuotientAlgebra { kind })
    }

    fn quotient_without_groebner_basis(&self) -> PyResult<PyQuotientAlgebra> {
        let kind = map_ref!(
            &self.kind,
            IdealKind => QuotientKind,
            |ideal| QuotientPathAlgebra::by_ideal(ideal).map_err(invalid)?
        );
        Ok(PyQuotientAlgebra { kind })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            IdealKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __contains__(&self, element: &PyAlgebraElement) -> PyResult<bool> {
        self.contains(element)
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .view_obj_with(&ideal.algebra().field().python_kind().display()))
    }
}

#[pymethods]
impl PyGroebnerBasis {
    #[getter]
    fn algebra(&self) -> PyPathAlgebra {
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => AlgebraKind,
            |basis| basis.algebra().clone()
        );
        PyPathAlgebra { kind }
    }

    #[getter]
    fn relations(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| wrap_algebra_elements(
            basis.relations().to_vec()
        ))
    }

    #[getter]
    fn kind(&self) -> &'static str {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| match basis.kind() {
            GroebnerBasisKind::TwoSided => "two-sided",
            GroebnerBasisKind::Right => "right",
        })
    }

    #[getter]
    fn is_complete(&self) -> bool {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis.is_complete())
    }

    #[getter]
    fn tip_reduced_property(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis
            .tip_reduced_property())
    }

    #[getter]
    fn completely_reduced_property(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis
            .completely_reduced_property())
    }

    #[getter]
    fn homogeneous_property(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis
            .homogeneous_property())
    }

    fn tip_reduced_basis(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => GroebnerKind,
            |basis| basis.tip_reduced_basis().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn completely_reduced_basis(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => GroebnerKind,
            |basis| basis.completely_reduced_basis().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn right_basis(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => GroebnerKind,
            |basis| basis.right_basis().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn tip_reduce(&self, element: &PyAlgebraElement) -> PyResult<PyAlgebraElement> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => AlgebraElementKind,
            |basis| basis
                .tip_reduce(runtime::Extract::extract(element)?)
                .map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn completely_reduce(&self, element: &PyAlgebraElement) -> PyResult<PyAlgebraElement> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => AlgebraElementKind,
            |basis| basis
                .completely_reduce(runtime::Extract::extract(element)?)
                .map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn remainder(&self, element: &PyAlgebraElement) -> PyResult<PyAlgebraElement> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => AlgebraElementKind,
            |basis| basis
                .remainder(runtime::Extract::extract(element)?)
                .map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn is_prefix_of_tip_in_tip_ideal(&self, element: &PyAlgebraElement) -> PyResult<bool> {
        let element = element.kind.clone();
        dispatch_ref!(&self.kind, GroebnerKind, |basis| {
            let element = runtime::Extract::extract(element)?;
            basis
                .is_prefix_of_tip_in_tip_ideal(&element)
                .map_err(invalid)
        })
    }

    #[getter]
    fn admits_finitely_many_nontips(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis
            .admits_finitely_many_nontips()
            .map_err(invalid))
    }

    #[getter]
    fn nontip_size(&self) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| Ok(basis
            .nontip_size()
            .map_err(invalid)?
            .map(|value| value.to_string())))
    }

    fn nontips(&self) -> PyResult<Option<Vec<PyPath>>> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| Ok(basis
            .nontips()
            .map_err(invalid)?
            .map(|paths| paths
                .into_iter()
                .map(|inner| PyPath { inner })
                .collect())))
    }

    fn __len__(&self) -> usize {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis.relations().len())
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis.view_string())
    }
}

#[pymethods]
impl PyQuotientAlgebra {
    #[staticmethod]
    fn truncated(
        field: PyRef<'_, PyField>,
        quiver: PyRef<'_, PyQuiver>,
        exponent: i64,
    ) -> PyResult<Self> {
        let kind = match field.kind {
            FieldKind::Rational => QuotientKind::Rational(
                qpa_gbnp::truncated_path_algebra(Q, quiver.inner.clone(), exponent)
                    .map_err(invalid)?,
            ),
            FieldKind::Binary => QuotientKind::Binary(
                qpa_gbnp::truncated_path_algebra(Z2, quiver.inner.clone(), exponent)
                    .map_err(invalid)?,
            ),
            FieldKind::Prime(characteristic) => QuotientKind::Prime(
                qpa_gbnp::truncated_path_algebra(
                    Zp::new(characteristic),
                    quiver.inner.clone(),
                    exponent,
                )
                .map_err(invalid)?,
            ),
        };
        Ok(Self { kind })
    }

    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra
            .original_path_algebra()
            .field()
            .python_kind());
        PyField { kind }
    }

    #[getter]
    fn original_path_algebra(&self) -> PyPathAlgebra {
        let kind = map_ref!(
            &self.kind,
            QuotientKind => AlgebraKind,
            |algebra| algebra.original_path_algebra().clone()
        );
        PyPathAlgebra { kind }
    }

    #[getter]
    fn quiver(&self) -> PyQuiver {
        let inner = dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra.quiver().clone());
        PyQuiver { inner }
    }

    #[getter]
    fn ideal(&self) -> PyPathAlgebraIdeal {
        let kind = map_ref!(
            &self.kind,
            QuotientKind => IdealKind,
            |algebra| algebra.ideal().clone()
        );
        PyPathAlgebraIdeal { kind }
    }

    #[getter]
    fn relations(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| wrap_algebra_elements(
            algebra.relations().to_vec()
        ))
    }

    #[getter]
    fn groebner_basis(&self) -> Option<PyGroebnerBasis> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra
            .groebner_basis()
            .cloned()
            .map(wrap_groebner_basis))
    }

    #[getter]
    fn has_normal_form_function(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra
            .has_normal_form_function())
    }

    fn generator(&self, name: &str) -> PyResult<PyQuotientElement> {
        let path = dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra
            .quiver()
            .generator_by_name(name)
            .ok_or_else(|| invalid(format!("unknown quiver generator {name:?}")))?
            .path());
        let kind = map_ref!(
            &self.kind,
            QuotientKind => QuotientElementKind,
            |algebra| algebra.path(path).map_err(invalid)?
        );
        Ok(PyQuotientElement { kind })
    }

    #[getter]
    fn generators(&self) -> PyResult<Vec<PyQuotientElement>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| Ok(
            wrap_quotient_elements(algebra.generators().map_err(invalid)?)
        ))
    }

    fn element(&self, element: &PyAlgebraElement) -> PyResult<PyQuotientElement> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            QuotientKind => QuotientElementKind,
            |algebra| algebra
                .element(runtime::Extract::extract(element)?)
                .map_err(invalid)?
        );
        Ok(PyQuotientElement { kind })
    }

    fn path(&self, path: &PyPath) -> PyResult<PyQuotientElement> {
        let path = path.inner.clone();
        let kind = map_ref!(
            &self.kind,
            QuotientKind => QuotientElementKind,
            |algebra| algebra.path(path).map_err(invalid)?
        );
        Ok(PyQuotientElement { kind })
    }

    fn zero(&self) -> PyQuotientElement {
        let kind = map_ref!(
            &self.kind,
            QuotientKind => QuotientElementKind,
            |algebra| algebra.zero()
        );
        PyQuotientElement { kind }
    }

    fn one(&self) -> PyResult<PyQuotientElement> {
        let kind = map_ref!(
            &self.kind,
            QuotientKind => QuotientElementKind,
            |algebra| algebra.one().map_err(invalid)?
        );
        Ok(PyQuotientElement { kind })
    }

    fn dimension(&self) -> PyResult<Option<String>> {
        let dimension = dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::quotient_dimension(algebra).map_err(invalid)?
        });
        Ok(match dimension {
            AlgebraDimension::Finite(value) => Some(value.to_string()),
            AlgebraDimension::Infinite => None,
        })
    }

    #[getter]
    fn is_finite_dimensional(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::quotient_is_finite_dimensional(algebra).map_err(invalid)
        })
    }

    fn canonical_basis(&self) -> PyResult<Option<Vec<PyQuotientElement>>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| Ok(algebra
            .canonical_basis()
            .map_err(invalid)?
            .map(wrap_quotient_elements)))
    }

    #[getter]
    fn is_admissible(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::quotient_is_admissible(algebra).map_err(invalid)
        })
    }

    #[getter]
    fn is_monomial(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::quotient_is_monomial(algebra).map_err(invalid)
        })
    }

    #[getter]
    fn is_special_biserial(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::quotient_is_special_biserial(algebra).map_err(invalid)
        })
    }

    #[getter]
    fn is_string(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::quotient_is_string(algebra).map_err(invalid)
        })
    }

    #[getter]
    fn is_gentle(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::quotient_is_gentle(algebra).map_err(invalid)
        })
    }

    #[getter]
    fn is_triangular_reduced(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::quotient_is_triangular_reduced(algebra).map_err(invalid)
        })
    }

    #[getter]
    fn is_semisimple(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::admissible_quotient_is_semisimple(algebra).map_err(invalid)
        })
    }

    #[getter]
    fn is_elementary(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::admissible_quotient_is_elementary(algebra).map_err(invalid)
        })
    }

    #[getter]
    fn is_basic(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::admissible_quotient_is_basic(algebra).map_err(invalid)
        })
    }

    fn opposite(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientKind => QuotientKind,
            |algebra| qpa_gbnp::opposite_path_algebra(algebra).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            QuotientKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra.view_obj_with(
            &algebra
                .original_path_algebra()
                .field()
                .python_kind()
                .display()
        ))
    }
}

#[pymethods]
impl PyQuotientElement {
    #[getter]
    fn algebra(&self) -> PyQuotientAlgebra {
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientKind,
            |element| element.algebra().clone()
        );
        PyQuotientAlgebra { kind }
    }

    #[getter]
    fn representative(&self) -> PyAlgebraElement {
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => AlgebraElementKind,
            |element| element.representative().clone()
        );
        PyAlgebraElement { kind }
    }

    #[getter]
    fn external_representation(&self) -> Vec<(Vec<usize>, PyScalar)> {
        dispatch_ref!(&self.kind, QuotientElementKind, |element| element
            .external_representation()
            .terms()
            .iter()
            .map(|(path, coefficient)| (
                path.clone(),
                element
                    .algebra()
                    .original_path_algebra()
                    .field()
                    .wrap_scalar(coefficient),
            ))
            .collect())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientElementKind, |element| element.is_zero())
    }

    #[getter]
    fn is_normal_form(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientElementKind, |element| element
            .is_normal_form())
    }

    #[getter]
    fn is_left_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientElementKind, |element| element
            .is_left_uniform())
    }

    #[getter]
    fn is_right_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientElementKind, |element| element
            .is_right_uniform())
    }

    #[getter]
    fn is_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientElementKind, |element| element
            .is_uniform())
    }

    #[getter]
    fn leading_term(&self) -> Self {
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| element.leading_term()
        );
        Self { kind }
    }

    #[getter]
    fn leading_coefficient(&self) -> PyScalar {
        dispatch_ref!(&self.kind, QuotientElementKind, |element| element
            .algebra()
            .original_path_algebra()
            .field()
            .wrap_scalar(&element.leading_coefficient()))
    }

    #[getter]
    fn leading_monomial(&self) -> PyPath {
        let inner = dispatch_ref!(&self.kind, QuotientElementKind, |element| element
            .leading_monomial());
        PyPath { inner }
    }

    #[getter]
    fn vertex_position(&self) -> PyResult<usize> {
        dispatch_ref!(&self.kind, QuotientElementKind, |element| element
            .vertex_position()
            .map_err(invalid))
    }

    fn normalized(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| element.normalized().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            QuotientElementKind => QuotientElementKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("quotient elements belong to different fields")
        )?;
        Ok(Self { kind })
    }

    fn negate(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| element.neg().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn subtract(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            QuotientElementKind => QuotientElementKind,
            |left, right| left.sub(right).map_err(invalid)?,
            invalid("quotient elements belong to different fields")
        )?;
        Ok(Self { kind })
    }

    fn multiply(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            QuotientElementKind => QuotientElementKind,
            |left, right| left.multiply(right).map_err(invalid)?,
            invalid("quotient elements belong to different fields")
        )?;
        Ok(Self { kind })
    }

    fn multiply_path_left(&self, path: &PyPath) -> PyResult<Self> {
        let path = path.inner.clone();
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| element.multiply_path_left(&path).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn multiply_path_right(&self, path: &PyPath) -> PyResult<Self> {
        let path = path.inner.clone();
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| element.multiply_path_right(&path).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_left(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| element
                .scale_left(&scalar.value_for(element.algebra().original_path_algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_right(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| element
                .scale_right(&scalar.value_for(element.algebra().original_path_algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn pow(&self, exponent: u64) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| element.pow(exponent).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn pow_integer(&self, exponent: i64) -> PyResult<Option<Self>> {
        dispatch_ref!(&self.kind, QuotientElementKind, |element| Ok(
            qpa_gbnp::power_quotient_element(element, exponent)
                .map_err(invalid)?
                .map(runtime::Wrap::wrap)
                .map(|kind| Self { kind })
        ))
    }

    fn inverse(&self) -> PyResult<Option<Self>> {
        dispatch_ref!(&self.kind, QuotientElementKind, |element| Ok(
            qpa_gbnp::inverse_quotient_element(element)
                .map_err(invalid)?
                .map(runtime::Wrap::wrap)
                .map(|kind| Self { kind })
        ))
    }

    fn compare(&self, right: &Self) -> PyResult<i8> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            QuotientElementKind,
            |left, right| Ok(ordering_value(left.compare(right).map_err(invalid)?)),
            Err(invalid("quotient elements belong to different fields"))
        )
    }

    fn opposite(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| qpa_gbnp::opposite_quotient_element(element).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn __add__(&self, right: &Self) -> PyResult<Self> {
        self.add(right)
    }

    fn __sub__(&self, right: &Self) -> PyResult<Self> {
        self.subtract(right)
    }

    fn __neg__(&self) -> PyResult<Self> {
        self.negate()
    }

    fn __mul__(&self, right: &Self) -> PyResult<Self> {
        self.multiply(right)
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, QuotientElementKind, |element| element
            .print_obj_with(|value| element
                .algebra()
                .original_path_algebra()
                .field()
                .scalar_text(value)))
    }

    fn __eq__(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            QuotientElementKind,
            |left, right| left == right,
            false
        )
    }
}

#[pymethods]
impl PyMatrixModule {
    #[getter]
    fn dimension_vector(&self) -> Vec<usize> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .dimension_vector()
            .to_vec())
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module.dimension())
    }

    #[getter]
    fn arrow_matrix_shapes(&self) -> Vec<(usize, usize)> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| matrix_shapes(
            module.arrow_matrices()
        ))
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            MatrixModuleKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .qpa_print_string())
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
    module.add_class::<PyPathAlgebraIdeal>()?;
    module.add_class::<PyGroebnerBasis>()?;
    module.add_class::<PyQuotientAlgebra>()?;
    module.add_class::<PyQuotientElement>()?;
    module.add_class::<PyMatrixModule>()?;
    Ok(())
}
