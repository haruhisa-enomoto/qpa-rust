//! Python wrappers for QPA's finite-dimensional algebra conversions.

mod quotient;

use quotient::structure_constants_quotient_for;
pub(crate) use quotient::{PyStructureConstantsHomomorphism, PyStructureConstantsQuotient};

use numerica::domains::{
    finite_field::{Z2, Zp},
    rational::{Q, RationalField},
};
use pyo3::prelude::*;
use qpa_core::{
    AlgebraGlobalDimension, OppositeStructureConstantsAlgebra,
    OppositeStructureConstantsAlgebraElement, PathAlgebraModule, QuiverAlgebra,
    StructureConstantsAlgebra, StructureConstantsAlgebraElement, StructureConstantsSubalgebra,
};

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

#[derive(Clone)]
pub(super) enum StructureConstantsAlgebraKind {
    Rational(StructureConstantsAlgebra<RationalField>),
    Binary(StructureConstantsAlgebra<Z2>),
    Prime(StructureConstantsAlgebra<Zp>),
}

/// A finite-dimensional algebra in QPA's ordered structure-constant basis.
#[pyclass(name = "StructureConstantsAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyStructureConstantsAlgebra {
    pub(super) kind: StructureConstantsAlgebraKind,
}

#[derive(Clone)]
enum StructureConstantsElementKind {
    Rational(StructureConstantsAlgebraElement<RationalField>),
    Binary(StructureConstantsAlgebraElement<Z2>),
    Prime(StructureConstantsAlgebraElement<Zp>),
}

/// An element in the ordered basis of a structure-constants algebra.
#[pyclass(name = "StructureConstantsElement", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyStructureConstantsElement {
    kind: StructureConstantsElementKind,
}

#[derive(Clone)]
enum StructureConstantsSubalgebraKind {
    Rational(StructureConstantsSubalgebra<RationalField>),
    Binary(StructureConstantsSubalgebra<Z2>),
    Prime(StructureConstantsSubalgebra<Zp>),
}

/// A QPA subalgebra with an explicit ordered basis in its ambient algebra.
#[pyclass(name = "StructureConstantsSubalgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyStructureConstantsSubalgebra {
    kind: StructureConstantsSubalgebraKind,
}

#[derive(Clone)]
enum OppositeStructureConstantsAlgebraKind {
    Rational(OppositeStructureConstantsAlgebra<RationalField>),
    Binary(OppositeStructureConstantsAlgebra<Z2>),
    Prime(OppositeStructureConstantsAlgebra<Zp>),
}

/// QPA's generic opposite of a structure-constants algebra.
#[pyclass(
    name = "OppositeStructureConstantsAlgebra",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PyOppositeStructureConstantsAlgebra {
    kind: OppositeStructureConstantsAlgebraKind,
}

#[derive(Clone)]
enum OppositeStructureConstantsElementKind {
    Rational(OppositeStructureConstantsAlgebraElement<RationalField>),
    Binary(OppositeStructureConstantsAlgebraElement<Z2>),
    Prime(OppositeStructureConstantsAlgebraElement<Zp>),
}

/// An element of QPA's generic opposite-algebra family.
#[pyclass(
    name = "OppositeStructureConstantsElement",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PyOppositeStructureConstantsElement {
    kind: OppositeStructureConstantsElementKind,
}

#[derive(Clone)]
pub(super) enum QuiverAlgebraKind {
    Rational(QuiverAlgebra<RationalField>),
    Binary(QuiverAlgebra<Z2>),
    Prime(QuiverAlgebra<Zp>),
}

/// A QPA path algebra or quotient path algebra returned by a conversion.
#[pyclass(name = "QuiverAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuiverAlgebra {
    pub(super) kind: QuiverAlgebraKind,
}

#[derive(Clone)]
pub(super) enum AlgebraAsQuiverAlgebraResultKind {
    Rational(qpa_gbnp::AlgebraAsQuiverAlgebraResult<RationalField>),
    Binary(qpa_gbnp::AlgebraAsQuiverAlgebraResult<Z2>),
    Prime(qpa_gbnp::AlgebraAsQuiverAlgebraResult<Zp>),
}

/// QPA's converted quiver algebra and ordered generator images.
#[pyclass(name = "AlgebraAsQuiverAlgebraResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAlgebraAsQuiverAlgebraResult {
    pub(super) kind: AlgebraAsQuiverAlgebraResultKind,
}

#[derive(Clone)]
pub(super) enum AdrAlgebraOfAlgebraResultKind {
    Rational(qpa_gbnp::AdrAlgebraOfAlgebraResult<RationalField>),
    Binary(qpa_gbnp::AdrAlgebraOfAlgebraResult<Z2>),
    Prime(qpa_gbnp::AdrAlgebraOfAlgebraResult<Zp>),
}

/// QPA's ADR algebra and installed finite-global-dimension fact.
#[pyclass(name = "AdrAlgebraOfAlgebraResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAdrAlgebraOfAlgebraResult {
    pub(super) kind: AdrAlgebraOfAlgebraResultKind,
}

macro_rules! impl_conversion_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<StructureConstantsAlgebra<$field>> for StructureConstantsAlgebraKind {
            fn wrap(value: StructureConstantsAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<StructureConstantsAlgebraElement<$field>>
            for StructureConstantsElementKind
        {
            fn wrap(value: StructureConstantsAlgebraElement<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Extract<StructureConstantsAlgebraElement<$field>>
            for StructureConstantsElementKind
        {
            fn extract(self) -> PyResult<StructureConstantsAlgebraElement<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "structure-constants element belongs to another coefficient field",
                    )),
                }
            }
        }

        impl runtime::Wrap<StructureConstantsSubalgebra<$field>>
            for StructureConstantsSubalgebraKind
        {
            fn wrap(value: StructureConstantsSubalgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<OppositeStructureConstantsAlgebra<$field>>
            for OppositeStructureConstantsAlgebraKind
        {
            fn wrap(value: OppositeStructureConstantsAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<OppositeStructureConstantsAlgebraElement<$field>>
            for OppositeStructureConstantsElementKind
        {
            fn wrap(value: OppositeStructureConstantsAlgebraElement<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<QuiverAlgebra<$field>> for QuiverAlgebraKind {
            fn wrap(value: QuiverAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::AlgebraAsQuiverAlgebraResult<$field>>
            for AlgebraAsQuiverAlgebraResultKind
        {
            fn wrap(value: qpa_gbnp::AlgebraAsQuiverAlgebraResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::AdrAlgebraOfAlgebraResult<$field>>
            for AdrAlgebraOfAlgebraResultKind
        {
            fn wrap(value: qpa_gbnp::AdrAlgebraOfAlgebraResult<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_conversion_runtime!(RationalField, Rational);
impl_conversion_runtime!(Z2, Binary);
impl_conversion_runtime!(Zp, Prime);

fn structure_constants_for<F: runtime::PythonField>(
    py: Python<'_>,
    field: F,
    products: Vec<Vec<Vec<Py<PyScalar>>>>,
) -> PyResult<PyStructureConstantsAlgebra>
where
    StructureConstantsAlgebraKind: runtime::Wrap<StructureConstantsAlgebra<F>>,
{
    let dimension = products.len();
    if products.iter().any(|left| {
        left.len() != dimension || left.iter().any(|product| product.len() != dimension)
    }) {
        return Err(invalid(format!(
            "structure constants must form a {dimension} x {dimension} x {dimension} tensor"
        )));
    }
    let products = products
        .into_iter()
        .flatten()
        .flatten()
        .map(|coefficient| field.scalar_value(&coefficient.borrow(py)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyStructureConstantsAlgebra {
        kind: runtime::Wrap::wrap(
            StructureConstantsAlgebra::new(field, dimension, products).map_err(invalid)?,
        ),
    })
}

fn structure_element_for<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: &StructureConstantsAlgebra<F>,
    coordinates: Vec<Py<PyScalar>>,
) -> PyResult<PyStructureConstantsElement>
where
    StructureConstantsElementKind: runtime::Wrap<StructureConstantsAlgebraElement<F>>,
{
    let coordinates = coordinates
        .into_iter()
        .map(|coefficient| algebra.field().scalar_value(&coefficient.borrow(py)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyStructureConstantsElement {
        kind: runtime::Wrap::wrap(algebra.element(coordinates).map_err(invalid)?),
    })
}

fn opposite_element_for<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: &OppositeStructureConstantsAlgebra<F>,
    coordinates: Vec<Py<PyScalar>>,
) -> PyResult<PyOppositeStructureConstantsElement>
where
    OppositeStructureConstantsElementKind:
        runtime::Wrap<OppositeStructureConstantsAlgebraElement<F>>,
{
    let coordinates = coordinates
        .into_iter()
        .map(|coefficient| algebra.field().scalar_value(&coefficient.borrow(py)))
        .collect::<PyResult<Vec<_>>>()?;
    Ok(wrap_opposite_element(
        algebra.element(coordinates).map_err(invalid)?,
    ))
}

fn wrap_structure_element<F: runtime::PythonField>(
    element: StructureConstantsAlgebraElement<F>,
) -> PyStructureConstantsElement
where
    StructureConstantsElementKind: runtime::Wrap<StructureConstantsAlgebraElement<F>>,
{
    PyStructureConstantsElement {
        kind: runtime::Wrap::wrap(element),
    }
}

fn wrap_structure_elements<F: runtime::PythonField>(
    elements: Vec<StructureConstantsAlgebraElement<F>>,
) -> Vec<PyStructureConstantsElement>
where
    StructureConstantsElementKind: runtime::Wrap<StructureConstantsAlgebraElement<F>>,
{
    elements.into_iter().map(wrap_structure_element).collect()
}

fn wrap_structure_subalgebra<F: runtime::PythonField>(
    subalgebra: StructureConstantsSubalgebra<F>,
) -> PyStructureConstantsSubalgebra
where
    StructureConstantsSubalgebraKind: runtime::Wrap<StructureConstantsSubalgebra<F>>,
{
    PyStructureConstantsSubalgebra {
        kind: runtime::Wrap::wrap(subalgebra),
    }
}

fn wrap_opposite_algebra<F: runtime::PythonField>(
    algebra: OppositeStructureConstantsAlgebra<F>,
) -> PyOppositeStructureConstantsAlgebra
where
    OppositeStructureConstantsAlgebraKind: runtime::Wrap<OppositeStructureConstantsAlgebra<F>>,
{
    PyOppositeStructureConstantsAlgebra {
        kind: runtime::Wrap::wrap(algebra),
    }
}

fn wrap_opposite_element<F: runtime::PythonField>(
    element: OppositeStructureConstantsAlgebraElement<F>,
) -> PyOppositeStructureConstantsElement
where
    OppositeStructureConstantsElementKind:
        runtime::Wrap<OppositeStructureConstantsAlgebraElement<F>>,
{
    PyOppositeStructureConstantsElement {
        kind: runtime::Wrap::wrap(element),
    }
}

fn wrap_opposite_elements<F: runtime::PythonField>(
    elements: Vec<OppositeStructureConstantsAlgebraElement<F>>,
) -> Vec<PyOppositeStructureConstantsElement>
where
    OppositeStructureConstantsElementKind:
        runtime::Wrap<OppositeStructureConstantsAlgebraElement<F>>,
{
    elements.into_iter().map(wrap_opposite_element).collect()
}

fn structure_matrix_rows<F: runtime::PythonField>(matrix: &DenseMatrix<F>) -> Vec<Vec<PyScalar>> {
    (0..matrix.nrows())
        .map(|row| {
            matrix
                .row(row)
                .iter()
                .map(|coefficient| matrix.field().wrap_scalar(coefficient))
                .collect()
        })
        .collect()
}

pub(super) fn wrap_quiver_algebra<F: runtime::PythonField>(
    algebra: QuiverAlgebra<F>,
) -> PyQuiverAlgebra
where
    QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    PyQuiverAlgebra {
        kind: runtime::Wrap::wrap(algebra),
    }
}

pub(super) fn adr_algebra_of_basic_module_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
) -> PyResult<PyQuiverAlgebra>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>>,
    QuiverAlgebraKind: runtime::Wrap<QuiverAlgebra<F>>,
{
    Ok(wrap_quiver_algebra(
        qpa_gbnp::adr_algebra_of_basic_module(&runtime::extract_all(modules)?).map_err(invalid)?,
    ))
}

#[pymethods]
impl PyStructureConstantsAlgebra {
    #[new]
    fn new(
        py: Python<'_>,
        field: PyRef<'_, PyField>,
        products: Vec<Vec<Vec<Py<PyScalar>>>>,
    ) -> PyResult<Self> {
        match field.kind {
            FieldKind::Rational => structure_constants_for(py, Q, products),
            FieldKind::Binary => structure_constants_for(py, Z2, products),
            FieldKind::Prime(characteristic) => {
                structure_constants_for(py, Zp::new(characteristic), products)
            }
        }
    }

    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            algebra.field().python_kind()
        });
        PyField { kind }
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            algebra.dimension()
        })
    }

    #[getter]
    fn multiplication_tensor(&self) -> Vec<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            (0..algebra.dimension())
                .map(|left| {
                    (0..algebra.dimension())
                        .map(|right| {
                            algebra
                                .basis_product(left, right)
                                .iter()
                                .map(|coefficient| algebra.field().wrap_scalar(coefficient))
                                .collect()
                        })
                        .collect()
                })
                .collect()
        })
    }

    fn zero(&self) -> PyStructureConstantsElement {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            wrap_structure_element(algebra.zero())
        })
    }

    #[getter]
    fn basis_elements(&self) -> Vec<PyStructureConstantsElement> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            wrap_structure_elements(algebra.basis_elements())
        })
    }

    fn element(
        &self,
        py: Python<'_>,
        coordinates: Vec<Py<PyScalar>>,
    ) -> PyResult<PyStructureConstantsElement> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            structure_element_for(py, algebra, coordinates)
        })
    }

    fn bind_named_element(
        &self,
        name: &str,
        element: &PyStructureConstantsElement,
    ) -> PyResult<()> {
        match (&self.kind, &element.kind) {
            (
                StructureConstantsAlgebraKind::Rational(algebra),
                StructureConstantsElementKind::Rational(element),
            ) => algebra.bind_named_element(name, element).map_err(invalid),
            (
                StructureConstantsAlgebraKind::Binary(algebra),
                StructureConstantsElementKind::Binary(element),
            ) => algebra.bind_named_element(name, element).map_err(invalid),
            (
                StructureConstantsAlgebraKind::Prime(algebra),
                StructureConstantsElementKind::Prime(element),
            ) => algebra.bind_named_element(name, element).map_err(invalid),
            _ => Err(invalid(
                "structure-constants algebra and element use different coefficient fields",
            )),
        }
    }

    fn opposite(&self) -> PyResult<PyOppositeStructureConstantsAlgebra> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            Ok(wrap_opposite_algebra(algebra.opposite().map_err(invalid)?))
        })
    }

    fn quotient_by_ideal_basis(
        &self,
        py: Python<'_>,
        ideal_basis: Vec<Vec<Py<PyScalar>>>,
    ) -> PyResult<PyStructureConstantsQuotient> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            structure_constants_quotient_for(py, algebra, ideal_basis)
        })
    }

    fn find_multiplicative_identity(&self) -> PyResult<Option<PyStructureConstantsElement>> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            Ok(algebra
                .find_multiplicative_identity()
                .map_err(invalid)?
                .map(wrap_structure_element))
        })
    }

    fn center(&self) -> PyResult<PyStructureConstantsSubalgebra> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            Ok(wrap_structure_subalgebra(
                algebra.centre().map_err(invalid)?,
            ))
        })
    }

    fn radical(&self) -> PyResult<PyStructureConstantsSubalgebra> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            Ok(wrap_structure_subalgebra(
                algebra.radical().map_err(invalid)?,
            ))
        })
    }

    #[getter]
    fn radical_series(&self) -> PyResult<Vec<Vec<Vec<PyScalar>>>> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            Ok(algebra
                .radical_series()
                .map_err(invalid)?
                .iter()
                .map(structure_matrix_rows)
                .collect())
        })
    }

    fn power_subalgebra(&self, degree: usize) -> PyResult<PyStructureConstantsSubalgebra> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            Ok(wrap_structure_subalgebra(
                algebra.power_subalgebra(degree).map_err(invalid)?,
            ))
        })
    }

    fn blocks(&self) -> PyResult<Vec<PyStructureConstantsSubalgebra>> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            Ok(algebra
                .blocks()
                .map_err(invalid)?
                .into_iter()
                .map(wrap_structure_subalgebra)
                .collect())
        })
    }

    #[getter]
    fn is_semisimple(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            algebra.is_semisimple_algebra().map_err(invalid)
        })
    }

    #[getter]
    fn is_radical_square_zero(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            algebra.is_radical_square_zero().map_err(invalid)
        })
    }

    #[getter]
    fn is_basic(&self) -> PyResult<Option<bool>> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            algebra.is_basic_algebra().map_err(invalid)
        })
    }

    #[getter]
    fn is_elementary(&self) -> PyResult<Option<bool>> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            algebra.is_elementary_algebra().map_err(invalid)
        })
    }

    #[getter]
    fn cached_global_dimension(&self) -> Option<String> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            algebra
                .cached_global_dimension()
                .map(|dimension| match dimension {
                    AlgebraGlobalDimension::Finite(dimension) => dimension.to_string(),
                    AlgebraGlobalDimension::Infinite => "infinity".into(),
                })
        })
    }

    fn global_dimension_from_known_properties(&self) -> Option<String> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            algebra
                .global_dimension_from_known_properties()
                .map(|dimension| match dimension {
                    AlgebraGlobalDimension::Finite(dimension) => dimension.to_string(),
                    AlgebraGlobalDimension::Infinite => "infinity".into(),
                })
        })
    }

    fn is_finite_global_dimension_from_known_properties(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            algebra.is_finite_global_dimension_from_known_properties()
        })
    }

    fn primitive_idempotents(&self) -> PyResult<Vec<PyStructureConstantsElement>> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            Ok(wrap_structure_elements(
                algebra.primitive_idempotents().map_err(invalid)?,
            ))
        })
    }

    fn idempotents_for_decomposition(&self) -> PyResult<Vec<PyStructureConstantsElement>> {
        dispatch_ref!(&self.kind, StructureConstantsAlgebraKind, |algebra| {
            Ok(wrap_structure_elements(
                algebra.idempotents_for_decomposition().map_err(invalid)?,
            ))
        })
    }

    fn as_quiver_algebra(&self) -> PyResult<PyAlgebraAsQuiverAlgebraResult> {
        let kind = map_ref!(
            &self.kind,
            StructureConstantsAlgebraKind => AlgebraAsQuiverAlgebraResultKind,
            |algebra| qpa_gbnp::algebra_as_quiver_algebra(algebra).map_err(invalid)?
        );
        Ok(PyAlgebraAsQuiverAlgebraResult { kind })
    }
}

#[pymethods]
impl PyOppositeStructureConstantsAlgebra {
    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| algebra.field().python_kind()
        );
        PyField { kind }
    }

    #[getter]
    fn underlying_algebra(&self) -> PyStructureConstantsAlgebra {
        let kind = map_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind => StructureConstantsAlgebraKind,
            |algebra| algebra.underlying_algebra().clone()
        );
        PyStructureConstantsAlgebra { kind }
    }

    #[getter]
    fn multiplication_algebra(&self) -> PyStructureConstantsAlgebra {
        let kind = map_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind => StructureConstantsAlgebraKind,
            |algebra| algebra.as_structure_constants_algebra().clone()
        );
        PyStructureConstantsAlgebra { kind }
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| algebra.dimension()
        )
    }

    #[getter]
    fn is_opposite_algebra(&self) -> bool {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| algebra.is_opposite_algebra()
        )
    }

    #[getter]
    fn is_identical_to_underlying(&self) -> bool {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| algebra.is_identical_to_underlying()
        )
    }

    #[getter]
    fn is_finite_dimensional(&self) -> bool {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| algebra.is_finite_dimensional()
        )
    }

    #[getter]
    fn basis_elements(&self) -> Vec<PyOppositeStructureConstantsElement> {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| wrap_opposite_elements(algebra.basis_elements())
        )
    }

    #[getter]
    fn canonical_basis_elements(&self) -> Vec<PyOppositeStructureConstantsElement> {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| wrap_opposite_elements(algebra.canonical_basis_elements())
        )
    }

    fn element(
        &self,
        py: Python<'_>,
        coordinates: Vec<Py<PyScalar>>,
    ) -> PyResult<PyOppositeStructureConstantsElement> {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| opposite_element_for(py, algebra, coordinates)
        )
    }

    fn element_from_underlying(
        &self,
        element: &PyStructureConstantsElement,
    ) -> PyResult<PyOppositeStructureConstantsElement> {
        match (&self.kind, &element.kind) {
            (
                OppositeStructureConstantsAlgebraKind::Rational(algebra),
                StructureConstantsElementKind::Rational(element),
            ) => Ok(wrap_opposite_element(
                algebra
                    .element_from_underlying(element.clone())
                    .map_err(invalid)?,
            )),
            (
                OppositeStructureConstantsAlgebraKind::Binary(algebra),
                StructureConstantsElementKind::Binary(element),
            ) => Ok(wrap_opposite_element(
                algebra
                    .element_from_underlying(element.clone())
                    .map_err(invalid)?,
            )),
            (
                OppositeStructureConstantsAlgebraKind::Prime(algebra),
                StructureConstantsElementKind::Prime(element),
            ) => Ok(wrap_opposite_element(
                algebra
                    .element_from_underlying(element.clone())
                    .map_err(invalid)?,
            )),
            _ => Err(invalid(
                "opposite algebra and underlying element use different coefficient fields",
            )),
        }
    }

    fn element_by_name(&self, name: &str) -> PyResult<Option<PyOppositeStructureConstantsElement>> {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| Ok(algebra
                .element_by_name(name)
                .map_err(invalid)?
                .map(wrap_opposite_element))
        )
    }

    fn coefficients(
        &self,
        element: &PyOppositeStructureConstantsElement,
    ) -> PyResult<Vec<PyScalar>> {
        match (&self.kind, &element.kind) {
            (
                OppositeStructureConstantsAlgebraKind::Rational(algebra),
                OppositeStructureConstantsElementKind::Rational(element),
            ) => Ok(algebra
                .coefficients(element)
                .map_err(invalid)?
                .iter()
                .map(|coefficient| algebra.field().wrap_scalar(coefficient))
                .collect()),
            (
                OppositeStructureConstantsAlgebraKind::Binary(algebra),
                OppositeStructureConstantsElementKind::Binary(element),
            ) => Ok(algebra
                .coefficients(element)
                .map_err(invalid)?
                .iter()
                .map(|coefficient| algebra.field().wrap_scalar(coefficient))
                .collect()),
            (
                OppositeStructureConstantsAlgebraKind::Prime(algebra),
                OppositeStructureConstantsElementKind::Prime(element),
            ) => Ok(algebra
                .coefficients(element)
                .map_err(invalid)?
                .iter()
                .map(|coefficient| algebra.field().wrap_scalar(coefficient))
                .collect()),
            _ => Err(invalid(
                "opposite algebra and element use different coefficient fields",
            )),
        }
    }

    fn zero(&self) -> PyOppositeStructureConstantsElement {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| wrap_opposite_element(algebra.zero())
        )
    }

    fn one(&self) -> PyResult<Option<PyOppositeStructureConstantsElement>> {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| Ok(algebra.one().map_err(invalid)?.map(wrap_opposite_element))
        )
    }

    fn contains(&self, element: &PyOppositeStructureConstantsElement) -> bool {
        match (&self.kind, &element.kind) {
            (
                OppositeStructureConstantsAlgebraKind::Rational(algebra),
                OppositeStructureConstantsElementKind::Rational(element),
            ) => algebra.contains(element),
            (
                OppositeStructureConstantsAlgebraKind::Binary(algebra),
                OppositeStructureConstantsElementKind::Binary(element),
            ) => algebra.contains(element),
            (
                OppositeStructureConstantsAlgebraKind::Prime(algebra),
                OppositeStructureConstantsElementKind::Prime(element),
            ) => algebra.contains(element),
            _ => false,
        }
    }

    fn __contains__(&self, element: &PyOppositeStructureConstantsElement) -> bool {
        self.contains(element)
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsAlgebraKind,
            |algebra| algebra.to_string()
        )
    }
}

#[pymethods]
impl PyOppositeStructureConstantsElement {
    #[getter]
    fn algebra(&self) -> PyOppositeStructureConstantsAlgebra {
        let kind = map_ref!(
            &self.kind,
            OppositeStructureConstantsElementKind => OppositeStructureConstantsAlgebraKind,
            |element| element.algebra().clone()
        );
        PyOppositeStructureConstantsAlgebra { kind }
    }

    #[getter]
    fn underlying_element(&self) -> PyStructureConstantsElement {
        let kind = map_ref!(
            &self.kind,
            OppositeStructureConstantsElementKind => StructureConstantsElementKind,
            |element| element.underlying_element().clone()
        );
        PyStructureConstantsElement { kind }
    }

    #[getter]
    fn coordinates(&self) -> Vec<PyScalar> {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsElementKind,
            |element| element
                .coordinates()
                .iter()
                .map(|coefficient| element.algebra().field().wrap_scalar(coefficient))
                .collect()
        )
    }

    #[getter]
    fn external_representation(&self) -> Vec<PyScalar> {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsElementKind,
            |element| element
                .external_representation()
                .iter()
                .map(|coefficient| element.algebra().field().wrap_scalar(coefficient))
                .collect()
        )
    }

    #[getter]
    fn is_opposite_algebra_element(&self) -> bool {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsElementKind,
            |element| element.is_opposite_algebra_element()
        )
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            OppositeStructureConstantsElementKind => OppositeStructureConstantsElementKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("opposite elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn multiply(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            OppositeStructureConstantsElementKind => OppositeStructureConstantsElementKind,
            |left, right| left.multiply(right).map_err(invalid)?,
            invalid("opposite elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn scale(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            OppositeStructureConstantsElementKind => OppositeStructureConstantsElementKind,
            |element| element
                .scale(&element.algebra().field().scalar_value(scalar)?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn negate(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            OppositeStructureConstantsElementKind => OppositeStructureConstantsElementKind,
            |element| element.negate().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn zero_same_parent(&self) -> Self {
        let kind = map_ref!(
            &self.kind,
            OppositeStructureConstantsElementKind => OppositeStructureConstantsElementKind,
            |element| element.zero_same_parent()
        );
        Self { kind }
    }

    fn compare(&self, right: &Self) -> PyResult<i8> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            OppositeStructureConstantsElementKind,
            |left, right| Ok(ordering_value(left.compare(right).map_err(invalid)?)),
            Err(invalid(
                "opposite elements belong to different coefficient fields"
            ))
        )
    }

    fn __add__(&self, right: &Self) -> PyResult<Self> {
        self.add(right)
    }

    fn __mul__(&self, right: &Self) -> PyResult<Self> {
        self.multiply(right)
    }

    fn __neg__(&self) -> PyResult<Self> {
        self.negate()
    }

    fn __lt__(&self, right: &Self) -> PyResult<bool> {
        Ok(self.compare(right)? < 0)
    }

    fn __eq__(&self, right: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            OppositeStructureConstantsElementKind,
            |left, right| left == right,
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(
            &self.kind,
            OppositeStructureConstantsElementKind,
            |element| element.to_string()
        )
    }
}

#[pymethods]
impl PyStructureConstantsElement {
    #[getter]
    fn algebra(&self) -> PyStructureConstantsAlgebra {
        let kind = map_ref!(
            &self.kind,
            StructureConstantsElementKind => StructureConstantsAlgebraKind,
            |element| element.algebra().clone()
        );
        PyStructureConstantsAlgebra { kind }
    }

    #[getter]
    fn coordinates(&self) -> Vec<PyScalar> {
        dispatch_ref!(&self.kind, StructureConstantsElementKind, |element| {
            element
                .coordinates()
                .iter()
                .map(|coefficient| element.algebra().field().wrap_scalar(coefficient))
                .collect()
        })
    }

    fn multiply(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            StructureConstantsElementKind => StructureConstantsElementKind,
            |left, right| left.multiply(right).map_err(invalid)?,
            invalid("structure-constants elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn __mul__(&self, right: &Self) -> PyResult<Self> {
        self.multiply(right)
    }

    fn __eq__(&self, right: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            StructureConstantsElementKind,
            |left, right| left == right,
            false
        )
    }
}

#[pymethods]
impl PyStructureConstantsSubalgebra {
    #[getter]
    fn ambient(&self) -> PyStructureConstantsAlgebra {
        let kind = map_ref!(
            &self.kind,
            StructureConstantsSubalgebraKind => StructureConstantsAlgebraKind,
            |subalgebra| subalgebra.ambient().clone()
        );
        PyStructureConstantsAlgebra { kind }
    }

    #[getter]
    fn basis(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, StructureConstantsSubalgebraKind, |subalgebra| {
            structure_matrix_rows(subalgebra.basis())
        })
    }

    #[getter]
    fn algebra(&self) -> PyStructureConstantsAlgebra {
        let kind = map_ref!(
            &self.kind,
            StructureConstantsSubalgebraKind => StructureConstantsAlgebraKind,
            |subalgebra| subalgebra.algebra().clone()
        );
        PyStructureConstantsAlgebra { kind }
    }

    fn lift_coordinates(
        &self,
        py: Python<'_>,
        coordinates: Vec<Py<PyScalar>>,
    ) -> PyResult<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, StructureConstantsSubalgebraKind, |subalgebra| {
            let coordinates = coordinates
                .into_iter()
                .map(|coefficient| {
                    subalgebra
                        .ambient()
                        .field()
                        .scalar_value(&coefficient.borrow(py))
                })
                .collect::<PyResult<Vec<_>>>()?;
            Ok(subalgebra
                .lift_coordinates(&coordinates)
                .map_err(invalid)?
                .iter()
                .map(|coefficient| subalgebra.ambient().field().wrap_scalar(coefficient))
                .collect())
        })
    }
}

#[pymethods]
impl PyQuiverAlgebra {
    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            algebra.field().python_kind()
        });
        PyField { kind }
    }

    #[getter]
    fn quiver(&self) -> PyQuiver {
        let inner = dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            algebra.quiver().clone()
        });
        PyQuiver { inner }
    }

    #[getter]
    fn is_quotient(&self) -> bool {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| {
            matches!(algebra, QuiverAlgebra::Quotient(_))
        })
    }

    #[getter]
    fn path_algebra(&self) -> Option<PyPathAlgebra> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| match algebra {
            QuiverAlgebra::Path(algebra) => Some(PyPathAlgebra {
                kind: runtime::Wrap::wrap(algebra.clone()),
            }),
            QuiverAlgebra::Quotient(_) => None,
        })
    }

    #[getter]
    fn quotient_algebra(&self) -> Option<PyQuotientAlgebra> {
        dispatch_ref!(&self.kind, QuiverAlgebraKind, |algebra| match algebra {
            QuiverAlgebra::Path(_) => None,
            QuiverAlgebra::Quotient(algebra) => Some(PyQuotientAlgebra {
                kind: runtime::Wrap::wrap(algebra.clone()),
            }),
        })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            QuiverAlgebraKind,
            |left, right| left.is_identical(right),
            false
        )
    }
}

#[pymethods]
impl PyAlgebraAsQuiverAlgebraResult {
    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, AlgebraAsQuiverAlgebraResultKind, |result| {
            wrap_quiver_algebra(result.algebra().clone())
        })
    }

    #[getter]
    fn generator_images(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, AlgebraAsQuiverAlgebraResultKind, |result| {
            result
                .generator_images()
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|coefficient| result.algebra().field().wrap_scalar(coefficient))
                        .collect()
                })
                .collect()
        })
    }

    #[getter]
    fn loewy_length(&self) -> usize {
        dispatch_ref!(&self.kind, AlgebraAsQuiverAlgebraResultKind, |result| {
            result.loewy_length()
        })
    }
}

#[pymethods]
impl PyAdrAlgebraOfAlgebraResult {
    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, AdrAlgebraOfAlgebraResultKind, |result| {
            wrap_quiver_algebra(result.algebra().clone())
        })
    }

    #[getter]
    fn is_finite_global_dimension_algebra(&self) -> bool {
        dispatch_ref!(&self.kind, AdrAlgebraOfAlgebraResultKind, |result| {
            result.is_finite_global_dimension_algebra()
        })
    }
}
