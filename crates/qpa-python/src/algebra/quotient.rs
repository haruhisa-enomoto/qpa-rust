//! Methods on quotient path algebras and their elements.

use pyo3::pymethods;

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

fn homomorphism_of_projectives_from_quotient_matrix<F: runtime::PythonField>(
    algebra: &QuotientPathAlgebra<F>,
    matrix: Vec<Vec<QuotientElementKind>>,
    source_vertices: &[usize],
    range_vertices: &[usize],
) -> PyResult<PyModuleHomomorphism>
where
    QuotientElementKind: runtime::Extract<QuotientPathAlgebraElement<F>>,
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    let matrix = matrix
        .into_iter()
        .map(|row| {
            Ok(runtime::extract_all(row)?
                .into_iter()
                .map(QuiverAlgebraElement::Quotient)
                .collect::<Vec<_>>())
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyModuleHomomorphism {
        kind: runtime::Wrap::wrap(
            qpa_gbnp::from_matrix_to_homomorphism_of_projectives(
                &QuiverAlgebra::Quotient(algebra.clone()),
                &matrix,
                source_vertices,
                range_vertices,
            )
            .map_err(invalid)?,
        ),
    })
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
                qpa_gbnp::truncated_path_algebra(Q, quiver.rust_quiver().clone(), exponent)
                    .map_err(invalid)?,
            ),
            FieldKind::Binary => QuotientKind::Binary(
                qpa_gbnp::truncated_path_algebra(Z2, quiver.rust_quiver().clone(), exponent)
                    .map_err(invalid)?,
            ),
            FieldKind::Prime(characteristic) => QuotientKind::Prime(
                qpa_gbnp::truncated_path_algebra(
                    Zp::new(characteristic),
                    quiver.rust_quiver().clone(),
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
        let ordering = dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra
            .ordering()
            .clone());
        PyQuiver::from_inner_with_ordering(inner, ordering)
    }

    #[getter]
    fn ordering(&self) -> PyQuiverOrdering {
        let inner = dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra
            .ordering()
            .clone());
        PyQuiverOrdering::from_inner(inner)
    }

    fn adr_algebra(&self) -> PyResult<PyAdrAlgebraOfAlgebraResult> {
        let kind = map_ref!(
            &self.kind,
            QuotientKind => AdrAlgebraOfAlgebraResultKind,
            |algebra| qpa_gbnp::adr_algebra_of_algebra(&QuiverAlgebra::Quotient(algebra.clone()))
                .map_err(invalid)?
        );
        Ok(PyAdrAlgebraOfAlgebraResult { kind })
    }

    fn global_dimension_bounded(&self, bound: usize) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| Ok(
            global_dimension_text(
                qpa_gbnp::global_dimension_of_algebra(
                    &QuiverAlgebra::Quotient(algebra.clone()),
                    bound,
                )
                .map_err(invalid)?,
            )
        ))
    }

    #[getter]
    fn cached_gorenstein_dimension(&self) -> Option<usize> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::gorenstein_dimension(&QuiverAlgebra::Quotient(algebra.clone()))
        })
    }

    fn gorenstein_dimension_bounded(&self, bound: i64) -> PyResult<Option<usize>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::gorenstein_dimension_of_algebra(
                &QuiverAlgebra::Quotient(algebra.clone()),
                &bound,
            )
            .map_err(invalid)
        })
    }

    fn dominant_dimension_bounded(&self, bound: i64) -> PyResult<Option<String>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, QuotientKind, |algebra| Ok(
            dominant_dimension_text(
                qpa_gbnp::dominant_dimension_of_algebra(
                    &QuiverAlgebra::Quotient(algebra.clone()),
                    &bound,
                )
                .map_err(invalid)?,
            )
        ))
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

    #[getter]
    fn is_gorenstein(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::gentle_algebra_is_gorenstein(&QuiverAlgebra::Quotient(algebra.clone()))
                .map_err(invalid)
        })
    }

    #[getter]
    fn is_hereditary(&self) -> PyResult<Option<bool>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::is_hereditary_algebra(&QuiverAlgebra::Quotient(algebra.clone()))
                .map_err(invalid)
        })
    }

    #[getter]
    fn is_finite_global_dimension(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::is_finite_global_dimension_from_known_properties(&QuiverAlgebra::Quotient(
                algebra.clone(),
            ))
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

    fn module_from_integer_matrices(
        &self,
        dimension_vector: Vec<usize>,
        actions: Vec<(String, Vec<Vec<i64>>)>,
    ) -> PyResult<PyMatrixModule> {
        let kind = map_ref!(
            &self.kind,
            QuotientKind => MatrixModuleKind,
            |algebra| build_matrix_module(
                QuiverAlgebra::Quotient(algebra.clone()),
                dimension_vector,
                actions,
            )?
        );
        Ok(PyMatrixModule { kind })
    }

    fn homomorphism_of_projectives_from_matrix(
        &self,
        py: Python<'_>,
        matrix: Vec<Vec<Py<PyQuotientElement>>>,
        source_vertices: Vec<usize>,
        range_vertices: Vec<usize>,
    ) -> PyResult<PyModuleHomomorphism> {
        let matrix = matrix
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|element| element.borrow(py).kind.clone())
                    .collect()
            })
            .collect();
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            homomorphism_of_projectives_from_quotient_matrix(
                algebra,
                matrix,
                &source_vertices,
                &range_vertices,
            )
        })
    }

    fn zero_module(&self) -> PyResult<PyMatrixModule> {
        let kind = map_ref!(
            &self.kind,
            QuotientKind => MatrixModuleKind,
            |algebra| qpa_core::zero_module(&QuiverAlgebra::Quotient(algebra.clone()))
                .map_err(invalid)?
        );
        Ok(PyMatrixModule { kind })
    }

    fn simple_modules(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| Ok(wrap_matrix_modules(
            qpa_gbnp::simple_modules(&QuiverAlgebra::Quotient(algebra.clone())).map_err(invalid)?,
        )))
    }

    fn indecomposable_projectives(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| Ok(wrap_matrix_modules(
            qpa_gbnp::indec_projective_modules(&QuiverAlgebra::Quotient(algebra.clone()))
                .map_err(invalid)?,
        )))
    }

    fn indecomposable_injectives(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| Ok(wrap_matrix_modules(
            qpa_gbnp::indec_injective_modules(&QuiverAlgebra::Quotient(algebra.clone()))
                .map_err(invalid)?,
        )))
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
