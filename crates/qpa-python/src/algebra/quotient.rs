//! Methods on quotient path algebras and their elements.

use pyo3::pymethods;

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

#[derive(Clone)]
enum QuotientCanonicalBasisKind {
    Rational(QuotientCanonicalBasis<RationalField>),
    Binary(QuotientCanonicalBasis<Z2>),
    Prime(QuotientCanonicalBasis<Zp>),
    Extension(QuotientCanonicalBasis<GaloisField>),
}

/// QPA's cached canonical basis object for a quotient path algebra.
#[pyclass(name = "QuotientCanonicalBasis", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuotientCanonicalBasis {
    kind: QuotientCanonicalBasisKind,
}

macro_rules! impl_quotient_basis_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<QuotientCanonicalBasis<$field>> for QuotientCanonicalBasisKind {
            fn wrap(value: QuotientCanonicalBasis<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_quotient_basis_runtime!(RationalField, Rational);
impl_quotient_basis_runtime!(Z2, Binary);
impl_quotient_basis_runtime!(Zp, Prime);
impl_quotient_basis_runtime!(GaloisField, Extension);

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
    #[pyo3(signature = (ideal, presentation_map, term=None))]
    fn projective_resolution_fp_path_algebra_module(
        &self,
        py: Python<'_>,
        ideal: &PyPathAlgebraIdeal,
        presentation_map: Vec<Vec<Py<PyQuotientElement>>>,
        term: Option<usize>,
    ) -> PyResult<PyGszProjectiveResolution> {
        module::explicit_gsz_resolution(self, py, ideal, presentation_map, term)
    }

    fn right_projective_module(
        &self,
        vertices: Vec<usize>,
    ) -> PyResult<PyRightProjectiveQuotientModule> {
        quotient_projective::right_projective_module(self, vertices)
    }

    fn complexity(&self, degree: i64) -> PyResult<usize> {
        let degree = Integer::from(degree);
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::complexity_of_algebra(&QuiverAlgebra::Quotient(algebra.clone()), &degree)
                .map_err(invalid)
        })
    }

    fn centre(&self) -> PyResult<PyQuotientPathAlgebraCentre> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            finite_structures::quotient_centre_for(algebra)
        })
    }

    fn units(&self) -> PyResult<PyQuiverAlgebraUnitGroup> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            finite_structures::units_for(&QuiverAlgebra::Quotient(algebra.clone()))
        })
    }

    #[getter]
    fn has_elements_of_units(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            QuiverAlgebra::Quotient(algebra.clone()).has_elements_of_units()
        })
    }

    fn elements_of_units(&self) -> PyResult<Vec<PyQuiverAlgebraElement>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            finite_structures::elements_of_units_for(&QuiverAlgebra::Quotient(algebra.clone()))
        })
    }

    /// Construct QPA's nice-basis conversion data for a quotient-element space.
    fn nice_basis(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyQuotientElement>>,
    ) -> PyResult<PyQuotientNiceBasis> {
        let generators = generators
            .into_iter()
            .map(|element| element.borrow(py).kind.clone())
            .collect();
        dispatch_ref!(&self.kind, QuotientKind, |quotient| {
            nice_basis::nice_basis_for(quotient, generators)
        })
    }

    fn left_support_of_quiver_algebra_element(
        &self,
        element: &PyQuotientElement,
    ) -> PyResult<Vec<usize>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            Ok(super::element_utilities::quotient_support_for(
                algebra,
                element.kind.clone(),
                super::element_utilities::SupportSide::Left,
            )?
            .pop()
            .unwrap_or_default())
        })
    }

    fn right_support_of_quiver_algebra_element(
        &self,
        element: &PyQuotientElement,
    ) -> PyResult<Vec<usize>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            Ok(super::element_utilities::quotient_support_for(
                algebra,
                element.kind.clone(),
                super::element_utilities::SupportSide::Right,
            )?
            .pop()
            .unwrap_or_default())
        })
    }

    fn support_of_quiver_algebra_element(
        &self,
        element: &PyQuotientElement,
    ) -> PyResult<Vec<Vec<usize>>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            super::element_utilities::quotient_support_for(
                algebra,
                element.kind.clone(),
                super::element_utilities::SupportSide::Both,
            )
        })
    }

    fn make_uniform_on_right(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyQuotientElement>>,
    ) -> PyResult<Vec<PyQuotientElement>> {
        let generators = generators
            .into_iter()
            .map(|element| element.borrow(py).kind.clone())
            .collect();
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            let generators = runtime::extract_all(generators)?;
            Ok(wrap_quotient_elements(
                algebra
                    .make_uniform_on_right(&generators)
                    .map_err(invalid)?,
            ))
        })
    }

    fn generators_times_arrows_on_right(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyQuotientElement>>,
    ) -> PyResult<Vec<PyQuotientElement>> {
        let generators = generators
            .into_iter()
            .map(|element| element.borrow(py).kind.clone())
            .collect();
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            let generators = runtime::extract_all(generators)?;
            Ok(wrap_quotient_elements(
                algebra
                    .generators_times_arrows_on_right(&generators)
                    .map_err(invalid)?,
            ))
        })
    }

    fn tip_monomial_and_coefficient_of_vector(
        &self,
        py: Python<'_>,
        vector: Vec<Py<PyQuotientElement>>,
    ) -> PyResult<(usize, PyPath, PyScalar)> {
        let vector = vector
            .into_iter()
            .map(|element| element.borrow(py).kind.clone())
            .collect();
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            vector_reduction::quotient_tip_for(algebra, vector)
        })
    }

    /// Mutate and return the identical outer Python list, matching QPA.
    fn tip_reduce_vectors(&self, vectors: &Bound<'_, PyList>) -> PyResult<Option<Py<PyList>>> {
        let vector_kinds = vector_reduction::quotient_rows(vectors)?;
        let reduced = dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            vector_reduction::quotient_reduce_for(algebra, vector_kinds)
        })?;
        let Some(reduced) = reduced else {
            return Ok(None);
        };
        vector_reduction::replace_quotient_rows(vectors, reduced)?;
        Ok(Some(vectors.clone().unbind()))
    }

    fn coefficients_of_vectors(
        &self,
        py: Python<'_>,
        basis: Vec<Vec<Py<PyQuotientElement>>>,
        vector: Vec<Py<PyQuotientElement>>,
    ) -> PyResult<PyQuotientVectorCoefficientsResult> {
        let basis = basis
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|element| element.borrow(py).kind.clone())
                    .collect()
            })
            .collect();
        let vector = vector
            .into_iter()
            .map(|element| element.borrow(py).kind.clone())
            .collect();
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            vector_reduction::quotient_coefficients_for(algebra, basis, vector)
        })
    }

    fn ordered_by(&self, ordering: &PyQuiverOrdering) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientKind => QuotientKind,
            |algebra| algebra.ordered_by(&ordering.inner).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[pyo3(signature = (read_only_names=Vec::new()))]
    fn generator_variable_assignments(
        &self,
        read_only_names: Vec<String>,
    ) -> PyResult<Vec<(String, PyQuiverAlgebraElement)>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            mapping::generator_variable_assignments_for(
                &QuiverAlgebra::Quotient(algebra.clone()),
                &read_only_names,
            )
        })
    }

    /// Construct QPA's quotient by the indicated power of the arrow ideal.
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
            FieldKind::Extension { .. } => QuotientKind::Extension(
                qpa_gbnp::truncated_path_algebra(
                    field
                        .kind
                        .galois_field()
                        .ok_or_else(|| invalid("expected an extension field"))?,
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

    #[getter]
    fn loewy_length(&self) -> PyResult<usize> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            qpa_gbnp::loewy_length_of_algebra(&QuiverAlgebra::Quotient(algebra.clone()))
                .map_err(invalid)
        })
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

    /// QPA's `RelatorsOfFpAlgebra` is the same ordered defining list.
    #[getter]
    fn relators(&self) -> Vec<PyAlgebraElement> {
        self.relations()
    }

    /// Return QPA's `IsPathRing` property.
    #[getter]
    fn is_path_ring(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra.is_path_ring())
    }

    /// Return the property installed by QPA's f.p. path-algebra constructor.
    #[getter]
    fn is_full_fp_path_algebra(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra
            .is_full_fp_path_algebra())
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

    /// Return the residue class of the named vertex or arrow generator.
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

    /// Invoke QPA's explicit three-argument quotient-element constructor.
    fn element_with_normality(
        &self,
        element: &PyAlgebraElement,
        known_normal: bool,
    ) -> PyResult<PyQuotientElement> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            QuotientKind => QuotientElementKind,
            |algebra| algebra
                .element_with_normality(runtime::Extract::extract(element)?, known_normal)
                .map_err(invalid)?
        );
        Ok(PyQuotientElement { kind })
    }

    /// Apply the normal-form function installed on this quotient family.
    fn normal_form(&self, element: &PyAlgebraElement) -> PyResult<PyQuotientElement> {
        if !self.has_normal_form_function() {
            return Err(invalid("the quotient family has no normal-form function"));
        }
        self.element_with_normality(element, false)
    }

    /// Rebuild a residue class through QPA's quotient `ObjByExtRep` method.
    fn element_from_external_representation(
        &self,
        py: Python<'_>,
        representation: Vec<(Vec<usize>, Py<PyScalar>)>,
    ) -> PyResult<PyQuotientElement> {
        let kind = map_ref!(
            &self.kind,
            QuotientKind => QuotientElementKind,
            |algebra| algebra
                .element_from_external_representation(&path_external_representation(
                    py,
                    algebra.original_path_algebra(),
                    representation,
                )?)
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

    /// Return the finite dimension as a decimal string, or `None`.
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

    fn canonical_basis_object(&self) -> PyResult<PyQuotientCanonicalBasis> {
        let kind = map_ref!(
            &self.kind,
            QuotientKind => QuotientCanonicalBasisKind,
            |algebra| algebra.canonical_basis_object().map_err(invalid)?
        );
        Ok(PyQuotientCanonicalBasis { kind })
    }

    /// Return QPA's Frobenius bilinear-form closure, or its cached `false`.
    fn frobenius_form(&self) -> PyResult<Option<super::PyFrobeniusForm>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            super::nakayama::form_for(algebra)
        })
    }

    /// Return QPA's Frobenius functional closure, or its cached `false`.
    fn frobenius_linear_functional(&self) -> PyResult<Option<super::PyFrobeniusLinearFunctional>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| {
            super::nakayama::functional_for(algebra)
        })
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

    /// Construct a right module and verify all quotient relations.
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

    /// Return standard simple modules in QPA vertex order.
    fn simple_modules(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| Ok(wrap_matrix_modules(
            qpa_gbnp::simple_modules(&QuiverAlgebra::Quotient(algebra.clone())).map_err(invalid)?,
        )))
    }

    /// Return indecomposable projectives in QPA vertex order.
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

    fn view(&self) -> String {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra.view_obj_with(
            &algebra
                .original_path_algebra()
                .field()
                .python_kind()
                .display()
        ))
    }

    fn print_obj(&self) -> String {
        dispatch_ref!(&self.kind, QuotientKind, |algebra| algebra.print_obj_with(
            &algebra
                .original_path_algebra()
                .field()
                .python_kind()
                .display()
        ))
    }

    fn __repr__(&self) -> String {
        self.view()
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
    fn tip(&self) -> Self {
        self.leading_term()
    }

    #[getter]
    fn tip_coefficient(&self) -> PyScalar {
        self.leading_coefficient()
    }

    #[getter]
    fn tip_monomial(&self) -> PyPath {
        self.leading_monomial()
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

    /// Return QPA's `OneOp` result in this element's quotient family.
    fn one_same_parent(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| element.algebra().one().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    /// Return QPA's `ZeroOp` result in this element's quotient family.
    fn zero_same_parent(&self) -> Self {
        let kind = map_ref!(
            &self.kind,
            QuotientElementKind => QuotientElementKind,
            |element| element.algebra().zero()
        );
        Self { kind }
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

    fn print_obj(&self) -> String {
        self.__repr__()
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
impl PyQuotientCanonicalBasis {
    #[getter]
    fn quotient(&self) -> PyQuotientAlgebra {
        let kind = map_ref!(
            &self.kind,
            QuotientCanonicalBasisKind => QuotientKind,
            |basis| basis.quotient().clone()
        );
        PyQuotientAlgebra { kind }
    }

    #[getter]
    fn has_basis_vectors(&self) -> bool {
        dispatch_ref!(&self.kind, QuotientCanonicalBasisKind, |basis| basis
            .has_basis_vectors())
    }

    #[getter]
    fn basis_vectors(&self) -> Option<Vec<PyQuotientElement>> {
        dispatch_ref!(&self.kind, QuotientCanonicalBasisKind, |basis| basis
            .basis_vectors()
            .map(wrap_quotient_elements))
    }

    fn coefficients(&self, element: &PyQuotientElement) -> PyResult<Vec<PyScalar>> {
        let element = element.kind.clone();
        dispatch_ref!(&self.kind, QuotientCanonicalBasisKind, |basis| {
            let coefficients = basis
                .coefficients(&runtime::Extract::extract(element)?)
                .map_err(invalid)?;
            Ok(coefficients
                .iter()
                .map(|value| {
                    basis
                        .quotient()
                        .original_path_algebra()
                        .field()
                        .wrap_scalar(value)
                })
                .collect())
        })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            QuotientCanonicalBasisKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}
