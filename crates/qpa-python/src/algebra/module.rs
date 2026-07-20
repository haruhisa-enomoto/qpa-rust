//! Matrix modules, module elements, and module homomorphisms.

use pyo3::pymethods;

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

fn wrap_module<F: runtime::PythonField>(module: PathAlgebraModule<F>) -> PyMatrixModule
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    PyMatrixModule {
        kind: runtime::Wrap::wrap(module),
    }
}

fn wrap_module_element<F: runtime::PythonField>(
    element: PathAlgebraModuleElement<F>,
) -> PyModuleElement
where
    ModuleElementKind: runtime::Wrap<PathAlgebraModuleElement<F>>,
{
    PyModuleElement {
        kind: runtime::Wrap::wrap(element),
    }
}

fn wrap_module_elements<F: runtime::PythonField>(
    elements: Vec<PathAlgebraModuleElement<F>>,
) -> Vec<PyModuleElement>
where
    ModuleElementKind: runtime::Wrap<PathAlgebraModuleElement<F>>,
{
    elements.into_iter().map(wrap_module_element).collect()
}

fn wrap_module_homomorphism<F: runtime::PythonField>(
    homomorphism: PathAlgebraModuleHomomorphism<F>,
) -> PyModuleHomomorphism
where
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    PyModuleHomomorphism {
        kind: runtime::Wrap::wrap(homomorphism),
    }
}

fn wrap_module_homomorphisms<F: runtime::PythonField>(
    homomorphisms: Vec<PathAlgebraModuleHomomorphism<F>>,
) -> Vec<PyModuleHomomorphism>
where
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    homomorphisms
        .into_iter()
        .map(wrap_module_homomorphism)
        .collect()
}

fn scalar_rows<F: runtime::PythonField>(field: &F, rows: &[Vec<F::Element>]) -> Vec<Vec<PyScalar>> {
    rows.iter()
        .map(|row| row.iter().map(|value| field.wrap_scalar(value)).collect())
        .collect()
}

fn matrix_rows<F: runtime::PythonField>(matrix: &DenseMatrix<F>) -> Vec<Vec<PyScalar>> {
    scalar_rows(
        matrix.field(),
        &(0..matrix.nrows())
            .map(|row| matrix.row(row).to_vec())
            .collect::<Vec<_>>(),
    )
}

fn extract_scalar_rows<F: runtime::PythonField>(
    py: Python<'_>,
    field: &F,
    rows: Vec<Vec<Py<PyScalar>>>,
) -> PyResult<Vec<Vec<F::Element>>> {
    rows.into_iter()
        .map(|row| {
            row.into_iter()
                .map(|value| field.scalar_value(&value.borrow(py)))
                .collect()
        })
        .collect()
}

fn build_module_element<F: runtime::PythonField>(
    py: Python<'_>,
    module: &PathAlgebraModule<F>,
    coordinates: Vec<Vec<Py<PyScalar>>>,
) -> PyResult<PathAlgebraModuleElement<F>> {
    let field = module.right_acting_algebra().field();
    module
        .element_from_vertex_coordinates(extract_scalar_rows(py, field, coordinates)?)
        .map_err(invalid)
}

fn build_module_homomorphism<F: runtime::PythonField>(
    py: Python<'_>,
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
    matrices: Vec<Vec<Vec<Py<PyScalar>>>>,
) -> PyResult<PathAlgebraModuleHomomorphism<F>> {
    let field = source.right_acting_algebra().field();
    let matrices = matrices
        .into_iter()
        .map(|rows| {
            DenseMatrix::from_rows(extract_scalar_rows(py, field, rows)?, field.clone())
                .map_err(invalid)
        })
        .collect::<PyResult<Vec<_>>>()?;
    qpa_core::right_module_hom_over_algebra(source.clone(), range.clone(), matrices)
        .map_err(invalid)
}

fn direct_sum_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
) -> PyResult<Option<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>> + runtime::Wrap<PathAlgebraModule<F>>,
{
    Ok(
        qpa_core::direct_sum_of_qpa_modules(&runtime::extract_all(modules)?)
            .map_err(invalid)?
            .map(wrap_module),
    )
}

#[pymethods]
impl PyMatrixModule {
    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .right_acting_algebra()
            .field()
            .python_kind());
        PyField { kind }
    }

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
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module.is_zero())
    }

    #[getter]
    fn arrow_matrix_shapes(&self) -> Vec<(usize, usize)> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| matrix_shapes(
            module.arrow_matrices()
        ))
    }

    #[getter]
    fn arrow_matrices(&self) -> Vec<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .arrow_matrices()
            .iter()
            .map(matrix_rows)
            .collect())
    }

    fn zero_element(&self) -> PyModuleElement {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleElementKind,
            |module| module.zero_element()
        );
        PyModuleElement { kind }
    }

    fn element(
        &self,
        py: Python<'_>,
        coordinates: Vec<Vec<Py<PyScalar>>>,
    ) -> PyResult<PyModuleElement> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleElementKind,
            |module| build_module_element(py, module, coordinates)?
        );
        Ok(PyModuleElement { kind })
    }

    fn basis_elements(&self) -> Vec<PyModuleElement> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| wrap_module_elements(
            module.basis_elements()
        ))
    }

    #[staticmethod]
    fn direct_sum(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
    ) -> PyResult<Option<PyMatrixModule>> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        match modules.first() {
            None => Ok(None),
            Some(MatrixModuleKind::Rational(_)) => direct_sum_for::<RationalField>(modules),
            Some(MatrixModuleKind::Binary(_)) => direct_sum_for::<Z2>(modules),
            Some(MatrixModuleKind::Prime(_)) => direct_sum_for::<Zp>(modules),
        }
    }

    #[getter]
    fn is_direct_sum(&self) -> bool {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .is_direct_sum_of_modules())
    }

    fn direct_sum_summands(&self) -> Option<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .direct_sum_summands()
            .map(|summands| summands.iter().cloned().map(wrap_module).collect()))
    }

    fn direct_sum_projections(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(qpa_core::direct_sum_projections(module).map_err(invalid)?,)
        ))
    }

    fn direct_sum_inclusions(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(qpa_core::direct_sum_inclusions(module).map_err(invalid)?,)
        ))
    }

    fn identity(&self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module| qpa_core::identity_module_mapping(module.clone()).map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    fn homomorphism_to(
        &self,
        py: Python<'_>,
        range: &Self,
        vertex_matrices: Vec<Vec<Vec<Py<PyScalar>>>>,
    ) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, range| build_module_homomorphism(py, source, range, vertex_matrices)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn radical(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_core::radical_of_module(module.clone()).map_err(invalid)?,
        )))
    }

    fn top(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_core::top_of_module(module.clone()).map_err(invalid)?,
        )))
    }

    fn socle(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_core::socle_of_module(module.clone()).map_err(invalid)?,
        )))
    }

    #[getter]
    fn is_simple(&self) -> bool {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_core::is_simple_qpa_module(module)
        })
    }

    #[getter]
    fn is_semisimple(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_core::is_semisimple_module(module).map_err(invalid)
        })
    }

    #[getter]
    fn loewy_length(&self) -> PyResult<usize> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_core::loewy_length_of_module(module.clone()).map_err(invalid)
        })
    }

    fn minimal_generating_set(&self) -> PyResult<Vec<PyModuleElement>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_elements(
                qpa_core::minimal_generating_set_of_module(module.clone()).map_err(invalid)?,
            )
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

    fn __eq__(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            MatrixModuleKind,
            |left, right| left == right,
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .qpa_print_string())
    }
}

#[pymethods]
impl PyModuleElement {
    #[getter]
    fn module(&self) -> PyMatrixModule {
        let kind = map_ref!(
            &self.kind,
            ModuleElementKind => MatrixModuleKind,
            |element| element.module().clone()
        );
        PyMatrixModule { kind }
    }

    #[getter]
    fn external_representation(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, ModuleElementKind, |element| scalar_rows(
            element.module().right_acting_algebra().field(),
            element.external_representation(),
        ))
    }

    #[getter]
    fn basis_coefficients(&self) -> Vec<PyScalar> {
        dispatch_ref!(&self.kind, ModuleElementKind, |element| {
            let field = element.module().right_acting_algebra().field();
            element
                .basis_coefficients()
                .iter()
                .map(|value| field.wrap_scalar(value))
                .collect()
        })
    }

    #[getter]
    fn support_vertices(&self) -> Vec<usize> {
        dispatch_ref!(&self.kind, ModuleElementKind, |element| element
            .support_vertices())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, ModuleElementKind, |element| element.is_zero())
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleElementKind => ModuleElementKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("module elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn subtract(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleElementKind => ModuleElementKind,
            |left, right| left.subtract(right).map_err(invalid)?,
            invalid("module elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn negate(&self) -> Self {
        let kind = map_ref!(
            &self.kind,
            ModuleElementKind => ModuleElementKind,
            |element| element.additive_inverse()
        );
        Self { kind }
    }

    fn scale_left(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleElementKind => ModuleElementKind,
            |element| element.scale_left(
                &element
                    .module()
                    .right_acting_algebra()
                    .field()
                    .scalar_value(scalar)?,
            )
        );
        Ok(Self { kind })
    }

    fn scale_right(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleElementKind => ModuleElementKind,
            |element| element.scale_right(
                &element
                    .module()
                    .right_acting_algebra()
                    .field()
                    .scalar_value(scalar)?,
            )
        );
        Ok(Self { kind })
    }

    fn act_by_algebra_element(&self, element: &PyAlgebraElement) -> PyResult<Self> {
        let kind = match (&self.kind, &element.kind) {
            (ModuleElementKind::Rational(left), AlgebraElementKind::Rational(right)) => {
                ModuleElementKind::Rational(
                    left.act(&QuiverAlgebraElement::Path(right.clone()))
                        .map_err(invalid)?,
                )
            }
            (ModuleElementKind::Binary(left), AlgebraElementKind::Binary(right)) => {
                ModuleElementKind::Binary(
                    left.act(&QuiverAlgebraElement::Path(right.clone()))
                        .map_err(invalid)?,
                )
            }
            (ModuleElementKind::Prime(left), AlgebraElementKind::Prime(right)) => {
                ModuleElementKind::Prime(
                    left.act(&QuiverAlgebraElement::Path(right.clone()))
                        .map_err(invalid)?,
                )
            }
            _ => {
                return Err(invalid(
                    "acting element belongs to another coefficient field",
                ));
            }
        };
        Ok(Self { kind })
    }

    fn act_by_quotient_element(&self, element: &PyQuotientElement) -> PyResult<Self> {
        let kind = match (&self.kind, &element.kind) {
            (ModuleElementKind::Rational(left), QuotientElementKind::Rational(right)) => {
                ModuleElementKind::Rational(
                    left.act(&QuiverAlgebraElement::Quotient(right.clone()))
                        .map_err(invalid)?,
                )
            }
            (ModuleElementKind::Binary(left), QuotientElementKind::Binary(right)) => {
                ModuleElementKind::Binary(
                    left.act(&QuiverAlgebraElement::Quotient(right.clone()))
                        .map_err(invalid)?,
                )
            }
            (ModuleElementKind::Prime(left), QuotientElementKind::Prime(right)) => {
                ModuleElementKind::Prime(
                    left.act(&QuiverAlgebraElement::Quotient(right.clone()))
                        .map_err(invalid)?,
                )
            }
            _ => {
                return Err(invalid(
                    "acting element belongs to another coefficient field",
                ));
            }
        };
        Ok(Self { kind })
    }

    fn __add__(&self, right: &Self) -> PyResult<Self> {
        self.add(right)
    }

    fn __sub__(&self, right: &Self) -> PyResult<Self> {
        self.subtract(right)
    }

    fn __neg__(&self) -> Self {
        self.negate()
    }

    fn __eq__(&self, right: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleElementKind,
            |left, right| left == right,
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, ModuleElementKind, |element| {
            let field = element.module().right_acting_algebra().field();
            element.qpa_print_string_with(|value| field.scalar_text(value))
        })
    }
}

#[pymethods]
impl PyModuleHomomorphism {
    #[getter]
    fn source(&self) -> PyMatrixModule {
        let kind = map_ref!(
            &self.kind,
            ModuleHomomorphismKind => MatrixModuleKind,
            |homomorphism| homomorphism.source().clone()
        );
        PyMatrixModule { kind }
    }

    #[getter]
    fn range(&self) -> PyMatrixModule {
        let kind = map_ref!(
            &self.kind,
            ModuleHomomorphismKind => MatrixModuleKind,
            |homomorphism| homomorphism.range().clone()
        );
        PyMatrixModule { kind }
    }

    #[getter]
    fn vertex_matrices(&self) -> Vec<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism
                .vertex_matrices()
                .iter()
                .map(matrix_rows)
                .collect()
        })
    }

    fn image(&self, element: &PyModuleElement) -> PyResult<PyModuleElement> {
        let kind = match (&self.kind, &element.kind) {
            (ModuleHomomorphismKind::Rational(map), ModuleElementKind::Rational(element)) => {
                ModuleElementKind::Rational(map.image(element).map_err(invalid)?)
            }
            (ModuleHomomorphismKind::Binary(map), ModuleElementKind::Binary(element)) => {
                ModuleElementKind::Binary(map.image(element).map_err(invalid)?)
            }
            (ModuleHomomorphismKind::Prime(map), ModuleElementKind::Prime(element)) => {
                ModuleElementKind::Prime(map.image(element).map_err(invalid)?)
            }
            _ => {
                return Err(invalid(
                    "map and element belong to different coefficient fields",
                ));
            }
        };
        Ok(PyModuleElement { kind })
    }

    fn preimage_representative(
        &self,
        element: &PyModuleElement,
    ) -> PyResult<Option<PyModuleElement>> {
        match (&self.kind, &element.kind) {
            (ModuleHomomorphismKind::Rational(map), ModuleElementKind::Rational(element)) => {
                Ok(map
                    .preimage_representative(element)
                    .map_err(invalid)?
                    .map(wrap_module_element))
            }
            (ModuleHomomorphismKind::Binary(map), ModuleElementKind::Binary(element)) => Ok(map
                .preimage_representative(element)
                .map_err(invalid)?
                .map(wrap_module_element)),
            (ModuleHomomorphismKind::Prime(map), ModuleElementKind::Prime(element)) => Ok(map
                .preimage_representative(element)
                .map_err(invalid)?
                .map(wrap_module_element)),
            _ => Err(invalid(
                "map and element belong to different coefficient fields",
            )),
        }
    }

    fn zero(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleHomomorphismKind => ModuleHomomorphismKind,
            |homomorphism| homomorphism.zero().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleHomomorphismKind => ModuleHomomorphismKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("homomorphisms belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn compose(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleHomomorphismKind => ModuleHomomorphismKind,
            |left, right| left.compose(right).map_err(invalid)?,
            invalid("homomorphisms belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn scale(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleHomomorphismKind => ModuleHomomorphismKind,
            |homomorphism| homomorphism
                .scale(
                    &homomorphism
                        .source()
                        .right_acting_algebra()
                        .field()
                        .scalar_value(scalar)?,
                )
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn negate(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleHomomorphismKind => ModuleHomomorphismKind,
            |homomorphism| homomorphism.negate().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism.is_zero()
        })
    }

    #[getter]
    fn is_injective(&self) -> bool {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism.is_injective()
        })
    }

    #[getter]
    fn is_surjective(&self) -> bool {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism.is_surjective()
        })
    }

    #[getter]
    fn is_isomorphism(&self) -> bool {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism.is_isomorphism()
        })
    }

    #[getter]
    fn is_split_monomorphism(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            qpa_core::is_split_monomorphism(homomorphism).map_err(invalid)
        })
    }

    #[getter]
    fn is_split_epimorphism(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            qpa_core::is_split_epimorphism(homomorphism).map_err(invalid)
        })
    }

    fn inverse(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_core::inverse_of_isomorphism(homomorphism).map_err(invalid)?,
            )
        ))
    }

    fn right_inverse(&self) -> PyResult<Option<Self>> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            qpa_core::right_inverse_of_homomorphism(homomorphism)
                .map_err(invalid)?
                .map(wrap_module_homomorphism)
        ))
    }

    fn left_inverse(&self) -> PyResult<Option<Self>> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            qpa_core::left_inverse_of_homomorphism(homomorphism)
                .map_err(invalid)?
                .map(wrap_module_homomorphism)
        ))
    }

    fn kernel(&self) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module(qpa_core::kernel_of_module_homomorphism(homomorphism).map_err(invalid)?,)
        ))
    }

    fn image_module(&self) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module(qpa_core::image_of_module_homomorphism(homomorphism).map_err(invalid)?,)
        ))
    }

    fn cokernel(&self) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module(qpa_core::cokernel_of_module_homomorphism(homomorphism).map_err(invalid)?,)
        ))
    }

    fn __add__(&self, right: &Self) -> PyResult<Self> {
        self.add(right)
    }

    fn __neg__(&self) -> PyResult<Self> {
        self.negate()
    }

    fn __eq__(&self, right: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleHomomorphismKind,
            |left, right| left == right,
            false
        )
    }

    fn is_identical(&self, right: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleHomomorphismKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism.qpa_print_string()
        })
    }
}

fn matrix_shapes<F: numerica::domains::Field>(matrices: &[DenseMatrix<F>]) -> Vec<(usize, usize)> {
    matrices
        .iter()
        .map(|matrix| (matrix.nrows(), matrix.ncols()))
        .collect()
}
