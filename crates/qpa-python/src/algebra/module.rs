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

fn wrap_module_homomorphism_pair<F: runtime::PythonField>(
    pair: Option<(
        PathAlgebraModuleHomomorphism<F>,
        PathAlgebraModuleHomomorphism<F>,
    )>,
) -> Option<(PyModuleHomomorphism, PyModuleHomomorphism)>
where
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    pair.map(|(first, second)| {
        (
            wrap_module_homomorphism(first),
            wrap_module_homomorphism(second),
        )
    })
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

fn right_approximation_by_add_list_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
    target: &PathAlgebraModule<F>,
) -> PyResult<PyModuleHomomorphism>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>>,
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    Ok(wrap_module_homomorphism(
        qpa_gbnp::right_approximation_by_add_m_list(&runtime::extract_all(modules)?, target)
            .map_err(invalid)?,
    ))
}

fn make_basic_list_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
) -> PyResult<Vec<PyMatrixModule>>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>> + runtime::Wrap<PathAlgebraModule<F>>,
{
    Ok(
        qpa_core::make_basic_list_of_modules(&runtime::extract_all(modules)?)
            .map_err(invalid)?
            .into_iter()
            .map(wrap_module)
            .collect(),
    )
}

fn radical_right_approximation_by_add_for<F: runtime::PythonField>(
    modules: Vec<MatrixModuleKind>,
    index: &Integer,
) -> PyResult<PyModuleHomomorphism>
where
    MatrixModuleKind: runtime::Extract<PathAlgebraModule<F>>,
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    Ok(wrap_module_homomorphism(
        qpa_gbnp::radical_right_approximation_by_add_m(&runtime::extract_all(modules)?, index)
            .map_err(invalid)?,
    ))
}

fn faithful_dimension_result_text(result: qpa_gbnp::FaithfulDimensionResult) -> String {
    match result {
        qpa_gbnp::FaithfulDimensionResult::Finite(value) => value.to_string(),
        qpa_gbnp::FaithfulDimensionResult::Infinite => "infinity".into(),
    }
}

fn number_of_complements_result_text(result: qpa_gbnp::NumberOfComplementsResult) -> String {
    match result {
        qpa_gbnp::NumberOfComplementsResult::Finite(value) => value.to_string(),
        qpa_gbnp::NumberOfComplementsResult::Infinite => "infinity".into(),
    }
}

fn clone_module_elements(
    py: Python<'_>,
    elements: Vec<Py<PyModuleElement>>,
) -> Vec<ModuleElementKind> {
    elements
        .into_iter()
        .map(|element| element.borrow(py).kind.clone())
        .collect()
}

fn wrap_module_submodule<F: runtime::PythonField>(
    submodule: PathAlgebraModuleSubmodule<F>,
) -> PyModuleSubmodule
where
    ModuleSubmoduleKind: runtime::Wrap<PathAlgebraModuleSubmodule<F>>,
{
    PyModuleSubmodule {
        kind: runtime::Wrap::wrap(submodule),
    }
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

    fn decompose(&self) -> PyResult<Vec<Self>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_core::decompose_module(module,)
                .map_err(invalid)?
                .into_iter()
                .map(wrap_module)
                .collect()
        ))
    }

    #[getter]
    fn is_indecomposable(&self) -> PyResult<Option<bool>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_core::is_indecomposable_module(module).map_err(invalid)
        })
    }

    fn decompose_with_inclusions(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(
                qpa_core::decompose_module_with_inclusions(module).map_err(invalid)?,
            )
        ))
    }

    fn decompose_with_multiplicities(
        &self,
    ) -> PyResult<Option<PyModuleDecompositionWithMultiplicities>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_core::decompose_module_with_multiplicities(module)
                .map_err(invalid)?
                .map(|result| PyModuleDecompositionWithMultiplicities {
                    kind: runtime::Wrap::wrap(result),
                })
        ))
    }

    fn block_splitting_idempotents(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(
                qpa_core::block_splitting_idempotents(module).map_err(invalid)?,
            )
        ))
    }

    fn block_decomposition(&self) -> PyResult<Vec<Self>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_core::block_decomposition_of_module(module)
                .map_err(invalid)?
                .into_iter()
                .map(wrap_module)
                .collect()
        ))
    }

    fn basic_version(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_core::basic_version_of_module(module).map_err(invalid)?,
        )))
    }

    #[staticmethod]
    fn make_basic(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
    ) -> PyResult<Vec<PyMatrixModule>> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        match modules.first() {
            Some(MatrixModuleKind::Rational(_)) | None => {
                make_basic_list_for::<RationalField>(modules)
            }
            Some(MatrixModuleKind::Binary(_)) => make_basic_list_for::<Z2>(modules),
            Some(MatrixModuleKind::Prime(_)) => make_basic_list_for::<Zp>(modules),
        }
    }

    fn nonisomorphic_direct_summand_data(&self) -> PyResult<PyNonisomorphicDirectSummandData> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            PyNonisomorphicDirectSummandData {
                value: qpa_core::number_of_nonisomorphic_direct_summands(module)
                    .map_err(invalid)?,
            }
        ))
    }

    fn isomorphic_to(&self, other: &Self) -> PyResult<bool> {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            MatrixModuleKind,
            |source, range| qpa_core::isomorphic_modules(source, range).map_err(invalid),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn isomorphism_to(&self, other: &Self) -> PyResult<Option<PyModuleHomomorphism>> {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            MatrixModuleKind,
            |source, range| Ok(qpa_core::isomorphism_of_modules(source, range)
                .map_err(invalid)?
                .map(wrap_module_homomorphism)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn is_direct_summand_of(&self, module: &Self) -> PyResult<bool> {
        dispatch_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind,
            |candidate, module| qpa_core::is_direct_summand(candidate, module).map_err(invalid),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn is_in_additive_closure_of(&self, generator: &Self) -> PyResult<Option<bool>> {
        dispatch_pair_ref!(
            &self.kind,
            &generator.kind,
            MatrixModuleKind,
            |candidate, generator| qpa_core::is_in_additive_closure(candidate, generator)
                .map_err(invalid),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn identity(&self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module| qpa_core::identity_module_mapping(module.clone()).map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    fn zero_map_to(&self, range: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, range| qpa_core::zero_module_mapping(source.clone(), range.clone())
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn hom_basis_to(&self, range: &Self) -> PyResult<Option<Vec<PyModuleHomomorphism>>> {
        dispatch_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind,
            |source, range| Ok(qpa_core::hom_over_algebra(source, range)
                .map_err(invalid)?
                .map(wrap_module_homomorphisms)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn hom_factoring_through_projectives_to(
        &self,
        range: &Self,
    ) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind,
            |source, range| Ok(wrap_module_homomorphisms(
                qpa_gbnp::hom_factoring_through_projectives(source, range).map_err(invalid)?,
            )),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn homomorphism_from_images(
        &self,
        py: Python<'_>,
        range: &Self,
        images: Vec<Py<PyModuleElement>>,
    ) -> PyResult<PyModuleHomomorphism> {
        let images = clone_module_elements(py, images);
        dispatch_pair_ref!(
            &self.kind,
            &range.kind,
            MatrixModuleKind,
            |source, range| Ok(wrap_module_homomorphism(
                qpa_core::homomorphism_from_images(source, range, &runtime::extract_all(images)?,)
                    .map_err(invalid)?,
            )),
            Err(invalid("modules belong to different coefficient fields"))
        )
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

    fn submodule(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyModuleElement>>,
    ) -> PyResult<PyModuleSubmodule> {
        let generators = clone_module_elements(py, generators);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_submodule(
                qpa_core::subalgebra_module(module.clone(), &runtime::extract_all(generators)?,)
                    .map_err(invalid)?,
            )
        ))
    }

    fn subrepresentation(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyModuleElement>>,
    ) -> PyResult<Self> {
        let generators = clone_module_elements(py, generators);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_core::subrepresentation(module.clone(), &runtime::extract_all(generators)?,)
                .map_err(invalid)?,
        )))
    }

    fn subrepresentation_inclusion(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyModuleElement>>,
    ) -> PyResult<PyModuleHomomorphism> {
        let generators = clone_module_elements(py, generators);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(
                qpa_core::subrepresentation_inclusion(
                    module.clone(),
                    &runtime::extract_all(generators)?,
                )
                .map_err(invalid)?,
            )
        ))
    }

    fn radical_inclusion(&self) -> PyResult<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(
                qpa_core::radical_of_module_inclusion(module.clone()).map_err(invalid)?,
            )
        ))
    }

    fn top_projection(&self) -> PyResult<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(
                qpa_core::top_of_module_projection(module.clone()).map_err(invalid)?,
            )
        ))
    }

    fn socle_inclusion(&self) -> PyResult<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(
                qpa_core::socle_of_module_inclusion(module.clone()).map_err(invalid)?,
            )
        ))
    }

    #[getter]
    fn is_projective(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_projective_module(module).map_err(invalid)
        })
    }

    #[getter]
    fn is_injective(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_injective_module(module).map_err(invalid)
        })
    }

    fn compare_with_indecomposable_projective(&self) -> PyResult<Option<usize>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::compare_with_indec_projective(module).map_err(invalid)
        })
    }

    fn compare_with_indecomposable_injective(&self) -> PyResult<Option<usize>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::compare_with_indec_injective(module).map_err(invalid)
        })
    }

    fn projective_description(&self) -> PyResult<Vec<Option<usize>>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::description_of_projective_module(module).map_err(invalid)
        })
    }

    fn direct_sum_minus_summands(
        &self,
        removed: Vec<usize>,
    ) -> PyResult<(PyModuleHomomorphism, PyModuleHomomorphism)> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            let (inclusion, projection) =
                qpa_gbnp::direct_sum_minus_summands(module, &removed).map_err(invalid)?;
            Ok((
                wrap_module_homomorphism(inclusion),
                wrap_module_homomorphism(projection),
            ))
        })
    }

    fn projective_cover(&self) -> PyResult<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(qpa_gbnp::projective_cover(module).map_err(invalid)?,)
        ))
    }

    fn injective_envelope(&self) -> PyResult<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism(qpa_gbnp::injective_envelope(module).map_err(invalid)?,)
        ))
    }

    fn first_syzygy(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_gbnp::first_syzygy(module).map_err(invalid)?,
        )))
    }

    fn nth_syzygy(&self, degree: i64) -> PyResult<Self> {
        let degree = Integer::from(degree);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(wrap_module(
            qpa_gbnp::nth_syzygy(module, &degree).map_err(invalid)?,
        )))
    }

    fn projective_resolution(&self) -> PyResult<PyProjectiveResolution> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ProjectiveResolutionKind,
            |module| qpa_gbnp::projective_resolution(module).map_err(invalid)?
        );
        Ok(PyProjectiveResolution { kind })
    }

    fn injective_resolution(&self) -> PyResult<PyInjectiveResolution> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => InjectiveResolutionKind,
            |module| qpa_gbnp::injective_resolution(module).map_err(invalid)?
        );
        Ok(PyInjectiveResolution { kind })
    }

    fn zero_complex(&self) -> PyResult<PyModuleComplex> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleComplexKind,
            |module| qpa_gbnp::zero_complex(
                qpa_gbnp::cat_of_right_algebra_modules(module.right_acting_algebra())
                    .map_err(invalid)?,
            )
            .map_err(invalid)?
        );
        Ok(PyModuleComplex { kind })
    }

    fn stalk_complex(&self, degree: i64) -> PyResult<PyModuleComplex> {
        let degree = Integer::from(degree);
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleComplexKind,
            |module| qpa_gbnp::stalk_complex(
                qpa_gbnp::cat_of_right_algebra_modules(module.right_acting_algebra())
                    .map_err(invalid)?,
                module,
                degree,
            )
            .map_err(invalid)?
        );
        Ok(PyModuleComplex { kind })
    }

    fn finite_complex(
        &self,
        py: Python<'_>,
        base_degree: i64,
        differentials: Vec<Py<PyModuleHomomorphism>>,
    ) -> PyResult<PyModuleComplex> {
        let differentials = clone_module_homomorphisms(py, differentials);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            let category = qpa_gbnp::cat_of_right_algebra_modules(module.right_acting_algebra())
                .map_err(invalid)?;
            let kind = runtime::Wrap::wrap(
                qpa_gbnp::finite_complex(
                    category,
                    Integer::from(base_degree),
                    runtime::extract_all(differentials)?,
                )
                .map_err(invalid)?,
            );
            Ok(PyModuleComplex { kind })
        })
    }

    fn projective_resolution_complex(&self) -> PyResult<PyModuleComplex> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleComplexKind,
            |module| qpa_gbnp::projective_resolution_complex(module).map_err(invalid)?
        );
        Ok(PyModuleComplex { kind })
    }

    fn injective_resolution_complex(&self) -> PyResult<PyModuleComplex> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleComplexKind,
            |module| qpa_gbnp::injective_resolution_complex(module).map_err(invalid)?
        );
        Ok(PyModuleComplex { kind })
    }

    #[getter]
    fn cached_projective_dimension(&self) -> Option<usize> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::proj_dimension(module)
        })
    }

    fn projective_dimension(&self, bound: i64) -> PyResult<Option<usize>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::proj_dimension_of_module(module, &bound).map_err(invalid)
        })
    }

    #[getter]
    fn cached_injective_dimension(&self) -> Option<usize> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::inj_dimension(module)
        })
    }

    fn injective_dimension(&self, bound: i64) -> PyResult<Option<usize>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::inj_dimension_of_module(module, &bound).map_err(invalid)
        })
    }

    fn dominant_dimension(&self, bound: i64) -> PyResult<Option<String>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            dominant_dimension_text(
                qpa_gbnp::dominant_dimension_of_module(module, &bound).map_err(invalid)?,
            )
        ))
    }

    fn is_nth_syzygy(&self, degree: i64) -> PyResult<bool> {
        let degree = Integer::from(degree);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_nth_syzygy(module, &degree).map_err(invalid)
        })
    }

    fn omega_period(&self, bound: i64) -> PyResult<Option<String>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_gbnp::is_omega_periodic(module, &bound)
                .map_err(invalid)?
                .map(|degree| degree.to_string())
        ))
    }

    fn tau_period(&self, bound: i64) -> PyResult<Option<String>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_gbnp::is_tau_periodic(module, &bound)
                .map_err(invalid)?
                .map(|degree| degree.to_string())
        ))
    }

    fn ext1(&self, coefficients_in: &Self) -> PyResult<PyExtOverAlgebraResult> {
        let kind = map_pair_ref!(
            &self.kind,
            &coefficients_in.kind,
            MatrixModuleKind => ExtOverAlgebraKind,
            |module, coefficients_in| qpa_gbnp::ext_over_algebra(module, coefficients_in)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyExtOverAlgebraResult { kind })
    }

    fn ext_algebra_generators(
        &self,
        maximum_degree: usize,
    ) -> PyResult<PyExtAlgebraGeneratorsResult> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ExtAlgebraGeneratorsKind,
            |module| qpa_gbnp::ext_algebra_generators(module, maximum_degree).map_err(invalid)?
        );
        Ok(PyExtAlgebraGeneratorsResult { kind })
    }

    #[getter]
    fn is_rigid(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_rigid_module(module).map_err(invalid)
        })
    }

    #[getter]
    fn is_tau_rigid(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_tau_rigid_module(module).map_err(invalid)
        })
    }

    #[getter]
    fn is_exceptional(&self) -> PyResult<Option<bool>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_exceptional_module(module).map_err(invalid)
        })
    }

    fn is_n_rigid(&self, bound: i64) -> PyResult<bool> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::n_rigid_module(module, &bound).map_err(invalid)
        })
    }

    #[getter]
    fn cached_is_tilting(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_tilting_module(module)
        })
    }

    fn finite_coresolution_in_add(
        &self,
        module: &Self,
        bound: i64,
    ) -> PyResult<Option<PyAddMCoresolution>> {
        let bound = Integer::from(bound);
        dispatch_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind,
            |source, module| Ok(qpa_gbnp::have_finite_coresolution_in_add_m(
                source, module, &bound,
            )
            .map_err(invalid)?
            .map(|result| PyAddMCoresolution {
                kind: runtime::Wrap::wrap(result),
            })),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn tilting(&self, bound: i64) -> PyResult<Option<PyTiltingModuleResult>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_gbnp::tilting_module(module, &bound)
                .map_err(invalid)?
                .map(|result| PyTiltingModuleResult {
                    kind: runtime::Wrap::wrap(result),
                })
        ))
    }

    #[getter]
    fn cached_is_cotilting(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| {
            qpa_gbnp::is_cotilting_module(module)
        })
    }

    fn finite_resolution_in_add(
        &self,
        module: &Self,
        bound: i64,
    ) -> PyResult<Option<PyFiniteAddMResolutionResult>> {
        let bound = Integer::from(bound);
        dispatch_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind,
            |target, module| Ok(
                qpa_gbnp::have_finite_resolution_in_add_m(target, module, &bound,)
                    .map_err(invalid)?
                    .map(|result| PyFiniteAddMResolutionResult {
                        kind: runtime::Wrap::wrap(result),
                    })
            ),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn cotilting(&self, bound: i64) -> PyResult<Option<PyCotiltingModuleResult>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            qpa_gbnp::cotilting_module(module, &bound)
                .map_err(invalid)?
                .map(|result| PyCotiltingModuleResult {
                    kind: runtime::Wrap::wrap(result),
                })
        ))
    }

    fn all_complements(&self, complement: &Self) -> PyResult<PyAllComplementsResult> {
        let kind = map_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind => AllComplementsResultKind,
            |module, complement| qpa_gbnp::all_complements_of_almost_complete_tilting_module(
                module, complement,
            )
            .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyAllComplementsResult { kind })
    }

    fn all_cotilting_complements(&self, complement: &Self) -> PyResult<PyAllComplementsResult> {
        let kind = map_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind => AllComplementsResultKind,
            |module, complement| qpa_gbnp::all_complements_of_almost_complete_cotilting_module(
                module, complement,
            )
            .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyAllComplementsResult { kind })
    }

    #[getter]
    fn cached_faithful_dimension(&self) -> Option<String> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| module
            .cached_faithful_dimension()
            .map(
                |dimension| dimension.map_or_else(|| "infinity".into(), |value| value.to_string())
            ))
    }

    fn faithful_dimension(&self) -> PyResult<String> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            faithful_dimension_result_text(qpa_gbnp::faithful_dimension(module).map_err(invalid)?)
        ))
    }

    fn number_of_tilting_complements(&self) -> PyResult<String> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            number_of_complements_result_text(
                qpa_gbnp::number_of_complements_of_almost_complete_tilting_module(module)
                    .map_err(invalid)?,
            )
        ))
    }

    fn number_of_cotilting_complements(&self) -> PyResult<String> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            number_of_complements_result_text(
                qpa_gbnp::number_of_complements_of_almost_complete_cotilting_module(module)
                    .map_err(invalid)?,
            )
        ))
    }

    fn left_tilting_complement_mutation(&self, complement: &Self) -> PyResult<Option<Self>> {
        dispatch_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind,
            |module, complement| Ok(qpa_gbnp::left_mutation_of_tilting_module_complement(
                module, complement
            )
            .map_err(invalid)?
            .map(wrap_module)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn right_tilting_complement_mutation(&self, complement: &Self) -> PyResult<Option<Self>> {
        dispatch_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind,
            |module, complement| Ok(qpa_gbnp::right_mutation_of_tilting_module_complement(
                module, complement
            )
            .map_err(invalid)?
            .map(wrap_module)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn left_cotilting_complement_mutation(&self, complement: &Self) -> PyResult<Option<Self>> {
        dispatch_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind,
            |module, complement| Ok(qpa_gbnp::left_mutation_of_cotilting_module_complement(
                module, complement
            )
            .map_err(invalid)?
            .map(wrap_module)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn right_cotilting_complement_mutation(&self, complement: &Self) -> PyResult<Option<Self>> {
        dispatch_pair_ref!(
            &self.kind,
            &complement.kind,
            MatrixModuleKind,
            |module, complement| Ok(qpa_gbnp::right_mutation_of_cotilting_module_complement(
                module, complement
            )
            .map_err(invalid)?
            .map(wrap_module)),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn right_approximation_by_add(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &target.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module, target| qpa_gbnp::right_approximation_by_add_m(module, target)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    #[staticmethod]
    fn right_approximation_by_add_list(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
        target: &Self,
    ) -> PyResult<PyModuleHomomorphism> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        match &target.kind {
            MatrixModuleKind::Rational(target) => {
                right_approximation_by_add_list_for::<RationalField>(modules, target)
            }
            MatrixModuleKind::Binary(target) => {
                right_approximation_by_add_list_for::<Z2>(modules, target)
            }
            MatrixModuleKind::Prime(target) => {
                right_approximation_by_add_list_for::<Zp>(modules, target)
            }
        }
    }

    fn minimal_right_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &target.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module, target| qpa_gbnp::minimal_right_approximation(module, target)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_right_add_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        self.minimal_right_approximation(target)
    }

    fn left_approximation_by_add(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, module| qpa_gbnp::left_approximation_by_add_m(source, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_left_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, module| qpa_gbnp::minimal_left_approximation(source, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_left_add_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        self.minimal_left_approximation(module)
    }

    fn right_fac_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &target.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module, target| qpa_core::right_fac_approximation(module, target).map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn right_fac_m_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        self.right_fac_approximation(target)
    }

    fn minimal_right_fac_m_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        self.right_fac_approximation(target)
    }

    fn left_fac_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, module| qpa_gbnp::left_fac_m_approximation(source, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_left_fac_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, module| qpa_gbnp::minimal_left_fac_m_approximation(source, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn right_sub_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &target.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module, target| qpa_gbnp::right_sub_m_approximation(module, target)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_right_sub_approximation(&self, target: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &target.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module, target| qpa_gbnp::minimal_right_sub_m_approximation(module, target)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn left_sub_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |source, module| qpa_gbnp::left_sub_m_approximation(source, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn minimal_left_sub_approximation(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        self.left_sub_approximation(module)
    }

    fn right_approximation_by_perp(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |cotilting, module| qpa_gbnp::right_approximation_by_perp_t(cotilting, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    fn left_approximation_by_add_t_hat(&self, module: &Self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |cotilting, module| qpa_gbnp::left_approximation_by_add_t_hat(cotilting, module)
                .map_err(invalid)?,
            invalid("modules belong to different coefficient fields")
        )?;
        Ok(PyModuleHomomorphism { kind })
    }

    #[staticmethod]
    fn radical_right_approximation_by_add(
        py: Python<'_>,
        modules: Vec<Py<PyMatrixModule>>,
        one_based_index: i64,
    ) -> PyResult<PyModuleHomomorphism> {
        let modules = modules
            .into_iter()
            .map(|module| module.borrow(py).kind.clone())
            .collect::<Vec<_>>();
        let index = Integer::from(one_based_index);
        match modules.first() {
            Some(MatrixModuleKind::Rational(_)) => {
                radical_right_approximation_by_add_for::<RationalField>(modules, &index)
            }
            Some(MatrixModuleKind::Binary(_)) => {
                radical_right_approximation_by_add_for::<Z2>(modules, &index)
            }
            Some(MatrixModuleKind::Prime(_)) => {
                radical_right_approximation_by_add_for::<Zp>(modules, &index)
            }
            // The Rust translation rejects the one-based index before it needs
            // a coefficient field, so any monomorphization preserves QPA's
            // empty-list error instead of inventing a binding-only shortcut.
            None => radical_right_approximation_by_add_for::<RationalField>(modules, &index),
        }
    }

    #[pyo3(signature = (side = "ending"))]
    fn almost_split_sequence(
        &self,
        side: &str,
    ) -> PyResult<Option<(PyModuleHomomorphism, PyModuleHomomorphism)>> {
        let side = match side {
            "ending" | "r" => qpa_gbnp::AlmostSplitSide::Ending,
            "starting" | "l" => qpa_gbnp::AlmostSplitSide::Starting,
            _ => {
                return Err(invalid(
                    "almost-split side must be 'ending'/'r' or 'starting'/'l'",
                ));
            }
        };
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphism_pair(
                qpa_gbnp::almost_split_sequence_from_side(module, side).map_err(invalid)?,
            )
        ))
    }

    fn almost_split_sequence_in_perp(
        &self,
        module: &Self,
    ) -> PyResult<Option<(PyModuleHomomorphism, PyModuleHomomorphism)>> {
        dispatch_pair_ref!(
            &self.kind,
            &module.kind,
            MatrixModuleKind,
            |cotilting, module| Ok(wrap_module_homomorphism_pair(
                qpa_gbnp::almost_split_sequence_in_perp_t(cotilting, module).map_err(invalid)?,
            )),
            Err(invalid("modules belong to different coefficient fields"))
        )
    }

    fn irreducible_morphisms_ending_in(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(
                qpa_gbnp::irreducible_morphisms_ending_in(module).map_err(invalid)?,
            )
        ))
    }

    fn irreducible_morphisms_starting_in(&self) -> PyResult<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, MatrixModuleKind, |module| Ok(
            wrap_module_homomorphisms(
                qpa_gbnp::irreducible_morphisms_starting_in(module).map_err(invalid)?,
            )
        ))
    }

    fn predecessors(&self, bound: i64) -> PyResult<PyPredecessorsOfModuleResult> {
        let bound = Integer::from(bound);
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => PredecessorsOfModuleResultKind,
            |module| qpa_gbnp::predecessors_of_module(module, &bound).map_err(invalid)?
        );
        Ok(PyPredecessorsOfModuleResult { kind })
    }

    fn dual(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::dual_of_module(module).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn transpose(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::transpose_of_module(module).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[pyo3(signature = (exponent=1))]
    fn d_tr(&self, exponent: i64) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::d_tr_iterated(module, exponent).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[pyo3(signature = (exponent=1))]
    fn tr_d(&self, exponent: i64) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::tr_d_iterated(module, exponent).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn star(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::star_of_module(module).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn identity_to_double_star(&self) -> PyResult<PyModuleHomomorphism> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => ModuleHomomorphismKind,
            |module| qpa_gbnp::from_identity_to_double_star_homomorphism(module)
                .map_err(invalid)?
        );
        Ok(PyModuleHomomorphism { kind })
    }

    fn nakayama(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::nakayama_functor_of_module(module).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn opposite_nakayama(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            MatrixModuleKind => MatrixModuleKind,
            |module| qpa_gbnp::opposite_nakayama_functor_of_module(module).map_err(invalid)?
        );
        Ok(Self { kind })
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
impl PyPredecessorValuation {
    #[getter]
    fn predecessor_position(&self) -> usize {
        self.value.predecessor_position
    }

    #[getter]
    fn successor_position(&self) -> usize {
        self.value.successor_position
    }

    #[getter]
    fn scalar_multiplicity(&self) -> Option<usize> {
        match self.value.multiplicity {
            qpa_gbnp::PredecessorMultiplicity::Scalar(value) => Some(value),
            qpa_gbnp::PredecessorMultiplicity::Pair { .. } => None,
        }
    }

    #[getter]
    fn pair_multiplicity(&self) -> Option<(Option<usize>, Option<usize>)> {
        match self.value.multiplicity {
            qpa_gbnp::PredecessorMultiplicity::Pair { first, second } => Some((first, second)),
            qpa_gbnp::PredecessorMultiplicity::Scalar(_) => None,
        }
    }
}

#[pymethods]
impl PyPredecessorsOfModuleResult {
    #[getter]
    fn is_module(&self) -> bool {
        dispatch_ref!(
            &self.kind,
            PredecessorsOfModuleResultKind,
            |result| matches!(result, qpa_gbnp::PredecessorsOfModuleResult::Module(_))
        )
    }

    #[getter]
    fn module(&self) -> Option<PyMatrixModule> {
        dispatch_ref!(
            &self.kind,
            PredecessorsOfModuleResultKind,
            |result| match result {
                qpa_gbnp::PredecessorsOfModuleResult::Module(module) => {
                    Some(wrap_module(module.clone()))
                }
                qpa_gbnp::PredecessorsOfModuleResult::Layers(_) => None,
            }
        )
    }

    #[getter]
    fn layers(&self) -> Option<Vec<Vec<PyMatrixModule>>> {
        dispatch_ref!(
            &self.kind,
            PredecessorsOfModuleResultKind,
            |result| match result {
                qpa_gbnp::PredecessorsOfModuleResult::Module(_) => None,
                qpa_gbnp::PredecessorsOfModuleResult::Layers(layers) => Some(
                    layers
                        .layers
                        .iter()
                        .map(|layer| layer.iter().cloned().map(wrap_module).collect())
                        .collect(),
                ),
            }
        )
    }

    #[getter]
    fn valuations(&self) -> Option<Vec<Vec<PyPredecessorValuation>>> {
        dispatch_ref!(
            &self.kind,
            PredecessorsOfModuleResultKind,
            |result| match result {
                qpa_gbnp::PredecessorsOfModuleResult::Module(_) => None,
                qpa_gbnp::PredecessorsOfModuleResult::Layers(layers) => Some(
                    layers
                        .valuations
                        .iter()
                        .map(|layer| {
                            layer
                                .iter()
                                .cloned()
                                .map(|value| PyPredecessorValuation { value })
                                .collect()
                        })
                        .collect(),
                ),
            }
        )
    }
}

#[pymethods]
impl PyModuleDecompositionWithMultiplicities {
    #[getter]
    fn basic_summands(&self) -> Vec<PyMatrixModule> {
        dispatch_ref!(
            &self.kind,
            ModuleDecompositionWithMultiplicitiesKind,
            |result| result
                .basic_summands()
                .iter()
                .cloned()
                .map(wrap_module)
                .collect()
        )
    }

    #[getter]
    fn multiplicities(&self) -> Vec<usize> {
        dispatch_ref!(
            &self.kind,
            ModuleDecompositionWithMultiplicitiesKind,
            |result| result.multiplicities().to_vec()
        )
    }
}

#[pymethods]
impl PyNonisomorphicDirectSummandData {
    #[getter]
    fn count(&self) -> usize {
        self.value.count()
    }

    #[getter]
    fn block_dimensions(&self) -> Vec<usize> {
        self.value.block_dimensions().to_vec()
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

    #[getter]
    fn dual_vertex_matrices(&self) -> Vec<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            qpa_gbnp::matrices_of_dual_map(homomorphism)
                .iter()
                .map(matrix_rows)
                .collect()
        })
    }

    fn all_components(&self) -> PyResult<Vec<Vec<Self>>> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            qpa_gbnp::find_all_map_components(homomorphism)
                .map_err(invalid)?
                .into_iter()
                .map(wrap_module_homomorphisms)
                .collect()
        ))
    }

    fn dual(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_gbnp::dual_of_module_homomorphism(homomorphism).map_err(invalid)?,
            )
        ))
    }

    fn transpose(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_gbnp::transpose_of_module_homomorphism(homomorphism).map_err(invalid)?,
            )
        ))
    }

    #[pyo3(signature = (exponent=1))]
    fn d_tr(&self, exponent: i64) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_gbnp::d_tr_homomorphism_iterated(homomorphism, exponent).map_err(invalid)?,
            )
        ))
    }

    #[pyo3(signature = (exponent=1))]
    fn tr_d(&self, exponent: i64) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_gbnp::tr_d_homomorphism_iterated(homomorphism, exponent).map_err(invalid)?,
            )
        ))
    }

    fn star(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_gbnp::star_of_module_homomorphism(homomorphism).map_err(invalid)?,
            )
        ))
    }

    fn nakayama(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_gbnp::nakayama_functor_of_module_homomorphism(homomorphism).map_err(invalid)?,
            )
        ))
    }

    fn opposite_nakayama(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_gbnp::opposite_nakayama_functor_of_module_homomorphism(homomorphism)
                    .map_err(invalid)?,
            )
        ))
    }

    fn star_between_indecomposable_projectives(
        &self,
        source_vertex: usize,
        range_description: Vec<usize>,
    ) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_gbnp::star_of_map_between_indec_projectives(
                    homomorphism,
                    source_vertex,
                    &range_description,
                )
                .map_err(invalid)?,
            )
        ))
    }

    fn star_between_decomposable_projectives(
        &self,
        source_description: Vec<usize>,
        range_description: Vec<usize>,
    ) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_gbnp::star_of_map_between_decomp_projectives(
                    homomorphism,
                    &source_description,
                    &range_description,
                )
                .map_err(invalid)?,
            )
        ))
    }

    fn star_between_projectives(
        &self,
        source_description: Vec<usize>,
        range_description: Vec<usize>,
    ) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_gbnp::star_of_map_between_projectives(
                    homomorphism,
                    &source_description,
                    &range_description,
                )
                .map_err(invalid)?,
            )
        ))
    }

    fn check_for_projective_isomorphism(
        &self,
    ) -> PyResult<Option<PyProjectiveIsomorphismComponent>> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            qpa_gbnp::check_for_isomorphism(homomorphism)
                .map_err(invalid)?
                .map(|component| PyProjectiveIsomorphismComponent {
                    kind: runtime::Wrap::wrap(component),
                })
        ))
    }

    fn reduce_differential(
        &self,
        cancellation: &PyProjectiveIsomorphismComponent,
        reduced_source: &PyMatrixModule,
        reduced_range: &PyMatrixModule,
    ) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |differential| {
            let cancellation = runtime::Extract::extract(cancellation.kind.clone())?;
            let reduced_source = runtime::Extract::extract(reduced_source.kind.clone())?;
            let reduced_range = runtime::Extract::extract(reduced_range.kind.clone())?;
            Ok(wrap_module_homomorphism(
                qpa_gbnp::reduce_differential(
                    differential,
                    &cancellation,
                    &reduced_source,
                    &reduced_range,
                )
                .map_err(invalid)?,
            ))
        })
    }

    fn reduce_previous_differential(&self, inclusion: &Self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |differential| {
            let inclusion = runtime::Extract::extract(inclusion.kind.clone())?;
            Ok(wrap_module_homomorphism(
                qpa_gbnp::reduce_previous_differential(differential, &inclusion)
                    .map_err(invalid)?,
            ))
        })
    }

    fn reduce_t_map(
        &self,
        differential: &Self,
        cancellation: &PyProjectiveIsomorphismComponent,
        inclusion: &Self,
    ) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |map| {
            let differential = runtime::Extract::extract(differential.kind.clone())?;
            let cancellation = runtime::Extract::extract(cancellation.kind.clone())?;
            let inclusion = runtime::Extract::extract(inclusion.kind.clone())?;
            Ok(wrap_module_homomorphism(
                qpa_gbnp::reduce_t_map(map, &differential, &cancellation, &inclusion)
                    .map_err(invalid)?,
            ))
        })
    }

    fn reduce_previous_t_map(&self, inclusion: &Self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |map| {
            let inclusion = runtime::Extract::extract(inclusion.kind.clone())?;
            Ok(wrap_module_homomorphism(
                qpa_gbnp::reduce_previous_t_map(map, &inclusion).map_err(invalid)?,
            ))
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

    fn images(
        &self,
        py: Python<'_>,
        elements: Vec<Py<PyModuleElement>>,
    ) -> PyResult<Vec<PyModuleElement>> {
        let elements = clone_module_elements(py, elements);
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_elements(
                homomorphism
                    .images_set(&runtime::extract_all(elements)?)
                    .map_err(invalid)?,
            )
        ))
    }

    fn push_out_with(&self, second: &Self) -> PyResult<Option<(Self, Self)>> {
        dispatch_pair_ref!(
            &self.kind,
            &second.kind,
            ModuleHomomorphismKind,
            |first, second| Ok(wrap_module_homomorphism_pair(
                qpa_core::push_out(first, second).map_err(invalid)?,
            )),
            Err(invalid(
                "homomorphisms belong to different coefficient fields"
            ))
        )
    }

    fn pull_back_with(&self, second: &Self) -> PyResult<Option<(Self, Self)>> {
        dispatch_pair_ref!(
            &self.kind,
            &second.kind,
            ModuleHomomorphismKind,
            |first, second| Ok(wrap_module_homomorphism_pair(
                qpa_core::pull_back(first, second).map_err(invalid)?,
            )),
            Err(invalid(
                "homomorphisms belong to different coefficient fields"
            ))
        )
    }

    fn lift_inclusion(&self, second: &Self) -> PyResult<Option<Self>> {
        dispatch_pair_ref!(
            &self.kind,
            &second.kind,
            ModuleHomomorphismKind,
            |first, second| Ok(qpa_core::lifting_inclusion_morphisms(first, second)
                .map_err(invalid)?
                .map(wrap_module_homomorphism)),
            Err(invalid(
                "homomorphisms belong to different coefficient fields"
            ))
        )
    }

    fn lift_from_projective(&self, second: &Self) -> PyResult<Option<Self>> {
        dispatch_pair_ref!(
            &self.kind,
            &second.kind,
            ModuleHomomorphismKind,
            |first, second| Ok(qpa_gbnp::lifting_morphism_from_projective(first, second)
                .map_err(invalid)?
                .map(wrap_module_homomorphism)),
            Err(invalid(
                "homomorphisms belong to different coefficient fields"
            ))
        )
    }

    #[getter]
    fn is_right_minimal(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            qpa_core::is_right_minimal(homomorphism).map_err(invalid)
        })
    }

    #[getter]
    fn is_left_minimal(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            qpa_core::is_left_minimal(homomorphism).map_err(invalid)
        })
    }

    fn right_minimal_version(&self) -> PyResult<(Self, Vec<PyMatrixModule>)> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            let (minimal, removed) =
                qpa_core::right_minimal_version(homomorphism).map_err(invalid)?;
            Ok((
                wrap_module_homomorphism(minimal),
                removed.into_iter().map(wrap_module).collect(),
            ))
        })
    }

    fn left_minimal_version(&self) -> PyResult<(Self, Vec<PyMatrixModule>)> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            let (minimal, removed) =
                qpa_core::left_minimal_version(homomorphism).map_err(invalid)?;
            Ok((
                wrap_module_homomorphism(minimal),
                removed.into_iter().map(wrap_module).collect(),
            ))
        })
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

    fn kernel_inclusion(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(qpa_core::kernel_inclusion(homomorphism).map_err(invalid)?,)
        ))
    }

    fn image_projection(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(qpa_core::image_projection(homomorphism).map_err(invalid)?,)
        ))
    }

    fn image_inclusion(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(qpa_core::image_inclusion(homomorphism).map_err(invalid)?,)
        ))
    }

    fn image_projection_inclusion(&self) -> PyResult<(Self, Self)> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            let (projection, inclusion) =
                qpa_core::image_projection_inclusion(homomorphism).map_err(invalid)?;
            Ok((
                wrap_module_homomorphism(projection),
                wrap_module_homomorphism(inclusion),
            ))
        })
    }

    fn cokernel_projection(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(qpa_core::cokernel_projection(homomorphism).map_err(invalid)?,)
        ))
    }

    fn top(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module_homomorphism(
                qpa_core::top_of_module_homomorphism(homomorphism).map_err(invalid)?,
            )
        ))
    }

    fn intersection(&self, right: &Self) -> PyResult<(Self, Self, Self)> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleHomomorphismKind,
            |left, right| {
                let (to_right, to_left, inclusion) =
                    qpa_core::intersection_of_submodules(left, right).map_err(invalid)?;
                Ok((
                    wrap_module_homomorphism(to_right),
                    wrap_module_homomorphism(to_left),
                    wrap_module_homomorphism(inclusion),
                ))
            },
            Err(invalid(
                "submodule maps belong to different coefficient fields",
            ))
        )
    }

    fn sum(&self, right: &Self) -> PyResult<(Self, Self, Self)> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleHomomorphismKind,
            |left, right| {
                let (inclusion, from_left, from_right) =
                    qpa_core::sum_of_submodules(left, right).map_err(invalid)?;
                Ok((
                    wrap_module_homomorphism(inclusion),
                    wrap_module_homomorphism(from_left),
                    wrap_module_homomorphism(from_right),
                ))
            },
            Err(invalid(
                "submodule maps belong to different coefficient fields",
            ))
        )
    }

    fn short_exact_sequence_with(&self, second: &Self) -> PyResult<PyModuleComplex> {
        let kind = map_pair_ref!(
            &self.kind,
            &second.kind,
            ModuleHomomorphismKind => ModuleComplexKind,
            |first, second| qpa_gbnp::short_exact_sequence(
                qpa_gbnp::cat_of_right_algebra_modules(first.source().right_acting_algebra())
                    .map_err(invalid)?,
                first.clone(),
                second.clone(),
            )
            .map_err(invalid)?,
            invalid("homomorphisms belong to different coefficient fields")
        )?;
        Ok(PyModuleComplex { kind })
    }

    fn comparison_lifting(
        &self,
        projective_complex: &PyModuleComplex,
        exact_complex: &PyModuleComplex,
    ) -> PyResult<PyModuleChainMap> {
        let kind = match (&self.kind, &projective_complex.kind, &exact_complex.kind) {
            (
                ModuleHomomorphismKind::Rational(morphism),
                ModuleComplexKind::Rational(projective),
                ModuleComplexKind::Rational(exact),
            ) => ModuleChainMapKind::Rational(
                qpa_gbnp::comparison_lifting(morphism, projective, exact).map_err(invalid)?,
            ),
            (
                ModuleHomomorphismKind::Binary(morphism),
                ModuleComplexKind::Binary(projective),
                ModuleComplexKind::Binary(exact),
            ) => ModuleChainMapKind::Binary(
                qpa_gbnp::comparison_lifting(morphism, projective, exact).map_err(invalid)?,
            ),
            (
                ModuleHomomorphismKind::Prime(morphism),
                ModuleComplexKind::Prime(projective),
                ModuleComplexKind::Prime(exact),
            ) => ModuleChainMapKind::Prime(
                qpa_gbnp::comparison_lifting(morphism, projective, exact).map_err(invalid)?,
            ),
            _ => {
                return Err(invalid(
                    "morphism and complexes belong to different coefficient fields",
                ));
            }
        };
        Ok(PyModuleChainMap { kind })
    }

    fn comparison_lifting_to_projective_resolution(&self) -> PyResult<PyModuleChainMap> {
        let kind = map_ref!(
            &self.kind,
            ModuleHomomorphismKind => ModuleChainMapKind,
            |morphism| qpa_gbnp::comparison_lifting_to_projective_resolution(morphism)
                .map_err(invalid)?
        );
        Ok(PyModuleChainMap { kind })
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

#[pymethods]
impl PyModuleSubmodule {
    #[getter]
    fn parent(&self) -> PyMatrixModule {
        let kind = map_ref!(
            &self.kind,
            ModuleSubmoduleKind => MatrixModuleKind,
            |submodule| submodule.parent().clone()
        );
        PyMatrixModule { kind }
    }

    #[getter]
    fn generators(&self) -> Vec<PyModuleElement> {
        dispatch_ref!(&self.kind, ModuleSubmoduleKind, |submodule| {
            wrap_module_elements(submodule.generators().to_vec())
        })
    }

    fn canonical_basis_elements(&self) -> PyResult<Vec<PyModuleElement>> {
        dispatch_ref!(&self.kind, ModuleSubmoduleKind, |submodule| Ok(
            wrap_module_elements(submodule.canonical_basis_elements().map_err(invalid)?,)
        ))
    }

    fn as_module(&self) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, ModuleSubmoduleKind, |submodule| Ok(
            wrap_module(submodule.as_module().map_err(invalid)?,)
        ))
    }

    fn is_identical(&self, right: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleSubmoduleKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __eq__(&self, right: &Self) -> bool {
        self.is_identical(right)
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, ModuleSubmoduleKind, |submodule| format!(
            "<submodule generated by {} element(s) in {}>",
            submodule.generators().len(),
            submodule.parent().qpa_print_string(),
        ))
    }
}

fn matrix_shapes<F: numerica::domains::Field>(matrices: &[DenseMatrix<F>]) -> Vec<(usize, usize)> {
    matrices
        .iter()
        .map(|matrix| (matrix.nrows(), matrix.ncols()))
        .collect()
}
