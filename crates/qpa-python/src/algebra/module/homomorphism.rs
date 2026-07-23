//! Python methods for matrix-module homomorphisms.

use pyo3::pymethods;

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

fn intersection_list_for<F: runtime::PythonField>(
    homomorphisms: Vec<ModuleHomomorphismKind>,
) -> PyResult<PyModuleHomomorphism>
where
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>
        + runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    Ok(wrap_module_homomorphism(
        qpa_core::intersection_of_submodules_list(&runtime::extract_all(homomorphisms)?)
            .map_err(invalid)?,
    ))
}

fn sum_list_for<F: runtime::PythonField>(
    homomorphisms: Vec<ModuleHomomorphismKind>,
) -> PyResult<PyModuleHomomorphism>
where
    ModuleHomomorphismKind: runtime::Extract<PathAlgebraModuleHomomorphism<F>>
        + runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    Ok(wrap_module_homomorphism(
        qpa_core::sum_of_submodules_list(&runtime::extract_all(homomorphisms)?).map_err(invalid)?,
    ))
}

#[pymethods]
impl PyModuleHomomorphism {
    #[getter]
    #[allow(clippy::unused_self)]
    fn is_algebra_module_homomorphism(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_path_algebra_mat_module_homomorphism_rep(&self) -> bool {
        true
    }

    #[getter]
    fn path_algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            conversion::wrap_quiver_algebra(homomorphism.source().right_acting_algebra().clone())
        })
    }

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

    fn endomorphism_matrix(&self) -> PyResult<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            matrix_rows(&qpa_core::from_hom_mm_to_end_m(homomorphism).map_err(invalid)?,)
        ))
    }

    fn underlying_linear_matrix(&self) -> PyResult<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            matrix_rows(&qpa_core::underlying_linear_map(homomorphism).map_err(invalid)?,)
        ))
    }

    fn projective_matrix(&self) -> PyResult<Vec<Vec<PyAlgebraElement>>> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            qpa_gbnp::matrix_of_homomorphism_between_projectives(homomorphism)
                .map_err(invalid)?
                .into_iter()
                .map(wrap_algebra_elements)
                .collect()
        ))
    }

    fn restrict_scalars(
        &self,
        algebra_homomorphism: &PyAlgebraHomomorphism,
    ) -> PyResult<PyModuleHomomorphism> {
        let algebra_homomorphism = algebra_homomorphism.kind.clone();
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |module_homomorphism| {
            restrict_module_homomorphism_for(module_homomorphism, algebra_homomorphism)
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
            (ModuleHomomorphismKind::Extension(map), ModuleElementKind::Extension(element)) => {
                ModuleElementKind::Extension(map.image(element).map_err(invalid)?)
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

    fn more_right_minimal_version(&self) -> PyResult<PyMoreMinimalVersionResult> {
        let kind = map_ref!(
            &self.kind,
            ModuleHomomorphismKind => MoreMinimalVersionResultKind,
            |homomorphism| qpa_core::more_right_minimal_version(homomorphism).map_err(invalid)?
        );
        Ok(PyMoreMinimalVersionResult { kind })
    }

    fn more_left_minimal_version(&self) -> PyResult<PyMoreMinimalVersionResult> {
        let kind = map_ref!(
            &self.kind,
            ModuleHomomorphismKind => MoreMinimalVersionResultKind,
            |homomorphism| qpa_core::more_left_minimal_version(homomorphism).map_err(invalid)?
        );
        Ok(PyMoreMinimalVersionResult { kind })
    }

    fn morphism_on_kernel(
        &self,
        second: &Self,
        beta: &Self,
        alpha: &Self,
    ) -> PyResult<Option<Self>> {
        let second = second.kind.clone();
        let beta = beta.kind.clone();
        let alpha = alpha.kind.clone();
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |first| {
            induced_map_for(first, second, beta, alpha, qpa_core::morphism_on_kernel)
        })
    }

    fn morphism_on_image(
        &self,
        second: &Self,
        beta: &Self,
        alpha: &Self,
    ) -> PyResult<Option<Self>> {
        let second = second.kind.clone();
        let beta = beta.kind.clone();
        let alpha = alpha.kind.clone();
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |first| {
            induced_map_for(first, second, beta, alpha, qpa_core::morphism_on_image)
        })
    }

    fn morphism_on_cokernel(
        &self,
        second: &Self,
        beta: &Self,
        alpha: &Self,
    ) -> PyResult<Option<Self>> {
        let second = second.kind.clone();
        let beta = beta.kind.clone();
        let alpha = alpha.kind.clone();
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |first| {
            induced_map_for(first, second, beta, alpha, qpa_core::morphism_on_cokernel)
        })
    }

    #[staticmethod]
    fn multiply_lists(
        py: Python<'_>,
        projections: Vec<Py<PyModuleHomomorphism>>,
        matrix: Vec<Vec<Py<PyModuleHomomorphism>>>,
        inclusions: Vec<Py<PyModuleHomomorphism>>,
    ) -> PyResult<PyModuleHomomorphism> {
        let projections = clone_module_homomorphisms(py, projections);
        let matrix = matrix
            .into_iter()
            .map(|row| clone_module_homomorphisms(py, row))
            .collect::<Vec<_>>();
        let inclusions = clone_module_homomorphisms(py, inclusions);
        match projections.first() {
            Some(ModuleHomomorphismKind::Rational(_)) | None => {
                multiply_lists_of_maps_for::<RationalField>(projections, matrix, inclusions)
            }
            Some(ModuleHomomorphismKind::Binary(_)) => {
                multiply_lists_of_maps_for::<Z2>(projections, matrix, inclusions)
            }
            Some(ModuleHomomorphismKind::Prime(_)) => {
                multiply_lists_of_maps_for::<Zp>(projections, matrix, inclusions)
            }
            Some(ModuleHomomorphismKind::Extension(_)) => {
                multiply_lists_of_maps_for::<GaloisField>(projections, matrix, inclusions)
            }
        }
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
            (ModuleHomomorphismKind::Extension(map), ModuleElementKind::Extension(element)) => {
                Ok(map
                    .preimage_representative(element)
                    .map_err(invalid)?
                    .map(wrap_module_element))
            }
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

    /// Apply this map first and `right` second, following QPA's row-action order.
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

    /// Return the kernel module computed by QPA's map routine.
    fn kernel(&self) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module(qpa_core::kernel_of_module_homomorphism(homomorphism).map_err(invalid)?,)
        ))
    }

    /// Return the image module computed by QPA's map routine.
    fn image_module(&self) -> PyResult<PyMatrixModule> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| Ok(
            wrap_module(qpa_core::image_of_module_homomorphism(homomorphism).map_err(invalid)?,)
        ))
    }

    /// Return the cokernel module computed by QPA's map routine.
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

    #[getter]
    fn kernel_of(&self) -> Option<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism.kernel_of_what().map(wrap_module_homomorphism)
        })
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

    #[getter]
    fn image_of(&self) -> Option<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism.image_of_what().map(wrap_module_homomorphism)
        })
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

    #[getter]
    fn cokernel_of(&self) -> Option<Self> {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism
                .cokernel_of_what()
                .map(wrap_module_homomorphism)
        })
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

    #[staticmethod]
    fn intersection_list(
        py: Python<'_>,
        homomorphisms: Vec<Py<PyModuleHomomorphism>>,
    ) -> PyResult<Self> {
        let homomorphisms = clone_module_homomorphisms(py, homomorphisms);
        match homomorphisms.first() {
            Some(ModuleHomomorphismKind::Rational(_)) | None => {
                intersection_list_for::<RationalField>(homomorphisms)
            }
            Some(ModuleHomomorphismKind::Binary(_)) => intersection_list_for::<Z2>(homomorphisms),
            Some(ModuleHomomorphismKind::Prime(_)) => intersection_list_for::<Zp>(homomorphisms),
            Some(ModuleHomomorphismKind::Extension(_)) => {
                intersection_list_for::<GaloisField>(homomorphisms)
            }
        }
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

    #[staticmethod]
    fn sum_list(py: Python<'_>, homomorphisms: Vec<Py<PyModuleHomomorphism>>) -> PyResult<Self> {
        let homomorphisms = clone_module_homomorphisms(py, homomorphisms);
        match homomorphisms.first() {
            Some(ModuleHomomorphismKind::Rational(_)) | None => {
                sum_list_for::<RationalField>(homomorphisms)
            }
            Some(ModuleHomomorphismKind::Binary(_)) => sum_list_for::<Z2>(homomorphisms),
            Some(ModuleHomomorphismKind::Prime(_)) => sum_list_for::<Zp>(homomorphisms),
            Some(ModuleHomomorphismKind::Extension(_)) => {
                sum_list_for::<GaloisField>(homomorphisms)
            }
        }
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
            (
                ModuleHomomorphismKind::Extension(morphism),
                ModuleComplexKind::Extension(projective),
                ModuleComplexKind::Extension(exact),
            ) => ModuleChainMapKind::Extension(
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

    fn view(&self) -> String {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism.to_string()
        })
    }

    fn print_obj(&self) -> String {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism.qpa_print_string()
        })
    }

    fn display(&self) -> String {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism.qpa_display_string_with(|value| {
                homomorphism
                    .source()
                    .right_acting_algebra()
                    .field()
                    .scalar_text(value)
            })
        })
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, ModuleHomomorphismKind, |homomorphism| {
            homomorphism.qpa_print_string()
        })
    }
}
