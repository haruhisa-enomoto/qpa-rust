//! Structured Python results for module and endomorphism computations.

use pyo3::pymethods;

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

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
impl PyAllSubmodulesOfModuleResult {
    #[getter]
    fn is_flat(&self) -> bool {
        dispatch_ref!(&self.kind, AllSubmodulesOfModuleResultKind, |result| result
            .flat()
            .is_some())
    }

    #[getter]
    fn flat(&self) -> Option<Vec<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, AllSubmodulesOfModuleResultKind, |result| result
            .flat()
            .map(|maps| wrap_module_homomorphisms(maps.to_vec())))
    }

    #[getter]
    fn by_dimension(&self) -> Option<Vec<Vec<PyModuleHomomorphism>>> {
        dispatch_ref!(&self.kind, AllSubmodulesOfModuleResultKind, |result| result
            .by_dimension()
            .map(|levels| levels
                .iter()
                .cloned()
                .map(wrap_module_homomorphisms)
                .collect()))
    }

    #[getter]
    fn flattened(&self) -> Vec<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, AllSubmodulesOfModuleResultKind, |result| {
            wrap_module_homomorphisms(result.clone().into_flattened())
        })
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
impl PyMoreMinimalVersionResult {
    #[getter]
    fn is_unchanged(&self) -> bool {
        dispatch_ref!(&self.kind, MoreMinimalVersionResultKind, |result| matches!(
            result,
            qpa_core::MoreMinimalVersionResult::Unchanged(_)
        ))
    }

    #[getter]
    fn minimal_part(&self) -> PyModuleHomomorphism {
        dispatch_ref!(
            &self.kind,
            MoreMinimalVersionResultKind,
            |result| match result {
                qpa_core::MoreMinimalVersionResult::Unchanged(homomorphism) => {
                    wrap_module_homomorphism(homomorphism.clone())
                }
                qpa_core::MoreMinimalVersionResult::Reduced { minimal_part, .. } => {
                    wrap_module_homomorphism(minimal_part.clone())
                }
            }
        )
    }

    #[getter]
    fn removed_part(&self) -> Option<PyModuleHomomorphism> {
        dispatch_ref!(
            &self.kind,
            MoreMinimalVersionResultKind,
            |result| match result {
                qpa_core::MoreMinimalVersionResult::Unchanged(_) => None,
                qpa_core::MoreMinimalVersionResult::Reduced { removed_part, .. } => {
                    Some(wrap_module_homomorphism(removed_part.clone()))
                }
            }
        )
    }
}

#[pymethods]
impl PyCommonDirectSummandResult {
    #[getter]
    fn status(&self) -> &'static str {
        dispatch_ref!(
            &self.kind,
            CommonDirectSummandResultKind,
            |result| match result {
                qpa_core::CommonDirectSummandResult::DifferentAlgebras => "different_algebras",
                qpa_core::CommonDirectSummandResult::NoCommonSummand => "none",
                qpa_core::CommonDirectSummandResult::Found(_) => "found",
            }
        )
    }

    #[getter]
    fn source_common(&self) -> Option<PyMatrixModule> {
        dispatch_ref!(
            &self.kind,
            CommonDirectSummandResultKind,
            |result| match result {
                qpa_core::CommonDirectSummandResult::Found(common) => {
                    Some(wrap_module(common.source_common().clone()))
                }
                _ => None,
            }
        )
    }

    #[getter]
    fn source_complement(&self) -> Option<PyMatrixModule> {
        dispatch_ref!(
            &self.kind,
            CommonDirectSummandResultKind,
            |result| match result {
                qpa_core::CommonDirectSummandResult::Found(common) => {
                    Some(wrap_module(common.source_complement().clone()))
                }
                _ => None,
            }
        )
    }

    #[getter]
    fn range_common(&self) -> Option<PyMatrixModule> {
        dispatch_ref!(
            &self.kind,
            CommonDirectSummandResultKind,
            |result| match result {
                qpa_core::CommonDirectSummandResult::Found(common) => {
                    Some(wrap_module(common.range_common().clone()))
                }
                _ => None,
            }
        )
    }

    #[getter]
    fn range_complement(&self) -> Option<PyMatrixModule> {
        dispatch_ref!(
            &self.kind,
            CommonDirectSummandResultKind,
            |result| match result {
                qpa_core::CommonDirectSummandResult::Found(common) => {
                    Some(wrap_module(common.range_complement().clone()))
                }
                _ => None,
            }
        )
    }
}

#[pymethods]
impl PyMaximalCommonDirectSummand {
    #[getter]
    fn common_summands(&self) -> Vec<PyMatrixModule> {
        dispatch_ref!(&self.kind, MaximalCommonDirectSummandKind, |result| result
            .common_summands()
            .iter()
            .cloned()
            .map(wrap_module)
            .collect())
    }

    #[getter]
    fn source_remainder(&self) -> PyMatrixModule {
        dispatch_ref!(&self.kind, MaximalCommonDirectSummandKind, |result| {
            wrap_module(result.source_remainder().clone())
        })
    }

    #[getter]
    fn range_remainder(&self) -> PyMatrixModule {
        dispatch_ref!(&self.kind, MaximalCommonDirectSummandKind, |result| {
            wrap_module(result.range_remainder().clone())
        })
    }
}

#[pymethods]
impl PyModuleHomomorphismBasis {
    #[getter]
    fn basis(&self) -> Vec<PyModuleHomomorphism> {
        dispatch_ref!(&self.kind, ModuleHomomorphismBasisKind, |basis| basis
            .basis()
            .iter()
            .cloned()
            .map(wrap_module_homomorphism)
            .collect())
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, ModuleHomomorphismBasisKind, |basis| basis
            .basis()
            .len())
    }

    fn coefficients(&self, homomorphism: &PyModuleHomomorphism) -> PyResult<Option<Vec<PyScalar>>> {
        let homomorphism = homomorphism.kind.clone();
        dispatch_ref!(&self.kind, ModuleHomomorphismBasisKind, |basis| {
            hom_basis_coefficients_for(basis, homomorphism)
        })
    }
}

#[pymethods]
impl PyModuleEndomorphismAlgebra {
    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(&self.kind, ModuleEndomorphismAlgebraKind, |algebra| algebra
            .field()
            .python_kind());
        PyField { kind }
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(&self.kind, ModuleEndomorphismAlgebraKind, |algebra| algebra
            .dimension())
    }

    #[getter]
    fn basis_matrices(&self) -> Vec<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, ModuleEndomorphismAlgebraKind, |algebra| algebra
            .basis_matrices()
            .iter()
            .map(matrix_rows)
            .collect())
    }

    #[getter]
    fn one(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, ModuleEndomorphismAlgebraKind, |algebra| {
            matrix_rows(algebra.one())
        })
    }

    #[getter]
    fn zero(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, ModuleEndomorphismAlgebraKind, |algebra| {
            matrix_rows(&algebra.zero())
        })
    }

    #[getter]
    fn radical_basis_matrices(&self) -> PyResult<Vec<Vec<Vec<PyScalar>>>> {
        dispatch_ref!(&self.kind, ModuleEndomorphismAlgebraKind, |algebra| {
            Ok(algebra
                .radical_basis_matrices()
                .map_err(invalid)?
                .iter()
                .map(matrix_rows)
                .collect())
        })
    }

    fn coefficients(
        &self,
        py: Python<'_>,
        matrix: Vec<Vec<Py<PyScalar>>>,
    ) -> PyResult<Option<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, ModuleEndomorphismAlgebraKind, |algebra| {
            let matrix = extract_matrix(py, algebra.field(), matrix)?;
            Ok(algebra
                .coefficients(&matrix)
                .map_err(invalid)?
                .map(|values| wrap_scalars(algebra.field(), &values)))
        })
    }

    fn matrix_from_coefficients(
        &self,
        py: Python<'_>,
        coefficients: Vec<Py<PyScalar>>,
    ) -> PyResult<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, ModuleEndomorphismAlgebraKind, |algebra| {
            let coefficients = extract_scalars(py, algebra.field(), coefficients)?;
            Ok(matrix_rows(
                &algebra
                    .matrix_from_coefficients(&coefficients)
                    .map_err(invalid)?,
            ))
        })
    }

    fn multiply(
        &self,
        py: Python<'_>,
        left: Vec<Vec<Py<PyScalar>>>,
        right: Vec<Vec<Py<PyScalar>>>,
    ) -> PyResult<Vec<Vec<PyScalar>>> {
        dispatch_ref!(&self.kind, ModuleEndomorphismAlgebraKind, |algebra| {
            let left = extract_matrix(py, algebra.field(), left)?;
            let right = extract_matrix(py, algebra.field(), right)?;
            Ok(matrix_rows(
                &algebra.multiply(&left, &right).map_err(invalid)?,
            ))
        })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            ModuleEndomorphismAlgebraKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn structure_constants_algebra(&self) -> PyResult<PyStructureConstantsAlgebra> {
        let kind = map_ref!(
            &self.kind,
            ModuleEndomorphismAlgebraKind => StructureConstantsAlgebraKind,
            |algebra| algebra.structure_constants_algebra().map_err(invalid)?
        );
        Ok(PyStructureConstantsAlgebra { kind })
    }
}

#[pymethods]
impl PyEndOfModuleAsQuiverAlgebraResult {
    #[getter]
    fn endomorphism_algebra(&self) -> PyModuleEndomorphismAlgebra {
        let kind = map_ref!(
            &self.kind,
            EndOfModuleAsQuiverAlgebraResultKind => ModuleEndomorphismAlgebraKind,
            |result| result.endomorphism_algebra().clone()
        );
        PyModuleEndomorphismAlgebra { kind }
    }

    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, EndOfModuleAsQuiverAlgebraResultKind, |result| {
            wrap_quiver_algebra(result.quiver_algebra().clone())
        })
    }

    #[getter]
    fn adjacency_matrix(&self) -> Vec<Vec<usize>> {
        dispatch_ref!(&self.kind, EndOfModuleAsQuiverAlgebraResultKind, |result| {
            result
                .adjacency_matrix()
                .iter()
                .map(|row| row.to_vec())
                .collect()
        })
    }

    #[getter]
    fn quiver(&self) -> PyQuiver {
        let inner = dispatch_ref!(&self.kind, EndOfModuleAsQuiverAlgebraResultKind, |result| {
            result.quiver_algebra().quiver().clone()
        });
        PyQuiver { inner }
    }

    #[getter]
    fn is_quotient(&self) -> bool {
        dispatch_ref!(
            &self.kind,
            EndOfModuleAsQuiverAlgebraResultKind,
            |result| matches!(result.quiver_algebra(), QuiverAlgebra::Quotient(_))
        )
    }

    #[getter]
    fn path_algebra(&self) -> Option<PyPathAlgebra> {
        dispatch_ref!(
            &self.kind,
            EndOfModuleAsQuiverAlgebraResultKind,
            |result| match result.quiver_algebra() {
                QuiverAlgebra::Path(algebra) => Some(PyPathAlgebra {
                    kind: runtime::Wrap::wrap(algebra.clone()),
                }),
                QuiverAlgebra::Quotient(_) => None,
            }
        )
    }

    #[getter]
    fn quotient_algebra(&self) -> Option<PyQuotientAlgebra> {
        dispatch_ref!(
            &self.kind,
            EndOfModuleAsQuiverAlgebraResultKind,
            |result| match result.quiver_algebra() {
                QuiverAlgebra::Path(_) => None,
                QuiverAlgebra::Quotient(algebra) => Some(PyQuotientAlgebra {
                    kind: runtime::Wrap::wrap(algebra.clone()),
                }),
            }
        )
    }
}

#[pymethods]
impl PyEndOfBasicModuleAsQuiverAlgebraResult {
    #[getter]
    fn algebra(&self) -> PyQuiverAlgebra {
        dispatch_ref!(
            &self.kind,
            EndOfBasicModuleAsQuiverAlgebraResultKind,
            |result| wrap_quiver_algebra(result.algebra().clone())
        )
    }

    #[getter]
    fn search_bound(&self) -> usize {
        dispatch_ref!(
            &self.kind,
            EndOfBasicModuleAsQuiverAlgebraResultKind,
            |result| result.search_bound()
        )
    }

    #[getter]
    fn final_search_length(&self) -> usize {
        dispatch_ref!(
            &self.kind,
            EndOfBasicModuleAsQuiverAlgebraResultKind,
            |result| result.final_search_length()
        )
    }

    #[getter]
    fn relations_may_be_incomplete(&self) -> bool {
        dispatch_ref!(
            &self.kind,
            EndOfBasicModuleAsQuiverAlgebraResultKind,
            |result| result.relations_may_be_incomplete()
        )
    }

    #[getter]
    fn quiver(&self) -> PyQuiver {
        let inner = dispatch_ref!(
            &self.kind,
            EndOfBasicModuleAsQuiverAlgebraResultKind,
            |result| result.algebra().quiver().clone()
        );
        PyQuiver { inner }
    }

    #[getter]
    fn is_quotient(&self) -> bool {
        dispatch_ref!(
            &self.kind,
            EndOfBasicModuleAsQuiverAlgebraResultKind,
            |result| matches!(result.algebra(), QuiverAlgebra::Quotient(_))
        )
    }

    #[getter]
    fn path_algebra(&self) -> Option<PyPathAlgebra> {
        dispatch_ref!(
            &self.kind,
            EndOfBasicModuleAsQuiverAlgebraResultKind,
            |result| match result.algebra() {
                QuiverAlgebra::Path(algebra) => Some(PyPathAlgebra {
                    kind: runtime::Wrap::wrap(algebra.clone()),
                }),
                QuiverAlgebra::Quotient(_) => None,
            }
        )
    }

    #[getter]
    fn quotient_algebra(&self) -> Option<PyQuotientAlgebra> {
        dispatch_ref!(
            &self.kind,
            EndOfBasicModuleAsQuiverAlgebraResultKind,
            |result| match result.algebra() {
                QuiverAlgebra::Path(_) => None,
                QuiverAlgebra::Quotient(algebra) => Some(PyQuotientAlgebra {
                    kind: runtime::Wrap::wrap(algebra.clone()),
                }),
            }
        )
    }
}

#[pymethods]
impl PyModuleEndomorphismQuotientMap {
    #[getter]
    fn source(&self) -> PyModuleEndomorphismAlgebra {
        let kind = map_ref!(
            &self.kind,
            ModuleEndomorphismQuotientMapKind => ModuleEndomorphismAlgebraKind,
            |homomorphism| homomorphism.source().clone()
        );
        PyModuleEndomorphismAlgebra { kind }
    }

    #[getter]
    fn range(&self) -> PyModuleEndomorphismQuotientAlgebra {
        let kind = map_ref!(
            &self.kind,
            ModuleEndomorphismQuotientMapKind => ModuleEndomorphismQuotientAlgebraKind,
            |homomorphism| homomorphism.range().clone()
        );
        PyModuleEndomorphismQuotientAlgebra { kind }
    }

    #[getter]
    fn kernel_basis_matrices(&self) -> Vec<Vec<Vec<PyScalar>>> {
        dispatch_ref!(
            &self.kind,
            ModuleEndomorphismQuotientMapKind,
            |homomorphism| homomorphism
                .kernel_basis_matrices()
                .iter()
                .map(matrix_rows)
                .collect()
        )
    }

    fn image(
        &self,
        py: Python<'_>,
        matrix: Vec<Vec<Py<PyScalar>>>,
    ) -> PyResult<PyModuleEndomorphismQuotientElement> {
        dispatch_ref!(
            &self.kind,
            ModuleEndomorphismQuotientMapKind,
            |homomorphism| {
                let field = homomorphism.source().field();
                Ok(PyModuleEndomorphismQuotientElement {
                    kind: runtime::Wrap::wrap(
                        homomorphism
                            .image(&extract_matrix(py, field, matrix)?)
                            .map_err(invalid)?,
                    ),
                })
            }
        )
    }
}

#[pymethods]
impl PyModuleEndomorphismQuotientAlgebra {
    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(
            &self.kind,
            ModuleEndomorphismQuotientAlgebraKind,
            |algebra| algebra.field().python_kind()
        );
        PyField { kind }
    }

    #[getter]
    fn dimension(&self) -> usize {
        dispatch_ref!(
            &self.kind,
            ModuleEndomorphismQuotientAlgebraKind,
            |algebra| algebra.dimension()
        )
    }

    #[getter]
    fn basis_elements(&self) -> Vec<PyModuleEndomorphismQuotientElement> {
        dispatch_ref!(
            &self.kind,
            ModuleEndomorphismQuotientAlgebraKind,
            |algebra| algebra
                .basis_elements()
                .into_iter()
                .map(|element| PyModuleEndomorphismQuotientElement {
                    kind: runtime::Wrap::wrap(element),
                })
                .collect()
        )
    }

    fn zero(&self) -> PyModuleEndomorphismQuotientElement {
        let kind = map_ref!(
            &self.kind,
            ModuleEndomorphismQuotientAlgebraKind => ModuleEndomorphismQuotientElementKind,
            |algebra| algebra.zero()
        );
        PyModuleEndomorphismQuotientElement { kind }
    }

    fn one(&self) -> PyResult<Option<PyModuleEndomorphismQuotientElement>> {
        dispatch_ref!(
            &self.kind,
            ModuleEndomorphismQuotientAlgebraKind,
            |algebra| Ok(algebra.one().map_err(invalid)?.map(|element| {
                PyModuleEndomorphismQuotientElement {
                    kind: runtime::Wrap::wrap(element),
                }
            }))
        )
    }

    fn multiply(
        &self,
        left: &PyModuleEndomorphismQuotientElement,
        right: &PyModuleEndomorphismQuotientElement,
    ) -> PyResult<PyModuleEndomorphismQuotientElement> {
        match (&self.kind, &left.kind, &right.kind) {
            (
                ModuleEndomorphismQuotientAlgebraKind::Rational(algebra),
                ModuleEndomorphismQuotientElementKind::Rational(left),
                ModuleEndomorphismQuotientElementKind::Rational(right),
            ) => Ok(PyModuleEndomorphismQuotientElement {
                kind: ModuleEndomorphismQuotientElementKind::Rational(
                    algebra.multiply(left, right).map_err(invalid)?,
                ),
            }),
            (
                ModuleEndomorphismQuotientAlgebraKind::Binary(algebra),
                ModuleEndomorphismQuotientElementKind::Binary(left),
                ModuleEndomorphismQuotientElementKind::Binary(right),
            ) => Ok(PyModuleEndomorphismQuotientElement {
                kind: ModuleEndomorphismQuotientElementKind::Binary(
                    algebra.multiply(left, right).map_err(invalid)?,
                ),
            }),
            (
                ModuleEndomorphismQuotientAlgebraKind::Prime(algebra),
                ModuleEndomorphismQuotientElementKind::Prime(left),
                ModuleEndomorphismQuotientElementKind::Prime(right),
            ) => Ok(PyModuleEndomorphismQuotientElement {
                kind: ModuleEndomorphismQuotientElementKind::Prime(
                    algebra.multiply(left, right).map_err(invalid)?,
                ),
            }),
            _ => Err(invalid(
                "quotient algebra and elements belong to different coefficient fields",
            )),
        }
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            ModuleEndomorphismQuotientAlgebraKind,
            |left, right| left.is_identical(right),
            false
        )
    }
}

#[pymethods]
impl PyModuleEndomorphismQuotientElement {
    #[getter]
    fn algebra(&self) -> PyModuleEndomorphismQuotientAlgebra {
        let kind = map_ref!(
            &self.kind,
            ModuleEndomorphismQuotientElementKind => ModuleEndomorphismQuotientAlgebraKind,
            |element| element.algebra().clone()
        );
        PyModuleEndomorphismQuotientAlgebra { kind }
    }

    #[getter]
    fn coordinates(&self) -> Vec<PyScalar> {
        dispatch_ref!(
            &self.kind,
            ModuleEndomorphismQuotientElementKind,
            |element| wrap_scalars(element.algebra().field(), element.coordinates())
        )
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(
            &self.kind,
            ModuleEndomorphismQuotientElementKind,
            |element| element.is_zero()
        )
    }
}
