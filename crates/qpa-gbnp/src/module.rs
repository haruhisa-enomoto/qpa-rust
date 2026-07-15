//! QPA modules constructed from enveloping algebras.
//!
//! This begins the source-order translation of `lib/pathalgtensor.gi`'s
//! module section. The regular bimodule method deliberately uses QPA's
//! canonical-basis partition and its checked path-module constructor.

use numerica::domains::Field;
use qpa_core::{
    DenseMatrix, MatrixError, NamedArrowMatrix, PathAlgebraError, PathAlgebraModule,
    PathAlgebraModuleElement, PathAlgebraModuleError, PathAlgebraModuleHomomorphism,
    PathAlgebraTensorError, PathEndpoint, PathError, QuiverAlgebra, QuiverAlgebraElement,
    QuiverAlgebraElementError, QuiverProduct, QuiverProductError, QuotientError,
    direct_sum_of_qpa_modules, direct_sum_projections, hom_over_algebra, include_in_path_algebra,
    loewy_length_of_module, minimal_generating_set_of_module, project_from_product_quiver,
    right_module_hom_over_algebra, right_module_over_path_algebra,
    right_module_over_path_algebra_with_tensor, top_of_module, zero_module_mapping,
};
use thiserror::Error;

use crate::{
    PathalgGbnpError, TensorProductError, ensure_quotient_groebner_basis, enveloping_algebra,
    opposite_path_algebra, quotient_is_finite_dimensional,
};

/// Translate QPA's cached `AlgebraAsModuleOverEnvelopingAlgebra(A)`
/// attribute (and its `AlgebraAsModuleOfEnvelopingAlgebra` synonym).
///
/// The canonical basis is partitioned into the spaces `e_i A e_j`. For each
/// product-quiver arrow, the source method acts on those basis vectors by the
/// opposite of its first projection on the left and its second projection on
/// the right, then selects canonical-basis coordinates in the target space.
/// The resulting row-action matrices are passed through QPA's checked
/// `RightModuleOverPathAlgebra` translation.
///
/// `Ok(None)` is QPA's `fail` for an infinite-dimensional algebra (or the
/// checked constructor's theoretically possible relation failure).
///
/// # Errors
///
/// Propagates the configured GBNP finite-dimensionality/basis path,
/// enveloping-algebra construction, product projections, exact algebra and
/// matrix arithmetic, and checked module construction.
pub fn algebra_as_module_over_enveloping_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_regular_enveloping_module() {
        return Ok(Some(cached));
    }
    if algebra.regular_enveloping_module_attribute_is_fail() {
        return Ok(None);
    }

    let Some(basis_vectors) = finite_canonical_basis(algebra)? else {
        algebra.mark_regular_enveloping_module_attribute_fail();
        return Ok(None);
    };

    let enveloping = enveloping_algebra(algebra)?;
    let product = enveloping.product_quiver();
    let vertices = algebra.quiver().vertices().collect::<Vec<_>>();
    let mut vector_spaces = vec![vec![Vec::<usize>::new(); vertices.len()]; vertices.len()];

    let vertex_elements = vertices
        .iter()
        .map(|vertex| include_in_path_algebra(&vertex.path(), algebra))
        .collect::<Result<Vec<_>, PathAlgebraTensorError>>()?;
    for (left_index, left) in vertex_elements.iter().enumerate() {
        for (right_index, right) in vertex_elements.iter().enumerate() {
            for (basis_index, basis_vector) in basis_vectors.iter().enumerate() {
                let component = left.multiply(basis_vector)?.multiply(right)?;
                if !component.is_zero() {
                    vector_spaces[left_index][right_index].push(basis_index);
                }
            }
        }
    }
    let dimension_vector = vector_spaces
        .iter()
        .flat_map(|row| row.iter().map(Vec::len))
        .collect::<Vec<_>>();

    let mut module_specification = Vec::new();
    for arrow in enveloping.algebra().quiver().arrows() {
        let components = [
            project_from_product_quiver(1, &arrow.path(), product)?,
            project_from_product_quiver(2, &arrow.path(), product)?,
        ];
        let source_indices = components.each_ref().map(|component| {
            let PathEndpoint::Vertex(vertex) = component.source() else {
                unreachable!("a nonzero projected product arrow has a vertex source")
            };
            vertex.index().get() - 1
        });
        let target_indices = components.each_ref().map(|component| {
            let PathEndpoint::Vertex(vertex) = component.target() else {
                unreachable!("a nonzero projected product arrow has a vertex target")
            };
            vertex.index().get() - 1
        });
        let source_space = &vector_spaces[source_indices[0]][source_indices[1]];
        let target_space = &vector_spaces[target_indices[0]][target_indices[1]];
        if source_space.is_empty() || target_space.is_empty() {
            continue;
        }

        let left_path = components[0].quiver().opposite_path(&components[0])?;
        let left = include_in_path_algebra(&left_path, algebra)?;
        let right = include_in_path_algebra(&components[1], algebra)?;
        let rows = source_space
            .iter()
            .map(|&basis_index| {
                let image = left
                    .multiply(&basis_vectors[basis_index])?
                    .multiply(&right)?;
                let coefficients = canonical_basis_coefficients(algebra, &basis_vectors, &image)?;
                Ok(target_space
                    .iter()
                    .map(|&target_index| coefficients[target_index].clone())
                    .collect::<Vec<_>>())
            })
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
        let matrix = DenseMatrix::from_rows(rows, algebra.field().clone())?;
        module_specification.push(NamedArrowMatrix::new(arrow.name(), matrix));
    }

    let Some(module) = right_module_over_path_algebra_with_tensor(
        enveloping.algebra().clone(),
        dimension_vector,
        &module_specification,
        enveloping,
    )?
    else {
        algebra.mark_regular_enveloping_module_attribute_fail();
        return Ok(None);
    };
    Ok(Some(algebra.cache_regular_enveloping_module(&module)))
}

/// QPA's spelling synonym for
/// [`algebra_as_module_over_enveloping_algebra`].
///
/// # Errors
///
/// Propagates the identical source method failures.
pub fn algebra_as_module_of_enveloping_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    algebra_as_module_over_enveloping_algebra(algebra)
}

/// Translate QPA's algebra `LoewyLength` method through its indecomposable
/// projective representations and repeated `RadicalOfModule` calls.
///
/// # Errors
///
/// Requires a finite path algebra or finite admissible quotient and
/// propagates canonical-basis, exact module, and radical failures.
pub fn loewy_length_of_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<usize, EnvelopingModuleError> {
    let projectives = indec_projective_modules(algebra)?;
    let mut maximum = 0;
    for projective in projectives {
        maximum = maximum.max(loewy_length_of_module(projective)?);
    }
    Ok(maximum)
}

/// Vertexwise canonical bases returned by QPA's `BasisOfProjectives`.
///
/// The outer index is the generating vertex and the inner index is the right
/// support vertex; elements retain QPA's canonical right-ideal basis order.
pub type ProjectiveModuleBases<F> = Vec<Vec<Vec<QuiverAlgebraElement<F>>>>;

struct ProjectiveBasisData<F: Field> {
    basis: Vec<QuiverAlgebraElement<F>>,
    intervals: Vec<Vec<Vec<usize>>>,
}

/// Translate QPA's two `BasisOfProjectives(A)` methods.
///
/// `Ok(None)` is the printed `fail` result for an infinite path algebra or a
/// quotient without finitely many nontips.
///
/// # Errors
///
/// Propagates the configured GBNP finite-basis path.
pub fn basis_of_projectives<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<ProjectiveModuleBases<F>>, EnvelopingModuleError> {
    let Some(data) = projective_basis_data(algebra)? else {
        return Ok(None);
    };
    Ok(Some(
        data.intervals
            .iter()
            .map(|projective| {
                projective
                    .iter()
                    .map(|interval| {
                        interval
                            .iter()
                            .map(|&index| data.basis[index].clone())
                            .collect()
                    })
                    .collect()
            })
            .collect(),
    ))
}

/// Translate QPA's cached `IndecProjectiveModules(A)` attribute.
///
/// Modules are returned in vertex order, with the exact canonical-right-ideal
/// action matrices and QPA's distinguished vertex generator installed.
///
/// # Errors
///
/// Requires a finite path algebra or finite admissible quotient and propagates
/// canonical-basis, exact multiplication, coordinate, and checked-module
/// construction failures.
pub fn indec_projective_modules<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Vec<PathAlgebraModule<F>>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_indec_projective_modules() {
        return Ok(cached);
    }
    let Some(data) = projective_basis_data(algebra)? else {
        return Err(EnvelopingModuleError::InfiniteLoewyAlgebra);
    };
    if matches!(algebra, QuiverAlgebra::Quotient(quotient) if !quotient.is_admissible()?) {
        return Err(EnvelopingModuleError::NonadmissibleLoewyAlgebra);
    }
    let basis = data.basis;
    let vertex_count = algebra.quiver().number_of_vertices();
    let arrows = algebra
        .quiver()
        .arrows()
        .map(|arrow| include_in_path_algebra(&arrow.path(), algebra))
        .collect::<Result<Vec<_>, _>>()?;
    let mut projectives = Vec::with_capacity(vertex_count);

    for projective_vertex in 0..vertex_count {
        let intervals = &data.intervals[projective_vertex];
        let dimensions = intervals.iter().map(Vec::len).collect::<Vec<_>>();
        let mut specification = Vec::new();
        for (arrow, arrow_element) in algebra.quiver().arrows().zip(&arrows) {
            let source = arrow.source().index().get() - 1;
            let target = arrow.target().index().get() - 1;
            if dimensions[source] == 0 || dimensions[target] == 0 {
                continue;
            }
            let rows = intervals[source]
                .iter()
                .map(|&basis_index| {
                    let image = basis[basis_index].multiply(arrow_element)?;
                    let coefficients = canonical_basis_coefficients(algebra, &basis, &image)?;
                    Ok(intervals[target]
                        .iter()
                        .map(|&index| coefficients[index].clone())
                        .collect::<Vec<_>>())
                })
                .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
            specification.push(NamedArrowMatrix::new(
                arrow.name(),
                DenseMatrix::from_rows(rows, algebra.field().clone())?,
            ));
        }
        let projective =
            right_module_over_path_algebra(algebra.clone(), dimensions.clone(), &specification)?
                .ok_or(EnvelopingModuleError::GeneratedProjectiveViolatesRelations)?;
        let mut generator_coordinates = dimensions
            .iter()
            .map(|&dimension| vec![algebra.field().zero(); dimension])
            .collect::<Vec<_>>();
        let generator_coordinate = intervals[projective_vertex]
            .iter()
            .position(|&basis_index| basis[basis_index].leading_monomial().length() == Some(0))
            .ok_or(EnvelopingModuleError::MissingProjectiveVertexGenerator {
                vertex: projective_vertex + 1,
            })?;
        generator_coordinates[projective_vertex][generator_coordinate] = algebra.field().one();
        let generator = projective.element_from_vertex_coordinates(generator_coordinates)?;
        projective.cache_minimal_generating_set(std::slice::from_ref(&generator))?;
        let _ = projective.cache_is_indecomposable_module(true);
        let _ = projective.cache_is_projective_module(true);
        projectives.push(projective);
    }
    Ok(algebra.merge_indec_projective_modules(&projectives))
}

fn projective_basis_data<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<ProjectiveBasisData<F>>, EnvelopingModuleError> {
    let Some(basis) = finite_canonical_basis(algebra)? else {
        return Ok(None);
    };
    let vertex_count = algebra.quiver().number_of_vertices();
    let mut intervals = vec![vec![Vec::new(); vertex_count]; vertex_count];
    for (basis_index, element) in basis.iter().enumerate() {
        let monomial = element.leading_monomial();
        let PathEndpoint::Vertex(source) = monomial.source() else {
            continue;
        };
        let PathEndpoint::Vertex(target) = monomial.target() else {
            continue;
        };
        intervals[source.index().get() - 1][target.index().get() - 1].push(basis_index);
    }
    Ok(Some(ProjectiveBasisData { basis, intervals }))
}

/// Translate QPA's cached `IsProjectiveModule(M)` dimension test.
///
/// # Errors
///
/// Propagates top and indecomposable-projective construction failures and
/// reports address-size overflow in QPA's dimension sum.
pub fn is_projective_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<bool, EnvelopingModuleError> {
    if let Some(cached) = module.cached_is_projective_module() {
        return Ok(cached);
    }
    let top = top_of_module(module.clone())?;
    let projectives = indec_projective_modules(module.right_acting_algebra())?;
    let mut expected_dimension = 0_usize;
    for (projective, &multiplicity) in projectives.iter().zip(top.dimension_vector()) {
        expected_dimension = expected_dimension
            .checked_add(
                projective
                    .dimension()
                    .checked_mul(multiplicity)
                    .ok_or(EnvelopingModuleError::ProjectiveDimensionOverflow)?,
            )
            .ok_or(EnvelopingModuleError::ProjectiveDimensionOverflow)?;
    }
    Ok(module.cache_is_projective_module(expected_dimension == module.dimension()))
}

/// Translate QPA's `HomFromProjective(m, M)` construction.
///
/// `Ok(None)` is QPA's `fail` when `m` is not supported at exactly one
/// vertex. The source is the cached indecomposable projective at that vertex,
/// matching GAP's identity-preserving `ShallowCopy` on these modules.
///
/// # Errors
///
/// Requires `m` in the identical target module and propagates projective
/// basis, module action, exact matrix, and checked-map failures.
pub fn hom_from_projective<F: Field>(
    element: &PathAlgebraModuleElement<F>,
    target: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, EnvelopingModuleError> {
    if !element.module().is_identical(target) {
        return Err(EnvelopingModuleError::ProjectiveElementInWrongModule);
    }
    let support = element.support_vertices();
    let [vertex] = support.as_slice() else {
        return Ok(None);
    };
    let bases = basis_of_projectives(target.right_acting_algebra())?
        .ok_or(EnvelopingModuleError::ProjectiveBasisUnavailable)?;
    let projectives = indec_projective_modules(target.right_acting_algebra())?;
    let projective = projectives[*vertex - 1].clone();
    let field = target.right_acting_algebra().field().clone();
    let matrices = projective
        .dimension_vector()
        .iter()
        .zip(target.dimension_vector())
        .enumerate()
        .map(
            |(target_vertex, (&source_dimension, &target_dimension))| -> Result<
                DenseMatrix<F>,
                EnvelopingModuleError,
            > {
                if source_dimension == 0 {
                    return Ok(DenseMatrix::zero(
                        1,
                        target_dimension.max(1),
                        field.clone(),
                    ));
                }
                let rows = bases[*vertex - 1][target_vertex]
                    .iter()
                    .map(|basis_element| {
                        Ok(element.act(basis_element)?.external_representation()[target_vertex]
                            .clone())
                    })
                    .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
                Ok(DenseMatrix::from_rows(rows, field.clone())?)
            },
        )
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Some(right_module_hom_over_algebra(
        projective,
        target.clone(),
        matrices,
    )?))
}

/// Translate QPA's cached `ProjectiveCover(M)` attribute.
///
/// The nonzero branch constructs one indecomposable projective map per
/// ordered minimal generator, takes their QPA direct sum, and sums the
/// projection composites in that same order.
///
/// # Errors
///
/// Propagates minimal-generator, projective-map, direct-sum, exact map
/// arithmetic, and checked construction failures.
pub fn projective_cover<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_projective_cover() {
        return Ok(cached);
    }
    let cover = if module.dimension() == 0 {
        let zero = right_module_over_path_algebra(
            module.right_acting_algebra().clone(),
            vec![0; module.dimension_vector().len()],
            &[],
        )?
        .ok_or(EnvelopingModuleError::GeneratedZeroProjectiveViolatesRelations)?;
        zero_module_mapping(zero, module.clone())?
    } else {
        let generators = minimal_generating_set_of_module(module.clone())?;
        let maps = generators
            .iter()
            .map(|generator| {
                hom_from_projective(generator, module)?
                    .ok_or(EnvelopingModuleError::MinimalGeneratorNotVertexUniform)
            })
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
        let sources = maps
            .iter()
            .map(|homomorphism| homomorphism.source().clone())
            .collect::<Vec<_>>();
        let direct_sum = direct_sum_of_qpa_modules(&sources)?
            .ok_or(EnvelopingModuleError::EmptyNonzeroProjectiveCover)?;
        let projections = direct_sum_projections(&direct_sum)?;
        let mut composites = projections
            .iter()
            .zip(&maps)
            .map(|(projection, map)| Ok(projection.compose(map)?))
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?
            .into_iter();
        let mut sum = composites
            .next()
            .ok_or(EnvelopingModuleError::EmptyNonzeroProjectiveCover)?;
        for composite in composites {
            sum = sum.add(&composite)?;
        }
        sum
    };
    Ok(module.cache_projective_cover(&cover))
}

/// Translate QPA's `LiftingMorphismFromProjective(f, g)`.
///
/// `Ok(None)` is QPA's `fail` when the projective source was not constructed
/// with `DirectSumOfQPAModules`; validation errors retain the source method's
/// order.
///
/// # Errors
///
/// Requires a projective source for `g`, a surjective `f`, and a common range,
/// then propagates exact preimage, projective-map, and map arithmetic failures.
pub fn lifting_morphism_from_projective<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, EnvelopingModuleError> {
    let projective = second.source();
    if !is_projective_module(projective)? {
        return Err(EnvelopingModuleError::LiftingSourceNotProjective);
    }
    if !first.is_surjective() {
        return Err(EnvelopingModuleError::LiftingFirstMapNotSurjective);
    }
    if first.range() != second.range() {
        return Err(EnvelopingModuleError::LiftingMapsDoNotHaveCommonRange);
    }
    if projective.dimension() == 0 {
        return Ok(Some(zero_module_mapping(
            projective.clone(),
            first.source().clone(),
        )?));
    }
    if !projective.is_direct_sum_of_modules() {
        return Ok(None);
    }
    let inclusions = qpa_core::direct_sum_inclusions(projective)?;
    let projections = direct_sum_projections(projective)?;
    let mut partial_liftings = Vec::with_capacity(inclusions.len());
    for inclusion in &inclusions {
        let generator = minimal_generating_set_of_module(inclusion.source().clone())?
            .into_iter()
            .next()
            .ok_or(EnvelopingModuleError::EmptyProjectiveSummandGenerators)?;
        let projective_element = inclusion.image(&generator)?;
        let image = second.image(&projective_element)?;
        let generated = if image.is_zero() {
            zero_module_mapping(inclusion.source().clone(), first.source().clone())?
        } else {
            let preimage = first
                .preimage_representative(&image)?
                .ok_or(EnvelopingModuleError::MissingProjectiveLiftPreimage)?;
            let map = hom_from_projective(&preimage, first.source())?
                .ok_or(EnvelopingModuleError::ProjectiveLiftNotVertexUniform)?;
            right_module_hom_over_algebra(
                inclusion.source().clone(),
                first.source().clone(),
                map.vertex_matrices().to_vec(),
            )?
        };
        partial_liftings.push(generated);
    }
    let mut composites = projections
        .iter()
        .zip(&partial_liftings)
        .map(|(projection, partial)| Ok(projection.compose(partial)?))
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?
        .into_iter();
    let mut lifting = composites
        .next()
        .ok_or(EnvelopingModuleError::EmptyProjectiveSummandGenerators)?;
    for composite in composites {
        lifting = lifting.add(&composite)?;
    }
    Ok(Some(lifting))
}

/// Translate QPA's `HomFactoringThroughProjOverAlgebra(M, N)`.
///
/// The returned basis is GAP's insertion-ordered vector-space basis of raw
/// physical map lists after composing `Hom(M, P(N))` with the projective
/// cover of `N` and removing zero rows.
///
/// # Errors
///
/// Requires a common acting algebra and propagates projective-cover, Hom,
/// exact matrix, and checked-map construction failures.
pub fn hom_factoring_through_projectives<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, EnvelopingModuleError> {
    if source.right_acting_algebra() != range.right_acting_algebra() {
        return Err(EnvelopingModuleError::FactoringModulesOverDifferentAlgebras);
    }
    if source.dimension() == 0 || range.dimension() == 0 {
        return Ok(Vec::new());
    }
    let cover = projective_cover(range)?;
    let homomorphisms = hom_over_algebra(source, cover.source())?
        .ok_or(EnvelopingModuleError::InternalProjectiveHomRejected)?;
    let composites = homomorphisms
        .iter()
        .map(|homomorphism| homomorphism.compose(&cover))
        .collect::<Result<Vec<_>, _>>()?;
    let nonzero = composites
        .iter()
        .filter(|homomorphism| !homomorphism.is_zero())
        .collect::<Vec<_>>();
    if nonzero.is_empty() {
        return Ok(Vec::new());
    }
    let field = source.right_acting_algebra().field().clone();
    let row_length = nonzero[0].flat_vertex_matrix_entries().len();
    let rows = DenseMatrix::from_flat(
        nonzero.len(),
        row_length,
        nonzero
            .iter()
            .flat_map(|homomorphism| homomorphism.flat_vertex_matrix_entries())
            .collect(),
        field.clone(),
    )?;
    let basis = rows.semi_echelon_transformation().vectors;
    (0..basis.nrows())
        .map(|row| module_map_from_physical_flat(source, range, basis.row(row), &field))
        .collect()
}

fn module_map_from_physical_flat<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
    entries: &[F::Element],
    field: &F,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let expected = source
        .dimension_vector()
        .iter()
        .zip(range.dimension_vector())
        .try_fold(0_usize, |total, (&source_dimension, &range_dimension)| {
            total.checked_add(
                source_dimension
                    .max(1)
                    .checked_mul(range_dimension.max(1))?,
            )
        })
        .ok_or(EnvelopingModuleError::PhysicalMapSizeOverflow)?;
    if entries.len() != expected {
        return Err(EnvelopingModuleError::PhysicalMapLength {
            actual: entries.len(),
            expected,
        });
    }
    let mut offset = 0_usize;
    let matrices = source
        .dimension_vector()
        .iter()
        .zip(range.dimension_vector())
        .map(|(&source_dimension, &range_dimension)| {
            let rows = source_dimension.max(1);
            let columns = range_dimension.max(1);
            let length = rows
                .checked_mul(columns)
                .ok_or(EnvelopingModuleError::PhysicalMapSizeOverflow)?;
            let matrix = DenseMatrix::from_flat(
                rows,
                columns,
                entries[offset..offset + length].to_vec(),
                field.clone(),
            )?;
            offset += length;
            Ok(matrix)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    Ok(right_module_hom_over_algebra(
        source.clone(),
        range.clone(),
        matrices,
    )?)
}

/// Translate QPA's cached `DualOfModule(M)` for path-algebra matrix modules.
///
/// Arrow matrices are transposed without reordering, the dimension vector is
/// retained, and the result acts over QPA's opposite algebra. For a nonzero
/// module, QPA explicitly installs the reciprocal dual attribute on the
/// result; its zero-module early return does not.
///
/// # Errors
///
/// Propagates opposite-algebra construction and the checked right-module
/// constructor. A generated action failing the opposite relations is reported
/// as an internal source-translation failure.
pub fn dual_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, ModuleFunctorError> {
    if let Some(cached) = module.cached_dual_module() {
        return Ok(cached);
    }
    let opposite = match module.right_acting_algebra() {
        QuiverAlgebra::Path(algebra) => QuiverAlgebra::Path(algebra.opposite()?),
        QuiverAlgebra::Quotient(algebra) => {
            QuiverAlgebra::Quotient(opposite_path_algebra(algebra)?)
        }
    };
    let dimensions = module.dimension_vector().to_vec();
    let mut specification = Vec::new();
    if module.dimension() != 0 {
        for (position, arrow) in opposite.quiver().arrows().enumerate() {
            let source = dimensions[arrow.source().index().get() - 1];
            let target = dimensions[arrow.target().index().get() - 1];
            if source != 0 && target != 0 {
                specification.push(NamedArrowMatrix::new(
                    arrow.name(),
                    module.arrow_matrices()[position].transpose(),
                ));
            }
        }
    }
    let dual = right_module_over_path_algebra(opposite, dimensions, &specification)?
        .ok_or(ModuleFunctorError::GeneratedDualViolatesRelations)?;
    if module.dimension() == 0 {
        Ok(module.cache_dual_module(&dual))
    } else {
        Ok(module.install_dual_module(&dual)?)
    }
}

/// Translate QPA's cached `DualOfModuleHomomorphism(f)` attribute.
///
/// The source transposes every vertex map in place, dualizes the source and
/// range modules, reverses them, and passes the matrices through the checked
/// module-homomorphism constructor. QPA caches only the entered map's
/// attribute; it does not install a reciprocal cache on the returned map.
///
/// # Errors
///
/// Propagates module dualization and checked map construction.
pub fn dual_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, ModuleFunctorError> {
    if let Some(cached) = homomorphism.cached_dual_homomorphism() {
        return Ok(cached);
    }
    let matrices = homomorphism
        .vertex_matrices()
        .iter()
        .map(DenseMatrix::transpose)
        .collect();
    let dual_source = dual_of_module(homomorphism.source())?;
    let dual_range = dual_of_module(homomorphism.range())?;
    let dual = right_module_hom_over_algebra(dual_range, dual_source, matrices)?;
    Ok(homomorphism.cache_dual_homomorphism(&dual))
}

/// Translate QPA's cached
/// `DualOfAlgebraAsModuleOverEnvelopingAlgebra(A)` attribute.
///
/// QPA first dualizes the regular bimodule over the opposite enveloping
/// algebra. It then identifies that opposite with `A^op tensor A` by swapping
/// the two product-quiver coordinates, permutes vertex spaces and arrow
/// matrices through the generated names, and reconstructs a checked module
/// over the original enveloping algebra.
///
/// `Ok(None)` is the cached GAP `fail` value for an infinite-dimensional
/// algebra or an unexpected checked-relation failure.
///
/// # Errors
///
/// Propagates finite-dimensionality, enveloping/dual construction, product
/// projection and name-map failures, and checked module construction.
pub fn dual_of_algebra_as_module_over_enveloping_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_dual_enveloping_module() {
        return Ok(Some(cached));
    }
    if algebra.dual_enveloping_module_attribute_is_fail() {
        return Ok(None);
    }
    if !is_finite_dimensional(algebra)? {
        algebra.mark_dual_enveloping_module_attribute_fail();
        return Ok(None);
    }

    let enveloping = enveloping_algebra(algebra)?;
    let Some(regular) = algebra_as_module_over_enveloping_algebra(algebra)? else {
        algebra.mark_dual_enveloping_module_attribute_fail();
        return Ok(None);
    };
    let dual = dual_of_module(&regular)?;
    let product = enveloping.product_quiver();
    let vertices = enveloping.algebra().quiver().vertices().collect::<Vec<_>>();
    let arrows = enveloping.algebra().quiver().arrows().collect::<Vec<_>>();
    let vertex_names = vertices
        .iter()
        .map(qpa_core::Vertex::name)
        .collect::<Vec<_>>();
    let arrow_names = arrows.iter().map(qpa_core::Arrow::name).collect::<Vec<_>>();

    let vertex_positions = vertices
        .iter()
        .map(|vertex| {
            let name = swapped_product_generator_name(&vertex.path(), product, algebra)?;
            vertex_names
                .iter()
                .position(|candidate| *candidate == name.as_ref())
                .ok_or(EnvelopingModuleError::MissingSwappedProductGenerator { name })
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let arrow_positions = arrows
        .iter()
        .map(|arrow| {
            let name = swapped_product_generator_name(&arrow.path(), product, algebra)?;
            arrow_names
                .iter()
                .position(|candidate| *candidate == name.as_ref())
                .ok_or(EnvelopingModuleError::MissingSwappedProductGenerator { name })
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let dimensions = vertex_positions
        .iter()
        .map(|&position| dual.dimension_vector()[position])
        .collect::<Vec<_>>();
    let matrices = arrow_positions
        .iter()
        .map(|&position| dual.arrow_matrices()[position].clone())
        .collect::<Vec<_>>();
    let specification = arrows
        .iter()
        .enumerate()
        .filter(|(_, arrow)| {
            dimensions[arrow.source().index().get() - 1] != 0
                && dimensions[arrow.target().index().get() - 1] != 0
        })
        .map(|(position, arrow)| NamedArrowMatrix::new(arrow.name(), matrices[position].clone()))
        .collect::<Vec<_>>();

    let Some(module) = right_module_over_path_algebra_with_tensor(
        enveloping.algebra().clone(),
        dimensions,
        &specification,
        enveloping,
    )?
    else {
        algebra.mark_dual_enveloping_module_attribute_fail();
        return Ok(None);
    };
    Ok(Some(algebra.cache_dual_enveloping_module(&module)))
}

fn is_finite_dimensional<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<bool, EnvelopingModuleError> {
    match algebra {
        QuiverAlgebra::Path(algebra) => Ok(algebra.has_finite_canonical_basis()),
        QuiverAlgebra::Quotient(algebra) => Ok(quotient_is_finite_dimensional(algebra)?),
    }
}

fn swapped_product_generator_name<F: Field>(
    generator: &qpa_core::Path,
    product: &QuiverProduct,
    algebra: &QuiverAlgebra<F>,
) -> Result<Box<str>, EnvelopingModuleError> {
    let first = project_from_product_quiver(1, generator, product)?;
    let second = project_from_product_quiver(2, generator, product)?;
    let first_name = first.quiver().opposite_quiver_name(&first.to_string())?;
    let second_name = algebra.quiver().opposite_quiver_name(&second.to_string())?;
    Ok(format!("{second_name}_{first_name}").into_boxed_str())
}

pub(crate) fn finite_canonical_basis<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<Vec<QuiverAlgebraElement<F>>>, EnvelopingModuleError> {
    match algebra {
        QuiverAlgebra::Path(path_algebra) => Ok(path_algebra
            .canonical_basis()
            .map(|basis| basis.into_iter().map(QuiverAlgebraElement::Path).collect())),
        QuiverAlgebra::Quotient(quotient) => {
            if !quotient_is_finite_dimensional(quotient)? {
                return Ok(None);
            }
            ensure_quotient_groebner_basis(quotient)?;
            Ok(Some(
                quotient
                    .canonical_basis()?
                    .ok_or(EnvelopingModuleError::MissingFiniteCanonicalBasis)?
                    .into_iter()
                    .map(QuiverAlgebraElement::Quotient)
                    .collect(),
            ))
        }
    }
}

pub(crate) fn canonical_basis_coefficients<F: Field>(
    algebra: &QuiverAlgebra<F>,
    basis_vectors: &[QuiverAlgebraElement<F>],
    element: &QuiverAlgebraElement<F>,
) -> Result<Vec<F::Element>, EnvelopingModuleError> {
    match (algebra, element) {
        (QuiverAlgebra::Path(path_algebra), QuiverAlgebraElement::Path(element)) => {
            let mut coordinates = vec![path_algebra.field().zero(); basis_vectors.len()];
            for term in element.terms() {
                let position = basis_vectors
                    .iter()
                    .position(|basis_vector| {
                        matches!(basis_vector, QuiverAlgebraElement::Path(basis)
                            if basis.terms().len() == 1
                                && basis.terms()[0].path() == term.path())
                    })
                    .ok_or(EnvelopingModuleError::ElementOutsideCanonicalBasis)?;
                coordinates[position] = term.coefficient().clone();
            }
            Ok(coordinates)
        }
        (QuiverAlgebra::Quotient(quotient), QuiverAlgebraElement::Quotient(element)) => {
            Ok(quotient.canonical_basis_coefficients(element)?)
        }
        _ => Err(EnvelopingModuleError::UnexpectedElementVariant),
    }
}

/// Failure while constructing QPA's regular bimodule over its enveloping
/// algebra.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum EnvelopingModuleError {
    /// A finite quotient must expose the canonical basis prepared by GBNP.
    #[error("a finite-dimensional quotient has no canonical basis vectors")]
    MissingFiniteCanonicalBasis,
    /// A path-algebra image contained a path absent from its canonical basis.
    #[error("an element path is absent from the finite canonical basis")]
    ElementOutsideCanonicalBasis,
    /// The algebra and element variants must remain paired.
    #[error("canonical-basis arithmetic returned another quiver-algebra variant")]
    UnexpectedElementVariant,
    /// The coordinate-swap name must resolve in the original product quiver.
    #[error("cannot find swapped product-quiver generator {name:?}")]
    MissingSwappedProductGenerator {
        /// Generated QPA product name.
        name: Box<str>,
    },
    /// QPA defines algebra Loewy length only for finite-dimensional input.
    #[error("the entered algebra is not finite dimensional")]
    InfiniteLoewyAlgebra,
    /// QPA's method requires an admissible quotient when input is not a path algebra.
    #[error("LoewyLength requires an admissible quotient of a path algebra")]
    NonadmissibleLoewyAlgebra,
    /// A projective action generated from canonical basis multiplication must satisfy relations.
    #[error("a generated indecomposable projective does not satisfy the algebra relations")]
    GeneratedProjectiveViolatesRelations,
    /// Every projective right ideal contains its generating vertex idempotent.
    #[error("the canonical basis of projective {vertex} has no vertex generator")]
    MissingProjectiveVertexGenerator {
        /// One-based projective/vertex position.
        vertex: usize,
    },
    /// QPA's projectivity dimension sum exceeded Rust's addressable size.
    #[error("the projective dimension test overflows usize")]
    ProjectiveDimensionOverflow,
    /// `HomFromProjective` requires an element of the identical target module.
    #[error("the projective generator element is not in the target module")]
    ProjectiveElementInWrongModule,
    /// A projective basis cannot be extracted from infinite-dimensional input.
    #[error("BasisOfProjectives returned fail")]
    ProjectiveBasisUnavailable,
    /// A generated all-zero module must satisfy the acting relations.
    #[error("the zero projective module does not satisfy the algebra relations")]
    GeneratedZeroProjectiveViolatesRelations,
    /// Every nonzero module has at least one top/minimal generator.
    #[error("a nonzero projective cover has no summands")]
    EmptyNonzeroProjectiveCover,
    /// QPA's minimal generators are vertex-uniform.
    #[error("a minimal generator is not supported at exactly one vertex")]
    MinimalGeneratorNotVertexUniform,
    /// Projective lifting validates the source before its other inputs.
    #[error("the source of the second map is not projective")]
    LiftingSourceNotProjective,
    /// Projective lifting requires an epimorphism as its first map.
    #[error("the first map is not surjective")]
    LiftingFirstMapNotSurjective,
    /// Both maps in a projective lifting diagram must end at the same module.
    #[error("the entered maps do not have a common range")]
    LiftingMapsDoNotHaveCommonRange,
    /// Every nonzero direct projective summand has its installed generator.
    #[error("a nonzero projective summand has no minimal generator")]
    EmptyProjectiveSummandGenerators,
    /// Surjectivity guarantees a lift of every projective-generator image.
    #[error("a projective-generator image unexpectedly has no preimage")]
    MissingProjectiveLiftPreimage,
    /// A projective lift preimage must be vertex-uniform.
    #[error("a projective lifting preimage is not vertex-uniform")]
    ProjectiveLiftNotVertexUniform,
    /// Factoring-through computations require a common acting algebra.
    #[error("the factoring modules are not modules over the same algebra")]
    FactoringModulesOverDifferentAlgebras,
    /// A Hom call cannot reject the checked projective-cover algebra.
    #[error("HomOverAlgebra rejected the projective cover's common algebra")]
    InternalProjectiveHomRejected,
    /// Flattened physical map data exceeded Rust's addressable size.
    #[error("the physical module-map size overflows usize")]
    PhysicalMapSizeOverflow,
    /// A raw physical map vector must fill every vertex block exactly.
    #[error("the physical module-map vector has length {actual}, expected {expected}")]
    PhysicalMapLength {
        /// Entered flat-vector length.
        actual: usize,
        /// Required total physical block length.
        expected: usize,
    },
    /// The configured QPA/GBNP finite-basis path failed.
    #[error(transparent)]
    Backend(#[from] PathalgGbnpError),
    /// Enveloping-algebra construction failed.
    #[error(transparent)]
    Tensor(#[from] TensorProductError),
    /// Product-quiver projection failed.
    #[error(transparent)]
    Product(#[from] QuiverProductError),
    /// Opposite-path reconstruction failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// Path inclusion in the source algebra failed.
    #[error(transparent)]
    PathTensor(#[from] PathAlgebraTensorError),
    /// Source-algebra multiplication failed.
    #[error(transparent)]
    Element(#[from] QuiverAlgebraElementError),
    /// Quotient canonical-basis coordinates failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Exact matrix construction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// The checked QPA right-module constructor rejected its generated data.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// Generic QPA module dualization failed.
    #[error(transparent)]
    Functor(#[from] ModuleFunctorError),
}

/// Failure in QPA's generic module-dual functor.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ModuleFunctorError {
    /// The generated transposed actions should satisfy the opposite relations.
    #[error("the generated dual module does not satisfy the opposite-algebra relations")]
    GeneratedDualViolatesRelations,
    /// Parent path-algebra opposite construction failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Quotient opposite construction through QPA's GBNP backend failed.
    #[error(transparent)]
    Backend(#[from] PathalgGbnpError),
    /// The checked module constructor or reciprocal cache failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring,
        finite_field::{FiniteFieldCore, Z2},
        integer::Integer,
        rational::{Q, Rational},
    };
    use qpa_core::path_module::right_fac_approximation_via_all_basis_images;
    use qpa_core::{
        ArrowSpec, DenseMatrix, NamedArrowMatrix, NamedArrowSpec, PathAlgebra, PathAlgebraModule,
        PathAlgebraModuleElement, PathAlgebraModuleHomomorphism, Quiver, QuiverAlgebra,
        VertexIndex, direct_sum_of_qpa_modules, lifting_inclusion_morphisms,
        minimal_generating_set_of_module, right_fac_approximation, right_module_hom_over_algebra,
        right_module_over_path_algebra, socle_of_module, socle_of_module_inclusion, top_of_module,
        top_of_module_projection, trace_of_module, zero_module_mapping,
    };
    use serde_json::{Value, json};

    use super::{
        algebra_as_module_of_enveloping_algebra, algebra_as_module_over_enveloping_algebra,
        basis_of_projectives, dual_of_algebra_as_module_over_enveloping_algebra, dual_of_module,
        dual_of_module_homomorphism, hom_factoring_through_projectives, hom_from_projective,
        indec_projective_modules, is_projective_module, lifting_morphism_from_projective,
        projective_cover,
    };
    use crate::{enveloping_algebra, quotient_by_relations};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn q(value: i64) -> Rational {
        Q.nth(Integer::from(value))
    }

    #[test]
    fn regular_bimodule_for_dual_numbers_has_qpa_row_actions() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(Q, quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );

        let module = algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        let again = algebra_as_module_of_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        let enveloping = enveloping_algebra(&algebra).unwrap();

        assert_eq!(module.dimension_vector(), &[2]);
        assert_eq!(module.arrow_matrices().len(), 2);
        assert_eq!(
            module.arrow_matrices()[0].as_flat(),
            &[q(0), q(1), q(0), q(0)]
        );
        assert_eq!(
            module.arrow_matrices()[1].as_flat(),
            &[q(0), q(1), q(0), q(0)]
        );
        assert!(module.is_identical(&again));
        assert!(
            module
                .retained_tensor_product()
                .unwrap()
                .is_identical(&enveloping)
        );

        let generic_dual = dual_of_module(&module).unwrap();
        let generic_back = dual_of_module(&generic_dual).unwrap();
        let bimodule_dual = dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        assert_eq!(
            generic_dual.arrow_matrices()[0].as_flat(),
            &[q(0), q(0), q(1), q(0)]
        );
        assert!(generic_back.is_identical(&module));
        assert_eq!(
            bimodule_dual.arrow_matrices()[0].as_flat(),
            &[q(0), q(0), q(1), q(0)]
        );
    }

    #[test]
    fn dual_of_module_homomorphism_transposes_reverses_and_caches_like_qpa() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra,
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![q(1)], vec![q(0)]], Q).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let map = right_module_hom_over_algebra(
            module.clone(),
            module,
            vec![
                DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
                DenseMatrix::identity(1, Q),
            ],
        )
        .unwrap();

        let dual = dual_of_module_homomorphism(&map).unwrap();
        assert_eq!(dual.source().dimension_vector(), &[2, 1]);
        assert_eq!(dual.range().dimension_vector(), &[2, 1]);
        assert_eq!(
            dual.vertex_matrices(),
            &[
                DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
                DenseMatrix::identity(1, Q),
            ]
        );
        assert!(dual.is_identical(&dual_of_module_homomorphism(&map).unwrap()));
        let double_dual = dual_of_module_homomorphism(&dual).unwrap();
        assert_eq!(double_dual, map);
        assert!(
            !double_dual.is_identical(&map),
            "the QPA source does not install a reciprocal map attribute"
        );

        let e = |value: i64| Z2.nth(Integer::from(value));
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra,
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![e(1)], vec![e(0)]], Z2).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let map = right_module_hom_over_algebra(
            module.clone(),
            module,
            vec![
                DenseMatrix::from_rows(vec![vec![e(1), e(0)], vec![e(0), e(0)]], Z2).unwrap(),
                DenseMatrix::identity(1, Z2),
            ],
        )
        .unwrap();
        assert_eq!(
            dual_of_module_homomorphism(&map).unwrap().vertex_matrices(),
            map.vertex_matrices()
        );
    }

    #[test]
    fn regular_bimodule_supports_characteristic_two_and_zero_spaces() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let module = algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();

        assert_eq!(module.dimension_vector(), &[1, 1, 0, 1]);
        assert_eq!(module.dimension(), 3);
        assert_eq!(module.arrow_matrices().len(), 4);
        assert!(module.arrow_matrices().iter().all(|matrix| {
            matrix
                .as_flat()
                .iter()
                .all(|entry| Z2.to_integer(entry) == 0 || Z2.to_integer(entry) == 1)
        }));

        let dual = dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        assert_eq!(dual.dimension_vector(), &[1, 0, 1, 1]);
    }

    #[test]
    fn cyclic_path_algebra_returns_qpa_fail_before_enveloping_construction() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());

        assert!(
            algebra_as_module_over_enveloping_algebra(&algebra)
                .unwrap()
                .is_none()
        );
        assert!(algebra.cached_enveloping_algebra().is_none());
        assert!(
            dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
                .unwrap()
                .is_none()
        );
        assert!(algebra.has_dual_enveloping_module_attribute());
    }

    fn matrix_json<F: Field>(
        matrix: &DenseMatrix<F>,
        coefficient: impl Fn(&F::Element) -> String,
    ) -> Value {
        json!({
            "shape": [matrix.nrows(), matrix.ncols()],
            "rows": (0..matrix.nrows())
                .map(|row| matrix.row(row).iter().map(&coefficient).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
        })
    }

    fn module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "dimension_vector": module.dimension_vector(),
            "dimension": module.dimension(),
            "arrows": module
                .right_acting_algebra()
                .quiver()
                .arrows()
                .map(|arrow| arrow.name().to_owned())
                .collect::<Vec<_>>(),
            "matrices": module
                .arrow_matrices()
                .iter()
                .map(|matrix| matrix_json(matrix, coefficient))
                .collect::<Vec<_>>(),
        })
    }

    fn module_map_json<F: Field>(
        map: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "source_dimension": map.source().dimension_vector(),
            "range_dimension": map.range().dimension_vector(),
            "matrices": map
                .vertex_matrices()
                .iter()
                .map(|matrix| matrix_json(matrix, coefficient))
                .collect::<Vec<_>>(),
        })
    }

    fn projective_module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "dimension": module.dimension_vector(),
            "matrices": module
                .arrow_matrices()
                .iter()
                .map(|matrix| matrix_json(matrix, coefficient))
                .collect::<Vec<_>>(),
        })
    }

    fn module_elements_json<F: Field>(
        elements: &[PathAlgebraModuleElement<F>],
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!(
            elements
                .iter()
                .map(|element| element
                    .external_representation()
                    .iter()
                    .map(|row| row.iter().map(coefficient).collect::<Vec<_>>())
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>()
        )
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_functor_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one.clone()], vec![zero.clone()]], field.clone())
                    .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let map = right_module_hom_over_algebra(
            module.clone(),
            module.clone(),
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![one.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::identity(1, field.clone()),
            ],
        )
        .unwrap();
        let minimal = minimal_generating_set_of_module(module.clone()).unwrap();
        let top_projection = top_of_module_projection(module.clone()).unwrap();
        let top = top_of_module(module.clone()).unwrap();
        let socle_inclusion = socle_of_module_inclusion(module.clone()).unwrap();
        let dual_map = dual_of_module_homomorphism(&map).unwrap();
        let target = right_module_over_path_algebra(
            algebra.clone(),
            vec![1, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one.clone(), one]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let zero_module = right_module_over_path_algebra(algebra, vec![0, 0], &[])
            .unwrap()
            .unwrap();
        let zero_socle = socle_of_module_inclusion(zero_module.clone()).unwrap();
        let zero_trace = trace_of_module(&simple_one, &target).unwrap();
        let zero_source_trace = trace_of_module(&zero_module, &target).unwrap();
        let zero_target_trace = trace_of_module(&simple_two, &zero_module).unwrap();
        let nonzero_trace = trace_of_module(&simple_two, &target).unwrap();
        let legacy_trace =
            right_fac_approximation_via_all_basis_images(&simple_two, &target).unwrap();
        let ordered_trace = trace_of_module(&module, &module).unwrap();
        let ordered_legacy =
            right_fac_approximation_via_all_basis_images(&module, &module).unwrap();

        json!({
            "minimal_generators": {
                "count": minimal.len(),
                "elements": module_elements_json(&minimal, coefficient),
                "zero_count": minimal_generating_set_of_module(zero_module.clone()).unwrap().len(),
                "top_projection_cached_identity": top_projection.is_identical(
                    &top_of_module_projection(module.clone()).unwrap()),
                "top_cached_identity": top.is_identical(&top_of_module(module).unwrap()),
                "top_projection_range_identity": top_projection.range().is_identical(&top),
            },
            "socle": {
                "inclusion": module_map_json(&socle_inclusion, coefficient),
                "object_dimension": socle_of_module(socle_inclusion.range().clone())
                    .unwrap()
                    .dimension_vector(),
                "zero_inclusion": module_map_json(&zero_socle, coefficient),
                "zero_source_identity": zero_socle.source().is_identical(&zero_module),
            },
            "dual_map": {
                "cached_identity": dual_map.is_identical(
                    &dual_of_module_homomorphism(&map).unwrap()),
                "source_dimension": dual_map.source().dimension_vector(),
                "range_dimension": dual_map.range().dimension_vector(),
                "matrices": dual_map.vertex_matrices().iter()
                    .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
            },
            "trace": {
                "zero": module_map_json(&zero_trace, coefficient),
                "zero_source": module_map_json(&zero_source_trace, coefficient),
                "zero_target": module_map_json(&zero_target_trace, coefficient),
                "nonzero": module_map_json(&nonzero_trace, coefficient),
                "legacy": module_map_json(&legacy_trace, coefficient),
                "legacy_equal": legacy_trace == nonzero_trace,
                "ordered": module_map_json(&ordered_trace, coefficient),
                "ordered_legacy": module_map_json(&ordered_legacy, coefficient),
                "ordered_legacy_equal": ordered_legacy == ordered_trace,
                "right_fac_equal": right_fac_approximation(&simple_two, &target).unwrap()
                    == nonzero_trace,
            },
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_projective_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let projectives = indec_projective_modules(&algebra).unwrap();
        let bases = basis_of_projectives(&algebra).unwrap().unwrap();
        let cached_projectives = indec_projective_modules(&algebra).unwrap();
        let basis_paths = bases
            .iter()
            .map(|projective| {
                projective
                    .iter()
                    .map(|vertex_basis| {
                        vertex_basis
                            .iter()
                            .map(|element| element.leading_monomial().to_string())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(
                    vec![
                        vec![one.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let generators = minimal_generating_set_of_module(module.clone()).unwrap();
        let cover = projective_cover(&module).unwrap();
        let from_generator = hom_from_projective(&generators[0], &module)
            .unwrap()
            .unwrap();
        let projective_sum = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
        let lifting = lifting_morphism_from_projective(&cover, &cover)
            .unwrap()
            .unwrap();
        let nondirect_lifting = lifting_morphism_from_projective(&cover, &from_generator).unwrap();
        let zero_module = right_module_over_path_algebra(algebra.clone(), vec![0, 0], &[])
            .unwrap()
            .unwrap();
        let zero_to_module = zero_module_mapping(zero_module.clone(), module.clone()).unwrap();
        let zero_lifting = lifting_morphism_from_projective(&cover, &zero_to_module)
            .unwrap()
            .unwrap();
        let factoring = hom_factoring_through_projectives(&module, &module).unwrap();
        let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let layer = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
            .unwrap()
            .unwrap();
        let first = right_module_hom_over_algebra(
            layer.clone(),
            module.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone())
                    .unwrap(),
                DenseMatrix::from_rows(vec![vec![one.clone(), zero.clone()]], field.clone())
                    .unwrap(),
            ],
        )
        .unwrap();
        let second = right_module_hom_over_algebra(
            layer,
            module.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone())
                    .unwrap(),
                DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone())
                    .unwrap(),
            ],
        )
        .unwrap();
        let kernel_module = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let kernel_inclusion = right_module_hom_over_algebra(
            kernel_module,
            module.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone())
                    .unwrap(),
                DenseMatrix::zero(1, 2, field.clone()),
            ],
        )
        .unwrap();
        let inclusion_lift = lifting_inclusion_morphisms(&first, &kernel_inclusion)
            .unwrap()
            .unwrap();
        let inclusion_fail = lifting_inclusion_morphisms(&first, &second).unwrap();

        let loop_quiver =
            Quiver::from_named(vec!["w".into()], vec![NamedArrowSpec::named("w", "w", "x")])
                .unwrap();
        let path_algebra = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let loop_element = path_algebra
            .element(path_algebra.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(
                path_algebra,
                &[loop_element.multiply(&loop_element).unwrap()],
            )
            .unwrap(),
        );
        let dual_projective = indec_projective_modules(&dual_numbers).unwrap()[0].clone();
        let dual_basis = basis_of_projectives(&dual_numbers).unwrap().unwrap();
        let dual_basis_paths = dual_basis
            .iter()
            .map(|projective| {
                projective
                    .iter()
                    .map(|vertex_basis| {
                        vertex_basis
                            .iter()
                            .map(|element| element.leading_monomial().to_string())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let dual_simple = right_module_over_path_algebra(dual_numbers, vec![1], &[])
            .unwrap()
            .unwrap();
        let dual_cover = projective_cover(&dual_simple).unwrap();

        json!({
            "a2": {
                "basis_paths": basis_paths,
                "projectives": projectives
                    .iter()
                    .map(|projective| projective_module_json(projective, coefficient))
                    .collect::<Vec<_>>(),
                "cache_identity": projectives
                    .iter()
                    .zip(&cached_projectives)
                    .map(|(projective, cached)| projective.is_identical(cached))
                    .collect::<Vec<_>>(),
                "installed_properties": [
                    projectives[0].cached_is_projective_module() == Some(true),
                    projectives[1].cached_is_projective_module() == Some(true),
                    projectives[0].cached_is_indecomposable_module() == Some(true),
                    projectives[1].cached_is_indecomposable_module() == Some(true),
                ],
                "projectivity": [
                    is_projective_module(&projective_sum).unwrap(),
                    is_projective_module(&simple_one).unwrap(),
                    is_projective_module(&simple_two).unwrap(),
                    is_projective_module(&zero_module).unwrap(),
                ],
                "minimal_generator_count": generators.len(),
                "from_generator": module_map_json(&from_generator, coefficient),
                "cover": module_map_json(&cover, coefficient),
                "cover_cache_identity": cover.is_identical(&projective_cover(&module).unwrap()),
                "lifting": module_map_json(&lifting, coefficient),
                "lifting_triangle": lifting.compose(&cover).unwrap() == cover,
                "nondirect_lifting_is_fail": nondirect_lifting.is_none(),
                "zero_lifting": module_map_json(&zero_lifting, coefficient),
                "factoring_count": factoring.len(),
                "factoring": factoring
                    .iter()
                    .map(|map| module_map_json(map, coefficient))
                    .collect::<Vec<_>>(),
                "factoring_zero_counts": [
                    hom_factoring_through_projectives(&zero_module, &module).unwrap().len(),
                    hom_factoring_through_projectives(&module, &zero_module).unwrap().len(),
                ],
                "inclusion_lift": module_map_json(&inclusion_lift, coefficient),
                "inclusion_triangle": inclusion_lift.compose(&first).unwrap() == kernel_inclusion,
                "inclusion_fail": inclusion_fail.is_none(),
            },
            "dual_numbers": {
                "basis_paths": dual_basis_paths,
                "projective": projective_module_json(&dual_projective, coefficient),
                "generator_count": minimal_generating_set_of_module(dual_projective.clone())
                    .unwrap().len(),
                "is_projective": is_projective_module(&dual_projective).unwrap(),
                "simple_is_projective": is_projective_module(&dual_simple).unwrap(),
                "cover": module_map_json(&dual_cover, coefficient),
            },
        })
    }

    #[test]
    fn pinned_qpa_modulehom_projective_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_projectives.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_projective_case(Q, ToString::to_string);
        let binary = modulehom_projective_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_modulehom_functor_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_functors.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_functor_case(Q, ToString::to_string);
        let binary = modulehom_functor_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_regular_enveloping_module_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalgtensor_regular_module.schema-v1.json"
        ))
        .unwrap();

        let loop_quiver =
            Quiver::from_named(vec!["v".into()], vec![NamedArrowSpec::named("v", "v", "x")])
                .unwrap();
        let parent = PathAlgebra::new(Q, loop_quiver.clone()).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let rational_algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let rational_cached_before = rational_algebra.has_regular_enveloping_module_attribute();
        let rational_module = algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let rational_again = algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let rational_synonym = algebra_as_module_of_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let rational_enveloping = enveloping_algebra(&rational_algebra).unwrap();
        let rational_json = json!({
            "attribute_identity": [
                rational_cached_before,
                rational_algebra.has_regular_enveloping_module_attribute(),
                rational_module.is_identical(&rational_again),
                rational_module.is_identical(&rational_synonym),
                rational_module.right_acting_algebra() == rational_enveloping.algebra(),
            ],
            "module": module_json(&rational_module, ToString::to_string),
        });

        let path_quiver = Quiver::from_named(
            vec!["u1".into(), "u2".into()],
            vec![NamedArrowSpec::named("u1", "u2", "a")],
        )
        .unwrap();
        let binary_algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, path_quiver).unwrap());
        let binary_module = algebra_as_module_over_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_again = algebra_as_module_of_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_enveloping = enveloping_algebra(&binary_algebra).unwrap();
        let binary_json = json!({
            "attribute_identity": [
                binary_module.is_identical(&binary_again),
                binary_module.right_acting_algebra() == binary_enveloping.algebra(),
            ],
            "module": module_json(&binary_module, |entry| {
                if Z2.is_zero(entry) {
                    "0*Z(2)".to_owned()
                } else {
                    "Z(2)^0".to_owned()
                }
            }),
        });

        let infinite_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, loop_quiver).unwrap());
        let infinite_result = algebra_as_module_over_enveloping_algebra(&infinite_algebra).unwrap();
        let infinite_json = json!({
            "is_fail": infinite_result.is_none(),
            "attribute_bound": infinite_algebra.has_regular_enveloping_module_attribute(),
            "enveloping_bound": infinite_algebra.cached_enveloping_algebra().is_some(),
        });

        assert_eq!(rational_json, expected["rational_quotient"]);
        assert_eq!(binary_json, expected["binary_path"]);
        assert_eq!(infinite_json, expected["infinite_path"]);
    }

    #[test]
    fn pinned_qpa_dual_enveloping_module_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalgtensor_dual_module.schema-v1.json"
        ))
        .unwrap();

        let loop_quiver =
            Quiver::from_named(vec!["v".into()], vec![NamedArrowSpec::named("v", "v", "x")])
                .unwrap();
        let parent = PathAlgebra::new(Q, loop_quiver.clone()).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let rational_algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let cached_before = rational_algebra.has_dual_enveloping_module_attribute();
        let bimodule_dual = dual_of_algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let bimodule_again = dual_of_algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let regular = algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let generic_dual = dual_of_module(&regular).unwrap();
        let generic_again = dual_of_module(&regular).unwrap();
        let generic_back = dual_of_module(&generic_dual).unwrap();
        let rational_enveloping = enveloping_algebra(&rational_algebra).unwrap();
        let rational_json = json!({
            "attribute_identity": [
                cached_before,
                rational_algebra.has_dual_enveloping_module_attribute(),
                bimodule_dual.is_identical(&bimodule_again),
                bimodule_dual.right_acting_algebra() == rational_enveloping.algebra(),
            ],
            "generic_identity": [
                generic_dual.is_identical(&generic_again),
                generic_back.is_identical(&regular),
            ],
            "generic_dual": module_json(&generic_dual, ToString::to_string),
            "bimodule_dual": module_json(&bimodule_dual, ToString::to_string),
        });

        let path_quiver = Quiver::from_named(
            vec!["u1".into(), "u2".into()],
            vec![NamedArrowSpec::named("u1", "u2", "a")],
        )
        .unwrap();
        let binary_algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, path_quiver).unwrap());
        let binary_dual = dual_of_algebra_as_module_over_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_again = dual_of_algebra_as_module_over_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_enveloping = enveloping_algebra(&binary_algebra).unwrap();
        let binary_json = json!({
            "attribute_identity": [
                binary_dual.is_identical(&binary_again),
                binary_dual.right_acting_algebra() == binary_enveloping.algebra(),
            ],
            "bimodule_dual": module_json(&binary_dual, |entry| {
                if Z2.is_zero(entry) {
                    "0*Z(2)".to_owned()
                } else {
                    "Z(2)^0".to_owned()
                }
            }),
        });

        let infinite_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, loop_quiver).unwrap());
        let infinite_result =
            dual_of_algebra_as_module_over_enveloping_algebra(&infinite_algebra).unwrap();
        let infinite_json = json!({
            "is_fail": infinite_result.is_none(),
            "attribute_bound": infinite_algebra.has_dual_enveloping_module_attribute(),
            "regular_bound": infinite_algebra.has_regular_enveloping_module_attribute(),
            "enveloping_bound": infinite_algebra.cached_enveloping_algebra().is_some(),
        });

        assert_eq!(rational_json, expected["rational_quotient"]);
        assert_eq!(binary_json, expected["binary_path"]);
        assert_eq!(infinite_json, expected["infinite_path"]);
    }
}
