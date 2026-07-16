//! QPA modules constructed from enveloping algebras.
//!
//! This begins the source-order translation of `lib/pathalgtensor.gi`'s
//! module section. The regular bimodule method deliberately uses QPA's
//! canonical-basis partition and its checked path-module constructor.

use numerica::domains::Field;
use qpa_core::{
    DenseMatrix, MatrixError, ModuleDimensionSeries, ModuleEndomorphismQuotientHomomorphism,
    NamedArrowMatrix, PathAlgebraElement, PathAlgebraError, PathAlgebraModule,
    PathAlgebraModuleElement, PathAlgebraModuleError, PathAlgebraModuleHomomorphism,
    PathAlgebraTensorError, PathEndpoint, PathError, QuiverAlgebra, QuiverAlgebraElement,
    QuiverAlgebraElementError, QuiverProduct, QuiverProductError, QuotientError,
    cokernel_of_module_homomorphism, direct_sum_inclusions, direct_sum_of_qpa_modules,
    direct_sum_projections, end_over_algebra, from_hom_mm_to_end_m, hom_over_algebra,
    identity_module_mapping, include_in_path_algebra, kernel_inclusion,
    kernel_of_module_homomorphism, loewy_length_of_module, minimal_generating_set_of_module,
    morphism_on_cokernel, morphism_on_kernel, natural_endomorphism_quotient_by_ideal,
    project_from_product_quiver, radical_of_module, right_module_hom_over_algebra,
    right_module_over_path_algebra, right_module_over_path_algebra_with_tensor,
    subrepresentation_inclusion, top_of_module, zero_module, zero_module_mapping,
};
use thiserror::Error;

use crate::{
    PathalgGbnpError, TensorProductError, ensure_quotient_groebner_basis, enveloping_algebra,
    opposite_path_algebra, opposite_quotient_element, quotient_is_finite_dimensional,
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

/// Translate QPA's cached `IndecInjectiveModules(A)` attribute.
///
/// The method literally constructs the opposite algebra's ordered
/// indecomposable projectives and applies [`dual_of_module`] to each one.
/// Dual property propagation installs both indecomposability and injectivity
/// on every result, exactly as in QPA.
///
/// # Errors
///
/// Propagates opposite-algebra, projective-family, and module-dual failures.
pub fn indec_injective_modules<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Vec<PathAlgebraModule<F>>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_indec_injective_modules() {
        return Ok(cached);
    }
    let opposite = match algebra {
        QuiverAlgebra::Path(algebra) => QuiverAlgebra::Path(algebra.opposite()?),
        QuiverAlgebra::Quotient(algebra) => {
            QuiverAlgebra::Quotient(opposite_path_algebra(algebra)?)
        }
    };
    let opposite_projectives = indec_projective_modules(&opposite)?;
    let injectives = opposite_projectives
        .iter()
        .map(dual_of_module)
        .collect::<Result<Vec<_>, ModuleFunctorError>>()?;
    Ok(algebra.merge_indec_injective_modules(&injectives))
}

/// Translate QPA's cached `SimpleModules(A)` attribute.
///
/// One vertex-simple representation is constructed per vertex in quiver
/// order. Each has a single one in its dimension vector, no entered arrow
/// matrices, and QPA's `IsSimpleQPAModule` property installed together with
/// its true-method implications.
///
/// # Errors
///
/// Requires a finite path algebra or finite admissible quotient and
/// propagates the checked module constructor.
pub fn simple_modules<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Vec<PathAlgebraModule<F>>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_simple_modules() {
        return Ok(cached);
    }
    if finite_canonical_basis(algebra)?.is_none() {
        return Err(EnvelopingModuleError::InfiniteSimpleModuleAlgebra);
    }
    if matches!(algebra, QuiverAlgebra::Quotient(quotient) if !quotient.is_admissible()?) {
        return Err(EnvelopingModuleError::NonadmissibleSimpleModuleAlgebra);
    }
    let vertex_count = algebra.quiver().number_of_vertices();
    let mut simples = Vec::with_capacity(vertex_count);
    for vertex in 0..vertex_count {
        let mut dimensions = vec![0; vertex_count];
        dimensions[vertex] = 1;
        let simple = right_module_over_path_algebra(algebra.clone(), dimensions, &[])?
            .ok_or(EnvelopingModuleError::GeneratedSimpleViolatesRelations)?;
        let _ = simple.cache_is_simple_qpa_module(true);
        simples.push(simple);
    }
    Ok(algebra.merge_simple_modules(&simples))
}

/// Translate QPA's cached `AnnihilatorOfModule(M)` basis attribute.
///
/// For a nonzero module, rows use `BasisVectors(Basis(A))` for a path algebra
/// and GBNP nontips for a quotient, exactly as in QPA. Each row concatenates
/// the physical external coordinates of that basis element acting on every
/// standard module basis vector. GAP's `NullspaceMat` left-nullspace rows are
/// then used as coefficients of the returned algebra elements. The
/// zero-module branch returns `BasisVectors(Basis(A))` before constructing
/// this matrix.
///
/// # Errors
///
/// Requires a finite algebra basis and propagates module action, exact matrix,
/// and algebra-element arithmetic failures.
pub fn annihilator_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<QuiverAlgebraElement<F>>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_annihilator_of_module() {
        return Ok(cached);
    }
    let algebra = module.right_acting_algebra();
    let canonical_basis = finite_canonical_basis(algebra)?
        .ok_or(EnvelopingModuleError::InfiniteAnnihilatorAlgebra)?;
    let default_basis = finite_qpa_default_basis(algebra, &canonical_basis)?;
    if module.dimension() == 0 {
        return Ok(module.cache_annihilator_of_module(&default_basis));
    }
    let basis_vectors = match algebra {
        QuiverAlgebra::Path(_) => &default_basis,
        QuiverAlgebra::Quotient(_) => &canonical_basis,
    };
    let module_basis = module.basis_elements();
    let rows = basis_vectors
        .iter()
        .map(|algebra_element| {
            let mut row = Vec::new();
            for module_element in &module_basis {
                row.extend(
                    module_element
                        .act(algebra_element)?
                        .flat_external_representation(),
                );
            }
            Ok(row)
        })
        .collect::<Result<Vec<_>, PathAlgebraModuleError>>()?;
    let matrix = DenseMatrix::from_rows(rows, algebra.field().clone())?;
    let solutions = matrix.nullspace_rows();
    let mut annihilator = Vec::with_capacity(solutions.nrows());
    for solution in 0..solutions.nrows() {
        let mut element = algebra.zero();
        for (basis_element, coefficient) in basis_vectors.iter().zip(solutions.row(solution)) {
            element = element.add(&basis_element.scale_left(coefficient)?)?;
        }
        annihilator.push(element);
    }
    Ok(module.cache_annihilator_of_module(&annihilator))
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

/// Translate QPA's cached `IsInjectiveModule(M)` property.
///
/// This is the source expression `IsProjectiveModule(DualOfModule(M))`.
///
/// # Errors
///
/// Propagates module dualization and the projectivity dimension test.
pub fn is_injective_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<bool, EnvelopingModuleError> {
    if let Some(cached) = module.cached_is_injective_module() {
        return Ok(cached);
    }
    let dual = dual_of_module(module)?;
    let injective = is_projective_module(&dual)?;
    Ok(module.cache_is_injective_module(injective))
}

/// Translate QPA's `SocleSeries(M)` through the radical series of the dual.
///
/// QPA records the tops of successive radicals of `DualOfModule(M)` and then
/// reverses that list. Its zero branch returns the dimension vector directly,
/// represented by [`ModuleDimensionSeries::ZeroDimensionVector`].
///
/// # Errors
///
/// Propagates module dualization plus every top and radical construction.
pub fn socle_series<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<ModuleDimensionSeries, EnvelopingModuleError> {
    if module.dimension() == 0 {
        return Ok(ModuleDimensionSeries::ZeroDimensionVector(
            module.dimension_vector().to_vec(),
        ));
    }
    let mut current = dual_of_module(module)?;
    let mut layers = Vec::new();
    loop {
        layers.push(top_of_module(current.clone())?.dimension_vector().to_vec());
        current = radical_of_module(current)?;
        if current.dimension() == 0 {
            break;
        }
    }
    layers.reverse();
    Ok(ModuleDimensionSeries::Layers(layers))
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

/// Translate both QPA methods for `ElementInIndecProjective(A, m, s)`.
///
/// The result is deliberately an element of `OriginalPathAlgebra(A)`, even
/// when `A` is a quotient. QPA unwraps every quotient projective-basis element
/// to its stored parent representative before taking the linear combination;
/// this observable detail means the result of
/// `MatrixOfHomomorphismBetweenProjectives` over a quotient is not itself a
/// matrix of quotient-family elements.
///
/// # Errors
///
/// Requires a finite projective basis and a valid one-based projective
/// position. The physical external coordinates of `m`, including QPA's zero
/// space stand-ins, must have exactly the flattened projective-basis length.
pub fn element_in_indec_projective<F: Field>(
    algebra: &QuiverAlgebra<F>,
    element: &PathAlgebraModuleElement<F>,
    projective_vertex: usize,
) -> Result<PathAlgebraElement<F>, EnvelopingModuleError> {
    match algebra {
        QuiverAlgebra::Path(algebra) if !algebra.quiver().is_acyclic() => {
            return Err(EnvelopingModuleError::ElementInIndecProjectiveAlgebraNotFinite);
        }
        QuiverAlgebra::Quotient(algebra) if !algebra.is_admissible()? => {
            return Err(EnvelopingModuleError::ElementInIndecProjectiveAlgebraNotFinite);
        }
        QuiverAlgebra::Path(_) | QuiverAlgebra::Quotient(_) => {}
    }
    let bases =
        basis_of_projectives(algebra)?.ok_or(EnvelopingModuleError::ProjectiveBasisUnavailable)?;
    let projective = bases
        .get(projective_vertex.checked_sub(1).ok_or(
            EnvelopingModuleError::ProjectiveVertexOutOfRange {
                vertex: projective_vertex,
                vertex_count: bases.len(),
            },
        )?)
        .ok_or(EnvelopingModuleError::ProjectiveVertexOutOfRange {
            vertex: projective_vertex,
            vertex_count: bases.len(),
        })?;
    let zero = algebra.original_path_algebra().zero();
    let flattened_basis = projective
        .iter()
        .flat_map(|component| {
            if component.is_empty() {
                vec![zero.clone()]
            } else {
                component
                    .iter()
                    .map(|basis_element| match basis_element {
                        QuiverAlgebraElement::Path(element) => element.clone(),
                        QuiverAlgebraElement::Quotient(element) => element.representative().clone(),
                    })
                    .collect()
            }
        })
        .collect::<Vec<_>>();
    let coordinates = element.flat_external_representation();
    if coordinates.len() != flattened_basis.len() {
        return Err(
            EnvelopingModuleError::ProjectiveElementPhysicalCoordinateCount {
                actual: coordinates.len(),
                expected: flattened_basis.len(),
            },
        );
    }
    let mut result = zero;
    for (basis_element, coefficient) in flattened_basis.iter().zip(&coordinates) {
        result = result.add(&basis_element.scale_left(coefficient)?)?;
    }
    Ok(result)
}

/// Translate QPA's `ElementIn_vA_AsElementInIndecProj(A, m)`.
///
/// QPA splits `m` by right support, repeatedly removes its `Tip`, and looks up
/// that coefficient-bearing tip directly in `BasisOfProjectives`. This lookup
/// fails for a tip whose coefficient is not one; Rust preserves that source
/// bug as [`EnvelopingModuleError::ProjectiveLeadingTermMissing`].
/// The returned module is the cached indecomposable projective selected by the
/// unique left support of `m`.
///
/// # Errors
///
/// Rejects foreign, zero, or non-left-uniform elements and reports an internal
/// mismatch if a leading monomial is absent from QPA's projective basis.
pub fn element_in_v_a_as_element_in_indec_proj<F: Field>(
    algebra: &QuiverAlgebra<F>,
    element: &QuiverAlgebraElement<F>,
) -> Result<PathAlgebraModuleElement<F>, EnvelopingModuleError> {
    if !algebra.contains_element(element) {
        return Err(EnvelopingModuleError::ForeignProjectiveAlgebraElement);
    }
    if element.is_zero() {
        return Err(EnvelopingModuleError::ZeroProjectiveAlgebraElement);
    }
    let left_support = match (algebra, element) {
        (QuiverAlgebra::Path(algebra), QuiverAlgebraElement::Path(element)) => {
            algebra.left_support_of_quiver_algebra_element(element)?
        }
        (QuiverAlgebra::Quotient(algebra), QuiverAlgebraElement::Quotient(element)) => {
            algebra.left_support_of_quiver_algebra_element(element)?
        }
        _ => return Err(EnvelopingModuleError::ForeignProjectiveAlgebraElement),
    };
    let [left_vertex] = left_support.as_slice() else {
        return Err(EnvelopingModuleError::ProjectiveAlgebraElementNotLeftUniform);
    };
    let right_support = match (algebra, element) {
        (QuiverAlgebra::Path(algebra), QuiverAlgebraElement::Path(element)) => {
            algebra.right_support_of_quiver_algebra_element(element)?
        }
        (QuiverAlgebra::Quotient(algebra), QuiverAlgebraElement::Quotient(element)) => {
            algebra.right_support_of_quiver_algebra_element(element)?
        }
        _ => return Err(EnvelopingModuleError::ForeignProjectiveAlgebraElement),
    };
    let projectives = indec_projective_modules(algebra)?;
    let projective = projectives.get(*left_vertex - 1).ok_or(
        EnvelopingModuleError::ProjectiveVertexOutOfRange {
            vertex: *left_vertex,
            vertex_count: projectives.len(),
        },
    )?;
    let bases =
        basis_of_projectives(algebra)?.ok_or(EnvelopingModuleError::ProjectiveBasisUnavailable)?;
    let projective_basis = &bases[*left_vertex - 1];
    let field = algebra.field();
    let mut coordinates = projective
        .dimension_vector()
        .iter()
        .map(|&dimension| vec![field.zero(); dimension])
        .collect::<Vec<_>>();
    for &right_vertex in &right_support {
        let vertex = algebra.quiver().vertices().nth(right_vertex - 1).ok_or(
            EnvelopingModuleError::ProjectiveVertexOutOfRange {
                vertex: right_vertex,
                vertex_count: algebra.quiver().number_of_vertices(),
            },
        )?;
        let idempotent = include_in_path_algebra(&vertex.path(), algebra)?;
        let mut remainder = element.multiply(&idempotent)?;
        while !remainder.is_zero() {
            let tip = remainder.leading_term();
            let coefficient = tip.leading_coefficient();
            remainder = remainder.sub(&tip.scale_left(&coefficient)?)?;
            let position = projective_basis[right_vertex - 1]
                .iter()
                .position(|basis_element| basis_element == &tip)
                .ok_or(EnvelopingModuleError::ProjectiveLeadingTermMissing {
                    left_vertex: *left_vertex,
                    right_vertex,
                })?;
            coordinates[right_vertex - 1][position] =
                field.add(&coordinates[right_vertex - 1][position], &coefficient);
        }
    }
    Ok(projective.element_from_vertex_coordinates(coordinates)?)
}

fn projective_component_maps<F: Field>(
    module: &PathAlgebraModule<F>,
    projections: bool,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, EnvelopingModuleError> {
    if module.is_direct_sum_of_modules() {
        if projections {
            Ok(direct_sum_projections(module)?)
        } else {
            Ok(direct_sum_inclusions(module)?)
        }
    } else {
        Ok(vec![identity_module_mapping(module.clone())?])
    }
}

/// Translate QPA's `MatrixOfHomomorphismBetweenProjectives(f)`.
///
/// Rows index range summands and columns index source summands. Each entry is
/// the image of the corresponding distinguished source generator, projected
/// to a range summand and converted by [`element_in_indec_projective`]. As in
/// QPA, entries always belong to the original path algebra, including when
/// the modules act over a quotient.
///
/// # Errors
///
/// Propagates direct-sum attributes, minimal generators, map images, finite
/// projective bases, and parent path-algebra arithmetic failures.
pub fn matrix_of_homomorphism_between_projectives<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Vec<Vec<PathAlgebraElement<F>>>, EnvelopingModuleError> {
    let algebra = homomorphism.source().right_acting_algebra();
    let zero = algebra.original_path_algebra().zero();
    if homomorphism.range().dimension() == 0 {
        return Ok(vec![vec![zero]]);
    }
    let projections = projective_component_maps(homomorphism.range(), true)?;
    if homomorphism.source().dimension() == 0 {
        return Ok(projections.iter().map(|_| vec![zero.clone()]).collect());
    }
    let inclusions = projective_component_maps(homomorphism.source(), false)?;
    let target_vertices = projections
        .iter()
        .map(|projection| {
            let generator = minimal_generating_set_of_module(projection.range().clone())?
                .into_iter()
                .next()
                .ok_or(EnvelopingModuleError::EmptyProjectiveSummandGenerators)?;
            let support = generator.support_vertices();
            let [vertex] = support.as_slice() else {
                return Err(EnvelopingModuleError::ProjectiveGeneratorNotVertexUniform);
            };
            Ok(*vertex)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let mut columns = Vec::with_capacity(inclusions.len());
    for inclusion in &inclusions {
        let generator = minimal_generating_set_of_module(inclusion.source().clone())?
            .into_iter()
            .next()
            .ok_or(EnvelopingModuleError::EmptyProjectiveSummandGenerators)?;
        let source_element = inclusion.image(&generator)?;
        let image = homomorphism.image(&source_element)?;
        let mut column = Vec::with_capacity(projections.len());
        for (projection, &target_vertex) in projections.iter().zip(&target_vertices) {
            let component = projection.image(&image)?;
            column.push(element_in_indec_projective(
                algebra,
                &component,
                target_vertex,
            )?);
        }
        columns.push(column);
    }
    Ok((0..projections.len())
        .map(|row| columns.iter().map(|column| column[row].clone()).collect())
        .collect())
}

/// Translate QPA's `FromMatrixToHomomorphismOfProjectives(A, mat, vert1,
/// vert0)`.
///
/// `vert1` indexes source summands (matrix columns) and `vert0` indexes range
/// summands (matrix rows). Nonzero entries must be uniform, row-homogeneous on
/// left support, column-homogeneous on right support, and agree exactly with
/// the entered vertex lists; zero entries obtain both supports from those
/// lists, exactly as in the source method.
///
/// # Errors
///
/// Rejects malformed or support-incompatible matrices and propagates
/// projective construction, element conversion, `HomFromProjective`, and map
/// composition failures.
#[allow(clippy::too_many_lines)]
pub fn from_matrix_to_homomorphism_of_projectives<F: Field>(
    algebra: &QuiverAlgebra<F>,
    matrix: &[Vec<QuiverAlgebraElement<F>>],
    source_vertices: &[usize],
    range_vertices: &[usize],
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let Some(first_row) = matrix.first() else {
        return Err(EnvelopingModuleError::EmptyProjectiveElementMatrix);
    };
    if first_row.is_empty() {
        return Err(EnvelopingModuleError::EmptyProjectiveElementMatrix);
    }
    let column_count = first_row.len();
    if matrix.iter().any(|row| row.len() != column_count) {
        return Err(EnvelopingModuleError::RaggedProjectiveElementMatrix);
    }
    if matrix
        .iter()
        .flatten()
        .any(|element| !algebra.contains_element(element))
    {
        return Err(EnvelopingModuleError::ForeignProjectiveMatrixElement);
    }
    let mut detected_range_vertices = Vec::with_capacity(matrix.len());
    for (row_index, row) in matrix.iter().enumerate() {
        let mut support = Vec::new();
        for element in row {
            let left = match (algebra, element) {
                (QuiverAlgebra::Path(algebra), QuiverAlgebraElement::Path(element)) => {
                    algebra.left_support_of_quiver_algebra_element(element)?
                }
                (QuiverAlgebra::Quotient(algebra), QuiverAlgebraElement::Quotient(element)) => {
                    algebra.left_support_of_quiver_algebra_element(element)?
                }
                _ => return Err(EnvelopingModuleError::ForeignProjectiveMatrixElement),
            };
            if left.len() > 1 {
                return Err(EnvelopingModuleError::ProjectiveMatrixElementNotUniform);
            }
            support.extend(left);
        }
        if let Some(&first) = support.first() {
            if support.iter().any(|&vertex| vertex != first) {
                return Err(EnvelopingModuleError::ProjectiveMatrixRowNotHomogeneous {
                    row: row_index + 1,
                });
            }
            detected_range_vertices.push(first);
        } else {
            detected_range_vertices.push(*range_vertices.get(row_index).ok_or(
                EnvelopingModuleError::ProjectiveMatrixVertexListDimensions {
                    rows: matrix.len(),
                    columns: column_count,
                    range_vertices: range_vertices.len(),
                    source_vertices: source_vertices.len(),
                },
            )?);
        }
    }
    let mut detected_source_vertices = Vec::with_capacity(column_count);
    for column in 0..column_count {
        let mut support = Vec::new();
        for row in matrix {
            let element = &row[column];
            let right = match (algebra, element) {
                (QuiverAlgebra::Path(algebra), QuiverAlgebraElement::Path(element)) => {
                    algebra.right_support_of_quiver_algebra_element(element)?
                }
                (QuiverAlgebra::Quotient(algebra), QuiverAlgebraElement::Quotient(element)) => {
                    algebra.right_support_of_quiver_algebra_element(element)?
                }
                _ => return Err(EnvelopingModuleError::ForeignProjectiveMatrixElement),
            };
            if right.len() > 1 {
                return Err(EnvelopingModuleError::ProjectiveMatrixElementNotUniform);
            }
            support.extend(right);
        }
        if let Some(&first) = support.first() {
            if support.iter().any(|&vertex| vertex != first) {
                return Err(
                    EnvelopingModuleError::ProjectiveMatrixColumnNotHomogeneous {
                        column: column + 1,
                    },
                );
            }
            detected_source_vertices.push(first);
        } else {
            detected_source_vertices.push(*source_vertices.get(column).ok_or(
                EnvelopingModuleError::ProjectiveMatrixVertexListDimensions {
                    rows: matrix.len(),
                    columns: column_count,
                    range_vertices: range_vertices.len(),
                    source_vertices: source_vertices.len(),
                },
            )?);
        }
    }
    if detected_source_vertices != source_vertices || detected_range_vertices != range_vertices {
        return Err(EnvelopingModuleError::ProjectiveMatrixSupportDiscrepancy);
    }
    let vertex_count = algebra.quiver().number_of_vertices();
    for &vertex in source_vertices.iter().chain(range_vertices) {
        if !(1..=vertex_count).contains(&vertex) {
            return Err(EnvelopingModuleError::ProjectiveVertexOutOfRange {
                vertex,
                vertex_count,
            });
        }
    }

    let projectives = indec_projective_modules(algebra)?;
    let source_summands = source_vertices
        .iter()
        .map(|&vertex| projectives[vertex - 1].clone())
        .collect::<Vec<_>>();
    let range_summands = range_vertices
        .iter()
        .map(|&vertex| projectives[vertex - 1].clone())
        .collect::<Vec<_>>();
    let source = direct_sum_of_qpa_modules(&source_summands)?
        .ok_or(EnvelopingModuleError::InternalProjectiveDirectSumFailed)?;
    let range = direct_sum_of_qpa_modules(&range_summands)?
        .ok_or(EnvelopingModuleError::InternalProjectiveDirectSumFailed)?;
    let projections = direct_sum_projections(&source)?;
    let inclusions = direct_sum_inclusions(&range)?;
    let mut composites = Vec::with_capacity(matrix.len() * column_count);
    for (row_index, row) in matrix.iter().enumerate() {
        for (column, element) in row.iter().enumerate() {
            let component = if element.is_zero() {
                zero_module_mapping(
                    projections[column].range().clone(),
                    inclusions[row_index].source().clone(),
                )?
            } else {
                let target_element = element_in_v_a_as_element_in_indec_proj(algebra, element)?;
                hom_from_projective(&target_element, inclusions[row_index].source())?
                    .ok_or(EnvelopingModuleError::InternalProjectiveMatrixHomRejected)?
            };
            composites.push(
                projections[column]
                    .compose(&component)?
                    .compose(&inclusions[row_index])?,
            );
        }
    }
    let mut composites = composites.into_iter();
    let mut result = composites
        .next()
        .ok_or(EnvelopingModuleError::EmptyProjectiveElementMatrix)?;
    for composite in composites {
        result = result.add(&composite)?;
    }
    Ok(result)
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
        let zero = zero_module(module.right_acting_algebra())?;
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

/// Translate QPA's cached `1stSyzygy(M)` attribute.
///
/// The zero branch returns the acting algebra's identical cached
/// [`zero_module`]. Otherwise this is literally `Kernel(ProjectiveCover(M))`,
/// with the returned kernel module cached on `M` as GAP's attribute does.
///
/// # Errors
///
/// Propagates zero-module, projective-cover, and kernel construction failures.
pub fn first_syzygy<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_first_syzygy() {
        return Ok(cached);
    }
    let syzygy = if module.dimension() == 0 {
        zero_module(module.right_acting_algebra())?
    } else {
        let cover = projective_cover(module)?;
        kernel_of_module_homomorphism(&cover)?
    };
    Ok(module.cache_first_syzygy(&syzygy))
}

/// Translate QPA's cached `InjectiveEnvelope(M)` attribute.
///
/// The computation is the exact source expression
/// `DualOfModuleHomomorphism(ProjectiveCover(DualOfModule(M)))`. Reciprocal
/// dual-module caching makes the returned map start in the identical entered
/// module, including for the zero module.
///
/// # Errors
///
/// Propagates opposite-algebra dualization, projective-cover, and dual-map
/// construction failures.
pub fn injective_envelope<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_injective_envelope() {
        return Ok(cached);
    }
    let dual = dual_of_module(module)?;
    let dual_cover = projective_cover(&dual)?;
    let envelope = dual_of_module_homomorphism(&dual_cover)?;
    Ok(module.cache_injective_envelope(&envelope))
}

/// QPA's three-part `ExtOverAlgebra(M, N)` result.
///
/// The first component is the inclusion `Omega(M) -> P(M)`, and [`Self::basis`]
/// is QPA's ordered basis of `Ext^1(M, N)` represented inside
/// `Hom(Omega(M), N)`. When this basis is nonempty, [`Self::coefficients`]
/// translates the closure returned as QPA's third component.
#[derive(Clone, Debug)]
pub struct ExtOverAlgebraResult<F: Field> {
    syzygy_inclusion: PathAlgebraModuleHomomorphism<F>,
    basis: Vec<PathAlgebraModuleHomomorphism<F>>,
    coordinate_basis: Option<DenseMatrix<F>>,
}

impl<F: Field> ExtOverAlgebraResult<F> {
    /// Return QPA's kernel inclusion `Omega(M) -> P(M)`.
    #[must_use]
    pub const fn syzygy_inclusion(&self) -> &PathAlgebraModuleHomomorphism<F> {
        &self.syzygy_inclusion
    }

    /// Return QPA's ordered representatives of the `Ext^1(M, N)` basis.
    #[must_use]
    pub fn basis(&self) -> &[PathAlgebraModuleHomomorphism<F>] {
        &self.basis
    }

    /// Return the dimension of `Ext^1(M, N)`.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.basis.len()
    }

    /// Apply QPA's returned coefficient closure to a map in
    /// `Hom(Omega(M), N)`.
    ///
    /// `Ok(None)` corresponds to QPA returning `[]` instead of a closure when
    /// the Ext group is zero.
    ///
    /// # Errors
    ///
    /// The entered map's raw physical vector must use the coefficient field
    /// and lie in the Hom space used by the original computation. Exact row
    /// solving failures propagate.
    pub fn coefficients(
        &self,
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
    ) -> Result<Option<Vec<F::Element>>, EnvelopingModuleError> {
        let Some(coordinate_basis) = &self.coordinate_basis else {
            return Ok(None);
        };
        if homomorphism.source().right_acting_algebra().field() != coordinate_basis.field() {
            return Err(EnvelopingModuleError::ExtCoefficientFieldMismatch);
        }
        let coefficients = coordinate_basis
            .solution_row(&homomorphism.flat_vertex_matrix_entries())?
            .ok_or(EnvelopingModuleError::ExtMapOutsideHomSpace)?;
        Ok(Some(
            coefficients[coefficients.len() - self.basis.len()..].to_vec(),
        ))
    }
}

/// Translate QPA's `ExtOverAlgebra(M, N)` construction of `Ext^1(M, N)`.
///
/// This follows the source method through the projective cover, its kernel
/// inclusion, `Hom(P(M), N)` and `Hom(Omega(M), N)`, the image of the induced
/// restriction map, and GAP's `NaturalHomomorphismBySubspace`. Basis maps are
/// the exact complement representatives selected from `Basis(V)` after the
/// image basis, preserving QPA's raw physical map-vector order.
///
/// # Errors
///
/// Requires both modules over one algebra and propagates projective-cover,
/// kernel, Hom, composition, exact vector-space, and checked-map failures.
pub fn ext_over_algebra<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficients_in: &PathAlgebraModule<F>,
) -> Result<ExtOverAlgebraResult<F>, EnvelopingModuleError> {
    if module.right_acting_algebra() != coefficients_in.right_acting_algebra() {
        return Err(EnvelopingModuleError::ExtModulesOverDifferentAlgebras);
    }
    let cover = projective_cover(module)?;
    let syzygy_inclusion = kernel_inclusion(&cover)?;
    let syzygy = syzygy_inclusion.source();
    let projective_hom = hom_over_algebra(cover.source(), coefficients_in)?
        .ok_or(EnvelopingModuleError::InternalExtHomRejected)?;
    let syzygy_hom = hom_over_algebra(syzygy, coefficients_in)?
        .ok_or(EnvelopingModuleError::InternalExtHomRejected)?;
    let Some(first_homomorphism) = syzygy_hom.first() else {
        return Ok(ExtOverAlgebraResult {
            syzygy_inclusion,
            basis: Vec::new(),
            coordinate_basis: None,
        });
    };

    let field = module.right_acting_algebra().field().clone();
    let physical_length = first_homomorphism.flat_vertex_matrix_entries().len();
    let hom_rows = DenseMatrix::from_flat(
        syzygy_hom.len(),
        physical_length,
        syzygy_hom
            .iter()
            .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
            .collect(),
        field.clone(),
    )?;
    let hom_basis = hom_rows.semi_echelon_transformation().vectors;

    let image_maps = projective_hom
        .iter()
        .map(|homomorphism| syzygy_inclusion.compose(homomorphism))
        .collect::<Result<Vec<_>, _>>()?;
    let nonzero_image_maps = image_maps
        .iter()
        .filter(|homomorphism| !homomorphism.is_zero())
        .collect::<Vec<_>>();
    let image_basis = if nonzero_image_maps.is_empty() {
        DenseMatrix::from_flat(0, physical_length, Vec::new(), field.clone())?
    } else {
        DenseMatrix::from_flat(
            nonzero_image_maps.len(),
            physical_length,
            nonzero_image_maps
                .iter()
                .flat_map(|homomorphism| homomorphism.flat_vertex_matrix_entries())
                .collect(),
            field.clone(),
        )?
        .semi_echelon_transformation()
        .vectors
    };

    let quotient = quotient_complement_in_space(&hom_basis, &image_basis)?;
    if quotient.representatives.is_empty() {
        return Ok(ExtOverAlgebraResult {
            syzygy_inclusion,
            basis: Vec::new(),
            coordinate_basis: None,
        });
    }
    let basis = quotient
        .representatives
        .iter()
        .map(|row| module_map_from_physical_flat(syzygy, coefficients_in, row, &field))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(ExtOverAlgebraResult {
        syzygy_inclusion,
        basis,
        coordinate_basis: Some(quotient.coordinate_basis),
    })
}

/// Translate QPA's cached `IsRigidModule(M)` property.
///
/// This is exactly the source predicate
/// `Length(ExtOverAlgebra(M, M)[2]) = 0`.
///
/// # Errors
///
/// Propagates every projective-cover, Hom, and quotient-space failure from
/// [`ext_over_algebra`].
pub fn is_rigid_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<bool, EnvelopingModuleError> {
    if let Some(cached) = module.cached_is_rigid_module() {
        return Ok(cached);
    }
    let rigid = ext_over_algebra(module, module)?.dimension() == 0;
    Ok(module.cache_is_rigid_module(rigid))
}

struct ExtQuotientComplement<F: Field> {
    representatives: Vec<Vec<F::Element>>,
    coordinate_basis: DenseMatrix<F>,
}

fn quotient_complement_in_space<F: Field>(
    space_basis: &DenseMatrix<F>,
    subspace_basis: &DenseMatrix<F>,
) -> Result<ExtQuotientComplement<F>, EnvelopingModuleError> {
    let field = space_basis.field().clone();
    let columns = space_basis.ncols();
    let mut mutable_basis = subspace_basis.clone();
    let mut complement = Vec::new();
    for row in 0..space_basis.nrows() {
        let candidate_row = space_basis.row(row);
        let mut entries = mutable_basis.as_flat().to_vec();
        entries.extend_from_slice(candidate_row);
        let candidate =
            DenseMatrix::from_flat(mutable_basis.nrows() + 1, columns, entries, field.clone())?
                .semi_echelon_transformation()
                .vectors;
        if candidate.nrows() > mutable_basis.nrows() {
            complement.push(candidate_row.to_vec());
            mutable_basis = candidate;
        }
    }
    if subspace_basis.nrows() + complement.len() != space_basis.nrows() {
        return Err(EnvelopingModuleError::IncompleteExtQuotientBasis);
    }
    let mut entries = subspace_basis.as_flat().to_vec();
    entries.extend(complement.iter().flatten().cloned());
    let coordinate_basis = DenseMatrix::from_flat(space_basis.nrows(), columns, entries, field)?;
    Ok(ExtQuotientComplement {
        representatives: complement,
        coordinate_basis,
    })
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

/// Translate QPA's `EndModuloProjOverAlgebra(M)` natural algebra
/// homomorphism.
///
/// This constructs `EndOverAlgebra(M)`, embeds QPA's ordered basis of maps
/// factoring through projectives as block matrices, forms their two-sided
/// ideal, and returns the natural quotient map onto the stable endomorphism
/// algebra.
///
/// # Errors
///
/// QPA rejects the zero module. Otherwise projective factoring,
/// endomorphism embedding, ideal closure, and exact quotient failures
/// propagate.
pub fn end_modulo_projectives<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<ModuleEndomorphismQuotientHomomorphism<F>, EnvelopingModuleError> {
    if module.dimension() == 0 {
        return Err(EnvelopingModuleError::ZeroStableEndomorphismModule);
    }
    let endomorphisms = end_over_algebra(module)?;
    let factor_matrices = hom_factoring_through_projectives(module, module)?
        .iter()
        .map(from_hom_mm_to_end_m)
        .collect::<Result<Vec<_>, _>>()?;
    Ok(natural_endomorphism_quotient_by_ideal(
        &endomorphisms,
        &factor_matrices,
    )?)
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
/// retained, and the result acts over QPA's opposite algebra. QPA's zero
/// branch returns the opposite algebra's cached [`zero_module`]. For a
/// nonzero module, QPA explicitly installs the reciprocal dual attribute on
/// the result and propagates already-true indecomposable, projective, and
/// injective properties in the corresponding dual directions; its
/// zero-module early return does not.
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
    if module.dimension() == 0 {
        let dual = zero_module(&opposite)?;
        return Ok(module.cache_dual_module(&dual));
    }
    let dimensions = module.dimension_vector().to_vec();
    let mut specification = Vec::new();
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
    let dual = right_module_over_path_algebra(opposite, dimensions, &specification)?
        .ok_or(ModuleFunctorError::GeneratedDualViolatesRelations)?;
    let dual = module.install_dual_module(&dual)?;
    if module.cached_is_indecomposable_module() == Some(true) {
        let _ = dual.cache_is_indecomposable_module(true);
    }
    if module.cached_is_projective_module() == Some(true) {
        let _ = dual.cache_is_injective_module(true);
    }
    if module.cached_is_injective_module() == Some(true) {
        let _ = dual.cache_is_projective_module(true);
    }
    Ok(dual)
}

/// Translate QPA's cached `TransposeOfModule(M)` construction.
///
/// This follows `functors.gi` directly: expand the projective-cover map in
/// the global canonical algebra basis, take its GAP row nullspace, quotient
/// the syzygy by the arrow-generated radical with GAP's mutable-basis
/// complement order, transpose and opposite every resulting presentation
/// entry, and form the cokernel inside the corresponding opposite
/// projective cover.
///
/// # Errors
///
/// Requires the same finite projective/canonical-basis domain as QPA and
/// propagates exact algebra, nullspace, direct-sum, subrepresentation, and
/// cokernel failures.
#[allow(clippy::too_many_lines)]
pub fn transpose_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_transpose_of_module() {
        return Ok(cached);
    }
    let algebra = module.right_acting_algebra();
    let opposite = match algebra {
        QuiverAlgebra::Path(path_algebra) => QuiverAlgebra::Path(path_algebra.opposite()?),
        QuiverAlgebra::Quotient(quotient) => {
            QuiverAlgebra::Quotient(opposite_path_algebra(quotient)?)
        }
    };
    if is_projective_module(module)? {
        let zero = zero_module(&opposite)?;
        return Ok(module.cache_transpose_of_module(&zero));
    }

    let basis =
        finite_canonical_basis(algebra)?.ok_or(EnvelopingModuleError::InfiniteTransposeAlgebra)?;
    let field = algebra.field().clone();
    let vertex_elements = algebra
        .quiver()
        .vertices()
        .map(|vertex| include_in_path_algebra(&vertex.path(), algebra))
        .collect::<Result<Vec<_>, _>>()?;
    let arrow_elements = algebra
        .quiver()
        .arrows()
        .map(|arrow| include_in_path_algebra(&arrow.path(), algebra))
        .collect::<Result<Vec<_>, _>>()?;

    // QPA's local BB uses global CanonicalBasis(A) order, filtered by left
    // vertex multiplication; it is intentionally not the target-grouped
    // BasisOfProjectives order.
    let mut projective_basis_indices = vec![Vec::new(); vertex_elements.len()];
    for (vertex, vertex_element) in vertex_elements.iter().enumerate() {
        for (basis_index, basis_element) in basis.iter().enumerate() {
            if !vertex_element.multiply(basis_element)?.is_zero() {
                projective_basis_indices[vertex].push(basis_index);
            }
        }
    }

    // CanonicalBasis(M) is assigned in the source even though the value is
    // not subsequently read.
    let _module_basis = module.basis_elements();
    let generators = minimal_generating_set_of_module(module.clone())?;
    let mut cover_vertices = Vec::with_capacity(generators.len());
    let mut cover_basis_images = Vec::new();
    for generator in &generators {
        let mut supported_vertex = None;
        for (vertex, vertex_element) in vertex_elements.iter().enumerate() {
            if !generator.act(vertex_element)?.is_zero() {
                supported_vertex = Some(vertex);
                break;
            }
        }
        let vertex =
            supported_vertex.ok_or(EnvelopingModuleError::MissingTransposeGeneratorSupport)?;
        cover_vertices.push(vertex);
        for &basis_index in &projective_basis_indices[vertex] {
            cover_basis_images.push(generator.act(&basis[basis_index])?);
        }
    }

    let cover_matrix = DenseMatrix::from_rows(
        cover_basis_images
            .iter()
            .map(PathAlgebraModuleElement::flat_external_representation)
            .collect(),
        field.clone(),
    )?;
    let solution_rows = cover_matrix.nullspace_rows();
    if solution_rows.nrows() == 0 {
        let zero = zero_module(&opposite)?;
        return Ok(module.cache_transpose_of_module(&zero));
    }
    let solution_vectors = projective_algebra_vectors_from_rows(
        algebra,
        &basis,
        &projective_basis_indices,
        &cover_vertices,
        (0..solution_rows.nrows()).map(|row| solution_rows.row(row)),
    )?;

    let mut radical_vectors = Vec::new();
    for solution in &solution_vectors {
        for arrow in &arrow_elements {
            let product = solution
                .iter()
                .map(|entry| entry.multiply(arrow))
                .collect::<Result<Vec<_>, _>>()?;
            if product.iter().any(|entry| !entry.is_zero()) {
                radical_vectors.push(product);
            }
        }
    }
    let radical_rows = radical_vectors
        .iter()
        .map(|vector| {
            let mut row = Vec::with_capacity(solution_rows.ncols());
            for (entry, &cover_vertex) in vector.iter().zip(&cover_vertices) {
                let coefficients = canonical_basis_coefficients(algebra, &basis, entry)?;
                row.extend(
                    projective_basis_indices[cover_vertex]
                        .iter()
                        .map(|&basis_index| coefficients[basis_index].clone()),
                );
            }
            Ok(row)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;

    let syzygy_basis = solution_rows.semi_echelon_transformation().vectors;
    let radical_basis = if radical_rows.is_empty() {
        DenseMatrix::from_flat(0, solution_rows.ncols(), Vec::new(), field.clone())?
    } else {
        DenseMatrix::from_rows(radical_rows, field.clone())?
            .semi_echelon_transformation()
            .vectors
    };
    let top_rows = if radical_basis.nrows() == 0 {
        // GAP's high-rank trivial-subspace method returns IdentityMapping(V).
        // The source then asks for CanonicalBasis(Range(f)) = CanonicalBasis(V),
        // which uses full triangulation rather than the ordinary Basis(V)
        // used by the nontrivial quotient method.
        let canonical = solution_rows.canonical_row_space_basis();
        (0..canonical.nrows())
            .map(|row| canonical.row(row).to_vec())
            .collect()
    } else {
        quotient_complement_in_space(&syzygy_basis, &radical_basis)?.representatives
    };
    let top_vectors = projective_algebra_vectors_from_rows(
        algebra,
        &basis,
        &projective_basis_indices,
        &cover_vertices,
        top_rows.iter().map(Vec::as_slice),
    )?;

    let projectives = indec_projective_modules(algebra)?;
    let cover_summands = cover_vertices
        .iter()
        .map(|&vertex| projectives[vertex].clone())
        .collect::<Vec<_>>();
    let _projective_cover = direct_sum_of_qpa_modules(&cover_summands)?
        .ok_or(EnvelopingModuleError::InternalTransposeDirectSumFailed)?;

    let mut syzygy_generator_vertices = Vec::with_capacity(top_vectors.len());
    for generator in &top_vectors {
        let mut support = Vec::new();
        for (vertex, vertex_element) in vertex_elements.iter().enumerate() {
            let supported = generator.iter().try_fold(false, |supported, entry| {
                Ok::<_, EnvelopingModuleError>(
                    supported || !entry.multiply(vertex_element)?.is_zero(),
                )
            })?;
            if supported {
                support.push(vertex);
            }
        }
        let [vertex] = support.as_slice() else {
            return Err(EnvelopingModuleError::NonuniformTransposeSyzygyGenerator);
        };
        syzygy_generator_vertices.push(*vertex);
    }

    let opposite_projectives = indec_projective_modules(&opposite)?;
    let opposite_cover_summands = syzygy_generator_vertices
        .iter()
        .map(|&vertex| opposite_projectives[vertex].clone())
        .collect::<Vec<_>>();
    let opposite_cover = direct_sum_of_qpa_modules(&opposite_cover_summands)?
        .ok_or(EnvelopingModuleError::InternalTransposeDirectSumFailed)?;
    let inclusions = direct_sum_inclusions(&opposite_cover)?;
    let embedded_generators = inclusions
        .iter()
        .zip(&opposite_cover_summands)
        .map(|(inclusion, summand)| {
            let generator = minimal_generating_set_of_module(summand.clone())?
                .into_iter()
                .next()
                .ok_or(EnvelopingModuleError::EmptyTransposeProjectiveGenerator)?;
            Ok(inclusion.image(&generator)?)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;

    let mut image_generators = Vec::with_capacity(generators.len());
    for cover_generator in 0..generators.len() {
        let mut image = opposite_cover.zero_element();
        for (top_generator, embedded_generator) in top_vectors.iter().zip(&embedded_generators) {
            let opposite_entry = opposite_quiver_algebra_element(&top_generator[cover_generator])?;
            image = image.add(&embedded_generator.act(&opposite_entry)?)?;
        }
        image_generators.push(image);
    }
    let inclusion = subrepresentation_inclusion(opposite_cover, &image_generators)?;
    let transpose = cokernel_of_module_homomorphism(&inclusion)?;
    if transpose.dimension() != 0 && module.cached_is_indecomposable_module() == Some(true) {
        let _ = transpose.cache_is_indecomposable_module(true);
    }
    Ok(module.cache_transpose_of_module(&transpose))
}

fn projective_algebra_vectors_from_rows<'a, F: Field + 'a>(
    algebra: &QuiverAlgebra<F>,
    basis: &[QuiverAlgebraElement<F>],
    projective_basis_indices: &[Vec<usize>],
    cover_vertices: &[usize],
    rows: impl IntoIterator<Item = &'a [F::Element]>,
) -> Result<Vec<Vec<QuiverAlgebraElement<F>>>, EnvelopingModuleError>
where
    F::Element: 'a,
{
    rows.into_iter()
        .map(|row| {
            let mut offset = 0;
            cover_vertices
                .iter()
                .map(|&vertex| {
                    let mut element = algebra.zero();
                    for &basis_index in &projective_basis_indices[vertex] {
                        let term = basis[basis_index].scale_left(&row[offset])?;
                        element = element.add(&term)?;
                        offset += 1;
                    }
                    Ok(element)
                })
                .collect::<Result<Vec<_>, EnvelopingModuleError>>()
        })
        .collect()
}

fn opposite_quiver_algebra_element<F: Field>(
    element: &QuiverAlgebraElement<F>,
) -> Result<QuiverAlgebraElement<F>, EnvelopingModuleError> {
    match element {
        QuiverAlgebraElement::Path(element) => Ok(QuiverAlgebraElement::Path(element.opposite()?)),
        QuiverAlgebraElement::Quotient(element) => Ok(QuiverAlgebraElement::Quotient(
            opposite_quotient_element(element)?,
        )),
    }
}

/// Translate QPA's cached module `DTr(M) = DualOfModule(TransposeOfModule(M))`.
///
/// # Errors
///
/// Propagates transpose and dual construction.
pub fn d_tr<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_d_tr() {
        return Ok(cached);
    }
    let transpose = transpose_of_module(module)?;
    let result = dual_of_module(&transpose)?;
    Ok(module.cache_d_tr(&result))
}

/// Translate QPA's cached module `TrD(M) = TransposeOfModule(DualOfModule(M))`.
///
/// # Errors
///
/// Propagates dual and transpose construction.
pub fn tr_d<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_tr_d() {
        return Ok(cached);
    }
    let dual = dual_of_module(module)?;
    let result = transpose_of_module(&dual)?;
    Ok(module.cache_tr_d(&result))
}

/// Translate QPA's integer-power `DTr(M, n)` dispatch.
///
/// Negative powers apply [`tr_d`], while zero returns the identical entered
/// module.
///
/// # Errors
///
/// Propagates every iteration's functor construction.
pub fn d_tr_iterated<F: Field>(
    module: &PathAlgebraModule<F>,
    exponent: i64,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    let mut result = module.clone();
    if exponent < 0 {
        for _ in 0..exponent.unsigned_abs() {
            result = tr_d(&result)?;
        }
    } else {
        for _ in 0..exponent.unsigned_abs() {
            result = d_tr(&result)?;
        }
    }
    Ok(result)
}

/// Translate QPA's integer-power `TrD(M, n)` dispatch.
///
/// Negative powers apply [`d_tr`], while zero returns the identical entered
/// module.
///
/// # Errors
///
/// Propagates every iteration's functor construction.
pub fn tr_d_iterated<F: Field>(
    module: &PathAlgebraModule<F>,
    exponent: i64,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    let mut result = module.clone();
    if exponent < 0 {
        for _ in 0..exponent.unsigned_abs() {
            result = d_tr(&result)?;
        }
    } else {
        for _ in 0..exponent.unsigned_abs() {
            result = tr_d(&result)?;
        }
    }
    Ok(result)
}

/// Translate QPA's cached `IsTauRigidModule(M)` property.
///
/// This is exactly `Length(HomOverAlgebra(M, DTr(M))) = 0`.
///
/// # Errors
///
/// Propagates `DTr` and Hom-space construction.
pub fn is_tau_rigid_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<bool, EnvelopingModuleError> {
    if let Some(cached) = module.cached_is_tau_rigid_module() {
        return Ok(cached);
    }
    let translate = d_tr(module)?;
    let hom = hom_over_algebra(module, &translate)?
        .ok_or(EnvelopingModuleError::InternalTauRigidHomRejected)?;
    Ok(module.cache_is_tau_rigid_module(hom.is_empty()))
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

/// Translate QPA's cached `TransposeOfModuleHomomorphism(f)` attribute.
///
/// This follows QPA's projective-presentation construction literally. It
/// lifts `f` to the projective covers and their first syzygy covers, applies
/// the contravariant star functor to all four maps, and takes the induced map
/// on cokernels.
///
/// # Errors
///
/// Propagates projective-cover, kernel, lifting, star, and cokernel-map
/// construction failures. The internal `fail` variants can occur only if a
/// square produced by the source construction ceases to commute.
pub fn transpose_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_transpose_homomorphism() {
        return Ok(cached);
    }

    let d0 = projective_cover(homomorphism.source())?;
    let d0_prime = projective_cover(homomorphism.range())?;
    let ker_d0 = kernel_inclusion(&d0)?;
    let ker_d0_prime = kernel_inclusion(&d0_prime)?;
    let d1_temp = projective_cover(ker_d0.source())?;
    let d1_prime_temp = projective_cover(ker_d0_prime.source())?;
    let d1 = d1_temp.compose(&ker_d0)?;
    let d1_prime = d1_prime_temp.compose(&ker_d0_prime)?;

    let f0 = lifting_morphism_from_projective(&d0_prime, &d0.compose(homomorphism)?)?
        .ok_or(EnvelopingModuleError::InternalTransposeHomomorphismLiftFailed)?;
    let f1_prime = morphism_on_kernel(&d0, &d0_prime, &f0, homomorphism)?
        .ok_or(EnvelopingModuleError::InternalTransposeHomomorphismKernelSquareFailed)?;
    let f1 = lifting_morphism_from_projective(&d1_prime_temp, &d1_temp.compose(&f1_prime)?)?
        .ok_or(EnvelopingModuleError::InternalTransposeHomomorphismLiftFailed)?;

    let d1_star = star_of_module_homomorphism(&d1)?;
    let d1_prime_star = star_of_module_homomorphism(&d1_prime)?;
    let f1_star = star_of_module_homomorphism(&f1)?;
    let f0_star = star_of_module_homomorphism(&f0)?;
    let transpose = morphism_on_cokernel(&d1_prime_star, &d1_star, &f0_star, &f1_star)?
        .ok_or(EnvelopingModuleError::InternalTransposeHomomorphismCokernelSquareFailed)?;

    Ok(homomorphism.cache_transpose_homomorphism(&transpose))
}

/// Translate QPA's cached `DTr(f)` attribute on module homomorphisms.
///
/// QPA computes the dual of the transposed map, then reconstructs it over the
/// exact cached `DTr` source and range module objects.
///
/// # Errors
///
/// Propagates transpose, dual, module `DTr`, and checked-map construction.
pub fn d_tr_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_d_tr_homomorphism() {
        return Ok(cached);
    }
    let transpose = transpose_of_module_homomorphism(homomorphism)?;
    let generated = dual_of_module_homomorphism(&transpose)?;
    let source = d_tr(homomorphism.source())?;
    let range = d_tr(homomorphism.range())?;
    let result =
        right_module_hom_over_algebra(source, range, generated.vertex_matrices().to_vec())?;
    Ok(homomorphism.cache_d_tr_homomorphism(&result))
}

/// Translate QPA's cached `TrD(f)` attribute on module homomorphisms.
///
/// QPA transposes the dual map, then reconstructs it over the exact cached
/// `TrD` source and range module objects.
///
/// # Errors
///
/// Propagates dual, transpose, module `TrD`, and checked-map construction.
pub fn tr_d_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_tr_d_homomorphism() {
        return Ok(cached);
    }
    let dual = dual_of_module_homomorphism(homomorphism)?;
    let generated = transpose_of_module_homomorphism(&dual)?;
    let source = tr_d(homomorphism.source())?;
    let range = tr_d(homomorphism.range())?;
    let result =
        right_module_hom_over_algebra(source, range, generated.vertex_matrices().to_vec())?;
    Ok(homomorphism.cache_tr_d_homomorphism(&result))
}

/// Translate QPA's integer-power `DTr(f, n)` method.
///
/// Negative powers apply [`tr_d_homomorphism`], while zero returns the
/// identical entered mapping object.
///
/// # Errors
///
/// Propagates each iterated homomorphism-functor construction.
pub fn d_tr_homomorphism_iterated<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    exponent: i64,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let mut result = homomorphism.clone();
    if exponent < 0 {
        for _ in 0..exponent.unsigned_abs() {
            result = tr_d_homomorphism(&result)?;
        }
    } else {
        for _ in 0..exponent.unsigned_abs() {
            result = d_tr_homomorphism(&result)?;
        }
    }
    Ok(result)
}

/// Translate QPA's integer-power `TrD(f, n)` method.
///
/// Negative powers apply [`d_tr_homomorphism`], while zero returns the
/// identical entered mapping object.
///
/// # Errors
///
/// Propagates each iterated homomorphism-functor construction.
pub fn tr_d_homomorphism_iterated<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    exponent: i64,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let mut result = homomorphism.clone();
    if exponent < 0 {
        for _ in 0..exponent.unsigned_abs() {
            result = d_tr_homomorphism(&result)?;
        }
    } else {
        for _ in 0..exponent.unsigned_abs() {
            result = tr_d_homomorphism(&result)?;
        }
    }
    Ok(result)
}

/// Translate QPA's cached `StarOfModule(M) = Hom_A(M, A)` construction.
///
/// The vertex at `i` is represented by QPA's ordered basis of
/// `Hom_A(M, e_i A)`.  For every arrow `i -> j`, the opposite-arrow action is
/// obtained by composing `Hom_A(M, e_j A)` with QPA's map `e_j A -> e_i A`
/// induced by left multiplication by that arrow, then solving in the original
/// raw physical Hom basis.
///
/// # Errors
///
/// Requires QPA's finite projective bases and propagates opposite-algebra,
/// Hom, projective-map, exact coordinate, and checked-module failures.
pub fn star_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_star_of_module() {
        return Ok(cached);
    }
    let algebra = module.right_acting_algebra();
    let opposite = match algebra {
        QuiverAlgebra::Path(path_algebra) => QuiverAlgebra::Path(path_algebra.opposite()?),
        QuiverAlgebra::Quotient(quotient) => {
            QuiverAlgebra::Quotient(opposite_path_algebra(quotient)?)
        }
    };
    let field = algebra.field().clone();
    let projectives = indec_projective_modules(algebra)?;
    let hom_spaces = projectives
        .iter()
        .map(|projective| {
            hom_over_algebra(module, projective)?
                .ok_or(EnvelopingModuleError::InternalStarHomRejected)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let hom_bases = hom_spaces
        .iter()
        .map(|basis| {
            let Some(first) = basis.first() else {
                return Ok(None);
            };
            Ok(Some(DenseMatrix::from_flat(
                basis.len(),
                first.flat_vertex_matrix_entries().len(),
                basis
                    .iter()
                    .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
                    .collect(),
                field.clone(),
            )?))
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let dimensions = hom_spaces.iter().map(Vec::len).collect::<Vec<_>>();

    // The source computes these subspace bases even though the local `Bproj`
    // value is not subsequently read.
    let _projective_bases =
        basis_of_projectives(algebra)?.ok_or(EnvelopingModuleError::ProjectiveBasisUnavailable)?;

    let arrows = algebra.quiver().arrows().collect::<Vec<_>>();
    let opposite_arrows = opposite.quiver().arrows().collect::<Vec<_>>();
    let mut left_multiplications = Vec::with_capacity(arrows.len());
    for arrow in &arrows {
        let source = arrow.source().index().get() - 1;
        let generator = minimal_generating_set_of_module(projectives[source].clone())?
            .into_iter()
            .next()
            .ok_or(EnvelopingModuleError::EmptyStarProjectiveGenerator)?;
        let arrow_element = include_in_path_algebra(&arrow.path(), algebra)?;
        let image = generator.act(&arrow_element)?;
        let left_multiplication = hom_from_projective(&image, &projectives[source])?
            .ok_or(EnvelopingModuleError::StarArrowMapNotUniform)?;
        left_multiplications.push(left_multiplication);
    }

    let mut specification = Vec::new();
    for (position, arrow) in arrows.iter().enumerate() {
        let source = arrow.source().index().get() - 1;
        let target = arrow.target().index().get() - 1;
        if dimensions[source] == 0 || dimensions[target] == 0 {
            continue;
        }
        let coordinate_basis = hom_bases[source]
            .as_ref()
            .ok_or(EnvelopingModuleError::MissingStarCoordinateBasis)?;
        let rows = hom_spaces[target]
            .iter()
            .map(|homomorphism| {
                let composite = homomorphism.compose(&left_multiplications[position])?;
                coordinate_basis
                    .solution_row(&composite.flat_vertex_matrix_entries())?
                    .ok_or(EnvelopingModuleError::StarCompositionOutsideHomBasis)
            })
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
        specification.push(NamedArrowMatrix::new(
            opposite_arrows[position].name(),
            DenseMatrix::from_rows(rows, field.clone())?,
        ));
    }
    let star = right_module_over_path_algebra(opposite, dimensions, &specification)?
        .ok_or(EnvelopingModuleError::GeneratedStarViolatesRelations)?;
    Ok(module.cache_star_of_module(&star))
}

/// Translate QPA's cached contravariant `StarOfModuleHomomorphism(f)`.
///
/// # Errors
///
/// Propagates both star-module constructions, Hom composition, exact
/// coordinate solving, and the checked induced-map constructor.
pub fn star_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_star_homomorphism() {
        return Ok(cached);
    }
    let source = homomorphism.source();
    let range = homomorphism.range();
    let star_source = star_of_module(source)?;
    let star_range = star_of_module(range)?;
    let algebra = source.right_acting_algebra();
    let field = algebra.field().clone();
    let projectives = indec_projective_modules(algebra)?;

    let range_homs = projectives
        .iter()
        .map(|projective| {
            hom_over_algebra(range, projective)?
                .ok_or(EnvelopingModuleError::InternalStarHomRejected)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let source_homs = projectives
        .iter()
        .map(|projective| {
            hom_over_algebra(source, projective)?
                .ok_or(EnvelopingModuleError::InternalStarHomRejected)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;

    let mut maps = Vec::with_capacity(projectives.len());
    for (vertex, (range_basis, source_basis)) in range_homs.iter().zip(&source_homs).enumerate() {
        if !range_basis.is_empty() && !source_basis.is_empty() {
            let coordinate_basis = DenseMatrix::from_flat(
                source_basis.len(),
                source_basis[0].flat_vertex_matrix_entries().len(),
                source_basis
                    .iter()
                    .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
                    .collect(),
                field.clone(),
            )?;
            let rows = range_basis
                .iter()
                .map(|basis_map| {
                    let composite = homomorphism.compose(basis_map)?;
                    coordinate_basis
                        .solution_row(&composite.flat_vertex_matrix_entries())?
                        .ok_or(EnvelopingModuleError::StarCompositionOutsideHomBasis)
                })
                .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
            maps.push(DenseMatrix::from_rows(rows, field.clone())?);
        } else {
            maps.push(DenseMatrix::zero(
                star_range.dimension_vector()[vertex].max(1),
                star_source.dimension_vector()[vertex].max(1),
                field.clone(),
            ));
        }
    }
    let star = right_module_hom_over_algebra(star_range, star_source, maps)?;
    Ok(homomorphism.cache_star_homomorphism(&star))
}

/// Translate QPA's cached `NakayamaFunctorOfModule(M)` composition.
///
/// # Errors
///
/// Propagates [`star_of_module`] and [`dual_of_module`].
pub fn nakayama_functor_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_nakayama_functor_of_module() {
        return Ok(cached);
    }
    let star = star_of_module(module)?;
    let nakayama = dual_of_module(&star)?;
    Ok(module.cache_nakayama_functor_of_module(&nakayama))
}

/// Translate QPA's cached `OppositeNakayamaFunctorOfModule(M)` composition.
///
/// # Errors
///
/// Propagates [`dual_of_module`] and [`star_of_module`].
pub fn opposite_nakayama_functor_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_opposite_nakayama_functor_of_module() {
        return Ok(cached);
    }
    let dual = dual_of_module(module)?;
    let nakayama = star_of_module(&dual)?;
    Ok(module.cache_opposite_nakayama_functor_of_module(&nakayama))
}

/// Translate QPA's cached Nakayama functor on a module homomorphism.
///
/// # Errors
///
/// Propagates star and dual homomorphism construction.
pub fn nakayama_functor_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_nakayama_functor_homomorphism() {
        return Ok(cached);
    }
    let star = star_of_module_homomorphism(homomorphism)?;
    let nakayama = dual_of_module_homomorphism(&star)?;
    Ok(homomorphism.cache_nakayama_functor_homomorphism(&nakayama))
}

/// Translate QPA's cached opposite Nakayama functor on a homomorphism.
///
/// # Errors
///
/// Propagates dual and star homomorphism construction.
pub fn opposite_nakayama_functor_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_opposite_nakayama_functor_homomorphism() {
        return Ok(cached);
    }
    let dual = dual_of_module_homomorphism(homomorphism)?;
    let nakayama = star_of_module_homomorphism(&dual)?;
    Ok(homomorphism.cache_opposite_nakayama_functor_homomorphism(&nakayama))
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

fn finite_qpa_default_basis<F: Field>(
    algebra: &QuiverAlgebra<F>,
    canonical_basis: &[QuiverAlgebraElement<F>],
) -> Result<Vec<QuiverAlgebraElement<F>>, EnvelopingModuleError> {
    let vertex_count = algebra.quiver().number_of_vertices();
    if canonical_basis.len() < vertex_count {
        return Err(EnvelopingModuleError::MissingCanonicalVertexBasis);
    }
    for (position, vertex) in algebra.quiver().vertices().enumerate() {
        let expected = include_in_path_algebra(&vertex.path(), algebra)?;
        if canonical_basis[position] != expected {
            return Err(EnvelopingModuleError::MissingCanonicalVertexBasis);
        }
    }

    // GAP's generic `Basis(A)` for an algebra with one first inserts `One(A)`
    // and then semi-echelonizes the vertex generators in vertex order.  In
    // canonical path coordinates this replaces e_1,...,e_n by the suffix
    // sums e_1+...+e_n, e_2+...+e_n, ..., e_n; all nontrivial paths retain
    // canonical order.  This is observably different from `CanonicalBasis(A)`.
    let mut suffix = algebra.zero();
    let mut vertex_basis = Vec::with_capacity(vertex_count);
    for element in canonical_basis[..vertex_count].iter().rev() {
        suffix = element.add(&suffix)?;
        vertex_basis.push(suffix.clone());
    }
    vertex_basis.reverse();
    vertex_basis.extend_from_slice(&canonical_basis[vertex_count..]);
    Ok(vertex_basis)
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
    /// `SimpleModules` rejects infinite-dimensional input before constructing
    /// vertex representations.
    #[error("SimpleModules requires a finite-dimensional algebra")]
    InfiniteSimpleModuleAlgebra,
    /// An annihilator basis requires QPA's finite canonical algebra basis.
    #[error("AnnihilatorOfModule requires a finite-dimensional algebra")]
    InfiniteAnnihilatorAlgebra,
    /// QPA's finite canonical path basis starts with all vertex idempotents.
    #[error("the finite canonical algebra basis does not start with the quiver vertices")]
    MissingCanonicalVertexBasis,
    /// `ExtOverAlgebra` requires one common acting algebra.
    #[error("the two modules entered are not modules over the same algebra")]
    ExtModulesOverDifferentAlgebras,
    /// A Hom call cannot reject the already checked Ext modules.
    #[error("HomOverAlgebra rejected modules in the ExtOverAlgebra construction")]
    InternalExtHomRejected,
    /// GAP's Steinitz extension must complete the image basis to the Hom basis.
    #[error("NaturalHomomorphismBySubspace did not produce a complete Ext quotient basis")]
    IncompleteExtQuotientBasis,
    /// The coefficient closure only accepts raw vectors over its original field.
    #[error("the Ext representative map uses a different coefficient field")]
    ExtCoefficientFieldMismatch,
    /// The coefficient closure accepts only maps in `Hom(Omega(M), N)`.
    #[error("the entered map is outside the Hom space used by ExtOverAlgebra")]
    ExtMapOutsideHomSpace,
    /// Internal star Hom calls use modules over the already checked algebra.
    #[error("HomOverAlgebra rejected an internal StarOfModule Hom space")]
    InternalStarHomRejected,
    /// Every indecomposable projective has its distinguished vertex generator.
    #[error("an indecomposable projective has no generator in StarOfModule")]
    EmptyStarProjectiveGenerator,
    /// Arrow multiplication preserves one-vertex support.
    #[error("left multiplication by an arrow did not define a projective map")]
    StarArrowMapNotUniform,
    /// A nonzero star vertex has a raw physical Hom basis.
    #[error("a nonzero StarOfModule vertex has no coordinate basis")]
    MissingStarCoordinateBasis,
    /// Composition with left multiplication remains in the corresponding Hom space.
    #[error("a StarOfModule composite lies outside its Hom basis")]
    StarCompositionOutsideHomBasis,
    /// Source-generated star actions satisfy the opposite relations.
    #[error("the generated StarOfModule does not satisfy the opposite-algebra relations")]
    GeneratedStarViolatesRelations,
    /// The explicit transpose construction requires QPA's finite canonical basis.
    #[error("TransposeOfModule requires a finite-dimensional algebra")]
    InfiniteTransposeAlgebra,
    /// Every minimal module generator is supported at a quiver vertex.
    #[error("a minimal generator has no vertex support in TransposeOfModule")]
    MissingTransposeGeneratorSupport,
    /// The source's internal projective sums are nonempty in this branch.
    #[error("an internal TransposeOfModule projective direct sum returned fail")]
    InternalTransposeDirectSumFailed,
    /// Projective covers created by the transpose-on-maps construction retain
    /// the direct-sum metadata required by QPA's lifting method.
    #[error("an internal TransposeOfModuleHomomorphism projective lifting returned fail")]
    InternalTransposeHomomorphismLiftFailed,
    /// The lifted projective-cover square commutes by construction.
    #[error("the internal TransposeOfModuleHomomorphism kernel square does not commute")]
    InternalTransposeHomomorphismKernelSquareFailed,
    /// Applying the contravariant star functor preserves the presentation
    /// square used to induce the transpose map.
    #[error("the internal TransposeOfModuleHomomorphism cokernel square does not commute")]
    InternalTransposeHomomorphismCokernelSquareFailed,
    /// Minimal syzygy generators must have exactly one right-support vertex.
    #[error("a minimal syzygy generator is not uniform in TransposeOfModule")]
    NonuniformTransposeSyzygyGenerator,
    /// Every opposite projective summand has its distinguished generator.
    #[error("an opposite projective has no generator in TransposeOfModule")]
    EmptyTransposeProjectiveGenerator,
    /// `DTr(M)` remains over the identical acting algebra of `M`.
    #[error("HomOverAlgebra rejected the internal tau-rigidity Hom space")]
    InternalTauRigidHomRejected,
    /// `SimpleModules` accepts quotients only after QPA's admissibility filter
    /// has been established.
    #[error("SimpleModules requires an admissible quotient of a path algebra")]
    NonadmissibleSimpleModuleAlgebra,
    /// A vertex-simple representation must satisfy every admissible relation.
    #[error("a generated vertex-simple module does not satisfy the algebra relations")]
    GeneratedSimpleViolatesRelations,
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
    /// QPA explicitly rejects the zero module for its stable endomorphism
    /// quotient operation.
    #[error("zero module given as input for EndModuleProjOverAlgebra")]
    ZeroStableEndomorphismModule,
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
    /// Projective positions are one-based quiver-vertex positions.
    #[error("projective vertex {vertex} is outside 1..={vertex_count}")]
    ProjectiveVertexOutOfRange {
        /// Entered vertex position.
        vertex: usize,
        /// Number of vertices/projectives.
        vertex_count: usize,
    },
    /// QPA flattens physical module coordinates against a basis with a zero
    /// placeholder at every zero-dimensional vertex space.
    #[error("projective element has {actual} physical coordinates, expected {expected}")]
    ProjectiveElementPhysicalCoordinateCount {
        /// Entered flattened physical coordinate count.
        actual: usize,
        /// Flattened projective-basis count.
        expected: usize,
    },
    /// The two source methods require an acyclic path algebra or an admissible
    /// quotient, respectively.
    #[error("the entered algebra is not finite dimensional in QPA's projective-element domain")]
    ElementInIndecProjectiveAlgebraNotFinite,
    /// The algebra-element-to-projective conversion requires the identical
    /// path algebra or quotient family.
    #[error("the projective algebra element belongs to another algebra")]
    ForeignProjectiveAlgebraElement,
    /// QPA rejects zero in `ElementIn_vA_AsElementInIndecProj`.
    #[error("the entered projective algebra element is zero")]
    ZeroProjectiveAlgebraElement,
    /// A nonzero element must have exactly one left-support vertex.
    #[error("the projective algebra element is not left uniform")]
    ProjectiveAlgebraElementNotLeftUniform,
    /// QPA's coefficient-bearing successive tips must occur in the selected
    /// projective-basis component.
    #[error(
        "a coefficient-bearing tip is absent from projective basis component ({left_vertex}, {right_vertex})"
    )]
    ProjectiveLeadingTermMissing {
        /// Generating/left-support vertex.
        left_vertex: usize,
        /// Right-support component vertex.
        right_vertex: usize,
    },
    /// A distinguished generator of an indecomposable projective must have
    /// exactly one support vertex.
    #[error("a projective summand generator is not vertex uniform")]
    ProjectiveGeneratorNotVertexUniform,
    /// GAP's `IsMatrix` method domain excludes an empty element matrix here.
    #[error("the projective element matrix must be nonempty")]
    EmptyProjectiveElementMatrix,
    /// GAP matrices are rectangular homogeneous lists.
    #[error("the projective element matrix is ragged")]
    RaggedProjectiveElementMatrix,
    /// Matrix rows/columns and the two entered support lists must agree.
    #[error(
        "projective matrix shape {rows}x{columns} disagrees with {range_vertices} range and {source_vertices} source vertices"
    )]
    ProjectiveMatrixVertexListDimensions {
        /// Matrix row count.
        rows: usize,
        /// Matrix column count.
        columns: usize,
        /// Entered range-vertex count.
        range_vertices: usize,
        /// Entered source-vertex count.
        source_vertices: usize,
    },
    /// Every matrix entry must belong to the identical entered algebra.
    #[error("a projective matrix entry belongs to another algebra")]
    ForeignProjectiveMatrixElement,
    /// Every nonzero matrix entry must be both left and right uniform.
    #[error("a projective matrix entry is not uniform")]
    ProjectiveMatrixElementNotUniform,
    /// All nonzero entries in a row must have one left support.
    #[error("projective matrix row {row} is not homogeneous")]
    ProjectiveMatrixRowNotHomogeneous {
        /// One-based row position.
        row: usize,
    },
    /// All nonzero entries in a column must have one right support.
    #[error("projective matrix column {column} is not homogeneous")]
    ProjectiveMatrixColumnNotHomogeneous {
        /// One-based column position.
        column: usize,
    },
    /// Detected supports must equal the entered lists, including their order.
    #[error("the entered vertex lists disagree with projective matrix support")]
    ProjectiveMatrixSupportDiscrepancy,
    /// Nonempty source/range vertex lists make these internal direct sums
    /// infallible in QPA's method domain.
    #[error("an internal nonempty projective direct sum returned fail")]
    InternalProjectiveDirectSumFailed,
    /// A nonzero uniform matrix entry must produce a projective homomorphism.
    #[error("HomFromProjective rejected an internally uniform matrix entry")]
    InternalProjectiveMatrixHomRejected,
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
    /// Arithmetic in `OriginalPathAlgebra(A)` failed.
    #[error(transparent)]
    ParentPathAlgebra(#[from] PathAlgebraError),
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
        ArrowSpec, DenseMatrix, ModuleDimensionSeries, NamedArrowMatrix, NamedArrowSpec,
        PathAlgebra, PathAlgebraModule, PathAlgebraModuleElement, PathAlgebraModuleHomomorphism,
        PathAlgebraModuleMatrix, Quiver, QuiverAlgebra, QuiverAlgebraElement, VertexIndex,
        dimension_vector_partial_order, direct_sum_of_qpa_modules, end_over_algebra,
        from_end_m_to_hom_mm, from_hom_mm_to_end_m, identity_module_mapping,
        include_in_path_algebra, is_semisimple_module, is_simple_qpa_module, kernel_inclusion,
        lifting_inclusion_morphisms, minimal_generating_set_of_module, radical_series,
        right_fac_approximation, right_module_hom_over_algebra, right_module_over_path_algebra,
        socle_of_module, socle_of_module_inclusion, top_of_module, top_of_module_projection,
        trace_of_module, underlying_linear_map, zero_module, zero_module_mapping,
    };
    use serde_json::{Value, json};

    use super::{
        ExtOverAlgebraResult, algebra_as_module_of_enveloping_algebra,
        algebra_as_module_over_enveloping_algebra, annihilator_of_module, basis_of_projectives,
        canonical_basis_coefficients, d_tr, d_tr_homomorphism, d_tr_homomorphism_iterated,
        d_tr_iterated, dual_of_algebra_as_module_over_enveloping_algebra, dual_of_module,
        dual_of_module_homomorphism, element_in_indec_projective,
        element_in_v_a_as_element_in_indec_proj, end_modulo_projectives, ext_over_algebra,
        finite_canonical_basis, finite_qpa_default_basis, first_syzygy,
        from_matrix_to_homomorphism_of_projectives, hom_factoring_through_projectives,
        hom_from_projective, indec_injective_modules, indec_projective_modules, injective_envelope,
        is_injective_module, is_projective_module, is_rigid_module, is_tau_rigid_module,
        lifting_morphism_from_projective, matrix_of_homomorphism_between_projectives,
        nakayama_functor_of_module, nakayama_functor_of_module_homomorphism,
        opposite_nakayama_functor_of_module, opposite_nakayama_functor_of_module_homomorphism,
        projective_cover, simple_modules, socle_series, star_of_module,
        star_of_module_homomorphism, tr_d, tr_d_homomorphism, tr_d_homomorphism_iterated,
        tr_d_iterated, transpose_of_module, transpose_of_module_homomorphism,
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

    fn path_element_json<F: Field>(
        element: &qpa_core::PathAlgebraElement<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!(
            element
                .external_representation()
                .terms()
                .iter()
                .map(|(walk, scalar)| json!({
                    "walk": walk,
                    "coefficient": coefficient(scalar),
                }))
                .collect::<Vec<_>>()
        )
    }

    fn path_element_matrix_json<F: Field>(
        matrix: &[Vec<qpa_core::PathAlgebraElement<F>>],
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!(
            matrix
                .iter()
                .map(|row| row
                    .iter()
                    .map(|element| path_element_json(element, coefficient))
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>()
        )
    }

    fn physical_module_element_json<F: Field>(
        element: &PathAlgebraModuleElement<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "dimension": element.module().dimension_vector(),
            "support": element.support_vertices(),
            "coordinates": element
                .flat_external_representation()
                .iter()
                .map(coefficient)
                .collect::<Vec<_>>(),
        })
    }

    fn entered_path_matrix_equals<F: Field>(
        returned: &[Vec<qpa_core::PathAlgebraElement<F>>],
        entered: &[Vec<QuiverAlgebraElement<F>>],
    ) -> bool {
        returned.len() == entered.len()
            && returned
                .iter()
                .zip(entered)
                .all(|(returned_row, entered_row)| {
                    returned_row.len() == entered_row.len()
                        && returned_row.iter().zip(entered_row).all(
                            |(returned_element, entered_element)| {
                                matches!(entered_element, QuiverAlgebraElement::Path(element)
                                if returned_element == element)
                            },
                        )
                })
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

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_projective_matrix_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let path_algebra = PathAlgebra::new(field.clone(), quiver).unwrap();
        let algebra = QuiverAlgebra::Path(path_algebra);
        let vertices = algebra.quiver().vertices().collect::<Vec<_>>();
        let arrow = algebra.quiver().arrows().next().unwrap();
        let v1 = include_in_path_algebra(&vertices[0].path(), &algebra).unwrap();
        let v2 = include_in_path_algebra(&vertices[1].path(), &algebra).unwrap();
        let a = include_in_path_algebra(&arrow.path(), &algebra).unwrap();
        let entered_matrix = vec![vec![v1.clone(), a.clone()], vec![algebra.zero(), v2]];
        let map =
            from_matrix_to_homomorphism_of_projectives(&algebra, &entered_matrix, &[1, 2], &[1, 2])
                .unwrap();
        let returned_matrix = matrix_of_homomorphism_between_projectives(&map).unwrap();
        let uniform = v1.add(&a).unwrap();
        let projective_element =
            element_in_v_a_as_element_in_indec_proj(&algebra, &uniform).unwrap();
        let recovered = element_in_indec_projective(&algebra, &projective_element, 1).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let sum = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
        let zero_module = right_module_over_path_algebra(algebra.clone(), vec![0, 0], &[])
            .unwrap()
            .unwrap();
        let zero_to_sum = matrix_of_homomorphism_between_projectives(
            &zero_module_mapping(zero_module.clone(), sum.clone()).unwrap(),
        )
        .unwrap();
        let sum_to_zero = matrix_of_homomorphism_between_projectives(
            &zero_module_mapping(sum, zero_module.clone()).unwrap(),
        )
        .unwrap();
        let zero_to_zero = matrix_of_homomorphism_between_projectives(
            &zero_module_mapping(zero_module.clone(), zero_module).unwrap(),
        )
        .unwrap();
        let single =
            from_matrix_to_homomorphism_of_projectives(&algebra, &[vec![a]], &[2], &[1]).unwrap();

        let quotient_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(3), "y"),
            ],
        )
        .unwrap();
        let quotient_parent = PathAlgebra::new(field, quotient_quiver).unwrap();
        let quotient_arrows = quotient_parent.quiver().arrows().collect::<Vec<_>>();
        let x = quotient_parent.element(quotient_arrows[0].path()).unwrap();
        let y = quotient_parent.element(quotient_arrows[1].path()).unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(quotient_parent.clone(), &[x.multiply(&y).unwrap()]).unwrap(),
        );
        let quotient_vertices = quotient.quiver().vertices().collect::<Vec<_>>();
        let qv1 = include_in_path_algebra(&quotient_vertices[0].path(), &quotient).unwrap();
        let qv2 = include_in_path_algebra(&quotient_vertices[1].path(), &quotient).unwrap();
        let qx = include_in_path_algebra(&quotient_arrows[0].path(), &quotient).unwrap();
        let qy = include_in_path_algebra(&quotient_arrows[1].path(), &quotient).unwrap();
        let quotient_matrix = vec![vec![qx.clone(), quotient.zero()], vec![qv2, qy]];
        let quotient_map = from_matrix_to_homomorphism_of_projectives(
            &quotient,
            &quotient_matrix,
            &[2, 3],
            &[1, 2],
        )
        .unwrap();
        let quotient_back = matrix_of_homomorphism_between_projectives(&quotient_map).unwrap();
        let quotient_uniform = qv1.add(&qx).unwrap();
        let quotient_projective_element =
            element_in_v_a_as_element_in_indec_proj(&quotient, &quotient_uniform).unwrap();
        let quotient_recovered =
            element_in_indec_projective(&quotient, &quotient_projective_element, 1).unwrap();

        json!({
            "path": {
                "map": module_map_json(&map, coefficient),
                "matrix": path_element_matrix_json(&returned_matrix, coefficient),
                "roundtrip_equal": entered_path_matrix_equals(
                    &returned_matrix,
                    &entered_matrix,
                ),
                "underlying": matrix_json(&underlying_linear_map(&map).unwrap(), coefficient),
                "uniform_projective_element": physical_module_element_json(
                    &projective_element,
                    coefficient,
                ),
                "uniform_recovered": path_element_json(&recovered, coefficient),
                "single": module_map_json(&single, coefficient),
                "single_underlying": matrix_json(
                    &underlying_linear_map(&single).unwrap(),
                    coefficient,
                ),
                "zero_to_sum": path_element_matrix_json(&zero_to_sum, coefficient),
                "sum_to_zero": path_element_matrix_json(&sum_to_zero, coefficient),
                "zero_to_zero": path_element_matrix_json(&zero_to_zero, coefficient),
            },
            "quotient": {
                "map": module_map_json(&quotient_map, coefficient),
                "matrix_parent_representatives": path_element_matrix_json(
                    &quotient_back,
                    coefficient,
                ),
                "underlying": matrix_json(
                    &underlying_linear_map(&quotient_map).unwrap(),
                    coefficient,
                ),
                "uniform_projective_element": physical_module_element_json(
                    &quotient_projective_element,
                    coefficient,
                ),
                "uniform_recovered_parent": path_element_json(
                    &quotient_recovered,
                    coefficient,
                ),
            },
        })
    }

    fn attribute_module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "dimension": module.dimension_vector(),
            "matrices": module
                .matrices_of_path_algebra_module()
                .iter()
                .map(|entry| match entry {
                    PathAlgebraModuleMatrix::ZeroScalar => json!("0"),
                    PathAlgebraModuleMatrix::Matrix(matrix) => matrix_json(matrix, coefficient),
                })
                .collect::<Vec<_>>(),
        })
    }

    fn attribute_module_list_json<F: Field>(
        modules: &[PathAlgebraModule<F>],
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Vec<Value> {
        modules
            .iter()
            .map(|module| attribute_module_json(module, coefficient))
            .collect()
    }

    fn star_module_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module.cached_star_of_module().is_some();
        let value = star_of_module(module).unwrap();
        let has_after = module.cached_star_of_module().is_some();
        let repeated = star_of_module(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn nakayama_module_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module.cached_nakayama_functor_of_module().is_some();
        let value = nakayama_functor_of_module(module).unwrap();
        let has_after = module.cached_nakayama_functor_of_module().is_some();
        let repeated = nakayama_functor_of_module(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn opposite_nakayama_module_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module
            .cached_opposite_nakayama_functor_of_module()
            .is_some();
        let value = opposite_nakayama_functor_of_module(module).unwrap();
        let has_after = module
            .cached_opposite_nakayama_functor_of_module()
            .is_some();
        let repeated = opposite_nakayama_functor_of_module(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn module_functor_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let star = star_module_attribute_json(module, coefficient);
        let nakayama = nakayama_module_attribute_json(module, coefficient);
        let opposite_nakayama = opposite_nakayama_module_attribute_json(module, coefficient);
        json!({
            "star": star,
            "nakayama": nakayama,
            "opposite_nakayama": opposite_nakayama,
        })
    }

    fn star_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism.cached_star_homomorphism().is_some();
        let value = star_of_module_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism.cached_star_homomorphism().is_some();
        let repeated = star_of_module_homomorphism(homomorphism).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn nakayama_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism
            .cached_nakayama_functor_homomorphism()
            .is_some();
        let value = nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism
            .cached_nakayama_functor_homomorphism()
            .is_some();
        let repeated = nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn opposite_nakayama_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism
            .cached_opposite_nakayama_functor_homomorphism()
            .is_some();
        let value = opposite_nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism
            .cached_opposite_nakayama_functor_homomorphism()
            .is_some();
        let repeated = opposite_nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn map_functor_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let star = star_map_attribute_json(homomorphism, coefficient);
        let nakayama = nakayama_map_attribute_json(homomorphism, coefficient);
        let opposite_nakayama = opposite_nakayama_map_attribute_json(homomorphism, coefficient);
        json!({
            "star": star,
            "nakayama": nakayama,
            "opposite_nakayama": opposite_nakayama,
        })
    }

    fn transpose_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module.cached_transpose_of_module().is_some();
        let value = transpose_of_module(module).unwrap();
        let has_after = module.cached_transpose_of_module().is_some();
        let repeated = transpose_of_module(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
            "indecomposable_bound": value.cached_is_indecomposable_module().is_some(),
            "indecomposable_value": value.cached_is_indecomposable_module(),
        })
    }

    fn d_tr_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module.cached_d_tr().is_some();
        let value = d_tr(module).unwrap();
        let has_after = module.cached_d_tr().is_some();
        let repeated = d_tr(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn tr_d_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module.cached_tr_d().is_some();
        let value = tr_d(module).unwrap();
        let has_after = module.cached_tr_d().is_some();
        let repeated = tr_d(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn transpose_functor_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let transpose = transpose_attribute_json(module, coefficient);
        let dtr = d_tr_attribute_json(module, coefficient);
        let trd = tr_d_attribute_json(module, coefficient);
        let tau_has_before = module.cached_is_tau_rigid_module().is_some();
        let tau_value = is_tau_rigid_module(module).unwrap();
        let tau_has_after = module.cached_is_tau_rigid_module().is_some();
        json!({
            "transpose": transpose,
            "dtr": dtr,
            "trd": trd,
            "tau_rigid": {
                "has_before": tau_has_before,
                "value": tau_value,
                "has_after": tau_has_after,
            },
        })
    }

    fn functor_iterations_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let d_zero = d_tr_iterated(module, 0).unwrap();
        let d_one = d_tr_iterated(module, 1).unwrap();
        let d_one_generators = minimal_generating_set_of_module(d_one.clone()).unwrap();
        let d_minus_one = d_tr_iterated(module, -1).unwrap();
        let d_two = d_tr_iterated(module, 2).unwrap();
        let t_zero = tr_d_iterated(module, 0).unwrap();
        let t_one = tr_d_iterated(module, 1).unwrap();
        let t_minus_one = tr_d_iterated(module, -1).unwrap();
        let t_two = tr_d_iterated(module, 2).unwrap();
        json!({
            "dtr_zero_identity": d_zero.is_identical(module),
            "dtr_one_identity": d_one.is_identical(&d_tr(module).unwrap()),
            "dtr_one_generators": d_one_generators.iter().map(|element|
                element.external_representation().iter().map(|row|
                    row.iter().map(coefficient).collect::<Vec<_>>()).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
            "dtr_minus_one_identity": d_minus_one.is_identical(&tr_d(module).unwrap()),
            "dtr_two": attribute_module_json(&d_two, coefficient),
            "trd_zero_identity": t_zero.is_identical(module),
            "trd_one_identity": t_one.is_identical(&tr_d(module).unwrap()),
            "trd_minus_one_identity": t_minus_one.is_identical(&d_tr(module).unwrap()),
            "trd_two": attribute_module_json(&t_two, coefficient),
        })
    }

    fn transpose_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism.cached_transpose_homomorphism().is_some();
        let value = transpose_of_module_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism.cached_transpose_homomorphism().is_some();
        let repeated = transpose_of_module_homomorphism(homomorphism).unwrap();
        let source_attribute = transpose_of_module(homomorphism.range()).unwrap();
        let range_attribute = transpose_of_module(homomorphism.source()).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
            "source_attribute_identity": value.source().is_identical(&source_attribute),
            "range_attribute_identity": value.range().is_identical(&range_attribute),
        })
    }

    fn d_tr_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism.cached_d_tr_homomorphism().is_some();
        let value = d_tr_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism.cached_d_tr_homomorphism().is_some();
        let repeated = d_tr_homomorphism(homomorphism).unwrap();
        let source_attribute = d_tr(homomorphism.source()).unwrap();
        let range_attribute = d_tr(homomorphism.range()).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
            "source_attribute_identity": value.source().is_identical(&source_attribute),
            "range_attribute_identity": value.range().is_identical(&range_attribute),
        })
    }

    fn tr_d_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism.cached_tr_d_homomorphism().is_some();
        let value = tr_d_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism.cached_tr_d_homomorphism().is_some();
        let repeated = tr_d_homomorphism(homomorphism).unwrap();
        let source_attribute = tr_d(homomorphism.source()).unwrap();
        let range_attribute = tr_d(homomorphism.range()).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
            "source_attribute_identity": value.source().is_identical(&source_attribute),
            "range_attribute_identity": value.range().is_identical(&range_attribute),
        })
    }

    fn ar_map_functor_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "entered": module_map_json(homomorphism, coefficient),
            "transpose": transpose_map_attribute_json(homomorphism, coefficient),
            "dtr": d_tr_map_attribute_json(homomorphism, coefficient),
            "trd": tr_d_map_attribute_json(homomorphism, coefficient),
        })
    }

    fn ar_map_iterations_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let d_zero = d_tr_homomorphism_iterated(homomorphism, 0).unwrap();
        let d_one = d_tr_homomorphism_iterated(homomorphism, 1).unwrap();
        let d_minus_one = d_tr_homomorphism_iterated(homomorphism, -1).unwrap();
        let d_two = d_tr_homomorphism_iterated(homomorphism, 2).unwrap();
        let t_zero = tr_d_homomorphism_iterated(homomorphism, 0).unwrap();
        let t_one = tr_d_homomorphism_iterated(homomorphism, 1).unwrap();
        let t_minus_one = tr_d_homomorphism_iterated(homomorphism, -1).unwrap();
        let t_two = tr_d_homomorphism_iterated(homomorphism, 2).unwrap();
        json!({
            "dtr_zero_identity": d_zero.is_identical(homomorphism),
            "dtr_one_identity": d_one.is_identical(&d_tr_homomorphism(homomorphism).unwrap()),
            "dtr_minus_one_identity": d_minus_one
                .is_identical(&tr_d_homomorphism(homomorphism).unwrap()),
            "dtr_two": module_map_json(&d_two, coefficient),
            "trd_zero_identity": t_zero.is_identical(homomorphism),
            "trd_one_identity": t_one.is_identical(&tr_d_homomorphism(homomorphism).unwrap()),
            "trd_minus_one_identity": t_minus_one
                .is_identical(&d_tr_homomorphism(homomorphism).unwrap()),
            "trd_two": module_map_json(&t_two, coefficient),
        })
    }

    fn property_bound<F: Field>(module: &PathAlgebraModule<F>) -> Vec<bool> {
        vec![
            module.cached_is_injective_module().is_some(),
            module.cached_is_simple_qpa_module().is_some(),
            module.cached_is_semisimple_module().is_some(),
        ]
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn module_families_properties_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let injectives = indec_injective_modules(&algebra).unwrap();
        let injectives_again = indec_injective_modules(&algebra).unwrap();
        let simples = simple_modules(&algebra).unwrap();
        let simples_again = simple_modules(&algebra).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let zero = zero_module(&algebra).unwrap();
        let semisimple_sum = direct_sum_of_qpa_modules(&simples).unwrap().unwrap();
        let injective_installed = injectives
            .iter()
            .map(|module| {
                vec![
                    module.cached_is_injective_module() == Some(true),
                    module.cached_is_indecomposable_module() == Some(true),
                    module.cached_is_projective_module().is_some(),
                ]
            })
            .collect::<Vec<_>>();
        let simple_installed = simples
            .iter()
            .map(|module| {
                vec![
                    module.cached_is_simple_qpa_module() == Some(true),
                    module.cached_is_indecomposable_module() == Some(true),
                    module.cached_is_semisimple_module() == Some(true),
                ]
            })
            .collect::<Vec<_>>();
        let injective_values = injectives
            .iter()
            .chain(&simples)
            .chain(&projectives)
            .chain(std::iter::once(&zero))
            .map(|module| is_injective_module(module).unwrap())
            .collect::<Vec<_>>();
        let simple_values = simples
            .iter()
            .chain(&projectives)
            .chain([&zero, &semisimple_sum])
            .map(is_simple_qpa_module)
            .collect::<Vec<_>>();
        let semisimple_values = simples
            .iter()
            .chain(&projectives)
            .chain([&zero, &semisimple_sum])
            .map(|module| is_semisimple_module(module).unwrap())
            .collect::<Vec<_>>();
        let properties_bound = injectives
            .iter()
            .chain(&simples)
            .chain(&projectives)
            .chain([&zero, &semisimple_sum])
            .map(property_bound)
            .collect::<Vec<_>>();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field, loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_injectives = indec_injective_modules(&quotient).unwrap();
        let quotient_injectives_again = indec_injective_modules(&quotient).unwrap();
        let quotient_simples = simple_modules(&quotient).unwrap();
        let quotient_simples_again = simple_modules(&quotient).unwrap();
        let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
        let quotient_zero = zero_module(&quotient).unwrap();
        let quotient_injective_installed = quotient_injectives
            .iter()
            .map(|module| {
                vec![
                    module.cached_is_injective_module() == Some(true),
                    module.cached_is_indecomposable_module() == Some(true),
                    module.cached_is_projective_module().is_some(),
                ]
            })
            .collect::<Vec<_>>();
        let quotient_simple_installed = quotient_simples
            .iter()
            .map(|module| {
                vec![
                    module.cached_is_simple_qpa_module() == Some(true),
                    module.cached_is_indecomposable_module() == Some(true),
                    module.cached_is_semisimple_module() == Some(true),
                ]
            })
            .collect::<Vec<_>>();
        let quotient_injective_values = quotient_injectives
            .iter()
            .chain(&quotient_simples)
            .chain([&quotient_projective, &quotient_zero])
            .map(|module| is_injective_module(module).unwrap())
            .collect::<Vec<_>>();
        let quotient_simple_values = quotient_simples
            .iter()
            .chain([&quotient_projective, &quotient_zero])
            .map(is_simple_qpa_module)
            .collect::<Vec<_>>();
        let quotient_semisimple_values = quotient_simples
            .iter()
            .chain([&quotient_projective, &quotient_zero])
            .map(|module| is_semisimple_module(module).unwrap())
            .collect::<Vec<_>>();
        let quotient_properties_bound = quotient_injectives
            .iter()
            .chain(&quotient_simples)
            .chain([&quotient_projective, &quotient_zero])
            .map(property_bound)
            .collect::<Vec<_>>();

        json!({
            "path": {
                "injectives": attribute_module_list_json(&injectives, coefficient),
                "injective_cache_identity": injectives.iter().zip(&injectives_again)
                    .map(|(first, second)| first.is_identical(second)).collect::<Vec<_>>(),
                "injective_installed": injective_installed,
                "simples": attribute_module_list_json(&simples, coefficient),
                "simple_cache_identity": simples.iter().zip(&simples_again)
                    .map(|(first, second)| first.is_identical(second)).collect::<Vec<_>>(),
                "simple_installed": simple_installed,
                "injective_values": injective_values,
                "simple_values": simple_values,
                "semisimple_values": semisimple_values,
                "property_bound": properties_bound,
            },
            "quotient": {
                "injectives": attribute_module_list_json(&quotient_injectives, coefficient),
                "injective_cache_identity": quotient_injectives.iter()
                    .zip(&quotient_injectives_again).map(|(first, second)|
                        first.is_identical(second)).collect::<Vec<_>>(),
                "injective_installed": quotient_injective_installed,
                "simples": attribute_module_list_json(&quotient_simples, coefficient),
                "simple_cache_identity": quotient_simples.iter().zip(&quotient_simples_again)
                    .map(|(first, second)| first.is_identical(second)).collect::<Vec<_>>(),
                "simple_installed": quotient_simple_installed,
                "injective_values": quotient_injective_values,
                "simple_values": quotient_simple_values,
                "semisimple_values": quotient_semisimple_values,
                "property_bound": quotient_properties_bound,
            },
        })
    }

    fn series_json(series: ModuleDimensionSeries) -> Value {
        match series {
            ModuleDimensionSeries::ZeroDimensionVector(vector) => json!(vector),
            ModuleDimensionSeries::Layers(layers) => json!(layers),
        }
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn module_series_case<F: Field>(field: F) -> Value {
        let quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let projectives = indec_projective_modules(&algebra).unwrap();
        let path_simple = simple_modules(&algebra).unwrap()[0].clone();
        let zero = zero_module(&algebra).unwrap();
        let semisimple = right_module_over_path_algebra(algebra, vec![1, 1, 1], &[])
            .unwrap()
            .unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field, loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let x_squared = x.multiply(&x).unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x_squared.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
        let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
        let quotient_zero = zero_module(&quotient).unwrap();

        json!({
            "path": {
                "dimensions": projectives.iter().map(PathAlgebraModule::dimension)
                    .chain([zero.dimension(), semisimple.dimension()]).collect::<Vec<_>>(),
                "partial_order": [
                    dimension_vector_partial_order(&projectives[2], &projectives[0]),
                    dimension_vector_partial_order(&projectives[0], &projectives[1]),
                    dimension_vector_partial_order(&projectives[0], &projectives[0]),
                    dimension_vector_partial_order(&projectives[1], &projectives[0]),
                ],
                "radical_series": projectives.iter().chain(std::iter::once(&semisimple))
                    .map(|module| series_json(radical_series(module).unwrap()))
                    .collect::<Vec<_>>(),
                "socle_series": projectives.iter().chain(std::iter::once(&semisimple))
                    .map(|module| series_json(socle_series(module).unwrap()))
                    .collect::<Vec<_>>(),
                "zero_radical_series": series_json(radical_series(&zero).unwrap()),
                "zero_socle_series": series_json(socle_series(&zero).unwrap()),
            },
            "quotient": {
                "dimensions": [quotient_projective.dimension(), quotient_simple.dimension(),
                    quotient_zero.dimension()],
                "partial_order": [
                    dimension_vector_partial_order(&quotient_simple, &quotient_projective),
                    dimension_vector_partial_order(&quotient_projective, &quotient_simple),
                ],
                "cross_length_partial_order": [
                    dimension_vector_partial_order(&quotient_simple, &path_simple),
                    dimension_vector_partial_order(&path_simple, &quotient_simple),
                    dimension_vector_partial_order(&quotient_simple, &projectives[0]),
                    dimension_vector_partial_order(&projectives[0], &quotient_simple),
                    dimension_vector_partial_order(&quotient_simple, &projectives[1]),
                    dimension_vector_partial_order(&projectives[1], &quotient_simple),
                ],
                "radical_projective": series_json(radical_series(&quotient_projective).unwrap()),
                "socle_projective": series_json(socle_series(&quotient_projective).unwrap()),
                "radical_simple": series_json(radical_series(&quotient_simple).unwrap()),
                "socle_simple": series_json(socle_series(&quotient_simple).unwrap()),
                "zero_radical_series": series_json(radical_series(&quotient_zero).unwrap()),
                "zero_socle_series": series_json(socle_series(&quotient_zero).unwrap()),
            },
        })
    }

    fn annihilator_json<F: Field>(
        algebra: &QuiverAlgebra<F>,
        elements: &[QuiverAlgebraElement<F>],
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let canonical_basis = finite_canonical_basis(algebra).unwrap().unwrap();
        let default_basis = finite_qpa_default_basis(algebra, &canonical_basis).unwrap();
        let basis_change = DenseMatrix::from_rows(
            default_basis
                .iter()
                .map(|element| {
                    canonical_basis_coefficients(algebra, &canonical_basis, element).unwrap()
                })
                .collect(),
            algebra.field().clone(),
        )
        .unwrap();
        json!(
            elements
                .iter()
                .map(|element| {
                    let canonical_coordinates =
                        canonical_basis_coefficients(algebra, &canonical_basis, element).unwrap();
                    basis_change
                        .solution_row(&canonical_coordinates)
                        .unwrap()
                        .unwrap()
                        .iter()
                        .map(coefficient)
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        )
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn module_annihilator_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let projectives = indec_projective_modules(&algebra).unwrap();
        let faithful = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
        let simples = simple_modules(&algebra).unwrap();
        let zero = zero_module(&algebra).unwrap();
        let arrow_zero = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
            .unwrap()
            .unwrap();
        let zero_annihilator = annihilator_of_module(&zero).unwrap();
        let simple1_annihilator = annihilator_of_module(&simples[0]).unwrap();
        let simple2_annihilator = annihilator_of_module(&simples[1]).unwrap();
        let arrow_zero_annihilator = annihilator_of_module(&arrow_zero).unwrap();
        let faithful_annihilator = annihilator_of_module(&faithful).unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
        let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
        let quotient_zero = zero_module(&quotient).unwrap();
        let quotient_two_zero = right_module_over_path_algebra(quotient.clone(), vec![2], &[])
            .unwrap()
            .unwrap();
        let quotient_projective_annihilator = annihilator_of_module(&quotient_projective).unwrap();
        let quotient_simple_annihilator = annihilator_of_module(&quotient_simple).unwrap();
        let quotient_zero_annihilator = annihilator_of_module(&quotient_zero).unwrap();
        let quotient_two_zero_annihilator = annihilator_of_module(&quotient_two_zero).unwrap();

        let bound_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "c"),
                ArrowSpec::named(vertex(2), vertex(3), "d"),
            ],
        )
        .unwrap();
        let bound_parent = PathAlgebra::new(field, bound_quiver).unwrap();
        let bound_arrows = bound_parent.quiver().arrows().collect::<Vec<_>>();
        let c = bound_parent.element(bound_arrows[0].path()).unwrap();
        let d = bound_parent.element(bound_arrows[1].path()).unwrap();
        let bound_quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(bound_parent, &[c.multiply(&d).unwrap()]).unwrap(),
        );
        let bound_zero = zero_module(&bound_quotient).unwrap();
        let bound_simple = simple_modules(&bound_quotient).unwrap()[2].clone();
        let bound_zero_annihilator = annihilator_of_module(&bound_zero).unwrap();
        let bound_simple_annihilator = annihilator_of_module(&bound_simple).unwrap();

        json!({
            "path": {
                "basis_dimension": finite_canonical_basis(&algebra).unwrap().unwrap().len(),
                "zero": annihilator_json(&algebra, &zero_annihilator, coefficient),
                "simple1": annihilator_json(&algebra, &simple1_annihilator, coefficient),
                "simple2": annihilator_json(&algebra, &simple2_annihilator, coefficient),
                "arrow_zero": annihilator_json(&algebra, &arrow_zero_annihilator, coefficient),
                "faithful": annihilator_json(&algebra, &faithful_annihilator, coefficient),
                "cache_bound": [
                    zero.cached_annihilator_of_module().is_some(),
                    simples[0].cached_annihilator_of_module().is_some(),
                    simples[1].cached_annihilator_of_module().is_some(),
                    arrow_zero.cached_annihilator_of_module().is_some(),
                    faithful.cached_annihilator_of_module().is_some(),
                ],
            },
            "quotient": {
                "basis_dimension": finite_canonical_basis(&quotient).unwrap().unwrap().len(),
                "projective": annihilator_json(
                    &quotient,
                    &quotient_projective_annihilator,
                    coefficient,
                ),
                "simple": annihilator_json(
                    &quotient,
                    &quotient_simple_annihilator,
                    coefficient,
                ),
                "zero": annihilator_json(&quotient, &quotient_zero_annihilator, coefficient),
                "two_zero": annihilator_json(
                    &quotient,
                    &quotient_two_zero_annihilator,
                    coefficient,
                ),
                "cache_bound": [
                    quotient_projective.cached_annihilator_of_module().is_some(),
                    quotient_simple.cached_annihilator_of_module().is_some(),
                    quotient_zero.cached_annihilator_of_module().is_some(),
                    quotient_two_zero.cached_annihilator_of_module().is_some(),
                ],
            },
            "bound_quotient": {
                "basis_dimension": finite_canonical_basis(&bound_quotient)
                    .unwrap().unwrap().len(),
                "zero": annihilator_json(
                    &bound_quotient,
                    &bound_zero_annihilator,
                    coefficient,
                ),
                "simple3": annihilator_json(
                    &bound_quotient,
                    &bound_simple_annihilator,
                    coefficient,
                ),
                "cache_bound": [
                    bound_zero.cached_annihilator_of_module().is_some(),
                    bound_simple.cached_annihilator_of_module().is_some(),
                ],
            },
        })
    }

    fn ext_json<F: Field>(
        ext: &ExtOverAlgebraResult<F>,
        field: &F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let basis_coordinates = if ext.basis().is_empty() {
            Value::Null
        } else {
            json!(
                ext.basis()
                    .iter()
                    .map(|homomorphism| ext
                        .coefficients(homomorphism)
                        .unwrap()
                        .unwrap()
                        .iter()
                        .map(coefficient)
                        .collect::<Vec<_>>())
                    .collect::<Vec<_>>()
            )
        };
        let (entered_coefficients, combo_coordinates) = if ext.basis().is_empty() {
            (Value::Null, Value::Null)
        } else {
            let entered = (0..ext.dimension())
                .map(|index| field.nth(Integer::from(i64::try_from(index + 2).unwrap())))
                .collect::<Vec<_>>();
            let mut terms = ext
                .basis()
                .iter()
                .zip(&entered)
                .map(|(homomorphism, scalar)| homomorphism.scale(scalar).unwrap());
            let mut combo = terms.next().unwrap();
            for term in terms {
                combo = combo.add(&term).unwrap();
            }
            let coordinates = ext.coefficients(&combo).unwrap().unwrap();
            (
                json!(entered.iter().map(coefficient).collect::<Vec<_>>()),
                json!(coordinates.iter().map(coefficient).collect::<Vec<_>>()),
            )
        };
        json!({
            "syzygy_inclusion": module_map_json(ext.syzygy_inclusion(), coefficient),
            "dimension": ext.dimension(),
            "basis": ext.basis().iter()
                .map(|homomorphism| module_map_json(homomorphism, coefficient))
                .collect::<Vec<_>>(),
            "basis_coordinates": basis_coordinates,
            "entered_coefficients": entered_coefficients,
            "combo_coordinates": combo_coordinates,
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn functors_star_nakayama_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let simples = simple_modules(&algebra).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let zero = zero_module(&algebra).unwrap();
        let cover = projective_cover(&simples[0]).unwrap();
        let inclusion = kernel_inclusion(&cover).unwrap();

        let mixed = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![field.one()], vec![field.zero()]], field.clone())
                    .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let mixed_map = right_module_hom_over_algebra(
            mixed.clone(),
            mixed.clone(),
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![field.one(), field.zero()],
                        vec![field.zero(), field.zero()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::from_rows(vec![vec![field.one()]], field.clone()).unwrap(),
            ],
        )
        .unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
        let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
        let quotient_x = include_in_path_algebra(
            &quotient.quiver().arrows().next().unwrap().path(),
            &quotient,
        )
        .unwrap();
        let quotient_generator =
            minimal_generating_set_of_module(quotient_projective.clone()).unwrap()[0].clone();
        let radical_image = quotient_generator.act(&quotient_x).unwrap();
        let quotient_radical_map = hom_from_projective(&radical_image, &quotient_projective)
            .unwrap()
            .unwrap();

        json!({
            "modules": {
                "a2_s2": module_functor_attribute_json(&simples[1], coefficient),
                "a2_p1": module_functor_attribute_json(&projectives[0], coefficient),
                "a2_zero": module_functor_attribute_json(&zero, coefficient),
                "a2_mixed": module_functor_attribute_json(&mixed, coefficient),
                "dual_simple": module_functor_attribute_json(&quotient_simple, coefficient),
                "dual_projective": module_functor_attribute_json(
                    &quotient_projective,
                    coefficient,
                ),
            },
            "maps": {
                "a2_inclusion": map_functor_attribute_json(&inclusion, coefficient),
                "a2_endomorphism": map_functor_attribute_json(&mixed_map, coefficient),
                "dual_radical": map_functor_attribute_json(&quotient_radical_map, coefficient),
            },
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn functors_transpose_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let simples = simple_modules(&algebra).unwrap();
        let mixed = right_module_over_path_algebra(
            algebra,
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![field.one()], vec![field.zero()]], field.clone())
                    .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        let kronecker_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "p"),
                ArrowSpec::named(vertex(1), vertex(2), "q"),
            ],
        )
        .unwrap();
        let kronecker =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
        let kronecker_simple = simple_modules(&kronecker).unwrap()[0].clone();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();

        let bound_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let bound_parent = PathAlgebra::new(field, bound_quiver).unwrap();
        let bound_arrows = bound_parent.quiver().arrows().collect::<Vec<_>>();
        let a = bound_parent.element(bound_arrows[0].path()).unwrap();
        let b = bound_parent.element(bound_arrows[1].path()).unwrap();
        let bound = QuiverAlgebra::Quotient(
            quotient_by_relations(bound_parent, &[a.multiply(&b).unwrap()]).unwrap(),
        );
        let bound_simple = simple_modules(&bound).unwrap()[0].clone();

        json!({
            "a2_s1": transpose_functor_attribute_json(&simples[0], coefficient),
            "a2_s2_projective": transpose_functor_attribute_json(&simples[1], coefficient),
            "a2_mixed": transpose_functor_attribute_json(&mixed, coefficient),
            "kronecker_s1": transpose_functor_attribute_json(&kronecker_simple, coefficient),
            "dual_simple": transpose_functor_attribute_json(&dual_simple, coefficient),
            "bound_a3_s1": transpose_functor_attribute_json(&bound_simple, coefficient),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn functors_ar_maps_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let a2_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());
        let a2_module = right_module_over_path_algebra(a2, vec![2, 0], &[])
            .unwrap()
            .unwrap();
        let a2_map = right_module_hom_over_algebra(
            a2_module.clone(),
            a2_module,
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![field.one(), field.one()],
                        vec![field.zero(), field.one()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::zero(1, 1, field.clone()),
            ],
        )
        .unwrap();

        let kronecker_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "p"),
                ArrowSpec::named(vertex(1), vertex(2), "q"),
            ],
        )
        .unwrap();
        let kronecker =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
        let kronecker_module = right_module_over_path_algebra(kronecker, vec![2, 0], &[])
            .unwrap()
            .unwrap();
        let kronecker_map = right_module_hom_over_algebra(
            kronecker_module.clone(),
            kronecker_module,
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![field.one(), field.one()],
                        vec![field.zero(), field.one()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::zero(1, 1, field.clone()),
            ],
        )
        .unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();
        let dual_identity = identity_module_mapping(dual_simple).unwrap();

        let bound_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let bound_parent = PathAlgebra::new(field, bound_quiver).unwrap();
        let bound_arrows = bound_parent.quiver().arrows().collect::<Vec<_>>();
        let a = bound_parent.element(bound_arrows[0].path()).unwrap();
        let b = bound_parent.element(bound_arrows[1].path()).unwrap();
        let bound = QuiverAlgebra::Quotient(
            quotient_by_relations(bound_parent, &[a.multiply(&b).unwrap()]).unwrap(),
        );
        let bound_simple = simple_modules(&bound).unwrap()[0].clone();
        let bound_identity = identity_module_mapping(bound_simple).unwrap();

        json!({
            "a2_shear": ar_map_functor_attribute_json(&a2_map, coefficient),
            "kronecker_shear": ar_map_functor_attribute_json(&kronecker_map, coefficient),
            "dual_simple_identity": ar_map_functor_attribute_json(&dual_identity, coefficient),
            "bound_a3_identity": ar_map_functor_attribute_json(&bound_identity, coefficient),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn functors_ar_map_iterations_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let a2_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());
        let a2_module = right_module_over_path_algebra(a2, vec![2, 0], &[])
            .unwrap()
            .unwrap();
        let a2_map = right_module_hom_over_algebra(
            a2_module.clone(),
            a2_module,
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![field.one(), field.one()],
                        vec![field.zero(), field.one()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::zero(1, 1, field.clone()),
            ],
        )
        .unwrap();

        let kronecker_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "p"),
                ArrowSpec::named(vertex(1), vertex(2), "q"),
            ],
        )
        .unwrap();
        let kronecker =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
        let kronecker_module = right_module_over_path_algebra(kronecker, vec![2, 0], &[])
            .unwrap()
            .unwrap();
        let kronecker_map = right_module_hom_over_algebra(
            kronecker_module.clone(),
            kronecker_module,
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![field.one(), field.one()],
                        vec![field.zero(), field.one()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::zero(1, 1, field.clone()),
            ],
        )
        .unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_parent = PathAlgebra::new(field, loop_quiver).unwrap();
        let x = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();
        let dual_identity = identity_module_mapping(dual_simple).unwrap();

        json!({
            "a2_shear": ar_map_iterations_json(&a2_map, coefficient),
            "kronecker_shear": ar_map_iterations_json(&kronecker_map, coefficient),
            "dual_simple_identity": ar_map_iterations_json(&dual_identity, coefficient),
        })
    }

    #[allow(clippy::needless_pass_by_value)]
    fn functors_iterations_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let a2_simple = simple_modules(&algebra).unwrap()[0].clone();

        let kronecker_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "p"),
                ArrowSpec::named(vertex(1), vertex(2), "q"),
            ],
        )
        .unwrap();
        let kronecker =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
        let kronecker_simple = simple_modules(&kronecker).unwrap()[0].clone();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field, loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();

        json!({
            "a2_s1": functor_iterations_json(&a2_simple, coefficient),
            "kronecker_s1": functor_iterations_json(&kronecker_simple, coefficient),
            "dual_simple": functor_iterations_json(&dual_simple, coefficient),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_ext_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let simples = simple_modules(&algebra).unwrap();
        let zero = zero_module(&algebra).unwrap();
        let a2_12 = ext_over_algebra(&simples[0], &simples[1]).unwrap();
        let a2_11 = ext_over_algebra(&simples[0], &simples[0]).unwrap();
        let a2_22 = ext_over_algebra(&simples[1], &simples[1]).unwrap();
        let a2_02 = ext_over_algebra(&zero, &simples[1]).unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
        let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
        let dual_simple_ext = ext_over_algebra(&quotient_simple, &quotient_simple).unwrap();
        let dual_projective_ext = ext_over_algebra(&quotient_simple, &quotient_projective).unwrap();
        let quotient_sum =
            direct_sum_of_qpa_modules(&[quotient_simple.clone(), quotient_projective.clone()])
                .unwrap()
                .unwrap();
        let dual_sum = ext_over_algebra(&quotient_simple, &quotient_sum).unwrap();

        let kronecker_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "p"),
                ArrowSpec::named(vertex(1), vertex(2), "q"),
            ],
        )
        .unwrap();
        let kronecker =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
        let kronecker_simples = simple_modules(&kronecker).unwrap();
        let kronecker_ext = ext_over_algebra(&kronecker_simples[0], &kronecker_simples[1]).unwrap();

        let rigid_json = |module: &PathAlgebraModule<F>| {
            let has_before = module.cached_is_rigid_module().is_some();
            let value = is_rigid_module(module).unwrap();
            let has_after = module.cached_is_rigid_module().is_some();
            json!({
                "has_before": has_before,
                "value": value,
                "has_after": has_after,
            })
        };

        json!({
            "a2_s1_s2": ext_json(&a2_12, &field, coefficient),
            "a2_s1_s1": ext_json(&a2_11, &field, coefficient),
            "a2_s2_s2": ext_json(&a2_22, &field, coefficient),
            "a2_zero_s2": ext_json(&a2_02, &field, coefficient),
            "dual_simple_simple": ext_json(&dual_simple_ext, &field, coefficient),
            "dual_simple_projective": ext_json(&dual_projective_ext, &field, coefficient),
            "dual_simple_sum": ext_json(&dual_sum, &field, coefficient),
            "kronecker_s1_s2": ext_json(&kronecker_ext, &field, coefficient),
            "rigid": {
                "a2_s1": rigid_json(&simples[0]),
                "a2_s2": rigid_json(&simples[1]),
                "a2_zero": rigid_json(&zero),
                "dual_simple": rigid_json(&quotient_simple),
                "dual_projective": rigid_json(&quotient_projective),
                "dual_sum": rigid_json(&quotient_sum),
                "kronecker_s1": rigid_json(&kronecker_simples[0]),
                "kronecker_s2": rigid_json(&kronecker_simples[1]),
            },
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_syzygy_injective_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let zero = zero_module(&algebra).unwrap();
        let zero_again = zero_module(&algebra).unwrap();
        let simple1 = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple2 = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let syzygy1 = first_syzygy(&simple1).unwrap();
        let syzygy2 = first_syzygy(&simple2).unwrap();
        let syzygy0 = first_syzygy(&zero).unwrap();
        let envelope1 = injective_envelope(&simple1).unwrap();
        let envelope2 = injective_envelope(&simple2).unwrap();
        let envelope0 = injective_envelope(&zero).unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field, loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_zero = zero_module(&quotient).unwrap();
        let quotient_simple = right_module_over_path_algebra(quotient.clone(), vec![1], &[])
            .unwrap()
            .unwrap();
        let quotient_syzygy = first_syzygy(&quotient_simple).unwrap();
        let quotient_envelope = injective_envelope(&quotient_simple).unwrap();

        json!({
            "path": {
                "zero": attribute_module_json(&zero, coefficient),
                "zero_cache_identity": zero.is_identical(&zero_again),
                "syzygy_simple1": attribute_module_json(&syzygy1, coefficient),
                "syzygy_simple2": attribute_module_json(&syzygy2, coefficient),
                "syzygy_zero": attribute_module_json(&syzygy0, coefficient),
                "syzygy_cache_identity": [
                    syzygy1.is_identical(&first_syzygy(&simple1).unwrap()),
                    syzygy2.is_identical(&first_syzygy(&simple2).unwrap()),
                    syzygy0.is_identical(&zero),
                ],
                "injective_simple1": module_map_json(&envelope1, coefficient),
                "injective_simple2": module_map_json(&envelope2, coefficient),
                "injective_zero": module_map_json(&envelope0, coefficient),
                "injective_cache_identity": [
                    envelope1.is_identical(&injective_envelope(&simple1).unwrap()),
                    envelope2.is_identical(&injective_envelope(&simple2).unwrap()),
                    envelope0.is_identical(&injective_envelope(&zero).unwrap()),
                ],
                "injective_source_identity": [
                    envelope1.source().is_identical(&simple1),
                    envelope2.source().is_identical(&simple2),
                    envelope0.source().is_identical(&zero),
                ],
            },
            "quotient": {
                "zero": attribute_module_json(&quotient_zero, coefficient),
                "syzygy_simple": attribute_module_json(&quotient_syzygy, coefficient),
                "injective_simple": module_map_json(&quotient_envelope, coefficient),
                "cache_identity": [
                    quotient_zero.is_identical(&zero_module(&quotient).unwrap()),
                    quotient_syzygy.is_identical(&first_syzygy(&quotient_simple).unwrap()),
                    quotient_envelope
                        .is_identical(&injective_envelope(&quotient_simple).unwrap()),
                ],
            },
        })
    }

    fn stable_module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        print_zero_vertex_map: bool,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let natural = end_modulo_projectives(module).unwrap();
        let factoring = hom_factoring_through_projectives(module, module).unwrap();
        let factor_matrices = factoring
            .iter()
            .map(|map| from_hom_mm_to_end_m(map).unwrap())
            .collect::<Vec<_>>();
        let quotient_basis = natural.range().basis_elements();
        let images = natural
            .source()
            .basis_matrices()
            .iter()
            .map(|matrix| {
                natural
                    .image(matrix)
                    .unwrap()
                    .coordinates()
                    .iter()
                    .map(coefficient)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let one = natural.range().one().unwrap().map(|element| {
            element
                .coordinates()
                .iter()
                .map(coefficient)
                .collect::<Vec<_>>()
        });
        let multiplication = quotient_basis
            .iter()
            .map(|left| {
                quotient_basis
                    .iter()
                    .map(|right| {
                        natural
                            .range()
                            .multiply(left, right)
                            .unwrap()
                            .coordinates()
                            .iter()
                            .map(coefficient)
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let roundtrip_matrices = natural
            .source()
            .basis_matrices()
            .iter()
            .map(|matrix| {
                from_hom_mm_to_end_m(&from_end_m_to_hom_mm(module, matrix).unwrap()).unwrap()
            })
            .collect::<Vec<_>>();
        let mut value = json!({
            "end_basis": natural.source().basis_matrices().iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
            "factor_matrices": factor_matrices.iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
            "ideal_basis": natural.range().ideal_basis_matrices().iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
            "kernel_basis": natural.kernel_basis_matrices().iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
            "source_cache_identity": natural.source().is_identical(&end_over_algebra(module).unwrap()),
            "dimension": natural.range().dimension(),
            "images": images,
            "one": one,
            "multiplication": multiplication,
            "roundtrip_matrices": roundtrip_matrices.iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
        });
        if print_zero_vertex_map {
            value["zero_vertex_map"] = module_map_json(
                &from_end_m_to_hom_mm(module, &natural.source().basis_matrices()[0]).unwrap(),
                coefficient,
            );
        }
        value
    }

    #[allow(clippy::needless_pass_by_value)]
    fn modulehom_stable_case<F: Field>(
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
        let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(algebra, vec![0, 1], &[])
            .unwrap()
            .unwrap();

        let loop_quiver =
            Quiver::from_named(vec!["w".into()], vec![NamedArrowSpec::named("w", "w", "x")])
                .unwrap();
        let path_algebra = PathAlgebra::new(field, loop_quiver).unwrap();
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
        let dual_simple = right_module_over_path_algebra(dual_numbers.clone(), vec![1], &[])
            .unwrap()
            .unwrap();
        let dual_projective = indec_projective_modules(&dual_numbers).unwrap()[0].clone();

        json!({
            "a2_module": stable_module_json(&module, false, coefficient),
            "a2_simple_one": stable_module_json(&simple_one, true, coefficient),
            "a2_simple_two": stable_module_json(&simple_two, true, coefficient),
            "dual_simple": stable_module_json(&dual_simple, false, coefficient),
            "dual_projective": stable_module_json(&dual_projective, false, coefficient),
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
    fn pinned_qpa_modulehom_projective_matrix_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_projective_matrices.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_projective_matrix_case(Q, ToString::to_string);
        let binary = modulehom_projective_matrix_case(Z2, |entry| {
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
    fn pinned_qpa_modulehom_syzygy_injective_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_syzygy_injective.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_syzygy_injective_case(Q, ToString::to_string);
        let binary = modulehom_syzygy_injective_case(Z2, |entry| {
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
    fn pinned_qpa_module_families_properties_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/module_families_properties.schema-v1.json"
        ))
        .unwrap();
        let rational = module_families_properties_case(Q, ToString::to_string);
        let binary = module_families_properties_case(Z2, |entry| {
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
    fn pinned_qpa_module_series_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/module_series.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(module_series_case(Q), expected["rational"]);
        assert_eq!(module_series_case(Z2), expected["binary"]);
    }

    #[test]
    fn pinned_qpa_module_annihilator_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/module_annihilator.schema-v1.json"
        ))
        .unwrap();
        let rational = module_annihilator_case(Q, ToString::to_string);
        let binary = module_annihilator_case(Z2, |entry| {
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
    fn pinned_qpa_modulehom_ext_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_ext.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_ext_case(Q, ToString::to_string);
        let binary = modulehom_ext_case(Z2, |entry| {
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
    fn pinned_qpa_star_and_nakayama_functor_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/functors_star_nakayama.schema-v1.json"
        ))
        .unwrap();
        let rational = functors_star_nakayama_case(Q, ToString::to_string);
        let binary = functors_star_nakayama_case(Z2, |entry| {
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
    fn pinned_qpa_transpose_and_ar_functor_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/functors_transpose.schema-v1.json"
        ))
        .unwrap();
        let rational = functors_transpose_case(Q, ToString::to_string);
        let binary = functors_transpose_case(Z2, |entry| {
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
    fn pinned_qpa_ar_functors_on_homomorphisms_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/functors_ar_maps.schema-v1.json"
        ))
        .unwrap();
        let rational = functors_ar_maps_case(Q, ToString::to_string);
        let binary = functors_ar_maps_case(Z2, |entry| {
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
    fn pinned_qpa_ar_functor_homomorphism_iterations_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/functors_ar_map_iterations.schema-v1.json"
        ))
        .unwrap();
        let rational = functors_ar_map_iterations_case(Q, ToString::to_string);
        let binary = functors_ar_map_iterations_case(Z2, |entry| {
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
    fn pinned_qpa_ar_functor_iterations_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/functors_iterations.schema-v1.json"
        ))
        .unwrap();
        let rational = functors_iterations_case(Q, ToString::to_string);
        let binary = functors_iterations_case(Z2, |entry| {
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
    fn annihilator_of_module_reports_its_finite_basis_domain() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let module = right_module_over_path_algebra(algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        assert_eq!(
            annihilator_of_module(&module).unwrap_err(),
            super::EnvelopingModuleError::InfiniteAnnihilatorAlgebra
        );
    }

    #[test]
    fn ext_over_algebra_preserves_domain_and_coefficient_boundaries() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let simples = simple_modules(&algebra).unwrap();
        let ext = ext_over_algebra(&simples[0], &simples[1]).unwrap();
        let outside = identity_module_mapping(simples[0].clone()).unwrap();
        assert_eq!(
            ext.coefficients(&outside).unwrap_err(),
            super::EnvelopingModuleError::ExtMapOutsideHomSpace
        );
        let zero_ext = ext_over_algebra(&simples[0], &simples[0]).unwrap();
        assert_eq!(zero_ext.coefficients(&outside).unwrap(), None);

        let other_quiver = Quiver::from_indexed(1, Vec::new()).unwrap();
        let other_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, other_quiver).unwrap());
        let other_simple = simple_modules(&other_algebra).unwrap()[0].clone();
        assert_eq!(
            ext_over_algebra(&simples[0], &other_simple).unwrap_err(),
            super::EnvelopingModuleError::ExtModulesOverDifferentAlgebras
        );
    }

    #[test]
    fn simple_modules_preserves_qpa_source_domain_errors() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(Q, quiver).unwrap();
        assert_eq!(
            simple_modules(&QuiverAlgebra::Path(parent.clone())).unwrap_err(),
            super::EnvelopingModuleError::InfiniteSimpleModuleAlgebra
        );
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x]).unwrap());
        assert_eq!(
            simple_modules(&quotient).unwrap_err(),
            super::EnvelopingModuleError::NonadmissibleSimpleModuleAlgebra
        );
    }

    #[test]
    fn projective_matrix_conversions_report_qpa_domain_errors() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        assert_eq!(
            element_in_v_a_as_element_in_indec_proj(&algebra, &algebra.zero()),
            Err(super::EnvelopingModuleError::ZeroProjectiveAlgebraElement)
        );
        let v1 = include_in_path_algebra(
            &algebra.quiver().vertices().next().unwrap().path(),
            &algebra,
        )
        .unwrap();
        assert_eq!(
            from_matrix_to_homomorphism_of_projectives(&algebra, &[vec![v1]], &[2], &[1],),
            Err(super::EnvelopingModuleError::ProjectiveMatrixSupportDiscrepancy)
        );
        let scaled_arrow =
            include_in_path_algebra(&algebra.quiver().arrows().next().unwrap().path(), &algebra)
                .unwrap()
                .scale_left(&q(2))
                .unwrap();
        assert_eq!(
            element_in_v_a_as_element_in_indec_proj(&algebra, &scaled_arrow),
            Err(super::EnvelopingModuleError::ProjectiveLeadingTermMissing {
                left_vertex: 1,
                right_vertex: 2,
            })
        );
        assert_eq!(
            from_matrix_to_homomorphism_of_projectives(&algebra, &[], &[], &[]),
            Err(super::EnvelopingModuleError::EmptyProjectiveElementMatrix)
        );
    }

    #[test]
    fn pinned_qpa_modulehom_stable_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_stable.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_stable_case(Q, ToString::to_string);
        let binary = modulehom_stable_case(Z2, |entry| {
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
    fn stable_endomorphism_quotient_rejects_the_zero_module_like_qpa() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let zero = right_module_over_path_algebra(algebra, vec![0, 0], &[])
            .unwrap()
            .unwrap();

        assert_eq!(
            end_modulo_projectives(&zero).unwrap_err(),
            super::EnvelopingModuleError::ZeroStableEndomorphismModule
        );
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
