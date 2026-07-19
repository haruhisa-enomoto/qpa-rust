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

/// Translate QPA's cached `FromIdentityToDoubleStarHomomorphism(M)` attribute.
///
/// At a basis vector `m` in vertex `i`, the source constructs the evaluation
/// map `Hom_A(M,A) -> e_i A` sending `g` to `g(m)`, converts its projective
/// values through the opposite algebra, and expresses that map in QPA's
/// ordered basis of `Hom_{A^op}(M^*, e_i A^op)`. Those coefficients are the
/// corresponding row of the natural map `M -> M^{**}`.
///
/// # Errors
///
/// Requires finite projective and canonical bases and propagates opposite
/// algebra, Hom-space, module-action, exact-coordinate, and checked-map
/// failures.
#[allow(clippy::too_many_lines)]
pub fn from_identity_to_double_star_homomorphism<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_identity_to_double_star() {
        return Ok(cached);
    }
    let algebra = module.right_acting_algebra();
    let opposite = opposite_module_algebra(algebra)?;
    let field = algebra.field().clone();
    let projectives = indec_projective_modules(algebra)?;
    let vertex_count = projectives.len();
    let hom_spaces = projectives
        .iter()
        .map(|projective| {
            hom_over_algebra(module, projective)?
                .ok_or(EnvelopingModuleError::InternalDoubleStarHomRejected)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let _hom_coordinate_bases = module_hom_flat_coordinate_bases(&hom_spaces, &field)?;

    let star = star_of_module(module)?;
    let opposite_projectives = indec_projective_modules(&opposite)?;
    let opposite_hom_spaces = opposite_projectives
        .iter()
        .map(|projective| {
            hom_over_algebra(&star, projective)?
                .ok_or(EnvelopingModuleError::InternalDoubleStarHomRejected)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let opposite_hom_coordinate_bases =
        module_hom_flat_coordinate_bases(&opposite_hom_spaces, &field)?;
    let double_star = star_of_module(&star)?;

    let opposite_projective_bases = basis_of_projectives(&opposite)?
        .ok_or(EnvelopingModuleError::ProjectiveBasisUnavailable)?;
    let opposite_canonical_basis = finite_canonical_basis(&opposite)?
        .ok_or(EnvelopingModuleError::MissingDoubleStarCanonicalBasis)?;
    let opposite_projective_basis_positions = opposite_projective_bases
        .iter()
        .enumerate()
        .map(|(projective, vertex_bases)| {
            vertex_bases
                .iter()
                .enumerate()
                .map(|(vertex, basis)| {
                    basis
                        .iter()
                        .enumerate()
                        .map(|(position, element)| {
                            opposite_canonical_basis
                                .iter()
                                .position(|candidate| candidate == element)
                                .ok_or(
                                    EnvelopingModuleError::MissingDoubleStarProjectiveBasisElement {
                                        projective: projective + 1,
                                        vertex: vertex + 1,
                                        position: position + 1,
                                    },
                                )
                        })
                        .collect::<Result<Vec<_>, EnvelopingModuleError>>()
                })
                .collect::<Result<Vec<_>, EnvelopingModuleError>>()
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;

    let module_basis = module.basis_elements();
    let module_dimensions = module.dimension_vector();
    let double_star_dimensions = double_star.dimension_vector();
    let mut matrices = Vec::with_capacity(vertex_count);
    let mut module_basis_offset = 0_usize;
    for projective_vertex in 0..vertex_count {
        let source_dimension = module_dimensions[projective_vertex];
        let range_dimension = double_star_dimensions[projective_vertex];
        if source_dimension == 0 || range_dimension == 0 {
            matrices.push(DenseMatrix::zero(
                source_dimension.max(1),
                range_dimension.max(1),
                field.clone(),
            ));
            module_basis_offset += source_dimension;
            continue;
        }
        let coordinate_basis = opposite_hom_coordinate_bases[projective_vertex]
            .as_ref()
            .ok_or(EnvelopingModuleError::MissingDoubleStarHomCoordinateBasis {
                vertex: projective_vertex + 1,
            })?;
        let mut rows = Vec::with_capacity(source_dimension);
        for module_element in
            &module_basis[module_basis_offset..module_basis_offset + source_dimension]
        {
            let mut evaluation_matrices = Vec::with_capacity(vertex_count);
            for (target_vertex, hom_basis) in hom_spaces.iter().enumerate() {
                let projective_dimension =
                    opposite_projectives[projective_vertex].dimension_vector()[target_vertex];
                if hom_basis.is_empty() || projective_dimension == 0 {
                    evaluation_matrices.push(DenseMatrix::zero(
                        hom_basis.len().max(1),
                        projective_dimension.max(1),
                        field.clone(),
                    ));
                    continue;
                }
                let value_rows = hom_basis
                    .iter()
                    .map(|homomorphism| {
                        let image = homomorphism.image(module_element)?;
                        let parent_element =
                            element_in_indec_projective(algebra, &image, target_vertex + 1)?;
                        let opposite_element =
                            opposite_projective_path_element(algebra, &parent_element)?;
                        let canonical_coordinates = canonical_basis_coefficients(
                            &opposite,
                            &opposite_canonical_basis,
                            &opposite_element,
                        )?;
                        Ok(
                            opposite_projective_basis_positions[projective_vertex][target_vertex]
                                .iter()
                                .map(|&position| canonical_coordinates[position].clone())
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
                evaluation_matrices.push(DenseMatrix::from_rows(value_rows, field.clone())?);
            }
            let evaluation = right_module_hom_over_algebra(
                star.clone(),
                opposite_projectives[projective_vertex].clone(),
                evaluation_matrices,
            )?;
            rows.push(
                coordinate_basis
                    .solution_row(&evaluation.flat_vertex_matrix_entries())?
                    .ok_or(EnvelopingModuleError::DoubleStarEvaluationOutsideHomBasis {
                        vertex: projective_vertex + 1,
                    })?,
            );
        }
        matrices.push(DenseMatrix::from_rows(rows, field.clone())?);
        module_basis_offset += source_dimension;
    }
    let evaluation = right_module_hom_over_algebra(module.clone(), double_star, matrices)?;
    Ok(module.cache_identity_to_double_star(&evaluation))
}

fn module_hom_flat_coordinate_bases<F: Field>(
    hom_spaces: &[Vec<PathAlgebraModuleHomomorphism<F>>],
    field: &F,
) -> Result<Vec<Option<DenseMatrix<F>>>, EnvelopingModuleError> {
    hom_spaces
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
        .collect()
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

struct NaturalTensorQuotient<F: Field> {
    relation_basis: DenseMatrix<F>,
    basis_lifts: DenseMatrix<F>,
    projection: DenseMatrix<F>,
}

fn natural_tensor_quotient<F: Field>(
    relation_rows: &DenseMatrix<F>,
) -> Result<NaturalTensorQuotient<F>, EnvelopingModuleError> {
    let field = relation_rows.field().clone();
    let ambient_dimension = relation_rows.ncols();
    let relation_basis = relation_rows.semi_echelon_transformation().vectors;
    let mut mutable_basis = relation_basis.clone();
    let mut complement = Vec::new();
    for coordinate in 0..ambient_dimension {
        let mut standard = vec![field.zero(); ambient_dimension];
        standard[coordinate] = field.one();
        let mut entries = mutable_basis.as_flat().to_vec();
        entries.extend_from_slice(&standard);
        let candidate = DenseMatrix::from_flat(
            mutable_basis.nrows() + 1,
            ambient_dimension,
            entries,
            field.clone(),
        )?
        .semi_echelon_transformation()
        .vectors;
        if candidate.nrows() > mutable_basis.nrows() {
            complement.push(standard);
            mutable_basis = candidate;
        }
    }
    if relation_basis.nrows() + complement.len() != ambient_dimension {
        return Err(EnvelopingModuleError::IncompleteTensorProductQuotientBasis);
    }

    let quotient_dimension = complement.len();
    let basis_lifts = DenseMatrix::from_flat(
        quotient_dimension,
        ambient_dimension,
        complement.iter().flatten().cloned().collect(),
        field.clone(),
    )?;
    let mut full_basis_entries = relation_basis.as_flat().to_vec();
    full_basis_entries.extend(complement.iter().flatten().cloned());
    let full_basis = DenseMatrix::from_flat(
        ambient_dimension,
        ambient_dimension,
        full_basis_entries,
        field.clone(),
    )?;
    let mut projection_entries = Vec::with_capacity(ambient_dimension * quotient_dimension);
    for coordinate in 0..ambient_dimension {
        let mut standard = vec![field.zero(); ambient_dimension];
        standard[coordinate] = field.one();
        let coefficients = full_basis
            .solution_row(&standard)?
            .ok_or(EnvelopingModuleError::IncompleteTensorProductQuotientBasis)?;
        projection_entries.extend_from_slice(
            &coefficients[ambient_dimension - quotient_dimension..ambient_dimension],
        );
    }
    let projection = DenseMatrix::from_flat(
        ambient_dimension,
        quotient_dimension,
        projection_entries,
        field,
    )?;
    Ok(NaturalTensorQuotient {
        relation_basis,
        basis_lifts,
        projection,
    })
}

fn tensor_dimension_pairing(
    left: &[usize],
    right: &[usize],
) -> Result<usize, EnvelopingModuleError> {
    left.iter().zip(right).try_fold(0_usize, |total, (&a, &b)| {
        total
            .checked_add(
                a.checked_mul(b)
                    .ok_or(EnvelopingModuleError::TensorProductDimensionOverflow)?,
            )
            .ok_or(EnvelopingModuleError::TensorProductDimensionOverflow)
    })
}

pub(crate) fn opposite_module_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<QuiverAlgebra<F>, EnvelopingModuleError> {
    match algebra {
        QuiverAlgebra::Path(algebra) => Ok(QuiverAlgebra::Path(algebra.opposite()?)),
        QuiverAlgebra::Quotient(algebra) => {
            Ok(QuiverAlgebra::Quotient(opposite_path_algebra(algebra)?))
        }
    }
}

fn opposite_projective_path_element<F: Field>(
    algebra: &QuiverAlgebra<F>,
    element: &PathAlgebraElement<F>,
) -> Result<QuiverAlgebraElement<F>, EnvelopingModuleError> {
    let opposite = element.opposite()?;
    match algebra {
        QuiverAlgebra::Path(_) => Ok(QuiverAlgebraElement::Path(opposite)),
        QuiverAlgebra::Quotient(algebra) => Ok(QuiverAlgebraElement::Quotient(
            opposite_path_algebra(algebra)?.element(opposite)?,
        )),
    }
}

fn tensor_projective_vertex<F: Field>(
    projective: &PathAlgebraModule<F>,
) -> Result<usize, EnvelopingModuleError> {
    let generator = minimal_generating_set_of_module(projective.clone())?
        .into_iter()
        .next()
        .ok_or(EnvelopingModuleError::EmptyTensorProjectiveGenerator)?;
    let support = generator.support_vertices();
    let [vertex] = support.as_slice() else {
        return Err(EnvelopingModuleError::NonuniformTensorProjectiveGenerator);
    };
    Ok(*vertex)
}
