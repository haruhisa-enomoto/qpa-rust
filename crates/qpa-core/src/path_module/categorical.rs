/// Translate QPA's `DirectSumOfQPAModules` block-matrix construction.
///
/// Summands and every vertex block retain their entered order. The returned
/// module carries QPA's direct-sum filter together with lazily recoverable
/// projection and inclusion attributes. `Ok(None)` is QPA's `fail` for an
/// empty list.
///
/// # Errors
///
/// Rejects different acting algebras, reports dimension overflow, and
/// propagates checked module/map construction failures.
pub fn direct_sum_of_qpa_modules<F: Field>(
    modules: &[PathAlgebraModule<F>],
) -> Result<Option<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    let Some(first) = modules.first() else {
        return Ok(None);
    };
    if modules
        .iter()
        .any(|module| module.right_acting_algebra() != first.right_acting_algebra())
    {
        return Err(PathAlgebraModuleError::DirectSumDifferentAlgebras);
    }

    let vertex_count = first.dimension_vector().len();
    let mut dimensions = vec![0_usize; vertex_count];
    for module in modules {
        for (total, &dimension) in dimensions.iter_mut().zip(module.dimension_vector()) {
            *total = total
                .checked_add(dimension)
                .ok_or(PathAlgebraModuleError::DimensionOverflow)?;
        }
    }

    let field = first.right_acting_algebra().field().clone();
    let mut specification = Vec::new();
    for arrow in first.right_acting_algebra().quiver().arrows() {
        let source = arrow.source().index().zero_based();
        let target = arrow.target().index().zero_based();
        if dimensions[source] == 0 || dimensions[target] == 0 {
            continue;
        }
        let mut matrix = DenseMatrix::zero(dimensions[source], dimensions[target], field.clone());
        let mut row_offset = 0;
        let mut column_offset = 0;
        for module in modules {
            let source_dimension = module.dimension_vector()[source];
            let target_dimension = module.dimension_vector()[target];
            if source_dimension != 0 && target_dimension != 0 {
                let block = logical_submatrix(
                    &module.arrow_matrices()[arrow.index().zero_based()],
                    source_dimension,
                    target_dimension,
                )?;
                for row in 0..source_dimension {
                    for column in 0..target_dimension {
                        matrix[(row_offset + row, column_offset + column)] =
                            block[(row, column)].clone();
                    }
                }
            }
            row_offset += source_dimension;
            column_offset += target_dimension;
        }
        specification.push(NamedArrowMatrix::new(arrow.name(), matrix));
    }

    let direct_sum = right_module_over_path_algebra(
        first.right_acting_algebra().clone(),
        dimensions,
        &specification,
    )?
    .ok_or(PathAlgebraModuleError::GeneratedDirectSumViolatesRelations)?;
    let projections = construct_direct_sum_component_maps(&direct_sum, modules, true)?;
    let inclusions = construct_direct_sum_component_maps(&direct_sum, modules, false)?;
    direct_sum
        .0
        .direct_sum
        .set(DirectSumMetadata {
            summands: modules.to_vec(),
            projections: RwLock::new(
                projections
                    .iter()
                    .map(PathAlgebraModuleHomomorphism::downgrade)
                    .collect(),
            ),
            inclusions: RwLock::new(
                inclusions
                    .iter()
                    .map(PathAlgebraModuleHomomorphism::downgrade)
                    .collect(),
            ),
        })
        .map_err(|_| PathAlgebraModuleError::DirectSumMetadataAlreadySet)?;
    Ok(Some(direct_sum))
}
/// Return QPA's cached direct-sum projections in summand order.
///
/// # Errors
///
/// Rejects a module not built by [`direct_sum_of_qpa_modules`] and propagates
/// checked map reconstruction failures.
pub fn direct_sum_projections<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    direct_sum_component_maps(module, true)
}

/// Return QPA's cached direct-sum inclusions in summand order.
///
/// # Errors
///
/// Rejects a module not built by [`direct_sum_of_qpa_modules`] and propagates
/// checked map reconstruction failures.
pub fn direct_sum_inclusions<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    direct_sum_component_maps(module, false)
}

fn direct_sum_component_maps<F: Field>(
    module: &PathAlgebraModule<F>,
    projections: bool,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    let metadata = module
        .0
        .direct_sum
        .get()
        .ok_or(PathAlgebraModuleError::NotDirectSumModule)?;
    let cache = if projections {
        &metadata.projections
    } else {
        &metadata.inclusions
    };
    let mut slot = cache
        .write()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let mut result = Vec::with_capacity(metadata.summands.len());
    for summand_index in 0..metadata.summands.len() {
        if let Some(cached) = slot[summand_index].upgrade() {
            result.push(cached);
            continue;
        }
        let generated = construct_direct_sum_component_map(
            module,
            &metadata.summands,
            summand_index,
            projections,
        )?;
        slot[summand_index] = generated.downgrade();
        result.push(generated);
    }
    Ok(result)
}

fn construct_direct_sum_component_maps<F: Field>(
    direct_sum: &PathAlgebraModule<F>,
    summands: &[PathAlgebraModule<F>],
    projections: bool,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    let mut result = Vec::with_capacity(summands.len());
    for summand_index in 0..summands.len() {
        result.push(construct_direct_sum_component_map(
            direct_sum,
            summands,
            summand_index,
            projections,
        )?);
    }
    Ok(result)
}

fn construct_direct_sum_component_map<F: Field>(
    direct_sum: &PathAlgebraModule<F>,
    summands: &[PathAlgebraModule<F>],
    summand_index: usize,
    projections: bool,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let summand = &summands[summand_index];
    let field = direct_sum.right_acting_algebra().field().clone();
    let matrices = direct_sum
        .dimension_vector()
        .iter()
        .zip(summand.dimension_vector())
        .enumerate()
        .map(|(vertex, (&total, &component))| {
            let offset = summands[..summand_index]
                .iter()
                .map(|previous| previous.dimension_vector()[vertex])
                .sum::<usize>();
            let (rows, columns) = if projections {
                (physical_dimension(total), physical_dimension(component))
            } else {
                (physical_dimension(component), physical_dimension(total))
            };
            let mut matrix = DenseMatrix::zero(rows, columns, field.clone());
            if total != 0 && component != 0 {
                for coordinate in 0..component {
                    if projections {
                        matrix[(offset + coordinate, coordinate)] = field.one();
                    } else {
                        matrix[(coordinate, offset + coordinate)] = field.one();
                    }
                }
            }
            matrix
        })
        .collect::<Vec<_>>();
    if projections {
        let map = right_module_hom_over_algebra(direct_sum.clone(), summand.clone(), matrices)?;
        map.set_is_surjective(true);
        Ok(map)
    } else {
        let map = right_module_hom_over_algebra(summand.clone(), direct_sum.clone(), matrices)?;
        map.set_is_injective(true);
        Ok(map)
    }
}

/// Translate QPA's `PushOut(f, g)` construction.
///
/// The returned maps are ordered from `Range(g)` and `Range(f)` into the
/// pushout, matching QPA's `[f', g']`. `Ok(None)` is QPA's printed `fail`
/// result when the maps do not have structurally equal sources.
///
/// # Errors
///
/// Propagates direct-sum, map arithmetic, and cokernel construction failures.
pub fn push_out<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphismPair<F>>, PathAlgebraModuleError> {
    if first.source() != second.source() {
        return Ok(None);
    }
    let direct_sum = direct_sum_of_qpa_modules(&[second.range().clone(), first.range().clone()])?
        .ok_or(PathAlgebraModuleError::InternalNonemptyDirectSumFailed)?;
    let inclusions = direct_sum_inclusions(&direct_sum)?;
    let difference = first
        .compose(&inclusions[1])?
        .add(&second.compose(&inclusions[0])?.negate()?)?;
    let cokernel = cokernel_projection(&difference)?;
    Ok(Some((
        inclusions[0].compose(&cokernel)?,
        inclusions[1].compose(&cokernel)?,
    )))
}

/// Translate QPA's `PullBack(f, g)` construction.
///
/// The returned maps are ordered into `Source(g)` and `Source(f)`, matching
/// QPA's `[f', g']`. `Ok(None)` is QPA's printed `fail` result when the maps
/// do not have structurally equal ranges.
///
/// # Errors
///
/// Propagates direct-sum, map arithmetic, and kernel construction failures.
pub fn pull_back<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphismPair<F>>, PathAlgebraModuleError> {
    if first.range() != second.range() {
        return Ok(None);
    }
    let direct_sum = direct_sum_of_qpa_modules(&[first.source().clone(), second.source().clone()])?
        .ok_or(PathAlgebraModuleError::InternalNonemptyDirectSumFailed)?;
    let projections = direct_sum_projections(&direct_sum)?;
    let difference = projections[0]
        .compose(first)?
        .add(&projections[1].compose(second)?.negate()?)?;
    let kernel = kernel_inclusion(&difference)?;
    Ok(Some((
        kernel.compose(&projections[1])?,
        kernel.compose(&projections[0])?,
    )))
}

/// Translate QPA's `MorphismOnKernel(f, g, beta, alpha)`.
///
/// `Ok(None)` is QPA's `fail` for a noncommutative square.
///
/// # Errors
///
/// Rejects malformed diagram compositions and propagates exact kernel,
/// preimage, image, and checked-map construction failures.
pub fn morphism_on_kernel<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
    beta: &PathAlgebraModuleHomomorphism<F>,
    alpha: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if first.compose(alpha)? != beta.compose(second)? {
        return Ok(None);
    }
    let first_kernel = kernel_inclusion(first)?;
    let second_kernel = kernel_inclusion(second)?;
    let route = first_kernel.compose(beta)?;
    let mut images = Vec::with_capacity(first_kernel.source().dimension());
    for basis_element in first_kernel.source().basis_elements() {
        let image = route.image(&basis_element)?;
        images.push(
            second_kernel
                .preimage_representative(&image)?
                .ok_or(PathAlgebraModuleError::MissingInducedKernelPreimage)?,
        );
    }
    Ok(Some(homomorphism_from_images(
        first_kernel.source(),
        second_kernel.source(),
        &images,
    )?))
}

/// Translate QPA's `MorphismOnImage(f, g, beta, alpha)`.
///
/// `Ok(None)` is QPA's `fail` for a noncommutative square.
///
/// # Errors
///
/// Rejects malformed diagram compositions and propagates exact image,
/// preimage, and checked-map construction failures.
pub fn morphism_on_image<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
    beta: &PathAlgebraModuleHomomorphism<F>,
    alpha: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if first.compose(alpha)? != beta.compose(second)? {
        return Ok(None);
    }
    let first_image = image_projection(first)?;
    let second_image = image_projection(second)?;
    let route = beta.compose(&second_image)?;
    let mut images = Vec::with_capacity(first_image.range().dimension());
    for basis_element in first_image.range().basis_elements() {
        let lift = first_image
            .preimage_representative(&basis_element)?
            .ok_or(PathAlgebraModuleError::MissingInducedImagePreimage)?;
        images.push(route.image(&lift)?);
    }
    Ok(Some(homomorphism_from_images(
        first_image.range(),
        second_image.range(),
        &images,
    )?))
}

/// Translate QPA's `MorphismOnCoKernel(f, g, beta, alpha)`.
///
/// `Ok(None)` is QPA's `fail` for a noncommutative square.
///
/// # Errors
///
/// Rejects malformed diagram compositions and propagates exact cokernel,
/// preimage, image, and checked-map construction failures.
pub fn morphism_on_cokernel<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
    beta: &PathAlgebraModuleHomomorphism<F>,
    alpha: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if first.compose(alpha)? != beta.compose(second)? {
        return Ok(None);
    }
    let first_cokernel = cokernel_projection(first)?;
    let second_cokernel = cokernel_projection(second)?;
    let route = alpha.compose(&second_cokernel)?;
    let mut images = Vec::with_capacity(first_cokernel.range().dimension());
    for basis_element in first_cokernel.range().basis_elements() {
        let lift = first_cokernel
            .preimage_representative(&basis_element)?
            .ok_or(PathAlgebraModuleError::MissingInducedCokernelPreimage)?;
        images.push(route.image(&lift)?);
    }
    Ok(Some(homomorphism_from_images(
        first_cokernel.range(),
        second_cokernel.range(),
        &images,
    )?))
}

/// Translate QPA's `LiftingInclusionMorphisms(f, g)`.
///
/// QPA first constructs the cokernel of `f`, then returns `fail` unless both
/// maps are injective with a common structural range and the image of `g`
/// lies in the image of `f`. `Ok(None)` represents that `fail` result.
///
/// # Errors
///
/// Propagates cokernel, image, preimage, and checked-map construction
/// failures. A missing preimage after the containment test is reported as an
/// internal source-translation failure.
pub fn lifting_inclusion_morphisms<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    let cokernel = cokernel_projection(first)?;
    if !first.is_injective()
        || !second.is_injective()
        || first.range() != second.range()
        || !second.compose(&cokernel)?.is_zero()
    {
        return Ok(None);
    }
    let mut images = Vec::with_capacity(second.source().dimension());
    for basis_element in second.source().basis_elements() {
        let image = second.image(&basis_element)?;
        images.push(
            first
                .preimage_representative(&image)?
                .ok_or(PathAlgebraModuleError::MissingInclusionLiftPreimage)?,
        );
    }
    Ok(Some(homomorphism_from_images(
        second.source(),
        first.source(),
        &images,
    )?))
}

/// Translate QPA's `HomomorphismFromImages(M, N, genImages)` operation.
///
/// Images correspond to QPA's standard basis in vertex-major order. A
/// nonzero image must be supported at the same vertex as its basis vector;
/// additional images are ignored exactly as in the source loop.
///
/// # Errors
///
/// Reports a short image list or wrong support and propagates exact matrix and
/// checked module-map construction failures.
pub fn homomorphism_from_images<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
    images: &[PathAlgebraModuleElement<F>],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if images.len() < source.dimension() {
        return Err(PathAlgebraModuleError::BasisImageCount {
            actual: images.len(),
            expected: source.dimension(),
        });
    }
    let mut image_position = 0;
    for (vertex, &dimension) in source.dimension_vector().iter().enumerate() {
        for _ in 0..dimension {
            let image = &images[image_position];
            let support = image.support_vertices();
            if !support.is_empty() && support != [vertex + 1] {
                return Err(PathAlgebraModuleError::WrongBasisImageSupport {
                    position: image_position + 1,
                });
            }
            image_position += 1;
        }
    }
    let field = source.right_acting_algebra().field().clone();
    let mut offset = 0_usize;
    let matrices = source
        .dimension_vector()
        .iter()
        .zip(range.dimension_vector())
        .enumerate()
        .map(|(vertex, (&source_dimension, &range_dimension))| {
            if source_dimension == 0 {
                return Ok(DenseMatrix::zero(
                    1,
                    physical_dimension(range_dimension),
                    field.clone(),
                ));
            }
            let mut entries =
                Vec::with_capacity(source_dimension * physical_dimension(range_dimension));
            for image in &images[offset..offset + source_dimension] {
                entries.extend_from_slice(&image.external_representation()[vertex]);
            }
            offset += source_dimension;
            Ok(DenseMatrix::from_flat(
                source_dimension,
                physical_dimension(range_dimension),
                entries,
                field.clone(),
            )?)
        })
        .collect::<Result<Vec<_>, PathAlgebraModuleError>>()?;
    right_module_hom_over_algebra(source.clone(), range.clone(), matrices)
}

/// Translate QPA's `MultiplyListsOfMaps(projections, matrix, inclusions)`.
///
/// The computation is the literal row-by-matrix-by-column product in QPA's
/// left-to-right map-composition convention. Extra matrix rows are ignored by
/// the source loop.
///
/// # Errors
///
/// The intended QPA domain requires nonempty outer lists, at least one matrix
/// row per inclusion, and one entry per projection in each consumed row.
/// Composition, addition, and checked map failures propagate.
pub fn multiply_lists_of_maps<F: Field>(
    projections: &[PathAlgebraModuleHomomorphism<F>],
    matrix: &[Vec<PathAlgebraModuleHomomorphism<F>>],
    inclusions: &[PathAlgebraModuleHomomorphism<F>],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if projections.is_empty() {
        return Err(PathAlgebraModuleError::EmptyMapProjectionList);
    }
    if inclusions.is_empty() {
        return Err(PathAlgebraModuleError::EmptyMapInclusionList);
    }
    if matrix.len() < inclusions.len() {
        return Err(PathAlgebraModuleError::MapMatrixRowCount {
            actual: matrix.len(),
            expected: inclusions.len(),
        });
    }
    let mut outer_terms = Vec::with_capacity(inclusions.len());
    for (row_position, (row, inclusion)) in matrix.iter().zip(inclusions).enumerate() {
        if row.len() != projections.len() {
            return Err(PathAlgebraModuleError::MapMatrixColumnCount {
                row: row_position + 1,
                actual: row.len(),
                expected: projections.len(),
            });
        }
        let mut row_terms = projections
            .iter()
            .zip(row)
            .map(|(projection, entry)| projection.compose(entry));
        let Some(first_term) = row_terms.next() else {
            return Err(PathAlgebraModuleError::EmptyMapProjectionList);
        };
        let mut row_sum = first_term?;
        for term in row_terms {
            row_sum = row_sum.add(&term?)?;
        }
        outer_terms.push(row_sum.compose(inclusion)?);
    }
    let mut outer_terms = outer_terms.into_iter();
    let Some(mut sum) = outer_terms.next() else {
        return Err(PathAlgebraModuleError::EmptyMapInclusionList);
    };
    for term in outer_terms {
        sum = sum.add(&term)?;
    }
    Ok(sum)
}

/// Translate QPA's `UnderlyingLinearMap(f)` block-diagonal construction.
///
/// Vertex spaces occur in quiver order and the stored row-action matrix at
/// each vertex is copied into the corresponding diagonal block. QPA retains
/// one physical zero row or column for every zero-dimensional vertex space;
/// consequently the returned dimensions are the sums of the physical vertex
/// dimensions rather than the logical module dimensions.
///
/// # Errors
///
/// Reports address-size overflow while summing physical dimensions.
pub fn underlying_linear_map<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<DenseMatrix<F>, PathAlgebraModuleError> {
    let source_dimension =
        homomorphism
            .vertex_matrices()
            .iter()
            .try_fold(0_usize, |total, matrix| {
                total
                    .checked_add(matrix.nrows())
                    .ok_or(PathAlgebraModuleError::DimensionOverflow)
            })?;
    let range_dimension =
        homomorphism
            .vertex_matrices()
            .iter()
            .try_fold(0_usize, |total, matrix| {
                total
                    .checked_add(matrix.ncols())
                    .ok_or(PathAlgebraModuleError::DimensionOverflow)
            })?;
    let field = homomorphism.source().right_acting_algebra().field().clone();
    let mut result = DenseMatrix::zero(source_dimension, range_dimension, field);
    let mut source_offset = 0;
    let mut range_offset = 0;
    for matrix in homomorphism.vertex_matrices() {
        for row in 0..matrix.nrows() {
            for column in 0..matrix.ncols() {
                result[(source_offset + row, range_offset + column)] =
                    matrix[(row, column)].clone();
            }
        }
        source_offset += matrix.nrows();
        range_offset += matrix.ncols();
    }
    Ok(result)
}

/// Translate QPA's Bongartz--Mroz `CommonDirectSummand(M, N)` search.
///
/// Hom bases, loop order, the source's floating-point `Ceil(Log2(l))`
/// exponent, repeated squaring, and the returned image/kernel order are all
/// retained. The result enum distinguishes QPA's printed-then-`fail` parent
/// mismatch from its `false` exhausted-search result.
///
/// # Errors
///
/// Propagates exact Hom, composition, image, and kernel construction.
pub fn common_direct_summand<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<CommonDirectSummandResult<F>, PathAlgebraModuleError> {
    if source.right_acting_algebra() != range.right_acting_algebra() {
        return Ok(CommonDirectSummandResult::DifferentAlgebras);
    }
    let source_to_range = hom_over_algebra(source, range)?
        .ok_or(PathAlgebraModuleError::InternalCommonSummandHomRejected)?;
    let range_to_source = hom_over_algebra(range, source)?
        .ok_or(PathAlgebraModuleError::InternalCommonSummandHomRejected)?;
    if source_to_range.is_empty() || range_to_source.is_empty() {
        return Ok(CommonDirectSummandResult::NoCommonSummand);
    }

    let source_maximum = source.dimension_vector().iter().copied().max().unwrap_or(0);
    let range_maximum = range.dimension_vector().iter().copied().max().unwrap_or(0);
    let fitting_length = if source_maximum == range_maximum {
        range_maximum
    } else {
        source_maximum.min(range_maximum) + 1
    };
    if fitting_length == 0 {
        return Ok(CommonDirectSummandResult::NoCommonSummand);
    }
    let squarings = usize::BITS as usize - (fitting_length - 1).leading_zeros() as usize;

    for range_map in &range_to_source {
        for source_map in &source_to_range {
            let fitting_map = if fitting_length > 1 {
                let mut power = source_map.compose(range_map)?;
                for _ in 0..squarings {
                    power = power.compose(&power)?;
                }
                power.compose(source_map)?
            } else {
                source_map.clone()
            };
            let source_endomorphism = fitting_map.compose(range_map)?;
            if source_endomorphism.is_zero() {
                continue;
            }
            let range_endomorphism = range_map.compose(&fitting_map)?;
            return Ok(CommonDirectSummandResult::Found(CommonDirectSummand {
                source_common: image_of_module_homomorphism(&source_endomorphism)?,
                source_complement: kernel_of_module_homomorphism(&source_endomorphism)?,
                range_common: image_of_module_homomorphism(&range_endomorphism)?,
                range_complement: kernel_of_module_homomorphism(&range_endomorphism)?,
            }));
        }
    }
    Ok(CommonDirectSummandResult::NoCommonSummand)
}

/// Translate QPA's repeated `MaximalCommonDirectSummand(M, N)` loop.
///
/// `None` is QPA's `false`; as in the source this also collapses the first
/// operation's different-algebra `fail` after its diagnostic.
///
/// # Errors
///
/// Propagates [`common_direct_summand`] failures.
pub fn maximal_common_direct_summand<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Option<MaximalCommonDirectSummand<F>>, PathAlgebraModuleError> {
    let mut source_remainder = source.clone();
    let mut range_remainder = range.clone();
    let mut common_summands = Vec::new();

    loop {
        let CommonDirectSummandResult::Found(common) =
            common_direct_summand(&source_remainder, &range_remainder)?
        else {
            break;
        };
        common_summands.push(common.source_common);
        source_remainder = common.source_complement;
        range_remainder = common.range_complement;
        if source_remainder.dimension() == 0 || range_remainder.dimension() == 0 {
            break;
        }
    }

    Ok(
        (!common_summands.is_empty()).then_some(MaximalCommonDirectSummand {
            common_summands,
            source_remainder,
            range_remainder,
        }),
    )
}

/// Translate QPA's Boolean `IsomorphicModules(M, N)` operation.
///
/// This is intentionally distinct from [`isomorphism_of_modules`]: it follows
/// QPA's older maximal-common-summand decision algorithm and returns no map.
///
/// # Errors
///
/// Propagates [`maximal_common_direct_summand`] failures.
pub fn isomorphic_modules<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<bool, PathAlgebraModuleError> {
    if source.dimension_vector() != range.dimension_vector() {
        return Ok(false);
    }
    if source.dimension() == 0 && range.dimension() == 0 {
        return Ok(true);
    }
    let Some(common) = maximal_common_direct_summand(source, range)? else {
        return Ok(false);
    };
    Ok(common.source_remainder.dimension() == 0 && common.range_remainder.dimension() == 0)
}

/// Translate QPA's `IsDirectSummand(M, N)` decision operation.
///
/// # Errors
///
/// Propagates [`maximal_common_direct_summand`] failures.
pub fn is_direct_summand<F: Field>(
    candidate: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
) -> Result<bool, PathAlgebraModuleError> {
    if !dimension_vector_partial_order(candidate, module) {
        return Ok(false);
    }
    let Some(common) = maximal_common_direct_summand(candidate, module)? else {
        return Ok(false);
    };
    Ok(common.source_remainder.dimension() == 0)
}

/// Translate QPA's product-span `IsInAdditiveClosure(M, N)` test.
///
/// `None` is QPA's `fail` for modules over different algebras. The source
/// compares the span dimension of all `Hom(M,N) Hom(N,M)` products with the
/// dimension of `End(M)`; Rust retains that test rather than substituting an
/// isomorphism/decomposition criterion.
///
/// # Errors
///
/// Propagates Hom, composition, and exact row-space construction.
pub fn is_in_additive_closure<F: Field>(
    candidate: &PathAlgebraModule<F>,
    generator: &PathAlgebraModule<F>,
) -> Result<Option<bool>, PathAlgebraModuleError> {
    if candidate.right_acting_algebra() != generator.right_acting_algebra() {
        return Ok(None);
    }
    let candidate_to_generator = hom_over_algebra(candidate, generator)?
        .ok_or(PathAlgebraModuleError::InternalAdditiveClosureHomRejected)?;
    let generator_to_candidate = hom_over_algebra(generator, candidate)?
        .ok_or(PathAlgebraModuleError::InternalAdditiveClosureHomRejected)?;
    let endomorphism_dimension = hom_dimension_over_algebra(candidate, candidate)?
        .ok_or(PathAlgebraModuleError::InternalAdditiveClosureHomRejected)?;
    let mut products = Vec::with_capacity(
        candidate_to_generator
            .len()
            .saturating_mul(generator_to_candidate.len()),
    );
    for forward in &candidate_to_generator {
        for backward in &generator_to_candidate {
            products.push(forward.compose(backward)?);
        }
    }
    let span_dimension = if products.is_empty() {
        0
    } else {
        let width = products[0].flat_vertex_matrix_entries().len();
        DenseMatrix::from_flat(
            products.len(),
            width,
            products
                .iter()
                .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
                .collect(),
            candidate.right_acting_algebra().field().clone(),
        )?
        .semi_echelon_basis()
        .nrows()
    };
    Ok(Some(span_dimension == endomorphism_dimension))
}

/// Translate QPA's recursive `IsomorphismOfModules(M, N)` operation.
///
/// `Ok(None)` is QPA's `false`. After the same algebra, dimension-vector,
/// zero, and structural-equality branches, this follows QPA's ordered Hom-basis
/// search for a common Fitting summand and recursively identifies the two
/// complements.
///
/// # Errors
///
/// Propagates exact Hom, composition, kernel/image, split-inverse, preimage,
/// basis-image, and checked map construction failures. Failures of the
/// theoretically split internal inclusions are reported explicitly.
#[allow(clippy::too_many_lines)]
pub fn isomorphism_of_modules<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if source.right_acting_algebra() != range.right_acting_algebra()
        || source.dimension_vector() != range.dimension_vector()
    {
        return Ok(None);
    }
    if source.dimension() == 0 && range.dimension() == 0 {
        return Ok(Some(zero_module_mapping(source.clone(), range.clone())?));
    }
    if source == range {
        let field = source.right_acting_algebra().field().clone();
        let matrices = source
            .dimension_vector()
            .iter()
            .map(|&dimension| {
                if dimension == 0 {
                    DenseMatrix::zero(1, 1, field.clone())
                } else {
                    DenseMatrix::identity(dimension, field.clone())
                }
            })
            .collect();
        return Ok(Some(right_module_hom_over_algebra(
            source.clone(),
            range.clone(),
            matrices,
        )?));
    }

    let Some((source_endomorphism, range_endomorphism, source_to_range)) =
        common_splitting_maps(source, range)?
    else {
        return Ok(None);
    };
    let source_summand_inclusion = image_inclusion(&source_endomorphism)?;
    let source_summand_retraction = right_inverse_of_homomorphism(&source_summand_inclusion)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismSummandNotSplit)?;
    let range_summand_inclusion = image_inclusion(&range_endomorphism)?;
    let range_summand_retraction = right_inverse_of_homomorphism(&range_summand_inclusion)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismSummandNotSplit)?;

    let source_idempotent = source_summand_retraction.compose(&source_summand_inclusion)?;
    let source_common_inclusion = image_inclusion(&source_idempotent)?;
    let source_complement_inclusion = kernel_inclusion(&source_idempotent)?;
    let range_idempotent = range_summand_retraction.compose(&range_summand_inclusion)?;
    let range_common_inclusion = image_inclusion(&range_idempotent)?;
    let range_complement_inclusion = kernel_inclusion(&range_idempotent)?;

    if source_complement_inclusion.source().dimension() == 0
        && range_complement_inclusion.source().dimension() == 0
    {
        return Ok(Some(source_endomorphism.compose(&source_to_range)?));
    }

    let common_route = source_common_inclusion
        .compose(&source_endomorphism)?
        .compose(&source_to_range)?;
    let mut common_images = Vec::with_capacity(common_route.source().dimension());
    for basis_element in common_route.source().basis_elements() {
        let image = common_route.image(&basis_element)?;
        common_images.push(
            range_common_inclusion
                .preimage_representative(&image)?
                .ok_or(PathAlgebraModuleError::MissingIsomorphismSummandPreimage)?,
        );
    }
    let common_isomorphism = homomorphism_from_images(
        source_common_inclusion.source(),
        range_common_inclusion.source(),
        &common_images,
    )?;
    let Some(complement_isomorphism) = isomorphism_of_modules(
        source_complement_inclusion.source(),
        range_complement_inclusion.source(),
    )?
    else {
        return Ok(None);
    };

    let source_common_retraction = right_inverse_of_homomorphism(&source_common_inclusion)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismSummandNotSplit)?;
    let source_complement_retraction = right_inverse_of_homomorphism(&source_complement_inclusion)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismSummandNotSplit)?;
    let common_map = source_common_retraction
        .compose(&common_isomorphism)?
        .compose(&range_common_inclusion)?;
    let complement_map = source_complement_retraction
        .compose(&complement_isomorphism)?
        .compose(&range_complement_inclusion)?;
    Ok(Some(common_map.add(&complement_map)?))
}

fn common_splitting_maps<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModuleHomomorphismTriple<F>>, PathAlgebraModuleError> {
    let source_to_range = hom_over_algebra(source, range)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismHomRejected)?;
    let range_to_source = hom_over_algebra(range, source)?
        .ok_or(PathAlgebraModuleError::InternalIsomorphismHomRejected)?;
    if source_to_range.is_empty() || range_to_source.is_empty() {
        return Ok(None);
    }
    let source_maximum = source.dimension_vector().iter().copied().max().unwrap_or(0);
    let range_maximum = range.dimension_vector().iter().copied().max().unwrap_or(0);
    let fitting_length = if source_maximum == range_maximum {
        range_maximum
    } else {
        source_maximum.min(range_maximum) + 1
    };

    for range_map in &range_to_source {
        for source_map in &source_to_range {
            let fitting_map = if fitting_length > 1 {
                homomorphism_power(&source_map.compose(range_map)?, fitting_length - 1)?
                    .compose(source_map)?
            } else {
                source_map.clone()
            };
            let source_endomorphism = fitting_map.compose(range_map)?;
            if !source_endomorphism.is_zero() {
                return Ok(Some((
                    source_endomorphism,
                    range_map.compose(&fitting_map)?,
                    source_map.clone(),
                )));
            }
        }
    }
    Ok(None)
}

/// Translate QPA's binary `IntersectionOfSubmodules(f, g)`.
///
/// The tuple contains the maps into the second and first source followed by
/// the canonical inclusion of their intersection into the common range.
///
/// # Errors
///
/// Requires two monomorphisms with structurally equal ranges and propagates
/// QPA's pullback construction failures.
pub fn intersection_of_submodules<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphismTriple<F>, PathAlgebraModuleError> {
    validate_submodule_pair(first, second)?;
    let (to_second, to_first) = pull_back(first, second)?
        .ok_or(PathAlgebraModuleError::InternalCommonRangePullbackFailed)?;
    let inclusion = to_first.compose(first)?;
    Ok((to_second, to_first, inclusion))
}

/// Translate QPA's list `IntersectionOfSubmodules` method.
///
/// A singleton list returns the identical entered map, while longer lists
/// fold binary pullbacks in entered order.
///
/// # Errors
///
/// Preserves QPA's validation order for an empty list, source algebras,
/// injectivity, and common ranges, then propagates binary construction errors.
pub fn intersection_of_submodules_list<F: Field>(
    homomorphisms: &[PathAlgebraModuleHomomorphism<F>],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    validate_submodule_list(homomorphisms, true)?;
    let mut intersection = homomorphisms[0].clone();
    for homomorphism in &homomorphisms[1..] {
        intersection = intersection_of_submodules(&intersection, homomorphism)?.2;
    }
    Ok(intersection)
}

/// Translate QPA's binary `SumOfSubmodules(f, g)`.
///
/// The tuple is the inclusion of the sum into the common range followed by
/// the maps from the first and second source into that sum.
///
/// # Errors
///
/// Requires two monomorphisms with structurally equal ranges and propagates
/// QPA's direct-sum and image-factorization construction failures.
pub fn sum_of_submodules<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphismTriple<F>, PathAlgebraModuleError> {
    validate_submodule_pair(first, second)?;
    let direct_sum = direct_sum_of_qpa_modules(&[first.source().clone(), second.source().clone()])?
        .ok_or(PathAlgebraModuleError::InternalNonemptyDirectSumFailed)?;
    let inclusions = direct_sum_inclusions(&direct_sum)?;
    let projections = direct_sum_projections(&direct_sum)?;
    let combined = projections[0]
        .compose(first)?
        .add(&projections[1].compose(second)?)?;
    let (projection, inclusion) = image_projection_inclusion(&combined)?;
    Ok((
        inclusion,
        inclusions[0].compose(&projection)?,
        inclusions[1].compose(&projection)?,
    ))
}

/// Translate QPA's list `SumOfSubmodules` method.
///
/// A singleton list returns the identical entered map, while longer lists
/// fold binary image sums in entered order.
///
/// # Errors
///
/// Preserves QPA's validation order for an empty list, source algebras,
/// common ranges, and injectivity, then propagates binary construction errors.
pub fn sum_of_submodules_list<F: Field>(
    homomorphisms: &[PathAlgebraModuleHomomorphism<F>],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    validate_submodule_list(homomorphisms, false)?;
    let mut sum = homomorphisms[0].clone();
    for homomorphism in &homomorphisms[1..] {
        sum = sum_of_submodules(&sum, homomorphism)?.0;
    }
    Ok(sum)
}

fn validate_submodule_pair<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<(), PathAlgebraModuleError> {
    if !first.is_injective() {
        return Err(PathAlgebraModuleError::FirstSubmoduleMapNotInjective);
    }
    if !second.is_injective() {
        return Err(PathAlgebraModuleError::SecondSubmoduleMapNotInjective);
    }
    if first.range() != second.range() {
        return Err(PathAlgebraModuleError::SubmoduleRangesDiffer);
    }
    Ok(())
}

fn validate_submodule_list<F: Field>(
    homomorphisms: &[PathAlgebraModuleHomomorphism<F>],
    injectivity_before_range: bool,
) -> Result<(), PathAlgebraModuleError> {
    let Some(first) = homomorphisms.first() else {
        return Err(PathAlgebraModuleError::EmptySubmoduleList);
    };
    if homomorphisms.iter().any(|homomorphism| {
        homomorphism.source().right_acting_algebra() != first.source().right_acting_algebra()
    }) {
        return Err(PathAlgebraModuleError::SubmoduleListDifferentAlgebras);
    }
    let injective = || {
        homomorphisms
            .iter()
            .all(PathAlgebraModuleHomomorphism::is_injective)
    };
    let common_range = || {
        homomorphisms
            .iter()
            .all(|homomorphism| homomorphism.range() == first.range())
    };
    if injectivity_before_range {
        if !injective() {
            return Err(PathAlgebraModuleError::SomeSubmoduleMapNotInjective);
        }
        if !common_range() {
            return Err(PathAlgebraModuleError::SubmoduleRangesDiffer);
        }
    } else {
        if !common_range() {
            return Err(PathAlgebraModuleError::SubmoduleRangesDiffer);
        }
        if !injective() {
            return Err(PathAlgebraModuleError::SomeSubmoduleMapNotInjective);
        }
    }
    Ok(())
}

/// Translate QPA's `RejectOfModule(N, M)` intersection-of-kernels method.
///
/// # Errors
///
/// Requires both modules over the same algebra and propagates exact Hom,
/// kernel, and submodule-intersection failures.
pub fn reject_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
    rejected_by: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if module.right_acting_algebra() != rejected_by.right_acting_algebra() {
        return Err(PathAlgebraModuleError::RejectModulesOverDifferentAlgebras);
    }
    let homomorphisms = hom_over_algebra(module, rejected_by)?
        .ok_or(PathAlgebraModuleError::InternalCommonAlgebraHomFailed)?;
    match homomorphisms.as_slice() {
        [] => identity_module_mapping(module.clone()),
        [homomorphism] => kernel_inclusion(homomorphism),
        _ => {
            let kernels = homomorphisms
                .iter()
                .map(kernel_inclusion)
                .collect::<Result<Vec<_>, _>>()?;
            intersection_of_submodules_list(&kernels)
        }
    }
}

#[derive(Clone, Debug)]
struct NaturalRowQuotient<F: Field> {
    projection: DenseMatrix<F>,
    lifts: DenseMatrix<F>,
}

fn natural_quotient_by_row_space<F: Field>(
    image_rows: &DenseMatrix<F>,
) -> Result<NaturalRowQuotient<F>, MatrixError> {
    let field = image_rows.field().clone();
    let dimension = image_rows.ncols();
    let kernel_basis = image_rows.semi_echelon_basis();
    let mut mutable_basis = kernel_basis.clone();
    let mut complement = Vec::new();
    for coordinate in 0..dimension {
        let mut standard = vec![field.zero(); dimension];
        standard[coordinate] = field.one();
        let mut candidate = mutable_basis.as_flat().to_vec();
        candidate.extend_from_slice(&standard);
        let candidate = DenseMatrix::from_flat(
            mutable_basis.nrows() + 1,
            dimension,
            candidate,
            field.clone(),
        )?;
        let semi = candidate.semi_echelon_basis();
        if semi.nrows() > mutable_basis.nrows() {
            complement.push(standard);
            mutable_basis = semi;
        }
    }

    let quotient_dimension = complement.len();
    let lifts = DenseMatrix::from_flat(
        quotient_dimension,
        dimension,
        complement.iter().flatten().cloned().collect(),
        field.clone(),
    )?;
    let mut basis_entries = kernel_basis.as_flat().to_vec();
    basis_entries.extend(complement.iter().flatten().cloned());
    let full_basis = DenseMatrix::from_flat(dimension, dimension, basis_entries, field.clone())?;
    let mut projection_entries = Vec::with_capacity(dimension * quotient_dimension);
    for coordinate in 0..dimension {
        let mut standard = vec![field.zero(); dimension];
        standard[coordinate] = field.one();
        let coefficients = full_basis
            .solution_row(&standard)?
            .expect("a kernel basis extended by independent vectors spans the ambient space");
        projection_entries
            .extend_from_slice(&coefficients[dimension - quotient_dimension..dimension]);
    }
    Ok(NaturalRowQuotient {
        projection: DenseMatrix::from_flat(
            dimension,
            quotient_dimension,
            projection_entries,
            field,
        )?,
        lifts,
    })
}

fn logical_submatrix<F: Field>(
    matrix: &DenseMatrix<F>,
    rows: usize,
    columns: usize,
) -> Result<DenseMatrix<F>, MatrixError> {
    let mut entries = Vec::with_capacity(rows * columns);
    for row in 0..rows {
        entries.extend_from_slice(&matrix.row(row)[..columns]);
    }
    DenseMatrix::from_flat(rows, columns, entries, matrix.field().clone())
}

fn logical_vertex_map<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    vertex: usize,
) -> DenseMatrix<F> {
    logical_submatrix(
        &homomorphism.vertex_matrices()[vertex],
        homomorphism.source().dimension_vector()[vertex],
        homomorphism.range().dimension_vector()[vertex],
    )
    .expect("a checked module homomorphism retains every logical vertex submatrix")
}

fn dimension_offsets(dimensions: &[usize]) -> Vec<usize> {
    let mut offset = 0;
    dimensions
        .iter()
        .map(|&dimension| {
            let result = offset;
            offset += dimension;
            result
        })
        .collect()
}

fn quiver_algebra_path_element<F: Field>(
    algebra: &QuiverAlgebra<F>,
    path: Path,
) -> Result<QuiverAlgebraElement<F>, PathAlgebraModuleError> {
    match algebra {
        QuiverAlgebra::Path(algebra) => Ok(QuiverAlgebraElement::Path(algebra.element(path)?)),
        QuiverAlgebra::Quotient(algebra) => Ok(QuiverAlgebraElement::Quotient(algebra.path(path)?)),
    }
}
