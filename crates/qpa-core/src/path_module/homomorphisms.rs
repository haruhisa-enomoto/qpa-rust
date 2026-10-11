/// Translate QPA's checked `RightModuleHomOverAlgebra` constructor.
///
/// # Errors
///
/// Requires structurally equal acting algebras, one correctly shaped matrix
/// over the common field per vertex, and every arrow square to commute.
pub fn right_module_hom_over_algebra<F: Field>(
    source: PathAlgebraModule<F>,
    range: PathAlgebraModule<F>,
    vertex_matrices: Vec<DenseMatrix<F>>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if source.right_acting_algebra() != range.right_acting_algebra() {
        return Err(PathAlgebraModuleError::DifferentHomomorphismAlgebras);
    }
    if vertex_matrices.len() != source.dimension_vector().len() {
        return Err(PathAlgebraModuleError::HomomorphismMatrixCount {
            actual: vertex_matrices.len(),
            expected: source.dimension_vector().len(),
        });
    }
    for (vertex, matrix) in vertex_matrices.iter().enumerate() {
        if matrix.field() != source.right_acting_algebra().field() {
            return Err(PathAlgebraModuleError::HomomorphismMatrixField { vertex: vertex + 1 });
        }
        let expected = (
            physical_dimension(source.dimension_vector()[vertex]),
            physical_dimension(range.dimension_vector()[vertex]),
        );
        let actual = (matrix.nrows(), matrix.ncols());
        if actual != expected
            || (source.dimension_vector()[vertex] == 0 || range.dimension_vector()[vertex] == 0)
                && !matrix.is_zero()
        {
            return Err(PathAlgebraModuleError::HomomorphismMatrixShape {
                vertex: vertex + 1,
                actual,
                expected,
            });
        }
    }

    for arrow in source.right_acting_algebra().quiver().arrows() {
        let source_vertex = arrow.source().index().zero_based();
        let target_vertex = arrow.target().index().zero_based();
        let arrow_index = arrow.index().zero_based();
        let left =
            source.arrow_matrices()[arrow_index].multiply(&vertex_matrices[target_vertex])?;
        let right =
            vertex_matrices[source_vertex].multiply(&range.arrow_matrices()[arrow_index])?;
        if left != right {
            return Err(PathAlgebraModuleError::NoncommutingHomomorphismSquare {
                arrow: arrow.name().into(),
            });
        }
    }

    Ok(module_homomorphism_from_verified_matrices(
        source,
        range,
        vertex_matrices,
    ))
}

// Only the checked constructor and operations that preserve its invariants
// may call this. Module representations and map matrices are immutable.
fn module_homomorphism_from_verified_matrices<F: Field>(
    source: PathAlgebraModule<F>,
    range: PathAlgebraModule<F>,
    vertex_matrices: Vec<DenseMatrix<F>>,
) -> PathAlgebraModuleHomomorphism<F> {
    PathAlgebraModuleHomomorphism(Arc::new(PathAlgebraModuleHomomorphismData {
        source,
        range,
        vertex_matrices,
        kernel_inclusion: OnceLock::new(),
        image_projection_inclusion: OnceLock::new(),
        cokernel_projection: OnceLock::new(),
        dual_homomorphism: OnceLock::new(),
        star_homomorphism: OnceLock::new(),
        transpose_homomorphism: OnceLock::new(),
        d_tr_homomorphism: OnceLock::new(),
        tr_d_homomorphism: OnceLock::new(),
        nakayama_homomorphism: OnceLock::new(),
        opposite_nakayama_homomorphism: OnceLock::new(),
        is_injective: OnceLock::new(),
        is_surjective: OnceLock::new(),
        is_split_monomorphism: OnceLock::new(),
        is_split_epimorphism: OnceLock::new(),
        is_right_minimal: OnceLock::new(),
        is_left_minimal: OnceLock::new(),
        left_inverse_weak: RwLock::new(None),
        left_inverse_fail: OnceLock::new(),
        right_inverse_weak: RwLock::new(None),
        right_inverse_fail: OnceLock::new(),
        right_minimal_version: OnceLock::new(),
        left_minimal_version: OnceLock::new(),
        kernel_of_what: RwLock::new(None),
        image_of_what: RwLock::new(None),
        cokernel_of_what: RwLock::new(None),
    }))
}
/// Translate QPA's cached `ZeroModule(A)` attribute.
///
/// The module has one zero dimension for every vertex and no entered arrow
/// matrices. A weak algebra-side cache preserves GAP object identity while a
/// returned zero module remains live without leaking an algebra/module cycle.
///
/// # Errors
///
/// Reports the theoretically impossible rejection of the all-zero
/// representation by the checked module constructor.
pub fn zero_module<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    if let Some(cached) = algebra.cached_zero_module() {
        return Ok(cached);
    }
    let generated = right_module_over_path_algebra(
        algebra.clone(),
        vec![0; algebra.quiver().number_of_vertices()],
        &[],
    )?
    .ok_or(PathAlgebraModuleError::GeneratedZeroModuleViolatesRelations)?;
    Ok(algebra.cache_zero_module(&generated))
}

/// Translate QPA's zero mapping between two modules.
///
/// # Errors
///
/// Rejects modules over different algebras.
pub fn zero_module_mapping<F: Field>(
    source: PathAlgebraModule<F>,
    range: PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if source.right_acting_algebra() != range.right_acting_algebra() {
        return Err(PathAlgebraModuleError::DifferentHomomorphismAlgebras);
    }
    let field = source.right_acting_algebra().field().clone();
    let matrices = source
        .dimension_vector()
        .iter()
        .zip(range.dimension_vector())
        .map(|(&source_dimension, &range_dimension)| {
            DenseMatrix::zero(
                physical_dimension(source_dimension),
                physical_dimension(range_dimension),
                field.clone(),
            )
        })
        .collect();
    right_module_hom_over_algebra(source, range, matrices)
}

/// Translate QPA's identity mapping on a module.
///
/// At a zero-dimensional vertex QPA stores a `1 x 1` zero matrix rather than
/// an identity matrix.
///
/// # Errors
///
/// Propagates the checked module-map constructor.
pub fn identity_module_mapping<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let field = module.right_acting_algebra().field().clone();
    let matrices = module
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
    right_module_hom_over_algebra(module.clone(), module, matrices)
}

/// Translate QPA's `InverseOfIsomorphism(f)` operation.
///
/// QPA retains physical zero matrices at zero-dimensional vertices and takes
/// the ordinary exact inverse of every other vertex matrix.
///
/// # Errors
///
/// Rejects a map that is not an isomorphism and propagates exact matrix or
/// checked module-map construction failures.
pub fn inverse_of_isomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if !homomorphism.is_isomorphism() {
        return Err(PathAlgebraModuleError::EnteredMapNotIsomorphism);
    }
    let inverse_matrices = homomorphism
        .vertex_matrices()
        .iter()
        .map(|matrix| {
            if matrix.is_zero() {
                Ok(matrix.clone())
            } else {
                matrix
                    .inverse()?
                    .ok_or(PathAlgebraModuleError::InternalIsomorphismMatrixSingular)
            }
        })
        .collect::<Result<Vec<_>, PathAlgebraModuleError>>()?;
    right_module_hom_over_algebra(
        homomorphism.range().clone(),
        homomorphism.source().clone(),
        inverse_matrices,
    )
}

/// Translate QPA's cached `RightInverseOfHomomorphism(f)` attribute.
///
/// `Ok(None)` is QPA's returned `false`: either `f` is not injective, the
/// reverse Hom basis is empty, or GAP's row solver finds no retraction.
/// Successful calls install the same reciprocal inverse attributes and split
/// properties as the source method.
///
/// # Errors
///
/// Propagates Hom-space construction, composition, exact row solving, and
/// checked linear-combination failures.
pub fn right_inverse_of_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_right_inverse() {
        return Ok(Some(cached));
    }
    if homomorphism.right_inverse_attribute_is_fail() {
        return Ok(None);
    }
    if !homomorphism.is_injective() {
        homomorphism.mark_right_inverse_attribute_fail();
        return Ok(None);
    }
    let candidates = hom_over_algebra(homomorphism.range(), homomorphism.source())?
        .ok_or(PathAlgebraModuleError::InternalSplittingHomRejected)?;
    if candidates.is_empty() {
        homomorphism.mark_right_inverse_attribute_fail();
        return Ok(None);
    }
    let composites = candidates
        .iter()
        .map(|candidate| homomorphism.compose(candidate))
        .collect::<Result<Vec<_>, _>>()?;
    let identity = identity_module_mapping(homomorphism.source().clone())?;
    let Some(coefficients) = solve_map_linear_combination(&composites, &identity)? else {
        let _ = homomorphism.cache_is_split_monomorphism(false);
        homomorphism.mark_right_inverse_attribute_fail();
        return Ok(None);
    };
    let inverse = linear_combination_of_homomorphisms(&candidates, &coefficients)?;
    let _ = inverse.cache_is_split_epimorphism(true);
    let _ = homomorphism.cache_is_split_monomorphism(true);
    let _ = inverse.cache_left_inverse(homomorphism);
    Ok(Some(homomorphism.cache_right_inverse(&inverse)))
}

/// Translate QPA's cached `LeftInverseOfHomomorphism(f)` attribute.
///
/// `Ok(None)` is QPA's returned `false`: either `f` is not surjective, the
/// reverse Hom basis is empty, or GAP's row solver finds no section.
/// Successful calls install the source method's reciprocal attributes and
/// split properties.
///
/// # Errors
///
/// Propagates Hom-space construction, composition, exact row solving, and
/// checked linear-combination failures.
pub fn left_inverse_of_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_left_inverse() {
        return Ok(Some(cached));
    }
    if homomorphism.left_inverse_attribute_is_fail() {
        return Ok(None);
    }
    if !homomorphism.is_surjective() {
        homomorphism.mark_left_inverse_attribute_fail();
        return Ok(None);
    }
    let candidates = hom_over_algebra(homomorphism.range(), homomorphism.source())?
        .ok_or(PathAlgebraModuleError::InternalSplittingHomRejected)?;
    if candidates.is_empty() {
        homomorphism.mark_left_inverse_attribute_fail();
        return Ok(None);
    }
    let composites = candidates
        .iter()
        .map(|candidate| candidate.compose(homomorphism))
        .collect::<Result<Vec<_>, _>>()?;
    let identity = identity_module_mapping(homomorphism.range().clone())?;
    let Some(coefficients) = solve_map_linear_combination(&composites, &identity)? else {
        let _ = homomorphism.cache_is_split_epimorphism(false);
        homomorphism.mark_left_inverse_attribute_fail();
        return Ok(None);
    };
    let inverse = linear_combination_of_homomorphisms(&candidates, &coefficients)?;
    let _ = inverse.cache_is_split_monomorphism(true);
    let _ = homomorphism.cache_is_split_epimorphism(true);
    let _ = inverse.cache_right_inverse(homomorphism);
    Ok(Some(homomorphism.cache_left_inverse(&inverse)))
}

/// Translate QPA's `IsSplitMonomorphism(f)` property method.
///
/// # Errors
///
/// Propagates [`right_inverse_of_homomorphism`].
pub fn is_split_monomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<bool, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_is_split_monomorphism() {
        return Ok(cached);
    }
    let value = right_inverse_of_homomorphism(homomorphism)?.is_some();
    Ok(homomorphism.cache_is_split_monomorphism(value))
}

/// Translate QPA's `IsSplitEpimorphism(f)` property method.
///
/// # Errors
///
/// Propagates [`left_inverse_of_homomorphism`].
pub fn is_split_epimorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<bool, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_is_split_epimorphism() {
        return Ok(cached);
    }
    let value = left_inverse_of_homomorphism(homomorphism)?.is_some();
    Ok(homomorphism.cache_is_split_epimorphism(value))
}

/// Translate QPA's annihilator-radical `IsRightMinimal(f)` property method.
///
/// QPA takes the left nullspace of the ordered maps `x * f`, rewrites those
/// relations in its ordered `End(Source(f))` basis, and checks that every
/// resulting annihilator element lies in the Jacobson radical. The computed
/// property is cached on the identical map object.
///
/// # Errors
///
/// Propagates Hom, composition, exact nullspace, endomorphism conversion, and
/// radical computation failures.
pub fn is_right_minimal<F: QpaCoefficientField>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<bool, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_is_right_minimal() {
        return Ok(cached);
    }
    let source = homomorphism.source();
    let endomorphisms = hom_over_algebra(source, source)?
        .ok_or(PathAlgebraModuleError::InternalMinimalityHomRejected)?;
    let composites = endomorphisms
        .iter()
        .map(|endomorphism| endomorphism.compose(homomorphism))
        .collect::<Result<Vec<_>, _>>()?;
    let result = annihilator_is_radical(source, &endomorphisms, &composites)?;
    Ok(homomorphism.cache_is_right_minimal(result))
}

/// Translate QPA's annihilator-radical `IsLeftMinimal(f)` property method.
///
/// This is the source-dual ordered computation using `f * x` for the
/// `End(Range(f))` basis. The computed property is cached on the map.
///
/// # Errors
///
/// Propagates Hom, composition, exact nullspace, endomorphism conversion, and
/// radical computation failures.
pub fn is_left_minimal<F: QpaCoefficientField>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<bool, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_is_left_minimal() {
        return Ok(cached);
    }
    let range = homomorphism.range();
    let endomorphisms = hom_over_algebra(range, range)?
        .ok_or(PathAlgebraModuleError::InternalMinimalityHomRejected)?;
    let composites = endomorphisms
        .iter()
        .map(|endomorphism| homomorphism.compose(endomorphism))
        .collect::<Result<Vec<_>, _>>()?;
    let result = annihilator_is_radical(range, &endomorphisms, &composites)?;
    Ok(homomorphism.cache_is_left_minimal(result))
}

fn annihilator_is_radical<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    endomorphisms: &[PathAlgebraModuleHomomorphism<F>],
    composites: &[PathAlgebraModuleHomomorphism<F>],
) -> Result<bool, PathAlgebraModuleError> {
    let first = composites
        .first()
        .ok_or(PathAlgebraModuleError::EmptyMinimalityEndomorphismBasis)?;
    let width = first.flat_vertex_matrix_entries().len();
    let composite_rows = DenseMatrix::from_flat(
        composites.len(),
        width,
        composites
            .iter()
            .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
            .collect(),
        module.right_acting_algebra().field().clone(),
    )?;
    let annihilator = composite_rows.nullspace_rows();
    let endomorphism_algebra = end_over_algebra(module)?;
    let radical = endomorphism_algebra.radical_basis_matrices()?;
    let matrix_width = module
        .dimension()
        .checked_mul(module.dimension())
        .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
    let radical_rows = DenseMatrix::from_flat(
        radical.len(),
        matrix_width,
        radical.iter().flat_map(flatten_matrix).collect(),
        module.right_acting_algebra().field().clone(),
    )?;

    for row in 0..annihilator.nrows() {
        let annihilating_map =
            linear_combination_of_homomorphisms(endomorphisms, annihilator.row(row))?;
        let matrix = from_hom_mm_to_end_m(&annihilating_map)?;
        if radical_rows.solution_row(matrix.as_flat())?.is_none() {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Translate QPA's one-step `MoreRightMinimalVersion(f)` operation.
///
/// Hom bases, the nested loop order, the Fitting exponent, and the returned
/// kernel/image restrictions follow `modulehom.gi` literally. `Unchanged`
/// contains the identical entered map, matching QPA's heterogeneous map-or-list
/// result.
///
/// # Errors
///
/// Propagates kernel, Hom-space, composition, power, and image construction
/// failures.
pub fn more_right_minimal_version<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<MoreMinimalVersionResult<F>, PathAlgebraModuleError> {
    let source = homomorphism.source();
    let kernel = kernel_inclusion(homomorphism)?;
    let kernel_module = kernel.source();
    let hom_source_kernel = hom_over_algebra(source, kernel_module)?
        .ok_or(PathAlgebraModuleError::InternalMinimalHomRejected)?;
    if hom_source_kernel.is_empty() {
        return Ok(MoreMinimalVersionResult::Unchanged(homomorphism.clone()));
    }
    let hom_kernel_kernel = hom_over_algebra(kernel_module, kernel_module)?
        .ok_or(PathAlgebraModuleError::InternalMinimalHomRejected)?;
    let exponent = kernel_module
        .dimension_vector()
        .iter()
        .chain(source.dimension_vector())
        .copied()
        .max()
        .unwrap_or(0);

    for source_to_kernel in &hom_source_kernel {
        for kernel_endomorphism in &hom_kernel_kernel {
            let kernel_cycle = kernel
                .compose(source_to_kernel)?
                .compose(kernel_endomorphism)?;
            let fitting_kernel = homomorphism_power(&kernel_cycle, exponent)?;
            if !fitting_kernel.is_zero() {
                let source_cycle = source_to_kernel
                    .compose(kernel_endomorphism)?
                    .compose(&kernel)?;
                let fitting_source = homomorphism_power(&source_cycle, exponent)?;
                return Ok(MoreMinimalVersionResult::Reduced {
                    minimal_part: kernel_inclusion(&fitting_source)?.compose(homomorphism)?,
                    removed_part: image_inclusion(&fitting_source)?.compose(homomorphism)?,
                });
            }
        }
    }
    Ok(MoreMinimalVersionResult::Unchanged(homomorphism.clone()))
}

/// Translate QPA's one-step `MoreLeftMinimalVersion(f)` operation.
///
/// This preserves QPA's cokernel Hom-basis order, its Fitting powers, and the
/// right inverse used to corestrict to the retained kernel summand.
///
/// # Errors
///
/// Propagates cokernel, Hom-space, composition, power, split-inverse, kernel,
/// and image construction failures. An absent internal Fitting retraction is
/// QPA's subsequent map-times-`false` error and is reported explicitly.
pub fn more_left_minimal_version<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<MoreMinimalVersionResult<F>, PathAlgebraModuleError> {
    let range = homomorphism.range();
    let cokernel = cokernel_projection(homomorphism)?;
    let cokernel_module = cokernel.range();
    let hom_cokernel_range = hom_over_algebra(cokernel_module, range)?
        .ok_or(PathAlgebraModuleError::InternalMinimalHomRejected)?;
    if hom_cokernel_range.is_empty() {
        return Ok(MoreMinimalVersionResult::Unchanged(homomorphism.clone()));
    }
    let hom_cokernel_cokernel = hom_over_algebra(cokernel_module, cokernel_module)?
        .ok_or(PathAlgebraModuleError::InternalMinimalHomRejected)?;
    let exponent = range
        .dimension_vector()
        .iter()
        .chain(cokernel_module.dimension_vector())
        .copied()
        .max()
        .unwrap_or(0);

    for cokernel_endomorphism in &hom_cokernel_cokernel {
        for cokernel_to_range in &hom_cokernel_range {
            let cokernel_cycle = cokernel_endomorphism
                .compose(cokernel_to_range)?
                .compose(&cokernel)?;
            let fitting_cokernel = homomorphism_power(&cokernel_cycle, exponent)?;
            if !fitting_cokernel.is_zero() {
                let range_cycle = cokernel
                    .compose(cokernel_endomorphism)?
                    .compose(cokernel_to_range)?;
                let fitting_range = homomorphism_power(&range_cycle, exponent)?;
                let retained_inclusion = kernel_inclusion(&fitting_range)?;
                let retraction = right_inverse_of_homomorphism(&retained_inclusion)?
                    .ok_or(PathAlgebraModuleError::InternalMinimalKernelNotSplit)?;
                return Ok(MoreMinimalVersionResult::Reduced {
                    minimal_part: homomorphism.compose(&retraction)?,
                    removed_part: homomorphism.compose(&image_projection(&fitting_range)?)?,
                });
            }
        }
    }
    Ok(MoreMinimalVersionResult::Unchanged(homomorphism.clone()))
}

/// Translate QPA's cached `RightMinimalVersion(f)` attribute.
///
/// One source summand is removed per successful one-step call and recorded in
/// order. The final map has QPA's `IsRightMinimal` property installed as true.
///
/// # Errors
///
/// Propagates [`more_right_minimal_version`].
pub fn right_minimal_version<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleMinimalVersion<F>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_right_minimal_version() {
        return Ok(cached);
    }
    let mut removed_modules = Vec::new();
    let mut minimal = homomorphism.clone();
    loop {
        match more_right_minimal_version(&minimal)? {
            MoreMinimalVersionResult::Unchanged(_) => break,
            MoreMinimalVersionResult::Reduced {
                minimal_part,
                removed_part,
            } => {
                removed_modules.push(removed_part.source().clone());
                minimal = minimal_part;
            }
        }
    }
    let _ = minimal.cache_is_right_minimal(true);
    Ok(homomorphism.cache_right_minimal_version(&minimal, &removed_modules))
}

/// Translate QPA's cached `LeftMinimalVersion(f)` attribute.
///
/// One range summand is removed per successful one-step call and recorded in
/// order. The final map has QPA's `IsLeftMinimal` property installed as true.
///
/// # Errors
///
/// Propagates [`more_left_minimal_version`].
pub fn left_minimal_version<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleMinimalVersion<F>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_left_minimal_version() {
        return Ok(cached);
    }
    let mut removed_modules = Vec::new();
    let mut minimal = homomorphism.clone();
    loop {
        match more_left_minimal_version(&minimal)? {
            MoreMinimalVersionResult::Unchanged(_) => break,
            MoreMinimalVersionResult::Reduced {
                minimal_part,
                removed_part,
            } => {
                removed_modules.push(removed_part.range().clone());
                minimal = minimal_part;
            }
        }
    }
    let _ = minimal.cache_is_left_minimal(true);
    Ok(homomorphism.cache_left_minimal_version(&minimal, &removed_modules))
}

fn homomorphism_power<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    exponent: usize,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if exponent == 0 {
        return identity_module_mapping(homomorphism.source().clone());
    }
    let mut power = homomorphism.clone();
    for _ in 1..exponent {
        power = power.compose(homomorphism)?;
    }
    Ok(power)
}

fn solve_map_linear_combination<F: Field>(
    basis: &[PathAlgebraModuleHomomorphism<F>],
    target: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<Vec<F::Element>>, PathAlgebraModuleError> {
    let field = target.source().right_acting_algebra().field().clone();
    let row_length = target.flat_vertex_matrix_entries().len();
    let row_count = basis.len();
    let entry_count = row_count
        .checked_mul(row_length)
        .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
    let mut entries = Vec::with_capacity(entry_count);
    for homomorphism in basis {
        entries.extend(homomorphism.flat_vertex_matrix_entries());
    }
    let rows = DenseMatrix::from_flat(row_count, row_length, entries, field)?;
    Ok(rows.solution_row(&target.flat_vertex_matrix_entries())?)
}

fn linear_combination_of_homomorphisms<F: Field>(
    basis: &[PathAlgebraModuleHomomorphism<F>],
    coefficients: &[F::Element],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if basis.len() != coefficients.len() {
        return Err(PathAlgebraModuleError::SplittingCoefficientCount {
            actual: coefficients.len(),
            expected: basis.len(),
        });
    }
    let first = basis
        .first()
        .ok_or(PathAlgebraModuleError::EmptySplittingHomBasis)?;
    let field = first.source().right_acting_algebra().field();
    let mut combination = zero_module_mapping(first.source().clone(), first.range().clone())?;
    for (homomorphism, coefficient) in basis.iter().zip(coefficients) {
        if !field.is_zero(coefficient) {
            combination = combination.add(&homomorphism.scale(coefficient)?)?;
        }
    }
    Ok(combination)
}

fn zero_module_for_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    right_module_over_path_algebra(
        algebra.clone(),
        vec![0; algebra.quiver().number_of_vertices()],
        &[],
    )?
    .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)
}

/// Translate QPA's cached `KernelInclusion(f)` attribute.
///
/// Component kernels use the exact GAP `NullspaceMat` row convention. Their
/// arrow actions are reconstructed by solving in those ordered bases, and the
/// zero-kernel branch returns QPA's physical zero mapping.
///
/// # Errors
///
/// Propagates exact matrix operations and checked module/map construction.
pub fn kernel_inclusion<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_kernel_inclusion() {
        return Ok(cached);
    }

    let source = homomorphism.source();
    let algebra = source.right_acting_algebra().clone();
    let field = algebra.field().clone();
    let kernel_bases = source
        .dimension_vector()
        .iter()
        .enumerate()
        .map(|(vertex, &dimension)| {
            if dimension == 0 {
                DenseMatrix::from_flat(0, 0, Vec::new(), field.clone())
            } else {
                Ok(logical_vertex_map(homomorphism, vertex).nullspace_rows())
            }
        })
        .collect::<Result<Vec<_>, MatrixError>>()?;
    let dimensions = kernel_bases
        .iter()
        .map(DenseMatrix::nrows)
        .collect::<Vec<_>>();

    let inclusion = if dimensions.iter().all(|&dimension| dimension == 0) {
        homomorphism.set_is_injective(true);
        zero_module_mapping(zero_module_for_algebra(&algebra)?, source.clone())?
    } else {
        let coordinate_solvers = kernel_bases
            .iter()
            .map(|_| std::cell::OnceCell::new())
            .collect::<Vec<_>>();
        let mut arrow_matrices = Vec::new();
        for arrow in algebra.quiver().arrows() {
            let source_vertex = arrow.source().index().zero_based();
            let target_vertex = arrow.target().index().zero_based();
            if dimensions[source_vertex] == 0 || dimensions[target_vertex] == 0 {
                continue;
            }
            let source_action = logical_submatrix(
                &source.arrow_matrices()[arrow.index().zero_based()],
                source.dimension_vector()[source_vertex],
                source.dimension_vector()[target_vertex],
            )?;
            let acted = kernel_bases[source_vertex].multiply(&source_action)?;
            let mut entries =
                Vec::with_capacity(dimensions[source_vertex] * dimensions[target_vertex]);
            let solver = coordinate_solvers[target_vertex]
                .get_or_init(|| kernel_bases[target_vertex].semi_echelon_transformation());
            for row in 0..acted.nrows() {
                let coordinates = solver
                    .solution_row(acted.row(row))?
                    .ok_or(PathAlgebraModuleError::InternalKernelCoordinates)?;
                entries.extend(coordinates);
            }
            arrow_matrices.push(NamedArrowMatrix::new(
                arrow.name(),
                DenseMatrix::from_flat(
                    dimensions[source_vertex],
                    dimensions[target_vertex],
                    entries,
                    field.clone(),
                )?,
            ));
        }
        let kernel = right_module_over_path_algebra(algebra, dimensions.clone(), &arrow_matrices)?
            .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)?;
        let inclusion_matrices = kernel_bases
            .iter()
            .enumerate()
            .map(|(vertex, basis)| {
                if dimensions[vertex] == 0 {
                    DenseMatrix::zero(
                        1,
                        physical_dimension(source.dimension_vector()[vertex]),
                        field.clone(),
                    )
                } else {
                    basis.clone()
                }
            })
            .collect();
        right_module_hom_over_algebra(kernel, source.clone(), inclusion_matrices)?
    };
    inclusion.set_kernel_of_what(homomorphism);
    inclusion.set_is_injective(true);
    Ok(homomorphism.cache_kernel_inclusion(&inclusion))
}

/// Return the source of QPA's canonical kernel inclusion.
///
/// # Errors
///
/// Propagates [`kernel_inclusion`].
pub fn kernel_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(kernel_inclusion(homomorphism)?.source().clone())
}

/// Translate QPA's cached image projection/inclusion factorization.
///
/// At each vertex, QPA takes the fully reduced canonical basis of the row
/// image, reconstructs the induced range action, and solves each source-basis
/// image in that basis. The returned pair is projection first, inclusion
/// second, exactly as in `ImageProjectionInclusion`.
///
/// # Errors
///
/// Propagates exact row-space/coordinate operations and checked module/map
/// construction.
#[allow(clippy::too_many_lines)]
pub fn image_projection_inclusion<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<
    (
        PathAlgebraModuleHomomorphism<F>,
        PathAlgebraModuleHomomorphism<F>,
    ),
    PathAlgebraModuleError,
> {
    if let Some(cached) = homomorphism.cached_image_projection_inclusion() {
        return Ok(cached);
    }

    let source = homomorphism.source();
    let range = homomorphism.range();
    let algebra = source.right_acting_algebra().clone();
    let field = algebra.field().clone();
    let image_bases = source
        .dimension_vector()
        .iter()
        .enumerate()
        .map(|(vertex, _)| logical_vertex_map(homomorphism, vertex).canonical_row_space_basis())
        .collect::<Vec<_>>();
    let dimensions = image_bases
        .iter()
        .map(DenseMatrix::nrows)
        .collect::<Vec<_>>();

    if dimensions.iter().all(|&dimension| dimension == 0) {
        let zero = zero_module_for_algebra(&algebra)?;
        let projection = zero_module_mapping(source.clone(), zero.clone())?;
        let inclusion = zero_module_mapping(zero, range.clone())?;
        return Ok(homomorphism.cache_image_projection_inclusion(&projection, &inclusion));
    }

    let coordinate_solvers = image_bases
        .iter()
        .map(|_| std::cell::OnceCell::new())
        .collect::<Vec<_>>();
    let mut arrow_matrices = Vec::new();
    for arrow in algebra.quiver().arrows() {
        let source_vertex = arrow.source().index().zero_based();
        let target_vertex = arrow.target().index().zero_based();
        if dimensions[source_vertex] == 0 || dimensions[target_vertex] == 0 {
            continue;
        }
        let range_action = logical_submatrix(
            &range.arrow_matrices()[arrow.index().zero_based()],
            range.dimension_vector()[source_vertex],
            range.dimension_vector()[target_vertex],
        )?;
        let acted = image_bases[source_vertex].multiply(&range_action)?;
        let mut entries = Vec::with_capacity(dimensions[source_vertex] * dimensions[target_vertex]);
        let solver = coordinate_solvers[target_vertex]
            .get_or_init(|| image_bases[target_vertex].semi_echelon_transformation());
        for row in 0..acted.nrows() {
            let coordinates = solver
                .solution_row(acted.row(row))?
                .ok_or(PathAlgebraModuleError::InternalImageCoordinates)?;
            entries.extend(coordinates);
        }
        arrow_matrices.push(NamedArrowMatrix::new(
            arrow.name(),
            DenseMatrix::from_flat(
                dimensions[source_vertex],
                dimensions[target_vertex],
                entries,
                field.clone(),
            )?,
        ));
    }
    let image = right_module_over_path_algebra(algebra, dimensions.clone(), &arrow_matrices)?
        .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)?;

    let inclusion_matrices = image_bases
        .iter()
        .enumerate()
        .map(|(vertex, basis)| {
            if dimensions[vertex] == 0 {
                DenseMatrix::zero(
                    1,
                    physical_dimension(range.dimension_vector()[vertex]),
                    field.clone(),
                )
            } else {
                basis.clone()
            }
        })
        .collect();
    let inclusion =
        right_module_hom_over_algebra(image.clone(), range.clone(), inclusion_matrices)?;
    inclusion.set_image_of_what(homomorphism);
    inclusion.set_is_injective(true);

    let projection_matrices = source
        .dimension_vector()
        .iter()
        .enumerate()
        .map(
            |(vertex, &source_dimension)| -> Result<DenseMatrix<F>, PathAlgebraModuleError> {
                if source_dimension == 0 || dimensions[vertex] == 0 {
                    return Ok(DenseMatrix::zero(
                        physical_dimension(source_dimension),
                        physical_dimension(dimensions[vertex]),
                        field.clone(),
                    ));
                }
                let logical_map = logical_vertex_map(homomorphism, vertex);
                let mut entries = Vec::with_capacity(source_dimension * dimensions[vertex]);
                let solver = coordinate_solvers[vertex]
                    .get_or_init(|| image_bases[vertex].semi_echelon_transformation());
                for row in 0..source_dimension {
                    let coordinates = solver
                        .solution_row(logical_map.row(row))?
                        .ok_or(PathAlgebraModuleError::InternalImageCoordinates)?;
                    entries.extend(coordinates);
                }
                Ok(DenseMatrix::from_flat(
                    source_dimension,
                    dimensions[vertex],
                    entries,
                    field.clone(),
                )?)
            },
        )
        .collect::<Result<Vec<_>, _>>()?;
    let projection = right_module_hom_over_algebra(source.clone(), image, projection_matrices)?;
    projection.set_image_of_what(homomorphism);
    projection.set_is_surjective(true);

    Ok(homomorphism.cache_image_projection_inclusion(&projection, &inclusion))
}

/// Return QPA's canonical projection onto the image.
///
/// # Errors
///
/// Propagates [`image_projection_inclusion`].
pub fn image_projection<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    Ok(image_projection_inclusion(homomorphism)?.0)
}

/// Return QPA's canonical inclusion of the image.
///
/// # Errors
///
/// Propagates [`image_projection_inclusion`].
pub fn image_inclusion<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    Ok(image_projection_inclusion(homomorphism)?.1)
}

/// Return QPA's image module object.
///
/// # Errors
///
/// Propagates [`image_projection_inclusion`].
pub fn image_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(image_projection(homomorphism)?.range().clone())
}

struct RestrictedModuleData<F: Field> {
    module: PathAlgebraModule<F>,
    original_basis_by_vertex: Vec<Vec<PathAlgebraModuleElement<F>>>,
    original_basis_rows: Vec<DenseMatrix<F>>,
}

/// Translate QPA's `RestrictionViaAlgebraHomomorphism(f, M)` operation.
///
/// For `f: A -> B`, QPA applies the images of the vertex idempotents of `A`
/// to the standard basis of the `B`-module, takes each resulting canonical
/// vector-space basis, and computes the `A`-arrow action in those bases.
///
/// # Errors
///
/// Requires the module to act over the range of `f` and propagates algebra-map
/// evaluation, exact basis, action, and checked restricted-module failures.
pub fn restriction_via_algebra_homomorphism<F: Field>(
    homomorphism: &QuiverAlgebraHomomorphism<F>,
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, RestrictionViaAlgebraHomomorphismError> {
    Ok(restricted_module_data(homomorphism, module)?.module)
}

/// Translate QPA's `RestrictionViaAlgebraHomomorphismMap(f, h)` operation.
///
/// The source and range are independently restricted by QPA's module method;
/// the entered map is then expressed between their vertex-idempotent bases.
///
/// # Errors
///
/// Requires `h` to be over the range of `f` and propagates both module
/// restrictions, map application, exact coordinates, and checked-map errors.
pub fn restriction_via_algebra_homomorphism_map<F: Field>(
    algebra_homomorphism: &QuiverAlgebraHomomorphism<F>,
    module_homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, RestrictionViaAlgebraHomomorphismError> {
    if module_homomorphism.source().right_acting_algebra() != algebra_homomorphism.range() {
        return Err(RestrictionViaAlgebraHomomorphismError::HomomorphismNotOverRangeAlgebra);
    }
    let source = restricted_module_data(algebra_homomorphism, module_homomorphism.source())?;
    let range = restricted_module_data(algebra_homomorphism, module_homomorphism.range())?;
    let field = algebra_homomorphism.source().field().clone();
    let matrices = source
        .original_basis_by_vertex
        .iter()
        .enumerate()
        .map(
            |(vertex, basis)| -> Result<DenseMatrix<F>, RestrictionViaAlgebraHomomorphismError> {
                let source_dimension = basis.len();
                let range_dimension = range.original_basis_by_vertex[vertex].len();
                if source_dimension == 0 || range_dimension == 0 {
                    return Ok(DenseMatrix::zero(
                        physical_dimension(source_dimension),
                        physical_dimension(range_dimension),
                        field.clone(),
                    ));
                }
                let mut entries = Vec::with_capacity(source_dimension * range_dimension);
                for element in basis {
                    let image = module_homomorphism.image(element)?;
                    entries.extend(
                        range.original_basis_rows[vertex]
                            .solution_row(&image.basis_coefficients())?
                            .ok_or(
                                RestrictionViaAlgebraHomomorphismError::MapImageOutsideRestrictedVertex,
                            )?,
                    );
                }
                Ok(DenseMatrix::from_flat(
                    source_dimension,
                    range_dimension,
                    entries,
                    field.clone(),
                )?)
            },
        )
        .collect::<Result<Vec<_>, _>>()?;
    Ok(right_module_hom_over_algebra(
        source.module,
        range.module,
        matrices,
    )?)
}

fn restricted_module_data<F: Field>(
    homomorphism: &QuiverAlgebraHomomorphism<F>,
    module: &PathAlgebraModule<F>,
) -> Result<RestrictedModuleData<F>, RestrictionViaAlgebraHomomorphismError> {
    if module.right_acting_algebra() != homomorphism.range() {
        return Err(RestrictionViaAlgebraHomomorphismError::ModuleNotOverRangeAlgebra);
    }
    let source = homomorphism.source();
    let field = source.field().clone();
    let vertex_count = source.quiver().number_of_vertices();
    let source_generators = source.generators_with_one()?;
    let module_basis = module.basis_elements();
    let mut basis_by_vertex = Vec::with_capacity(vertex_count);
    let mut basis_rows = Vec::with_capacity(vertex_count);
    for vertex in &source_generators[..vertex_count] {
        let vertex_image = homomorphism.image(vertex)?;
        let spanning_rows = module_basis
            .iter()
            .map(|basis| basis.act(&vertex_image))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .filter(|element| !element.is_zero())
            .map(|element| element.basis_coefficients())
            .collect::<Vec<_>>();
        let basis = DenseMatrix::from_flat(
            spanning_rows.len(),
            module.dimension(),
            spanning_rows.into_iter().flatten().collect(),
            field.clone(),
        )?
        .canonical_row_space_basis();
        let elements = (0..basis.nrows())
            .map(|row| PathAlgebraModuleElement::from_logical_flat(module.clone(), basis.row(row)))
            .collect::<Result<Vec<_>, _>>()?;
        basis_rows.push(basis);
        basis_by_vertex.push(elements);
    }

    let dimensions = basis_by_vertex.iter().map(Vec::len).collect::<Vec<_>>();
    let mut arrow_matrices = Vec::new();
    for (arrow, arrow_generator) in source
        .quiver()
        .arrows()
        .zip(source_generators.iter().skip(vertex_count))
    {
        let arrow_source = arrow.source().index().zero_based();
        let arrow_target = arrow.target().index().zero_based();
        if dimensions[arrow_source] == 0 || dimensions[arrow_target] == 0 {
            continue;
        }
        let arrow_image = homomorphism.image(arrow_generator)?;
        let mut entries =
            Vec::with_capacity(dimensions[arrow_source].saturating_mul(dimensions[arrow_target]));
        for element in &basis_by_vertex[arrow_source] {
            let image = element.act(&arrow_image)?;
            entries.extend(
                basis_rows[arrow_target]
                    .solution_row(&image.basis_coefficients())?
                    .ok_or(
                        RestrictionViaAlgebraHomomorphismError::ArrowImageOutsideRestrictedVertex,
                    )?,
            );
        }
        arrow_matrices.push(NamedArrowMatrix::new(
            arrow.name(),
            DenseMatrix::from_flat(
                dimensions[arrow_source],
                dimensions[arrow_target],
                entries,
                field.clone(),
            )?,
        ));
    }
    let restricted = right_module_over_path_algebra(source.clone(), dimensions, &arrow_matrices)?
        .ok_or(RestrictionViaAlgebraHomomorphismError::RestrictedModuleRejected)?;
    Ok(RestrictedModuleData {
        module: restricted,
        original_basis_by_vertex: basis_by_vertex,
        original_basis_rows: basis_rows,
    })
}

/// Translate QPA's linear-equation `HomOverAlgebra(M, N)` method.
///
/// Unknown vertex-map entries are ordered vertex-major, then row-major. The
/// equation columns follow supported source vertices, their outgoing arrows,
/// source rows, and target coordinates exactly as in `modulehom.gi`.
/// `Ok(None)` is QPA's printed-then-`fail` result for modules over different
/// algebras.
///
/// # Errors
///
/// Reports address-size overflow and propagates exact matrix or checked-map
/// construction failures.
pub fn hom_over_algebra<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Option<Vec<PathAlgebraModuleHomomorphism<F>>>, PathAlgebraModuleError> {
    let Some(equations) = hom_equations(source, range)? else {
        return Ok(None);
    };
    if equations.ncols() == 0 {
        // Construct elementary maps one at a time, without retaining a dense
        // identity nullspace for an unconstrained Hom space.
        return (0..equations.nrows())
            .map(|index| {
                let mut entries = vec![equations.field().zero(); equations.nrows()];
                entries[index] = equations.field().one();
                module_homomorphism_from_logical_flat(source, range, &entries)
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Some);
    }
    let nullspace = equations.nullspace_rows();
    let homomorphisms = (0..nullspace.nrows())
        .map(|row| module_homomorphism_from_logical_flat(source, range, nullspace.row(row)))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Some(homomorphisms))
}

/// Compute the dimension of `Hom(source, range)` without constructing basis maps.
///
/// Uses the same equations as [`hom_over_algebra`], returning their left nullity.
/// Returns `Ok(None)` for different acting algebras, including zero modules.
///
/// # Errors
///
/// Propagates equation-size overflow and exact matrix construction failures.
pub fn hom_dimension_over_algebra<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Option<usize>, PathAlgebraModuleError> {
    Ok(hom_equations(source, range)?.map(|equations| {
        if equations.ncols() == 0 {
            equations.nrows()
        } else {
            equations.nrows() - equations.rank()
        }
    }))
}

#[allow(clippy::too_many_lines)]
fn hom_equations<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Option<crate::matrix::sparse::SparseMatrix<F>>, PathAlgebraModuleError> {
    if source.right_acting_algebra() != range.right_acting_algebra() {
        return Ok(None);
    }
    let algebra = source.right_acting_algebra();
    let field = algebra.field().clone();
    let source_dimensions = source.dimension_vector();
    let range_dimensions = range.dimension_vector();
    let vertex_count = source_dimensions.len();
    let source_support = source_dimensions
        .iter()
        .map(|&dimension| dimension != 0)
        .collect::<Vec<_>>();
    let range_support = range_dimensions
        .iter()
        .map(|&dimension| dimension != 0)
        .collect::<Vec<_>>();

    let mut row_count = 0_usize;
    let mut block_rows = vec![None; vertex_count];
    for vertex in 0..vertex_count {
        if source_support[vertex] {
            block_rows[vertex] = Some(row_count);
            row_count = row_count
                .checked_add(
                    source_dimensions[vertex]
                        .checked_mul(range_dimensions[vertex])
                        .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?,
                )
                .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        }
    }

    let arrows = algebra.quiver().arrows().collect::<Vec<_>>();
    let mut column_count = 0_usize;
    for vertex in 0..vertex_count {
        if !source_support[vertex] {
            continue;
        }
        for arrow in arrows
            .iter()
            .filter(|arrow| arrow.source().index().zero_based() == vertex)
        {
            let arrow_source = arrow.source().index().zero_based();
            let arrow_target = arrow.target().index().zero_based();
            if range_support[arrow_target]
                && (range_support[arrow_source] || source_support[arrow_target])
            {
                column_count = column_count
                    .checked_add(
                        source_dimensions[arrow_source]
                            .checked_mul(range_dimensions[arrow_target])
                            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?,
                    )
                    .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
            }
        }
    }

    // Keep QPA's unknown/equation order, but do not allocate their usually
    // almost entirely zero Cartesian product. Sparse elimination preserves
    // the dense algorithm's ordered nullspace representatives exactly.
    let mut equations = crate::matrix::sparse::SparseMatrix::zero(row_count, column_count, field.clone());
    if row_count != 0 && column_count != 0 {
        let mut previous_column = 0_usize;
        for vertex in 0..vertex_count {
            if !source_support[vertex] {
                continue;
            }
            for arrow in arrows
                .iter()
                .filter(|arrow| arrow.source().index().zero_based() == vertex)
            {
                let arrow_source = arrow.source().index().zero_based();
                let arrow_target = arrow.target().index().zero_based();
                if !range_support[arrow_target]
                    || !(range_support[arrow_source] || source_support[arrow_target])
                {
                    continue;
                }
                let source_action = logical_submatrix(
                    &source.arrow_matrices()[arrow.index().zero_based()],
                    source_dimensions[arrow_source],
                    source_dimensions[arrow_target],
                )?;
                let range_action = logical_submatrix(
                    &range.arrow_matrices()[arrow.index().zero_based()],
                    range_dimensions[arrow_source],
                    range_dimensions[arrow_target],
                )?;
                for source_row in 0..source_dimensions[arrow_source] {
                    let column_start =
                        previous_column + source_row * range_dimensions[arrow_target];
                    if range_support[arrow_source] {
                        let row_start = block_rows[arrow_source].ok_or(
                            PathAlgebraModuleError::MissingHomUnknownBlock {
                                vertex: arrow_source + 1,
                            },
                        )? + source_row * range_dimensions[arrow_source];
                        for range_row in 0..range_dimensions[arrow_source] {
                            for target_column in 0..range_dimensions[arrow_target] {
                                equations.set(row_start + range_row, column_start + target_column,
                                    range_action[(range_row, target_column)].clone());
                            }
                        }
                    }
                    if source_support[arrow_target] {
                        let target_row_start = block_rows[arrow_target].ok_or(
                            PathAlgebraModuleError::MissingHomUnknownBlock {
                                vertex: arrow_target + 1,
                            },
                        )?;
                        for target_row in 0..source_dimensions[arrow_target] {
                            let negative = field.neg(&source_action[(source_row, target_row)]);
                            for target_column in 0..range_dimensions[arrow_target] {
                                equations.add(
                                    target_row_start
                                        + target_row * range_dimensions[arrow_target]
                                        + target_column,
                                    column_start + target_column,
                                    &negative,
                                );
                            }
                        }
                    }
                }
                previous_column += source_dimensions[arrow_source] * range_dimensions[arrow_target];
            }
        }
    }
    Ok(Some(equations))
}

/// Translate QPA's `HomOverAlgebraWithBasisFunction` result pair.
///
/// `Ok(None)` is the underlying `HomOverAlgebra` failure for different acting
/// algebras.
///
/// # Errors
///
/// Propagates [`hom_over_algebra`] and exact flattened-basis construction.
pub fn hom_over_algebra_with_basis_function<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModuleHomomorphismBasis<F>>, PathAlgebraModuleError> {
    let Some(basis) = hom_over_algebra(source, range)? else {
        return Ok(None);
    };
    let flattened_basis = if let Some(first) = basis.first() {
        let column_count = flatten_module_homomorphism(first).len();
        Some(DenseMatrix::from_flat(
            basis.len(),
            column_count,
            basis.iter().flat_map(flatten_module_homomorphism).collect(),
            source.right_acting_algebra().field().clone(),
        )?)
    } else {
        None
    };
    Ok(Some(PathAlgebraModuleHomomorphismBasis {
        basis,
        flattened_basis,
        field: source.right_acting_algebra().field().clone(),
    }))
}

/// Translate QPA's cached `EndOverAlgebra(M)` matrix-algebra construction.
///
/// QPA embeds each vertex map as a diagonal block in the underlying
/// vertex-major vector space, forms the direct sum of full matrix algebras,
/// and then asks GAP for the basis of the generated subalgebra with one. The
/// latter is GAP's insertion-ordered semi-echelon basis of the flattened block
/// matrices.
///
/// # Errors
///
/// Propagates [`hom_over_algebra`], address-size overflow, and exact matrix
/// construction failures.
pub fn end_over_algebra<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<ModuleEndomorphismAlgebra<F>, PathAlgebraModuleError> {
    if let Some(cached) = module.cached_endomorphism_algebra() {
        return Ok(cached);
    }
    if module.dimension() == 0 {
        return Err(PathAlgebraModuleError::EmptyEndomorphismDirectSum);
    }
    let homomorphisms = hom_over_algebra(module, module)?
        .ok_or(PathAlgebraModuleError::EndomorphismAlgebraMismatch)?;
    let dimension = module.dimension();
    let flattened_dimension = dimension
        .checked_mul(dimension)
        .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
    let field = module.right_acting_algebra().field().clone();
    let embedded = homomorphisms
        .iter()
        .map(from_hom_mm_to_end_m)
        .collect::<Result<Vec<_>, _>>()?;
    let generator_rows = DenseMatrix::from_flat(
        embedded.len(),
        flattened_dimension,
        embedded.iter().flat_map(flatten_matrix).collect(),
        field.clone(),
    )?;
    let flattened_basis = generator_rows.semi_echelon_basis();
    let basis_matrices = (0..flattened_basis.nrows())
        .map(|row| {
            DenseMatrix::from_flat(
                dimension,
                dimension,
                flattened_basis.row(row).to_vec(),
                field.clone(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let algebra = ModuleEndomorphismAlgebra(Arc::new(ModuleEndomorphismAlgebraData {
        field: field.clone(),
        basis_matrices,
        flattened_basis,
        one: DenseMatrix::identity(dimension, field),
    }));
    Ok(module.cache_endomorphism_algebra(&algebra))
}

/// Translate QPA's `NumberOfNonIsoDirSummands(M)` operation.
///
/// The source constructs `End(M)`, quotients by its algebra radical, obtains
/// the central-idempotent direct-sum decomposition, and returns the number and
/// dimensions of those semisimple blocks. This implementation follows the
/// same structure-constant route on the currently supported split prime-field
/// coefficient tier.
///
/// # Errors
///
/// Propagates endomorphism-algebra, radical, quotient, central-idempotent, and
/// block construction failures. In particular, the present central splitting
/// backend reports a typed finite-field error for rational coefficients.
pub fn number_of_nonisomorphic_direct_summands<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<NonisomorphicDirectSummandData, PathAlgebraModuleError> {
    let endomorphism_algebra = end_over_algebra(module)?.structure_constants_algebra()?;
    let radical = endomorphism_algebra.radical()?;
    let semisimple = endomorphism_algebra.quotient_by_ideal_basis(radical.basis().clone())?;
    let blocks = semisimple.algebra().blocks()?;
    Ok(NonisomorphicDirectSummandData {
        count: blocks.len(),
        block_dimensions: blocks
            .iter()
            .map(|block| block.algebra().dimension())
            .collect(),
    })
}

/// Represent QPA's declared but uninstalled `LiftIdempotents(M)` operation.
///
/// # Errors
///
/// Pinned QPA declares this operation for matrix modules but installs no
/// method, so every ordinary call reaches the missing-method boundary.
pub const fn lift_idempotents<F: Field>(
    _module: &PathAlgebraModule<F>,
) -> Result<(), PathAlgebraModuleError> {
    Err(PathAlgebraModuleError::UninstalledLiftIdempotents)
}

/// Translate QPA's `FromHomMMToEndM(f)` block-diagonal embedding.
///
/// # Errors
///
/// The source and range must be structurally equal modules; exact matrix
/// construction failures propagate.
pub fn from_hom_mm_to_end_m<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<DenseMatrix<F>, PathAlgebraModuleError> {
    if homomorphism.source() != homomorphism.range() {
        return Err(PathAlgebraModuleError::NonEndomorphismConversion);
    }
    let dimensions = homomorphism.source().dimension_vector();
    let total = homomorphism.source().dimension();
    let field = homomorphism.source().right_acting_algebra().field().clone();
    let offsets = dimension_offsets(dimensions);
    let mut matrix = DenseMatrix::zero(total, total, field);
    for (vertex, &dimension) in dimensions.iter().enumerate() {
        if dimension == 0 {
            continue;
        }
        let block = logical_vertex_map(homomorphism, vertex);
        for row in 0..dimension {
            for column in 0..dimension {
                matrix[(offsets[vertex] + row, offsets[vertex] + column)] =
                    block[(row, column)].clone();
            }
        }
    }
    Ok(matrix)
}

/// Translate QPA's `FromEndMToHomMM(M, mat)` diagonal-block extraction.
///
/// Entries outside the vertex-diagonal blocks are ignored, just as in QPA.
/// A matrix larger than the underlying module dimension is therefore
/// accepted, while a smaller or wrong-field matrix is rejected explicitly.
///
/// # Errors
///
/// Propagates field/size validation and QPA's checked module-map constructor.
pub fn from_end_m_to_hom_mm<F: Field>(
    module: &PathAlgebraModule<F>,
    matrix: &DenseMatrix<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let field = module.right_acting_algebra().field().clone();
    if matrix.field() != &field {
        return Err(PathAlgebraModuleError::EndomorphismMatrixField);
    }
    let total = module.dimension();
    if matrix.nrows() < total || matrix.ncols() < total {
        return Err(PathAlgebraModuleError::EndomorphismMatrixTooSmall {
            actual: (matrix.nrows(), matrix.ncols()),
            minimum: (total, total),
        });
    }
    let offsets = dimension_offsets(module.dimension_vector());
    let matrices = module
        .dimension_vector()
        .iter()
        .enumerate()
        .map(|(vertex, &dimension)| {
            if dimension == 0 {
                return Ok(DenseMatrix::zero(1, 1, field.clone()));
            }
            let offset = offsets[vertex];
            let mut entries = Vec::with_capacity(dimension * dimension);
            for row in offset..offset + dimension {
                entries.extend_from_slice(&matrix.row(row)[offset..offset + dimension]);
            }
            Ok(DenseMatrix::from_flat(
                dimension,
                dimension,
                entries,
                field.clone(),
            )?)
        })
        .collect::<Result<Vec<_>, PathAlgebraModuleError>>()?;
    right_module_hom_over_algebra(module.clone(), module.clone(), matrices)
}

/// Construct QPA's natural quotient homomorphism by the two-sided ideal
/// generated by represented endomorphism matrices.
///
/// Generator order is retained through GAP's insertion-ordered semi-echelon
/// basis algorithm. The ideal is closed under left and right multiplication
/// by the parent basis before the natural row-space quotient is formed.
///
/// # Errors
///
/// Every generator and generated product must lie in the parent algebra;
/// exact multiplication, coordinate, ideal-closure, and quotient failures
/// propagate.
pub fn natural_endomorphism_quotient_by_ideal<F: Field>(
    parent: &ModuleEndomorphismAlgebra<F>,
    generators: &[DenseMatrix<F>],
) -> Result<ModuleEndomorphismQuotientHomomorphism<F>, PathAlgebraModuleError> {
    let ideal_coordinates = endomorphism_ideal_basis_coordinates(parent, generators)?;
    let ideal_basis_matrices = (0..ideal_coordinates.nrows())
        .map(|row| parent.matrix_from_coefficients(ideal_coordinates.row(row)))
        .collect::<Result<Vec<_>, _>>()?;
    let quotient = natural_quotient_by_row_space(&ideal_coordinates)?;
    let range =
        ModuleEndomorphismQuotientAlgebra(Arc::new(ModuleEndomorphismQuotientAlgebraData {
            parent: parent.clone(),
            ideal_basis_matrices,
            projection: quotient.projection,
            lifts: quotient.lifts,
        }));
    Ok(ModuleEndomorphismQuotientHomomorphism {
        source: parent.clone(),
        range,
    })
}

fn endomorphism_ideal_basis_coordinates<F: Field>(
    parent: &ModuleEndomorphismAlgebra<F>,
    generators: &[DenseMatrix<F>],
) -> Result<DenseMatrix<F>, PathAlgebraModuleError> {
    let dimension = parent.dimension();
    let capacity = generators
        .len()
        .checked_mul(dimension)
        .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
    let mut generator_coordinates = Vec::with_capacity(capacity);
    for generator in generators {
        generator_coordinates.extend(
            parent
                .coefficients(generator)?
                .ok_or(PathAlgebraModuleError::EndomorphismIdealGeneratorOutsideAlgebra)?,
        );
    }
    let mut basis = DenseMatrix::from_flat(
        generators.len(),
        dimension,
        generator_coordinates,
        parent.field().clone(),
    )?
    .semi_echelon_basis();

    loop {
        let ideal_matrices = (0..basis.nrows())
            .map(|row| parent.matrix_from_coefficients(basis.row(row)))
            .collect::<Result<Vec<_>, _>>()?;
        let product_rows = ideal_matrices
            .len()
            .checked_mul(parent.basis_matrices().len())
            .and_then(|count| count.checked_mul(2))
            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        let expanded_rows = basis
            .nrows()
            .checked_add(product_rows)
            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        let expanded_capacity = expanded_rows
            .checked_mul(dimension)
            .ok_or(PathAlgebraModuleError::HomEquationSizeOverflow)?;
        let mut expanded_coordinates = Vec::with_capacity(expanded_capacity);
        expanded_coordinates.extend_from_slice(basis.as_flat());
        for ideal_element in &ideal_matrices {
            for algebra_element in parent.basis_matrices() {
                for product in [
                    parent.multiply(algebra_element, ideal_element)?,
                    parent.multiply(ideal_element, algebra_element)?,
                ] {
                    expanded_coordinates.extend(
                        parent.coefficients(&product)?.ok_or(
                            PathAlgebraModuleError::EndomorphismIdealProductOutsideAlgebra,
                        )?,
                    );
                }
            }
        }
        let expanded = DenseMatrix::from_flat(
            expanded_rows,
            dimension,
            expanded_coordinates,
            parent.field().clone(),
        )?
        .semi_echelon_basis();
        if expanded.nrows() == basis.nrows() {
            return Ok(basis);
        }
        basis = expanded;
    }
}

fn module_homomorphism_from_logical_flat<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
    entries: &[F::Element],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let field = source.right_acting_algebra().field().clone();
    let mut offset = 0_usize;
    let matrices = source
        .dimension_vector()
        .iter()
        .zip(range.dimension_vector())
        .map(
            |(&source_dimension, &range_dimension)| -> Result<DenseMatrix<F>, MatrixError> {
                if source_dimension == 0 || range_dimension == 0 {
                    return Ok(DenseMatrix::zero(
                        physical_dimension(source_dimension),
                        physical_dimension(range_dimension),
                        field.clone(),
                    ));
                }
                let length = source_dimension * range_dimension;
                let matrix = DenseMatrix::from_flat(
                    source_dimension,
                    range_dimension,
                    entries[offset..offset + length].to_vec(),
                    field.clone(),
                )?;
                offset += length;
                Ok(matrix)
            },
        )
        .collect::<Result<Vec<_>, _>>()?;
    right_module_hom_over_algebra(source.clone(), range.clone(), matrices)
}

fn flatten_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Vec<F::Element> {
    homomorphism.flat_vertex_matrix_entries()
}

fn flatten_matrix<F: Field>(matrix: &DenseMatrix<F>) -> Vec<F::Element> {
    (0..matrix.nrows())
        .flat_map(|row| matrix.row(row).iter().cloned())
        .collect()
}

#[allow(clippy::too_many_lines)]
fn hopf_matrix_module_hom<F: Field>(
    algebra: &QuiverAlgebra<F>,
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Vec<DenseMatrix<F>>, PathAlgebraModuleError> {
    if !algebra.is_identical(source.right_acting_algebra())
        || !algebra.is_identical(range.right_acting_algebra())
    {
        return Err(PathAlgebraModuleError::GenericHomActingAlgebraMismatch);
    }
    if source.dimension() == 0 || range.dimension() == 0 {
        return Err(PathAlgebraModuleError::ZeroGenericHomModuleFailure);
    }

    let field = algebra.field().clone();
    let source_dimensions = source.dimension_vector();
    let range_dimensions = range.dimension_vector();
    let vertex_count = source_dimensions.len();
    let mut column_count = 0_usize;
    let mut column_starts = Vec::with_capacity(vertex_count);
    for vertex in 0..vertex_count {
        column_starts.push(column_count);
        column_count = column_count
            .checked_add(
                source_dimensions[vertex]
                    .checked_mul(range_dimensions[vertex])
                    .ok_or(PathAlgebraModuleError::GenericHomEquationSizeOverflow)?,
            )
            .ok_or(PathAlgebraModuleError::GenericHomEquationSizeOverflow)?;
    }

    let mut row_count = 0_usize;
    for arrow in algebra.quiver().arrows() {
        let source_vertex = arrow.source().index().zero_based();
        let target_vertex = arrow.target().index().zero_based();
        row_count = row_count
            .checked_add(
                source_dimensions[source_vertex]
                    .checked_mul(range_dimensions[target_vertex])
                    .ok_or(PathAlgebraModuleError::GenericHomEquationSizeOverflow)?,
            )
            .ok_or(PathAlgebraModuleError::GenericHomEquationSizeOverflow)?;
    }
    if row_count == 0 {
        return Err(PathAlgebraModuleError::DegenerateGenericHomNullspaceFailure);
    }

    let mut equations = DenseMatrix::zero(row_count, column_count, field.clone());
    let mut row_start = 0_usize;
    for arrow in algebra.quiver().arrows() {
        let arrow_index = arrow.index().zero_based();
        let source_vertex = arrow.source().index().zero_based();
        let target_vertex = arrow.target().index().zero_based();
        let source_action = logical_submatrix(
            &source.arrow_matrices()[arrow_index],
            source_dimensions[source_vertex],
            source_dimensions[target_vertex],
        )?;
        let transposed_range_action = logical_submatrix(
            &range.arrow_matrices()[arrow_index],
            range_dimensions[source_vertex],
            range_dimensions[target_vertex],
        )?
        .transpose();
        let range_target_dimension = range_dimensions[target_vertex];
        let range_source_dimension = range_dimensions[source_vertex];

        for source_row in 0..source_dimensions[source_vertex] {
            let block_row = row_start + source_row * range_target_dimension;
            let block_column = column_starts[source_vertex] + source_row * range_source_dimension;
            for row in 0..range_target_dimension {
                for column in 0..range_source_dimension {
                    equations[(block_row + row, block_column + column)] =
                        transposed_range_action[(row, column)].clone();
                }
            }
        }

        let mut equation_row = row_start;
        for source_row in 0..source_dimensions[source_vertex] {
            for range_column in 0..range_target_dimension {
                let mut column = column_starts[target_vertex] + range_column;
                for source_target_row in 0..source_dimensions[target_vertex] {
                    let negative = field.neg(&source_action[(source_row, source_target_row)]);
                    field.add_assign(&mut equations[(equation_row, column)], &negative);
                    column += range_target_dimension;
                }
                equation_row += 1;
            }
        }
        row_start += source_dimensions[source_vertex] * range_target_dimension;
    }

    let nullspace = equations.transpose().nullspace_rows();
    let mut basis = Vec::with_capacity(nullspace.nrows());
    for basis_row in 0..nullspace.nrows() {
        let coefficients = nullspace.row(basis_row);
        let mut matrix = DenseMatrix::zero(source.dimension(), range.dimension(), field.clone());
        let mut coefficient = 0_usize;
        let mut row_offset = 0_usize;
        let mut column_offset = 0_usize;
        for vertex in 0..vertex_count {
            for row in 0..source_dimensions[vertex] {
                for column in 0..range_dimensions[vertex] {
                    matrix[(row_offset + row, column_offset + column)] =
                        coefficients[coefficient].clone();
                    coefficient += 1;
                }
            }
            row_offset += source_dimensions[vertex];
            column_offset += range_dimensions[vertex];
        }
        basis.push(matrix);
    }
    Ok(basis)
}
