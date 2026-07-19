/// Translate QPA's `RightApproximationByPerpT(T, M)` operation.
///
/// QPA requires the already-installed `IsCotiltingModule(T)` property, then
/// returns the zero map immediately for zero `M`. For nonzero `M`, this
/// source-faithful path takes the projective resolution through
/// `InjDimension(T)`, splices its short exact sequences by the same ordered
/// minimal left add(T)-approximations and pushouts, and finally applies
/// `RightMinimalVersion`.
///
/// # Errors
///
/// Preserves the unbound/false cotilting-property distinction, reports QPA's
/// unbound injective-dimension and zero-dimension indexing failures as typed
/// errors, and propagates resolution, approximation, isomorphism, pushout,
/// cokernel, composition, and right-minimal-version failures.
pub fn right_approximation_by_perp_t<F: QpaCoefficientField>(
    cotilting: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    require_cotilting_for_approximation(cotilting)?;
    if module.dimension() == 0 {
        return Ok(zero_module_mapping(
            zero_module(module.right_acting_algebra())?,
            module.clone(),
        )?);
    }
    let dimension = inj_dimension(cotilting)
        .ok_or(EnvelopingModuleError::MissingCotiltingApproximationInjectiveDimensionCache)?;
    let exact_sequences = cotilting_approximation_exact_sequences(module, dimension)?;
    let mut approximation =
        minimal_left_approximation(exact_sequences[dimension - 1].0.source(), cotilting)?;
    let mut result = None;
    for iteration in 1..=dimension {
        let exact_sequence = &exact_sequences[dimension - iteration];
        let comparison =
            cotilting_approximation_isomorphism(exact_sequence.0.source(), approximation.source())?;
        approximation = comparison.compose(&approximation)?;
        let pushout = cotilting_approximation_push_out(&exact_sequence.0, &approximation)?;
        let mut projection = cokernel_projection(&pushout)?;
        let comparison =
            cotilting_approximation_isomorphism(projection.range(), exact_sequence.1.range())?;
        projection = projection.compose(&comparison)?;
        if iteration < dimension {
            let next = minimal_left_approximation(projection.source(), cotilting)?;
            approximation = cotilting_approximation_push_out(&next, &projection)?;
        }
        result = Some(projection);
    }
    let result = result.ok_or(EnvelopingModuleError::EmptyCotiltingApproximationExactSequences)?;
    Ok(right_minimal_version(&result)?.0)
}
/// Translate QPA's `LeftApproximationByAddTHat(T, M)` operation.
///
/// This follows the same installed-cotilting check and projective-resolution
/// splice as [`right_approximation_by_perp_t`]. QPA omits the intermediate
/// source isomorphism here and performs one final minimal left
/// add(T)-approximation/pushout before taking `LeftMinimalVersion`.
///
/// # Errors
///
/// Preserves QPA's property and attribute ordering and propagates the exact
/// resolution, approximation, isomorphism, pushout, cokernel, composition,
/// and left-minimal-version failure paths.
pub fn left_approximation_by_add_t_hat<F: QpaCoefficientField>(
    cotilting: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    require_cotilting_for_approximation(cotilting)?;
    if module.dimension() == 0 {
        return Ok(zero_module_mapping(
            module.clone(),
            zero_module(module.right_acting_algebra())?,
        )?);
    }
    let dimension = inj_dimension(cotilting)
        .ok_or(EnvelopingModuleError::MissingCotiltingApproximationInjectiveDimensionCache)?;
    let exact_sequences = cotilting_approximation_exact_sequences(module, dimension)?;
    let mut approximation =
        minimal_left_approximation(exact_sequences[dimension - 1].0.source(), cotilting)?;
    for iteration in 1..=dimension {
        let exact_sequence = &exact_sequences[dimension - iteration];
        let pushout = cotilting_approximation_push_out(&exact_sequence.0, &approximation)?;
        let mut projection = cokernel_projection(&pushout)?;
        let comparison =
            cotilting_approximation_isomorphism(projection.range(), exact_sequence.1.range())?;
        projection = projection.compose(&comparison)?;
        let next = minimal_left_approximation(projection.source(), cotilting)?;
        approximation = cotilting_approximation_push_out(&next, &projection)?;
    }
    Ok(left_minimal_version(&approximation)?.0)
}

/// Translate QPA's `IsNthSyzygy(M, n)` operation.
///
/// The source takes the `n`th syzygy after dualizing, dualizes and takes the
/// `n`th syzygy again, adds the source of `ProjectiveCover(M)`, and asks QPA's
/// literal `IsDirectSummand(M, N)` decision. In particular, this preserves
/// QPA's perhaps surprising `false` result for the zero module.
///
/// # Errors
///
/// Propagates both dualizations, both `NthSyzygy` calls, the projective cover,
/// direct sum, and direct-summand decision. Negative degrees retain
/// [`nth_syzygy`]'s typed source failure.
pub fn is_nth_syzygy<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    degree: &Integer,
) -> Result<bool, EnvelopingModuleError> {
    let dual = dual_of_module(module)?;
    let dual_syzygy = nth_syzygy(&dual, degree)?;
    let double_dual_syzygy = dual_of_module(&dual_syzygy)?;
    let candidate = nth_syzygy(&double_dual_syzygy, degree)?;
    let projective = projective_cover(module)?;
    let ambient = direct_sum_of_qpa_modules(&[candidate, projective.source().clone()])?
        .ok_or(EnvelopingModuleError::InternalNthSyzygyDirectSumFailed)?;
    Ok(is_direct_summand(module, &ambient)?)
}

/// Translate QPA's list overload `RightApproximationByAddM(L, C)`.
///
/// QPA scans `L` in order. For each `L_i`, it starts with maps to `C` that
/// factor through already selected approximation maps and maps obtained by
/// precomposing `Hom(L_i,C)` with `rad End(L_i)`. It then scans QPA's ordered
/// Hom basis and retains precisely those flattened map rows that extend this
/// span. The selected sources are direct-summed in selection order and their
/// projection composites are added. This order-sensitive algorithm is kept
/// separate from the one-module quotient-basis implementation in
/// [`right_approximation_by_add_m`].
///
/// # Errors
///
/// An empty list or zero target takes QPA's zero-map shortcut before checking
/// acting algebras. Otherwise this rejects a foreign list entry, preserves the
/// source's failure when no map is selected from a nonempty list, and
/// propagates Hom, radical, exact row-span, direct-sum, and map failures.
pub fn right_approximation_by_add_m_list<F: QpaCoefficientField>(
    modules: &[PathAlgebraModule<F>],
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let algebra = target.right_acting_algebra();
    if modules.is_empty() || target.dimension() == 0 {
        return Ok(zero_module_mapping(zero_module(algebra)?, target.clone())?);
    }
    if modules
        .iter()
        .any(|module| module.right_acting_algebra() != algebra)
    {
        return Err(EnvelopingModuleError::ListApproximationModulesOverDifferentAlgebras);
    }

    let field = algebra.field().clone();
    let mut approximation = Vec::<PathAlgebraModuleHomomorphism<F>>::new();
    for module in modules {
        let homomorphisms = hom_over_algebra(module, target)?
            .ok_or(EnvelopingModuleError::ApproximationHomUnavailable)?;
        let Some(first_homomorphism) = homomorphisms.first() else {
            continue;
        };
        let row_length = first_homomorphism.flat_vertex_matrix_entries().len();
        let mut generator_rows = Vec::<Vec<F::Element>>::new();

        for selected in &approximation {
            let factoring_maps = hom_over_algebra(module, selected.source())?
                .ok_or(EnvelopingModuleError::ApproximationHomUnavailable)?;
            for factoring in factoring_maps {
                let composite = factoring.compose(selected)?;
                if !composite.is_zero() {
                    generator_rows.push(composite.flat_vertex_matrix_entries());
                }
            }
        }

        let radical_maps = end_over_algebra(module)?
            .radical_basis_matrices()?
            .iter()
            .map(|matrix| from_end_m_to_hom_mm(module, matrix))
            .collect::<Result<Vec<_>, _>>()?;
        for homomorphism in &homomorphisms {
            for radical in &radical_maps {
                generator_rows.push(radical.compose(homomorphism)?.flat_vertex_matrix_entries());
            }
        }

        let homomorphism_rows = homomorphisms
            .iter()
            .map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
            .collect::<Vec<_>>();
        if generator_rows.is_empty() {
            generator_rows.push(
                zero_module_mapping(module.clone(), target.clone())?.flat_vertex_matrix_entries(),
            );
        }
        for candidate in &homomorphism_rows {
            let span = DenseMatrix::from_flat(
                generator_rows.len(),
                row_length,
                generator_rows.iter().flatten().cloned().collect(),
                field.clone(),
            )?;
            if span.solution_row(candidate)?.is_none() {
                generator_rows.push(candidate.clone());
                let position = homomorphism_rows
                    .iter()
                    .position(|row| row == candidate)
                    .ok_or(EnvelopingModuleError::MissingListApproximationHomPosition)?;
                approximation.push(homomorphisms[position].clone());
            }
        }
    }

    let sources = approximation
        .iter()
        .map(|homomorphism| homomorphism.source().clone())
        .collect::<Vec<_>>();
    let direct_sum = direct_sum_of_qpa_modules(&sources)?
        .ok_or(EnvelopingModuleError::EmptyListApproximationDirectSum)?;
    let projections = direct_sum_projections(&direct_sum)?;
    if projections.len() != approximation.len() {
        return Err(EnvelopingModuleError::ApproximationMapCountMismatch);
    }
    sum_nonempty_module_maps(
        projections
            .iter()
            .zip(&approximation)
            .map(|(projection, homomorphism)| Ok(projection.compose(homomorphism)?))
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?,
    )
}

/// Translate QPA's `RadicalRightApproximationByAddM(modulelist, t)`.
///
/// The one-based entry `N = modulelist[t]` is first approximated by all other
/// list entries using [`right_approximation_by_add_m_list`]. QPA then maps
/// `Hom(N, Source(f))` through `f`, quotients `rad End(N)` by the span of those
/// composites, and appends one copy of `N` for each radical basis map with a
/// nonzero quotient image. Projection composites are summed in the resulting
/// source order.
///
/// # Errors
///
/// Rejects a GAP integer outside the one-based module-list interval and
/// propagates list approximation, Hom, radical, row-span, direct-sum, and map
/// failures.
pub fn radical_right_approximation_by_add_m<F: QpaCoefficientField>(
    modules: &[PathAlgebraModule<F>],
    index: &Integer,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let length = modules.len();
    if index < &Integer::one() || index > &Integer::from(length) {
        return Err(EnvelopingModuleError::RadicalApproximationIndexOutOfRange {
            index: index.clone(),
            length,
        });
    }
    let position = usize::try_from(index.clone()).map_err(|_| {
        EnvelopingModuleError::RadicalApproximationIndexOutOfRange {
            index: index.clone(),
            length,
        }
    })? - 1;
    let module = &modules[position];
    let other_modules = modules
        .iter()
        .enumerate()
        .filter(|(candidate, _)| *candidate != position)
        .map(|(_, module)| module.clone())
        .collect::<Vec<_>>();
    let mut approximation = right_approximation_by_add_m_list(&other_modules, module)?;
    let approximation_source = approximation.source().clone();

    let radical_maps = end_over_algebra(module)?
        .radical_basis_matrices()?
        .iter()
        .map(|matrix| from_end_m_to_hom_mm(module, matrix))
        .collect::<Result<Vec<_>, _>>()?;
    if radical_maps.is_empty() {
        return Ok(approximation);
    }

    let factoring_composites = hom_over_algebra(module, &approximation_source)?
        .ok_or(EnvelopingModuleError::ApproximationHomUnavailable)?
        .iter()
        .map(|homomorphism| homomorphism.compose(&approximation))
        .collect::<Result<Vec<_>, _>>()?;
    let row_length = radical_maps[0].flat_vertex_matrix_entries().len();
    let factoring_span = if factoring_composites.is_empty() {
        None
    } else {
        Some(DenseMatrix::from_flat(
            factoring_composites.len(),
            row_length,
            factoring_composites
                .iter()
                .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
                .collect(),
            module.right_acting_algebra().field().clone(),
        )?)
    };
    let selected_radicals = radical_maps
        .iter()
        .filter_map(|radical| {
            let row = radical.flat_vertex_matrix_entries();
            match flat_row_in_span(
                factoring_span.as_ref(),
                &row,
                module.right_acting_algebra().field(),
            ) {
                Ok(true) => None,
                Ok(false) => Some(Ok(radical.clone())),
                Err(error) => Some(Err(error)),
            }
        })
        .collect::<Result<Vec<_>, MatrixError>>()?;

    let mut sources = Vec::with_capacity(selected_radicals.len() + 1);
    sources.push(approximation_source);
    sources.extend((0..selected_radicals.len()).map(|_| module.clone()));
    let direct_sum = direct_sum_of_qpa_modules(&sources)?
        .ok_or(EnvelopingModuleError::InternalRadicalApproximationDirectSumFailed)?;
    let projections = direct_sum_projections(&direct_sum)?;
    let mut homomorphisms = Vec::with_capacity(selected_radicals.len() + 1);
    homomorphisms.push(approximation);
    homomorphisms.extend(selected_radicals);
    if projections.len() != homomorphisms.len() {
        return Err(EnvelopingModuleError::ApproximationMapCountMismatch);
    }
    approximation = sum_nonempty_module_maps(
        projections
            .iter()
            .zip(&homomorphisms)
            .map(|(projection, homomorphism)| Ok(projection.compose(homomorphism)?))
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?,
    )?;
    Ok(approximation)
}

/// QPA's heterogeneous first result component from
/// `ProjectiveResolutionOfSimpleModuleOverEndo`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ProjectiveResolutionOfSimpleModuleOverEndoStatus {
    /// QPA found this finite projective dimension.
    Finite(Integer),
    /// QPA exhausted the entered loop and returned `"projdim > t"`.
    ///
    /// The source interpolates the one-based module index `t`, rather than
    /// the entered resolution length. The index is retained here so bindings
    /// can reproduce that literal result without making the Rust API stringly
    /// typed.
    GreaterThanModuleIndex(Integer),
}

/// QPA's pair returned by `ProjectiveResolutionOfSimpleModuleOverEndo`.
#[derive(Clone, Debug)]
pub struct ProjectiveResolutionOfSimpleModuleOverEndoResult<F: Field> {
    projective_dimension: ProjectiveResolutionOfSimpleModuleOverEndoStatus,
    syzygies: Vec<PathAlgebraModule<F>>,
}

impl<F: Field> ProjectiveResolutionOfSimpleModuleOverEndoResult<F> {
    /// Return QPA's integer-or-string first result component in typed form.
    #[must_use]
    pub const fn projective_dimension(&self) -> &ProjectiveResolutionOfSimpleModuleOverEndoStatus {
        &self.projective_dimension
    }

    /// Return the nonzero stripped kernels accumulated in QPA loop order.
    #[must_use]
    pub fn syzygies(&self) -> &[PathAlgebraModule<F>] {
        &self.syzygies
    }
}

/// Translate QPA's `ProjectiveResolutionOfSimpleModuleOverEndo(L, t, length)`.
///
/// The source begins with [`radical_right_approximation_by_add_m`] at the
/// one-based list position `t`. In each degree it checks a zero source before
/// injectivity, takes the kernel, and repeatedly removes every common direct
/// summand with each entry of `L` in list order. A nonzero stripped kernel is
/// recorded and approximated through the list overload
/// [`right_approximation_by_add_m_list`]. The GAP range `[1..length-2]` and
/// all of its early-return offsets are preserved with arbitrary-precision
/// integers.
///
/// # Errors
///
/// Rejects a negative length before inspecting the list or index, then
/// propagates radical/list approximation, kernel, and common-summand
/// failures. QPA's invalid attempt to index the `fail` returned by
/// `CommonDirectSummand` for different algebras is exposed as a typed error.
pub fn projective_resolution_of_simple_module_over_endo<F: QpaCoefficientField>(
    modules: &[PathAlgebraModule<F>],
    index: &Integer,
    length: &Integer,
) -> Result<ProjectiveResolutionOfSimpleModuleOverEndoResult<F>, EnvelopingModuleError> {
    if length < &Integer::from(0) {
        return Err(
            EnvelopingModuleError::NegativeSimpleEndoProjectiveResolutionLength {
                length: length.clone(),
            },
        );
    }

    let mut approximation = radical_right_approximation_by_add_m(modules, index)?;
    let mut syzygies = Vec::new();
    let cutoff = length - 2;
    let mut degree = Integer::one();
    while degree <= cutoff {
        if approximation.source().dimension() == 0 {
            return Ok(ProjectiveResolutionOfSimpleModuleOverEndoResult {
                projective_dimension: ProjectiveResolutionOfSimpleModuleOverEndoStatus::Finite(
                    &degree - 1,
                ),
                syzygies,
            });
        }
        if approximation.is_injective() {
            return Ok(ProjectiveResolutionOfSimpleModuleOverEndoResult {
                projective_dimension: ProjectiveResolutionOfSimpleModuleOverEndoStatus::Finite(
                    degree,
                ),
                syzygies,
            });
        }

        let mut kernel = kernel_of_module_homomorphism(&approximation)?;
        for module in modules {
            loop {
                match common_direct_summand(&kernel, module)? {
                    CommonDirectSummandResult::DifferentAlgebras => {
                        return Err(EnvelopingModuleError::SimpleEndoResolutionCommonSummandFail);
                    }
                    CommonDirectSummandResult::NoCommonSummand => break,
                    CommonDirectSummandResult::Found(common) => {
                        kernel = common.source_complement().clone();
                    }
                }
            }
        }
        if kernel.dimension() == 0 {
            return Ok(ProjectiveResolutionOfSimpleModuleOverEndoResult {
                projective_dimension: ProjectiveResolutionOfSimpleModuleOverEndoStatus::Finite(
                    &degree + 1,
                ),
                syzygies,
            });
        }
        syzygies.push(kernel.clone());
        approximation = right_approximation_by_add_m_list(modules, &kernel)?;
        degree += 1;
    }

    Ok(ProjectiveResolutionOfSimpleModuleOverEndoResult {
        projective_dimension:
            ProjectiveResolutionOfSimpleModuleOverEndoStatus::GreaterThanModuleIndex(index.clone()),
        syzygies,
    })
}

#[derive(Clone, Copy)]
enum ApproximationSide {
    Right,
    Left,
}

fn require_common_approximation_algebra<F: Field>(
    module: &PathAlgebraModule<F>,
    other: &PathAlgebraModule<F>,
) -> Result<(), EnvelopingModuleError> {
    if module.right_acting_algebra() != other.right_acting_algebra() {
        return Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras);
    }
    Ok(())
}

fn add_m_approximation_generators<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    other: &PathAlgebraModule<F>,
    side: ApproximationSide,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, EnvelopingModuleError> {
    let homomorphisms = match side {
        ApproximationSide::Right => hom_over_algebra(module, other)?,
        ApproximationSide::Left => hom_over_algebra(other, module)?,
    }
    .ok_or(EnvelopingModuleError::ApproximationHomUnavailable)?;
    let Some(first) = homomorphisms.first() else {
        return Ok(Vec::new());
    };
    let field = module.right_acting_algebra().field().clone();
    let physical_length = first.flat_vertex_matrix_entries().len();
    let hom_basis = DenseMatrix::from_flat(
        homomorphisms.len(),
        physical_length,
        homomorphisms
            .iter()
            .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
            .collect(),
        field.clone(),
    )?;
    let radical_maps = end_over_algebra(module)?
        .radical_basis_matrices()?
        .iter()
        .map(|matrix| from_end_m_to_hom_mm(module, matrix))
        .collect::<Result<Vec<_>, _>>()?;
    let mut radical_composites = Vec::with_capacity(
        homomorphisms
            .len()
            .checked_mul(radical_maps.len())
            .ok_or(EnvelopingModuleError::ApproximationMapCountOverflow)?,
    );
    for homomorphism in &homomorphisms {
        for radical in &radical_maps {
            radical_composites.push(match side {
                ApproximationSide::Right => radical.compose(homomorphism)?,
                ApproximationSide::Left => homomorphism.compose(radical)?,
            });
        }
    }
    let radical_basis = if radical_composites.is_empty() {
        DenseMatrix::from_flat(0, physical_length, Vec::new(), field.clone())?
    } else {
        DenseMatrix::from_flat(
            radical_composites.len(),
            physical_length,
            radical_composites
                .iter()
                .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
                .collect(),
            field.clone(),
        )?
        .semi_echelon_transformation()
        .vectors
    };
    let quotient = quotient_complement_in_space(&hom_basis, &radical_basis)?;
    let (source, range) = match side {
        ApproximationSide::Right => (module, other),
        ApproximationSide::Left => (other, module),
    };
    quotient
        .representatives
        .iter()
        .map(|row| module_map_from_physical_flat(source, range, row, &field))
        .collect()
}

fn sum_nonempty_module_maps<F: Field>(
    maps: Vec<PathAlgebraModuleHomomorphism<F>>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let mut maps = maps.into_iter();
    let mut sum = maps
        .next()
        .ok_or(EnvelopingModuleError::EmptyApproximationGeneratorList)?;
    for homomorphism in maps {
        sum = sum.add(&homomorphism)?;
    }
    Ok(sum)
}
