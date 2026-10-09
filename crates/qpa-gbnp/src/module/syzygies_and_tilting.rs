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
        // HomFromProjective requests the complete family for each generator.
        // Its algebra-side cache is weak: retaining only the selected sources
        // lets unused members expire and forces the next call to rebuild the
        // family. Keep all members live until every generator map is made.
        // Initialize on the first map, preserving HomFromProjective's checks
        // and error order for unsupported algebras.
        let mut projectives = None;
        let maps = generators
            .iter()
            .map(|generator| {
                hom_from_projective_with_family(generator, module, &mut projectives)?
                    .ok_or(EnvelopingModuleError::MinimalGeneratorNotVertexUniform)
            })
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
        drop(projectives);
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

/// Translate QPA's `NthSyzygy(M, n)` operation.
///
/// Negative degrees are rejected, degree zero returns the identical entered
/// module, and a positive degree returns the kernel of differential `n - 1`
/// in QPA's lazy [`crate::ProjectiveResolution`].
///
/// # Errors
///
/// Rejects a negative GAP integer and propagates projective-resolution and
/// kernel construction failures.
pub fn nth_syzygy<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    degree: &Integer,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if degree < &Integer::from(0) {
        return Err(EnvelopingModuleError::NegativeSyzygyDegree {
            degree: degree.clone(),
        });
    }
    if degree == &Integer::from(0) {
        return Ok(module.clone());
    }
    let resolution = projective_resolution(module)?;
    let differential = resolution.differential(&(degree - &Integer::one()))?;
    Ok(kernel_of_module_homomorphism(&differential)?)
}

/// Translate QPA's `IsOmegaPeriodic(M, n)` search.
///
/// Starting with the entered module, the source computes successive first
/// syzygies and returns the first one-based degree at which the result is
/// isomorphic to `M`. `None` is QPA's `false`, including for zero and negative
/// bounds. The unbounded GAP integer is retained rather than narrowed to a
/// machine integer.
///
/// QPA prints `Computing syzygy number: i` before every iteration. The normal
/// Rust entry point is quiet; [`is_omega_periodic_with_progress`] exposes the
/// same one-based iteration values without imposing process-global output.
///
/// # Errors
///
/// Propagates first-syzygy construction and QPA's module-isomorphism test.
pub fn is_omega_periodic<F: Field>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<Option<Integer>, EnvelopingModuleError> {
    is_omega_periodic_with_progress(module, bound, |_| {})
}

/// Run [`is_omega_periodic`] while reporting each iteration at the exact point
/// where QPA emits its progress line.
///
/// # Errors
///
/// Propagates first-syzygy construction and QPA's module-isomorphism test.
pub fn is_omega_periodic_with_progress<F: Field>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
    mut progress: impl FnMut(&Integer),
) -> Result<Option<Integer>, EnvelopingModuleError> {
    let mut current = module.clone();
    let mut degree = Integer::one();
    while &degree <= bound {
        progress(&degree);
        let next = first_syzygy(&current)?;
        if isomorphic_modules(module, &next)? {
            return Ok(Some(degree));
        }
        current = next;
        degree += 1;
    }
    Ok(None)
}

/// Translate QPA's `IsTauPeriodic(M, n)` search.
///
/// Starting with the entered module, the source repeatedly applies [`d_tr`]
/// and returns the first one-based degree at which the result is isomorphic to
/// `M`. `None` is QPA's `false`, including for zero and negative bounds. The
/// unbounded GAP integer is retained rather than narrowed to a machine integer.
///
/// # Errors
///
/// Propagates `DTr` construction and QPA's module-isomorphism test.
pub fn is_tau_periodic<F: Field>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<Option<Integer>, EnvelopingModuleError> {
    let mut current = module.clone();
    let mut degree = Integer::one();
    while &degree <= bound {
        let next = d_tr(&current)?;
        if isomorphic_modules(module, &next)? {
            return Ok(Some(degree));
        }
        current = next;
        degree += 1;
    }
    Ok(None)
}

/// Translate QPA's `RightApproximationByAddM(M, C)` operation.
///
/// The source computes the ordered basis of `Hom(M,C)`, quotients it by
/// `rad(End(M)) Hom(M,C)`, chooses GAP's natural quotient-basis lifts, and
/// sums the corresponding maps from a direct sum of copies of `M`. This map
/// is not necessarily minimal.
///
/// # Errors
///
/// Requires identical acting algebras and propagates Hom, endomorphism
/// radical, exact quotient-space, direct-sum, and map-arithmetic failures.
pub fn right_approximation_by_add_m<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    require_common_approximation_algebra(module, target)?;
    if target.dimension() == 0 {
        return Ok(zero_module_mapping(
            zero_module(module.right_acting_algebra())?,
            target.clone(),
        )?);
    }
    let generators = add_m_approximation_generators(module, target, ApproximationSide::Right)?;
    if generators.is_empty() {
        return Ok(zero_module_mapping(
            zero_module(module.right_acting_algebra())?,
            target.clone(),
        )?);
    }
    let copies = vec![module.clone(); generators.len()];
    let approximation = direct_sum_of_qpa_modules(&copies)?
        .ok_or(EnvelopingModuleError::EmptyApproximationDirectSum)?;
    let projections = direct_sum_projections(&approximation)?;
    if projections.len() != generators.len() {
        return Err(EnvelopingModuleError::ApproximationMapCountMismatch);
    }
    sum_nonempty_module_maps(
        projections
            .iter()
            .zip(&generators)
            .map(|(projection, generator)| Ok(projection.compose(generator)?))
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?,
    )
}

/// Translate QPA's `MinimalRightApproximation(M, C)` operation.
///
/// This is the exact source composition of [`right_approximation_by_add_m`]
/// with QPA's cached `RightMinimalVersion`, returning its first component.
///
/// # Errors
///
/// Propagates approximation and right-minimal-version failures.
pub fn minimal_right_approximation<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let approximation = right_approximation_by_add_m(module, target)?;
    Ok(right_minimal_version(&approximation)?.0)
}

/// Translate QPA's `LeftApproximationByAddM(C, M)` operation.
///
/// The source computes the ordered basis of `Hom(C,M)`, quotients it by
/// `Hom(C,M) rad(End(M))`, chooses GAP's natural quotient-basis lifts, and
/// sums the corresponding maps into a direct sum of copies of `M`. This map
/// is not necessarily minimal.
///
/// # Errors
///
/// Requires identical acting algebras and propagates Hom, endomorphism
/// radical, exact quotient-space, direct-sum, and map-arithmetic failures.
pub fn left_approximation_by_add_m<F: QpaCoefficientField>(
    source: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    require_common_approximation_algebra(module, source)?;
    if source.dimension() == 0 {
        return Ok(zero_module_mapping(
            source.clone(),
            zero_module(module.right_acting_algebra())?,
        )?);
    }
    let generators = add_m_approximation_generators(module, source, ApproximationSide::Left)?;
    if generators.is_empty() {
        return Ok(zero_module_mapping(
            source.clone(),
            zero_module(module.right_acting_algebra())?,
        )?);
    }
    let copies = vec![module.clone(); generators.len()];
    let approximation = direct_sum_of_qpa_modules(&copies)?
        .ok_or(EnvelopingModuleError::EmptyApproximationDirectSum)?;
    let inclusions = direct_sum_inclusions(&approximation)?;
    if inclusions.len() != generators.len() {
        return Err(EnvelopingModuleError::ApproximationMapCountMismatch);
    }
    sum_nonempty_module_maps(
        generators
            .iter()
            .zip(&inclusions)
            .map(|(generator, inclusion)| Ok(generator.compose(inclusion)?))
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?,
    )
}

/// Translate QPA's `MinimalLeftApproximation(C, M)` operation.
///
/// This is the exact source composition of [`left_approximation_by_add_m`]
/// with QPA's cached `LeftMinimalVersion`, returning its first component.
///
/// # Errors
///
/// Propagates approximation and left-minimal-version failures.
pub fn minimal_left_approximation<F: QpaCoefficientField>(
    source: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let approximation = left_approximation_by_add_m(source, module)?;
    Ok(left_minimal_version(&approximation)?.0)
}

/// Translate QPA's `LeftFacMApproximation(C, M)` operation.
///
/// The source takes the projective cover `p: P(C) -> C`, computes the minimal
/// left add(M)-approximation `f: P(C) -> M'`, and returns `PushOut(p, f)[2]`,
/// the induced map from `C` into the pushout.
///
/// # Errors
///
/// Rejects modules over different algebras and propagates projective-cover,
/// add(M)-approximation, and pushout failures.
pub fn left_fac_m_approximation<F: QpaCoefficientField>(
    source: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    require_common_approximation_algebra(source, module)?;
    let cover = projective_cover(source)?;
    let approximation = minimal_left_approximation(cover.source(), module)?;
    let pushout = push_out(&cover, &approximation)?
        .ok_or(EnvelopingModuleError::InternalLeftFacPushOutFailed)?;
    Ok(pushout.1)
}

/// Translate QPA's `MinimalLeftFacMApproximation(C, M)` operation.
///
/// This is the exact source composition of [`left_fac_m_approximation`] with
/// QPA's cached `LeftMinimalVersion`, returning its first component.
///
/// # Errors
///
/// Propagates left-Fac(M)-approximation and left-minimal-version failures.
pub fn minimal_left_fac_m_approximation<F: QpaCoefficientField>(
    source: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let approximation = left_fac_m_approximation(source, module)?;
    Ok(left_minimal_version(&approximation)?.0)
}

/// Translate QPA's `RightSubMApproximation(M, C)` operation.
///
/// QPA constructs the injective-envelope map of `C` as the dual of the
/// projective cover of `DualOfModule(C)`, takes the minimal right
/// add(M)-approximation of its range, and returns `PullBack(i, f)[2]`, the map
/// from the pullback into `C`.
///
/// # Errors
///
/// Rejects modules over different algebras and propagates dualization,
/// projective-cover, add(M)-approximation, and pullback failures.
pub fn right_sub_m_approximation<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    require_common_approximation_algebra(module, target)?;
    let dual_target = dual_of_module(target)?;
    let opposite_cover = projective_cover(&dual_target)?;
    let envelope = dual_of_module_homomorphism(&opposite_cover)?;
    let approximation = minimal_right_approximation(module, envelope.range())?;
    let pullback = pull_back(&envelope, &approximation)?
        .ok_or(EnvelopingModuleError::InternalRightSubPullBackFailed)?;
    Ok(pullback.1)
}

/// Translate QPA's `MinimalRightSubMApproximation(M, C)` operation.
///
/// This is the exact source composition of [`right_sub_m_approximation`] with
/// QPA's cached `RightMinimalVersion`, returning its first component.
///
/// # Errors
///
/// Propagates right-Sub(M)-approximation and right-minimal-version failures.
pub fn minimal_right_sub_m_approximation<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let approximation = right_sub_m_approximation(module, target)?;
    Ok(right_minimal_version(&approximation)?.0)
}

/// Translate QPA's `LeftSubMApproximation(C, M)` operation.
///
/// The source computes the minimal left add(M)-approximation of `C` and
/// returns its natural projection onto its image. QPA documents this map as
/// already minimal.
///
/// # Errors
///
/// Rejects modules over different algebras and propagates add(M)-approximation
/// and image-projection failures.
pub fn left_sub_m_approximation<F: QpaCoefficientField>(
    source: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    require_common_approximation_algebra(source, module)?;
    let approximation = minimal_left_approximation(source, module)?;
    Ok(image_projection(&approximation)?)
}

/// QPA synonym of [`left_sub_m_approximation`].
///
/// # Errors
///
/// Propagates [`left_sub_m_approximation`].
pub fn minimal_left_sub_m_approximation<F: QpaCoefficientField>(
    source: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    left_sub_m_approximation(source, module)
}

/// The finite complex returned by QPA's
/// `HaveFiniteCoresolutionInAddM(N, M, n)` operation.
///
/// `base_position` and `differentials` are exactly the arguments passed to
/// QPA's `FiniteComplex`: the first stored map is the differential in the
/// base degree and maps to the object one degree below. Outside the stored
/// interval, [`Self::differential`] reproduces `FiniteComplex`'s zero tails.
#[derive(Clone, Debug)]
pub struct AddMCoresolution<F: Field> {
    base_position: Integer,
    differentials: Vec<PathAlgebraModuleHomomorphism<F>>,
    zero: PathAlgebraModule<F>,
}

impl<F: Field> AddMCoresolution<F> {
    /// Return the base degree passed to QPA's `FiniteComplex` constructor.
    #[must_use]
    pub const fn base_position(&self) -> &Integer {
        &self.base_position
    }

    /// Return QPA's reversed list passed to `FiniteComplex`.
    #[must_use]
    pub fn differentials(&self) -> &[PathAlgebraModuleHomomorphism<F>] {
        &self.differentials
    }

    /// Return QPA's `DifferentialOfComplex` map in an arbitrary degree.
    ///
    /// # Errors
    ///
    /// Propagates zero-map construction at either finite boundary.
    pub fn differential(
        &self,
        degree: &Integer,
    ) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
        let offset = degree - &self.base_position;
        if offset >= 0
            && let Ok(position) = usize::try_from(offset)
            && let Some(differential) = self.differentials.get(position)
        {
            return Ok(differential.clone());
        }
        if degree == &(&self.base_position - Integer::one()) {
            return zero_module_mapping(self.differentials[0].range().clone(), self.zero.clone());
        }
        if degree == &(&self.base_position + Integer::from(self.differentials.len())) {
            let Some(last) = self.differentials.last() else {
                return zero_module_mapping(self.zero.clone(), self.zero.clone());
            };
            return zero_module_mapping(self.zero.clone(), last.source().clone());
        }
        zero_module_mapping(self.zero.clone(), self.zero.clone())
    }

    /// Return QPA's `ObjectOfComplex` value in an arbitrary degree.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::differential`].
    pub fn object(&self, degree: &Integer) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
        Ok(self.differential(degree)?.source().clone())
    }
}

/// Translate QPA's `HaveFiniteCoresolutionInAddM(N, M, n)` operation.
///
/// `None` is QPA's `false`. The zero bound uses QPA's separate isomorphism
/// test. A negative bound makes GAP's literal range `[0..n]` empty, after
/// which `FiniteComplex` fails by indexing the empty differential list; Rust
/// reports that source failure after computing the initial approximation. For
/// a positive bound, every iteration requires the current minimal left
/// add(M)-approximation to be injective, appends `g*f`, and continues from its
/// cokernel. QPA also computes the next approximation after a nonzero final
/// cokernel before leaving the range. The maps are then reversed and passed
/// to `FiniteComplex` in degree `-Length + 1`.
///
/// # Errors
///
/// Rejects modules over different algebras and propagates approximation,
/// composition, cokernel, identity-map, and zero-module failures.
pub fn have_finite_coresolution_in_add_m<F: QpaCoefficientField>(
    source: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<Option<AddMCoresolution<F>>, EnvelopingModuleError> {
    require_common_approximation_algebra(source, module)?;
    let mut differentials = Vec::new();
    let mut connecting = identity_module_mapping(source.clone())?;
    let mut approximation = minimal_left_approximation(source, module)?;

    if bound == &Integer::from(0) {
        if !approximation.is_isomorphism() {
            return Ok(None);
        }
        return Ok(Some(AddMCoresolution {
            base_position: Integer::from(0),
            differentials: vec![approximation],
            zero: zero_module(module.right_acting_algebra())?,
        }));
    }

    if bound < &Integer::from(0) {
        return Err(EnvelopingModuleError::EmptyCoresolutionDifferentialList {
            bound: bound.clone(),
        });
    }

    let mut degree = Integer::from(0);
    loop {
        if !approximation.is_injective() {
            return Ok(None);
        }
        differentials.push(connecting.compose(&approximation)?);
        connecting = cokernel_projection(&approximation)?;
        if connecting.range().dimension() == 0 {
            break;
        }
        approximation = minimal_left_approximation(connecting.range(), module)?;
        if &degree == bound {
            break;
        }
        degree += 1;
    }

    differentials.reverse();
    let base_position = Integer::one() - Integer::from(differentials.len());
    Ok(Some(AddMCoresolution {
        base_position,
        differentials,
        zero: zero_module(module.right_acting_algebra())?,
    }))
}

/// Return QPA's installed `IsTiltingModule(M)` attribute, if bound.
///
/// QPA's bounded operation installs only the successful value `true`; Rust
/// exposes the unbound attribute state as `None`.
#[must_use]
pub fn is_tilting_module<F: Field>(module: &PathAlgebraModule<F>) -> Option<bool> {
    module.cached_is_tilting_module()
}

/// QPA's successful `[m, coresolution]` value from `TiltingModule(M, n)`.
#[derive(Clone, Debug)]
pub struct TiltingModuleResult<F: Field> {
    projective_dimension: usize,
    coresolutions: Vec<AddMCoresolution<F>>,
}

impl<F: Field> TiltingModuleResult<F> {
    /// Return QPA's first result component, the projective dimension of `M`.
    #[must_use]
    pub const fn projective_dimension(&self) -> usize {
        self.projective_dimension
    }

    /// Return the projective-indexed add(M)-coresolutions in QPA order.
    #[must_use]
    pub fn coresolutions(&self) -> &[AddMCoresolution<F>] {
        &self.coresolutions
    }
}

/// Translate QPA's bounded `TiltingModule(M, n)` operation.
///
/// `None` is QPA's `false`. The source first establishes projective dimension
/// at most `n`, then performs its literal self-orthogonality loop. In
/// particular, for projective dimension `m`, QPA makes `m + 1` Ext calls but
/// advances the syzygy only while the one-based loop counter is strictly less
/// than `m`, so its final Ext degree is repeated. Finally, every
/// indecomposable projective must have the source's finite add(M)-coresolution
/// of bound `m`. Only success installs `IsTiltingModule(M) := true`.
///
/// # Errors
///
/// Propagates bounded projective-dimension, Ext, syzygy, indecomposable
/// projective, and finite-coresolution failures.
pub fn tilting_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<Option<TiltingModuleResult<F>>, EnvelopingModuleError> {
    match proj_dimension(module) {
        Some(dimension) => {
            if &Integer::from(dimension) > bound {
                return Ok(None);
            }
        }
        None => {
            if proj_dimension_of_module(module, bound)?.is_none() {
                return Ok(None);
            }
        }
    }
    let dimension = proj_dimension(module)
        .ok_or(EnvelopingModuleError::MissingTiltingProjectiveDimensionCache)?;

    let mut syzygy = module.clone();
    let syzygy_cutoff = dimension.saturating_sub(1);
    for iteration in 0..=dimension {
        if ext_dimension_over_algebra(&syzygy, module)? != 0 {
            return Ok(None);
        }
        if iteration < syzygy_cutoff {
            syzygy = nth_syzygy(&syzygy, &Integer::one())?;
        }
    }

    let coresolution_bound = Integer::from(dimension);
    let mut coresolutions = Vec::new();
    for projective in indec_projective_modules(module.right_acting_algebra())? {
        let Some(coresolution) =
            have_finite_coresolution_in_add_m(&projective, module, &coresolution_bound)?
        else {
            return Ok(None);
        };
        coresolutions.push(coresolution);
    }
    let _ = module.cache_is_tilting_module(true);
    Ok(Some(TiltingModuleResult {
        projective_dimension: dimension,
        coresolutions,
    }))
}

/// The finite complex returned by a successful positive-bound
/// `HaveFiniteResolutionInAddM(N, M, n)` call.
///
/// The first stored map is QPA's differential in degree one. The remaining
/// maps retain the source order, and [`Self::differential`] reproduces
/// `FiniteComplex`'s zero tails.
#[derive(Clone, Debug)]
pub struct AddMResolution<F: Field> {
    differentials: Vec<PathAlgebraModuleHomomorphism<F>>,
    zero: PathAlgebraModule<F>,
}

impl<F: Field> AddMResolution<F> {
    /// Return the source list passed to `FiniteComplex(cat, 1, ...)`.
    #[must_use]
    pub fn differentials(&self) -> &[PathAlgebraModuleHomomorphism<F>] {
        &self.differentials
    }

    /// Return QPA's `DifferentialOfComplex` map in an arbitrary degree.
    ///
    /// # Errors
    ///
    /// Propagates zero-map construction at either finite boundary.
    pub fn differential(
        &self,
        degree: &Integer,
    ) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
        let offset = degree - Integer::one();
        if offset >= 0
            && let Ok(position) = usize::try_from(offset)
            && let Some(differential) = self.differentials.get(position)
        {
            return Ok(differential.clone());
        }
        if degree == &Integer::from(0) {
            return zero_module_mapping(self.differentials[0].range().clone(), self.zero.clone());
        }
        if degree == &Integer::from(self.differentials.len() + 1) {
            let Some(last) = self.differentials.last() else {
                return zero_module_mapping(self.zero.clone(), self.zero.clone());
            };
            return zero_module_mapping(self.zero.clone(), last.source().clone());
        }
        zero_module_mapping(self.zero.clone(), self.zero.clone())
    }

    /// Return QPA's `ObjectOfComplex` value in an arbitrary degree.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::differential`].
    pub fn object(&self, degree: &Integer) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
        Ok(self.differential(degree)?.source().clone())
    }
}

/// QPA's three-valued successful result for
/// `HaveFiniteResolutionInAddM(N, M, n)`.
#[derive(Clone, Debug)]
pub enum FiniteAddMResolutionResult<F: Field> {
    /// QPA's literal `true` result for the zero-bound isomorphism branch.
    InAddM,
    /// QPA's finite complex result for a nonzero bound.
    Resolution(AddMResolution<F>),
}

/// Translate QPA's `HaveFiniteResolutionInAddM(N, M, n)` operation.
///
/// `None` is QPA's `false`. At bound zero, QPA returns the boolean `true`
/// exactly when the minimal right add(M)-approximation is an isomorphism. A
/// negative bound computes that initial approximation and then fails when
/// `FiniteComplex` indexes the empty GAP range `[0..n]`. For a positive bound,
/// QPA requires each approximation to be surjective, appends `f*g`, and
/// continues from its kernel. It eagerly computes the next approximation
/// after a nonzero final kernel before returning `FiniteComplex(cat, 1, ...)`.
///
/// # Errors
///
/// Rejects modules over different algebras and propagates approximation,
/// composition, kernel, identity-map, and zero-module failures.
pub fn have_finite_resolution_in_add_m<F: QpaCoefficientField>(
    target: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<Option<FiniteAddMResolutionResult<F>>, EnvelopingModuleError> {
    require_common_approximation_algebra(target, module)?;
    let mut differentials = Vec::new();
    let mut connecting = identity_module_mapping(target.clone())?;
    let mut approximation = minimal_right_approximation(module, target)?;

    if bound == &Integer::from(0) {
        return Ok(approximation
            .is_isomorphism()
            .then_some(FiniteAddMResolutionResult::InAddM));
    }
    if bound < &Integer::from(0) {
        return Err(EnvelopingModuleError::EmptyResolutionDifferentialList {
            bound: bound.clone(),
        });
    }

    let mut degree = Integer::from(0);
    loop {
        if !approximation.is_surjective() {
            return Ok(None);
        }
        differentials.push(approximation.compose(&connecting)?);
        connecting = kernel_inclusion(&approximation)?;
        if connecting.source().dimension() == 0 {
            break;
        }
        approximation = minimal_right_approximation(module, connecting.source())?;
        if &degree == bound {
            break;
        }
        degree += 1;
    }

    Ok(Some(FiniteAddMResolutionResult::Resolution(
        AddMResolution {
            differentials,
            zero: zero_module(module.right_acting_algebra())?,
        },
    )))
}

/// Return QPA's installed `IsCotiltingModule(M)` attribute, if bound.
///
/// QPA installs `false` only on the self-orthogonality failure branch and
/// installs `true` after every indecomposable injective has been resolved.
#[must_use]
pub fn is_cotilting_module<F: Field>(module: &PathAlgebraModule<F>) -> Option<bool> {
    module.cached_is_cotilting_module()
}

/// QPA's successful `[m, resolution]` value from `CotiltingModule(M, n)`.
#[derive(Clone, Debug)]
pub struct CotiltingModuleResult<F: Field> {
    injective_dimension: usize,
    resolutions: Vec<FiniteAddMResolutionResult<F>>,
}

impl<F: Field> CotiltingModuleResult<F> {
    /// Return QPA's first result component, the injective dimension of `M`.
    #[must_use]
    pub const fn injective_dimension(&self) -> usize {
        self.injective_dimension
    }

    /// Return the injective-indexed add(M)-resolutions in QPA order.
    #[must_use]
    pub fn resolutions(&self) -> &[FiniteAddMResolutionResult<F>] {
        &self.resolutions
    }
}

/// Translate QPA's bounded `CotiltingModule(M, n)` operation.
///
/// `None` is QPA's `false`. The source first establishes injective dimension
/// at most `n`, then performs the same literal syzygy-based self-orthogonality
/// loop as `TiltingModule`: for dimension `m`, it makes `m + 1` Ext calls but
/// advances the syzygy only while the one-based counter is strictly less than
/// `m`. This intentionally repeats the final Ext degree. A nonzero Ext binds
/// `IsCotiltingModule(M) := false`. Finally, every indecomposable injective
/// must have a finite add(M)-resolution of bound `m`; success binds the
/// attribute to `true`.
///
/// # Errors
///
/// Propagates bounded injective-dimension, Ext, syzygy, indecomposable
/// injective, and finite-resolution failures. QPA also assumes that a
/// successful bounded dimension call installed `InjDimension`; the source's
/// unbound-attribute failure is represented explicitly.
pub fn cotilting_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<Option<CotiltingModuleResult<F>>, EnvelopingModuleError> {
    match inj_dimension(module) {
        Some(dimension) => {
            if &Integer::from(dimension) > bound {
                return Ok(None);
            }
        }
        None => {
            if inj_dimension_of_module(module, bound)?.is_none() {
                return Ok(None);
            }
        }
    }
    let dimension = inj_dimension(module)
        .ok_or(EnvelopingModuleError::MissingCotiltingInjectiveDimensionCache)?;

    let mut syzygy = module.clone();
    let syzygy_cutoff = dimension.saturating_sub(1);
    for iteration in 0..=dimension {
        if ext_dimension_over_algebra(&syzygy, module)? != 0 {
            let _ = module.cache_is_cotilting_module(false);
            return Ok(None);
        }
        if iteration < syzygy_cutoff {
            syzygy = nth_syzygy(&syzygy, &Integer::one())?;
        }
    }

    let resolution_bound = Integer::from(dimension);
    let mut resolutions = Vec::new();
    for injective in indec_injective_modules(module.right_acting_algebra())? {
        let Some(resolution) =
            have_finite_resolution_in_add_m(&injective, module, &resolution_bound)?
        else {
            return Ok(None);
        };
        resolutions.push(resolution);
    }
    let _ = module.cache_is_cotilting_module(true);
    Ok(Some(CotiltingModuleResult {
        injective_dimension: dimension,
        resolutions,
    }))
}

/// QPA's pair of possibly empty complement sequences.
///
/// `None` preserves QPA's literal empty list for that side. A nonempty value
/// is the corresponding `FiniteComplex` constructed by the source method.
#[derive(Clone, Debug)]
pub struct AllComplementsResult<F: Field> {
    resolution: Option<AddMResolution<F>>,
    coresolution: Option<AddMCoresolution<F>>,
}

impl<F: Field> AllComplementsResult<F> {
    /// Return QPA's first result component.
    #[must_use]
    pub const fn resolution(&self) -> Option<&AddMResolution<F>> {
        self.resolution.as_ref()
    }

    /// Return QPA's second result component.
    #[must_use]
    pub const fn coresolution(&self) -> Option<&AddMCoresolution<F>> {
        self.coresolution.as_ref()
    }
}

/// Translate `AllComplementsOfAlmostCompleteTiltingModule(M, X)`.
///
/// QPA first follows successive surjective minimal right add(M)
/// approximations of `X`. If at least one exists, it appends the kernel
/// inclusion of the final composite and creates a degree-one finite complex.
/// It independently follows injective minimal left approximations, appends
/// the cokernel projection of the final composite, reverses the maps, and
/// creates the source's negative-base finite complex. A side with no accepted
/// approximation remains QPA's empty list.
///
/// # Errors
///
/// Rejects modules over different algebras and propagates approximation,
/// composition, kernel, cokernel, identity-map, and zero-module failures.
pub fn all_complements_of_almost_complete_tilting_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    complement: &PathAlgebraModule<F>,
) -> Result<AllComplementsResult<F>, EnvelopingModuleError> {
    require_common_approximation_algebra(module, complement)?;

    let mut resolution_differentials = Vec::new();
    let mut connecting = identity_module_mapping(complement.clone())?;
    let mut approximation = minimal_right_approximation(module, complement)?;
    while approximation.is_surjective() {
        resolution_differentials.push(approximation.compose(&connecting)?);
        connecting = kernel_inclusion(&approximation)?;
        if connecting.source().dimension() == 0 {
            break;
        }
        approximation = minimal_right_approximation(module, connecting.source())?;
    }
    let resolution = if let Some(last) = resolution_differentials.last().cloned() {
        let terminal = kernel_inclusion(&last)?;
        resolution_differentials.push(terminal);
        Some(AddMResolution {
            differentials: resolution_differentials,
            zero: zero_module(module.right_acting_algebra())?,
        })
    } else {
        None
    };

    let mut coresolution_differentials = Vec::new();
    let mut connecting = identity_module_mapping(complement.clone())?;
    let mut approximation = minimal_left_approximation(complement, module)?;
    while approximation.is_injective() {
        coresolution_differentials.push(connecting.compose(&approximation)?);
        connecting = cokernel_projection(&approximation)?;
        approximation = minimal_left_approximation(connecting.range(), module)?;
    }
    let coresolution = if let Some(last) = coresolution_differentials.last().cloned() {
        let terminal = cokernel_projection(&last)?;
        coresolution_differentials.push(terminal);
        coresolution_differentials.reverse();
        Some(AddMCoresolution {
            base_position: Integer::one() - Integer::from(coresolution_differentials.len()),
            differentials: coresolution_differentials,
            zero: zero_module(module.right_acting_algebra())?,
        })
    } else {
        None
    };

    Ok(AllComplementsResult {
        resolution,
        coresolution,
    })
}

/// QPA's synonym `AllComplementsOfAlmostCompleteCotiltingModule`.
///
/// # Errors
///
/// Propagates [`all_complements_of_almost_complete_tilting_module`].
pub fn all_complements_of_almost_complete_cotilting_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    complement: &PathAlgebraModule<F>,
) -> Result<AllComplementsResult<F>, EnvelopingModuleError> {
    all_complements_of_almost_complete_tilting_module(module, complement)
}

/// QPA's integer-or-infinity `FaithfulDimension(M)` value.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FaithfulDimensionResult {
    /// A finite faithful dimension.
    Finite(usize),
    /// QPA's `infinity` value.
    Infinite,
}

impl FaithfulDimensionResult {
    fn from_cache(value: Option<usize>) -> Self {
        value.map_or(Self::Infinite, Self::Finite)
    }
}

/// Translate QPA's cached `FaithfulDimension(M)` attribute.
///
/// For every indecomposable projective, the source follows successive
/// injective minimal left add(M)-approximations. A failed initial injection
/// returns zero immediately. Reaching a zero cokernel gives that projective
/// infinite length; otherwise the first failed successor records the number
/// of accepted injections. The result is the minimum of those values, with
/// QPA's usual ordering in which every finite integer is below `infinity`.
/// The attribute result is retained on the identical module object.
///
/// # Errors
///
/// Propagates indecomposable-projective, approximation, and cokernel failures.
/// Reports QPA's `Minimum([])` failure explicitly for an algebra with no
/// indecomposable projectives, and reports a non-addressable chain length.
pub fn faithful_dimension<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<FaithfulDimensionResult, EnvelopingModuleError> {
    if let Some(cached) = module.cached_faithful_dimension() {
        return Ok(FaithfulDimensionResult::from_cache(cached));
    }

    let projectives = indec_projective_modules(module.right_acting_algebra())?;
    if projectives.is_empty() {
        return Err(EnvelopingModuleError::EmptyFaithfulDimensionProjectiveList);
    }
    let mut finite_minimum = None;
    for projective in projectives {
        let mut length = 0_usize;
        let mut approximation = minimal_left_approximation(&projective, module)?;
        if !approximation.is_injective() {
            let cached = module.cache_faithful_dimension(Some(0));
            return Ok(FaithfulDimensionResult::from_cache(cached));
        }
        loop {
            length = length
                .checked_add(1)
                .ok_or(EnvelopingModuleError::FaithfulDimensionOverflow)?;
            let cokernel = cokernel_projection(&approximation)?;
            if cokernel.range().dimension() == 0 {
                break;
            }
            approximation = minimal_left_approximation(cokernel.range(), module)?;
            if !approximation.is_injective() {
                finite_minimum =
                    Some(finite_minimum.map_or(length, |value: usize| value.min(length)));
                break;
            }
        }
    }

    let cached = module.cache_faithful_dimension(finite_minimum);
    Ok(FaithfulDimensionResult::from_cache(cached))
}

/// QPA's integer-or-infinity complement count.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NumberOfComplementsResult {
    /// A finite number of complements.
    Finite(usize),
    /// QPA's `infinity` value.
    Infinite,
}

/// Translate `NumberOfComplementsOfAlmostCompleteTiltingModule(M)`.
///
/// The source returns `FaithfulDimension(M) + 1`, including GAP's arithmetic
/// rule `infinity + 1 = infinity`.
///
/// # Errors
///
/// Propagates [`faithful_dimension`] and reports a finite result beyond
/// Rust's addressable range.
pub fn number_of_complements_of_almost_complete_tilting_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<NumberOfComplementsResult, EnvelopingModuleError> {
    match faithful_dimension(module)? {
        FaithfulDimensionResult::Finite(value) => value
            .checked_add(1)
            .map(NumberOfComplementsResult::Finite)
            .ok_or(EnvelopingModuleError::NumberOfComplementsOverflow),
        FaithfulDimensionResult::Infinite => Ok(NumberOfComplementsResult::Infinite),
    }
}

/// QPA's synonym `NumberOfComplementsOfAlmostCompleteCotiltingModule`.
///
/// # Errors
///
/// Propagates [`number_of_complements_of_almost_complete_tilting_module`].
pub fn number_of_complements_of_almost_complete_cotilting_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<NumberOfComplementsResult, EnvelopingModuleError> {
    number_of_complements_of_almost_complete_tilting_module(module)
}

/// Translate `LeftMutationOfTiltingModuleComplement(M, N)`.
///
/// QPA constructs the minimal left add(M)-approximation of `N` and returns
/// its cokernel exactly when that approximation is injective. `None` is the
/// source's `false` result.
///
/// # Errors
///
/// Propagates minimal-approximation and cokernel construction failures,
/// including the source's different-algebra rejection.
pub fn left_mutation_of_tilting_module_complement<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    complement: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    let approximation = minimal_left_approximation(complement, module)?;
    if approximation.is_injective() {
        return Ok(Some(cokernel_of_module_homomorphism(&approximation)?));
    }
    Ok(None)
}

/// Translate `RightMutationOfTiltingModuleComplement(M, N)`.
///
/// QPA constructs the minimal right add(M)-approximation of `N` and returns
/// its kernel exactly when that approximation is surjective. `None` is the
/// source's `false` result.
///
/// # Errors
///
/// Propagates minimal-approximation and kernel construction failures,
/// including the source's different-algebra rejection.
pub fn right_mutation_of_tilting_module_complement<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    complement: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    let approximation = minimal_right_approximation(module, complement)?;
    if approximation.is_surjective() {
        return Ok(Some(kernel_of_module_homomorphism(&approximation)?));
    }
    Ok(None)
}

/// QPA's synonym `LeftMutationOfCotiltingModuleComplement`.
///
/// # Errors
///
/// Propagates [`left_mutation_of_tilting_module_complement`].
pub fn left_mutation_of_cotilting_module_complement<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    complement: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    left_mutation_of_tilting_module_complement(module, complement)
}

/// QPA's synonym `RightMutationOfCotiltingModuleComplement`.
///
/// # Errors
///
/// Propagates [`right_mutation_of_tilting_module_complement`].
pub fn right_mutation_of_cotilting_module_complement<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    complement: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    right_mutation_of_tilting_module_complement(module, complement)
}

fn require_cotilting_for_approximation<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<(), EnvelopingModuleError> {
    match is_cotilting_module(module) {
        Some(true) => Ok(()),
        Some(false) => Err(EnvelopingModuleError::FirstApproximationArgumentNotCotilting),
        None => Err(EnvelopingModuleError::CotiltingApproximationPropertyUnbound),
    }
}

fn cotilting_approximation_exact_sequences<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    dimension: usize,
) -> Result<Vec<PathAlgebraModuleHomomorphismPair<F>>, EnvelopingModuleError> {
    let resolution = projective_resolution(module)?;
    let _ = resolution.object(&Integer::from(dimension))?;
    let mut exact_sequences = Vec::with_capacity(dimension);
    for degree in 0..dimension {
        let differential = resolution.differential(&Integer::from(degree))?;
        exact_sequences.push((
            kernel_inclusion(&differential)?,
            image_projection(&differential)?,
        ));
    }
    if exact_sequences.is_empty() {
        return Err(EnvelopingModuleError::EmptyCotiltingApproximationExactSequences);
    }
    Ok(exact_sequences)
}

fn cotilting_approximation_isomorphism<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    isomorphism_of_modules(source, range)?
        .ok_or(EnvelopingModuleError::InternalCotiltingApproximationIsomorphismFailed)
}

fn cotilting_approximation_push_out<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    Ok(push_out(first, second)?
        .ok_or(EnvelopingModuleError::InternalCotiltingApproximationPushOutFailed)?
        .0)
}
