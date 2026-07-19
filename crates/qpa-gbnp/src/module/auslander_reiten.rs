/// Failure while evaluating QPA's `IsExceptionalModule` property.
#[derive(Debug, Error)]
pub enum ExceptionalModuleError {
    /// Rigidity and Ext construction failed.
    #[error(transparent)]
    Enveloping(#[from] EnvelopingModuleError),
    /// Indecomposable decomposition failed.
    #[error(transparent)]
    Decomposition(#[from] ModuleDecompositionError),
}
/// Choose whether QPA's two-argument `AlmostSplitSequence(M, e)` dispatch
/// starts at or ends in the entered module.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlmostSplitSide {
    /// QPA's `"l"`: construct a sequence starting at the module by duality.
    Starting,
    /// QPA's `"r"`: construct the cached sequence ending in the module.
    Ending,
}

/// Translate QPA's cached `AlmostSplitSequence(M)` attribute.
///
/// For a projective module, `Ok(None)` is QPA's `fail`. Otherwise the source
/// chooses the first basis class outside the image of
/// `Hom(P(M), DTr(M))`, repeatedly acts by the radical of `End(DTr(M))`
/// until that Ext class reaches the socle, and realizes the class as a
/// pushout. The returned pair is the monomorphism followed by the epimorphism.
///
/// # Errors
///
/// Propagates projectivity, `DTr`, projective-cover, Hom-space,
/// endomorphism-radical, exact row-space, pushout, and cokernel failures.
pub fn almost_split_sequence<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModuleHomomorphismPair<F>>, AuslanderReitenError> {
    match module.cached_almost_split_sequence() {
        CachedAttribute::Unbound => {}
        CachedAttribute::Fail => return Ok(None),
        CachedAttribute::False => {
            return Err(AuslanderReitenError::InvalidAlmostSplitCacheState);
        }
        CachedAttribute::Value(sequence) => return Ok(Some(sequence)),
    }
    if is_projective_module(module)? {
        return Ok(module.cache_almost_split_sequence(None));
    }

    let translate = d_tr(module)?;
    let cover = projective_cover(module)?;
    let syzygy_inclusion = kernel_inclusion(&cover)?;
    let projective_homs = hom_over_algebra(cover.source(), &translate)?
        .ok_or(AuslanderReitenError::InternalProjectiveHomRejected)?;
    let syzygy_homs = hom_over_algebra(syzygy_inclusion.source(), &translate)?
        .ok_or(AuslanderReitenError::InternalSyzygyHomRejected)?;
    let Some(first_syzygy_hom) = syzygy_homs.first() else {
        return Err(AuslanderReitenError::NoNonzeroExtClass);
    };
    let flat_length = first_syzygy_hom.flat_vertex_matrix_entries().len();
    let image_generators = projective_homs
        .iter()
        .map(|homomorphism| syzygy_inclusion.compose(homomorphism))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .filter(|homomorphism| !homomorphism.is_zero())
        .collect::<Vec<_>>();
    let image_span = if image_generators.is_empty() {
        None
    } else {
        Some(DenseMatrix::from_flat(
            image_generators.len(),
            flat_length,
            image_generators
                .iter()
                .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
                .collect(),
            module.right_acting_algebra().field().clone(),
        )?)
    };
    let mut extension = syzygy_homs
        .into_iter()
        .find_map(|homomorphism| {
            match flat_row_in_span(
                image_span.as_ref(),
                &homomorphism.flat_vertex_matrix_entries(),
                module.right_acting_algebra().field(),
            ) {
                Ok(true) => None,
                Ok(false) => Some(Ok(homomorphism)),
                Err(error) => Some(Err(error)),
            }
        })
        .transpose()?
        .ok_or(AuslanderReitenError::NoNonzeroExtClass)?;

    let endomorphism_algebra = end_over_algebra(&translate)?;
    let structure_constants = endomorphism_algebra.structure_constants_algebra()?;
    let radical = structure_constants.radical()?;
    let radical_homomorphisms = (0..radical.basis().nrows())
        .map(|row| {
            let matrix = endomorphism_algebra.matrix_from_coefficients(radical.basis().row(row))?;
            Ok(from_end_m_to_hom_mm(&translate, &matrix)?)
        })
        .collect::<Result<Vec<_>, AuslanderReitenError>>()?;

    loop {
        let composites = radical_homomorphisms
            .iter()
            .map(|radical_map| extension.compose(radical_map))
            .collect::<Result<Vec<_>, _>>()?;
        let mut first_outside = None;
        for (position, composite) in composites.iter().enumerate() {
            if !flat_row_in_span(
                image_span.as_ref(),
                &composite.flat_vertex_matrix_entries(),
                module.right_acting_algebra().field(),
            )? {
                first_outside = Some(position);
                break;
            }
        }
        let Some(position) = first_outside else {
            break;
        };
        extension = composites[position].clone();
    }

    let pushout = push_out(&syzygy_inclusion, &extension)?
        .ok_or(AuslanderReitenError::InternalAlmostSplitPushOutFailed)?;
    let monomorphism = pushout.0;
    let epimorphism = cokernel_projection(&monomorphism)?;
    let sequence = (monomorphism, epimorphism);
    Ok(module.cache_almost_split_sequence(Some(&sequence)))
}

/// Translate QPA's two-argument `AlmostSplitSequence(M, "l"|"r")`
/// overload with a typed side selector.
///
/// # Errors
///
/// The ending case delegates to [`almost_split_sequence`]. The starting case
/// additionally propagates module and map dualization failures.
pub fn almost_split_sequence_from_side<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    side: AlmostSplitSide,
) -> Result<Option<PathAlgebraModuleHomomorphismPair<F>>, AuslanderReitenError> {
    if side == AlmostSplitSide::Ending {
        return almost_split_sequence(module);
    }
    let dual = dual_of_module(module)?;
    let Some((monomorphism, epimorphism)) = almost_split_sequence(&dual)? else {
        return Ok(None);
    };
    Ok(Some((
        dual_of_module_homomorphism(&epimorphism)?,
        dual_of_module_homomorphism(&monomorphism)?,
    )))
}

/// Translate QPA's `AlmostSplitSequenceInPerpT(T, M)` operation.
///
/// The source first requires only that `IsCotiltingModule(T)` has already
/// been installed. It returns `fail` for projective `M`; otherwise it takes
/// the almost-split epimorphism ending in `M`, precomposes it with
/// `RightApproximationByPerpT(T, Source(epi))`, applies QPA's right-minimal
/// version, and returns that epimorphism with its kernel inclusion.
///
/// # Errors
///
/// Preserves the unbound cotilting-property error and propagates
/// almost-split, perpendicular-approximation, composition, minimal-version,
/// and kernel failures.
pub fn almost_split_sequence_in_perp_t<F: QpaCoefficientField>(
    cotilting: &PathAlgebraModule<F>,
    module: &PathAlgebraModule<F>,
) -> Result<Option<PathAlgebraModuleHomomorphismPair<F>>, AuslanderReitenError> {
    if cotilting.cached_is_cotilting_module().is_none() {
        return Err(AuslanderReitenError::CotiltingPerpPropertyUnbound);
    }
    if is_projective_module(module)? {
        return Ok(None);
    }
    let (_, epimorphism) = almost_split_sequence(module)?
        .ok_or(AuslanderReitenError::UnexpectedProjectiveAlmostSplitFail)?;
    let approximation = right_approximation_by_perp_t(cotilting, epimorphism.source())?;
    let composite = approximation.compose(&epimorphism)?;
    let minimal = right_minimal_version(&composite)?.0;
    Ok(Some((kernel_inclusion(&minimal)?, minimal)))
}

/// Translate QPA's cached `IrreducibleMorphismsEndingIn(M)` attribute.
///
/// Over a finite field, the source uses the radical inclusion for a
/// projective module and the almost-split epimorphism otherwise. It decomposes
/// that map's source with inclusions and composes every inclusion into the
/// selected map, preserving decomposition order and repeated isomorphic
/// summands.
///
/// # Errors
///
/// Requires a supported finite field and propagates projectivity,
/// almost-split, radical, decomposition, and map-composition failures.
pub fn irreducible_morphisms_ending_in<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, AuslanderReitenError> {
    if let Some(cached) = module.cached_irreducible_morphisms_ending_in() {
        return Ok(cached);
    }
    if module
        .right_acting_algebra()
        .field()
        .qpa_elements()?
        .is_none()
    {
        return Err(AuslanderReitenError::FiniteFieldRequired);
    }
    let right_almost_split = if is_projective_module(module)? {
        radical_of_module_inclusion(module.clone())?
    } else {
        almost_split_sequence(module)?
            .ok_or(AuslanderReitenError::UnexpectedProjectiveAlmostSplitFail)?
            .1
    };
    if right_almost_split.source().dimension() == 0 {
        return Err(AuslanderReitenError::ZeroIrreducibleSourceDecomposition);
    }
    let inclusions = decompose_module_with_inclusions(right_almost_split.source())?;
    let homomorphisms = inclusions
        .iter()
        .map(|inclusion| inclusion.compose(&right_almost_split))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(module.cache_irreducible_morphisms_ending_in(&homomorphisms))
}

/// Translate QPA's cached `IrreducibleMorphismsStartingIn(M)` attribute by
/// dualizing every irreducible morphism ending in `DualOfModule(M)`.
///
/// # Errors
///
/// Requires a supported finite field and propagates dualization and ending-map
/// construction failures.
pub fn irreducible_morphisms_starting_in<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, AuslanderReitenError> {
    if let Some(cached) = module.cached_irreducible_morphisms_starting_in() {
        return Ok(cached);
    }
    if module
        .right_acting_algebra()
        .field()
        .qpa_elements()?
        .is_none()
    {
        return Err(AuslanderReitenError::FiniteFieldRequired);
    }
    let dual = dual_of_module(module)?;
    let homomorphisms = irreducible_morphisms_ending_in(&dual)?
        .iter()
        .map(dual_of_module_homomorphism)
        .collect::<Result<Vec<_>, _>>()?;
    Ok(module.cache_irreducible_morphisms_starting_in(&homomorphisms))
}

/// QPA's heterogeneous third component in a predecessor-arrow valuation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum PredecessorMultiplicity {
    /// The usual `[integer|false, integer|false]` value.
    Pair {
        /// QPA's first list position.
        first: Option<usize>,
        /// QPA's second list position.
        second: Option<usize>,
    },
    /// The projective/new-summand branch literally stores an integer instead
    /// of the usual pair.
    Scalar(usize),
}

/// One QPA predecessor valuation entry `[source_position, range_position, v]`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PredecessorValuation {
    /// One-based position in the farther predecessor layer.
    pub predecessor_position: usize,
    /// One-based position in the layer nearer the entered module.
    pub successor_position: usize,
    /// QPA's pair-or-scalar valuation payload.
    pub multiplicity: PredecessorMultiplicity,
}

/// Layer and valuation lists returned by QPA's general predecessor branch.
#[derive(Clone, Debug)]
pub struct PredecessorLayers<F: Field> {
    /// The `n + 1` module layers in QPA order.
    pub layers: Vec<Vec<PathAlgebraModule<F>>>,
    /// The `n` arrow-valuation layers in QPA order.
    pub valuations: Vec<Vec<PredecessorValuation>>,
}

/// QPA's heterogeneous `PredecessorsOfModule(M,n)` return value.
#[derive(Clone, Debug)]
pub enum PredecessorsOfModuleResult<F: Field> {
    /// The literal `n = 1` early return is the entered module itself.
    Module(PathAlgebraModule<F>),
    /// Every supported `n >= 2` result contains layers and valuations.
    Layers(PredecessorLayers<F>),
}

/// Translate QPA's `PredecessorsOfModule(M, n)` layer algorithm.
///
/// The source assumes indecomposability, rejects a projective module before
/// inspecting `n`, returns `M` itself at `n = 1`, and for larger bounds walks
/// almost-split middle terms. Isomorphic occurrences are detected through
/// `CommonDirectSummand`, while every stored position remains one-based.
///
/// # Errors
///
/// Rejects projective input and nonpositive or unrepresentable bounds, then
/// propagates almost-split, radical, decomposition, and common-summand
/// failures. Internal valuation lookup/shape failures mirror source indexing
/// errors explicitly.
#[allow(clippy::if_not_else, clippy::too_many_lines)]
pub fn predecessors_of_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<PredecessorsOfModuleResult<F>, AuslanderReitenError> {
    if is_projective_module(module)? {
        return Err(AuslanderReitenError::PredecessorModuleProjective);
    }
    if bound == &Integer::from(1) {
        return Ok(PredecessorsOfModuleResult::Module(module.clone()));
    }
    if bound < &Integer::from(2) {
        return Err(AuslanderReitenError::PredecessorBoundBelowOne {
            bound: bound.clone(),
        });
    }
    let bound = usize::try_from(bound.clone()).map_err(|_| {
        AuslanderReitenError::PredecessorBoundNotRepresentable {
            bound: bound.clone(),
        }
    })?;
    let layer_count = bound
        .checked_add(1)
        .ok_or(AuslanderReitenError::PredecessorLayerCountOverflow { bound })?;
    let mut layers = vec![Vec::new(); layer_count];
    let mut valuations = vec![Vec::new(); bound];
    layers[0].push(module.clone());

    let first_sequence = almost_split_sequence(module)?
        .ok_or(AuslanderReitenError::UnexpectedProjectiveAlmostSplitFail)?;
    let first_middle = decompose_module_with_multiplicities(first_sequence.0.range())?
        .ok_or(AuslanderReitenError::EmptyAlmostSplitMiddleTerm)?;
    layers[1].extend(first_middle.basic_summands().iter().cloned());
    layers[2].push(first_sequence.0.source().clone());
    for (position, &multiplicity) in first_middle.multiplicities().iter().enumerate() {
        valuations[0].push(PredecessorValuation {
            predecessor_position: position + 1,
            successor_position: 1,
            multiplicity: PredecessorMultiplicity::Pair {
                first: Some(multiplicity),
                second: None,
            },
        });
        valuations[1].push(PredecessorValuation {
            predecessor_position: 1,
            successor_position: position + 1,
            multiplicity: PredecessorMultiplicity::Pair {
                first: None,
                second: Some(multiplicity),
            },
        });
    }

    for layer_number in 2..bound {
        let current_modules = layers[layer_number - 1].clone();
        for current in current_modules {
            let current_position = layers[layer_number - 1]
                .iter()
                .position(|candidate| candidate == &current)
                .ok_or(AuslanderReitenError::MissingCurrentPredecessorModule)?
                + 1;
            if !is_projective_module(&current)? {
                let sequence = almost_split_sequence(&current)?
                    .ok_or(AuslanderReitenError::UnexpectedProjectiveAlmostSplitFail)?;
                let middle = decompose_module_with_multiplicities(sequence.0.range())?
                    .ok_or(AuslanderReitenError::EmptyAlmostSplitMiddleTerm)?;
                layers[layer_number + 1].push(sequence.0.source().clone());
                let translate_position = layers[layer_number + 1].len();

                for (summand, &multiplicity) in
                    middle.basic_summands().iter().zip(middle.multiplicities())
                {
                    let mut found = false;
                    for next_position in 0..layers[layer_number].len() {
                        if matches!(
                            common_direct_summand(summand, &layers[layer_number][next_position])?,
                            CommonDirectSummandResult::Found(_)
                        ) {
                            update_predecessor_forward_multiplicity(
                                &mut valuations[layer_number - 1],
                                next_position + 1,
                                current_position,
                                multiplicity,
                            )?;
                            valuations[layer_number].push(PredecessorValuation {
                                predecessor_position: translate_position,
                                successor_position: next_position + 1,
                                multiplicity: PredecessorMultiplicity::Pair {
                                    first: None,
                                    second: Some(multiplicity),
                                },
                            });
                            found = true;
                        }
                    }
                    if !found {
                        layers[layer_number].push(summand.clone());
                        let next_position = layers[layer_number].len();
                        valuations[layer_number - 1].push(PredecessorValuation {
                            predecessor_position: next_position,
                            successor_position: current_position,
                            multiplicity: PredecessorMultiplicity::Pair {
                                first: Some(multiplicity),
                                second: None,
                            },
                        });
                        valuations[layer_number].push(PredecessorValuation {
                            predecessor_position: translate_position,
                            successor_position: next_position,
                            multiplicity: PredecessorMultiplicity::Pair {
                                first: None,
                                second: Some(multiplicity),
                            },
                        });
                    }
                }
            } else {
                let radical = radical_of_module(current.clone())?;
                if radical.dimension() == 0 {
                    continue;
                }
                let middle = decompose_module_with_multiplicities(&radical)?
                    .ok_or(AuslanderReitenError::EmptyProjectiveRadicalTerm)?;
                for (summand, &multiplicity) in
                    middle.basic_summands().iter().zip(middle.multiplicities())
                {
                    let mut found = false;
                    for next_position in 0..layers[layer_number].len() {
                        if matches!(
                            common_direct_summand(summand, &layers[layer_number][next_position])?,
                            CommonDirectSummandResult::Found(_)
                        ) {
                            update_predecessor_forward_multiplicity(
                                &mut valuations[layer_number - 1],
                                next_position + 1,
                                current_position,
                                multiplicity,
                            )?;
                            found = true;
                        }
                    }
                    if !found {
                        layers[layer_number].push(summand.clone());
                        valuations[layer_number - 1].push(PredecessorValuation {
                            predecessor_position: layers[layer_number].len(),
                            successor_position: current_position,
                            multiplicity: PredecessorMultiplicity::Scalar(multiplicity),
                        });
                    }
                }
            }
        }
    }

    Ok(PredecessorsOfModuleResult::Layers(PredecessorLayers {
        layers,
        valuations,
    }))
}

fn update_predecessor_forward_multiplicity(
    valuations: &mut [PredecessorValuation],
    predecessor_position: usize,
    successor_position: usize,
    multiplicity: usize,
) -> Result<(), AuslanderReitenError> {
    let valuation = valuations
        .iter_mut()
        .find(|valuation| {
            valuation.predecessor_position == predecessor_position
                && valuation.successor_position == successor_position
        })
        .ok_or(AuslanderReitenError::MissingPredecessorValuation {
            predecessor_position,
            successor_position,
        })?;
    let PredecessorMultiplicity::Pair { first, .. } = &mut valuation.multiplicity else {
        return Err(AuslanderReitenError::ScalarPredecessorValuationUpdate);
    };
    *first = Some(multiplicity);
    Ok(())
}

fn flat_row_in_span<F: Field>(
    span: Option<&DenseMatrix<F>>,
    row: &[F::Element],
    field: &F,
) -> Result<bool, MatrixError> {
    match span {
        Some(span) => Ok(span.solution_row(row)?.is_some()),
        None => Ok(row.iter().all(|entry| field.is_zero(entry))),
    }
}

/// Failure while constructing QPA's almost-split sequences and neighboring
/// Auslander--Reiten data.
#[derive(Debug, Error)]
pub enum AuslanderReitenError {
    /// QPA requires users to run `CotiltingModule(T, n)` before this method.
    #[error("the first argument has no installed cotilting-module property")]
    CotiltingPerpPropertyUnbound,
    /// The attribute accepts only unbound, `fail`, or a two-map sequence.
    #[error("the AlmostSplitSequence attribute cache contains boolean false")]
    InvalidAlmostSplitCacheState,
    /// The projective-cover Hom call uses modules over one acting algebra.
    #[error("HomOverAlgebra rejected the projective-to-DTr Hom space")]
    InternalProjectiveHomRejected,
    /// The syzygy Hom call uses modules over one acting algebra.
    #[error("HomOverAlgebra rejected the syzygy-to-DTr Hom space")]
    InternalSyzygyHomRejected,
    /// An indecomposable nonprojective module must have nonzero
    /// `Ext^1(M,DTr(M))`.
    #[error("AlmostSplitSequence found no nonzero Ext class")]
    NoNonzeroExtClass,
    /// Both source maps in the almost-split pushout share the syzygy source.
    #[error("the internal AlmostSplitSequence pushout returned fail")]
    InternalAlmostSplitPushOutFailed,
    /// QPA restricts irreducible-morphism enumeration to finite fields.
    #[error("module is not over a quiver algebra with a finite field as ground ring")]
    FiniteFieldRequired,
    /// The nonprojective branch cannot receive the projective `fail` result.
    #[error("AlmostSplitSequence unexpectedly returned fail for a nonprojective module")]
    UnexpectedProjectiveAlmostSplitFail,
    /// QPA's decomposition route errors when a projective simple has zero
    /// radical (and dually for a simple injective starting map).
    #[error("IrreducibleMorphisms attempts to decompose a zero map source")]
    ZeroIrreducibleSourceDecomposition,
    /// Supported-field enumeration failed.
    #[error(transparent)]
    Coefficient(#[from] CoefficientDomainError),
    /// Enveloping-module and AR-functor construction failed.
    #[error(transparent)]
    Enveloping(#[from] EnvelopingModuleError),
    /// Generic module or map dualization failed.
    #[error(transparent)]
    Functor(#[from] ModuleFunctorError),
    /// Module maps, projective covers, pushouts, or cokernels failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// Exact row-space or endomorphism-matrix arithmetic failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// The radical of the exact endomorphism algebra failed.
    #[error(transparent)]
    StructureConstants(#[from] StructureConstantsError),
    /// Finite-field indecomposable decomposition failed.
    #[error(transparent)]
    Decomposition(#[from] ModuleDecompositionError),
    /// The source rejects projective input before its `n = 1` early return.
    #[error("entered module is projective")]
    PredecessorModuleProjective,
    /// Nonpositive bounds make QPA's subsequent layer indexing invalid.
    #[error("PredecessorsOfModule bound must be positive: {bound}")]
    PredecessorBoundBelowOne {
        /// Rejected arbitrary-precision bound.
        bound: Integer,
    },
    /// A positive layer bound must fit addressable list lengths.
    #[error("PredecessorsOfModule bound is not representable by usize: {bound}")]
    PredecessorBoundNotRepresentable {
        /// Rejected arbitrary-precision bound.
        bound: Integer,
    },
    /// QPA requests `n + 1` layer lists.
    #[error("PredecessorsOfModule layer count overflows usize for bound {bound}")]
    PredecessorLayerCountOverflow {
        /// Entered representable bound whose successor is not representable.
        bound: usize,
    },
    /// A nonprojective almost-split sequence has a nonzero middle term.
    #[error("an almost-split sequence has an empty middle-term decomposition")]
    EmptyAlmostSplitMiddleTerm,
    /// A nonzero projective radical has a nonempty decomposition.
    #[error("a nonzero projective radical has an empty decomposition")]
    EmptyProjectiveRadicalTerm,
    /// The iterated layer must still contain the module currently being read.
    #[error("the current predecessor module disappeared from its layer")]
    MissingCurrentPredecessorModule,
    /// QPA's `Position` lookup must find the already-created arrow valuation.
    #[error(
        "missing predecessor valuation at positions ({predecessor_position},{successor_position})"
    )]
    MissingPredecessorValuation {
        /// One-based farther-layer position.
        predecessor_position: usize,
        /// One-based nearer-layer position.
        successor_position: usize,
    },
    /// QPA indexes the first slot of a valuation here, which fails if an
    /// earlier projective branch stored its scalar payload.
    #[error("QPA attempted to update the first slot of a scalar predecessor valuation")]
    ScalarPredecessorValuationUpdate,
}

/// Translate QPA's `AllModulesOfLengthPlusOne(list)` operation.
///
/// The input is a nonempty list of equal-dimensional modules over one
/// finite-field algebra. QPA constructs extensions on both sides by every
/// vertex simple, retains first occurrences under structural equality, and
/// finally removes later isomorphic modules. The source's first `s = 1`
/// branch literally takes `list[u]`, the last entered module, rather than
/// `list[j]`; this observable indexing behavior is preserved here.
///
/// # Errors
///
/// Preserves QPA's three input errors and finite-field restriction, then
/// propagates simple-module, projective-cover, Ext, pushout, direct-sum, and
/// isomorphism failures.
pub fn all_modules_of_length_plus_one<F: QpaCoefficientField>(
    list: &[PathAlgebraModule<F>],
) -> Result<Vec<PathAlgebraModule<F>>, ModuleLengthError> {
    let Some(first) = list.first() else {
        return Err(ModuleLengthError::EmptyModuleList);
    };
    let algebra = first.right_acting_algebra();
    let field_elements = algebra
        .field()
        .qpa_elements()?
        .ok_or(ModuleLengthError::FiniteFieldRequired)?;
    if list
        .iter()
        .any(|module| module.right_acting_algebra() != algebra)
    {
        return Err(ModuleLengthError::ModulesOverDifferentAlgebras);
    }
    let dimension = first.dimension();
    if list.iter().any(|module| module.dimension() != dimension) {
        return Err(ModuleLengthError::ModulesOfDifferentLengths);
    }

    let simples = simple_modules(algebra)?;
    let simple_syzygies = simples
        .iter()
        .map(|simple| Ok(kernel_inclusion(&projective_cover(simple)?)?))
        .collect::<Result<Vec<_>, ModuleLengthError>>()?;
    let mut extensions = Vec::with_capacity(simples.len());
    for simple in &simples {
        let mut row = Vec::with_capacity(list.len());
        for module in list {
            row.push(ext_over_algebra(simple, module)?.basis().to_vec());
        }
        extensions.push(row);
    }

    let mut modules = Vec::new();
    let Some(source_first_split_module) = list.last() else {
        return Err(ModuleLengthError::EmptyModuleList);
    };
    for (simple_index, extension_row) in extensions.iter().enumerate() {
        for basis in extension_row {
            append_extension_modules(
                &mut modules,
                basis,
                source_first_split_module,
                &simples[simple_index],
                &simple_syzygies[simple_index],
                &field_elements,
            )?;
        }
    }
    retain_first_equal_modules(&mut modules);

    let module_syzygies = list
        .iter()
        .map(|module| Ok(kernel_inclusion(&projective_cover(module)?)?))
        .collect::<Result<Vec<_>, ModuleLengthError>>()?;
    let mut extensions = Vec::with_capacity(list.len());
    for module in list {
        let mut row = Vec::with_capacity(simples.len());
        for simple in &simples {
            row.push(ext_over_algebra(module, simple)?.basis().to_vec());
        }
        extensions.push(row);
    }
    for (module_index, extension_row) in extensions.iter().enumerate() {
        for (simple_index, basis) in extension_row.iter().enumerate() {
            append_extension_modules(
                &mut modules,
                basis,
                &list[module_index],
                &simples[simple_index],
                &module_syzygies[module_index],
                &field_elements,
            )?;
        }
    }
    retain_first_equal_modules(&mut modules);

    let Some(first) = modules.first().cloned() else {
        return Err(ModuleLengthError::EmptyConstructedModuleList);
    };
    let mut nonisomorphic = vec![first];
    for module in modules.into_iter().skip(1) {
        let mut new_isomorphism_class = true;
        for representative in &nonisomorphic {
            if isomorphic_modules(representative, &module)? {
                new_isomorphism_class = false;
                break;
            }
        }
        if new_isomorphism_class {
            nonisomorphic.push(module);
        }
    }
    Ok(nonisomorphic)
}

/// Translate QPA's `AllModulesOfLengthAtMost(A, n)` operation.
///
/// Despite its name, the source tests total vector-space dimension. Negative
/// bounds return an empty list; zero returns the zero module; later layers
/// are obtained by repeated [`all_modules_of_length_plus_one`] calls.
///
/// # Errors
///
/// Requires a supported finite field, reports a nonnegative bound too large
/// for Rust's address space, and propagates layer construction failures.
pub fn all_modules_of_length_at_most<F: QpaCoefficientField>(
    algebra: &QuiverAlgebra<F>,
    bound: &Integer,
) -> Result<Vec<PathAlgebraModule<F>>, ModuleLengthError> {
    if algebra.field().qpa_elements()?.is_none() {
        return Err(ModuleLengthError::FiniteFieldRequired);
    }
    if bound < &Integer::from(0) {
        return Ok(Vec::new());
    }
    let bound = usize::try_from(bound.clone()).map_err(|_| {
        ModuleLengthError::LengthBoundNotRepresentable {
            value: bound.clone(),
        }
    })?;
    let zero = zero_module(algebra)?;
    if bound == 0 {
        return Ok(vec![zero]);
    }
    let simples = simple_modules(algebra)?;
    let mut all_modules = Vec::with_capacity(simples.len() + 1);
    all_modules.push(zero);
    all_modules.extend(simples.iter().cloned());
    if bound == 1 {
        return Ok(all_modules);
    }

    let mut previous = simples;
    for _ in 2..=bound {
        previous = all_modules_of_length_plus_one(&previous)?;
        all_modules.extend(previous.iter().cloned());
    }
    Ok(all_modules)
}

/// Translate QPA's `AllIndecModulesOfLengthAtMost(A, n)` operation.
///
/// QPA removes the leading zero module from
/// [`all_modules_of_length_at_most`] and applies its finite-field
/// `IsIndecomposableModule` method to every remaining module.
///
/// # Errors
///
/// Propagates module enumeration and deterministic decomposition failures.
pub fn all_indec_modules_of_length_at_most<F: QpaCoefficientField>(
    algebra: &QuiverAlgebra<F>,
    bound: &Integer,
) -> Result<Vec<PathAlgebraModule<F>>, ModuleLengthError> {
    let all_modules = all_modules_of_length_at_most(algebra, bound)?;
    let mut indecomposables = Vec::new();
    for module in all_modules.into_iter().skip(1) {
        if is_indecomposable_module(&module)?
            .ok_or(ModuleLengthError::IndecomposabilityMethodUnavailable)?
        {
            indecomposables.push(module);
        }
    }
    Ok(indecomposables)
}

fn append_extension_modules<F: Field>(
    modules: &mut Vec<PathAlgebraModule<F>>,
    basis: &[PathAlgebraModuleHomomorphism<F>],
    split_module: &PathAlgebraModule<F>,
    simple: &PathAlgebraModule<F>,
    syzygy_inclusion: &PathAlgebraModuleHomomorphism<F>,
    field_elements: &[F::Element],
) -> Result<(), ModuleLengthError> {
    match basis.len() {
        0 => {}
        1 => {
            modules.push(
                direct_sum_of_qpa_modules(&[split_module.clone(), simple.clone()])?
                    .ok_or(ModuleLengthError::InternalDirectSumFailed)?,
            );
            modules.push(extension_middle_term(&basis[0], syzygy_inclusion)?);
        }
        dimension => {
            for coefficients in cartesian_field_power(field_elements, dimension) {
                let homomorphism = linear_combination_of_module_maps(basis, &coefficients)?;
                modules.push(extension_middle_term(&homomorphism, syzygy_inclusion)?);
            }
        }
    }
    Ok(())
}

fn extension_middle_term<F: Field>(
    extension: &PathAlgebraModuleHomomorphism<F>,
    syzygy_inclusion: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModule<F>, ModuleLengthError> {
    Ok(push_out(extension, syzygy_inclusion)?
        .ok_or(ModuleLengthError::InternalPushOutFailed)?
        .1
        .range()
        .clone())
}

fn linear_combination_of_module_maps<F: Field>(
    basis: &[PathAlgebraModuleHomomorphism<F>],
    coefficients: &[F::Element],
) -> Result<PathAlgebraModuleHomomorphism<F>, ModuleLengthError> {
    let first = basis
        .first()
        .ok_or(ModuleLengthError::EmptyLinearCombinationBasis)?;
    let field = first.source().right_acting_algebra().field();
    let mut result = first.scale(&field.zero())?;
    for (homomorphism, coefficient) in basis.iter().zip(coefficients) {
        result = result.add(&homomorphism.scale(coefficient)?)?;
    }
    Ok(result)
}

fn cartesian_field_power<E: Clone>(elements: &[E], exponent: usize) -> Vec<Vec<E>> {
    let mut product = vec![Vec::new()];
    for _ in 0..exponent {
        let mut next = Vec::new();
        for prefix in product {
            for element in elements {
                let mut tuple = prefix.clone();
                tuple.push(element.clone());
                next.push(tuple);
            }
        }
        product = next;
    }
    product
}

fn retain_first_equal_modules<F: Field>(modules: &mut Vec<PathAlgebraModule<F>>) {
    let mut unique = Vec::new();
    for module in modules.drain(..) {
        if !unique.contains(&module) {
            unique.push(module);
        }
    }
    *modules = unique;
}

/// Failure in QPA's finite-field module length enumeration operations.
#[derive(Debug, Error)]
pub enum ModuleLengthError {
    /// QPA rejects an empty list before inspecting an acting algebra.
    #[error("list entered into AllModulesOfLengthPlusOne is empty")]
    EmptyModuleList,
    /// All three source methods are finite-field-only.
    #[error("the algebra is not over a finite field")]
    FiniteFieldRequired,
    /// Every entered module must act over the identical algebra.
    #[error("the entered list does not consist of modules over the same algebra")]
    ModulesOverDifferentAlgebras,
    /// QPA calls total vector-space dimension the module length here.
    #[error("the entered list does not have the same length")]
    ModulesOfDifferentLengths,
    /// A nonnegative GAP integer bound cannot be represented by `usize`.
    #[error("module length bound is not representable by usize: {value}")]
    LengthBoundNotRepresentable {
        /// Rejected bound.
        value: Integer,
    },
    /// A nonempty two-module direct sum cannot return `fail` in this branch.
    #[error("an internal nonempty extension direct sum returned fail")]
    InternalDirectSumFailed,
    /// Both maps supplied by the source have the same syzygy source.
    #[error("an internal extension pushout returned fail")]
    InternalPushOutFailed,
    /// QPA subsequently indexes `modules[1]` after the two Ext loops.
    #[error("the extension loops constructed no modules")]
    EmptyConstructedModuleList,
    /// The helper is entered only for Ext dimension greater than one.
    #[error("an internal linear-combination basis was empty")]
    EmptyLinearCombinationBasis,
    /// Finite fields always dispatch QPA's decomposition-based method.
    #[error("IsIndecomposableModule unexpectedly invoked TryNextMethod")]
    IndecomposabilityMethodUnavailable,
    /// Supported-field enumeration failed.
    #[error(transparent)]
    Coefficient(#[from] CoefficientDomainError),
    /// Simple modules, projective covers, and Ext construction failed.
    #[error(transparent)]
    Enveloping(#[from] EnvelopingModuleError),
    /// Direct sums, pushouts, map arithmetic, or isomorphism failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// Deterministic indecomposable decomposition failed.
    #[error(transparent)]
    Decomposition(#[from] ModuleDecompositionError),
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
