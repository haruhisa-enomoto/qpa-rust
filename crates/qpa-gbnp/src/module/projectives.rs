/// Invoke QPA's immediate `GlobalDimension` methods from already-known
/// algebra properties.
///
/// Path algebras return the constructor-installed value. A quotient first
/// returns an explicitly cached attribute; otherwise known semisimplicity
/// gives zero, while known selfinjectivity together with known
/// nonsemisimplicity gives `infinity`. `None` means neither immediate method
/// is applicable yet. A successful inferred value is cached like GAP's
/// attribute call.
#[must_use]
pub fn global_dimension_from_known_properties<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Option<AlgebraGlobalDimension> {
    let QuiverAlgebra::Quotient(quotient) = algebra else {
        let QuiverAlgebra::Path(path_algebra) = algebra else {
            unreachable!("every quiver algebra is a path algebra or quotient")
        };
        return Some(AlgebraGlobalDimension::Finite(
            path_algebra.qpa_global_dimension(),
        ));
    };
    if let Some(cached) = quotient.cached_global_dimension() {
        return Some(cached);
    }
    let inferred = if quotient.cached_is_semisimple_algebra() == Some(true) {
        AlgebraGlobalDimension::Finite(0)
    } else if quotient.cached_is_selfinjective_algebra() == Some(true)
        && quotient.cached_is_semisimple_algebra() == Some(false)
    {
        AlgebraGlobalDimension::Infinite
    } else {
        return None;
    };
    Some(quotient.cache_global_dimension(inferred))
}
/// Translate QPA's `IsHereditaryAlgebra` dispatch for quiver algebras.
///
/// Every path algebra is hereditary by QPA's true-method implication. For an
/// admissible quotient, an already cached global dimension is compared with
/// two; otherwise `GlobalDimensionOfAlgebra(A, 1)` supplies the exact bounded
/// test. `None` preserves the absence of this method for a nonadmissible
/// quotient.
///
/// # Errors
///
/// Propagates the default admissibility and bounded-global-dimension paths.
pub fn is_hereditary_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<bool>, EnvelopingModuleError> {
    let QuiverAlgebra::Quotient(quotient) = algebra else {
        return Ok(Some(true));
    };
    if let Some(cached) = quotient.cached_is_hereditary_algebra() {
        return Ok(Some(cached));
    }
    if !quotient_is_admissible(quotient)? {
        return Ok(None);
    }
    let result = if let Some(cached) = quotient.cached_global_dimension() {
        matches!(cached, AlgebraGlobalDimension::Finite(dimension) if dimension < 2)
    } else {
        match global_dimension_of_algebra(algebra, 1)? {
            Some(GlobalDimensionOfAlgebraResult::Finite(dimension)) => dimension < 2,
            Some(
                GlobalDimensionOfAlgebraResult::Infinite
                | GlobalDimensionOfAlgebraResult::ExceedsBound,
            ) => false,
            None => return Ok(None),
        }
    };
    Ok(Some(quotient.cache_is_hereditary_algebra(result)))
}

/// Return QPA's true-method implications for
/// `IsFiniteGlobalDimensionAlgebra`.
///
/// Every path algebra carries the property. A quotient acquires it when
/// either semisimplicity or heredity is already known true. `None` means QPA
/// has installed no conclusion from these filters; it does not mean false.
#[must_use]
pub fn is_finite_global_dimension_from_known_properties<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Option<bool> {
    match algebra {
        QuiverAlgebra::Path(_) => Some(true),
        QuiverAlgebra::Quotient(quotient)
            if quotient.cached_is_semisimple_algebra() == Some(true)
                || quotient.cached_is_hereditary_algebra() == Some(true) =>
        {
            Some(true)
        }
        QuiverAlgebra::Quotient(_) => None,
    }
}

/// Translate QPA's enveloping-module test for `IsSymmetricAlgebra(A)`.
///
/// Infinite-dimensional input is false. A finite path algebra is symmetric
/// exactly when its quiver has no arrows. For a quotient, QPA first requires
/// selfinjectivity and then compares the regular bimodule with its dual over
/// the enveloping algebra. `None` preserves the selfinjectivity method's
/// nonadmissible quotient fall-through.
///
/// # Errors
///
/// Propagates the finite-basis, selfinjectivity, enveloping-module, dual, and
/// module-isomorphism algorithms.
pub fn is_symmetric_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<bool>, EnvelopingModuleError> {
    if let QuiverAlgebra::Quotient(quotient) = algebra
        && let Some(cached) = quotient.cached_is_symmetric_algebra()
    {
        return Ok(Some(cached));
    }
    if finite_canonical_basis(algebra)?.is_none() {
        if let QuiverAlgebra::Quotient(quotient) = algebra {
            return Ok(Some(quotient.cache_is_symmetric_algebra(false)));
        }
        return Ok(Some(false));
    }
    if let QuiverAlgebra::Path(path_algebra) = algebra {
        return Ok(Some(path_algebra.quiver().number_of_arrows() == 0));
    }
    let QuiverAlgebra::Quotient(quotient) = algebra else {
        unreachable!("every quiver algebra is a path algebra or quotient")
    };
    match is_selfinjective_algebra(algebra)? {
        Some(false) => return Ok(Some(quotient.cache_is_symmetric_algebra(false))),
        Some(true) => {}
        None => return Ok(None),
    }
    let regular = algebra_as_module_over_enveloping_algebra(algebra)?
        .ok_or(EnvelopingModuleError::InternalSymmetricRegularModuleFailed)?;
    let dual = dual_of_algebra_as_module_over_enveloping_algebra(algebra)?
        .ok_or(EnvelopingModuleError::InternalSymmetricDualModuleFailed)?;
    let result = isomorphism_of_modules(&regular, &dual)?.is_some();
    Ok(Some(quotient.cache_is_symmetric_algebra(result)))
}

/// Translate QPA's projective top/socle test for
/// `IsWeaklySymmetricAlgebra(A)`.
///
/// A selfinjective algebra is weakly symmetric precisely when every
/// indecomposable projective has the same dimension vector in its socle and
/// top. `None` preserves the selfinjectivity method's nonadmissible quotient
/// fall-through.
///
/// # Errors
///
/// Propagates selfinjectivity, projective, socle, and top computations.
pub fn is_weakly_symmetric_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<bool>, EnvelopingModuleError> {
    if let QuiverAlgebra::Quotient(quotient) = algebra
        && let Some(cached) = quotient.cached_is_weakly_symmetric_algebra()
    {
        return Ok(Some(cached));
    }
    let Some(selfinjective) = is_selfinjective_algebra(algebra)? else {
        return Ok(None);
    };
    if !selfinjective {
        if let QuiverAlgebra::Quotient(quotient) = algebra {
            return Ok(Some(quotient.cache_is_weakly_symmetric_algebra(false)));
        }
        return Ok(Some(false));
    }
    let mut result = true;
    for projective in indec_projective_modules(algebra)? {
        let socle = socle_of_module(projective.clone())?;
        let top = top_of_module(projective)?;
        if socle.dimension_vector() != top.dimension_vector() {
            result = false;
            break;
        }
    }
    if let QuiverAlgebra::Quotient(quotient) = algebra {
        result = quotient.cache_is_weakly_symmetric_algebra(result);
    }
    Ok(Some(result))
}

/// Translate QPA's `IsSchurianAlgebra(A)` dimension-vector method.
///
/// The method applies to path algebras and admissible quotients. It requires
/// finite dimension and flattens the indecomposable-projective dimension
/// vectors, accepting exactly when every entry is at most one. `None`
/// preserves `TryNextMethod()` for a nonadmissible quotient.
///
/// # Errors
///
/// Reports QPA's explicit infinite-dimensional error and propagates the
/// default admissibility, basis, and projective computations.
pub fn is_schurian_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<bool>, EnvelopingModuleError> {
    if let QuiverAlgebra::Quotient(quotient) = algebra {
        if let Some(cached) = quotient.cached_is_schurian_algebra() {
            return Ok(Some(cached));
        }
        if !quotient_is_admissible(quotient)? {
            return Ok(None);
        }
    }
    if finite_canonical_basis(algebra)?.is_none() {
        return Err(EnvelopingModuleError::InfiniteSchurianAlgebra);
    }
    let result = indec_projective_modules(algebra)?
        .iter()
        .flat_map(PathAlgebraModule::dimension_vector)
        .all(|&dimension| dimension <= 1);
    if let QuiverAlgebra::Quotient(quotient) = algebra {
        return Ok(Some(quotient.cache_is_schurian_algebra(result)));
    }
    Ok(Some(result))
}

/// Translate both QPA `IsSemicommutativeAlgebra` methods.
///
/// A path algebra must be schurian and acyclic. A quotient first rejects a
/// cyclic original quiver, then requires the specialized schurian property,
/// and finally checks for each ordered vertex pair that either every
/// nontrivial path lies in the defining ideal or none does. `None` preserves
/// the nested schurian method's nonadmissible fall-through.
///
/// # Errors
///
/// Propagates schurian, path-embedding, and exact ideal-membership failures.
pub fn is_semicommutative_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<bool>, EnvelopingModuleError> {
    if let QuiverAlgebra::Path(path_algebra) = algebra {
        if !is_schurian_algebra(algebra)?.unwrap_or(false) {
            return Ok(Some(false));
        }
        return Ok(Some(path_algebra.quiver().is_acyclic()));
    }
    let QuiverAlgebra::Quotient(quotient) = algebra else {
        unreachable!("every quiver algebra is a path algebra or quotient")
    };
    if let Some(cached) = quotient.cached_is_semicommutative_algebra() {
        return Ok(Some(cached));
    }
    let original = quotient.original_path_algebra();
    if !original.quiver().is_acyclic() {
        return Ok(Some(quotient.cache_is_semicommutative_algebra(false)));
    }
    match is_schurian_algebra(algebra)? {
        Some(false) => return Ok(Some(quotient.cache_is_semicommutative_algebra(false))),
        Some(true) => {}
        None => return Ok(None),
    }

    let vertex_count = original.quiver().number_of_vertices();
    let pair_count = vertex_count
        .checked_mul(vertex_count)
        .ok_or(EnvelopingModuleError::SemicommutativePairCountOverflow { vertex_count })?;
    let mut membership = vec![(false, false); pair_count];
    for path in original
        .quiver()
        .paths()
        .filter(|path| path.length().is_some_and(|length| length > 0))
    {
        let (PathEndpoint::Vertex(source), PathEndpoint::Vertex(target)) =
            (path.source(), path.target())
        else {
            unreachable!("a nontrivial path has vertex endpoints")
        };
        let position = (source.index().get() - 1) * vertex_count + target.index().get() - 1;
        let lies_in_ideal = quotient.ideal().contains(original.element(path)?)?;
        if lies_in_ideal {
            membership[position].0 = true;
        } else {
            membership[position].1 = true;
        }
        if membership[position].0 && membership[position].1 {
            return Ok(Some(quotient.cache_is_semicommutative_algebra(false)));
        }
    }
    Ok(Some(quotient.cache_is_semicommutative_algebra(true)))
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
    hom_from_projective_with_family(element, target, &mut None)
}

fn hom_from_projective_with_family<F: Field>(
    element: &PathAlgebraModuleElement<F>,
    target: &PathAlgebraModule<F>,
    projectives: &mut Option<Vec<PathAlgebraModule<F>>>,
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
    let projectives = match projectives {
        Some(family) => family,
        None => projectives.insert(indec_projective_modules(target.right_acting_algebra())?),
    };
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
