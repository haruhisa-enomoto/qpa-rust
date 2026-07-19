/// Translate QPA's generic `Hom(R, M, N)` method for path-ring matrix modules.
///
/// This is the `HOPF_MatrixModuleHom` implementation from `module.gi`, not
/// the separately installed `HomOverAlgebra` algorithm. For submodules QPA
/// computes their ordered canonical bases, converts them with
/// `SubmoduleAsModule`, and then reuses the standalone full matrices without
/// changing their row or column order.
///
/// # Errors
///
/// Preserves QPA's pure-path-ring dispatch, identical acting-algebra check,
/// zero-module and degenerate-nullspace failures, and every submodule
/// conversion failure. Exact matrix failures and address-size overflow also
/// propagate.
pub fn generic_hom_over_path_ring<F: Field>(
    algebra: &QuiverAlgebra<F>,
    source: PathAlgebraMatModuleObject<F>,
    range: PathAlgebraMatModuleObject<F>,
) -> Result<GenericPathAlgebraModuleHomSpace<F>, PathAlgebraModuleError> {
    if !algebra.is_path_algebra() {
        return Err(PathAlgebraModuleError::GenericHomPathRingDispatchFailure);
    }
    if !algebra.is_identical(source.right_acting_algebra())
        || !algebra.is_identical(range.right_acting_algebra())
    {
        return Err(PathAlgebraModuleError::GenericHomActingAlgebraMismatch);
    }

    let source_basis = source.canonical_basis_elements()?;
    let range_basis = range.canonical_basis_elements()?;
    let source_standalone = match &source {
        PathAlgebraMatModuleObject::Module(module) => module.clone(),
        PathAlgebraMatModuleObject::Submodule(submodule) => submodule_as_module(submodule)?,
    };
    let range_standalone = match &range {
        PathAlgebraMatModuleObject::Module(module) => module.clone(),
        PathAlgebraMatModuleObject::Submodule(submodule) => submodule_as_module(submodule)?,
    };
    let matrices = hopf_matrix_module_hom(algebra, &source_standalone, &range_standalone)?;
    let flattened_width = source_basis
        .len()
        .checked_mul(range_basis.len())
        .ok_or(PathAlgebraModuleError::GenericHomEquationSizeOverflow)?;
    let flattened_generators = DenseMatrix::from_flat(
        matrices.len(),
        flattened_width,
        matrices.iter().flat_map(flatten_matrix).collect(),
        algebra.field().clone(),
    )?;
    let generators = matrices
        .into_iter()
        .map(|matrix| {
            GenericPathAlgebraModuleHomomorphism(Arc::new(
                GenericPathAlgebraModuleHomomorphismData {
                    source: source.clone(),
                    range: range.clone(),
                    source_basis: source_basis.clone().into_boxed_slice(),
                    range_basis: range_basis.clone().into_boxed_slice(),
                    matrix,
                },
            ))
        })
        .collect();
    Ok(GenericPathAlgebraModuleHomSpace {
        field: algebra.field().clone(),
        source,
        range,
        generators,
        flattened_generators,
    })
}
/// Translate QPA's generic `End(R, M)` method for a path-ring matrix module.
///
/// The returned basis is exactly the generic Hom basis supplied to
/// `AlgebraWithOne(..., "basis")`; unlike [`end_over_algebra`], QPA performs
/// no generated-subalgebra semi-echelon conversion here.
///
/// # Errors
///
/// Propagates [`generic_hom_over_path_ring`] and exact identity construction.
pub fn generic_end_over_path_ring<F: Field>(
    algebra: &QuiverAlgebra<F>,
    module: PathAlgebraMatModuleObject<F>,
) -> Result<GenericPathAlgebraModuleEndomorphismAlgebra<F>, PathAlgebraModuleError> {
    let hom = generic_hom_over_path_ring(algebra, module.clone(), module.clone())?;
    let basis = module.canonical_basis_elements()?;
    let dimension = basis.len();
    let one =
        GenericPathAlgebraModuleHomomorphism(Arc::new(GenericPathAlgebraModuleHomomorphismData {
            source: module.clone(),
            range: module.clone(),
            source_basis: basis.clone().into_boxed_slice(),
            range_basis: basis.into_boxed_slice(),
            matrix: DenseMatrix::identity(dimension, algebra.field().clone()),
        }));
    Ok(GenericPathAlgebraModuleEndomorphismAlgebra {
        field: hom.field,
        module,
        generators: hom.generators,
        flattened_generators: hom.flattened_generators,
        one,
    })
}

/// Translate QPA's `RightAlgebraModuleToPathAlgebraMatModule(M)` method.
///
/// QPA projects each entered basis vector through every vertex idempotent,
/// takes `Basis(Subspace(...))` in insertion-preserving semi-echelon order,
/// acts those vertex-basis rows by each arrow, and solves coordinates in the
/// target vertex basis. This function follows that loop literally on the
/// generator-action matrices of [`BasedRightAlgebraModule`].
///
/// # Errors
///
/// Reports an action that leaves the computed target vertex space, a converted
/// representation that violates the entered quotient relations, and exact
/// matrix failures.
pub fn right_algebra_module_to_path_algebra_mat_module<F: Field>(
    module: &BasedRightAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    let algebra = module.algebra();
    let field = algebra.field().clone();
    let mut vertex_bases = Vec::with_capacity(algebra.quiver().number_of_vertices());
    for action in module.vertex_actions() {
        let projected = (0..module.dimension())
            .filter(|&row| action.row(row).iter().any(|entry| !field.is_zero(entry)))
            .map(|row| action.row(row).to_vec())
            .collect::<Vec<_>>();
        let basis = if projected.is_empty() {
            DenseMatrix::from_flat(0, module.dimension(), Vec::new(), field.clone())?
        } else {
            DenseMatrix::from_rows(projected, field.clone())?
                .semi_echelon_transformation()
                .vectors
        };
        vertex_bases.push(basis);
    }
    let dimensions = vertex_bases
        .iter()
        .map(DenseMatrix::nrows)
        .collect::<Vec<_>>();
    let mut arrow_matrices = Vec::new();
    for arrow in algebra.quiver().arrows() {
        let source = arrow.source().index().zero_based();
        let target = arrow.target().index().zero_based();
        if dimensions[source] == 0 || dimensions[target] == 0 {
            continue;
        }
        let acted =
            vertex_bases[source].multiply(&module.arrow_actions()[arrow.index().zero_based()])?;
        let mut entries = Vec::with_capacity(dimensions[source] * dimensions[target]);
        for row in 0..acted.nrows() {
            entries.extend(vertex_bases[target].solution_row(acted.row(row))?.ok_or(
                PathAlgebraModuleError::BasedModuleArrowImageOutsideTarget {
                    arrow: arrow.name().into(),
                },
            )?);
        }
        arrow_matrices.push(NamedArrowMatrix::new(
            arrow.name(),
            DenseMatrix::from_flat(
                dimensions[source],
                dimensions[target],
                entries,
                field.clone(),
            )?,
        ));
    }
    right_module_over_path_algebra(algebra.clone(), dimensions, &arrow_matrices)?
        .ok_or(PathAlgebraModuleError::GeneratedBasedModuleConversionRejectedRelations)
}

/// Translate QPA's matrix-module `SubAlgebraModule(V, generators)` method.
///
/// Entered elements are split by vertex idempotents in generator-major,
/// vertex-major order, and zero components are removed before storage.
///
/// An empty entered list follows GAP's higher-ranking generic dispatch and
/// returns the zero submodule. Nonempty lists use QPA's installed method.
///
/// # Errors
///
/// Rejects elements outside the identical ambient module family and preserves
/// QPA's quotient-module action failure.
pub fn subalgebra_module<F: Field>(
    parent: PathAlgebraModule<F>,
    generators: &[PathAlgebraModuleElement<F>],
) -> Result<PathAlgebraModuleSubmodule<F>, PathAlgebraModuleError> {
    if generators.is_empty() {
        return Ok(PathAlgebraModuleSubmodule(Arc::new(
            PathAlgebraModuleSubmoduleData {
                parent,
                generators: Box::new([]),
            },
        )));
    }
    if !parent.right_acting_algebra().is_path_algebra() {
        return Err(PathAlgebraModuleError::QuotientSubalgebraModuleFailure);
    }
    if generators
        .iter()
        .any(|generator| !parent.is_identical(generator.module()))
    {
        return Err(PathAlgebraModuleError::ForeignSubalgebraModuleGenerator);
    }
    let mut uniform = Vec::new();
    for generator in generators {
        for vertex in 0..parent.dimension_vector().len() {
            let mut component = parent.zero_element();
            component.coordinates[vertex].clone_from(&generator.coordinates[vertex]);
            if !component.is_zero() {
                uniform.push(component);
            }
        }
    }
    Ok(PathAlgebraModuleSubmodule(Arc::new(
        PathAlgebraModuleSubmoduleData {
            parent,
            generators: uniform.into_boxed_slice(),
        },
    )))
}

fn canonical_basis_of_subalgebra_module<F: Field>(
    submodule: &PathAlgebraModuleSubmodule<F>,
) -> Result<Vec<PathAlgebraModuleElement<F>>, PathAlgebraModuleError> {
    if submodule.generators().is_empty() {
        return Ok(Vec::new());
    }
    let parent = submodule.parent();
    let field = parent.right_acting_algebra().field().clone();
    let mut basis = submodule.generators().to_vec();
    let entered_rows = DenseMatrix::from_rows(
        basis
            .iter()
            .map(PathAlgebraModuleElement::basis_coefficients)
            .collect(),
        field.clone(),
    )?;
    if entered_rows.canonical_row_space_basis().nrows() != basis.len() {
        return Err(PathAlgebraModuleError::DependentSubalgebraModuleGenerators);
    }

    let arrow_elements = parent
        .right_acting_algebra()
        .quiver()
        .arrows()
        .map(|arrow| quiver_algebra_path_element(parent.right_acting_algebra(), arrow.path()))
        .collect::<Result<Vec<_>, _>>()?;
    let mut position = 0;
    while position < basis.len() {
        let element = basis[position].clone();
        for arrow in &arrow_elements {
            let image = element.act(arrow)?;
            if image.is_zero() {
                continue;
            }
            let rows = DenseMatrix::from_rows(
                basis
                    .iter()
                    .map(PathAlgebraModuleElement::basis_coefficients)
                    .collect(),
                field.clone(),
            )?;
            if rows.solution_row(&image.basis_coefficients())?.is_none() {
                basis.push(image);
            }
        }
        position += 1;
    }
    Ok(basis)
}

fn canonical_vertex_bases_of_submodule<F: Field>(
    submodule: &PathAlgebraModuleSubmodule<F>,
) -> Result<Vec<DenseMatrix<F>>, PathAlgebraModuleError> {
    let parent = submodule.parent();
    let field = parent.right_acting_algebra().field().clone();
    let canonical = canonical_basis_of_subalgebra_module(submodule)?;
    let mut rows = vec![Vec::new(); parent.dimension_vector().len()];
    for element in canonical {
        for (vertex, &dimension) in parent.dimension_vector().iter().enumerate() {
            if element.coordinates[vertex][..dimension]
                .iter()
                .any(|entry| !field.is_zero(entry))
            {
                rows[vertex].push(element.coordinates[vertex][..dimension].to_vec());
            }
        }
    }
    rows.into_iter()
        .zip(parent.dimension_vector())
        .map(|(rows, &dimension)| {
            if rows.is_empty() {
                DenseMatrix::from_flat(0, dimension, Vec::new(), field.clone())
            } else {
                Ok(DenseMatrix::from_rows(rows, field.clone())?.canonical_row_space_basis())
            }
        })
        .collect::<Result<Vec<_>, MatrixError>>()
        .map_err(Into::into)
}

/// Translate QPA's `SubmoduleAsModule(N)` operation.
///
/// The returned representation acts over QPA's stored original path algebra,
/// including when the ambient module acts over a quotient.
///
/// # Errors
///
/// Preserves QPA's canonical-basis and two-argument-constructor failures and
/// propagates exact coordinate solving.
pub fn submodule_as_module<F: Field>(
    submodule: &PathAlgebraModuleSubmodule<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    if submodule.generators().is_empty() {
        return Err(PathAlgebraModuleError::ZeroSubmoduleAsModuleFailure);
    }
    let parent = submodule.parent();
    let field = parent.right_acting_algebra().field().clone();
    let bases = canonical_vertex_bases_of_submodule(submodule)?;
    let dimensions = bases.iter().map(DenseMatrix::nrows).collect::<Vec<_>>();
    let original = QuiverAlgebra::Path(
        parent
            .right_acting_algebra()
            .original_path_algebra()
            .clone(),
    );
    if original.quiver().number_of_arrows() == 0 {
        return Err(PathAlgebraModuleError::ArrowlessSubmoduleAsModuleFailure);
    }
    let mut matrices = Vec::with_capacity(original.quiver().number_of_arrows());
    for arrow in original.quiver().arrows() {
        let source = arrow.source().index().zero_based();
        let target = arrow.target().index().zero_based();
        if dimensions[source] == 0 || dimensions[target] == 0 {
            return Err(PathAlgebraModuleError::ZeroSpaceSubmoduleAsModuleFailure {
                arrow: arrow.name().into(),
            });
        }
        let ambient_action = logical_submatrix(
            &parent.arrow_matrices()[arrow.index().zero_based()],
            parent.dimension_vector()[source],
            parent.dimension_vector()[target],
        )?;
        let acted = bases[source].multiply(&ambient_action)?;
        let mut entries = Vec::with_capacity(dimensions[source] * dimensions[target]);
        for row in 0..acted.nrows() {
            entries.extend(
                bases[target]
                    .solution_row(acted.row(row))?
                    .ok_or(PathAlgebraModuleError::InternalSubmoduleCoordinates)?,
            );
        }
        matrices.push(DenseMatrix::from_flat(
            dimensions[source],
            dimensions[target],
            entries,
            field.clone(),
        )?);
    }
    right_module_over_path_algebra_infer_dimension(
        original,
        &PathAlgebraModuleArrowActions::Positional(matrices),
    )?
    .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)
}

/// Translate QPA's `NaturalHomomorphismBySubAlgebraModule(M, N)` method.
///
/// Despite its name and comment, QPA returns a basis-image mapping from `M`
/// to `M/N`. Its mutable `base` list is not restored between basis images, so
/// coordinates from earlier vertices leak into later images. This function
/// intentionally preserves that observable mapping.
///
/// # Errors
///
/// Requires `N` to have the identical entered parent and preserves QPA's
/// canonical-basis, zero-space positional-constructor, and zero-quotient
/// failures.
#[allow(clippy::too_many_lines)]
pub fn natural_homomorphism_by_subalgebra_module<F: Field>(
    module: &PathAlgebraModule<F>,
    submodule: &PathAlgebraModuleSubmodule<F>,
) -> Result<PathAlgebraModuleBasisMapping<F>, PathAlgebraModuleError> {
    if !module.is_identical(submodule.parent()) {
        return Err(PathAlgebraModuleError::SubalgebraModuleNotInEnteredParent);
    }
    if submodule.generators().is_empty() {
        return Err(PathAlgebraModuleError::ZeroNaturalSubmoduleFailure);
    }
    let subspace_bases = canonical_vertex_bases_of_submodule(submodule)?;
    let quotients = subspace_bases
        .iter()
        .map(natural_quotient_by_row_space)
        .collect::<Result<Vec<_>, _>>()?;
    let quotient_dimensions = quotients
        .iter()
        .map(|quotient| quotient.lifts.nrows())
        .collect::<Vec<_>>();
    let original = QuiverAlgebra::Path(
        module
            .right_acting_algebra()
            .original_path_algebra()
            .clone(),
    );
    if original.quiver().number_of_arrows() == 0 {
        return Err(PathAlgebraModuleError::ArrowlessNaturalSubmoduleFailure);
    }
    let mut matrices = Vec::with_capacity(original.quiver().number_of_arrows());
    for arrow in original.quiver().arrows() {
        let source = arrow.source().index().zero_based();
        let target = arrow.target().index().zero_based();
        if quotient_dimensions[source] == 0 || quotient_dimensions[target] == 0 {
            return Err(PathAlgebraModuleError::ZeroSpaceNaturalSubmoduleFailure {
                arrow: arrow.name().into(),
            });
        }
        let action = logical_submatrix(
            &module.arrow_matrices()[arrow.index().zero_based()],
            module.dimension_vector()[source],
            module.dimension_vector()[target],
        )?;
        matrices.push(
            quotients[source]
                .lifts
                .multiply(&action)?
                .multiply(&quotients[target].projection)?,
        );
    }
    let quotient = right_module_over_path_algebra_infer_dimension(
        original,
        &PathAlgebraModuleArrowActions::Positional(matrices),
    )?
    .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)?;
    if quotient.dimension() == 0 {
        return Err(PathAlgebraModuleError::ZeroNaturalSubmoduleQuotientFailure);
    }

    let mut base = quotient.zero_element();
    let mut images = Vec::with_capacity(module.dimension());
    for (vertex, &dimension) in module.dimension_vector().iter().enumerate() {
        for row in 0..dimension {
            if quotient_dimensions[vertex] != 0 {
                base.coordinates[vertex] = quotients[vertex].projection.row(row).to_vec();
            }
            images.push(base.clone());
        }
    }
    Ok(PathAlgebraModuleBasisMapping(Arc::new(
        PathAlgebraModuleBasisMappingData {
            source: module.clone(),
            range: quotient,
            basis_images: images.into_boxed_slice(),
            surjective: true,
        },
    )))
}

/// Translate QPA's `SubRepresentationInclusion(M, generators)` method.
///
/// QPA first makes every generator vertex-uniform, repeatedly acts by every
/// arrow until no nonzero elements remain, takes GAP's fully reduced canonical
/// row-space basis, and then reconstructs both the subrepresentation and its
/// inclusion.  This function follows that sequence and ordering.
///
/// # Errors
///
/// Rejects generators outside the identical module family and propagates
/// algebra, matrix, and checked module-construction failures.
#[allow(clippy::too_many_lines)]
pub fn subrepresentation_inclusion<F: Field>(
    module: PathAlgebraModule<F>,
    generators: &[PathAlgebraModuleElement<F>],
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if generators
        .iter()
        .any(|generator| !module.is_identical(generator.module()))
    {
        return Err(PathAlgebraModuleError::ForeignSubrepresentationGenerator);
    }

    let algebra = module.right_acting_algebra().clone();
    let vertex_count = module.dimension_vector().len();
    if generators.is_empty() {
        let zero = right_module_over_path_algebra(algebra, vec![0; vertex_count], &[])?
            .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)?;
        return zero_module_mapping(zero, module);
    }

    let mut uniform_generators = Vec::new();
    for generator in generators {
        for vertex in 0..vertex_count {
            let mut uniform = module.zero_element();
            uniform.coordinates[vertex].clone_from(&generator.coordinates[vertex]);
            if !uniform.is_zero() {
                uniform_generators.push(uniform);
            }
        }
    }

    let arrow_elements = algebra
        .quiver()
        .arrows()
        .map(|arrow| quiver_algebra_path_element(&algebra, arrow.path()))
        .collect::<Result<Vec<_>, _>>()?;
    let mut spanning_elements = Vec::new();
    let mut new = uniform_generators;
    while !new.is_empty() {
        let mut next = Vec::new();
        for element in &new {
            for arrow in &arrow_elements {
                let image = element.act(arrow)?;
                if !image.is_zero() {
                    next.push(image);
                }
            }
        }
        spanning_elements.extend(new);
        new = next;
    }

    let field = algebra.field().clone();
    let spanning_rows = DenseMatrix::from_flat(
        spanning_elements.len(),
        module.dimension(),
        spanning_elements
            .iter()
            .flat_map(PathAlgebraModuleElement::basis_coefficients)
            .collect(),
        field.clone(),
    )?;
    let canonical_basis = spanning_rows.canonical_row_space_basis();
    let basis_elements = (0..canonical_basis.nrows())
        .map(|row| {
            PathAlgebraModuleElement::from_logical_flat(module.clone(), canonical_basis.row(row))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut sub_dimension = vec![0; vertex_count];
    for element in &basis_elements {
        for vertex in element.support_vertices() {
            sub_dimension[vertex - 1] += 1;
        }
    }
    if canonical_basis.nrows() != sub_dimension.iter().sum::<usize>() {
        return Err(PathAlgebraModuleError::NonuniformSubrepresentationBasis);
    }
    let sub_offsets = dimension_offsets(&sub_dimension);
    let module_offsets = dimension_offsets(module.dimension_vector());

    let mut arrow_matrices = Vec::new();
    for (arrow, arrow_element) in algebra.quiver().arrows().zip(&arrow_elements) {
        let source = arrow.source().index().zero_based();
        let target = arrow.target().index().zero_based();
        if sub_dimension[source] == 0 || sub_dimension[target] == 0 {
            continue;
        }
        let mut rows = Vec::with_capacity(sub_dimension[source] * sub_dimension[target]);
        for basis_element in basis_elements
            .iter()
            .skip(sub_offsets[source])
            .take(sub_dimension[source])
        {
            let image = basis_element.act(arrow_element)?;
            let coefficients = canonical_basis
                .solution_row(&image.basis_coefficients())?
                .ok_or(PathAlgebraModuleError::InternalSubrepresentationCoordinates)?;
            rows.extend_from_slice(
                &coefficients[sub_offsets[target]..sub_offsets[target] + sub_dimension[target]],
            );
        }
        arrow_matrices.push(NamedArrowMatrix::new(
            arrow.name(),
            DenseMatrix::from_flat(
                sub_dimension[source],
                sub_dimension[target],
                rows,
                field.clone(),
            )?,
        ));
    }

    let submodule =
        right_module_over_path_algebra(algebra, sub_dimension.clone(), &arrow_matrices)?
            .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)?;

    let inclusion_matrices = (0..vertex_count)
        .map(|vertex| -> Result<DenseMatrix<F>, MatrixError> {
            if sub_dimension[vertex] == 0 || module.dimension_vector()[vertex] == 0 {
                return Ok(DenseMatrix::zero(
                    physical_dimension(sub_dimension[vertex]),
                    physical_dimension(module.dimension_vector()[vertex]),
                    field.clone(),
                ));
            }
            let mut entries =
                Vec::with_capacity(sub_dimension[vertex] * module.dimension_vector()[vertex]);
            for row in sub_offsets[vertex]..sub_offsets[vertex] + sub_dimension[vertex] {
                entries.extend_from_slice(
                    &canonical_basis.row(row)[module_offsets[vertex]
                        ..module_offsets[vertex] + module.dimension_vector()[vertex]],
                );
            }
            DenseMatrix::from_flat(
                sub_dimension[vertex],
                module.dimension_vector()[vertex],
                entries,
                field.clone(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    right_module_hom_over_algebra(submodule, module, inclusion_matrices)
}

/// Translate QPA's `SubRepresentation(M, generators)` convenience method.
///
/// # Errors
///
/// Propagates [`subrepresentation_inclusion`].
pub fn subrepresentation<F: Field>(
    module: PathAlgebraModule<F>,
    generators: &[PathAlgebraModuleElement<F>],
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(subrepresentation_inclusion(module, generators)?
        .source()
        .clone())
}

/// Translate QPA's `RadicalOfModuleInclusion` method for finite path algebras
/// and finite admissible quotients.
///
/// # Errors
///
/// Reports the same source-domain restriction and propagates generated
/// subrepresentation failures.
pub fn radical_of_module_inclusion<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let algebra = module.right_acting_algebra();
    let supported = match algebra {
        QuiverAlgebra::Path(algebra) => algebra.has_finite_canonical_basis(),
        QuiverAlgebra::Quotient(algebra) => {
            algebra.is_finite_dimensional()? && algebra.is_admissible()?
        }
    };
    if !supported {
        return Err(PathAlgebraModuleError::UnsupportedRadicalAlgebra);
    }

    let basis = module.basis_elements();
    let arrows = algebra
        .quiver()
        .arrows()
        .map(|arrow| quiver_algebra_path_element(algebra, arrow.path()))
        .collect::<Result<Vec<_>, _>>()?;
    let mut generators = Vec::new();
    for arrow in &arrows {
        for element in &basis {
            let image = element.act(arrow)?;
            if !image.is_zero() && !generators.contains(&image) {
                generators.push(image);
            }
        }
    }
    subrepresentation_inclusion(module, &generators)
}

/// Translate QPA's `RadicalOfModule(M)` convenience method.
///
/// # Errors
///
/// Propagates [`radical_of_module_inclusion`].
pub fn radical_of_module<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(radical_of_module_inclusion(module)?.source().clone())
}

/// Translate QPA's componentwise `DimensionVectorPartialOrder(M, N)`.
///
/// GAP additive-list arithmetic pads the shorter dimension vector with
/// zeros, so modules over quivers with different vertex counts remain in the
/// method domain.
#[must_use]
pub fn dimension_vector_partial_order<F: Field>(
    first: &PathAlgebraModule<F>,
    second: &PathAlgebraModule<F>,
) -> bool {
    let length = first
        .dimension_vector()
        .len()
        .max(second.dimension_vector().len());
    (0..length).all(|index| {
        first.dimension_vector().get(index).copied().unwrap_or(0)
            <= second.dimension_vector().get(index).copied().unwrap_or(0)
    })
}

/// Translate QPA's `RadicalSeries(M)` layer loop.
///
/// The zero branch deliberately returns [`ModuleDimensionSeries::ZeroDimensionVector`],
/// matching QPA's dimension vector rather than a list containing no layers.
///
/// # Errors
///
/// Propagates top and radical construction at every layer.
pub fn radical_series<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<ModuleDimensionSeries, PathAlgebraModuleError> {
    if module.dimension() == 0 {
        return Ok(ModuleDimensionSeries::ZeroDimensionVector(
            module.dimension_vector().to_vec(),
        ));
    }
    let mut current = module.clone();
    let mut layers = Vec::new();
    loop {
        layers.push(top_of_module(current.clone())?.dimension_vector().to_vec());
        current = radical_of_module(current)?;
        if current.dimension() == 0 {
            break;
        }
    }
    Ok(ModuleDimensionSeries::Layers(layers))
}

/// Translate QPA's cached `IsSimpleQPAModule(M)` dimension test.
///
/// QPA defines this property to be true exactly in total dimension one. A
/// true result also installs its declared implications to
/// `IsIndecomposableModule` and `IsSemisimpleModule`.
#[must_use]
pub fn is_simple_qpa_module<F: Field>(module: &PathAlgebraModule<F>) -> bool {
    if let Some(cached) = module.cached_is_simple_qpa_module() {
        return cached;
    }
    module.cache_is_simple_qpa_module(module.dimension() == 1)
}

/// Translate QPA's cached `IsSemisimpleModule(M)` radical test.
///
/// # Errors
///
/// Propagates the exact radical-module construction.
pub fn is_semisimple_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<bool, PathAlgebraModuleError> {
    if let Some(cached) = module.cached_is_semisimple_module() {
        return Ok(cached);
    }
    let radical = radical_of_module(module.clone())?;
    Ok(module.cache_is_semisimple_module(radical.dimension() == 0))
}

/// Translate QPA's module `LoewyLength` method by repeated radicals.
///
/// # Errors
///
/// Propagates the source-domain and generated-radical failures of
/// [`radical_of_module`].
pub fn loewy_length_of_module<F: Field>(
    mut module: PathAlgebraModule<F>,
) -> Result<usize, PathAlgebraModuleError> {
    let mut length = 0;
    while module.dimension() != 0 {
        module = radical_of_module(module)?;
        length += 1;
    }
    Ok(length)
}

/// Translate QPA's `CoKernelProjection` for a path-algebra module map.
///
/// At each vertex this uses GAP's `NaturalHomomorphismBySubspace`: a
/// semi-echelon image basis is greedily extended by ambient standard vectors,
/// and those selected standard vectors are retained as canonical quotient
/// basis lifts.
///
/// # Errors
///
/// Propagates exact quotient-space, module-construction, and commutativity
/// failures.
pub fn cokernel_projection<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if let Some(cached) = homomorphism.cached_cokernel_projection() {
        return Ok(cached);
    }
    let range = homomorphism.range();
    let algebra = range.right_acting_algebra().clone();
    let field = algebra.field().clone();
    let mut quotients = Vec::with_capacity(range.dimension_vector().len());
    for (vertex, &range_dimension) in range.dimension_vector().iter().enumerate() {
        let source_dimension = homomorphism.source().dimension_vector()[vertex];
        let mut image_entries = Vec::with_capacity(source_dimension * range_dimension);
        for row in 0..source_dimension {
            image_entries.extend_from_slice(
                &homomorphism.vertex_matrices()[vertex].row(row)[..range_dimension],
            );
        }
        let image_rows = DenseMatrix::from_flat(
            source_dimension,
            range_dimension,
            image_entries,
            field.clone(),
        )?;
        quotients.push(natural_quotient_by_row_space(&image_rows)?);
    }
    let quotient_dimension = quotients
        .iter()
        .map(|quotient| quotient.lifts.nrows())
        .collect::<Vec<_>>();

    let mut arrow_matrices = Vec::new();
    for arrow in algebra.quiver().arrows() {
        let source = arrow.source().index().zero_based();
        let target = arrow.target().index().zero_based();
        if quotient_dimension[source] == 0 || quotient_dimension[target] == 0 {
            continue;
        }
        let range_action = logical_submatrix(
            &range.arrow_matrices()[arrow.index().zero_based()],
            range.dimension_vector()[source],
            range.dimension_vector()[target],
        )?;
        let action = quotients[source]
            .lifts
            .multiply(&range_action)?
            .multiply(&quotients[target].projection)?;
        arrow_matrices.push(NamedArrowMatrix::new(arrow.name(), action));
    }

    let cokernel =
        right_module_over_path_algebra(algebra, quotient_dimension.clone(), &arrow_matrices)?
            .ok_or(PathAlgebraModuleError::GeneratedModuleRejectedRelations)?;
    let projection_matrices = range
        .dimension_vector()
        .iter()
        .enumerate()
        .map(|(vertex, &range_dimension)| {
            let quotient_dimension = quotient_dimension[vertex];
            if range_dimension == 0 || quotient_dimension == 0 {
                DenseMatrix::zero(
                    physical_dimension(range_dimension),
                    physical_dimension(quotient_dimension),
                    field.clone(),
                )
            } else {
                quotients[vertex].projection.clone()
            }
        })
        .collect();
    let projection = right_module_hom_over_algebra(range.clone(), cokernel, projection_matrices)?;
    projection.set_cokernel_of_what(homomorphism);
    projection.set_is_surjective(true);
    Ok(homomorphism.cache_cokernel_projection(&projection))
}

/// Translate QPA's cokernel-object convenience method as the range of its
/// canonical projection.
///
/// # Errors
///
/// Propagates [`cokernel_projection`].
pub fn cokernel_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(cokernel_projection(homomorphism)?.range().clone())
}

/// Translate QPA's `TopOfModuleProjection` attribute computation.
///
/// The total-zero branch returns the zero endomorphism of the identical input
/// module, exactly as the source method does.
///
/// # Errors
///
/// Propagates radical and cokernel construction failures.
#[allow(clippy::needless_pass_by_value)]
pub fn top_of_module_projection<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if let Some(cached) = module.cached_top_projection() {
        return Ok(cached);
    }
    let projection = if module.dimension() == 0 {
        zero_module_mapping(module.clone(), module.clone())?
    } else {
        let generated = cokernel_projection(&radical_of_module_inclusion(module.clone())?)?;
        let top = module.cache_top_module(generated.range());
        if top.is_identical(generated.range()) {
            generated
        } else {
            let projection = right_module_hom_over_algebra(
                module.clone(),
                top,
                generated.vertex_matrices().to_vec(),
            )?;
            projection.set_is_surjective(true);
            projection
        }
    };
    Ok(module.cache_top_projection(&projection))
}

/// Translate QPA's `TopOfModule(M)` convenience method.
///
/// # Errors
///
/// Propagates [`top_of_module_projection`].
pub fn top_of_module<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    if module.dimension() == 0 {
        return Ok(module);
    }
    if let Some(cached) = module.cached_top_module() {
        return Ok(cached);
    }
    Ok(top_of_module_projection(module)?.range().clone())
}

/// Translate QPA's `TopOfModule(f)` overload for a module homomorphism.
///
/// QPA lifts each vertex-major basis vector of the source top through its
/// canonical top projection, applies `f`, and projects to the range top. This
/// preserves the source method's chosen quotient representatives rather than
/// constructing an abstract induced linear map.
///
/// # Errors
///
/// Propagates top construction, preimage, image, and checked-map failures.
pub fn top_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let source_projection = top_of_module_projection(homomorphism.source().clone())?;
    let range_projection = top_of_module_projection(homomorphism.range().clone())?;
    let source_top = source_projection.range();
    let range_top = range_projection.range();
    let field = source_top.right_acting_algebra().field().clone();
    let basis = source_top.basis_elements();
    let offsets = dimension_offsets(source_top.dimension_vector());
    let matrices = source_top
        .dimension_vector()
        .iter()
        .enumerate()
        .map(
            |(vertex, &source_dimension)| -> Result<DenseMatrix<F>, PathAlgebraModuleError> {
                let range_dimension = range_top.dimension_vector()[vertex];
                if source_dimension == 0 || range_dimension == 0 {
                    return Ok(DenseMatrix::zero(
                        physical_dimension(source_dimension),
                        physical_dimension(range_dimension),
                        field.clone(),
                    ));
                }
                let mut entries = Vec::with_capacity(source_dimension * range_dimension);
                for basis_element in basis.iter().skip(offsets[vertex]).take(source_dimension) {
                    let lift = source_projection
                        .preimage_representative(basis_element)?
                        .ok_or(PathAlgebraModuleError::MissingTopBasisPreimage)?;
                    let image = homomorphism.image(&lift)?;
                    let top_image = range_projection.image(&image)?;
                    entries.extend_from_slice(
                        &top_image.external_representation()[vertex][..range_dimension],
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
    right_module_hom_over_algebra(source_top.clone(), range_top.clone(), matrices)
}

/// Translate QPA's cached `MinimalGeneratingSetOfModule(M)` construction.
///
/// QPA takes the vertex-major basis of the canonical top and lifts every
/// basis element through `TopOfModuleProjection` using its chosen natural
/// quotient representatives.
///
/// # Errors
///
/// Propagates top construction and exact preimage solving. A missing lift is
/// reported as an internal source-translation failure.
#[allow(clippy::needless_pass_by_value)]
pub fn minimal_generating_set_of_module<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleElement<F>>, PathAlgebraModuleError> {
    if let Some(cached) = module.cached_minimal_generating_set()? {
        return Ok(cached);
    }
    let projection = top_of_module_projection(module.clone())?;
    let generators = projection
        .range()
        .basis_elements()
        .iter()
        .map(|basis_element| {
            projection
                .preimage_representative(basis_element)?
                .ok_or(PathAlgebraModuleError::MissingTopBasisPreimage)
        })
        .collect::<Result<Vec<_>, _>>()?;
    module.cache_minimal_generating_set(&generators)
}

/// Translate QPA's arrow-kernel `SocleOfModuleInclusion(M)` method.
///
/// At each vertex with outgoing arrows, the source method concatenates their
/// action rows and takes GAP's left nullspace. Sink spaces are retained in
/// full, and zero spaces preserve QPA's physical `1 x 1` zero matrix.
///
/// # Errors
///
/// Propagates exact matrix arithmetic and checked module/map construction.
pub fn socle_of_module_inclusion<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if module.dimension() == 0 {
        return zero_module_mapping(module.clone(), module);
    }

    let algebra = module.right_acting_algebra();
    let field = algebra.field().clone();
    let arrows = algebra.quiver().arrows().collect::<Vec<_>>();
    let mut socle_dimensions = Vec::with_capacity(module.dimension_vector().len());
    let mut inclusion_matrices = Vec::with_capacity(module.dimension_vector().len());

    for (vertex, &dimension) in module.dimension_vector().iter().enumerate() {
        let outgoing = arrows
            .iter()
            .filter(|arrow| arrow.source().index().zero_based() == vertex)
            .collect::<Vec<_>>();
        if dimension == 0 {
            socle_dimensions.push(0);
            inclusion_matrices.push(DenseMatrix::zero(1, 1, field.clone()));
            continue;
        }
        if outgoing.is_empty() {
            socle_dimensions.push(dimension);
            inclusion_matrices.push(DenseMatrix::identity(dimension, field.clone()));
            continue;
        }

        let concatenated_columns = outgoing
            .iter()
            .map(|arrow| module.dimension_vector()[arrow.target().index().zero_based()])
            .sum::<usize>();
        let mut entries = Vec::with_capacity(dimension * concatenated_columns);
        for row in 0..dimension {
            for arrow in &outgoing {
                let target = arrow.target().index().zero_based();
                let action = logical_submatrix(
                    &module.arrow_matrices()[arrow.index().zero_based()],
                    dimension,
                    module.dimension_vector()[target],
                )?;
                entries.extend_from_slice(action.row(row));
            }
        }
        let concatenated =
            DenseMatrix::from_flat(dimension, concatenated_columns, entries, field.clone())?;
        let kernel = concatenated.nullspace_rows();
        socle_dimensions.push(kernel.nrows());
        if kernel.nrows() == 0 {
            inclusion_matrices.push(DenseMatrix::zero(1, dimension, field.clone()));
        } else {
            inclusion_matrices.push(kernel);
        }
    }

    let socle = right_module_over_path_algebra(
        module.right_acting_algebra().clone(),
        socle_dimensions,
        &[],
    )?
    .ok_or(PathAlgebraModuleError::GeneratedSocleViolatesRelations)?;
    right_module_hom_over_algebra(socle, module, inclusion_matrices)
}

/// Return QPA's `SocleOfModule(M)` object as the source of its inclusion.
///
/// # Errors
///
/// Propagates [`socle_of_module_inclusion`].
pub fn socle_of_module<F: Field>(
    module: PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, PathAlgebraModuleError> {
    Ok(socle_of_module_inclusion(module)?.source().clone())
}

/// Translate the earlier, runtime-shadowed QPA `RightFacApproximation`
/// method from `modulehom.gi:1735-1750`.
///
/// QPA later installs a higher-ranked method on the synonymous
/// `TraceOfModule` operation; user-facing aliases therefore call
/// [`trace_of_module`]. This helper retains the earlier source algorithm for
/// the operation ledger and method-level differential tests.
///
/// # Errors
///
/// Rejects different acting algebras and propagates Hom, image, and
/// subrepresentation construction failures.
#[doc(hidden)]
pub fn right_fac_approximation_via_all_basis_images<F: Field>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let homomorphisms = hom_over_algebra(module, target)?
        .ok_or(PathAlgebraModuleError::DifferentHomomorphismAlgebras)?;
    let basis = module.basis_elements();
    let mut generators = Vec::new();
    for homomorphism in homomorphisms {
        generators.extend(homomorphism.images_set(&basis)?);
    }
    subrepresentation_inclusion(target.clone(), &generators)
}

/// Translate QPA's selected `TraceOfModule(M, N)` method.
///
/// This is also the operation behind `RightFacMApproximation`,
/// `RightFacApproximation`, and `MinimalRightFacMApproximation`. QPA orders
/// the generators first by a minimal top lift of `M`, then by its ordered Hom
/// basis, before closing their images inside `N`.
///
/// # Errors
///
/// Rejects modules over different algebras and propagates top, Hom, image,
/// zero-module, and subrepresentation construction failures.
pub fn trace_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    if module.right_acting_algebra() != target.right_acting_algebra() {
        return Err(PathAlgebraModuleError::DifferentHomomorphismAlgebras);
    }
    if module.dimension() == 0 || target.dimension() == 0 {
        return zero_submodule_inclusion(target);
    }
    let homomorphisms = hom_over_algebra(module, target)?
        .ok_or(PathAlgebraModuleError::DifferentHomomorphismAlgebras)?;
    if homomorphisms.is_empty() {
        return zero_submodule_inclusion(target);
    }
    let generators = minimal_generating_set_of_module(module.clone())?;
    let mut trace_generators = Vec::with_capacity(generators.len() * homomorphisms.len());
    for generator in &generators {
        for homomorphism in &homomorphisms {
            trace_generators.push(homomorphism.image(generator)?);
        }
    }
    subrepresentation_inclusion(target.clone(), &trace_generators)
}

/// QPA synonym of [`trace_of_module`].
///
/// # Errors
///
/// Propagates [`trace_of_module`].
pub fn right_fac_approximation<F: Field>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    trace_of_module(module, target)
}

/// QPA synonym of [`trace_of_module`].
///
/// # Errors
///
/// Propagates [`trace_of_module`].
pub fn right_fac_m_approximation<F: Field>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    trace_of_module(module, target)
}

/// QPA synonym of [`trace_of_module`].
///
/// # Errors
///
/// Propagates [`trace_of_module`].
pub fn minimal_right_fac_m_approximation<F: Field>(
    module: &PathAlgebraModule<F>,
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    trace_of_module(module, target)
}

fn zero_submodule_inclusion<F: Field>(
    target: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, PathAlgebraModuleError> {
    let vertex_count = target.dimension_vector().len();
    let zero = right_module_over_path_algebra(
        target.right_acting_algebra().clone(),
        vec![0; vertex_count],
        &[],
    )?
    .ok_or(PathAlgebraModuleError::GeneratedZeroModuleViolatesRelations)?;
    zero_module_mapping(zero, target.clone())
}
