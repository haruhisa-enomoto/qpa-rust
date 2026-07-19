/// Translate QPA's `TensorProductOfModules(M, N)` operation.
///
/// The method uses exactly the first two differentials of QPA's minimal
/// projective resolution of `M`: its projective cover and the projective cover
/// of that map's kernel. It tensors the resulting presentation with `N`,
/// forms GAP's natural quotient by the row space of `d1 tensor 1`, and retains
/// the returned elementary-tensor closure data.
///
/// The source's two broken branches remain observable as typed errors from
/// [`TensorProductOfModulesResult::elementary_tensor`] or this constructor:
/// an unassigned `size_matrix` when `P1 tensor N` is zero, and quotient-only
/// projective-generator extraction in the nontrivial path-algebra branch.
///
/// # Errors
///
/// Requires the second module to act over the opposite of the first module's
/// algebra. Propagates projective-presentation, opposite-algebra, exact module
/// action, relation-space, and matrix failures, plus typed versions of QPA's
/// source-level failures.
#[allow(clippy::too_many_lines)]
pub fn tensor_product_of_modules<F: Field>(
    left: &PathAlgebraModule<F>,
    right: &PathAlgebraModule<F>,
) -> Result<TensorProductOfModulesResult<F>, EnvelopingModuleError> {
    let algebra = left.right_acting_algebra();
    if algebra != &opposite_module_algebra(right.right_acting_algebra())? {
        return Err(EnvelopingModuleError::TensorProductModulesOverWrongAlgebras);
    }
    let field = algebra.field().clone();
    let top_left = top_of_module(left.clone())?.dimension_vector().to_vec();
    let right_dimensions = right.dimension_vector();
    let top_pairing = tensor_dimension_pairing(&top_left, right_dimensions)?;
    if top_pairing == 0 {
        return Ok(TensorProductOfModulesResult {
            dimension: 0,
            ambient_dimension: 1,
            relation_basis: DenseMatrix::from_flat(0, 1, Vec::new(), field.clone())?,
            basis_lifts: DenseMatrix::from_flat(0, 1, Vec::new(), field.clone())?,
            projection: DenseMatrix::from_flat(1, 0, Vec::new(), field)?,
            elementary: ModuleTensorProductElementary::AlwaysAmbientZero,
        });
    }

    // `ProjectiveResolution(left)` is forced through degree one before QPA
    // reads its degree-zero and degree-one differentials.
    let cover = projective_cover(left)?;
    let kernel = kernel_inclusion(&cover)?;
    let first_syzygy_cover = projective_cover(kernel.source())?;
    let first_differential = first_syzygy_cover.compose(&kernel)?;

    let projections = direct_sum_projections(first_differential.range())?;
    let projective_zero_vertices = projections
        .iter()
        .map(|projection| tensor_projective_vertex(projection.range()))
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;

    if tensor_dimension_pairing(
        first_differential.source().dimension_vector(),
        right_dimensions,
    )? == 0
    {
        return Ok(TensorProductOfModulesResult {
            dimension: top_pairing,
            ambient_dimension: top_pairing,
            relation_basis: DenseMatrix::from_flat(0, top_pairing, Vec::new(), field.clone())?,
            basis_lifts: DenseMatrix::identity(top_pairing, field.clone()),
            projection: DenseMatrix::identity(top_pairing, field),
            elementary: ModuleTensorProductElementary::UnboundSizeMatrix,
        });
    }

    let inclusions = direct_sum_inclusions(first_differential.source())?;
    let mut images_of_generators = Vec::with_capacity(inclusions.len());
    let mut projective_one_vertices = Vec::with_capacity(inclusions.len());
    for inclusion in &inclusions {
        let generator = minimal_generating_set_of_module(inclusion.source().clone())?
            .into_iter()
            .next()
            .ok_or(EnvelopingModuleError::EmptyTensorProjectiveGenerator)?;
        let included = inclusion.image(&generator)?;
        images_of_generators.push(first_differential.image(&included)?);
        projective_one_vertices.push(tensor_projective_vertex(inclusion.source())?);
    }
    if matches!(algebra, QuiverAlgebra::Path(_)) {
        return Err(EnvelopingModuleError::TensorProductPathAlgebraP1Extraction);
    }

    let mut opposite_matrix = Vec::with_capacity(images_of_generators.len());
    for image in &images_of_generators {
        let mut row = Vec::with_capacity(projections.len());
        for (projection, &vertex) in projections.iter().zip(&projective_zero_vertices) {
            let component = projection.image(image)?;
            let element = element_in_indec_projective(algebra, &component, vertex)?;
            row.push(opposite_projective_path_element(algebra, &element)?);
        }
        opposite_matrix.push(row);
    }

    // QPA computes these three values even though it never subsequently reads
    // them in this method.
    let _top_of_p1 = top_of_module(first_differential.source().clone())?;
    let _top_of_p0 = top_of_module(first_differential.range().clone())?;
    let _left_dimensions = left.dimension_vector();

    let right_basis = right.basis_elements();
    let mut relation_rows = Vec::new();
    for (row_index, &source_vertex) in projective_one_vertices.iter().enumerate() {
        let source_dimension = right_dimensions[source_vertex - 1];
        let source_offset = right_dimensions[..source_vertex - 1].iter().sum::<usize>();
        let mut blocks = Vec::with_capacity(projective_zero_vertices.len());
        for (column_index, &target_vertex) in projective_zero_vertices.iter().enumerate() {
            let block = if source_dimension == 0 {
                vec![vec![
                    field.zero();
                    right_dimensions[target_vertex - 1].max(1)
                ]]
            } else {
                right_basis[source_offset..source_offset + source_dimension]
                    .iter()
                    .map(|basis_element| {
                        Ok(basis_element
                            .act(&opposite_matrix[row_index][column_index])?
                            .external_representation()[target_vertex - 1]
                            .clone())
                    })
                    .collect::<Result<Vec<_>, EnvelopingModuleError>>()?
            };
            blocks.push(block);
        }
        let row_count = blocks
            .first()
            .ok_or(EnvelopingModuleError::EmptyTensorProjectivePresentation)?
            .len();
        for physical_row in 0..row_count {
            relation_rows.push(
                blocks
                    .iter()
                    .flat_map(|block| block[physical_row].iter().cloned())
                    .collect::<Vec<_>>(),
            );
        }
    }
    let ambient_dimension = relation_rows
        .first()
        .ok_or(EnvelopingModuleError::EmptyTensorRelationMatrix)?
        .len();
    let relation_matrix = DenseMatrix::from_rows(relation_rows, field)?;
    let quotient = natural_tensor_quotient(&relation_matrix)?;
    let dimension = quotient.basis_lifts.nrows();
    Ok(TensorProductOfModulesResult {
        dimension,
        ambient_dimension,
        relation_basis: quotient.relation_basis,
        basis_lifts: quotient.basis_lifts,
        projection: quotient.projection,
        elementary: ModuleTensorProductElementary::Presentation {
            cover,
            projections,
            projective_vertices: projective_zero_vertices,
            closure_projection_count: *projective_one_vertices
                .last()
                .ok_or(EnvelopingModuleError::EmptyTensorProjectivePresentation)?,
        },
    })
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
