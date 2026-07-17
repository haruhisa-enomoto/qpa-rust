//! QPA's conversion of a finite-dimensional algebra to a quiver algebra.
//!
//! The implementation follows `AlgebraAsQuiverAlgebra` in `lib/algebra.gi`.
//! The first field contract is the split elementary prime-field branch: this
//! is exactly the branch used by QPA once its constructed field `K` equals the
//! entered prime field. Extension-field scalar restriction remains explicit.

use numerica::domains::Field;
use qpa_core::{
    CoefficientDomainError, DenseMatrix, MatrixError, ModuleDecompositionError,
    ModuleEndomorphismAlgebra, PathAlgebra, PathAlgebraElement, PathAlgebraError,
    PathAlgebraModule, PathAlgebraModuleError, PathAlgebraModuleHomomorphism, PathEndpoint,
    QpaCoefficientField, Quiver, QuiverAlgebra, QuiverError, QuotientError,
    StructureConstantsAlgebra, StructureConstantsError, cokernel_of_module_homomorphism,
    direct_sum_inclusions, direct_sum_of_qpa_modules, direct_sum_projections, end_over_algebra,
    from_end_m_to_hom_mm, from_hom_mm_to_end_m, hom_over_algebra, identity_module_mapping,
    loewy_length_of_module, make_basic_list_of_modules, radical_of_module_inclusion,
};
use thiserror::Error;

use crate::{
    EnvelopingModuleError, PathalgGbnpError, indec_projective_modules, loewy_length_of_algebra,
    quotient_by_relations, quotient_is_finite_dimensional,
};

/// QPA's `[KQ/I, images]` return value from `AlgebraAsQuiverAlgebra`.
#[derive(Clone, Debug)]
pub struct AlgebraAsQuiverAlgebraResult<F: Field> {
    algebra: QuiverAlgebra<F>,
    generator_images: Vec<Vec<F::Element>>,
    loewy_length: usize,
}

impl<F: Field> AlgebraAsQuiverAlgebraResult<F> {
    /// Return the path algebra or quotient path algebra constructed by QPA.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Return images of vertices followed by arrows in the entered algebra's
    /// standard coordinates.
    #[must_use]
    pub fn generator_images(&self) -> &[Vec<F::Element>] {
        &self.generator_images
    }

    /// Return QPA's `Length(RadicalSeriesOfAlgebra(A)) - 1` exponent.
    #[must_use]
    pub const fn loewy_length(&self) -> usize {
        self.loewy_length
    }

    /// Consume the result and return QPA's two public list components.
    #[must_use]
    pub fn into_parts(self) -> (QuiverAlgebra<F>, Vec<Vec<F::Element>>) {
        (self.algebra, self.generator_images)
    }
}

/// QPA's successful `EndOfModuleAsQuiverAlgebra` return value.
#[derive(Clone, Debug)]
pub struct EndOfModuleAsQuiverAlgebraResult<F: Field> {
    endomorphism_algebra: ModuleEndomorphismAlgebra<F>,
    quiver_algebra: QuiverAlgebra<F>,
}

impl<F: Field> EndOfModuleAsQuiverAlgebraResult<F> {
    /// Return the entered module's QPA endomorphism algebra.
    #[must_use]
    pub const fn endomorphism_algebra(&self) -> &ModuleEndomorphismAlgebra<F> {
        &self.endomorphism_algebra
    }

    /// Return the adjacency matrix of the converted endomorphism quiver.
    #[must_use]
    pub fn adjacency_matrix(&self) -> &[Box<[usize]>] {
        self.quiver_algebra.quiver().adjacency_matrix()
    }

    /// Return the endomorphism ring as QPA's quiver algebra.
    #[must_use]
    pub const fn quiver_algebra(&self) -> &QuiverAlgebra<F> {
        &self.quiver_algebra
    }

    /// Consume QPA's three-component result in the same order.
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        ModuleEndomorphismAlgebra<F>,
        Vec<Vec<usize>>,
        QuiverAlgebra<F>,
    ) {
        let adjacency = self
            .quiver_algebra
            .quiver()
            .adjacency_matrix()
            .iter()
            .map(|row| row.to_vec())
            .collect();
        (self.endomorphism_algebra, adjacency, self.quiver_algebra)
    }
}

/// QPA's quiver-algebra presentation of the endomorphism ring of a basic sum.
#[derive(Clone, Debug)]
pub struct EndOfBasicModuleAsQuiverAlgebraResult<F: Field> {
    algebra: QuiverAlgebra<F>,
    search_bound: usize,
    final_search_length: usize,
    relations_may_be_incomplete: bool,
}

impl<F: Field> EndOfBasicModuleAsQuiverAlgebraResult<F> {
    /// Return QPA's quotient of the constructed endomorphism quiver.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Return the computed Loewy length or the explicitly entered bound.
    #[must_use]
    pub const fn search_bound(&self) -> usize {
        self.search_bound
    }

    /// Return QPA's final `len` value after its relation-search loop.
    #[must_use]
    pub const fn final_search_length(&self) -> usize {
        self.final_search_length
    }

    /// Return whether QPA emits its bounded-search incompleteness warning.
    #[must_use]
    pub const fn relations_may_be_incomplete(&self) -> bool {
        self.relations_may_be_incomplete
    }

    /// Consume the metadata wrapper and return QPA's public algebra value.
    #[must_use]
    pub fn into_algebra(self) -> QuiverAlgebra<F> {
        self.algebra
    }
}

/// Translate QPA's unbounded `EndOfBasicModuleAsQuiverAlgebra(LL)` call.
///
/// The entered list is assumed basic, exactly as in QPA. The implementation
/// computes radical Hom spaces, arrows from `rad/rad^2`, embeds them in the
/// endomorphism ring of the ordered direct sum, and follows QPA's path
/// dictionary relation search through the default GBNP quotient backend.
///
/// # Errors
///
/// Reproduces QPA's nonempty, common-algebra, and K-elementary checks and
/// propagates exact Hom, radical, matrix, quiver, path, and GBNP failures.
pub fn end_of_basic_module_as_quiver_algebra<F: QpaCoefficientField>(
    modules: &[PathAlgebraModule<F>],
) -> Result<EndOfBasicModuleAsQuiverAlgebraResult<F>, EndOfBasicModuleAsQuiverAlgebraError> {
    end_of_basic_module_as_quiver_algebra_impl(modules, None)
}

/// Translate QPA's bounded `EndOfBasicModuleAsQuiverAlgebra([LL, bound])` call.
///
/// QPA's warning is represented by
/// [`EndOfBasicModuleAsQuiverAlgebraResult::relations_may_be_incomplete`]
/// instead of being written to a process-global GAP info stream.
///
/// # Errors
///
/// In addition to the unbounded errors, rejects zero because GAP's second
/// list entry must satisfy `IsPosInt`.
pub fn end_of_basic_module_as_quiver_algebra_with_bound<F: QpaCoefficientField>(
    modules: &[PathAlgebraModule<F>],
    bound: usize,
) -> Result<EndOfBasicModuleAsQuiverAlgebraResult<F>, EndOfBasicModuleAsQuiverAlgebraError> {
    if bound == 0 {
        return Err(EndOfBasicModuleAsQuiverAlgebraError::NonPositiveSearchBound);
    }
    end_of_basic_module_as_quiver_algebra_impl(modules, Some(bound))
}

#[allow(clippy::too_many_lines)]
fn end_of_basic_module_as_quiver_algebra_impl<F: QpaCoefficientField>(
    modules: &[PathAlgebraModule<F>],
    explicit_bound: Option<usize>,
) -> Result<EndOfBasicModuleAsQuiverAlgebraResult<F>, EndOfBasicModuleAsQuiverAlgebraError> {
    let Some(first) = modules.first() else {
        return Err(EndOfBasicModuleAsQuiverAlgebraError::EmptyModuleList);
    };
    if modules
        .iter()
        .any(|module| module.right_acting_algebra() != first.right_acting_algebra())
    {
        return Err(EndOfBasicModuleAsQuiverAlgebraError::DifferentActingAlgebras);
    }
    let field = first.right_acting_algebra().field().clone();
    let module_count = modules.len();

    // `homs`, `endofdiagonal`, and the diagonal radical replacement in
    // QPA's `radhoms := ShallowCopy(homs)`.
    let mut homs = Vec::with_capacity(module_count);
    for source in modules {
        let mut row = Vec::with_capacity(module_count);
        for range in modules {
            row.push(
                hom_over_algebra(source, range)?.ok_or(
                    EndOfBasicModuleAsQuiverAlgebraError::InternalCommonAlgebraHomRejected,
                )?,
            );
        }
        homs.push(row);
    }
    let mut radical_homs = homs.clone();
    for (index, module) in modules.iter().enumerate() {
        let endomorphisms = end_over_algebra(module)?;
        let structure_constants = endomorphisms.structure_constants_algebra()?;
        let radical = structure_constants.radical()?;
        if endomorphisms
            .dimension()
            .saturating_sub(radical.basis().nrows())
            > 1
        {
            return Err(
                EndOfBasicModuleAsQuiverAlgebraError::EndomorphismRingNotKElementary {
                    module: index + 1,
                },
            );
        }
        let mut diagonal = Vec::with_capacity(radical.basis().nrows());
        for row in 0..radical.basis().nrows() {
            let matrix = endomorphisms.matrix_from_coefficients(radical.basis().row(row))?;
            diagonal.push(from_end_m_to_hom_mm(module, &matrix)?);
        }
        radical_homs[index][index] = diagonal;
    }

    // Preserve the source's outer `radhoms[i][j] <> []` guard and its
    // source/intermediate/target composition order.
    let mut radical_square_homs = (0..module_count)
        .map(|_| (0..module_count).map(|_| Vec::new()).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    for source in 0..module_count {
        for range in 0..module_count {
            if radical_homs[source][range].is_empty() {
                continue;
            }
            for intermediate in 0..module_count {
                for left in &radical_homs[source][intermediate] {
                    for right in &radical_homs[intermediate][range] {
                        radical_square_homs[source][range].push(left.compose(right)?);
                    }
                }
            }
        }
    }

    // QPA's quotient spaces R1/R2 determine both the row-major adjacency
    // matrix and the representative module maps corresponding to arrows.
    let mut adjacency = vec![vec![0usize; module_count]; module_count];
    let mut arrow_maps = Vec::new();
    for source in 0..module_count {
        for range in 0..module_count {
            let Some(first_hom) = homs[source][range].first() else {
                continue;
            };
            let width = first_hom.flat_vertex_matrix_entries().len();
            let radical_rows = unique_rows(
                radical_homs[source][range]
                    .iter()
                    .map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries),
            );
            let radical_square_rows = unique_rows(
                radical_square_homs[source][range]
                    .iter()
                    .map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries),
            );
            let radical_basis = span_basis(&field, width, radical_rows)?;
            let radical_square_basis = span_basis(&field, width, radical_square_rows)?;
            let representatives =
                quotient_basis_representatives(&field, &radical_basis, &radical_square_basis)?;
            adjacency[source][range] = representatives.nrows();

            let hom_basis = DenseMatrix::from_flat(
                homs[source][range].len(),
                width,
                homs[source][range]
                    .iter()
                    .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
                    .collect(),
                field.clone(),
            )?;
            for row in 0..representatives.nrows() {
                let coefficients = hom_basis.solution_row(representatives.row(row))?.ok_or(
                    EndOfBasicModuleAsQuiverAlgebraError::ArrowRepresentativeOutsideHomSpace,
                )?;
                arrow_maps.push(linear_combination_of_homomorphisms(
                    &homs[source][range],
                    &coefficients,
                )?);
            }
        }
    }

    let quiver = Quiver::from_adjacency(adjacency)?;
    let path_algebra = PathAlgebra::new(field.clone(), quiver)?;
    let direct_sum = direct_sum_of_qpa_modules(modules)?
        .ok_or(EndOfBasicModuleAsQuiverAlgebraError::InternalNonemptyDirectSumFailed)?;
    let inclusions = direct_sum_inclusions(&direct_sum)?;
    let projections = direct_sum_projections(&direct_sum)?;
    let endomorphisms = end_over_algebra(&direct_sum)?;
    let endomorphism_constants = endomorphisms.structure_constants_algebra()?;

    let mut primitive_idempotents = Vec::with_capacity(module_count);
    for (projection, inclusion) in projections.iter().zip(&inclusions) {
        let endomorphism = projection.compose(inclusion)?;
        primitive_idempotents.push(endomorphism_coordinates(
            &endomorphisms,
            &from_hom_mm_to_end_m(&endomorphism)?,
        )?);
    }

    let mut endomorphism_arrows = Vec::with_capacity(arrow_maps.len());
    let mut arrow_index = 0;
    for (source, projection) in projections.iter().enumerate() {
        for (range, inclusion) in inclusions.iter().enumerate() {
            for _ in 0..path_algebra.quiver().adjacency_matrix()[source][range] {
                let embedded = projection
                    .compose(&arrow_maps[arrow_index])?
                    .compose(inclusion)?;
                endomorphism_arrows.push(endomorphism_coordinates(
                    &endomorphisms,
                    &from_hom_mm_to_end_m(&embedded)?,
                )?);
                arrow_index += 1;
            }
        }
    }

    let mut radical_endomorphisms = Vec::new();
    for (source, projection) in projections.iter().enumerate() {
        for (range, inclusion) in inclusions.iter().enumerate() {
            for radical_map in &radical_homs[source][range] {
                let embedded = projection.compose(radical_map)?.compose(inclusion)?;
                radical_endomorphisms.push(endomorphism_coordinates(
                    &endomorphisms,
                    &from_hom_mm_to_end_m(&embedded)?,
                )?);
            }
        }
    }
    let radical_basis = span_basis(&field, endomorphisms.dimension(), radical_endomorphisms)?;
    let search_bound = match explicit_bound {
        Some(bound) => bound,
        None => radical_nilpotency_index_from_basis(&endomorphism_constants, &radical_basis)?,
    };

    let vertex_paths = path_algebra
        .quiver()
        .vertices()
        .map(|vertex| path_algebra.element(vertex.path()))
        .collect::<Result<Vec<_>, _>>()?;
    let arrow_paths = path_algebra
        .quiver()
        .arrows()
        .map(|arrow| path_algebra.element(arrow.path()))
        .collect::<Result<Vec<_>, _>>()?;
    let mut layers = vec![
        primitive_idempotents
            .into_iter()
            .zip(vertex_paths)
            .collect::<Vec<_>>(),
        endomorphism_arrows
            .iter()
            .cloned()
            .zip(arrow_paths)
            .collect::<Vec<_>>(),
    ];
    let mut relations = Vec::new();
    let mut length = 1usize;
    while layers.last().is_some_and(|layer| !layer.is_empty()) && length <= search_bound {
        let current_layer = layers.last().cloned().unwrap_or_default();
        let mut new_layer: Vec<(Vec<F::Element>, PathAlgebraElement<F>)> = Vec::new();
        for (element, element_path) in current_layer {
            for (arrow, arrow_path) in endomorphism_arrows
                .iter()
                .zip(layers.get(1).into_iter().flatten().map(|(_, path)| path))
            {
                let product_path = element_path.multiply(arrow_path)?;
                if product_path.is_zero() {
                    continue;
                }
                let product = endomorphism_constants.multiply_coordinates(&element, arrow)?;
                let generators = layers
                    .iter()
                    .flatten()
                    .chain(&new_layer)
                    .collect::<Vec<_>>();
                let generator_matrix = DenseMatrix::from_flat(
                    generators.len(),
                    endomorphisms.dimension(),
                    generators
                        .iter()
                        .flat_map(|(coordinates, _)| coordinates.iter().cloned())
                        .collect(),
                    field.clone(),
                )?;
                if let Some(coefficients) = generator_matrix.solution_row(&product)? {
                    relations.push(path_relation(product_path, &generators, &coefficients)?);
                } else {
                    new_layer.push((product, product_path));
                }
            }
        }
        layers.push(new_layer);
        length = length
            .checked_add(1)
            .ok_or(EndOfBasicModuleAsQuiverAlgebraError::SearchLengthOverflow)?;
    }

    let relations_may_be_incomplete = explicit_bound.is_some() && length > search_bound;
    let algebra = QuiverAlgebra::Quotient(quotient_by_relations(path_algebra, &relations)?);
    Ok(EndOfBasicModuleAsQuiverAlgebraResult {
        algebra,
        search_bound,
        final_search_length: length,
        relations_may_be_incomplete,
    })
}

/// QPA's `ADRAlgebraOfAlgebra` result and installed finite-global-dimension fact.
#[derive(Clone, Debug)]
pub struct AdrAlgebraOfAlgebraResult<F: Field> {
    algebra: QuiverAlgebra<F>,
}

impl<F: Field> AdrAlgebraOfAlgebraResult<F> {
    /// Return the ADR quiver algebra.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Return QPA's explicit `SetIsFiniteGlobalDimensionAlgebra(B, true)` fact.
    #[must_use]
    pub const fn is_finite_global_dimension_algebra(&self) -> bool {
        true
    }

    /// Consume the metadata wrapper and return QPA's public algebra value.
    #[must_use]
    pub fn into_algebra(self) -> QuiverAlgebra<F> {
        self.algebra
    }
}

/// Translate QPA's `ADRAlgebraOfAlgebra(A)` construction.
///
/// Starting from QPA's ordered indecomposable projectives, this follows the
/// successive radical-inclusion composites and cokernels, reverses the full
/// list of proper radical quotients, and delegates the presentation to
/// [`end_of_basic_module_as_quiver_algebra`].
///
/// # Errors
///
/// Reproduces the finite-dimensional precondition and propagates projective,
/// Loewy-length, radical, cokernel, and basic-endomorphism failures.
pub fn adr_algebra_of_algebra<F: QpaCoefficientField>(
    algebra: &QuiverAlgebra<F>,
) -> Result<AdrAlgebraOfAlgebraResult<F>, AdrAlgebraError> {
    if !quiver_algebra_is_finite_dimensional(algebra)? {
        return Err(AdrAlgebraError::InfiniteDimensionalAlgebra);
    }
    let loewy_length = loewy_length_of_algebra(algebra)?;
    let mut modules = indec_projective_modules(algebra)?;
    let mut new_modules = proper_radical_quotients(&modules, loewy_length)?;
    new_modules.reverse();
    modules.extend(new_modules);
    let algebra = end_of_basic_module_as_quiver_algebra(&modules)?.into_algebra();
    Ok(AdrAlgebraOfAlgebraResult { algebra })
}

/// Translate QPA's `ADRAlgebraOfBasicModule(L)` construction.
///
/// QPA forms every proper quotient `M_i/rad^j(M_i)`, reverses the generated
/// quotient list, applies `MakeBasicListOfModules`, and then invokes
/// `EndOfBasicModuleAsQuiverAlgebra`. Rust preserves that exact dependency
/// chain and its finite-field restriction.
///
/// # Errors
///
/// Reproduces QPA's nonempty, common-algebra, finite-field, and
/// finite-dimensional checks and propagates radical, cokernel, deterministic
/// decomposition, and basic-endomorphism failures.
pub fn adr_algebra_of_basic_module<F: QpaCoefficientField>(
    modules: &[PathAlgebraModule<F>],
) -> Result<QuiverAlgebra<F>, AdrAlgebraError> {
    let Some(first) = modules.first() else {
        return Err(AdrAlgebraError::EmptyModuleList);
    };
    let algebra = first.right_acting_algebra();
    if modules
        .iter()
        .any(|module| module.right_acting_algebra() != algebra)
    {
        return Err(AdrAlgebraError::DifferentActingAlgebras);
    }
    if algebra.field().qpa_elements()?.is_none() {
        return Err(AdrAlgebraError::FiniteFieldRequired);
    }
    if !quiver_algebra_is_finite_dimensional(algebra)? {
        return Err(AdrAlgebraError::InfiniteDimensionalAlgebra);
    }
    let mut loewy_length = 0usize;
    for module in modules {
        loewy_length = loewy_length.max(loewy_length_of_module(module.clone())?);
    }
    let mut new_modules = proper_radical_quotients(modules, loewy_length)?;
    new_modules.reverse();
    let mut all_modules = modules.to_vec();
    all_modules.extend(new_modules);
    let basic_modules = make_basic_list_of_modules(&all_modules)?;
    Ok(end_of_basic_module_as_quiver_algebra(&basic_modules)?.into_algebra())
}

fn proper_radical_quotients<F: Field>(
    modules: &[PathAlgebraModule<F>],
    loewy_length: usize,
) -> Result<Vec<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    let mut inclusions = modules
        .iter()
        .cloned()
        .map(identity_module_mapping)
        .collect::<Result<Vec<_>, _>>()?;
    let mut quotients = Vec::new();
    for _ in 0..loewy_length {
        let mut new_inclusions = Vec::new();
        for inclusion in inclusions {
            let radical = radical_of_module_inclusion(inclusion.source().clone())?;
            if radical.source().dimension() > 0 {
                let composite = radical.compose(&inclusion)?;
                quotients.push(cokernel_of_module_homomorphism(&composite)?);
                new_inclusions.push(composite);
            }
        }
        inclusions = new_inclusions;
    }
    Ok(quotients)
}

fn quiver_algebra_is_finite_dimensional<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<bool, PathalgGbnpError> {
    match algebra {
        QuiverAlgebra::Path(algebra) => Ok(algebra.has_finite_canonical_basis()),
        QuiverAlgebra::Quotient(algebra) => quotient_is_finite_dimensional(algebra),
    }
}

/// Translate QPA's `EndOfModuleAsQuiverAlgebra(M)` adapter.
///
/// QPA constructs `EndOverAlgebra(M)`, delegates the actual computation to
/// `AlgebraAsQuiverAlgebra`, and returns the endomorphism algebra, converted
/// quiver adjacency matrix, and converted quiver algebra. Rust preserves that
/// exact dependency chain rather than reconstructing a quiver from a module
/// decomposition.
///
/// # Errors
///
/// Propagates endomorphism structure-constant conversion and QPA's finite,
/// indecomposable, basic, elementary prime-field conversion boundaries.
pub fn end_of_module_as_quiver_algebra<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<EndOfModuleAsQuiverAlgebraResult<F>, EndOfModuleAsQuiverAlgebraError> {
    let endomorphism_algebra = end_over_algebra(module)?;
    let structure_constants = endomorphism_algebra.structure_constants_algebra()?;
    let converted = algebra_as_quiver_algebra(&structure_constants)?;
    let (quiver_algebra, _) = converted.into_parts();
    Ok(EndOfModuleAsQuiverAlgebraResult {
        endomorphism_algebra,
        quiver_algebra,
    })
}

/// Translate QPA's `AlgebraAsQuiverAlgebra` on its split prime-field branch.
///
/// QPA checks indecomposability, basicness, and elementarity; lifts primitive
/// idempotents across the radical; obtains arrows from
/// `e_i rad(A)/rad(A)^2 e_j`; constructs `KQ/J^t`; and finds the final
/// relations as the left nullspace of the resulting surjection to `A`. This
/// function preserves that sequence and uses the translated default GBNP
/// backend for both path-algebra quotients.
///
/// The current Numerica contract supports the case in which QPA's common
/// elementary component field is the entered prime field. A common component
/// of larger dimension would make QPA construct `GF(p^r)` and is reported
/// explicitly until an extension-field backend exists.
///
/// # Errors
///
/// Reproduces QPA's finite-field, indecomposable, basic, and elementary
/// rejections; reports the deferred extension-field branch; and propagates
/// exact structure-constant, matrix, quiver, path, quotient, and GBNP errors.
#[allow(clippy::too_many_lines)]
pub fn algebra_as_quiver_algebra<F: QpaCoefficientField>(
    algebra: &StructureConstantsAlgebra<F>,
) -> Result<AlgebraAsQuiverAlgebraResult<F>, AlgebraAsQuiverAlgebraError> {
    if algebra.field().qpa_elements()?.is_none() {
        return Err(AlgebraAsQuiverAlgebraError::InfiniteCoefficientField);
    }
    let identity = algebra
        .find_multiplicative_identity()?
        .ok_or(AlgebraAsQuiverAlgebraError::MissingMultiplicativeIdentity)?;

    // QPA first decomposes the center and rejects more than one block.
    let central_blocks = algebra.central_idempotents_split_prime()?;
    if central_blocks.len() > 1 {
        return Err(AlgebraAsQuiverAlgebraError::DecomposableAlgebra);
    }

    let radical = algebra.radical()?;
    let semisimple = algebra.quotient_by_ideal_basis(radical.basis().clone())?;
    let primitive_central = semisimple.algebra().central_idempotents_split_prime()?;
    if primitive_central.is_empty() {
        return Err(AlgebraAsQuiverAlgebraError::MissingSemisimpleComponent);
    }

    // `factoralgebradecomp`: c * BasisVectors(Basis(A/rad A)).
    let mut component_dimensions = Vec::with_capacity(primitive_central.len());
    for idempotent in &primitive_central {
        let mut generators = Vec::with_capacity(semisimple.algebra().dimension());
        for basis_index in 0..semisimple.algebra().dimension() {
            let mut basis = vec![algebra.field().zero(); semisimple.algebra().dimension()];
            basis[basis_index] = algebra.field().one();
            let product = semisimple
                .algebra()
                .multiply_coordinates(idempotent, &basis)?;
            if !is_zero_row(algebra.field(), &product) {
                generators.push(product);
            }
        }
        let component = span_basis(
            algebra.field(),
            semisimple.algebra().dimension(),
            generators,
        )?;
        if !subspace_is_commutative(semisimple.algebra(), &component)? {
            return Err(AlgebraAsQuiverAlgebraError::NonBasicAlgebra);
        }
        component_dimensions.push(component.nrows());
    }
    let component_dimension = component_dimensions[0];
    if component_dimensions
        .iter()
        .any(|&dimension| dimension != component_dimension)
    {
        return Err(AlgebraAsQuiverAlgebraError::NonElementaryAlgebra);
    }
    if component_dimension != 1 {
        return Err(AlgebraAsQuiverAlgebraError::ExtensionFieldRequired {
            degree: component_dimension,
        });
    }

    let vertices = semisimple.lift_complete_set_of_orthogonal_idempotents(&primitive_central)?;

    // ProductSpace(rad A, rad A), followed by A/rad(A)^2.
    let radical_square = product_space_basis(algebra, radical.basis(), radical.basis())?;
    let modulo_square = algebra.quotient_by_ideal_basis(radical_square)?;
    let projected_radical = (0..radical.basis().nrows())
        .map(|row| modulo_square.project_coordinates(radical.basis().row(row)))
        .collect::<Result<Vec<_>, _>>()?;
    let radical_modulo_square = span_basis(
        algebra.field(),
        modulo_square.algebra().dimension(),
        projected_radical,
    )?;
    let quotient_vertices = vertices
        .iter()
        .map(|vertex| modulo_square.project_coordinates(vertex))
        .collect::<Result<Vec<_>, _>>()?;

    // QPA's row-major arrows[i][j] list and adjacency matrix.
    let vertex_count = vertices.len();
    let mut adjacency = vec![vec![0usize; vertex_count]; vertex_count];
    let mut arrows = Vec::new();
    for left_vertex in 0..vertex_count {
        for right_vertex in 0..vertex_count {
            let mut candidates = Vec::new();
            for row in 0..radical_modulo_square.nrows() {
                let left_product = modulo_square.algebra().multiply_coordinates(
                    &quotient_vertices[left_vertex],
                    radical_modulo_square.row(row),
                )?;
                let component = modulo_square
                    .algebra()
                    .multiply_coordinates(&left_product, &quotient_vertices[right_vertex])?;
                if !is_zero_row(algebra.field(), &component) {
                    candidates.push(component);
                }
            }
            let arrow_space = span_basis(
                algebra.field(),
                modulo_square.algebra().dimension(),
                candidates,
            )?;
            adjacency[left_vertex][right_vertex] = arrow_space.nrows();
            for row in 0..arrow_space.nrows() {
                let representative = modulo_square.lift_coordinates(arrow_space.row(row))?;
                let left_product =
                    algebra.multiply_coordinates(&vertices[left_vertex], &representative)?;
                arrows.push(algebra.multiply_coordinates(&left_product, &vertices[right_vertex])?);
            }
        }
    }

    let loewy_length = radical_nilpotency_index(algebra, radical.basis())?;
    let exponent = i64::try_from(loewy_length)
        .map_err(|_| AlgebraAsQuiverAlgebraError::LoewyLengthOverflow)?;
    let quiver = Quiver::from_adjacency(adjacency)?;
    let path_algebra = PathAlgebra::new(algebra.field().clone(), quiver)?;
    let arrow_power = path_algebra.nth_power_of_arrow_ideal(exponent)?;
    let mut generator_images = vertices;
    generator_images.extend(arrows);

    // QPA computes a basis of AA = KQ/J^t, applies the generator map to each
    // basis vector, and takes the left nullspace of that image matrix.
    let (basis_representatives, initial_algebra) = if arrow_power.is_empty() {
        let basis = path_algebra
            .canonical_basis()
            .ok_or(AlgebraAsQuiverAlgebraError::MissingTruncatedBasis)?;
        (basis, QuiverAlgebra::Path(path_algebra.clone()))
    } else {
        let truncated = quotient_by_relations(path_algebra.clone(), &arrow_power)?;
        let basis = truncated
            .canonical_basis()?
            .ok_or(AlgebraAsQuiverAlgebraError::MissingTruncatedBasis)?;
        (
            basis
                .into_iter()
                .map(|element| element.representative().clone())
                .collect(),
            QuiverAlgebra::Quotient(truncated),
        )
    };
    let image_rows = basis_representatives
        .iter()
        .map(|basis| {
            image_of_path_element(
                algebra,
                identity.coordinates(),
                basis,
                vertex_count,
                &generator_images,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let image_matrix = DenseMatrix::from_flat(
        image_rows.len(),
        algebra.dimension(),
        image_rows.into_iter().flatten().collect(),
        algebra.field().clone(),
    )?;
    let kernel = image_matrix.nullspace_rows();
    let mut ideal_generators = Vec::with_capacity(kernel.nrows() + arrow_power.len());
    for row in 0..kernel.nrows() {
        ideal_generators.push(linear_combination_of_path_basis(
            &path_algebra,
            &basis_representatives,
            kernel.row(row),
        )?);
    }
    if !arrow_power.is_empty() {
        ideal_generators.extend(arrow_power);
    }

    let result_algebra = if ideal_generators.is_empty() {
        initial_algebra
    } else {
        QuiverAlgebra::Quotient(quotient_by_relations(path_algebra, &ideal_generators)?)
    };
    Ok(AlgebraAsQuiverAlgebraResult {
        algebra: result_algebra,
        generator_images,
        loewy_length,
    })
}

fn unique_rows<T>(rows: impl IntoIterator<Item = Vec<T>>) -> Vec<Vec<T>>
where
    T: PartialEq,
{
    let mut unique = Vec::new();
    for row in rows {
        if !unique.contains(&row) {
            unique.push(row);
        }
    }
    unique
}

fn quotient_basis_representatives<F: Field>(
    field: &F,
    space: &DenseMatrix<F>,
    subspace: &DenseMatrix<F>,
) -> Result<DenseMatrix<F>, EndOfBasicModuleAsQuiverAlgebraError> {
    let mut subspace_coordinates = Vec::with_capacity(subspace.nrows() * space.nrows());
    for row in 0..subspace.nrows() {
        subspace_coordinates.extend(
            space
                .solution_row(subspace.row(row))?
                .ok_or(EndOfBasicModuleAsQuiverAlgebraError::RadicalSquareOutsideRadical)?,
        );
    }
    let kernel = DenseMatrix::from_flat(
        subspace.nrows(),
        space.nrows(),
        subspace_coordinates,
        field.clone(),
    )?
    .semi_echelon_transformation()
    .vectors;
    let mut basis = kernel;
    let mut lifts = Vec::new();
    for coordinate in 0..space.nrows() {
        let mut standard = vec![field.zero(); space.nrows()];
        standard[coordinate] = field.one();
        let mut candidate_entries = basis.as_flat().to_vec();
        candidate_entries.extend_from_slice(&standard);
        let candidate = DenseMatrix::from_flat(
            basis.nrows() + 1,
            space.nrows(),
            candidate_entries,
            field.clone(),
        )?
        .semi_echelon_transformation()
        .vectors;
        if candidate.nrows() > basis.nrows() {
            lifts.extend(standard);
            basis = candidate;
        }
    }
    DenseMatrix::from_flat(
        lifts.len().checked_div(space.nrows()).unwrap_or(0),
        space.nrows(),
        lifts,
        field.clone(),
    )?
    .multiply(space)
    .map_err(Into::into)
}

fn linear_combination_of_homomorphisms<F: Field>(
    basis: &[PathAlgebraModuleHomomorphism<F>],
    coefficients: &[F::Element],
) -> Result<PathAlgebraModuleHomomorphism<F>, EndOfBasicModuleAsQuiverAlgebraError> {
    let first = basis
        .first()
        .ok_or(EndOfBasicModuleAsQuiverAlgebraError::EmptyArrowHomBasis)?;
    let mut result = first.zero()?;
    for (homomorphism, coefficient) in basis.iter().zip(coefficients) {
        result = result.add(&homomorphism.scale(coefficient)?)?;
    }
    Ok(result)
}

fn endomorphism_coordinates<F: Field>(
    algebra: &ModuleEndomorphismAlgebra<F>,
    matrix: &DenseMatrix<F>,
) -> Result<Vec<F::Element>, EndOfBasicModuleAsQuiverAlgebraError> {
    algebra
        .coefficients(matrix)?
        .ok_or(EndOfBasicModuleAsQuiverAlgebraError::EmbeddedMapOutsideEndomorphismAlgebra)
}

fn radical_nilpotency_index_from_basis<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    radical: &DenseMatrix<F>,
) -> Result<usize, EndOfBasicModuleAsQuiverAlgebraError> {
    let mut power = radical.clone();
    let mut exponent = 1usize;
    while power.nrows() > 0 {
        let mut products = Vec::with_capacity(power.nrows().saturating_mul(radical.nrows()));
        for left in 0..power.nrows() {
            for right in 0..radical.nrows() {
                products.push(algebra.multiply_coordinates(power.row(left), radical.row(right))?);
            }
        }
        power = span_basis(algebra.field(), algebra.dimension(), products)?;
        exponent = exponent
            .checked_add(1)
            .ok_or(EndOfBasicModuleAsQuiverAlgebraError::SearchLengthOverflow)?;
        if exponent > algebra.dimension().saturating_add(1) {
            return Err(EndOfBasicModuleAsQuiverAlgebraError::RadicalNotNilpotent);
        }
    }
    Ok(exponent)
}

fn path_relation<F: Field>(
    product_path: PathAlgebraElement<F>,
    generators: &[&(Vec<F::Element>, PathAlgebraElement<F>)],
    coefficients: &[F::Element],
) -> Result<PathAlgebraElement<F>, PathAlgebraError> {
    let field = product_path.algebra().field().clone();
    let mut relation = product_path;
    for ((_, path), coefficient) in generators.iter().copied().zip(coefficients) {
        if !field.is_zero(coefficient) {
            relation = relation.sub(&path.scale_left(coefficient)?)?;
        }
    }
    Ok(relation)
}

fn span_basis<F: Field>(
    field: &F,
    width: usize,
    rows: Vec<Vec<F::Element>>,
) -> Result<DenseMatrix<F>, MatrixError> {
    DenseMatrix::from_flat(
        rows.len(),
        width,
        rows.into_iter().flatten().collect(),
        field.clone(),
    )
    .map(|matrix| matrix.semi_echelon_transformation().vectors)
}

fn product_space_basis<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    left: &DenseMatrix<F>,
    right: &DenseMatrix<F>,
) -> Result<DenseMatrix<F>, AlgebraAsQuiverAlgebraError> {
    let mut products = Vec::with_capacity(left.nrows().saturating_mul(right.nrows()));
    for left_row in 0..left.nrows() {
        for right_row in 0..right.nrows() {
            products.push(algebra.multiply_coordinates(left.row(left_row), right.row(right_row))?);
        }
    }
    Ok(span_basis(algebra.field(), algebra.dimension(), products)?)
}

fn subspace_is_commutative<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    basis: &DenseMatrix<F>,
) -> Result<bool, StructureConstantsError> {
    for left in 0..basis.nrows() {
        for right in left + 1..basis.nrows() {
            if algebra.multiply_coordinates(basis.row(left), basis.row(right))?
                != algebra.multiply_coordinates(basis.row(right), basis.row(left))?
            {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn radical_nilpotency_index<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    radical: &DenseMatrix<F>,
) -> Result<usize, AlgebraAsQuiverAlgebraError> {
    let mut power = radical.clone();
    let mut exponent = 1usize;
    while power.nrows() > 0 {
        power = product_space_basis(algebra, &power, radical)?;
        exponent = exponent
            .checked_add(1)
            .ok_or(AlgebraAsQuiverAlgebraError::LoewyLengthOverflow)?;
        if exponent > algebra.dimension().saturating_add(1) {
            return Err(AlgebraAsQuiverAlgebraError::RadicalNotNilpotent);
        }
    }
    Ok(exponent)
}

fn image_of_path_element<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    identity: &[F::Element],
    element: &PathAlgebraElement<F>,
    vertex_count: usize,
    images: &[Vec<F::Element>],
) -> Result<Vec<F::Element>, AlgebraAsQuiverAlgebraError> {
    let mut total = vec![algebra.field().zero(); algebra.dimension()];
    for term in element.terms() {
        let walk = term
            .path()
            .walk()
            .ok_or(AlgebraAsQuiverAlgebraError::ZeroPathInTruncatedBasis)?;
        let mut image = if walk.is_empty() {
            let PathEndpoint::Vertex(vertex) = term.path().source() else {
                return Err(AlgebraAsQuiverAlgebraError::InvalidVertexBasisPath);
            };
            images
                .get(vertex.index().get() - 1)
                .cloned()
                .ok_or(AlgebraAsQuiverAlgebraError::MissingGeneratorImage)?
        } else {
            let mut product = identity.to_vec();
            for arrow in walk {
                let arrow_image = images
                    .get(vertex_count + arrow.get() - 1)
                    .ok_or(AlgebraAsQuiverAlgebraError::MissingGeneratorImage)?;
                product = algebra.multiply_coordinates(&product, arrow_image)?;
            }
            product
        };
        for coefficient in &mut image {
            algebra.field().mul_assign(coefficient, term.coefficient());
        }
        add_row_assign(algebra.field(), &mut total, &image);
    }
    Ok(total)
}

fn linear_combination_of_path_basis<F: Field>(
    algebra: &PathAlgebra<F>,
    basis: &[PathAlgebraElement<F>],
    coefficients: &[F::Element],
) -> Result<PathAlgebraElement<F>, PathAlgebraError> {
    algebra.from_terms(
        basis
            .iter()
            .zip(coefficients)
            .flat_map(|(element, scalar)| {
                element.terms().iter().map(move |term| {
                    (
                        term.path().clone(),
                        algebra.field().mul(scalar, term.coefficient()),
                    )
                })
            }),
    )
}

fn is_zero_row<F: Field>(field: &F, row: &[F::Element]) -> bool {
    row.iter().all(|entry| field.is_zero(entry))
}

fn add_row_assign<F: Field>(field: &F, target: &mut [F::Element], source: &[F::Element]) {
    for (target, source) in target.iter_mut().zip(source) {
        field.add_assign(target, source);
    }
}

/// Failure in QPA's two ADR-algebra constructors.
#[derive(Debug, Error)]
pub enum AdrAlgebraError {
    /// `ADRAlgebraOfBasicModule` requires a nonempty list.
    #[error("the entered module list is empty")]
    EmptyModuleList,
    /// Every basic-module input must use one structurally equal algebra.
    #[error("not all modules are over the same algebra")]
    DifferentActingAlgebras,
    /// QPA's basic-module overload is restricted to finite fields.
    #[error("the entered module is not over an algebra over a finite field")]
    FiniteFieldRequired,
    /// Both source methods explicitly require finite-dimensional algebras.
    #[error("the entered algebra is not a finite-dimensional algebra")]
    InfiniteDimensionalAlgebra,
    /// Coefficient-field enumeration failed.
    #[error(transparent)]
    Coefficient(#[from] CoefficientDomainError),
    /// The configured QPA/GBNP finite-dimensionality path failed.
    #[error(transparent)]
    Backend(#[from] PathalgGbnpError),
    /// QPA's projective or algebra-Loewy construction failed.
    #[error(transparent)]
    Enveloping(#[from] EnvelopingModuleError),
    /// Exact radical, cokernel, module-Loewy, or map construction failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// `MakeBasicListOfModules` failed in the basic-module overload.
    #[error(transparent)]
    Decomposition(#[from] ModuleDecompositionError),
    /// The final QPA basic-endomorphism presentation failed.
    #[error(transparent)]
    BasicEndomorphism(#[from] EndOfBasicModuleAsQuiverAlgebraError),
}

/// Failure in QPA's basic-module endomorphism presentation algorithm.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum EndOfBasicModuleAsQuiverAlgebraError {
    /// QPA requires a nonempty dense list.
    #[error("the entered module list must be nonempty")]
    EmptyModuleList,
    /// QPA requires all modules to have one structurally equal acting algebra.
    #[error("not all entered modules are modules over the same algebra")]
    DifferentActingAlgebras,
    /// GAP's optional second list entry must satisfy `IsPosInt`.
    #[error("the maximum relation-search length must be positive")]
    NonPositiveSearchBound,
    /// A common-algebra Hom query cannot take its printed-then-fail branch.
    #[error("an internal HomOverAlgebra call rejected the checked common algebra")]
    InternalCommonAlgebraHomRejected,
    /// QPA checks every diagonal endomorphism-ring semisimple quotient.
    #[error("the endomorphism ring of module {module} is not K-elementary")]
    EndomorphismRingNotKElementary {
        /// One-based module position in the entered list.
        module: usize,
    },
    /// `rad^2 Hom` must be contained in `rad Hom`.
    #[error("an internal radical-square map lies outside the radical Hom space")]
    RadicalSquareOutsideRadical,
    /// Every quotient representative was constructed inside its full Hom space.
    #[error("an arrow representative lies outside its full Hom space")]
    ArrowRepresentativeOutsideHomSpace,
    /// An arrow branch necessarily has a nonempty full Hom basis.
    #[error("an internal arrow linear combination has an empty Hom basis")]
    EmptyArrowHomBasis,
    /// A nonempty input must produce QPA's ordered direct sum.
    #[error("an internal nonempty DirectSumOfQPAModules call returned fail")]
    InternalNonemptyDirectSumFailed,
    /// Every embedded projection, arrow, and radical map belongs to `End(M)`.
    #[error("an embedded module map lies outside the represented endomorphism algebra")]
    EmbeddedMapOutsideEndomorphismAlgebra,
    /// The radical of a finite-dimensional endomorphism ring must be nilpotent.
    #[error("the computed radical of the endomorphism ring is not nilpotent")]
    RadicalNotNilpotent,
    /// QPA's relation-search layer counter overflowed Rust's address space.
    #[error("the relation-search length overflows usize")]
    SearchLengthOverflow,
    /// Exact module and module-map construction failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// Exact vector-space construction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// Exact structure-constant arithmetic failed.
    #[error(transparent)]
    Structure(#[from] StructureConstantsError),
    /// The radical-Hom adjacency matrix did not define a QPA quiver.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// Path-algebra construction or relation arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// QPA's default GBNP-backed quotient construction failed.
    #[error(transparent)]
    Pathalg(#[from] PathalgGbnpError),
}

/// Failure in QPA's finite-dimensional-algebra conversion.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum AlgebraAsQuiverAlgebraError {
    /// QPA's method is restricted to finite coefficient fields.
    #[error("the entered algebra is not an algebra over a finite field")]
    InfiniteCoefficientField,
    /// QPA could not find an identity in the entered algebra.
    #[error("the entered algebra does not have a multiplicative identity")]
    MissingMultiplicativeIdentity,
    /// QPA rejects more than one central block.
    #[error("the entered algebra is not indecomposable")]
    DecomposableAlgebra,
    /// QPA's semisimple factor unexpectedly supplied no central component.
    #[error("the semisimple factor has no central component")]
    MissingSemisimpleComponent,
    /// One simple factor algebra is noncommutative.
    #[error("the entered algebra is not basic")]
    NonBasicAlgebra,
    /// The simple factor fields do not all have equal dimension.
    #[error("the entered algebra is not elementary, hence not a quotient of a path algebra")]
    NonElementaryAlgebra,
    /// QPA would construct an extension field not yet supplied by Numerica.
    #[error("the conversion requires a degree-{degree} extension of the prime coefficient field")]
    ExtensionFieldRequired {
        /// Common dimension of the elementary division-field components.
        degree: usize,
    },
    /// `Length(RadicalSeriesOfAlgebra(A)) - 1` did not fit QPA's exponent API.
    #[error("the radical-series length overflows the supported path exponent")]
    LoewyLengthOverflow,
    /// A purported Jacobson radical failed to become zero.
    #[error("the computed radical is not nilpotent")]
    RadicalNotNilpotent,
    /// QPA's `KQ/J^t` basis was unavailable.
    #[error("the truncated path algebra has no finite canonical basis")]
    MissingTruncatedBasis,
    /// A canonical basis must never contain QPA's distinguished zero path.
    #[error("the truncated canonical basis contains the zero path")]
    ZeroPathInTruncatedBasis,
    /// A length-zero canonical path did not identify a vertex.
    #[error("a truncated canonical basis path has no vertex endpoint")]
    InvalidVertexBasisPath,
    /// A basis path referred outside the QPA generator-image list.
    #[error("a truncated canonical basis path has no generator image")]
    MissingGeneratorImage,
    /// Structure-constant arithmetic or decomposition failed.
    #[error(transparent)]
    Structure(#[from] StructureConstantsError),
    /// Exact vector-space construction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// The adjacency matrix did not define a QPA quiver.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// Parent path-algebra construction or arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Quotient canonical-basis construction failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// QPA's default GBNP-backed quotient construction failed.
    #[error(transparent)]
    Pathalg(#[from] PathalgGbnpError),
    /// Coefficient-field enumeration failed.
    #[error(transparent)]
    Coefficient(#[from] qpa_core::CoefficientDomainError),
}

/// Failure while applying QPA's module-endomorphism conversion adapter.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum EndOfModuleAsQuiverAlgebraError {
    /// Constructing or converting `EndOverAlgebra(M)` failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// QPA's delegated `AlgebraAsQuiverAlgebra` computation failed.
    #[error(transparent)]
    Conversion(#[from] AlgebraAsQuiverAlgebraError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, Rational},
    };
    use qpa_core::{
        ArrowSpec, DenseMatrix, NamedArrowMatrix, PathAlgebra, PathAlgebraModule,
        QpaCoefficientField, Quiver, QuiverAlgebra, StructureConstantsAlgebra, VertexIndex,
        direct_sum_of_qpa_modules, right_module_over_path_algebra,
    };
    use serde_json::{Value, json};

    use crate::quotient_by_relations;

    use super::{
        AdrAlgebraError, AlgebraAsQuiverAlgebraError, EndOfBasicModuleAsQuiverAlgebraError,
        EndOfModuleAsQuiverAlgebraError, adr_algebra_of_algebra, adr_algebra_of_basic_module,
        algebra_as_quiver_algebra, end_of_basic_module_as_quiver_algebra,
        end_of_basic_module_as_quiver_algebra_with_bound, end_of_module_as_quiver_algebra,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn set_basis_product<F: Field>(
        tensor: &mut [F::Element],
        field: &F,
        dimension: usize,
        left: usize,
        right: usize,
        output: usize,
    ) {
        tensor[(left * dimension + right) * dimension + output] = field.one();
    }

    fn upper_triangular_two_by_two<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        let dimension = 3;
        let mut tensor = vec![field.zero(); dimension * dimension * dimension];
        set_basis_product(&mut tensor, &field, dimension, 0, 0, 0);
        set_basis_product(&mut tensor, &field, dimension, 0, 1, 1);
        set_basis_product(&mut tensor, &field, dimension, 1, 2, 1);
        set_basis_product(&mut tensor, &field, dimension, 2, 2, 2);
        StructureConstantsAlgebra::new(field, dimension, tensor).unwrap()
    }

    fn dual_numbers<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        let dimension = 2;
        let mut tensor = vec![field.zero(); dimension * dimension * dimension];
        set_basis_product(&mut tensor, &field, dimension, 0, 0, 0);
        set_basis_product(&mut tensor, &field, dimension, 0, 1, 1);
        set_basis_product(&mut tensor, &field, dimension, 1, 0, 1);
        StructureConstantsAlgebra::new(field, dimension, tensor).unwrap()
    }

    fn product_of_two_fields<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        let dimension = 2;
        let mut tensor = vec![field.zero(); dimension * dimension * dimension];
        set_basis_product(&mut tensor, &field, dimension, 0, 0, 0);
        set_basis_product(&mut tensor, &field, dimension, 1, 1, 1);
        StructureConstantsAlgebra::new(field, dimension, tensor).unwrap()
    }

    fn two_by_two_matrix_algebra<F: Field>(field: F) -> StructureConstantsAlgebra<F> {
        // e11, e12, e21, e22 in row-major matrix-unit order.
        let dimension = 4;
        let mut tensor = vec![field.zero(); dimension * dimension * dimension];
        let units = [(0usize, 0usize), (0, 1), (1, 0), (1, 1)];
        for (left, &(i, j)) in units.iter().enumerate() {
            for (right, &(k, l)) in units.iter().enumerate() {
                if j == k {
                    let output = units
                        .iter()
                        .position(|&(row, column)| row == i && column == l)
                        .unwrap();
                    set_basis_product(&mut tensor, &field, dimension, left, right, output);
                }
            }
        }
        StructureConstantsAlgebra::new(field, dimension, tensor).unwrap()
    }

    fn conversion_json<F: qpa_core::QpaCoefficientField>(
        input: &StructureConstantsAlgebra<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        let converted = algebra_as_quiver_algebra(input).unwrap();
        let result = converted.algebra();
        let dimension = match result {
            QuiverAlgebra::Path(algebra) => algebra.canonical_basis().unwrap().len(),
            QuiverAlgebra::Quotient(algebra) => algebra.canonical_basis().unwrap().unwrap().len(),
        };
        let arrows = result
            .quiver()
            .arrows()
            .map(|arrow| json!([arrow.source().index().get(), arrow.target().index().get()]))
            .collect::<Vec<_>>();
        let images = converted
            .generator_images()
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|(_, coefficient)| !input.field().is_zero(coefficient))
                    .map(|(coordinate, coefficient)| {
                        json!([coordinate + 1, coefficient_json(coefficient)])
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let relations = result
            .relations()
            .iter()
            .map(|relation| {
                relation
                    .terms()
                    .iter()
                    .map(|term| {
                        json!({
                            "walk": term.path().walk().unwrap().iter().map(|arrow| arrow.get()).collect::<Vec<_>>(),
                            "coefficient": coefficient_json(term.coefficient())
                        })
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        json!({
            "dimension": dimension,
            "vertices": result.quiver().number_of_vertices(),
            "arrows": arrows,
            "generator_images": images,
            "loewy_length": converted.loewy_length(),
            "relations": relations
        })
    }

    fn small_integer(value: &Integer) -> i64 {
        value.to_string().parse().unwrap()
    }

    fn end_of_module_json<F: QpaCoefficientField>(field: F) -> Value {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field.clone(), quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let regular = right_module_over_path_algebra(
            algebra,
            vec![2],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(
                    vec![
                        vec![field.zero(), field.one()],
                        vec![field.zero(), field.zero()],
                    ],
                    field,
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let converted = end_of_module_as_quiver_algebra(&regular).unwrap();
        let quiver_dimension = match converted.quiver_algebra() {
            QuiverAlgebra::Path(algebra) => algebra.canonical_basis().unwrap().len(),
            QuiverAlgebra::Quotient(algebra) => algebra.canonical_basis().unwrap().unwrap().len(),
        };
        json!({
            "endomorphism_dimension": converted.endomorphism_algebra().dimension(),
            "adjacency": converted.adjacency_matrix().iter()
                .map(|row| row.to_vec()).collect::<Vec<_>>(),
            "quiver_algebra_dimension": quiver_dimension,
        })
    }

    fn basic_end_modules<F: QpaCoefficientField>(field: F) -> Vec<PathAlgebraModule<F>> {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let projective = right_module_over_path_algebra(
            algebra,
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![field.one()]], field).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        vec![simple_one, projective, simple_two]
    }

    fn dual_numbers_generator<F: QpaCoefficientField>(field: F) -> Vec<PathAlgebraModule<F>> {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field.clone(), quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let regular = right_module_over_path_algebra(
            algebra.clone(),
            vec![2],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(
                    vec![
                        vec![field.zero(), field.one()],
                        vec![field.zero(), field.zero()],
                    ],
                    field,
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let simple = right_module_over_path_algebra(algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        vec![regular, simple]
    }

    fn basic_end_json<F: QpaCoefficientField>(
        modules: &[PathAlgebraModule<F>],
        coefficient_json: &impl Fn(&F::Element) -> Value,
        bound: Option<usize>,
    ) -> Value {
        let result = match bound {
            Some(bound) => end_of_basic_module_as_quiver_algebra_with_bound(modules, bound),
            None => end_of_basic_module_as_quiver_algebra(modules),
        }
        .unwrap();
        quiver_algebra_json(result.algebra(), coefficient_json)
    }

    fn quiver_algebra_json<F: Field>(
        algebra: &QuiverAlgebra<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        let dimension = match algebra {
            QuiverAlgebra::Path(algebra) => algebra.canonical_basis().unwrap().len(),
            QuiverAlgebra::Quotient(algebra) => algebra.canonical_basis().unwrap().unwrap().len(),
        };
        json!({
            "dimension": dimension,
            "adjacency": algebra.quiver().adjacency_matrix().iter()
                .map(|row| row.to_vec()).collect::<Vec<_>>(),
            "arrows": algebra.quiver().arrows().map(|arrow| json!([
                arrow.source().index().get(), arrow.target().index().get()
            ])).collect::<Vec<_>>(),
            "relations": algebra.relations().iter().map(|relation| {
                relation.terms().iter().map(|term| json!({
                    "walk": term.path().walk().unwrap().iter()
                        .map(|arrow| arrow.get()).collect::<Vec<_>>(),
                    "coefficient": coefficient_json(term.coefficient())
                })).collect::<Vec<_>>()
            }).collect::<Vec<_>>()
        })
    }

    fn adr_algebra_json<F: QpaCoefficientField>(
        field: F,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        let a2_modules = basic_end_modules(field.clone());
        let a2 = adr_algebra_of_algebra(a2_modules[0].right_acting_algebra()).unwrap();
        let dual_modules = dual_numbers_generator(field);
        let dual = adr_algebra_of_algebra(dual_modules[0].right_acting_algebra()).unwrap();
        json!({
            "a2": quiver_algebra_json(a2.algebra(), coefficient_json),
            "a2_finite_global_dimension_cached":
                a2.is_finite_global_dimension_algebra(),
            "dual_numbers": quiver_algebra_json(dual.algebra(), coefficient_json),
            "dual_numbers_finite_global_dimension_cached":
                dual.is_finite_global_dimension_algebra()
        })
    }

    #[test]
    fn pinned_qpa_conversion_branches_match() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_module_preprojective_structure.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            conversion_json(&upper_triangular_two_by_two(Z2), &|coefficient| {
                small_integer(&Z2.to_integer(coefficient)).into()
            }),
            fixture["conversion_cases"]["upper_triangular_f2"]
        );
        let field = Zp::new(5);
        assert_eq!(
            conversion_json(&dual_numbers(field.clone()), &|coefficient| {
                small_integer(&field.to_integer(coefficient)).into()
            }),
            fixture["conversion_cases"]["dual_numbers_f5"]
        );
    }

    #[test]
    fn qpa_conversion_rejects_infinite_decomposable_and_nonbasic_inputs() {
        assert!(matches!(
            algebra_as_quiver_algebra(&dual_numbers(Q)),
            Err(AlgebraAsQuiverAlgebraError::InfiniteCoefficientField)
        ));
        assert!(matches!(
            algebra_as_quiver_algebra(&product_of_two_fields(Z2)),
            Err(AlgebraAsQuiverAlgebraError::DecomposableAlgebra)
        ));
        assert!(matches!(
            algebra_as_quiver_algebra(&two_by_two_matrix_algebra(Z2)),
            Err(AlgebraAsQuiverAlgebraError::NonBasicAlgebra)
        ));
    }

    #[test]
    fn pinned_qpa_end_of_module_as_quiver_algebra_matches() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_end_quiver.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            end_of_module_json(Z2),
            fixture["binary"]["regular_dual_numbers"]
        );
        assert_eq!(
            end_of_module_json(Zp::new(5)),
            fixture["quinary"]["regular_dual_numbers"]
        );

        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let repeated = right_module_over_path_algebra(algebra, vec![2], &[])
            .unwrap()
            .unwrap();
        assert!(matches!(
            end_of_module_as_quiver_algebra(&repeated),
            Err(EndOfModuleAsQuiverAlgebraError::Conversion(
                AlgebraAsQuiverAlgebraError::NonBasicAlgebra
            ))
        ));

        // QPA reaches the same zero radical here, but its multiplication of
        // the empty radical-basis list aborts before it can return K.  The
        // translated converter handles that zero space explicitly.
        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let simple = right_module_over_path_algebra(algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        let converted = end_of_module_as_quiver_algebra(&simple).unwrap();
        assert_eq!(converted.endomorphism_algebra().dimension(), 1);
        assert_eq!(converted.adjacency_matrix(), &[vec![0].into_boxed_slice()]);
        assert!(matches!(
            converted.quiver_algebra(),
            QuiverAlgebra::Path(algebra) if algebra.canonical_basis().unwrap().len() == 1
        ));
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_end_of_basic_module_as_quiver_algebra_matches() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_basic_end.schema-v1.json"
        ))
        .unwrap();

        let binary = basic_end_modules(Z2);
        assert_eq!(
            basic_end_json(
                &binary,
                &|coefficient| { small_integer(&Z2.to_integer(coefficient)).into() },
                None
            ),
            fixture["binary"]["auslander_a2"]
        );
        assert_eq!(
            basic_end_json(
                &[binary[0].clone(), binary[2].clone()],
                &|coefficient| { small_integer(&Z2.to_integer(coefficient)).into() },
                None
            ),
            fixture["binary"]["semisimple"]
        );
        let bounded = end_of_basic_module_as_quiver_algebra_with_bound(&binary, 1).unwrap();
        assert_eq!(bounded.search_bound(), 1);
        assert_eq!(bounded.final_search_length(), 2);
        assert!(bounded.relations_may_be_incomplete());
        assert_eq!(
            basic_end_json(
                &binary,
                &|coefficient| { small_integer(&Z2.to_integer(coefficient)).into() },
                Some(1)
            ),
            fixture["binary"]["bounded_one"]
        );
        assert_eq!(
            basic_end_json(
                &dual_numbers_generator(Z2),
                &|coefficient| { small_integer(&Z2.to_integer(coefficient)).into() },
                None
            ),
            fixture["binary"]["dual_numbers_generator"]
        );

        let quinary_field = Zp::new(5);
        let quinary = basic_end_modules(quinary_field.clone());
        assert_eq!(
            basic_end_json(
                &quinary,
                &|coefficient| { small_integer(&quinary_field.to_integer(coefficient)).into() },
                None
            ),
            fixture["quinary"]["auslander_a2"]
        );
        assert_eq!(
            basic_end_json(
                &[quinary[0].clone(), quinary[2].clone()],
                &|coefficient| { small_integer(&quinary_field.to_integer(coefficient)).into() },
                None
            ),
            fixture["quinary"]["semisimple"]
        );
        assert_eq!(
            basic_end_json(
                &quinary,
                &|coefficient| { small_integer(&quinary_field.to_integer(coefficient)).into() },
                Some(1)
            ),
            fixture["quinary"]["bounded_one"]
        );
        assert_eq!(
            basic_end_json(
                &dual_numbers_generator(quinary_field.clone()),
                &|coefficient| { small_integer(&quinary_field.to_integer(coefficient)).into() },
                None
            ),
            fixture["quinary"]["dual_numbers_generator"]
        );

        let rational = basic_end_modules(Q);
        let rational_json = |coefficient: &Rational| {
            json!([
                small_integer(coefficient.numerator_ref()),
                small_integer(coefficient.denominator_ref())
            ])
        };
        assert_eq!(
            basic_end_json(&rational, &rational_json, None),
            fixture["rational"]["auslander_a2"]
        );
        assert_eq!(
            basic_end_json(
                &[rational[0].clone(), rational[2].clone()],
                &rational_json,
                None
            ),
            fixture["rational"]["semisimple"]
        );
        assert_eq!(
            basic_end_json(&rational, &rational_json, Some(1)),
            fixture["rational"]["bounded_one"]
        );
        assert_eq!(
            basic_end_json(&dual_numbers_generator(Q), &rational_json, None),
            fixture["rational"]["dual_numbers_generator"]
        );
    }

    #[test]
    fn qpa_end_of_basic_module_validation_branches_match() {
        assert!(matches!(
            end_of_basic_module_as_quiver_algebra::<Z2>(&[]),
            Err(EndOfBasicModuleAsQuiverAlgebraError::EmptyModuleList)
        ));
        let modules = basic_end_modules(Z2);
        assert!(matches!(
            end_of_basic_module_as_quiver_algebra_with_bound(&modules, 0),
            Err(EndOfBasicModuleAsQuiverAlgebraError::NonPositiveSearchBound)
        ));

        let foreign_parent =
            PathAlgebra::new(Z2, Quiver::from_indexed(1, vec![]).unwrap()).unwrap();
        let foreign =
            right_module_over_path_algebra(QuiverAlgebra::Path(foreign_parent), vec![1], &[])
                .unwrap()
                .unwrap();
        assert!(matches!(
            end_of_basic_module_as_quiver_algebra(&[modules[0].clone(), foreign]),
            Err(EndOfBasicModuleAsQuiverAlgebraError::DifferentActingAlgebras)
        ));

        let repeated = direct_sum_of_qpa_modules(&[modules[0].clone(), modules[0].clone()])
            .unwrap()
            .unwrap();
        assert!(matches!(
            end_of_basic_module_as_quiver_algebra(&[repeated]),
            Err(EndOfBasicModuleAsQuiverAlgebraError::EndomorphismRingNotKElementary { module: 1 })
        ));
    }

    #[test]
    fn pinned_qpa_adr_algebra_constructions_match() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_adr.schema-v1.json"
        ))
        .unwrap();

        assert_eq!(
            adr_algebra_json(Z2, &|coefficient| {
                small_integer(&Z2.to_integer(coefficient)).into()
            }),
            fixture["binary"]["algebra"]
        );
        let binary = basic_end_modules(Z2);
        let binary_basic = adr_algebra_of_basic_module(&[binary[1].clone()]).unwrap();
        assert_eq!(
            quiver_algebra_json(&binary_basic, &|coefficient| {
                small_integer(&Z2.to_integer(coefficient)).into()
            }),
            fixture["binary"]["basic_module"]
        );

        let quinary_field = Zp::new(5);
        assert_eq!(
            adr_algebra_json(quinary_field.clone(), &|coefficient| {
                small_integer(&quinary_field.to_integer(coefficient)).into()
            }),
            fixture["quinary"]["algebra"]
        );
        let quinary = basic_end_modules(quinary_field.clone());
        let quinary_basic = adr_algebra_of_basic_module(&[quinary[1].clone()]).unwrap();
        assert_eq!(
            quiver_algebra_json(&quinary_basic, &|coefficient| {
                small_integer(&quinary_field.to_integer(coefficient)).into()
            }),
            fixture["quinary"]["basic_module"]
        );

        let rational_json = |coefficient: &Rational| {
            json!([
                small_integer(coefficient.numerator_ref()),
                small_integer(coefficient.denominator_ref())
            ])
        };
        assert_eq!(
            adr_algebra_json(Q, &rational_json),
            fixture["rational"]["algebra"]
        );
    }

    #[test]
    fn qpa_adr_validation_branches_match() {
        assert!(matches!(
            adr_algebra_of_basic_module::<Z2>(&[]),
            Err(AdrAlgebraError::EmptyModuleList)
        ));

        let binary = basic_end_modules(Z2);
        let foreign_parent =
            PathAlgebra::new(Z2, Quiver::from_indexed(1, vec![]).unwrap()).unwrap();
        let foreign =
            right_module_over_path_algebra(QuiverAlgebra::Path(foreign_parent), vec![1], &[])
                .unwrap()
                .unwrap();
        assert!(matches!(
            adr_algebra_of_basic_module(&[binary[0].clone(), foreign]),
            Err(AdrAlgebraError::DifferentActingAlgebras)
        ));

        let rational = basic_end_modules(Q);
        assert!(matches!(
            adr_algebra_of_basic_module(&[rational[0].clone()]),
            Err(AdrAlgebraError::FiniteFieldRequired)
        ));

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "a")]).unwrap();
        let loop_algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, loop_quiver).unwrap());
        assert!(matches!(
            adr_algebra_of_algebra(&loop_algebra),
            Err(AdrAlgebraError::InfiniteDimensionalAlgebra)
        ));
        let loop_module = right_module_over_path_algebra(loop_algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        assert!(matches!(
            adr_algebra_of_basic_module(&[loop_module]),
            Err(AdrAlgebraError::InfiniteDimensionalAlgebra)
        ));
    }
}
