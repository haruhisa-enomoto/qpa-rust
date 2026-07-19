/// Translate QPA's cached `InjectiveEnvelope(M)` attribute.
///
/// The computation is the exact source expression
/// `DualOfModuleHomomorphism(ProjectiveCover(DualOfModule(M)))`. Reciprocal
/// dual-module caching makes the returned map start in the identical entered
/// module, including for the zero module.
///
/// # Errors
///
/// Propagates opposite-algebra dualization, projective-cover, and dual-map
/// construction failures.
pub fn injective_envelope<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_injective_envelope() {
        return Ok(cached);
    }
    let dual = dual_of_module(module)?;
    let dual_cover = projective_cover(&dual)?;
    let envelope = dual_of_module_homomorphism(&dual_cover)?;
    Ok(module.cache_injective_envelope(&envelope))
}
/// QPA's three-part `ExtOverAlgebra(M, N)` result.
///
/// The first component is the inclusion `Omega(M) -> P(M)`, and [`Self::basis`]
/// is QPA's ordered basis of `Ext^1(M, N)` represented inside
/// `Hom(Omega(M), N)`. When this basis is nonempty, [`Self::coefficients`]
/// translates the closure returned as QPA's third component.
#[derive(Clone, Debug)]
pub struct ExtOverAlgebraResult<F: Field> {
    syzygy_inclusion: PathAlgebraModuleHomomorphism<F>,
    basis: Vec<PathAlgebraModuleHomomorphism<F>>,
    coordinate_basis: Option<DenseMatrix<F>>,
}

/// QPA's three-part `ExtAlgebraGenerators(M, n)` result.
#[derive(Clone, Debug)]
pub struct ExtAlgebraGeneratorsResult<F: Field> {
    dimensions: Vec<usize>,
    generator_counts: Vec<usize>,
    generators: Vec<Vec<PathAlgebraModuleHomomorphism<F>>>,
}

impl<F: Field> ExtAlgebraGeneratorsResult<F> {
    /// Dimensions of `Ext^0(M,M), ..., Ext^n(M,M)` in degree order.
    #[must_use]
    pub fn dimensions(&self) -> &[usize] {
        &self.dimensions
    }

    /// Number of retained algebra generators in each degree.
    #[must_use]
    pub fn generator_counts(&self) -> &[usize] {
        &self.generator_counts
    }

    /// Retained module-map representatives in each degree.
    #[must_use]
    pub fn generators(&self) -> &[Vec<PathAlgebraModuleHomomorphism<F>>] {
        &self.generators
    }
}

/// Failure in QPA's bounded Ext-algebra generator construction.
#[derive(Debug, Error)]
pub enum ExtAlgebraGeneratorsError {
    /// Projective-cover, Ext, lifting, or exact matrix computation failed.
    #[error(transparent)]
    Enveloping(#[from] EnvelopingModuleError),
    /// Endomorphism-radical or module-map arithmetic failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// Exact product-coordinate or quotient-space matrix arithmetic failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// A source lifting unexpectedly returned QPA's `fail` value.
    #[error("an Ext-product projective lifting returned fail")]
    MissingProjectiveLifting,
    /// A source kernel-square construction unexpectedly returned `fail`.
    #[error("an Ext-product induced kernel morphism returned fail")]
    MissingInducedKernelMorphism,
    /// A nonzero Ext group must expose its coefficient closure.
    #[error("a nonzero Ext group has no coefficient function")]
    MissingExtCoefficients,
    /// Every represented endomorphism product belongs to `End(M)`.
    #[error("an endomorphism ideal product is outside the End(M) basis")]
    EndomorphismProductOutsideBasis,
    /// A nonempty map basis is required for a returned quotient representative.
    #[error("an Ext quotient representative has an empty map basis")]
    EmptyMapBasis,
}

/// Failure in QPA's partial or complete Iyama-generator construction.
#[derive(Debug, Error)]
pub enum IyamaGeneratorError {
    /// Endomorphism-radical, image, subrepresentation, or direct-sum failure.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// QPA's finite-field decomposition path failed.
    #[error(transparent)]
    Decomposition(#[from] ModuleDecompositionError),
    /// Detecting whether the Numerica coefficient field is finite failed.
    #[error(transparent)]
    Coefficient(#[from] CoefficientDomainError),
    /// A nonempty QPA direct-sum input unexpectedly returned `fail`.
    #[error("an internal nonempty Iyama-generator direct sum returned fail")]
    InternalDirectSumFailed,
    /// QPA indexes the failed multiplicity decomposition of a finite-field zero module.
    #[error("the finite-field Iyama generator has no multiplicity decomposition")]
    MissingMultiplicityDecomposition,
    /// A nonempty basic-summand list unexpectedly returned `fail` on direct sum.
    #[error("the finite-field Iyama basic-summand direct sum returned fail")]
    InternalBasicDirectSumFailed,
}

#[derive(Clone, Debug)]
enum ModuleTensorProductElementary<F: Field> {
    AlwaysAmbientZero,
    UnboundSizeMatrix,
    Presentation {
        cover: PathAlgebraModuleHomomorphism<F>,
        projections: Vec<PathAlgebraModuleHomomorphism<F>>,
        projective_vertices: Vec<usize>,
        closure_projection_count: usize,
    },
}

/// QPA's two-part `TensorProductOfModules(M, N)` result.
///
/// QPA returns a vector space together with an elementary-tensor closure. The
/// Rust value retains the same natural-quotient presentation: relation rows
/// inside the ambient free vector space, the selected standard-vector basis
/// lifts, and the projection into quotient coordinates.
#[derive(Clone, Debug)]
pub struct TensorProductOfModulesResult<F: Field> {
    dimension: usize,
    ambient_dimension: usize,
    relation_basis: DenseMatrix<F>,
    basis_lifts: DenseMatrix<F>,
    projection: DenseMatrix<F>,
    elementary: ModuleTensorProductElementary<F>,
}

impl<F: Field> TensorProductOfModulesResult<F> {
    /// Return the dimension of QPA's returned vector space.
    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    /// Return the dimension of the ambient free vector space before taking
    /// QPA's relation subspace.
    #[must_use]
    pub const fn ambient_dimension(&self) -> usize {
        self.ambient_dimension
    }

    /// Return GAP's semi-echelon basis of the relation subspace.
    #[must_use]
    pub const fn relation_basis(&self) -> &DenseMatrix<F> {
        &self.relation_basis
    }

    /// Return the ambient representatives of QPA's selected quotient basis.
    #[must_use]
    pub const fn basis_lifts(&self) -> &DenseMatrix<F> {
        &self.basis_lifts
    }

    /// Return the matrix sending ambient row vectors to quotient coordinates.
    #[must_use]
    pub const fn projection(&self) -> &DenseMatrix<F> {
        &self.projection
    }

    /// Return QPA's standard coordinate basis of the resulting vector space.
    #[must_use]
    pub fn basis_vectors(&self) -> Vec<Vec<F::Element>> {
        (0..self.dimension)
            .map(|coordinate| {
                let mut row = vec![self.projection.field().zero(); self.dimension];
                row[coordinate] = self.projection.field().one();
                row
            })
            .collect()
    }

    /// Apply the elementary-tensor closure returned by QPA.
    ///
    /// The easy zero branch deliberately ignores both arguments and returns
    /// QPA's one-coordinate ambient zero. The projective-presentation branch
    /// preserves the source bug where `size_matrix` is unassigned as a typed
    /// error.
    ///
    /// # Errors
    ///
    /// Propagates projective-cover preimages, opposite-element conversion,
    /// module action, and exact projection failures. It also reports QPA's
    /// unassigned-local failure in the corresponding source branch.
    pub fn elementary_tensor(
        &self,
        left: &PathAlgebraModuleElement<F>,
        right: &PathAlgebraModuleElement<F>,
    ) -> Result<Vec<F::Element>, EnvelopingModuleError> {
        match &self.elementary {
            ModuleTensorProductElementary::AlwaysAmbientZero => {
                Ok(vec![self.projection.field().zero()])
            }
            ModuleTensorProductElementary::UnboundSizeMatrix => {
                Err(EnvelopingModuleError::TensorProductUnboundSizeMatrix)
            }
            ModuleTensorProductElementary::Presentation {
                cover,
                projections,
                projective_vertices,
                closure_projection_count,
            } => {
                let projective_element = cover
                    .preimage_representative(left)?
                    .ok_or(EnvelopingModuleError::MissingTensorProjectiveCoverPreimage)?;
                let algebra = cover.range().right_acting_algebra();
                let mut acted_summands = Vec::with_capacity(*closure_projection_count);
                if *closure_projection_count > projections.len() {
                    return Err(
                        EnvelopingModuleError::TensorProductClosureProjectionOutOfRange {
                            index: projections.len() + 1,
                            count: projections.len(),
                        },
                    );
                }
                for (projection, &vertex) in projections
                    .iter()
                    .zip(projective_vertices)
                    .take(*closure_projection_count)
                {
                    let summand_element = projection.image(&projective_element)?;
                    let path_element =
                        element_in_indec_projective(algebra, &summand_element, vertex)?;
                    let opposite = opposite_projective_path_element(algebra, &path_element)?;
                    acted_summands.push(right.act(&opposite)?);
                }
                let mut ambient = Vec::with_capacity(self.ambient_dimension);
                for (index, &vertex) in projective_vertices.iter().enumerate() {
                    let acted = acted_summands.get(index).ok_or(
                        EnvelopingModuleError::TensorProductClosureActedIndexOutOfRange {
                            index: index + 1,
                            count: acted_summands.len(),
                        },
                    )?;
                    ambient.extend_from_slice(&acted.external_representation()[vertex - 1]);
                }
                let ambient_row = DenseMatrix::from_flat(
                    1,
                    self.ambient_dimension,
                    ambient,
                    self.projection.field().clone(),
                )?;
                Ok(ambient_row.multiply(&self.projection)?.row(0).to_vec())
            }
        }
    }
}

impl<F: Field> ExtOverAlgebraResult<F> {
    /// Return QPA's kernel inclusion `Omega(M) -> P(M)`.
    #[must_use]
    pub const fn syzygy_inclusion(&self) -> &PathAlgebraModuleHomomorphism<F> {
        &self.syzygy_inclusion
    }

    /// Return QPA's ordered representatives of the `Ext^1(M, N)` basis.
    #[must_use]
    pub fn basis(&self) -> &[PathAlgebraModuleHomomorphism<F>] {
        &self.basis
    }

    /// Return the dimension of `Ext^1(M, N)`.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.basis.len()
    }

    /// Apply QPA's returned coefficient closure to a map in
    /// `Hom(Omega(M), N)`.
    ///
    /// `Ok(None)` corresponds to QPA returning `[]` instead of a closure when
    /// the Ext group is zero.
    ///
    /// # Errors
    ///
    /// The entered map's raw physical vector must use the coefficient field
    /// and lie in the Hom space used by the original computation. Exact row
    /// solving failures propagate.
    pub fn coefficients(
        &self,
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
    ) -> Result<Option<Vec<F::Element>>, EnvelopingModuleError> {
        let Some(coordinate_basis) = &self.coordinate_basis else {
            return Ok(None);
        };
        if homomorphism.source().right_acting_algebra().field() != coordinate_basis.field() {
            return Err(EnvelopingModuleError::ExtCoefficientFieldMismatch);
        }
        let coefficients = coordinate_basis
            .solution_row(&homomorphism.flat_vertex_matrix_entries())?
            .ok_or(EnvelopingModuleError::ExtMapOutsideHomSpace)?;
        Ok(Some(
            coefficients[coefficients.len() - self.basis.len()..].to_vec(),
        ))
    }
}

/// Translate QPA's `ExtOverAlgebra(M, N)` construction of `Ext^1(M, N)`.
///
/// This follows the source method through the projective cover, its kernel
/// inclusion, `Hom(P(M), N)` and `Hom(Omega(M), N)`, the image of the induced
/// restriction map, and GAP's `NaturalHomomorphismBySubspace`. Basis maps are
/// the exact complement representatives selected from `Basis(V)` after the
/// image basis, preserving QPA's raw physical map-vector order.
///
/// # Errors
///
/// Requires both modules over one algebra and propagates projective-cover,
/// kernel, Hom, composition, exact vector-space, and checked-map failures.
pub fn ext_over_algebra<F: Field>(
    module: &PathAlgebraModule<F>,
    coefficients_in: &PathAlgebraModule<F>,
) -> Result<ExtOverAlgebraResult<F>, EnvelopingModuleError> {
    if module.right_acting_algebra() != coefficients_in.right_acting_algebra() {
        return Err(EnvelopingModuleError::ExtModulesOverDifferentAlgebras);
    }
    let cover = projective_cover(module)?;
    let syzygy_inclusion = kernel_inclusion(&cover)?;
    let syzygy = syzygy_inclusion.source();
    let projective_hom = hom_over_algebra(cover.source(), coefficients_in)?
        .ok_or(EnvelopingModuleError::InternalExtHomRejected)?;
    let syzygy_hom = hom_over_algebra(syzygy, coefficients_in)?
        .ok_or(EnvelopingModuleError::InternalExtHomRejected)?;
    let Some(first_homomorphism) = syzygy_hom.first() else {
        return Ok(ExtOverAlgebraResult {
            syzygy_inclusion,
            basis: Vec::new(),
            coordinate_basis: None,
        });
    };

    let field = module.right_acting_algebra().field().clone();
    let physical_length = first_homomorphism.flat_vertex_matrix_entries().len();
    let hom_rows = DenseMatrix::from_flat(
        syzygy_hom.len(),
        physical_length,
        syzygy_hom
            .iter()
            .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
            .collect(),
        field.clone(),
    )?;
    let hom_basis = hom_rows.semi_echelon_transformation().vectors;

    let image_maps = projective_hom
        .iter()
        .map(|homomorphism| syzygy_inclusion.compose(homomorphism))
        .collect::<Result<Vec<_>, _>>()?;
    let nonzero_image_maps = image_maps
        .iter()
        .filter(|homomorphism| !homomorphism.is_zero())
        .collect::<Vec<_>>();
    let image_basis = if nonzero_image_maps.is_empty() {
        DenseMatrix::from_flat(0, physical_length, Vec::new(), field.clone())?
    } else {
        DenseMatrix::from_flat(
            nonzero_image_maps.len(),
            physical_length,
            nonzero_image_maps
                .iter()
                .flat_map(|homomorphism| homomorphism.flat_vertex_matrix_entries())
                .collect(),
            field.clone(),
        )?
        .semi_echelon_transformation()
        .vectors
    };

    let quotient = quotient_complement_in_space(&hom_basis, &image_basis)?;
    if quotient.representatives.is_empty() {
        return Ok(ExtOverAlgebraResult {
            syzygy_inclusion,
            basis: Vec::new(),
            coordinate_basis: None,
        });
    }
    let basis = quotient
        .representatives
        .iter()
        .map(|row| module_map_from_physical_flat(syzygy, coefficients_in, row, &field))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(ExtOverAlgebraResult {
        syzygy_inclusion,
        basis,
        coordinate_basis: Some(quotient.coordinate_basis),
    })
}

/// Translate QPA's `ExtAlgebraGenerators(M, n)` method through degree `n`.
///
/// This is the source algorithm: construct successive projective covers,
/// compute the ordered Ext bases, lift every product of lower degrees through
/// those covers and their kernels, and retain GAP's natural-quotient basis
/// complements. Degree zero is `End(M)` modulo the ideal generated by QPA's
/// `rad(End(M)) * End(M) * rad(End(M))` product list.
///
/// # Errors
///
/// Propagates exact projective-cover, Ext, radical, lifting, kernel, map,
/// coordinate, and quotient-complement failures.
#[allow(clippy::too_many_lines)]
pub fn ext_algebra_generators<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    maximum_degree: usize,
) -> Result<ExtAlgebraGeneratorsResult<F>, ExtAlgebraGeneratorsError> {
    let field = module.right_acting_algebra().field().clone();
    let mut syzygy = module.clone();
    let mut projective_covers = Vec::with_capacity(maximum_degree);
    for _ in 0..maximum_degree {
        let cover = projective_cover(&syzygy)?;
        syzygy = kernel_of_module_homomorphism(&cover)?;
        projective_covers.push(cover);
    }

    let endomorphism_algebra = end_over_algebra(module)?;
    let radical_matrices = endomorphism_algebra.radical_basis_matrices()?;
    let radical_maps = radical_matrices
        .iter()
        .map(|matrix| from_end_m_to_hom_mm(module, matrix))
        .collect::<Result<Vec<_>, _>>()?;
    let mut ext_groups = Vec::with_capacity(maximum_degree + 1);
    ext_groups.push(ExtDegreeData::DegreeZero(radical_maps));
    for degree in 1..=maximum_degree {
        let result = if degree == 1 {
            ext_over_algebra(module, module)?
        } else {
            let previous_syzygy = kernel_of_module_homomorphism(&projective_covers[degree - 2])?;
            ext_over_algebra(&previous_syzygy, module)?
        };
        ext_groups.push(ExtDegreeData::Positive(result));
    }

    let mut dimensions = ext_groups
        .iter()
        .map(ExtDegreeData::basis_len)
        .collect::<Vec<_>>();
    let mut generator_counts = vec![0; maximum_degree + 1];
    let mut generators = vec![Vec::new(); maximum_degree + 1];

    for degree in 1..=maximum_degree {
        let mut product_coordinates = Vec::new();
        for left_degree in 0..=degree {
            let right_degree = degree - left_degree;
            if dimensions[degree] == 0
                || dimensions[left_degree] == 0
                || dimensions[right_degree] == 0
            {
                continue;
            }
            let mut induced = ext_groups[right_degree].basis().to_vec();
            for lift_degree in 1..=left_degree {
                let upper_cover = &projective_covers[right_degree + lift_degree - 1];
                let lower_cover = &projective_covers[lift_degree - 1];
                let mut liftings = Vec::with_capacity(induced.len());
                for map in &induced {
                    let target = upper_cover.compose(map)?;
                    liftings.push(
                        lifting_morphism_from_projective(lower_cover, &target)?
                            .ok_or(ExtAlgebraGeneratorsError::MissingProjectiveLifting)?,
                    );
                }
                induced = liftings
                    .iter()
                    .zip(&induced)
                    .map(|(lifting, map)| {
                        morphism_on_kernel(upper_cover, lower_cover, lifting, map)?
                            .ok_or(ExtAlgebraGeneratorsError::MissingInducedKernelMorphism)
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
            for left in ext_groups[left_degree].basis() {
                for right in &induced {
                    let product = right.compose(left)?;
                    let coordinates = ext_groups[degree]
                        .coefficients(&product)?
                        .ok_or(ExtAlgebraGeneratorsError::MissingExtCoefficients)?;
                    if coordinates.iter().any(|entry| !field.is_zero(entry)) {
                        product_coordinates.push(coordinates);
                    }
                }
            }
        }

        if dimensions[degree] == 0 {
            continue;
        }
        if product_coordinates.is_empty() {
            generators[degree] = ext_groups[degree].basis().to_vec();
        } else {
            let products = DenseMatrix::from_rows(product_coordinates, field.clone())?
                .canonical_row_space_basis();
            if products.nrows() < dimensions[degree] {
                let ambient = DenseMatrix::identity(dimensions[degree], field.clone());
                let complement = quotient_complement_in_space(&ambient, &products)?;
                generators[degree] = complement
                    .representatives
                    .iter()
                    .map(|coordinates| {
                        linear_combination_of_ext_maps(
                            ext_groups[degree].basis(),
                            coordinates,
                            &field,
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()?;
            }
        }
        generator_counts[degree] = generators[degree].len();
    }

    dimensions[0] = endomorphism_algebra.dimension();
    let endomorphism_basis = endomorphism_algebra.basis_matrices();
    let mut middle_products = Vec::new();
    for radical in &radical_matrices {
        for basis in endomorphism_basis {
            let product = endomorphism_algebra.multiply(basis, radical)?;
            if !product.is_zero() {
                middle_products.push(product);
            }
        }
    }
    let mut ideal_square_coordinates = Vec::new();
    for middle in &middle_products {
        for radical in &radical_matrices {
            let product = endomorphism_algebra.multiply(radical, middle)?;
            let coordinates = endomorphism_algebra
                .coefficients(&product)?
                .ok_or(ExtAlgebraGeneratorsError::EndomorphismProductOutsideBasis)?;
            ideal_square_coordinates.push(coordinates);
        }
    }
    let ideal_square = if ideal_square_coordinates.is_empty() {
        DenseMatrix::zero(0, dimensions[0], field.clone())
    } else {
        DenseMatrix::from_rows(ideal_square_coordinates, field.clone())?.canonical_row_space_basis()
    };
    let ambient = DenseMatrix::identity(dimensions[0], field.clone());
    let quotient = quotient_complement_in_space(&ambient, &ideal_square)?;
    generators[0] = quotient
        .representatives
        .iter()
        .map(|coordinates| {
            let matrix = endomorphism_algebra.matrix_from_coefficients(coordinates)?;
            Ok(from_end_m_to_hom_mm(module, &matrix)?)
        })
        .collect::<Result<Vec<_>, ExtAlgebraGeneratorsError>>()?;
    generator_counts[0] = generators[0].len();

    Ok(ExtAlgebraGeneratorsResult {
        dimensions,
        generator_counts,
        generators,
    })
}

enum ExtDegreeData<F: Field> {
    DegreeZero(Vec<PathAlgebraModuleHomomorphism<F>>),
    Positive(ExtOverAlgebraResult<F>),
}

impl<F: Field> ExtDegreeData<F> {
    fn basis(&self) -> &[PathAlgebraModuleHomomorphism<F>] {
        match self {
            Self::DegreeZero(basis) => basis,
            Self::Positive(result) => result.basis(),
        }
    }

    fn basis_len(&self) -> usize {
        self.basis().len()
    }

    fn coefficients(
        &self,
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
    ) -> Result<Option<Vec<F::Element>>, EnvelopingModuleError> {
        match self {
            Self::DegreeZero(_) => Ok(None),
            Self::Positive(result) => result.coefficients(homomorphism),
        }
    }
}

fn linear_combination_of_ext_maps<F: Field>(
    basis: &[PathAlgebraModuleHomomorphism<F>],
    coefficients: &[F::Element],
    field: &F,
) -> Result<PathAlgebraModuleHomomorphism<F>, ExtAlgebraGeneratorsError> {
    if basis.len() != coefficients.len() {
        return Err(MatrixError::TargetLength {
            actual: coefficients.len(),
            expected: basis.len(),
        }
        .into());
    }
    let first = basis
        .first()
        .ok_or(ExtAlgebraGeneratorsError::EmptyMapBasis)?;
    let mut combination = zero_module_mapping(first.source().clone(), first.range().clone())?;
    for (map, coefficient) in basis.iter().zip(coefficients) {
        if !field.is_zero(coefficient) {
            combination = combination.add(&map.scale(coefficient)?)?;
        }
    }
    Ok(combination)
}

/// Translate QPA's `PartialIyamaGenerator(M)` construction.
///
/// The source computes `rad End(M)`, converts its ordered basis matrices back
/// to endomorphisms, applies them in module-basis-major order, and generates
/// the resulting subrepresentation. The zero module is returned identically.
///
/// # Errors
///
/// Propagates the exact endomorphism-radical, matrix-to-map, element-image,
/// and subrepresentation algorithms.
pub fn partial_iyama_generator<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, IyamaGeneratorError> {
    if module.dimension() == 0 {
        return Ok(module.clone());
    }

    let basis = module.basis_elements();
    let endomorphism_algebra = end_over_algebra(module)?;
    let radical_matrices = endomorphism_algebra.radical_basis_matrices()?;
    let radical_maps = radical_matrices
        .iter()
        .map(|matrix| from_end_m_to_hom_mm(module, matrix))
        .collect::<Result<Vec<_>, _>>()?;
    let mut generators = Vec::with_capacity(
        basis
            .len()
            .checked_mul(radical_maps.len())
            .ok_or(PathAlgebraModuleError::DimensionOverflow)?,
    );
    for basis_element in &basis {
        for radical_map in &radical_maps {
            generators.push(radical_map.image(basis_element)?);
        }
    }
    Ok(subrepresentation(module.clone(), &generators)?)
}

/// Translate QPA's `IyamaGenerator(M)` radical-filtration direct sum.
///
/// QPA repeatedly applies [`partial_iyama_generator`], retaining every
/// nonzero term (or the entered zero module), then takes their direct sum. Over a finite
/// coefficient field it replaces that sum by one representative from every
/// isomorphism class returned by `DecomposeModuleWithMultiplicities`.
///
/// # Errors
///
/// Propagates partial-generator, direct-sum, finite-field detection, and
/// deterministic decomposition failures. The finite-field zero input exposes
/// QPA's subsequent indexing of a failed multiplicity decomposition as a
/// typed error.
pub fn iyama_generator<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, IyamaGeneratorError> {
    let mut filtration = Vec::new();
    let mut current = module.clone();
    loop {
        filtration.push(current.clone());
        current = partial_iyama_generator(&current)?;
        if current.dimension() == 0 {
            break;
        }
    }

    let mut generator = direct_sum_of_qpa_modules(&filtration)?
        .ok_or(IyamaGeneratorError::InternalDirectSumFailed)?;
    if module
        .right_acting_algebra()
        .field()
        .qpa_elements()?
        .is_some()
    {
        let decomposition = decompose_module_with_multiplicities(&generator)?
            .ok_or(IyamaGeneratorError::MissingMultiplicityDecomposition)?;
        generator = direct_sum_of_qpa_modules(decomposition.basic_summands())?
            .ok_or(IyamaGeneratorError::InternalBasicDirectSumFailed)?;
    }
    Ok(generator)
}

/// Translate QPA's cached `IsRigidModule(M)` property.
///
/// This is exactly the source predicate
/// `Length(ExtOverAlgebra(M, M)[2]) = 0`.
///
/// # Errors
///
/// Propagates every projective-cover, Hom, and quotient-space failure from
/// [`ext_over_algebra`].
pub fn is_rigid_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<bool, EnvelopingModuleError> {
    if let Some(cached) = module.cached_is_rigid_module() {
        return Ok(cached);
    }
    let rigid = ext_over_algebra(module, module)?.dimension() == 0;
    Ok(module.cache_is_rigid_module(rigid))
}

/// Translate QPA's cached `IsExceptionalModule(M)` property.
///
/// The source evaluates `IsRigidModule(M) and IsIndecomposableModule(M)` in
/// that order, including GAP's short-circuit. `Ok(None)` preserves the
/// infinite-field `TryNextMethod()` branch of the latter property.
///
/// # Errors
///
/// Propagates rigidity/Ext construction and indecomposable-decomposition
/// failures.
pub fn is_exceptional_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Option<bool>, ExceptionalModuleError> {
    if let Some(cached) = module.cached_is_exceptional_module() {
        return Ok(Some(cached));
    }
    if !is_rigid_module(module)? {
        return Ok(Some(module.cache_is_exceptional_module(false)));
    }
    let Some(indecomposable) = is_indecomposable_module(module)? else {
        return Ok(None);
    };
    Ok(Some(module.cache_is_exceptional_module(indecomposable)))
}

/// Translate QPA's estimate `ComplexityOfModule(M, n)`.
///
/// If finite global dimension is already known, QPA returns zero without
/// inspecting `n`. Otherwise it requires an admissible quotient, asks QPA's
/// lazy [`crate::ProjectiveResolution`] for its object in degree `n`, and returns
/// `Int(Round(log10(dim(P_n))/log10(n) + 1.0))`. Degrees below `-1` address
/// the zero tail of QPA's complex and therefore also return zero.
///
/// # Errors
///
/// Rejects a nonadmissible or nonquotient domain and makes QPA's undefined
/// nonzero estimates at degrees `-1`, `0`, and `1` explicit. Lazy
/// projective-resolution failures are propagated.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss
)]
pub fn complexity_of_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    degree: &Integer,
) -> Result<usize, EnvelopingModuleError> {
    let algebra = module.right_acting_algebra();
    if has_known_finite_global_dimension(algebra) {
        return Ok(0);
    }
    require_complexity_admissible_quotient(algebra)?;
    let projective_dimension = projective_resolution_object_dimension(module, degree)?;
    if projective_dimension == 0 {
        return Ok(0);
    }
    let degree_index = usize::try_from(degree.clone()).map_err(|_| {
        EnvelopingModuleError::UndefinedComplexityEstimate {
            degree: degree.clone(),
        }
    })?;
    if degree_index < 2 {
        return Err(EnvelopingModuleError::UndefinedComplexityEstimate {
            degree: degree.clone(),
        });
    }
    let estimate = (projective_dimension as f64).log10() / (degree_index as f64).log10();
    Ok((estimate + 1.0).round() as usize)
}

/// Translate QPA's maximum `ComplexityOfAlgebra(A, n)` estimate over all
/// vertex simple modules.
///
/// # Errors
///
/// Preserves the same cached-finite-global-dimension and admissibility order
/// as [`complexity_of_module`], then propagates simple-module and individual
/// complexity failures.
pub fn complexity_of_algebra<F: QpaCoefficientField>(
    algebra: &QuiverAlgebra<F>,
    degree: &Integer,
) -> Result<usize, EnvelopingModuleError> {
    if has_known_finite_global_dimension(algebra) {
        return Ok(0);
    }
    require_complexity_admissible_quotient(algebra)?;
    let mut estimates = simple_modules(algebra)?
        .into_iter()
        .map(|simple| complexity_of_module(&simple, degree));
    let mut maximum = estimates
        .next()
        .ok_or(EnvelopingModuleError::ComplexityAlgebraHasNoSimpleModules)??;
    for estimate in estimates {
        maximum = maximum.max(estimate?);
    }
    Ok(maximum)
}

fn has_known_finite_global_dimension<F: Field>(algebra: &QuiverAlgebra<F>) -> bool {
    match algebra {
        QuiverAlgebra::Path(_) => true,
        QuiverAlgebra::Quotient(quotient) => matches!(
            quotient.cached_global_dimension(),
            Some(AlgebraGlobalDimension::Finite(_))
        ),
    }
}

fn require_complexity_admissible_quotient<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<(), EnvelopingModuleError> {
    match algebra {
        QuiverAlgebra::Quotient(quotient) if quotient.is_admissible()? => Ok(()),
        QuiverAlgebra::Path(_) | QuiverAlgebra::Quotient(_) => {
            Err(EnvelopingModuleError::ComplexityRequiresAdmissibleQuotient)
        }
    }
}

fn projective_resolution_object_dimension<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    degree: &Integer,
) -> Result<usize, EnvelopingModuleError> {
    Ok(projective_resolution(module)?.object(degree)?.dimension())
}
