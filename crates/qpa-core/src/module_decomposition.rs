//! QPA's finite-dimensional module-decomposition foundations.
//!
//! This module begins the source-faithful translation of
//! `lib/moduledecomp.gi`: row-space complements and the idempotent algorithms
//! used by the deterministic module decomposition methods.

use std::cmp::Ordering;

use numerica::domains::{Field, integer::Integer};
use rand::{Rng, SeedableRng, rngs::StdRng};
use thiserror::Error;

use crate::{
    CoefficientDomainError, CommonDirectSummandResult, DenseMatrix, HopfError, MatrixError,
    PathAlgebraModule, PathAlgebraModuleError, PathAlgebraModuleHomomorphism, QpaCoefficientField,
    QpaFieldCharacteristic, QuiverAlgebra, QuotientError, StructureConstantsAlgebra,
    StructureConstantsAlgebraElement, StructureConstantsError, StructureConstantsQuotient,
    cokernel_projection, common_direct_summand, direct_sum_of_qpa_modules, end_over_algebra,
    from_end_m_to_hom_mm, from_hom_mm_to_end_m, hom_over_algebra,
    idempotent::{factor_finite_field_polynomial, polynomial_product},
    identity_module_mapping, image_inclusion, image_of_module_homomorphism, image_projection,
    is_simple_qpa_module, isomorphism_of_modules, kernel_inclusion, kernel_of_module_homomorphism,
    maximal_common_direct_summand,
    path_module::ModuleDecompositionMultiplicitiesCache,
    pull_back, right_module_hom_over_algebra, socle_of_module, socle_of_module_inclusion,
    structure_constants::{StructureConstantsSubalgebra, add_rows, span_rows, subtract_rows},
    subrepresentation_inclusion, top_of_module, top_of_module_homomorphism,
};

/// Translate QPA's `ComplementInFullRowSpace(V)` operation.
///
/// `basis` contains `BasisVectors(Basis(V))` as rows. The source triangulates
/// those rows, removes their pivot positions from `[1..n]`, and returns the
/// corresponding standard basis vectors in increasing coordinate order.
///
/// # Errors
///
/// QPA errors on the trivial space. Exact shaped-matrix failures propagate.
pub fn complement_in_full_row_space<F: Field>(
    basis: &DenseMatrix<F>,
) -> Result<DenseMatrix<F>, ModuleDecompositionError> {
    let triangulated = basis.semi_echelon_transformation().vectors;
    if triangulated.nrows() == 0 {
        return Err(ModuleDecompositionError::TrivialRowSpace);
    }
    if triangulated.nrows() == triangulated.ncols() {
        return Ok(DenseMatrix::zero(
            0,
            triangulated.ncols(),
            triangulated.field().clone(),
        ));
    }

    let one = triangulated.field().one();
    let mut needed = vec![true; triangulated.ncols()];
    for row in 0..triangulated.nrows() {
        let pivot = triangulated
            .row(row)
            .iter()
            .position(|entry| entry == &one)
            .ok_or(ModuleDecompositionError::TriangulatedRowWithoutPivot { row: row + 1 })?;
        needed[pivot] = false;
    }

    let complement_dimension = needed.iter().filter(|&&is_needed| is_needed).count();
    let mut entries = Vec::with_capacity(complement_dimension * triangulated.ncols());
    for (coordinate, &is_needed) in needed.iter().enumerate() {
        if !is_needed {
            continue;
        }
        for column in 0..triangulated.ncols() {
            entries.push(if column == coordinate {
                triangulated.field().one()
            } else {
                triangulated.field().zero()
            });
        }
    }
    Ok(DenseMatrix::from_flat(
        complement_dimension,
        triangulated.ncols(),
        entries,
        triangulated.field().clone(),
    )?)
}

/// A linear algebra homomorphism between finite-dimensional
/// structure-constant algebras.
///
/// Rows of `basis_images` are the range coordinates of the ordered source
/// basis. This is the supported exact analogue of GAP's
/// `IsAlgebraGeneralMapping` inputs in `moduledecomp.gi`.
#[derive(Clone, Debug)]
pub struct StructureConstantsAlgebraHomomorphism<F: Field> {
    source: StructureConstantsAlgebra<F>,
    range: StructureConstantsAlgebra<F>,
    basis_images: DenseMatrix<F>,
}

impl<F: Field> StructureConstantsAlgebraHomomorphism<F> {
    /// Construct and validate an algebra homomorphism from ordered source-basis
    /// images.
    ///
    /// # Errors
    ///
    /// Requires a common field, the exact `dim(source) × dim(range)` shape,
    /// and preservation of every ordered basis product.
    pub fn by_basis_images(
        source: StructureConstantsAlgebra<F>,
        range: StructureConstantsAlgebra<F>,
        basis_images: DenseMatrix<F>,
    ) -> Result<Self, ModuleDecompositionError> {
        if source.field() != range.field() || basis_images.field() != source.field() {
            return Err(ModuleDecompositionError::MappingFieldsDiffer);
        }
        let expected = (source.dimension(), range.dimension());
        let actual = (basis_images.nrows(), basis_images.ncols());
        if actual != expected {
            return Err(ModuleDecompositionError::MappingShape { actual, expected });
        }
        let mapping = Self {
            source,
            range,
            basis_images,
        };
        for left in 0..mapping.source.dimension() {
            for right in 0..mapping.source.dimension() {
                let source_product = mapping.source.basis_product(left, right);
                let image_product = mapping.image_coordinates(source_product)?;
                let range_product = mapping.range.multiply_coordinates(
                    mapping.basis_images.row(left),
                    mapping.basis_images.row(right),
                )?;
                if image_product != range_product {
                    return Err(ModuleDecompositionError::MappingNotMultiplicative {
                        left: left + 1,
                        right: right + 1,
                    });
                }
            }
        }
        Ok(mapping)
    }

    /// Return the source algebra.
    #[must_use]
    pub const fn source(&self) -> &StructureConstantsAlgebra<F> {
        &self.source
    }

    /// Return the range algebra.
    #[must_use]
    pub const fn range(&self) -> &StructureConstantsAlgebra<F> {
        &self.range
    }

    /// Return the ordered source-basis image matrix.
    #[must_use]
    pub const fn basis_images(&self) -> &DenseMatrix<F> {
        &self.basis_images
    }

    /// Apply the homomorphism to source coordinates.
    ///
    /// # Errors
    ///
    /// Rejects a row whose width differs from the source dimension.
    pub fn image_coordinates(
        &self,
        coordinates: &[F::Element],
    ) -> Result<Vec<F::Element>, ModuleDecompositionError> {
        if coordinates.len() != self.source.dimension() {
            return Err(ModuleDecompositionError::MappingSourceLength {
                actual: coordinates.len(),
                expected: self.source.dimension(),
            });
        }
        Ok(multiply_row_matrix(coordinates, &self.basis_images))
    }

    /// Return GAP's linear preimage representative when the range element is
    /// in the image.
    ///
    /// # Errors
    ///
    /// Rejects a row whose width differs from the range dimension and
    /// propagates exact row solving.
    pub fn preimage_representative(
        &self,
        coordinates: &[F::Element],
    ) -> Result<Option<Vec<F::Element>>, ModuleDecompositionError> {
        if coordinates.len() != self.range.dimension() {
            return Err(ModuleDecompositionError::MappingRangeLength {
                actual: coordinates.len(),
                expected: self.range.dimension(),
            });
        }
        Ok(self.basis_images.solution_row(coordinates)?)
    }
}

/// Translate QPA's `LiftIdempotentsForDecomposition(map, ids)` operation.
///
/// Preimages are adjusted successively to be orthogonal to the accumulated
/// lift, corrected by the literal `3 e^2 - 2 e^3` iteration, and returned in
/// GAP sorted-set order.
///
/// # Errors
///
/// QPA indexes the first entered idempotent, so an empty list is rejected.
/// Missing preimages, exact products, coefficient ordering, and a correction
/// that exceeds the nilpotent-kernel bound are reported explicitly.
pub fn lift_idempotents_for_decomposition<F: QpaCoefficientField>(
    mapping: &StructureConstantsAlgebraHomomorphism<F>,
    idempotents: &[Vec<F::Element>],
) -> Result<Vec<Vec<F::Element>>, ModuleDecompositionError> {
    if idempotents.is_empty() {
        return Err(ModuleDecompositionError::EmptyIdempotentSet);
    }
    let mut lifts = idempotents
        .iter()
        .enumerate()
        .map(|(position, idempotent)| {
            mapping.preimage_representative(idempotent)?.ok_or(
                ModuleDecompositionError::MissingPreimage {
                    position: position + 1,
                },
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let source = mapping.source();
    let mut accumulated = vec![source.field().zero(); source.dimension()];
    let two = source.field().nth(Integer::from(2));
    let three = source.field().nth(Integer::from(3));

    for lift in &mut lifts {
        let accumulated_left = source.multiply_coordinates(&accumulated, lift)?;
        let accumulated_right = source.multiply_coordinates(lift, &accumulated)?;
        let sandwich = source.multiply_coordinates(&accumulated_left, &accumulated)?;
        let mut adjusted = subtract_rows(source.field(), lift, &accumulated_left);
        adjusted = subtract_rows(source.field(), &adjusted, &accumulated_right);
        adjusted = add_rows(source.field(), &adjusted, &sandwich);

        let mut iterations = 0usize;
        loop {
            let square = source.multiply_coordinates(&adjusted, &adjusted)?;
            if square == adjusted {
                break;
            }
            let cube = source.multiply_coordinates(&square, &adjusted)?;
            adjusted = square
                .iter()
                .zip(cube)
                .map(|(square, cube)| {
                    source.field().sub(
                        &source.field().mul(&three, square),
                        &source.field().mul(&two, &cube),
                    )
                })
                .collect();
            iterations = iterations
                .checked_add(1)
                .ok_or(ModuleDecompositionError::IdempotentCorrectionDidNotConverge)?;
            if iterations > source.dimension().saturating_add(1) {
                return Err(ModuleDecompositionError::IdempotentCorrectionDidNotConverge);
            }
        }
        *lift = adjusted;
        accumulated = add_rows(source.field(), &accumulated, lift);
    }

    gap_sort_and_deduplicate(source.field(), &mut lifts)?;
    Ok(lifts)
}

/// Translate QPA's Anderson--Fuller `LiftIdempotent(f, v)` operation.
///
/// # Errors
///
/// The entered range element must be idempotent and lie in the mapping image.
/// A missing source identity or nonnilpotent defect is reported explicitly.
pub fn lift_idempotent<F: Field>(
    mapping: &StructureConstantsAlgebraHomomorphism<F>,
    idempotent: &[F::Element],
) -> Result<Vec<F::Element>, ModuleDecompositionError> {
    let representative = mapping
        .preimage_representative(idempotent)?
        .ok_or(ModuleDecompositionError::EnteredElementOutsideImage)?;
    if mapping
        .range()
        .multiply_coordinates(idempotent, idempotent)?
        != idempotent
    {
        return Err(ModuleDecompositionError::EnteredElementNotIdempotent);
    }
    let source = mapping.source();
    let square = source.multiply_coordinates(&representative, &representative)?;
    let defect = subtract_rows(source.field(), &representative, &square);
    let mut power = defect.clone();
    let mut nilindex = 1usize;
    while !is_zero_row(source.field(), &power) {
        nilindex = nilindex
            .checked_add(1)
            .ok_or(ModuleDecompositionError::NilpotentSeriesOverflow)?;
        power = source.multiply_coordinates(&defect, &power)?;
        if nilindex > source.dimension().saturating_add(1) {
            return Err(ModuleDecompositionError::MappingKernelNotNilpotent);
        }
    }
    if nilindex == 1 {
        return Ok(representative);
    }
    let identity = source
        .find_multiplicative_identity()?
        .ok_or(ModuleDecompositionError::MissingSourceIdentity)?;
    let mut polynomial = vec![source.field().zero(); source.dimension()];
    for exponent in 0..nilindex {
        let coefficient = binomial_integer(nilindex, exponent + 1);
        let coefficient = source.field().nth(if exponent.is_multiple_of(2) {
            coefficient
        } else {
            -coefficient
        });
        let term = algebra_power(source, &representative, exponent, identity.coordinates())?;
        for (entry, term) in polynomial.iter_mut().zip(term) {
            source.field().add_mul_assign(entry, &coefficient, &term);
        }
    }
    let representative_power =
        algebra_power(source, &representative, nilindex, identity.coordinates())?;
    let polynomial_power = algebra_power(source, &polynomial, nilindex, identity.coordinates())?;
    Ok(source.multiply_coordinates(&representative_power, &polynomial_power)?)
}

/// Translate QPA's Lam-series `LiftTwoOrthogonalIdempotents(f, v, w)`.
///
/// # Errors
///
/// Rejects nonorthogonal entered idempotents and propagates idempotent lifting,
/// identity, exact-product, and nilpotence failures.
pub fn lift_two_orthogonal_idempotents<F: Field>(
    mapping: &StructureConstantsAlgebraHomomorphism<F>,
    first: &[F::Element],
    second_range: &[F::Element],
) -> Result<[Vec<F::Element>; 2], ModuleDecompositionError> {
    let first_image = mapping.image_coordinates(first)?;
    let left = mapping
        .range()
        .multiply_coordinates(&first_image, second_range)?;
    let right = mapping
        .range()
        .multiply_coordinates(second_range, &first_image)?;
    if !is_zero_row(mapping.range().field(), &left) || !is_zero_row(mapping.range().field(), &right)
    {
        return Err(ModuleDecompositionError::EnteredIdempotentsNotOrthogonal);
    }
    let second = lift_idempotent(mapping, second_range)?;
    let source = mapping.source();
    let product = source.multiply_coordinates(&second, first)?;
    let identity = source
        .find_multiplicative_identity()?
        .ok_or(ModuleDecompositionError::MissingSourceIdentity)?;
    let mut series = identity.coordinates().to_vec();
    let mut power = product.clone();
    let mut steps = 0usize;
    while !is_zero_row(source.field(), &power) {
        series = add_rows(source.field(), &series, &power);
        power = source.multiply_coordinates(&product, &power)?;
        steps = steps
            .checked_add(1)
            .ok_or(ModuleDecompositionError::NilpotentSeriesOverflow)?;
        if steps > source.dimension() {
            return Err(ModuleDecompositionError::MappingKernelNotNilpotent);
        }
    }
    let one_minus_first = subtract_rows(source.field(), identity.coordinates(), first);
    let one_minus_product = subtract_rows(source.field(), identity.coordinates(), &product);
    let corrected = source.multiply_coordinates(&one_minus_first, &series)?;
    let corrected = source.multiply_coordinates(&corrected, &second)?;
    let corrected = source.multiply_coordinates(&corrected, &one_minus_product)?;
    Ok([first.to_vec(), corrected])
}

/// Translate QPA's `LiftingCompleteSetOfOrthogonalIdempotents(map, ids)` for
/// a general finite-dimensional structure-constant algebra mapping.
///
/// The source validates idempotence, ordered orthogonality, and completeness,
/// then applies `LiftIdempotent` once and `LiftTwoOrthogonalIdempotents`
/// successively against the sum already lifted.
///
/// # Errors
///
/// Rejects an empty, non-idempotent, nonorthogonal, or incomplete entered set
/// and propagates exact mapping and lifting failures.
pub fn lifting_complete_set_of_orthogonal_idempotents<F: Field>(
    mapping: &StructureConstantsAlgebraHomomorphism<F>,
    idempotents: &[Vec<F::Element>],
) -> Result<Vec<Vec<F::Element>>, ModuleDecompositionError> {
    if idempotents.is_empty() {
        return Err(ModuleDecompositionError::EmptyIdempotentSet);
    }
    let range = mapping.range();
    let range_identity = range
        .find_multiplicative_identity()?
        .ok_or(ModuleDecompositionError::MissingRangeIdentity)?;
    let mut sum = vec![range.field().zero(); range.dimension()];
    for (left_position, left) in idempotents.iter().enumerate() {
        if range.multiply_coordinates(left, left)?.as_slice() != left {
            return Err(ModuleDecompositionError::EnteredElementNotIdempotent);
        }
        for right in idempotents.iter().skip(left_position + 1) {
            let left_right = range.multiply_coordinates(left, right)?;
            let right_left = range.multiply_coordinates(right, left)?;
            if !is_zero_row(range.field(), &left_right) || !is_zero_row(range.field(), &right_left)
            {
                return Err(ModuleDecompositionError::EnteredIdempotentsNotOrthogonal);
            }
        }
        sum = add_rows(range.field(), &sum, left);
    }
    if sum != range_identity.coordinates() {
        return Err(ModuleDecompositionError::EnteredIdempotentsNotComplete);
    }

    let mut lifted = vec![lift_idempotent(mapping, &idempotents[0])?];
    for idempotent in &idempotents[1..] {
        let accumulated = lifted.iter().fold(
            vec![mapping.source().field().zero(); mapping.source().dimension()],
            |sum, idempotent| add_rows(mapping.source().field(), &sum, idempotent),
        );
        let pair = lift_two_orthogonal_idempotents(mapping, &accumulated, idempotent)?;
        lifted.push(pair[1].clone());
    }
    Ok(lifted)
}

impl<F: QpaCoefficientField> StructureConstantsQuotient<F> {
    /// Return this natural quotient as a validated general algebra mapping.
    ///
    /// # Errors
    ///
    /// Propagates projection of the ordered original basis, exact matrix
    /// construction, and multiplicativity validation.
    pub fn natural_homomorphism(
        &self,
    ) -> Result<StructureConstantsAlgebraHomomorphism<F>, ModuleDecompositionError> {
        let mut entries = Vec::with_capacity(
            self.original()
                .dimension()
                .saturating_mul(self.algebra().dimension()),
        );
        for basis_index in 0..self.original().dimension() {
            let mut basis = vec![self.original().field().zero(); self.original().dimension()];
            basis[basis_index] = self.original().field().one();
            entries.extend(self.project_coordinates(&basis)?);
        }
        StructureConstantsAlgebraHomomorphism::by_basis_images(
            self.original().clone(),
            self.algebra().clone(),
            DenseMatrix::from_flat(
                self.original().dimension(),
                self.algebra().dimension(),
                entries,
                self.original().field().clone(),
            )?,
        )
    }

    /// Translate QPA's `LiftIdempotentsForDecomposition(map, ids)` for a
    /// natural algebra quotient.
    ///
    /// The stored GAP quotient representatives are adjusted successively to
    /// be orthogonal to the accumulated lift, corrected by the literal
    /// `3 e^2 - 2 e^3` iteration, and returned in GAP sorted-set order.
    ///
    /// # Errors
    ///
    /// QPA indexes the first entered idempotent, so an empty list is rejected.
    /// Invalid coordinate rows, exact products, coefficient ordering, and a
    /// correction that does not terminate within the nilpotent-kernel bound
    /// are reported explicitly.
    pub fn lift_idempotents_for_decomposition(
        &self,
        idempotents: &[Vec<F::Element>],
    ) -> Result<Vec<Vec<F::Element>>, ModuleDecompositionError> {
        lift_idempotents_for_decomposition(&self.natural_homomorphism()?, idempotents)
    }
}

impl<F: QpaCoefficientField> StructureConstantsAlgebra<F> {
    /// Translate QPA's `IdempotentsForDecomposition(A)` operation.
    ///
    /// The source passes to `A/rad(A)`, splits it by central idempotents,
    /// computes primitive idempotents inside each simple central component,
    /// lifts their component coordinates to the semisimple quotient, and then
    /// invokes `LiftIdempotentsForDecomposition`.
    ///
    /// # Errors
    ///
    /// Propagates radical/quotient construction, the QPA primitive-idempotent
    /// algorithm, exact coordinate lifting, and decomposition lifting errors.
    pub fn idempotents_for_decomposition(
        &self,
    ) -> Result<Vec<StructureConstantsAlgebraElement<F>>, ModuleDecompositionError> {
        let radical = self.radical()?;
        let quotient = self.quotient_by_ideal_basis(radical.basis().clone())?;
        let central_idempotents = quotient.algebra().central_idempotents_split_prime()?;
        let mut primitive_quotient_coordinates = Vec::new();

        for central_idempotent in central_idempotents {
            let generators = (0..quotient.algebra().dimension())
                .map(|basis_index| {
                    let mut basis = vec![self.field().zero(); quotient.algebra().dimension()];
                    basis[basis_index] = self.field().one();
                    quotient
                        .algebra()
                        .multiply_coordinates(&basis, &central_idempotent)
                })
                .collect::<Result<Vec<_>, _>>()?;
            let basis = span_rows(self.field(), quotient.algebra().dimension(), generators)?;
            let component =
                StructureConstantsSubalgebra::from_basis(quotient.algebra().clone(), basis)?;
            for primitive in component.algebra().primitive_idempotents()? {
                primitive_quotient_coordinates
                    .push(component.lift_coordinates(primitive.coordinates())?);
            }
        }

        quotient
            .lift_idempotents_for_decomposition(&primitive_quotient_coordinates)?
            .into_iter()
            .map(|coordinates| self.element(coordinates).map_err(Into::into))
            .collect()
    }
}

/// Translate QPA's cached `DecomposeModule(M)` attribute.
///
/// For a module built by `DirectSumOfQPAModules`, QPA recursively decomposes
/// the retained summands and flattens the result. Otherwise it forms
/// `EndOverAlgebra(M)`, reverses multiplication through the source's
/// transpose construction, lifts a complete primitive-idempotent set, and
/// returns the images of the corresponding endomorphisms.
///
/// # Errors
///
/// Propagates endomorphism-algebra construction, primitive-idempotent
/// computation, conversion back to module maps, and image construction.
pub fn decompose_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModule<F>>, ModuleDecompositionError> {
    if let Some(cached) = module.cached_decomposition() {
        return Ok(cached);
    }
    if let Some(summands) = module.direct_sum_summands() {
        if module.dimension() == 0 {
            return Ok(module.cache_decomposition(&[]));
        }
        let mut decomposition = Vec::new();
        for summand in summands {
            decomposition.extend(decompose_module(summand)?);
        }
        return Ok(module.cache_decomposition(&decomposition));
    }

    let endomorphisms = decomposition_endomorphisms(module)?;
    let decomposition = endomorphisms
        .iter()
        .map(image_of_module_homomorphism)
        .collect::<Result<Vec<_>, _>>()?;
    Ok(module.cache_decomposition(&decomposition))
}

/// Translate QPA's `IsIndecomposableModule(M)` method.
///
/// The zero module is false. Over a finite field QPA decides by whether
/// `DecomposeModule(M)` has one summand. Over an infinite field it returns
/// true when either the top or socle is simple, and otherwise invokes
/// `TryNextMethod`; that final branch is represented by `Ok(None)`.
///
/// # Errors
///
/// Propagates coefficient-domain enumeration, decomposition, top, and socle
/// failures.
pub fn is_indecomposable_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Option<bool>, ModuleDecompositionError> {
    if let Some(cached) = module.cached_is_indecomposable_module() {
        return Ok(Some(cached));
    }
    if module.dimension() == 0 {
        return Ok(Some(module.cache_is_indecomposable_module(false)));
    }
    if module
        .right_acting_algebra()
        .field()
        .qpa_elements()?
        .is_some()
    {
        let indecomposable = decompose_module(module)?.len() == 1;
        return Ok(Some(module.cache_is_indecomposable_module(indecomposable)));
    }
    if is_simple_qpa_module(&top_of_module(module.clone())?)
        || is_simple_qpa_module(&socle_of_module(module.clone())?)
    {
        return Ok(Some(module.cache_is_indecomposable_module(true)));
    }
    Ok(None)
}

/// Translate QPA's cached `DecomposeModuleWithInclusions(M)` attribute.
///
/// QPA repeats the deterministic endomorphism/idempotent construction used by
/// `DecomposeModule` and returns the canonical image inclusion for each
/// idempotent endomorphism.
///
/// # Errors
///
/// Propagates endomorphism-algebra construction, primitive-idempotent
/// computation, conversion back to module maps, and image-inclusion
/// construction.
pub fn decompose_module_with_inclusions<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, ModuleDecompositionError> {
    if let Some(cached) = module.cached_decomposition_inclusions() {
        return Ok(cached);
    }
    let inclusions = decomposition_endomorphisms(module)?
        .iter()
        .map(image_inclusion)
        .collect::<Result<Vec<_>, _>>()?;
    Ok(module.cache_decomposition_inclusions(&inclusions))
}

/// QPA's nonzero `DecomposeModuleWithMultiplicities(M)` result.
#[derive(Clone, Debug)]
pub struct ModuleDecompositionWithMultiplicities<F: Field> {
    basic_summands: Vec<PathAlgebraModule<F>>,
    multiplicities: Vec<usize>,
}

/// QPA's heterogeneous `AllSubmodulesOfModule` return value.
///
/// The zero and simple branches return one flat list, while the general
/// branch returns a list indexed by vector-space dimension.
#[derive(Clone, Debug)]
pub enum AllSubmodulesOfModuleResult<F: Field> {
    /// QPA's zero- or simple-module early return.
    Flat(Vec<PathAlgebraModuleHomomorphism<F>>),
    /// General result, with level `d` containing the `d`-dimensional
    /// submodules (and therefore stored at zero-based index `d`).
    ByDimension(Vec<Vec<PathAlgebraModuleHomomorphism<F>>>),
}

impl<F: Field> AllSubmodulesOfModuleResult<F> {
    /// Return the early-return list, if QPA selected that branch.
    #[must_use]
    pub fn flat(&self) -> Option<&[PathAlgebraModuleHomomorphism<F>]> {
        match self {
            Self::Flat(maps) => Some(maps),
            Self::ByDimension(_) => None,
        }
    }

    /// Return the dimension-indexed levels, if QPA selected its general
    /// branch.
    #[must_use]
    pub fn by_dimension(&self) -> Option<&[Vec<PathAlgebraModuleHomomorphism<F>>]> {
        match self {
            Self::Flat(_) => None,
            Self::ByDimension(levels) => Some(levels),
        }
    }

    /// Consume the source-shaped result and flatten its maps in QPA order.
    #[must_use]
    pub fn into_flattened(self) -> Vec<PathAlgebraModuleHomomorphism<F>> {
        match self {
            Self::Flat(maps) => maps,
            Self::ByDimension(levels) => levels.into_iter().flatten().collect(),
        }
    }
}

impl<F: Field> ModuleDecompositionWithMultiplicities<F> {
    /// Return the first summand retained from each isomorphism class.
    #[must_use]
    pub fn basic_summands(&self) -> &[PathAlgebraModule<F>] {
        &self.basic_summands
    }

    /// Return multiplicities in the corresponding basic-summand order.
    #[must_use]
    pub fn multiplicities(&self) -> &[usize] {
        &self.multiplicities
    }
}

/// Translate QPA's cached `DecomposeModuleWithMultiplicities(M)` attribute.
///
/// QPA returns `fail` for the zero module, represented by `Ok(None)`. For a
/// nonzero module it walks `DecomposeModule(M)` in order, retaining the first
/// member of each isomorphism class and detecting equal classes through the
/// source's `CommonDirectSummand` test after a dimension-vector check.
///
/// # Errors
///
/// Propagates deterministic decomposition and common-direct-summand search.
pub fn decompose_module_with_multiplicities<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Option<ModuleDecompositionWithMultiplicities<F>>, ModuleDecompositionError> {
    if let Some(cached) = module.cached_decomposition_multiplicities() {
        return Ok(match cached {
            ModuleDecompositionMultiplicitiesCache::Fail => None,
            ModuleDecompositionMultiplicitiesCache::Value {
                basic_summands,
                multiplicities,
            } => Some(ModuleDecompositionWithMultiplicities {
                basic_summands,
                multiplicities,
            }),
        });
    }
    if module.dimension() == 0 {
        let _ =
            module.cache_decomposition_multiplicities(ModuleDecompositionMultiplicitiesCache::Fail);
        return Ok(None);
    }

    let decomposition = decompose_module(module)?;
    if decomposition.is_empty() {
        return Err(ModuleDecompositionError::EmptyNonzeroDecomposition);
    }
    let mut current = 0usize;
    let mut rest = (1..decomposition.len()).collect::<Vec<_>>();
    let mut basic_summands = vec![decomposition[0].clone()];
    let mut multiplicities = vec![1usize];

    while !rest.is_empty() {
        let positions = rest.clone();
        for position in positions {
            if decomposition[current].dimension_vector()
                == decomposition[position].dimension_vector()
                && matches!(
                    common_direct_summand(&decomposition[current], &decomposition[position])?,
                    CommonDirectSummandResult::Found(_)
                )
            {
                let Some(multiplicity) = multiplicities.last_mut() else {
                    return Err(ModuleDecompositionError::EmptyNonzeroDecomposition);
                };
                *multiplicity += 1;
                rest.retain(|&remaining| remaining != position);
            }
        }
        if let Some(&next) = rest.first() {
            current = next;
            basic_summands.push(decomposition[current].clone());
            rest.remove(0);
            multiplicities.push(1);
        }
    }

    match module.cache_decomposition_multiplicities(ModuleDecompositionMultiplicitiesCache::Value {
        basic_summands,
        multiplicities,
    }) {
        ModuleDecompositionMultiplicitiesCache::Fail => Ok(None),
        ModuleDecompositionMultiplicitiesCache::Value {
            basic_summands,
            multiplicities,
        } => Ok(Some(ModuleDecompositionWithMultiplicities {
            basic_summands,
            multiplicities,
        })),
    }
}

/// Translate QPA's `BlockSplittingIdempotents(M)` operation.
///
/// The source quotients `EndOverAlgebra(M)` by its radical, takes the central
/// idempotents of the semisimple quotient, and lifts them successively with
/// `LiftIdempotent` and `LiftTwoOrthogonalIdempotents` before converting the
/// resulting matrices back to module endomorphisms.
///
/// # Errors
///
/// Propagates endomorphism-algebra, radical, quotient, central-idempotent,
/// lifting, and module-map conversion failures.
pub fn block_splitting_idempotents<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, ModuleDecompositionError> {
    if module.dimension() == 0 {
        return Ok(Vec::new());
    }
    let endomorphism_algebra = end_over_algebra(module)?;
    let structure_constants = endomorphism_algebra.structure_constants_algebra()?;
    let radical = structure_constants.radical()?;
    let quotient = structure_constants.quotient_by_ideal_basis(radical.basis().clone())?;
    let top = quotient.algebra().central_idempotents_split_prime()?;
    let Some(first_top) = top.first() else {
        return Err(ModuleDecompositionError::EmptySemisimpleEndomorphismTop);
    };
    let mapping = quotient.natural_homomorphism()?;
    let mut orthogonal_idempotents = vec![lift_idempotent(&mapping, first_top)?];
    for next_top in top.iter().skip(1) {
        let accumulated = orthogonal_idempotents.iter().fold(
            vec![structure_constants.field().zero(); structure_constants.dimension()],
            |sum, idempotent| add_rows(structure_constants.field(), &sum, idempotent),
        );
        let pair = lift_two_orthogonal_idempotents(&mapping, &accumulated, next_top)?;
        orthogonal_idempotents.push(pair[1].clone());
    }
    orthogonal_idempotents
        .into_iter()
        .map(|coordinates| {
            let matrix = endomorphism_algebra.matrix_from_coefficients(&coordinates)?;
            Ok(from_end_m_to_hom_mm(module, &matrix)?)
        })
        .collect()
}

/// Translate QPA's `BlockDecompositionOfModule(M)` operation.
///
/// # Errors
///
/// Propagates [`block_splitting_idempotents`] and canonical image
/// construction.
pub fn block_decomposition_of_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModule<F>>, ModuleDecompositionError> {
    block_splitting_idempotents(module)?
        .iter()
        .map(image_of_module_homomorphism)
        .collect::<Result<Vec<_>, _>>()
        .map_err(Into::into)
}

/// Translate QPA's `AllSimpleSubmodulesOfModule(M)` operation.
///
/// QPA passes to the socle, splits its isotypic blocks with
/// `BlockSplittingIdempotents`, enumerates every one-dimensional field
/// subspace of each block in GAP order, and closes its canonical basis under
/// the module action. The resulting inclusions are composed back into `M`.
///
/// # Errors
///
/// Rejects infinite coefficient fields and propagates socle, block-splitting,
/// subrepresentation, image, and composition failures.
pub fn all_simple_submodules_of_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, ModuleDecompositionError> {
    let field = module.right_acting_algebra().field();
    let elements = field
        .qpa_elements()?
        .ok_or(ModuleDecompositionError::FiniteFieldRequired)?;
    if module.dimension() == 0 || is_simple_qpa_module(module) {
        return Ok(vec![identity_module_mapping(module.clone())?]);
    }

    let socle_inclusion = socle_of_module_inclusion(module.clone())?;
    let socle = socle_inclusion.source();
    let mut simple_inclusions = Vec::new();
    for block_idempotent in block_splitting_idempotents(socle)? {
        let block_inclusion = image_inclusion(&block_idempotent)?;
        let block = block_inclusion.source();
        let basis = block.basis_elements();
        for coordinates in one_dimensional_subspace_bases(field, &elements, block.dimension()) {
            let mut generator = block.zero_element();
            for (basis_element, coefficient) in basis.iter().zip(&coordinates) {
                generator = generator.add(&basis_element.scale(coefficient))?;
            }
            let simple = subrepresentation_inclusion(block.clone(), &[generator])?;
            simple_inclusions.push(
                simple
                    .compose(&block_inclusion)?
                    .compose(&socle_inclusion)?,
            );
        }
    }
    Ok(simple_inclusions)
}

/// Translate QPA's `AllSubmodulesOfModule(M)` operation.
///
/// Starting with all simple submodules, the source repeatedly quotients by a
/// submodule, enumerates the simple submodules of that quotient, and pulls
/// them back. Equal underlying row spaces are removed within each dimension
/// in their first-occurrence order.
///
/// # Errors
///
/// Rejects infinite coefficient fields and propagates simple-submodule,
/// cokernel, pullback, map-image, and exact row-space failures.
pub fn all_submodules_of_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<AllSubmodulesOfModuleResult<F>, ModuleDecompositionError> {
    if module
        .right_acting_algebra()
        .field()
        .qpa_elements()?
        .is_none()
    {
        return Err(ModuleDecompositionError::FiniteFieldRequired);
    }
    if module.dimension() == 0 {
        return Ok(AllSubmodulesOfModuleResult::Flat(vec![
            identity_module_mapping(module.clone())?,
        ]));
    }

    let zero = subrepresentation_inclusion(module.clone(), &[])?;
    if is_simple_qpa_module(module) {
        return Ok(AllSubmodulesOfModuleResult::Flat(vec![
            zero,
            identity_module_mapping(module.clone())?,
        ]));
    }

    let simples = all_simple_submodules_of_module(module)?;
    let mut levels = vec![vec![zero], simples.clone()];
    let mut previous = simples;
    let mut submodule_dimension = 1_usize;
    while submodule_dimension < module.dimension() {
        submodule_dimension += 1;
        let mut spaces = Vec::<DenseMatrix<F>>::new();
        let mut next = Vec::new();
        for submodule in &previous {
            let quotient = cokernel_projection(submodule)?;
            for simple in all_simple_submodules_of_module(quotient.range())? {
                let pullback = pull_back(&quotient, &simple)?
                    .ok_or(PathAlgebraModuleError::InternalCommonRangePullbackFailed)?;
                let inclusion = pullback.1;
                let space = underlying_submodule_row_space(&inclusion)?;
                if !spaces.contains(&space) {
                    spaces.push(space);
                    next.push(inclusion);
                }
            }
        }
        levels.push(next.clone());
        previous = next;
    }
    Ok(AllSubmodulesOfModuleResult::ByDimension(levels))
}

/// Translate QPA's `MakeBasicListOfModules(L)` operation.
///
/// QPA walks the entered list from left to right, removes from each later
/// module the maximal common direct summand with the current module, filters
/// removed entries, and finally decomposes every retained remainder. This
/// preserves the source's early termination when a previously removed entry
/// becomes the current outer-loop position.
///
/// # Errors
///
/// Rejects an empty list, different acting algebras, infinite fields, and
/// infinite-dimensional algebras; propagates maximal-common-summand and
/// deterministic-decomposition failures.
pub fn make_basic_list_of_modules<F: QpaCoefficientField>(
    modules: &[PathAlgebraModule<F>],
) -> Result<Vec<PathAlgebraModule<F>>, ModuleDecompositionError> {
    let Some(first) = modules.first() else {
        return Err(ModuleDecompositionError::EmptyBasicModuleList);
    };
    let algebra = first.right_acting_algebra();
    if modules
        .iter()
        .any(|module| module.right_acting_algebra() != algebra)
    {
        return Err(ModuleDecompositionError::BasicModulesOverDifferentAlgebras);
    }
    if algebra.field().qpa_elements()?.is_none() {
        return Err(ModuleDecompositionError::FiniteFieldRequired);
    }
    let finite_dimensional = match algebra {
        QuiverAlgebra::Path(path_algebra) => path_algebra.has_finite_canonical_basis(),
        QuiverAlgebra::Quotient(quotient) => quotient.is_finite_dimensional()?,
    };
    if !finite_dimensional {
        return Err(ModuleDecompositionError::InfiniteDimensionalActingAlgebra);
    }

    let mut remainders = modules.iter().cloned().map(Some).collect::<Vec<_>>();
    for current in 0..modules.len().saturating_sub(1) {
        let Some(current_module) = remainders[current].clone() else {
            break;
        };
        for remainder in remainders.iter_mut().skip(current + 1) {
            let Some(later_module) = remainder.clone() else {
                continue;
            };
            if let Some(common) = maximal_common_direct_summand(&current_module, &later_module)? {
                *remainder = (common.range_remainder().dimension() > 0)
                    .then(|| common.range_remainder().clone());
            }
        }
    }

    let mut basic = Vec::new();
    for module in remainders.into_iter().flatten() {
        if module.dimension() == 0 {
            return Err(ModuleDecompositionError::ZeroModuleInBasicList);
        }
        basic.extend(decompose_module(&module)?);
    }
    Ok(basic)
}

fn one_dimensional_subspace_bases<F: Field>(
    field: &F,
    elements: &[F::Element],
    dimension: usize,
) -> Vec<Vec<F::Element>> {
    debug_assert!(!elements.is_empty());
    let mut result = Vec::new();
    for pivot in 0..dimension {
        let suffix_length = dimension - pivot - 1;
        let mut indices = vec![0_usize; suffix_length];
        loop {
            let mut row = vec![field.zero(); dimension];
            row[pivot] = field.one();
            for (offset, &index) in indices.iter().enumerate() {
                row[pivot + 1 + offset] = elements[index].clone();
            }
            result.push(row);
            if suffix_length == 0 {
                break;
            }
            let mut position = suffix_length;
            let mut exhausted = false;
            loop {
                position -= 1;
                indices[position] += 1;
                if indices[position] < elements.len() {
                    break;
                }
                indices[position] = 0;
                if position == 0 {
                    exhausted = true;
                    break;
                }
            }
            if exhausted {
                break;
            }
        }
    }
    result
}

fn underlying_submodule_row_space<F: Field>(
    inclusion: &PathAlgebraModuleHomomorphism<F>,
) -> Result<DenseMatrix<F>, ModuleDecompositionError> {
    let rows = inclusion
        .source()
        .basis_elements()
        .iter()
        .map(|basis_element| {
            inclusion
                .image(basis_element)
                .map(|image| image.basis_coefficients())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let entries = rows.into_iter().flatten().collect();
    Ok(DenseMatrix::from_flat(
        inclusion.source().dimension(),
        inclusion.range().dimension(),
        entries,
        inclusion.range().right_acting_algebra().field().clone(),
    )?
    .canonical_row_space_basis())
}

/// Translate QPA's `BasicVersionOfModule(M)` operation.
///
/// The zero-module branch returns the identical entered object. Otherwise
/// QPA takes the basic summands from
/// `DecomposeModuleWithMultiplicities(M)` and builds their ordered direct sum.
///
/// # Errors
///
/// Propagates decomposition with multiplicities and direct-sum construction.
pub fn basic_version_of_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, ModuleDecompositionError> {
    if module.dimension() == 0 {
        return Ok(module.clone());
    }
    let decomposition = decompose_module_with_multiplicities(module)?
        .ok_or(ModuleDecompositionError::MissingNonzeroMultiplicityDecomposition)?;
    direct_sum_of_qpa_modules(decomposition.basic_summands())?
        .ok_or(ModuleDecompositionError::BasicDirectSumFailed)
}

/// Translate QPA's deterministic `DecomposeModuleViaTop(M)` operation.
///
/// QPA sends the ordered endomorphism basis through `TopOfModule`, constructs
/// the generated algebra of induced top matrices, decomposes that algebra,
/// and lifts its complete orthogonal idempotent set through the resulting
/// algebra homomorphism. The lifted endomorphisms are returned through their
/// canonical images.
///
/// # Errors
///
/// Requires a finite supported coefficient field and propagates top-functor,
/// generated-matrix-algebra, idempotent, lifting, and image failures.
pub fn decompose_module_via_top<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModule<F>>, ModuleDecompositionError> {
    if module.right_acting_algebra().field().qpa_characteristic() == QpaFieldCharacteristic::Zero {
        return Err(ModuleDecompositionError::FiniteFieldRequired);
    }
    let endomorphism_algebra = end_over_algebra(module)?;
    let source_algebra = endomorphism_algebra.structure_constants_algebra()?;
    let top_matrices = endomorphism_algebra
        .basis_matrices()
        .iter()
        .map(|matrix| {
            let homomorphism = from_end_m_to_hom_mm(module, matrix)?;
            let top = top_of_module_homomorphism(&homomorphism)?;
            from_hom_mm_to_end_m(&top)
        })
        .collect::<Result<Vec<_>, PathAlgebraModuleError>>()?;
    let top_algebra = GeneratedMatrixAlgebra::from_generators(
        module.right_acting_algebra().field().clone(),
        &top_matrices,
    )?;
    let mut basis_images = Vec::with_capacity(
        source_algebra
            .dimension()
            .checked_mul(top_algebra.algebra.dimension())
            .ok_or(ModuleDecompositionError::GeneratedMatrixAlgebraSizeOverflow)?,
    );
    for matrix in &top_matrices {
        basis_images.extend(
            top_algebra
                .coordinates(matrix)?
                .ok_or(ModuleDecompositionError::TopImageOutsideGeneratedAlgebra)?,
        );
    }
    let mapping = StructureConstantsAlgebraHomomorphism::by_basis_images(
        source_algebra,
        top_algebra.algebra.clone(),
        DenseMatrix::from_flat(
            endomorphism_algebra.dimension(),
            top_algebra.algebra.dimension(),
            basis_images,
            module.right_acting_algebra().field().clone(),
        )?,
    )?;
    let top_idempotents = top_algebra.idempotents_for_decomposition()?;
    lifting_complete_set_of_orthogonal_idempotents(&mapping, &top_idempotents)?
        .into_iter()
        .map(|coordinates| {
            let matrix = endomorphism_algebra.matrix_from_coefficients(&coordinates)?;
            let homomorphism = from_end_m_to_hom_mm(module, &matrix)?;
            Ok(image_of_module_homomorphism(&homomorphism)?)
        })
        .collect()
}

/// Reproducible default seed for QPA's experimental probabilistic module
/// decomposition.
pub const MODULE_DECOMPOSITION_DEFAULT_SEED: u64 = 0x5150_414d_4445_4301;

/// QPA's heterogeneous `DecomposeModuleProbabilistic` return contract.
#[derive(Clone, Debug)]
pub enum ProbabilisticModuleDecomposition<F: Field> {
    /// No splitting element was found and QPA returned the entered module.
    Unsplit(PathAlgebraModule<F>),
    /// At least one splitting succeeded and QPA returned the flattened list of
    /// recursively obtained summands.
    Split(Vec<PathAlgebraModule<F>>),
}

impl<F: Field> ProbabilisticModuleDecomposition<F> {
    /// Consume the heterogeneous result as a normalized summand list.
    #[must_use]
    pub fn into_modules(self) -> Vec<PathAlgebraModule<F>> {
        match self {
            Self::Unsplit(module) => vec![module],
            Self::Split(modules) => modules,
        }
    }
}

/// Translate QPA's experimental `DecomposeModuleProbabilistic(HomMM, M)`
/// operation with a reproducible default random source.
///
/// # Errors
///
/// Propagates [`decompose_module_probabilistic_with_rng`].
pub fn decompose_module_probabilistic<F: QpaCoefficientField>(
    endomorphism_basis: &[PathAlgebraModuleHomomorphism<F>],
    module: &PathAlgebraModule<F>,
) -> Result<ProbabilisticModuleDecomposition<F>, ModuleDecompositionError> {
    let mut rng = StdRng::seed_from_u64(MODULE_DECOMPOSITION_DEFAULT_SEED);
    decompose_module_probabilistic_with_rng(endomorphism_basis, module, &mut rng)
}

/// Translate QPA's experimental `DecomposeModuleProbabilistic(HomMM, M)`
/// operation using the caller's random source.
///
/// Rust draws each full-row-space coordinate uniformly from GAP's ordered
/// finite-field element list. Its engine state is intentionally explicit;
/// the square-power test, attempt bound, image/kernel splitting, restricted
/// endomorphism spans, and recursive flattening are the literal source
/// algorithm.
///
/// # Errors
///
/// Requires a nonempty endomorphism list over the entered nonzero module and
/// a finite supported coefficient field. Exact map, determinant,
/// factorization, isomorphism, and recursive failures propagate.
pub fn decompose_module_probabilistic_with_rng<F: QpaCoefficientField, R: Rng + ?Sized>(
    endomorphism_basis: &[PathAlgebraModuleHomomorphism<F>],
    module: &PathAlgebraModule<F>,
    rng: &mut R,
) -> Result<ProbabilisticModuleDecomposition<F>, ModuleDecompositionError> {
    let elements = module
        .right_acting_algebra()
        .field()
        .qpa_elements()?
        .ok_or(ModuleDecompositionError::FiniteFieldRequired)?;
    if endomorphism_basis.is_empty() {
        return Err(ModuleDecompositionError::EmptyProbabilisticEndomorphismList);
    }
    if module.dimension() == 0 {
        return Err(ModuleDecompositionError::ZeroProbabilisticModule);
    }
    if endomorphism_basis
        .iter()
        .any(|homomorphism| homomorphism.source() != module || homomorphism.range() != module)
    {
        return Err(ModuleDecompositionError::ProbabilisticEndomorphismModuleMismatch);
    }
    let maximum_dimension = *module
        .dimension_vector()
        .iter()
        .max()
        .ok_or(ModuleDecompositionError::ZeroProbabilisticModule)?;
    if maximum_dimension == 0 {
        return Err(ModuleDecompositionError::ZeroProbabilisticModule);
    }
    let logarithmic_bound = if maximum_dimension == 1 {
        0usize
    } else {
        usize::BITS as usize - (maximum_dimension - 1).leading_zeros() as usize
    };
    let attempt_limit = elements
        .len()
        .checked_mul(endomorphism_basis.len())
        .and_then(|product| product.checked_mul(2))
        .ok_or(ModuleDecompositionError::ProbabilisticAttemptBoundOverflow)?;

    for _ in 1..attempt_limit {
        let coefficients = (0..endomorphism_basis.len())
            .map(|_| elements[rng.random_range(0..elements.len())].clone())
            .collect::<Vec<_>>();
        let candidate =
            linear_combination_of_module_homomorphisms(endomorphism_basis, &coefficients)?;
        let candidate_matrix = from_hom_mm_to_end_m(&candidate)?;
        if !candidate_matrix
            .field()
            .is_zero(&candidate_matrix.determinant()?)
        {
            continue;
        }
        let mut power = candidate;
        for _ in 0..=logarithmic_bound {
            power = power.compose(&power)?;
        }
        if power.is_zero() {
            continue;
        }

        let image_projection = image_projection(&power)?;
        let image_inclusion = image_inclusion(&power)?;
        let kernel_inclusion = kernel_inclusion(&power)?;
        let cokernel_projection = cokernel_projection(&power)?;
        let kernel_isomorphism =
            isomorphism_of_modules(cokernel_projection.range(), kernel_inclusion.source())?
                .ok_or(ModuleDecompositionError::MissingProbabilisticKernelIsomorphism)?;
        let kernel_projection = cokernel_projection.compose(&kernel_isomorphism)?;

        let image_endomorphisms =
            restrict_endomorphism_basis(endomorphism_basis, &image_inclusion, &image_projection)?;
        let kernel_endomorphisms =
            restrict_endomorphism_basis(endomorphism_basis, &kernel_inclusion, &kernel_projection)?;
        let image_result = decompose_module_probabilistic_with_rng(
            &image_endomorphisms,
            image_projection.range(),
            rng,
        )?;
        let kernel_result = decompose_module_probabilistic_with_rng(
            &kernel_endomorphisms,
            kernel_projection.range(),
            rng,
        )?;
        let mut modules = image_result.into_modules();
        modules.extend(kernel_result.into_modules());
        return Ok(ProbabilisticModuleDecomposition::Split(modules));
    }
    Ok(ProbabilisticModuleDecomposition::Unsplit(module.clone()))
}

/// Reproducible default seed for QPA's experimental characteristic-polynomial
/// module decomposition.
pub const MODULE_DECOMPOSITION_CHAR_POLY_DEFAULT_SEED: u64 = 0x5150_414d_4445_4302;

/// Translate QPA's experimental `DecomposeModuleViaCharPoly(M)` operation
/// with a reproducible default random source.
///
/// # Errors
///
/// Propagates [`decompose_module_via_char_poly_with_rng`].
pub fn decompose_module_via_char_poly<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModule<F>>, ModuleDecompositionError> {
    let mut rng = StdRng::seed_from_u64(MODULE_DECOMPOSITION_CHAR_POLY_DEFAULT_SEED);
    decompose_module_via_char_poly_with_rng(module, &mut rng)
}

/// Translate QPA's experimental `DecomposeModuleViaCharPoly(M)` operation
/// using the caller's random source.
///
/// Each trial chooses a uniform random element of `End(M)`, factors the
/// characteristic polynomials of its nonzero vertex matrices, evaluates the
/// maximal primary powers on the endomorphism, and recursively decomposes the
/// resulting kernels. The trial count and the small-module deterministic
/// branch are exactly those in `lib/moduledecomp.gi`.
///
/// # Errors
///
/// Requires a supported finite coefficient field and propagates exact Hom,
/// characteristic-polynomial, finite-field factorization, polynomial-matrix
/// evaluation, kernel, and recursive decomposition failures.
pub fn decompose_module_via_char_poly_with_rng<F: QpaCoefficientField, R: Rng + ?Sized>(
    module: &PathAlgebraModule<F>,
    rng: &mut R,
) -> Result<Vec<PathAlgebraModule<F>>, ModuleDecompositionError> {
    let field = module.right_acting_algebra().field();
    let elements = field
        .qpa_elements()?
        .ok_or(ModuleDecompositionError::FiniteFieldRequired)?;
    let dimensions = module.dimension_vector();
    let total_dimension = module.dimension();
    if total_dimension < dimensions.len().saturating_add(1) {
        return decompose_module(module);
    }

    let endomorphism_basis = hom_over_algebra(module, module)?
        .ok_or(ModuleDecompositionError::MissingCharPolyEndomorphismSpace)?;
    if endomorphism_basis.is_empty() {
        return Err(ModuleDecompositionError::EmptyCharPolyEndomorphismList);
    }
    let maximum_runs = 40 / total_dimension + 10;
    for _ in 0..maximum_runs {
        let coefficients = (0..endomorphism_basis.len())
            .map(|_| elements[rng.random_range(0..elements.len())].clone())
            .collect::<Vec<_>>();
        let candidate =
            linear_combination_of_module_homomorphisms(&endomorphism_basis, &coefficients)?;

        let mut factor_lists = Vec::new();
        for matrix in candidate.vertex_matrices() {
            if matrix.is_zero() {
                continue;
            }
            let characteristic = matrix.characteristic_polynomial_coefficients()?;
            factor_lists.push(factor_finite_field_polynomial(
                field,
                &elements,
                &characteristic,
                rng,
            )?);
        }
        let mut occurring_primes = Vec::new();
        for factors in &factor_lists {
            for factor in factors {
                if !occurring_primes.contains(factor) {
                    occurring_primes.push(factor.clone());
                }
            }
        }
        if occurring_primes.len() <= 1 {
            continue;
        }

        let mut result = Vec::new();
        for prime in &occurring_primes {
            let maximum_multiplicity = factor_lists
                .iter()
                .map(|factors| factors.iter().filter(|factor| *factor == prime).count())
                .max()
                .unwrap_or(0);
            let primary_polynomial = polynomial_power(field, prime, maximum_multiplicity);
            let mut matrices = Vec::with_capacity(dimensions.len());
            for (&dimension, matrix) in dimensions.iter().zip(candidate.vertex_matrices()) {
                matrices.push(if dimension == 0 {
                    DenseMatrix::zero(1, 1, field.clone())
                } else {
                    evaluate_matrix_polynomial(matrix, &primary_polynomial)?
                });
            }
            let primary_map =
                right_module_hom_over_algebra(module.clone(), module.clone(), matrices)?;
            let kernel = kernel_of_module_homomorphism(&primary_map)?;
            result.extend(decompose_module_via_char_poly_with_rng(&kernel, rng)?);
        }
        return Ok(result);
    }

    Ok(vec![module.clone()])
}

fn polynomial_power<F: Field>(
    field: &F,
    polynomial: &[F::Element],
    exponent: usize,
) -> Vec<F::Element> {
    (0..exponent).fold(vec![field.one()], |power, _| {
        polynomial_product(field, &power, polynomial)
    })
}

fn evaluate_matrix_polynomial<F: Field>(
    matrix: &DenseMatrix<F>,
    polynomial: &[F::Element],
) -> Result<DenseMatrix<F>, MatrixError> {
    let identity = DenseMatrix::identity(matrix.nrows(), matrix.field().clone());
    let mut result = DenseMatrix::zero(matrix.nrows(), matrix.ncols(), matrix.field().clone());
    for coefficient in polynomial.iter().rev() {
        result = result.multiply(matrix)?;
        result = result.add(&identity.scale_left(coefficient))?;
    }
    Ok(result)
}

fn decomposition_endomorphisms<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, ModuleDecompositionError> {
    let endomorphism_algebra = end_over_algebra(module)?;
    let (transposed, endomorphism_basis) =
        endomorphism_algebra.transposed_structure_constants_model()?;
    transposed
        .idempotents_for_decomposition()?
        .into_iter()
        .map(|idempotent| {
            let mut matrix = endomorphism_algebra.zero();
            for (coefficient, basis) in idempotent.coordinates().iter().zip(&endomorphism_basis) {
                matrix = matrix.add(&basis.scale_left(coefficient))?;
            }
            Ok(from_end_m_to_hom_mm(module, &matrix)?)
        })
        .collect()
}

fn linear_combination_of_module_homomorphisms<F: Field>(
    basis: &[PathAlgebraModuleHomomorphism<F>],
    coefficients: &[F::Element],
) -> Result<PathAlgebraModuleHomomorphism<F>, ModuleDecompositionError> {
    if basis.len() != coefficients.len() {
        return Err(ModuleDecompositionError::ProbabilisticCoefficientCount {
            actual: coefficients.len(),
            expected: basis.len(),
        });
    }
    let first = basis
        .first()
        .ok_or(ModuleDecompositionError::EmptyProbabilisticEndomorphismList)?;
    let field = first.source().right_acting_algebra().field();
    let mut combination = first.scale(&field.zero())?;
    for (homomorphism, coefficient) in basis.iter().zip(coefficients) {
        if !field.is_zero(coefficient) {
            combination = combination.add(&homomorphism.scale(coefficient)?)?;
        }
    }
    Ok(combination)
}

fn restrict_endomorphism_basis<F: Field>(
    endomorphism_basis: &[PathAlgebraModuleHomomorphism<F>],
    inclusion: &PathAlgebraModuleHomomorphism<F>,
    projection: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, ModuleDecompositionError> {
    let mut unique = Vec::new();
    for endomorphism in endomorphism_basis {
        let restricted = inclusion.compose(endomorphism)?.compose(projection)?;
        if !unique.iter().any(|existing| existing == &restricted) {
            unique.push(restricted);
        }
    }
    let dimension = projection.range().dimension();
    if dimension == 0 {
        return Err(ModuleDecompositionError::ZeroProbabilisticRestrictedModule);
    }
    let flattened_dimension = dimension
        .checked_mul(dimension)
        .ok_or(ModuleDecompositionError::GeneratedMatrixAlgebraSizeOverflow)?;
    let field = projection.range().right_acting_algebra().field().clone();
    let rows = DenseMatrix::from_flat(
        unique.len(),
        flattened_dimension,
        unique
            .iter()
            .map(from_hom_mm_to_end_m)
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flat_map(|matrix| matrix.as_flat().to_vec())
            .collect(),
        field,
    )?
    .semi_echelon_transformation()
    .vectors;
    if rows.nrows() == 0 {
        return Err(ModuleDecompositionError::EmptyProbabilisticRestrictedBasis);
    }
    (0..rows.nrows())
        .map(|row| {
            let matrix = DenseMatrix::from_flat(
                dimension,
                dimension,
                rows.row(row).to_vec(),
                rows.field().clone(),
            )?;
            Ok(from_end_m_to_hom_mm(projection.range(), &matrix)?)
        })
        .collect()
}

struct GeneratedMatrixAlgebra<F: Field> {
    algebra: StructureConstantsAlgebra<F>,
    basis_matrices: Vec<DenseMatrix<F>>,
    flattened_basis: DenseMatrix<F>,
}

impl<F: Field> GeneratedMatrixAlgebra<F> {
    fn from_generators(
        field: F,
        generators: &[DenseMatrix<F>],
    ) -> Result<Self, ModuleDecompositionError> {
        let Some(first) = generators.first() else {
            return Err(ModuleDecompositionError::EmptyGeneratedMatrixAlgebra);
        };
        if first.nrows() == 0 || first.nrows() != first.ncols() || first.field() != &field {
            return Err(ModuleDecompositionError::InvalidGeneratedMatrixAlgebraGenerator);
        }
        let matrix_dimension = first.nrows();
        if generators.iter().any(|matrix| {
            matrix.field() != &field
                || matrix.nrows() != matrix_dimension
                || matrix.ncols() != matrix_dimension
        }) {
            return Err(ModuleDecompositionError::InvalidGeneratedMatrixAlgebraGenerator);
        }
        let flattened_dimension = matrix_dimension
            .checked_mul(matrix_dimension)
            .ok_or(ModuleDecompositionError::GeneratedMatrixAlgebraSizeOverflow)?;
        let flattened_basis = DenseMatrix::from_flat(
            generators.len(),
            flattened_dimension,
            generators
                .iter()
                .flat_map(|matrix| matrix.as_flat().iter().cloned())
                .collect(),
            field.clone(),
        )?
        .semi_echelon_transformation()
        .vectors;
        if flattened_basis.nrows() == 0 {
            return Err(ModuleDecompositionError::EmptyGeneratedMatrixAlgebra);
        }
        let basis = (0..flattened_basis.nrows())
            .map(|row| {
                DenseMatrix::from_flat(
                    matrix_dimension,
                    matrix_dimension,
                    flattened_basis.row(row).to_vec(),
                    field.clone(),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let dimension = basis.len();
        let capacity = dimension
            .checked_mul(dimension)
            .and_then(|square| square.checked_mul(dimension))
            .ok_or(ModuleDecompositionError::GeneratedMatrixAlgebraSizeOverflow)?;
        let mut products = Vec::with_capacity(capacity);
        for left in &basis {
            for right in &basis {
                let product = left.multiply(right)?;
                products.extend(
                    flattened_basis
                        .solution_row(product.as_flat())?
                        .ok_or(ModuleDecompositionError::GeneratedMatrixProductOutsideSpan)?,
                );
            }
        }
        Ok(Self {
            algebra: StructureConstantsAlgebra::new(field, dimension, products)?,
            basis_matrices: basis,
            flattened_basis,
        })
    }

    fn coordinates(&self, matrix: &DenseMatrix<F>) -> Result<Option<Vec<F::Element>>, MatrixError> {
        self.flattened_basis.solution_row(matrix.as_flat())
    }
}

impl<F: QpaCoefficientField> GeneratedMatrixAlgebra<F> {
    fn idempotents_for_decomposition(
        &self,
    ) -> Result<Vec<Vec<F::Element>>, ModuleDecompositionError> {
        let mut idempotents = self
            .algebra
            .idempotents_for_decomposition()?
            .into_iter()
            .map(|idempotent| idempotent.coordinates().to_vec())
            .collect::<Vec<_>>();
        for index in 1..idempotents.len() {
            let mut position = index;
            while position > 0
                && self.compare_matrix_coordinates(
                    &idempotents[position],
                    &idempotents[position - 1],
                )? == Ordering::Less
            {
                idempotents.swap(position, position - 1);
                position -= 1;
            }
        }
        Ok(idempotents)
    }

    fn compare_matrix_coordinates(
        &self,
        left: &[F::Element],
        right: &[F::Element],
    ) -> Result<Ordering, ModuleDecompositionError> {
        let left = self.matrix_from_coordinates(left)?;
        let right = self.matrix_from_coordinates(right)?;
        for (left, right) in left.as_flat().iter().zip(right.as_flat()) {
            match self.algebra.field().qpa_compare(left, right)? {
                Ordering::Equal => {}
                ordering => return Ok(ordering),
            }
        }
        Ok(Ordering::Equal)
    }

    fn matrix_from_coordinates(
        &self,
        coordinates: &[F::Element],
    ) -> Result<DenseMatrix<F>, ModuleDecompositionError> {
        if coordinates.len() != self.basis_matrices.len() {
            return Err(ModuleDecompositionError::GeneratedMatrixCoordinateCount {
                actual: coordinates.len(),
                expected: self.basis_matrices.len(),
            });
        }
        let dimension = self.basis_matrices[0].nrows();
        let mut matrix = DenseMatrix::zero(dimension, dimension, self.algebra.field().clone());
        for (coefficient, basis) in coordinates.iter().zip(&self.basis_matrices) {
            matrix = matrix.add(&basis.scale_left(coefficient))?;
        }
        Ok(matrix)
    }
}

fn multiply_row_matrix<F: Field>(row: &[F::Element], matrix: &DenseMatrix<F>) -> Vec<F::Element> {
    debug_assert_eq!(row.len(), matrix.nrows());
    let mut result = vec![matrix.field().zero(); matrix.ncols()];
    for (inner, coefficient) in row.iter().enumerate() {
        if matrix.field().is_zero(coefficient) {
            continue;
        }
        for column in 0..matrix.ncols() {
            matrix.field().add_mul_assign(
                &mut result[column],
                coefficient,
                &matrix[(inner, column)],
            );
        }
    }
    result
}

fn is_zero_row<F: Field>(field: &F, row: &[F::Element]) -> bool {
    row.iter().all(|entry| field.is_zero(entry))
}

fn binomial_integer(n: usize, k: usize) -> Integer {
    let k = k.min(n - k);
    let mut result = Integer::from(1);
    for step in 1..=k {
        result *= Integer::from(n - k + step);
        result /= Integer::from(step);
    }
    result
}

fn algebra_power<F: Field>(
    algebra: &StructureConstantsAlgebra<F>,
    element: &[F::Element],
    exponent: usize,
    identity: &[F::Element],
) -> Result<Vec<F::Element>, StructureConstantsError> {
    let mut result = identity.to_vec();
    let mut base = element.to_vec();
    let mut exponent = exponent;
    while exponent > 0 {
        if !exponent.is_multiple_of(2) {
            result = algebra.multiply_coordinates(&result, &base)?;
        }
        exponent /= 2;
        if exponent > 0 {
            base = algebra.multiply_coordinates(&base, &base)?;
        }
    }
    Ok(result)
}

fn gap_sort_and_deduplicate<F: QpaCoefficientField>(
    field: &F,
    rows: &mut Vec<Vec<F::Element>>,
) -> Result<(), ModuleDecompositionError> {
    for index in 1..rows.len() {
        let mut position = index;
        while position > 0
            && gap_compare_rows(field, &rows[position], &rows[position - 1])? == Ordering::Less
        {
            rows.swap(position, position - 1);
            position -= 1;
        }
    }
    rows.dedup();
    Ok(())
}

fn gap_compare_rows<F: QpaCoefficientField>(
    field: &F,
    left: &[F::Element],
    right: &[F::Element],
) -> Result<Ordering, CoefficientDomainError> {
    for (left, right) in left.iter().zip(right) {
        match field.qpa_compare(left, right)? {
            Ordering::Equal => {}
            ordering => return Ok(ordering),
        }
    }
    Ok(left.len().cmp(&right.len()))
}

/// Failure in QPA's module-decomposition foundations.
#[derive(Debug, Error)]
pub enum ModuleDecompositionError {
    /// `ComplementInFullRowSpace` rejects a zero-dimensional entered space.
    #[error("cannot find complement of trivial space")]
    TrivialRowSpace,
    /// A translated triangulated basis row unexpectedly had no normalized
    /// pivot.
    #[error("triangulated row {row} has no unit pivot")]
    TriangulatedRowWithoutPivot {
        /// One-based row position.
        row: usize,
    },
    /// The source, range, and image matrix of an algebra mapping must use one
    /// exact field.
    #[error("the algebra mapping fields differ")]
    MappingFieldsDiffer,
    /// The basis-image matrix has the wrong source/range dimensions.
    #[error("algebra mapping matrix has shape {actual:?}, expected {expected:?}")]
    MappingShape {
        /// Entered matrix shape.
        actual: (usize, usize),
        /// Required `dim(source) × dim(range)` shape.
        expected: (usize, usize),
    },
    /// Entered basis images fail one multiplicativity equation.
    #[error("algebra mapping does not preserve basis product ({left},{right})")]
    MappingNotMultiplicative {
        /// One-based left source basis position.
        left: usize,
        /// One-based right source basis position.
        right: usize,
    },
    /// A source coordinate row has the wrong length.
    #[error("mapping source row has length {actual}, expected {expected}")]
    MappingSourceLength {
        /// Entered length.
        actual: usize,
        /// Source algebra dimension.
        expected: usize,
    },
    /// A range coordinate row has the wrong length.
    #[error("mapping range row has length {actual}, expected {expected}")]
    MappingRangeLength {
        /// Entered length.
        actual: usize,
        /// Range algebra dimension.
        expected: usize,
    },
    /// QPA indexes `ids[1]` before entering its lifting loop.
    #[error("the entered idempotent list is empty")]
    EmptyIdempotentSet,
    /// `PreImagesRepresentativeNC` was requested outside the mapping image.
    #[error("entered element {position} has no mapping preimage")]
    MissingPreimage {
        /// One-based entered idempotent position.
        position: usize,
    },
    /// `LiftIdempotent` requires an idempotent range element.
    #[error("the entered range element is not idempotent")]
    EnteredElementNotIdempotent,
    /// `LiftIdempotent` requires the range element to lie in the image.
    #[error("the entered range element is outside the mapping image")]
    EnteredElementOutsideImage,
    /// `LiftTwoOrthogonalIdempotents` requires orthogonality in the range.
    #[error("the entered idempotents are not orthogonal in the mapping range")]
    EnteredIdempotentsNotOrthogonal,
    /// The lifting formulas require a multiplicative identity in the source.
    #[error("the source algebra has no multiplicative identity")]
    MissingSourceIdentity,
    /// Complete-set validation requires a multiplicative identity in the
    /// mapping range.
    #[error("the range algebra has no multiplicative identity")]
    MissingRangeIdentity,
    /// A complete orthogonal idempotent set must sum to the range identity.
    #[error("the entered idempotents are not a complete set")]
    EnteredIdempotentsNotComplete,
    /// A checked nilpotent-series index overflowed.
    #[error("the lifting nilpotence index overflows usize")]
    NilpotentSeriesOverflow,
    /// The source algorithms assume a nilpotent mapping kernel.
    #[error("the algebra mapping kernel is not nilpotent")]
    MappingKernelNotNilpotent,
    /// The source's Newton correction failed to reach an idempotent within the
    /// finite-dimensional nilpotent-kernel bound.
    #[error("the idempotent correction did not converge")]
    IdempotentCorrectionDidNotConverge,
    /// Exact shaped-matrix failure.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// Finite-dimensional structure-constant algebra failure.
    #[error(transparent)]
    StructureConstants(#[from] StructureConstantsError),
    /// QPA primitive-idempotent failure.
    #[error(transparent)]
    Hopf(#[from] HopfError),
    /// GAP-compatible coefficient ordering failure.
    #[error(transparent)]
    Coefficient(#[from] CoefficientDomainError),
    /// Path-algebra module or homomorphism failure.
    #[error(transparent)]
    PathAlgebraModule(#[from] PathAlgebraModuleError),
    /// Finite-dimensional quotient detection can require its Gröbner basis.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// A nonzero module must have at least one primitive-idempotent image.
    #[error("a nonzero module unexpectedly has an empty decomposition")]
    EmptyNonzeroDecomposition,
    /// A nonzero endomorphism algebra has a nonempty semisimple central
    /// idempotent set.
    #[error("the semisimple endomorphism top has no central idempotents")]
    EmptySemisimpleEndomorphismTop,
    /// The nonzero branch cannot receive QPA's failed multiplicity result.
    #[error("a nonzero module unexpectedly has no multiplicity decomposition")]
    MissingNonzeroMultiplicityDecomposition,
    /// A nonempty basic-summand list must have a direct sum.
    #[error("the basic-summand direct sum unexpectedly failed")]
    BasicDirectSumFailed,
    /// QPA's alternative top-decomposition method is finite-field-only.
    #[error("the entered module is not over a finite field")]
    FiniteFieldRequired,
    /// `MakeBasicListOfModules` indexes its first entered module.
    #[error("the entered list of modules is empty")]
    EmptyBasicModuleList,
    /// Every module in `MakeBasicListOfModules` must use one acting algebra.
    #[error("not all modules are over the same algebra")]
    BasicModulesOverDifferentAlgebras,
    /// QPA requires a finite-dimensional acting algebra before simplifying
    /// the entered module list.
    #[error("the entered algebra is not finite dimensional")]
    InfiniteDimensionalActingAlgebra,
    /// QPA's final `DecomposeModule(0)` step fails while constructing the
    /// zero module's basis or empty endomorphism algebra.
    #[error("MakeBasicListOfModules cannot decompose a zero module")]
    ZeroModuleInBasicList,
    /// A generated matrix algebra needs at least one nonzero generator.
    #[error("cannot construct the generated matrix algebra from an empty basis")]
    EmptyGeneratedMatrixAlgebra,
    /// Generated-algebra matrices must be equally sized square matrices over
    /// one field.
    #[error("invalid generated matrix-algebra generator")]
    InvalidGeneratedMatrixAlgebraGenerator,
    /// A generated matrix-algebra allocation exceeded addressable size.
    #[error("the generated matrix algebra exceeds addressable size")]
    GeneratedMatrixAlgebraSizeOverflow,
    /// The retained generator span was unexpectedly not closed under matrix
    /// multiplication.
    #[error("a generated matrix-algebra product left the retained span")]
    GeneratedMatrixProductOutsideSpan,
    /// Every induced top endomorphism belongs to its generated algebra.
    #[error("an induced top endomorphism is outside its generated algebra")]
    TopImageOutsideGeneratedAlgebra,
    /// Generated matrix-algebra coordinates carry one coefficient per basis
    /// matrix.
    #[error("generated matrix coordinates have length {actual}, expected {expected}")]
    GeneratedMatrixCoordinateCount {
        /// Entered coordinate count.
        actual: usize,
        /// Retained basis dimension.
        expected: usize,
    },
    /// QPA's probabilistic method cannot choose coordinates from an empty Hom
    /// list.
    #[error("the probabilistic endomorphism list is empty")]
    EmptyProbabilisticEndomorphismList,
    /// The source's logarithmic bound is undefined for the zero module.
    #[error("probabilistic decomposition does not accept the zero module")]
    ZeroProbabilisticModule,
    /// Every entered map must be an endomorphism of the entered module.
    #[error("the probabilistic list contains a map outside End(M)")]
    ProbabilisticEndomorphismModuleMismatch,
    /// The finite-field attempt bound exceeded addressable size.
    #[error("the probabilistic decomposition attempt bound overflows usize")]
    ProbabilisticAttemptBoundOverflow,
    /// A random linear combination needs one coefficient per basis map.
    #[error("the probabilistic coefficient row has length {actual}, expected {expected}")]
    ProbabilisticCoefficientCount {
        /// Entered coefficient count.
        actual: usize,
        /// Required endomorphism count.
        expected: usize,
    },
    /// QPA's Fitting cokernel and kernel must be isomorphic.
    #[error("the probabilistic Fitting factors are unexpectedly nonisomorphic")]
    MissingProbabilisticKernelIsomorphism,
    /// A successful split must have nonzero image and kernel factors.
    #[error("a probabilistic restricted factor is zero")]
    ZeroProbabilisticRestrictedModule,
    /// Restricting the entered endomorphism span must retain an identity.
    #[error("a probabilistic restricted endomorphism basis is empty")]
    EmptyProbabilisticRestrictedBasis,
    /// `HomOverAlgebra(M,M)` must be defined for a module over itself.
    #[error("the characteristic-polynomial endomorphism space is unavailable")]
    MissingCharPolyEndomorphismSpace,
    /// Every nonzero unital module has at least its identity endomorphism.
    #[error("the characteristic-polynomial endomorphism list is empty")]
    EmptyCharPolyEndomorphismList,
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring,
        finite_field::{FiniteFieldCore, FiniteFieldElement, Z2, Zp},
        rational::Q,
    };
    use serde_json::{Value, json};

    use crate::{
        ArrowSpec, DenseMatrix, GaloisField, NamedArrowMatrix, PathAlgebra, PathAlgebraModule,
        PathAlgebraModuleHomomorphism, QpaCoefficientField, Quiver, QuiverAlgebra,
        StructureConstantsAlgebra, VertexIndex, direct_sum_of_qpa_modules, end_over_algebra,
        hom_over_algebra, number_of_nonisomorphic_direct_summands, right_module_over_path_algebra,
    };

    use super::{
        AllSubmodulesOfModuleResult, ModuleDecompositionError, ProbabilisticModuleDecomposition,
        all_simple_submodules_of_module, all_submodules_of_module, basic_version_of_module,
        block_decomposition_of_module, block_splitting_idempotents, complement_in_full_row_space,
        decompose_module, decompose_module_probabilistic, decompose_module_via_char_poly,
        decompose_module_via_top, decompose_module_with_inclusions,
        decompose_module_with_multiplicities, lift_idempotent, lift_idempotents_for_decomposition,
        lift_two_orthogonal_idempotents, make_basic_list_of_modules,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn shifted_upper_triangular_two_by_two(field: Zp) -> StructureConstantsAlgebra<Zp> {
        let zero = field.zero();
        let one = field.one();
        let two = field.nth(2.into());
        let mut products = vec![zero; 27];
        let mut set = |left: usize, right: usize, coordinates: [FiniteFieldElement<u32>; 3]| {
            let start = (left * 3 + right) * 3;
            products[start..start + 3].clone_from_slice(&coordinates);
        };
        set(0, 0, [one, field.zero(), field.zero()]);
        set(0, 1, [field.zero(), field.zero(), two]);
        set(0, 2, [field.zero(), field.zero(), one]);
        set(1, 1, [field.zero(), one, field.zero()]);
        set(2, 1, [field.zero(), field.zero(), one]);
        StructureConstantsAlgebra::new(field, 3, products).unwrap()
    }

    fn integer_rows(field: &Zp, rows: &[Vec<FiniteFieldElement<u32>>]) -> Value {
        json!(
            rows.iter()
                .map(|row| row
                    .iter()
                    .map(|entry| field.from_element(entry))
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>()
        )
    }

    fn matrix_json<F: Field>(
        matrix: &DenseMatrix<F>,
        integer: impl Fn(&F::Element) -> i64 + Copy,
    ) -> Value {
        json!(
            matrix
                .as_flat()
                .chunks(matrix.ncols())
                .map(|row| row.iter().map(integer).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        )
    }

    fn module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        integer: impl Fn(&F::Element) -> i64 + Copy,
    ) -> Value {
        json!({
            "dimension": module.dimension_vector(),
            "matrices": module.arrow_matrices().iter()
                .map(|matrix| matrix_json(matrix, integer)).collect::<Vec<_>>(),
        })
    }

    fn inclusion_json<F: Field>(
        inclusion: &PathAlgebraModuleHomomorphism<F>,
        integer: impl Fn(&F::Element) -> i64 + Copy,
    ) -> Value {
        json!({
            "source_dimension": inclusion.source().dimension_vector(),
            "matrices": inclusion.vertex_matrices().iter()
                .map(|matrix| matrix_json(matrix, integer)).collect::<Vec<_>>(),
        })
    }

    fn all_submodules_json<F: Field>(
        result: &AllSubmodulesOfModuleResult<F>,
        integer: impl Fn(&F::Element) -> i64 + Copy,
    ) -> Value {
        match result {
            AllSubmodulesOfModuleResult::Flat(maps) => json!({
                "shape": "flat",
                "maps": maps.iter().map(|map| inclusion_json(map, integer)).collect::<Vec<_>>(),
            }),
            AllSubmodulesOfModuleResult::ByDimension(levels) => json!({
                "shape": "levels",
                "levels": levels.iter().map(|level| level.iter()
                    .map(|map| inclusion_json(map, integer)).collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
            }),
        }
    }

    fn submodules_case<F: QpaCoefficientField>(
        field: F,
        integer: impl Fn(&F::Element) -> i64 + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let zero = right_module_over_path_algebra(algebra.clone(), vec![0, 0], &[])
            .unwrap()
            .unwrap();
        let simple = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let semisimple = right_module_over_path_algebra(algebra.clone(), vec![2, 1], &[])
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

        let simple_json = |module: &PathAlgebraModule<F>| {
            all_simple_submodules_of_module(module)
                .unwrap()
                .iter()
                .map(|map| inclusion_json(map, integer))
                .collect::<Vec<_>>()
        };
        json!({
            "all_simple": {
                "zero": simple_json(&zero),
                "simple": simple_json(&simple),
                "semisimple": simple_json(&semisimple),
                "projective": simple_json(&projective),
            },
            "all_submodules": {
                "zero": all_submodules_json(&all_submodules_of_module(&zero).unwrap(), integer),
                "simple": all_submodules_json(&all_submodules_of_module(&simple).unwrap(), integer),
                "semisimple": all_submodules_json(
                    &all_submodules_of_module(&semisimple).unwrap(), integer),
                "projective": all_submodules_json(
                    &all_submodules_of_module(&projective).unwrap(), integer),
            },
        })
    }

    fn make_basic_case<F: QpaCoefficientField>(
        field: F,
        integer: impl Fn(&F::Element) -> i64 + Copy,
    ) -> Value {
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
        let projective_plus_simple =
            direct_sum_of_qpa_modules(&[projective.clone(), simple_one.clone()])
                .unwrap()
                .unwrap();
        let ordered = make_basic_list_of_modules(&[
            projective.clone(),
            simple_one.clone(),
            simple_two,
            projective.clone(),
            projective_plus_simple,
        ])
        .unwrap();
        let early_break = make_basic_list_of_modules(&[
            projective.clone(),
            projective,
            simple_one.clone(),
            simple_one,
        ])
        .unwrap();
        json!({
            "ordered_removal": ordered.iter()
                .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
            "early_break": early_break.iter()
                .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
        })
    }

    fn probabilistic_json<F: Field>(
        result: &ProbabilisticModuleDecomposition<F>,
        integer: impl Fn(&F::Element) -> i64 + Copy,
    ) -> Value {
        let (kind, modules) = match result {
            ProbabilisticModuleDecomposition::Unsplit(module) => {
                ("unsplit", std::slice::from_ref(module))
            }
            ProbabilisticModuleDecomposition::Split(modules) => ("split", modules.as_slice()),
        };
        json!({
            "kind": kind,
            "modules": modules.iter()
                .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn moduledecomp_case<F: QpaCoefficientField>(
        field: F,
        integer: impl Fn(&F::Element) -> i64 + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![field.one()], vec![field.zero()]], field.clone())
                    .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let endomorphism_algebra = end_over_algebra(&module).unwrap();
        let (transposed_algebra, endomorphism_basis) = endomorphism_algebra
            .transposed_structure_constants_model()
            .unwrap();
        let transposed_basis = endomorphism_basis
            .iter()
            .map(DenseMatrix::transpose)
            .collect::<Vec<_>>();
        let idempotent_matrices = transposed_algebra
            .idempotents_for_decomposition()
            .unwrap()
            .into_iter()
            .map(|idempotent| {
                let mut matrix = DenseMatrix::zero(3, 3, field.clone());
                for (coefficient, basis) in idempotent.coordinates().iter().zip(&transposed_basis) {
                    matrix = matrix.add(&basis.scale_left(coefficient)).unwrap();
                }
                matrix
            })
            .collect::<Vec<_>>();
        let decomposition = decompose_module(&module).unwrap();
        let inclusions = decompose_module_with_inclusions(&module).unwrap();
        let repeated_decomposition = decompose_module(&module).unwrap();
        let repeated_inclusions = decompose_module_with_inclusions(&module).unwrap();

        let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(algebra, vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let direct_sum = direct_sum_of_qpa_modules(&[simple_one, simple_two])
            .unwrap()
            .unwrap();
        let direct_decomposition = decompose_module(&direct_sum).unwrap();
        let mixed_multiplicities = decompose_module_with_multiplicities(&module)
            .unwrap()
            .unwrap();

        let vector_quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let vector_algebra =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), vector_quiver).unwrap());
        let repeated_module = right_module_over_path_algebra(vector_algebra.clone(), vec![2], &[])
            .unwrap()
            .unwrap();
        let repeated_multiplicities = decompose_module_with_multiplicities(&repeated_module)
            .unwrap()
            .unwrap();
        let repeated_again = decompose_module_with_multiplicities(&repeated_module)
            .unwrap()
            .unwrap();
        let zero_module = right_module_over_path_algebra(vector_algebra.clone(), vec![0], &[])
            .unwrap()
            .unwrap();
        let simple_module = right_module_over_path_algebra(vector_algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        let zero_is_fail = decompose_module_with_multiplicities(&zero_module)
            .unwrap()
            .is_none();
        let block_idempotents = block_splitting_idempotents(&module).unwrap();
        let block_decomposition = block_decomposition_of_module(&module).unwrap();
        let repeated_block_idempotents = block_splitting_idempotents(&repeated_module).unwrap();
        let repeated_block_decomposition = block_decomposition_of_module(&repeated_module).unwrap();
        let zero_idempotents_empty = block_splitting_idempotents(&zero_module)
            .unwrap()
            .is_empty();
        let basic_mixed = basic_version_of_module(&module).unwrap();
        let basic_repeated = basic_version_of_module(&repeated_module).unwrap();
        let basic_zero = basic_version_of_module(&zero_module).unwrap();
        let via_top = decompose_module_via_top(&module).unwrap();
        let repeated_via_top = decompose_module_via_top(&repeated_module).unwrap();
        let repeated_endomorphisms = hom_over_algebra(&repeated_module, &repeated_module)
            .unwrap()
            .unwrap();
        let simple_endomorphisms = hom_over_algebra(&simple_module, &simple_module)
            .unwrap()
            .unwrap();
        let probabilistic_repeated =
            decompose_module_probabilistic(&repeated_endomorphisms, &repeated_module).unwrap();
        let probabilistic_simple =
            decompose_module_probabilistic(&simple_endomorphisms, &simple_module).unwrap();
        let char_poly_repeated = decompose_module_via_char_poly(&repeated_module).unwrap();
        let char_poly_simple = decompose_module_via_char_poly(&simple_module).unwrap();
        let noniso_mixed = number_of_nonisomorphic_direct_summands(&module).unwrap();
        let noniso_repeated = number_of_nonisomorphic_direct_summands(&repeated_module).unwrap();

        json!({
            "algorithm": {
                "endomorphism_basis": endomorphism_algebra.basis_matrices().iter()
                    .map(|matrix| matrix_json(matrix, integer)).collect::<Vec<_>>(),
                "transposed_basis": transposed_basis.iter()
                    .map(|matrix| matrix_json(matrix, integer)).collect::<Vec<_>>(),
                "idempotents": idempotent_matrices.iter()
                    .map(|matrix| matrix_json(matrix, integer)).collect::<Vec<_>>(),
            },
            "decomposition": decomposition.iter()
                .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
            "inclusions": inclusions.iter()
                .map(|inclusion| inclusion_json(inclusion, integer)).collect::<Vec<_>>(),
            "cache": {
                "decomposition_first": decomposition[0]
                    .is_identical(&repeated_decomposition[0]),
                "inclusion_first": inclusions[0].is_identical(&repeated_inclusions[0]),
            },
            "direct_sum": direct_decomposition.iter()
                .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
            "multiplicities": {
                "mixed": {
                    "basic": mixed_multiplicities.basic_summands().iter()
                        .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
                    "values": mixed_multiplicities.multiplicities(),
                },
                "repeated": {
                    "basic": repeated_multiplicities.basic_summands().iter()
                        .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
                    "values": repeated_multiplicities.multiplicities(),
                    "cache_first": repeated_multiplicities.basic_summands()[0]
                        .is_identical(&repeated_again.basic_summands()[0]),
                },
                "zero_is_fail": zero_is_fail,
            },
            "blocks": {
                "mixed_idempotents": block_idempotents.iter().map(|homomorphism| {
                    homomorphism.vertex_matrices().iter()
                        .map(|matrix| matrix_json(matrix, integer)).collect::<Vec<_>>()
                }).collect::<Vec<_>>(),
                "mixed_modules": block_decomposition.iter()
                    .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
                "repeated_idempotents": repeated_block_idempotents.iter().map(|homomorphism| {
                    homomorphism.vertex_matrices().iter()
                        .map(|matrix| matrix_json(matrix, integer)).collect::<Vec<_>>()
                }).collect::<Vec<_>>(),
                "repeated_modules": repeated_block_decomposition.iter()
                    .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
                "zero_idempotents_empty": zero_idempotents_empty,
            },
            "basic": {
                "mixed": [module_json(&basic_mixed, integer)],
                "repeated": [module_json(&basic_repeated, integer)],
                "zero_identical": basic_zero.is_identical(&zero_module),
            },
            "via_top": {
                "mixed": via_top.iter()
                    .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
                "repeated": repeated_via_top.iter()
                    .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
            },
            "probabilistic": {
                "repeated": probabilistic_json(&probabilistic_repeated, integer),
                "simple": probabilistic_json(&probabilistic_simple, integer),
            },
            "char_poly": {
                "repeated": char_poly_repeated.iter()
                    .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
                "simple": char_poly_simple.iter()
                    .map(|module| module_json(module, integer)).collect::<Vec<_>>(),
            },
            "noniso": {
                "mixed": [json!(noniso_mixed.count()), json!(noniso_mixed.block_dimensions())],
                "repeated": [
                    json!(noniso_repeated.count()),
                    json!(noniso_repeated.block_dimensions()),
                ],
            },
        })
    }

    #[test]
    fn pinned_qpa_moduledecomp_foundation_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/moduledecomp_foundation.schema-v1.json"
        ))
        .unwrap();

        let rational_basis = DenseMatrix::from_rows(
            vec![
                vec![Q.nth(1.into()), Q.nth(2.into()), Q.zero(), Q.nth(1.into())],
                vec![Q.zero(), Q.zero(), Q.nth(1.into()), Q.nth(3.into())],
            ],
            Q,
        )
        .unwrap();
        let rational_complement = complement_in_full_row_space(&rational_basis).unwrap();
        let binary_basis = DenseMatrix::from_rows(
            vec![
                vec![Z2.one(), Z2.one(), Z2.zero(), Z2.one()],
                vec![Z2.zero(), Z2.zero(), Z2.one(), Z2.one()],
            ],
            Z2,
        )
        .unwrap();
        let binary_complement = complement_in_full_row_space(&binary_basis).unwrap();

        let f5 = Zp::new(5);
        let full_complement =
            complement_in_full_row_space(&DenseMatrix::identity(3, f5.clone())).unwrap();
        let algebra = shifted_upper_triangular_two_by_two(f5.clone());
        let ideal =
            DenseMatrix::from_rows(vec![vec![f5.zero(), f5.zero(), f5.one()]], f5.clone()).unwrap();
        let quotient = algebra.quotient_by_ideal_basis(ideal).unwrap();
        let quotient_ids = vec![vec![f5.zero(), f5.one()], vec![f5.one(), f5.zero()]];
        let mapping = quotient.natural_homomorphism().unwrap();
        let lifted = lift_idempotents_for_decomposition(&mapping, &quotient_ids).unwrap();
        assert_eq!(
            quotient
                .lift_idempotents_for_decomposition(&quotient_ids)
                .unwrap(),
            lifted
        );
        let primitive = algebra
            .idempotents_for_decomposition()
            .unwrap()
            .into_iter()
            .map(|element| element.coordinates().to_vec())
            .collect::<Vec<_>>();
        let first = lift_idempotent(&mapping, &quotient_ids[0]).unwrap();
        assert_eq!(quotient.lift_idempotent(&quotient_ids[0]).unwrap(), first);
        let pair = lift_two_orthogonal_idempotents(&mapping, &first, &quotient_ids[1]).unwrap();
        assert_eq!(
            quotient
                .lift_two_orthogonal_idempotents(&first, &quotient_ids[1])
                .unwrap(),
            pair
        );

        let actual = json!({
            "schema_version": 1,
            "kind": "moduledecomp_foundation",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "complement": {
                "rational": rational_complement.as_flat().chunks(4)
                    .map(|row| row.iter().map(ToString::to_string)
                        .map(|entry| entry.parse::<i64>().unwrap()).collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
                "f2": binary_complement.as_flat().chunks(4)
                    .map(<[u8]>::to_vec)
                    .collect::<Vec<_>>(),
                "full_f5": integer_rows(&f5, &[]),
            },
            "f5": {
                "quotient_ids": integer_rows(&f5, &quotient_ids),
                "lifted_for_decomposition": integer_rows(&f5, &lifted),
                "idempotents_for_decomposition": integer_rows(&f5, &primitive),
                "lift_idempotent": integer_rows(&f5, &[first]),
                "lift_two": integer_rows(&f5, &pair),
            }
        });
        assert_eq!(actual, expected);
        assert_eq!(full_complement.nrows(), 0);
        assert_eq!(full_complement.ncols(), 3);
    }

    #[test]
    fn complement_rejects_qpa_trivial_space() {
        assert!(matches!(
            complement_in_full_row_space(&DenseMatrix::zero(0, 3, Q)),
            Err(ModuleDecompositionError::TrivialRowSpace)
        ));
    }

    #[test]
    fn pinned_qpa_moduledecomp_modules_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/moduledecomp_modules.schema-v1.json"
        ))
        .unwrap();
        let f5 = Zp::new(5);
        let actual = json!({
            "schema_version": 1,
            "kind": "moduledecomp_modules",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654",
            },
            "binary": moduledecomp_case(Z2, |entry| i64::from(*entry)),
            "quinary": moduledecomp_case(f5.clone(), |entry| i64::from(f5.from_element(entry))),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn extension_fields_run_qpa_decomposition_and_enumeration_algorithms() {
        for (prime, degree, order) in [(2, 2, 4usize), (2, 3, 8), (3, 2, 9)] {
            let field = GaloisField::conway(prime, degree).unwrap();
            let quiver =
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
            let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
            let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
                .unwrap()
                .unwrap();
            let projective_one = right_module_over_path_algebra(
                algebra.clone(),
                vec![1, 1],
                &[NamedArrowMatrix::new(
                    "a",
                    DenseMatrix::from_rows(vec![vec![field.one()]], field.clone()).unwrap(),
                )],
            )
            .unwrap()
            .unwrap();
            let mixed = direct_sum_of_qpa_modules(&[simple_one.clone(), projective_one])
                .unwrap()
                .unwrap();
            assert_eq!(
                decompose_module(&mixed)
                    .unwrap()
                    .iter()
                    .map(PathAlgebraModule::dimension_vector)
                    .collect::<Vec<_>>(),
                [vec![1, 0], vec![1, 1]]
            );

            let repeated = direct_sum_of_qpa_modules(&[simple_one.clone(), simple_one])
                .unwrap()
                .unwrap();
            assert_eq!(
                decompose_module_with_multiplicities(&repeated)
                    .unwrap()
                    .unwrap()
                    .multiplicities(),
                &[2]
            );
            assert_eq!(
                all_simple_submodules_of_module(&repeated).unwrap().len(),
                order + 1
            );

            let semisimple = right_module_over_path_algebra(algebra, vec![2, 1], &[])
                .unwrap()
                .unwrap();
            let AllSubmodulesOfModuleResult::ByDimension(levels) =
                all_submodules_of_module(&semisimple).unwrap()
            else {
                panic!("a nonzero finite-field module must use dimension levels");
            };
            assert_eq!(
                levels.iter().map(Vec::len).collect::<Vec<_>>(),
                [1, order + 2, order + 2, 1]
            );
        }
    }

    #[test]
    fn alternative_decompositions_preserve_qpa_finite_field_precondition() {
        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let module = right_module_over_path_algebra(algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        assert!(matches!(
            decompose_module_via_top(&module),
            Err(ModuleDecompositionError::FiniteFieldRequired)
        ));
        assert!(matches!(
            decompose_module_via_char_poly(&module),
            Err(ModuleDecompositionError::FiniteFieldRequired)
        ));
    }

    #[test]
    fn pinned_qpa_submodule_enumeration_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_submodules.schema-v1.json"
        ))
        .unwrap();
        let f5 = Zp::new(5);
        let actual = json!({
            "schema_version": 1,
            "kind": "modulehom_submodules",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654",
            },
            "binary": submodules_case(Z2, |entry| i64::from(*entry)),
            "quinary": submodules_case(
                f5.clone(), |entry| i64::from(f5.from_element(entry))),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn submodule_enumeration_preserves_qpa_finite_field_precondition() {
        let quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let simple = right_module_over_path_algebra(algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        assert!(matches!(
            all_simple_submodules_of_module(&simple),
            Err(ModuleDecompositionError::FiniteFieldRequired)
        ));
        assert!(matches!(
            all_submodules_of_module(&simple),
            Err(ModuleDecompositionError::FiniteFieldRequired)
        ));
    }

    #[test]
    fn pinned_qpa_make_basic_list_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_make_basic.schema-v1.json"
        ))
        .unwrap();
        let f5 = Zp::new(5);
        let actual = json!({
            "schema_version": 1,
            "kind": "modulehom_make_basic",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654",
            },
            "binary": make_basic_case(Z2, |entry| i64::from(*entry)),
            "quinary": make_basic_case(
                f5.clone(), |entry| i64::from(f5.from_element(entry))),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn make_basic_list_preserves_qpa_input_failures() {
        assert!(matches!(
            make_basic_list_of_modules::<numerica::domains::finite_field::Z2>(&[]),
            Err(ModuleDecompositionError::EmptyBasicModuleList)
        ));

        let rational_quiver = Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap();
        let rational_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, rational_quiver).unwrap());
        let rational = right_module_over_path_algebra(rational_algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        assert!(matches!(
            make_basic_list_of_modules(&[rational]),
            Err(ModuleDecompositionError::FiniteFieldRequired)
        ));

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, loop_quiver).unwrap());
        let loop_module = right_module_over_path_algebra(loop_algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        assert!(matches!(
            make_basic_list_of_modules(&[loop_module]),
            Err(ModuleDecompositionError::InfiniteDimensionalActingAlgebra)
        ));

        let first_algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
            )
            .unwrap(),
        );
        let second_algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(2, Vec::<ArrowSpec>::new()).unwrap(),
            )
            .unwrap(),
        );
        let first = right_module_over_path_algebra(first_algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        let second = right_module_over_path_algebra(second_algebra, vec![1, 0], &[])
            .unwrap()
            .unwrap();
        assert!(matches!(
            make_basic_list_of_modules(&[first, second]),
            Err(ModuleDecompositionError::BasicModulesOverDifferentAlgebras)
        ));

        let zero_algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
            )
            .unwrap(),
        );
        let zero = right_module_over_path_algebra(zero_algebra, vec![0], &[])
            .unwrap()
            .unwrap();
        assert!(matches!(
            make_basic_list_of_modules(&[zero]),
            Err(ModuleDecompositionError::ZeroModuleInBasicList)
        ));
    }
}
