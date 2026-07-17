//! QPA's finite-dimensional module-decomposition foundations.
//!
//! This module begins the source-faithful translation of
//! `lib/moduledecomp.gi`: row-space complements and the idempotent algorithms
//! used by the deterministic module decomposition methods.

use std::cmp::Ordering;

use numerica::domains::{Field, integer::Integer};
use thiserror::Error;

use crate::{
    CoefficientDomainError, CommonDirectSummandResult, DenseMatrix, HopfError, MatrixError,
    PathAlgebraModule, PathAlgebraModuleError, PathAlgebraModuleHomomorphism, QpaCoefficientField,
    StructureConstantsAlgebra, StructureConstantsAlgebraElement, StructureConstantsError,
    StructureConstantsQuotient, common_direct_summand, direct_sum_of_qpa_modules, end_over_algebra,
    from_end_m_to_hom_mm, image_inclusion, image_of_module_homomorphism,
    path_module::ModuleDecompositionMultiplicitiesCache,
    structure_constants::{StructureConstantsSubalgebra, add_rows, span_rows, subtract_rows},
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
        ArrowSpec, DenseMatrix, NamedArrowMatrix, PathAlgebra, PathAlgebraModule,
        PathAlgebraModuleHomomorphism, QpaCoefficientField, Quiver, QuiverAlgebra,
        StructureConstantsAlgebra, VertexIndex, direct_sum_of_qpa_modules, end_over_algebra,
        right_module_over_path_algebra,
    };

    use super::{
        ModuleDecompositionError, basic_version_of_module, block_decomposition_of_module,
        block_splitting_idempotents, complement_in_full_row_space, decompose_module,
        decompose_module_with_inclusions, decompose_module_with_multiplicities, lift_idempotent,
        lift_idempotents_for_decomposition, lift_two_orthogonal_idempotents,
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
        let zero_module = right_module_over_path_algebra(vector_algebra, vec![0], &[])
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
}
