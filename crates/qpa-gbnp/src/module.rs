//! QPA modules constructed from enveloping algebras.
//!
//! This begins the source-order translation of `lib/pathalgtensor.gi`'s
//! module section. The regular bimodule method deliberately uses QPA's
//! canonical-basis partition and its checked path-module constructor.

use std::num::NonZeroUsize;

use numerica::domains::{
    Field, Ring, RingOps,
    finite_field::{Z2, Zp},
    integer::Integer,
    rational::{Q, Rational, RationalField},
};
use qpa_core::{
    AlgebraGlobalDimension, CachedAttribute, CachedNakayamaAutomorphismOrder,
    CoefficientDomainError, CommonDirectSummandResult, DenseMatrix, MatrixError,
    ModuleDecompositionError, ModuleDimensionSeries, ModuleEndomorphismQuotientHomomorphism,
    NamedArrowMatrix, PathAlgebraElement, PathAlgebraError, PathAlgebraModule,
    PathAlgebraModuleElement, PathAlgebraModuleError, PathAlgebraModuleHomomorphism,
    PathAlgebraModuleHomomorphismPair, PathAlgebraTensorError, PathEndpoint, PathError,
    QpaCoefficientField, QpaFieldCharacteristic, QuiverAlgebra, QuiverAlgebraElement,
    QuiverAlgebraElementError, QuiverProduct, QuiverProductError, QuotientError,
    QuotientPathAlgebra, StructureConstantsError, cokernel_of_module_homomorphism,
    cokernel_projection, common_direct_summand, decompose_module_with_inclusions,
    decompose_module_with_multiplicities, direct_sum_inclusions, direct_sum_of_qpa_modules,
    direct_sum_projections, end_over_algebra, from_end_m_to_hom_mm, from_hom_mm_to_end_m,
    hom_over_algebra, identity_module_mapping, include_in_path_algebra, is_indecomposable_module,
    isomorphic_modules, isomorphism_of_modules, kernel_inclusion, kernel_of_module_homomorphism,
    left_minimal_version, loewy_length_of_module, minimal_generating_set_of_module,
    morphism_on_cokernel, morphism_on_kernel, natural_endomorphism_quotient_by_ideal,
    project_from_product_quiver, push_out, radical_of_module, radical_of_module_inclusion,
    right_minimal_version, right_module_hom_over_algebra, right_module_over_path_algebra,
    right_module_over_path_algebra_with_tensor, socle_of_module, socle_of_module_inclusion,
    subrepresentation_inclusion, top_of_module, zero_module, zero_module_mapping,
};
use rand::Rng;
use thiserror::Error;

use crate::complex::projective_resolution;
use crate::{
    PathalgGbnpError, TensorProductError, ensure_quotient_groebner_basis, enveloping_algebra,
    opposite_path_algebra, opposite_quotient_element, quotient_is_admissible,
    quotient_is_finite_dimensional,
};

/// Translate QPA's cached `AlgebraAsModuleOverEnvelopingAlgebra(A)`
/// attribute (and its `AlgebraAsModuleOfEnvelopingAlgebra` synonym).
///
/// The canonical basis is partitioned into the spaces `e_i A e_j`. For each
/// product-quiver arrow, the source method acts on those basis vectors by the
/// opposite of its first projection on the left and its second projection on
/// the right, then selects canonical-basis coordinates in the target space.
/// The resulting row-action matrices are passed through QPA's checked
/// `RightModuleOverPathAlgebra` translation.
///
/// `Ok(None)` is QPA's `fail` for an infinite-dimensional algebra (or the
/// checked constructor's theoretically possible relation failure).
///
/// # Errors
///
/// Propagates the configured GBNP finite-dimensionality/basis path,
/// enveloping-algebra construction, product projections, exact algebra and
/// matrix arithmetic, and checked module construction.
pub fn algebra_as_module_over_enveloping_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_regular_enveloping_module() {
        return Ok(Some(cached));
    }
    if algebra.regular_enveloping_module_attribute_is_fail() {
        return Ok(None);
    }

    let Some(basis_vectors) = finite_canonical_basis(algebra)? else {
        algebra.mark_regular_enveloping_module_attribute_fail();
        return Ok(None);
    };

    let enveloping = enveloping_algebra(algebra)?;
    let product = enveloping.product_quiver();
    let vertices = algebra.quiver().vertices().collect::<Vec<_>>();
    let mut vector_spaces = vec![vec![Vec::<usize>::new(); vertices.len()]; vertices.len()];

    let vertex_elements = vertices
        .iter()
        .map(|vertex| include_in_path_algebra(&vertex.path(), algebra))
        .collect::<Result<Vec<_>, PathAlgebraTensorError>>()?;
    for (left_index, left) in vertex_elements.iter().enumerate() {
        for (right_index, right) in vertex_elements.iter().enumerate() {
            for (basis_index, basis_vector) in basis_vectors.iter().enumerate() {
                let component = left.multiply(basis_vector)?.multiply(right)?;
                if !component.is_zero() {
                    vector_spaces[left_index][right_index].push(basis_index);
                }
            }
        }
    }
    let dimension_vector = vector_spaces
        .iter()
        .flat_map(|row| row.iter().map(Vec::len))
        .collect::<Vec<_>>();

    let mut module_specification = Vec::new();
    for arrow in enveloping.algebra().quiver().arrows() {
        let components = [
            project_from_product_quiver(1, &arrow.path(), product)?,
            project_from_product_quiver(2, &arrow.path(), product)?,
        ];
        let source_indices = components.each_ref().map(|component| {
            let PathEndpoint::Vertex(vertex) = component.source() else {
                unreachable!("a nonzero projected product arrow has a vertex source")
            };
            vertex.index().get() - 1
        });
        let target_indices = components.each_ref().map(|component| {
            let PathEndpoint::Vertex(vertex) = component.target() else {
                unreachable!("a nonzero projected product arrow has a vertex target")
            };
            vertex.index().get() - 1
        });
        let source_space = &vector_spaces[source_indices[0]][source_indices[1]];
        let target_space = &vector_spaces[target_indices[0]][target_indices[1]];
        if source_space.is_empty() || target_space.is_empty() {
            continue;
        }

        let left_path = components[0].quiver().opposite_path(&components[0])?;
        let left = include_in_path_algebra(&left_path, algebra)?;
        let right = include_in_path_algebra(&components[1], algebra)?;
        let rows = source_space
            .iter()
            .map(|&basis_index| {
                let image = left
                    .multiply(&basis_vectors[basis_index])?
                    .multiply(&right)?;
                let coefficients = canonical_basis_coefficients(algebra, &basis_vectors, &image)?;
                Ok(target_space
                    .iter()
                    .map(|&target_index| coefficients[target_index].clone())
                    .collect::<Vec<_>>())
            })
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
        let matrix = DenseMatrix::from_rows(rows, algebra.field().clone())?;
        module_specification.push(NamedArrowMatrix::new(arrow.name(), matrix));
    }

    let Some(module) = right_module_over_path_algebra_with_tensor(
        enveloping.algebra().clone(),
        dimension_vector,
        &module_specification,
        enveloping,
    )?
    else {
        algebra.mark_regular_enveloping_module_attribute_fail();
        return Ok(None);
    };
    Ok(Some(algebra.cache_regular_enveloping_module(&module)))
}

/// QPA's spelling synonym for
/// [`algebra_as_module_over_enveloping_algebra`].
///
/// # Errors
///
/// Propagates the identical source method failures.
pub fn algebra_as_module_of_enveloping_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    algebra_as_module_over_enveloping_algebra(algebra)
}

/// Translate QPA's algebra `LoewyLength` method through its indecomposable
/// projective representations and repeated `RadicalOfModule` calls.
///
/// # Errors
///
/// Requires a finite path algebra or finite admissible quotient and
/// propagates canonical-basis, exact module, and radical failures.
pub fn loewy_length_of_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<usize, EnvelopingModuleError> {
    let projectives = indec_projective_modules(algebra)?;
    let mut maximum = 0;
    for projective in projectives {
        maximum = maximum.max(loewy_length_of_module(projective)?);
    }
    Ok(maximum)
}

/// Vertexwise canonical bases returned by QPA's `BasisOfProjectives`.
///
/// The outer index is the generating vertex and the inner index is the right
/// support vertex; elements retain QPA's canonical right-ideal basis order.
pub type ProjectiveModuleBases<F> = Vec<Vec<Vec<QuiverAlgebraElement<F>>>>;

struct ProjectiveBasisData<F: Field> {
    basis: Vec<QuiverAlgebraElement<F>>,
    intervals: Vec<Vec<Vec<usize>>>,
}

/// Translate QPA's two `BasisOfProjectives(A)` methods.
///
/// `Ok(None)` is the printed `fail` result for an infinite path algebra or a
/// quotient without finitely many nontips.
///
/// # Errors
///
/// Propagates the configured GBNP finite-basis path.
pub fn basis_of_projectives<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<ProjectiveModuleBases<F>>, EnvelopingModuleError> {
    let Some(data) = projective_basis_data(algebra)? else {
        return Ok(None);
    };
    Ok(Some(
        data.intervals
            .iter()
            .map(|projective| {
                projective
                    .iter()
                    .map(|interval| {
                        interval
                            .iter()
                            .map(|&index| data.basis[index].clone())
                            .collect()
                    })
                    .collect()
            })
            .collect(),
    ))
}

/// Translate QPA's cached `IndecProjectiveModules(A)` attribute.
///
/// Modules are returned in vertex order, with the exact canonical-right-ideal
/// action matrices and QPA's distinguished vertex generator installed.
///
/// # Errors
///
/// Requires a finite path algebra or finite admissible quotient and propagates
/// canonical-basis, exact multiplication, coordinate, and checked-module
/// construction failures.
pub fn indec_projective_modules<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Vec<PathAlgebraModule<F>>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_indec_projective_modules() {
        return Ok(cached);
    }
    let Some(data) = projective_basis_data(algebra)? else {
        return Err(EnvelopingModuleError::InfiniteLoewyAlgebra);
    };
    if matches!(algebra, QuiverAlgebra::Quotient(quotient) if !quotient.is_admissible()?) {
        return Err(EnvelopingModuleError::NonadmissibleLoewyAlgebra);
    }
    let basis = data.basis;
    let vertex_count = algebra.quiver().number_of_vertices();
    let arrows = algebra
        .quiver()
        .arrows()
        .map(|arrow| include_in_path_algebra(&arrow.path(), algebra))
        .collect::<Result<Vec<_>, _>>()?;
    let mut projectives = Vec::with_capacity(vertex_count);

    for projective_vertex in 0..vertex_count {
        let intervals = &data.intervals[projective_vertex];
        let dimensions = intervals.iter().map(Vec::len).collect::<Vec<_>>();
        let mut specification = Vec::new();
        for (arrow, arrow_element) in algebra.quiver().arrows().zip(&arrows) {
            let source = arrow.source().index().get() - 1;
            let target = arrow.target().index().get() - 1;
            if dimensions[source] == 0 || dimensions[target] == 0 {
                continue;
            }
            let rows = intervals[source]
                .iter()
                .map(|&basis_index| {
                    let image = basis[basis_index].multiply(arrow_element)?;
                    let coefficients = canonical_basis_coefficients(algebra, &basis, &image)?;
                    Ok(intervals[target]
                        .iter()
                        .map(|&index| coefficients[index].clone())
                        .collect::<Vec<_>>())
                })
                .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
            specification.push(NamedArrowMatrix::new(
                arrow.name(),
                DenseMatrix::from_rows(rows, algebra.field().clone())?,
            ));
        }
        let projective =
            right_module_over_path_algebra(algebra.clone(), dimensions.clone(), &specification)?
                .ok_or(EnvelopingModuleError::GeneratedProjectiveViolatesRelations)?;
        let mut generator_coordinates = dimensions
            .iter()
            .map(|&dimension| vec![algebra.field().zero(); dimension])
            .collect::<Vec<_>>();
        let generator_coordinate = intervals[projective_vertex]
            .iter()
            .position(|&basis_index| basis[basis_index].leading_monomial().length() == Some(0))
            .ok_or(EnvelopingModuleError::MissingProjectiveVertexGenerator {
                vertex: projective_vertex + 1,
            })?;
        generator_coordinates[projective_vertex][generator_coordinate] = algebra.field().one();
        let generator = projective.element_from_vertex_coordinates(generator_coordinates)?;
        projective.cache_minimal_generating_set(std::slice::from_ref(&generator))?;
        let _ = projective.cache_is_indecomposable_module(true);
        let _ = projective.cache_is_projective_module(true);
        projectives.push(projective);
    }
    Ok(algebra.merge_indec_projective_modules(&projectives))
}

/// Translate QPA's cached `IndecInjectiveModules(A)` attribute.
///
/// The method literally constructs the opposite algebra's ordered
/// indecomposable projectives and applies [`dual_of_module`] to each one.
/// Dual property propagation installs both indecomposability and injectivity
/// on every result, exactly as in QPA.
///
/// # Errors
///
/// Propagates opposite-algebra, projective-family, and module-dual failures.
pub fn indec_injective_modules<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Vec<PathAlgebraModule<F>>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_indec_injective_modules() {
        return Ok(cached);
    }
    let opposite = match algebra {
        QuiverAlgebra::Path(algebra) => QuiverAlgebra::Path(algebra.opposite()?),
        QuiverAlgebra::Quotient(algebra) => {
            QuiverAlgebra::Quotient(opposite_path_algebra(algebra)?)
        }
    };
    let opposite_projectives = indec_projective_modules(&opposite)?;
    let injectives = opposite_projectives
        .iter()
        .map(dual_of_module)
        .collect::<Result<Vec<_>, ModuleFunctorError>>()?;
    Ok(algebra.merge_indec_injective_modules(&injectives))
}

/// Numerica fields for which GAP's `Random(field)` scalar distribution is
/// translated by this port.
///
/// Implementations currently cover `QQ`, `GF(2)`, and Numerica's runtime
/// `GF(p)` backend. Extension-field backends can implement this contract when
/// they are added without changing QPA's `RandomModule` algorithm.
pub trait QpaRandomScalar: Field {
    /// Draw one scalar with GAP's `Random(field)` distribution.
    fn qpa_random_scalar<R: Rng>(&self, rng: &mut R) -> Self::Element;
}

fn qpa_random_integer<R: Rng>(rng: &mut R) -> i64 {
    i64::from(rng.random_range(0_u32..(1_u32 << 20)).count_ones()) - 10
}

impl QpaRandomScalar for RationalField {
    fn qpa_random_scalar<R: Rng>(&self, rng: &mut R) -> Self::Element {
        let denominator = loop {
            let value = qpa_random_integer(rng);
            if value != 0 {
                break value;
            }
        };
        let numerator = qpa_random_integer(rng);
        self.div(
            &self.nth(Integer::from(numerator)),
            &self.nth(Integer::from(denominator)),
        )
    }
}

impl QpaRandomScalar for Z2 {
    fn qpa_random_scalar<R: Rng>(&self, rng: &mut R) -> Self::Element {
        self.sample(rng, (0, 2))
    }
}

impl QpaRandomScalar for Zp {
    fn qpa_random_scalar<R: Rng>(&self, rng: &mut R) -> Self::Element {
        self.sample(rng, (0, i64::MAX))
    }
}

/// Translate QPA's `RandomModule(A, m, n)` method.
///
/// Each of at most ten attempts independently chooses `m` and `n`
/// indecomposable projective summands, computes QPA's ordered Hom basis,
/// scales every basis map by `Random(LeftActingDomain(A))`, sums the maps, and
/// returns its cokernel as soon as that Hom space is nonzero. `Ok(None)` is
/// QPA's `fail` after ten zero-Hom attempts.
///
/// Rust receives the random source explicitly. Its engine state is not meant
/// to be byte-identical to GAP's global Mersenne Twister; projective choices
/// and scalar distributions are otherwise the source method's ones.
///
/// # Errors
///
/// Propagates projective construction, direct sums, Hom, exact map arithmetic,
/// and cokernel construction. The two internal errors represent impossible
/// failures for the nonempty, common-algebra projective lists built here.
pub fn random_module<F: QpaRandomScalar, R: Rng>(
    algebra: &QuiverAlgebra<F>,
    source_summands: NonZeroUsize,
    range_summands: NonZeroUsize,
    rng: &mut R,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    let projectives = indec_projective_modules(algebra)?;
    for _ in 0..10 {
        let source_parts = (0..source_summands.get())
            .map(|_| projectives[rng.random_range(0..projectives.len())].clone())
            .collect::<Vec<_>>();
        let range_parts = (0..range_summands.get())
            .map(|_| projectives[rng.random_range(0..projectives.len())].clone())
            .collect::<Vec<_>>();
        let source = direct_sum_of_qpa_modules(&source_parts)?
            .ok_or(EnvelopingModuleError::InternalRandomModuleDirectSumFailed)?;
        let range = direct_sum_of_qpa_modules(&range_parts)?
            .ok_or(EnvelopingModuleError::InternalRandomModuleDirectSumFailed)?;
        let basis = hom_over_algebra(&source, &range)?
            .ok_or(EnvelopingModuleError::InternalRandomModuleHomRejected)?;
        if let Some((first, rest)) = basis.split_first() {
            let mut random_map = first.scale(&algebra.field().qpa_random_scalar(rng))?;
            for homomorphism in rest {
                random_map = random_map
                    .add(&homomorphism.scale(&algebra.field().qpa_random_scalar(rng))?)?;
            }
            return Ok(Some(cokernel_of_module_homomorphism(&random_map)?));
        }
    }
    Ok(None)
}

/// Translate QPA's cached `SimpleModules(A)` attribute.
///
/// One vertex-simple representation is constructed per vertex in quiver
/// order. Each has a single one in its dimension vector, no entered arrow
/// matrices, and QPA's `IsSimpleQPAModule` property installed together with
/// its true-method implications.
///
/// # Errors
///
/// Requires a finite path algebra or finite admissible quotient and
/// propagates the checked module constructor.
pub fn simple_modules<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Vec<PathAlgebraModule<F>>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_simple_modules() {
        return Ok(cached);
    }
    if finite_canonical_basis(algebra)?.is_none() {
        return Err(EnvelopingModuleError::InfiniteSimpleModuleAlgebra);
    }
    if matches!(algebra, QuiverAlgebra::Quotient(quotient) if !quotient.is_admissible()?) {
        return Err(EnvelopingModuleError::NonadmissibleSimpleModuleAlgebra);
    }
    let vertex_count = algebra.quiver().number_of_vertices();
    let mut simples = Vec::with_capacity(vertex_count);
    for vertex in 0..vertex_count {
        let mut dimensions = vec![0; vertex_count];
        dimensions[vertex] = 1;
        let simple = right_module_over_path_algebra(algebra.clone(), dimensions, &[])?
            .ok_or(EnvelopingModuleError::GeneratedSimpleViolatesRelations)?;
        let _ = simple.cache_is_simple_qpa_module(true);
        simples.push(simple);
    }
    Ok(algebra.merge_simple_modules(&simples))
}

/// QPA's result type for `GlobalDimensionOfAlgebra(A, n)`.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum GlobalDimensionOfAlgebraResult {
    /// QPA found and cached this finite global dimension.
    Finite(usize),
    /// QPA's loop test returned `infinity`.
    Infinite,
    /// QPA returned `false` because some simple still has a nonzero projective
    /// term in degree `n + 1`.
    ExceedsBound,
}

/// QPA's three-valued bounded dominant-dimension result.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum DominantDimensionResult {
    /// QPA found this finite dominant dimension within the entered bound.
    Finite(usize),
    /// The module is projective-injective, or the algebra is selfinjective.
    Infinite,
    /// QPA returned `false` because the dominant dimension exceeds the bound.
    ExceedsBound,
}

/// Translate QPA's bounded `GlobalDimensionOfAlgebra(A, n)` algorithm.
///
/// `None` preserves `TryNextMethod()` for infinite-dimensional input. A path
/// algebra already has QPA's constructor-installed global dimension (zero for
/// a discrete quiver and one otherwise), so that cached value is returned
/// before consulting the bound. For a finite quotient, a quiver loop gives
/// `infinity`; otherwise the method follows the projective resolutions of all
/// vertex simples through degree `n + 1` and returns QPA's integer or `false`.
///
/// # Errors
///
/// Propagates finite-basis, simple-module, projective-cover, syzygy, and exact
/// module failures, and rejects an overflowing `n + 1` bound.
pub fn global_dimension_of_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
    bound: usize,
) -> Result<Option<GlobalDimensionOfAlgebraResult>, EnvelopingModuleError> {
    if finite_canonical_basis(algebra)?.is_none() {
        return Ok(None);
    }
    if let QuiverAlgebra::Path(path_algebra) = algebra {
        return Ok(Some(GlobalDimensionOfAlgebraResult::Finite(
            path_algebra.qpa_global_dimension(),
        )));
    }
    let QuiverAlgebra::Quotient(quotient) = algebra else {
        unreachable!("every quiver algebra is a path algebra or quotient")
    };
    if let Some(cached) = quotient.cached_global_dimension() {
        return Ok(Some(match cached {
            AlgebraGlobalDimension::Finite(dimension) => {
                GlobalDimensionOfAlgebraResult::Finite(dimension)
            }
            AlgebraGlobalDimension::Infinite => GlobalDimensionOfAlgebraResult::Infinite,
        }));
    }
    if algebra
        .quiver()
        .arrows()
        .any(|arrow| arrow.source() == arrow.target())
    {
        let _ = quotient.cache_global_dimension(AlgebraGlobalDimension::Infinite);
        return Ok(Some(GlobalDimensionOfAlgebraResult::Infinite));
    }

    let cutoff = bound
        .checked_add(1)
        .ok_or(EnvelopingModuleError::GlobalDimensionBoundOverflow { bound })?;
    let mut dimension = 0usize;
    for simple in simple_modules(algebra)? {
        let mut syzygy = simple;
        let mut degree = 0usize;
        while degree < cutoff {
            let projective = projective_cover(&syzygy)?;
            if projective.source().dimension() == 0 {
                break;
            }
            syzygy = first_syzygy(&syzygy)?;
            degree += 1;
        }
        if degree < cutoff {
            dimension = dimension.max(degree.saturating_sub(1));
            continue;
        }
        let projective = projective_cover(&syzygy)?;
        if projective.source().dimension() != 0 {
            return Ok(Some(GlobalDimensionOfAlgebraResult::ExceedsBound));
        }
        dimension = dimension.max(bound);
    }
    let _ = quotient.cache_global_dimension(AlgebraGlobalDimension::Finite(dimension));
    Ok(Some(GlobalDimensionOfAlgebraResult::Finite(dimension)))
}

/// Return QPA's installed `ProjDimension(M)` attribute, if it is bound.
///
/// QPA declares this attribute without a computing method: bounded
/// [`proj_dimension_of_module`] calls install finite values when they find
/// them. Rust exposes the unbound state explicitly as `None`.
#[must_use]
pub fn proj_dimension<F: Field>(module: &PathAlgebraModule<F>) -> Option<usize> {
    module.cached_projective_dimension()
}

/// Translate QPA's bounded `ProjDimensionOfModule(M, n)` method.
///
/// An already-installed dimension is compared with the entered GAP integer.
/// Otherwise QPA scans its minimal projective resolution in degrees
/// `1, ..., n + 1`, installs `i - 1` at the first zero object, and returns
/// `false` when the scan does not find one. `None` represents that `false`.
/// For an uncached negative bound, QPA still constructs the resolution and
/// then returns `false` without evaluating a positive-degree object; that
/// eager constructor ordering is retained.
///
/// # Errors
///
/// Propagates lazy projective-resolution construction and evaluation
/// failures.
pub fn proj_dimension_of_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<Option<usize>, EnvelopingModuleError> {
    if let Some(dimension) = module.cached_projective_dimension() {
        return Ok((dimension <= *bound).then_some(dimension));
    }
    let cutoff = bound + &Integer::from(2);
    let resolution = projective_resolution(module)?;
    let mut degree = Integer::one();
    while degree < cutoff {
        if resolution.object(&degree)?.dimension() == 0 {
            let dimension = usize::try_from(degree.clone()).map_err(|_| {
                EnvelopingModuleError::ResolutionDegreeNotRepresentable {
                    degree: degree.clone(),
                }
            })? - 1;
            return Ok(Some(module.cache_projective_dimension(dimension)));
        }
        degree += &Integer::one();
    }
    Ok(None)
}

/// Return QPA's installed `InjDimension(M)` attribute, if it is bound.
///
/// As with [`proj_dimension`], QPA only binds this attribute as a side effect
/// of a successful bounded computation.
#[must_use]
pub fn inj_dimension<F: Field>(module: &PathAlgebraModule<F>) -> Option<usize> {
    module.cached_injective_dimension()
}

/// Translate QPA's bounded `InjDimensionOfModule(M, n)` method.
///
/// The source computes projective dimension on `DualOfModule(M)`. Its exact
/// cache ordering is retained: an already-installed injective dimension is
/// compared with `n`, while an already-installed projective dimension on the
/// dual is returned immediately even when it exceeds `n`. Otherwise the same
/// degree scan as [`proj_dimension_of_module`] installs only the injective
/// dimension on `M`. `None` represents QPA's `false`.
///
/// # Errors
///
/// Propagates dualization and lazy projective-resolution failures. As in the
/// projective method, an uncached negative bound still constructs the dual
/// resolution before its scan ends immediately.
pub fn inj_dimension_of_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<Option<usize>, EnvelopingModuleError> {
    if let Some(dimension) = module.cached_injective_dimension() {
        return Ok((dimension <= *bound).then_some(dimension));
    }
    let dual = dual_of_module(module)?;
    if let Some(dimension) = dual.cached_projective_dimension() {
        return Ok(Some(dimension));
    }
    let cutoff = bound + &Integer::from(2);
    let resolution = projective_resolution(&dual)?;
    let mut degree = Integer::one();
    while degree < cutoff {
        if resolution.object(&degree)?.dimension() == 0 {
            let dimension = usize::try_from(degree.clone()).map_err(|_| {
                EnvelopingModuleError::ResolutionDegreeNotRepresentable {
                    degree: degree.clone(),
                }
            })? - 1;
            return Ok(Some(module.cache_injective_dimension(dimension)));
        }
        degree += &Integer::one();
    }
    Ok(None)
}

/// Translate QPA's bounded `DominantDimensionOfModule(M, n)` method.
///
/// `None` preserves the method's `TryNextMethod()` on an infinite-dimensional
/// acting algebra. A projective-injective module returns `Infinite` before
/// the bound is used. Otherwise QPA dualizes `M` and scans the projective
/// resolution of the dual for the first noninjective projective term.
/// `ExceedsBound` is QPA's `false`. In particular, negative bounds retain the
/// source behavior: degree `n` is still inspected, so the result can be
/// finite zero or `false`.
///
/// # Errors
///
/// Propagates finite-dimensionality, projectivity/injectivity, dualization,
/// and lazy projective-resolution failures.
pub fn dominant_dimension_of_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<Option<DominantDimensionResult>, EnvelopingModuleError> {
    if finite_canonical_basis(module.right_acting_algebra())?.is_none() {
        return Ok(None);
    }
    if is_injective_module(module)? && is_projective_module(module)? {
        return Ok(Some(DominantDimensionResult::Infinite));
    }

    let dual = dual_of_module(module)?;
    let resolution = projective_resolution(&dual)?;
    let mut degree = Integer::from(0);
    while degree < *bound && is_injective_module(&resolution.object(&degree)?)? {
        degree += &Integer::one();
    }
    if degree < *bound {
        let dimension = usize::try_from(degree.clone()).map_err(|_| {
            EnvelopingModuleError::ResolutionDegreeNotRepresentable {
                degree: degree.clone(),
            }
        })?;
        return Ok(Some(DominantDimensionResult::Finite(dimension)));
    }
    if !is_injective_module(&resolution.object(bound)?)? {
        let dimension = if bound < &Integer::from(0) {
            0
        } else {
            usize::try_from(bound.clone()).map_err(|_| {
                EnvelopingModuleError::ResolutionDegreeNotRepresentable {
                    degree: bound.clone(),
                }
            })?
        };
        return Ok(Some(DominantDimensionResult::Finite(dimension)));
    }
    Ok(Some(DominantDimensionResult::ExceedsBound))
}

/// Translate QPA's bounded `DominantDimensionOfAlgebra(A, n)` method.
///
/// `None` preserves the method's explicit `TryNextMethod()` for
/// infinite-dimensional input. A selfinjective algebra returns `Infinite`.
/// Otherwise QPA computes the module result for every indecomposable
/// projective, discards `Infinite` and `ExceedsBound`, and returns the minimum
/// remaining finite value; if there is none, it returns `false`, represented
/// by `ExceedsBound`.
///
/// # Errors
///
/// Propagates finite-dimensionality, selfinjectivity, indecomposable
/// projective construction, and module dominant-dimension failures.
pub fn dominant_dimension_of_algebra<F: QpaCoefficientField>(
    algebra: &QuiverAlgebra<F>,
    bound: &Integer,
) -> Result<Option<DominantDimensionResult>, EnvelopingModuleError> {
    if finite_canonical_basis(algebra)?.is_none() {
        return Ok(None);
    }
    match is_selfinjective_algebra(algebra)? {
        Some(true) => return Ok(Some(DominantDimensionResult::Infinite)),
        Some(false) => {}
        None => {
            return Err(EnvelopingModuleError::DominantSelfinjectivityUnavailable);
        }
    }
    let mut minimum = None;
    for projective in indec_projective_modules(algebra)? {
        let result = dominant_dimension_of_module(&projective, bound)?
            .ok_or(EnvelopingModuleError::FiniteDominantModuleUnavailable)?;
        if let DominantDimensionResult::Finite(dimension) = result {
            minimum = Some(minimum.map_or(dimension, |current: usize| current.min(dimension)));
        }
    }
    Ok(Some(minimum.map_or(
        DominantDimensionResult::ExceedsBound,
        DominantDimensionResult::Finite,
    )))
}

/// Return QPA's installed `GorensteinDimension(A)` attribute, if bound.
///
/// The bounded operation installs finite values; Rust exposes QPA's unbound
/// attribute state as `None`.
#[must_use]
pub fn gorenstein_dimension<F: Field>(algebra: &QuiverAlgebra<F>) -> Option<usize> {
    match algebra {
        QuiverAlgebra::Path(algebra) => algebra.cached_gorenstein_dimension(),
        QuiverAlgebra::Quotient(algebra) => algebra.cached_gorenstein_dimension(),
    }
}

fn cache_gorenstein_dimension<F: Field>(algebra: &QuiverAlgebra<F>, value: usize) -> usize {
    match algebra {
        QuiverAlgebra::Path(algebra) => algebra.cache_gorenstein_dimension(value),
        QuiverAlgebra::Quotient(algebra) => algebra.cache_gorenstein_dimension(value),
    }
}

/// Translate QPA's bounded `GorensteinDimensionOfAlgebra(A, n)` method.
///
/// `None` is QPA's `false`. QPA first returns any known finite global
/// dimension, regardless of `n`, and then any cached Gorenstein dimension,
/// again without comparing the bound. Otherwise it takes the maximum
/// projective dimension of the indecomposable injectives on `A` and on
/// `OppositeAlgebra(A)`, stopping with `false` as soon as either side exceeds
/// `n`. The ordinary Rust entry point suppresses QPA's conjecture diagnostic;
/// use [`gorenstein_dimension_of_algebra_with_progress`] to receive it.
///
/// # Errors
///
/// Propagates indecomposable-injective, bounded projective-dimension, and
/// opposite-algebra construction failures.
pub fn gorenstein_dimension_of_algebra<F: QpaCoefficientField>(
    algebra: &QuiverAlgebra<F>,
    bound: &Integer,
) -> Result<Option<usize>, EnvelopingModuleError> {
    gorenstein_dimension_of_algebra_with_progress(algebra, bound, |_| {})
}

/// Evaluate [`gorenstein_dimension_of_algebra`] while exposing QPA's printed
/// warning when the left and right dimensions differ.
///
/// # Errors
///
/// Propagates the same failures as [`gorenstein_dimension_of_algebra`].
pub fn gorenstein_dimension_of_algebra_with_progress<F, R>(
    algebra: &QuiverAlgebra<F>,
    bound: &Integer,
    mut report: R,
) -> Result<Option<usize>, EnvelopingModuleError>
where
    F: QpaCoefficientField,
    R: FnMut(&'static str),
{
    let known_global_dimension = match algebra {
        QuiverAlgebra::Path(algebra) => Some(AlgebraGlobalDimension::Finite(
            algebra.qpa_global_dimension(),
        )),
        QuiverAlgebra::Quotient(algebra) => algebra.cached_global_dimension(),
    };
    if let Some(AlgebraGlobalDimension::Finite(dimension)) = known_global_dimension {
        return Ok(Some(cache_gorenstein_dimension(algebra, dimension)));
    }
    if let Some(dimension) = gorenstein_dimension(algebra) {
        return Ok(Some(dimension));
    }

    let mut dimension_right = 0;
    for injective in indec_injective_modules(algebra)? {
        let Some(dimension) = proj_dimension_of_module(&injective, bound)? else {
            return Ok(None);
        };
        dimension_right = dimension_right.max(dimension);
    }

    let opposite = opposite_module_algebra(algebra)?;
    let mut dimension_left = 0;
    for injective in indec_injective_modules(&opposite)? {
        let Some(dimension) = proj_dimension_of_module(&injective, bound)? else {
            return Ok(None);
        };
        dimension_left = dimension_left.max(dimension);
    }
    if dimension_left != dimension_right {
        report("You have a counterexample to the Gorenstein conjecture!\n\n");
    }
    Ok(Some(cache_gorenstein_dimension(
        algebra,
        dimension_left.max(dimension_right),
    )))
}

/// Translate QPA's Nakayama-socle test for `IsSelfinjectiveAlgebra(A)`.
///
/// Infinite-dimensional input returns `Some(false)` before path-algebra
/// dispatch. A finite path algebra is selfinjective exactly when its quiver
/// has no arrows. For a finite admissible quotient, QPA requires every
/// indecomposable projective socle to be one-dimensional and every simple to
/// occur exactly once across their dimension vectors. `None` preserves
/// `TryNextMethod()` for a finite nonadmissible quotient.
///
/// # Errors
///
/// Propagates default finite-basis and admissibility computations, projective
/// and socle construction, and reports an unrepresentable dimension-vector
/// sum.
pub fn is_selfinjective_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<bool>, EnvelopingModuleError> {
    if let QuiverAlgebra::Quotient(quotient) = algebra
        && let Some(cached) = quotient.cached_is_selfinjective_algebra()
    {
        return Ok(Some(cached));
    }
    if finite_canonical_basis(algebra)?.is_none() {
        if let QuiverAlgebra::Quotient(quotient) = algebra {
            return Ok(Some(quotient.cache_is_selfinjective_algebra(false)));
        }
        return Ok(Some(false));
    }
    if let QuiverAlgebra::Path(path_algebra) = algebra {
        return Ok(Some(path_algebra.quiver().number_of_arrows() == 0));
    }
    let QuiverAlgebra::Quotient(quotient) = algebra else {
        unreachable!("every quiver algebra is a path algebra or quotient")
    };
    if !quotient_is_admissible(quotient)? {
        return Ok(None);
    }

    let socles = indec_projective_modules(algebra)?
        .into_iter()
        .map(socle_of_module)
        .collect::<Result<Vec<_>, _>>()?;
    if socles.iter().any(|socle| socle.dimension() != 1) {
        return Ok(Some(quotient.cache_is_selfinjective_algebra(false)));
    }
    let mut multiplicities = vec![0usize; algebra.quiver().number_of_vertices()];
    for socle in socles {
        for (vertex, (&dimension, total)) in socle
            .dimension_vector()
            .iter()
            .zip(&mut multiplicities)
            .enumerate()
        {
            *total = total
                .checked_add(dimension)
                .ok_or(EnvelopingModuleError::SelfinjectiveDimensionVectorOverflow { vertex })?;
        }
    }
    let result = multiplicities.into_iter().all(|entry| entry == 1);
    let result = quotient.cache_is_selfinjective_algebra(result);
    if result && quotient.cached_is_semisimple_algebra() == Some(false) {
        let _ = quotient.cache_global_dimension(AlgebraGlobalDimension::Infinite);
    }
    Ok(Some(result))
}

/// Translate QPA's `CartanMatrix(A)` method for finite bound quiver algebras.
///
/// Rows are the dimension vectors of `IndecProjectiveModules(A)` in quiver
/// vertex order. The matrix is integral and independent of the coefficient
/// characteristic. `Ok(None)` is QPA's cached `fail` for an infinite algebra
/// or for a finite quotient that is not admissible.
///
/// # Errors
///
/// Propagates QPA's finite-dimensionality/admissibility backend and exact
/// indecomposable-projective construction.
pub fn cartan_matrix<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<Vec<Vec<usize>>>, EnvelopingModuleError> {
    match algebra.cached_cartan_matrix_attribute() {
        CachedAttribute::Value(value) => return Ok(Some(value)),
        CachedAttribute::Fail | CachedAttribute::False => return Ok(None),
        CachedAttribute::Unbound => {}
    }
    if !is_finite_dimensional(algebra)? {
        return Ok(algebra.cache_cartan_matrix_attribute(None));
    }
    if let QuiverAlgebra::Quotient(quotient) = algebra
        && !quotient_is_admissible(quotient)?
    {
        return Ok(algebra.cache_cartan_matrix_attribute(None));
    }
    let matrix = indec_projective_modules(algebra)?
        .into_iter()
        .map(|projective| projective.dimension_vector().to_vec())
        .collect();
    Ok(algebra.cache_cartan_matrix_attribute(Some(matrix)))
}

/// Translate QPA's `CoxeterMatrix(A) = -C^-1 C^T` method.
///
/// GAP constructs `C` as an integer matrix, so this result is always an exact
/// rational matrix, including when `A` is defined over a prime field. This is
/// deliberately not computed in the coefficient field. `Ok(None)` preserves
/// cached `fail` from `CartanMatrix` or a singular Cartan matrix.
///
/// # Errors
///
/// Propagates Cartan construction and checked exact rational matrix
/// construction, inversion, and multiplication.
pub fn coxeter_matrix<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<DenseMatrix<RationalField>>, EnvelopingModuleError> {
    match algebra.cached_coxeter_matrix_attribute() {
        CachedAttribute::Value(value) => return Ok(Some(value)),
        CachedAttribute::Fail | CachedAttribute::False => return Ok(None),
        CachedAttribute::Unbound => {}
    }
    let Some(cartan_rows) = cartan_matrix(algebra)? else {
        return Ok(algebra.cache_coxeter_matrix_attribute(None));
    };
    let rational_rows = cartan_rows
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|entry| Q.nth(Integer::from(entry)))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let cartan = DenseMatrix::from_rows(rational_rows, Q)?;
    if Q.is_zero(&cartan.determinant()?) {
        return Ok(algebra.cache_coxeter_matrix_attribute(None));
    }
    let inverse = cartan
        .inverse()?
        .ok_or(EnvelopingModuleError::NonzeroCartanDeterminantNotInvertible)?;
    let negative_one = Q.neg(&Q.one());
    let matrix = inverse
        .multiply(&cartan.transpose())?
        .scale_left(&negative_one);
    Ok(algebra.cache_coxeter_matrix_attribute(Some(matrix)))
}

/// Translate QPA's cached `CoxeterPolynomial(A)` attribute.
///
/// The returned exact rational coefficients are in ascending degree order,
/// matching `CoefficientsOfUnivariatePolynomial` in the oracle fixtures.
/// `Ok(None)` preserves a cached failed Coxeter-matrix computation.
///
/// # Errors
///
/// Propagates Coxeter-matrix construction and GAP's source-translated exact
/// characteristic-polynomial algorithm.
pub fn coxeter_polynomial<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<Vec<Rational>>, EnvelopingModuleError> {
    match algebra.cached_coxeter_polynomial_attribute() {
        CachedAttribute::Value(value) => return Ok(Some(value)),
        CachedAttribute::Fail | CachedAttribute::False => return Ok(None),
        CachedAttribute::Unbound => {}
    }
    let Some(matrix) = coxeter_matrix(algebra)? else {
        return Ok(algebra.cache_coxeter_polynomial_attribute(None));
    };
    let coefficients = matrix.characteristic_polynomial_coefficients()?;
    Ok(algebra.cache_coxeter_polynomial_attribute(Some(coefficients)))
}

/// The two closures returned by QPA's `NakayamaPermutation(A)` attribute.
///
/// `image_of_module` is the first closure on simple modules and
/// `image_of_index` is the second closure on the one-based indexing set. The
/// quotient case retains the exact socle modules captured by the source
/// method; the selfinjective path-algebra case is the identity closure.
#[derive(Clone, Debug)]
pub struct NakayamaPermutation<F: Field> {
    algebra: QuiverAlgebra<F>,
    index_images: Vec<usize>,
    quotient_socles: Option<Vec<PathAlgebraModule<F>>>,
}

impl<F: Field> NakayamaPermutation<F> {
    /// Return the algebra whose attribute produced these closures.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Return the one-based image of every simple-module index.
    #[must_use]
    pub fn index_images(&self) -> &[usize] {
        &self.index_images
    }

    /// Apply QPA's simple-module closure.
    ///
    /// The source checks only that the entered module has total dimension one;
    /// it does not require the module to act over the captured algebra. A path
    /// algebra returns the entered module itself. A quotient selects the
    /// captured projective socle at the position of the dimension-vector `1`.
    ///
    /// # Errors
    ///
    /// Preserves QPA's error for a nonsimple module and reports the source's
    /// eventual list-index failure if a foreign simple has too many vertices.
    pub fn image_of_module(
        &self,
        module: &PathAlgebraModule<F>,
    ) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
        if module.dimension() != 1 {
            return Err(EnvelopingModuleError::NakayamaPermutationModuleNotSimple);
        }
        let Some(socles) = &self.quotient_socles else {
            return Ok(module.clone());
        };
        let position = module
            .dimension_vector()
            .iter()
            .position(|&dimension| dimension == 1)
            .ok_or(EnvelopingModuleError::NakayamaPermutationModuleNotSimple)?;
        socles.get(position).cloned().ok_or(
            EnvelopingModuleError::NakayamaPermutationModulePositionOutOfRange {
                position: position + 1,
                vertex_count: socles.len(),
            },
        )
    }

    /// Apply QPA's one-based index closure.
    ///
    /// # Errors
    ///
    /// Preserves the explicit source error for every integer outside the
    /// algebra's simple-module indexing set.
    pub fn image_of_index(&self, index: &Integer) -> Result<usize, EnvelopingModuleError> {
        let position = (0..self.index_images.len())
            .find(|&position| index == &Integer::from(position + 1))
            .ok_or_else(
                || EnvelopingModuleError::IncorrectNakayamaPermutationIndex {
                    index: index.clone(),
                    vertex_count: self.index_images.len(),
                },
            )?;
        Ok(self.index_images[position])
    }
}

/// Translate QPA's cached `NakayamaPermutation(A)` attribute.
///
/// `Ok(None)` is the attribute's cached boolean `false` for a nonselfinjective
/// algebra. An admissible selfinjective quotient captures the ordered socles
/// of its indecomposable projectives and their simple indices. A selfinjective
/// path algebra has no arrows, so both closures are identities.
///
/// # Errors
///
/// Preserves the absence of a selfinjectivity method on a nonadmissible
/// quotient and propagates projective/socle construction.
pub fn nakayama_permutation<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<NakayamaPermutation<F>>, EnvelopingModuleError> {
    match algebra.cached_nakayama_permutation_attribute() {
        CachedAttribute::Value(index_images) => {
            return Ok(Some(rebuild_nakayama_permutation(algebra, index_images)?));
        }
        CachedAttribute::False => return Ok(None),
        CachedAttribute::Unbound => {}
        CachedAttribute::Fail => {
            return Err(EnvelopingModuleError::InvalidNakayamaPermutationCacheState);
        }
    }

    match is_selfinjective_algebra(algebra)? {
        Some(false) => {
            let _ = algebra.cache_nakayama_permutation_attribute(None);
            return Ok(None);
        }
        None => return Err(EnvelopingModuleError::NakayamaPermutationSelfinjectivityUnavailable),
        Some(true) => {}
    }

    let (index_images, quotient_socles) = match algebra {
        QuiverAlgebra::Path(_) => ((1..=algebra.quiver().number_of_vertices()).collect(), None),
        QuiverAlgebra::Quotient(_) => {
            let socles = indec_projective_modules(algebra)?
                .into_iter()
                .map(socle_of_module)
                .collect::<Result<Vec<_>, _>>()?;
            let indices = nakayama_indices_from_socles(&socles)?;
            (indices, Some(socles))
        }
    };
    let index_images = algebra
        .cache_nakayama_permutation_attribute(Some(index_images))
        .ok_or(EnvelopingModuleError::NakayamaPermutationCacheConflict)?;
    Ok(Some(NakayamaPermutation {
        algebra: algebra.clone(),
        index_images,
        quotient_socles,
    }))
}

fn rebuild_nakayama_permutation<F: Field>(
    algebra: &QuiverAlgebra<F>,
    index_images: Vec<usize>,
) -> Result<NakayamaPermutation<F>, EnvelopingModuleError> {
    let quotient_socles = if matches!(algebra, QuiverAlgebra::Quotient(_)) {
        Some(
            indec_projective_modules(algebra)?
                .into_iter()
                .map(socle_of_module)
                .collect::<Result<Vec<_>, _>>()?,
        )
    } else {
        None
    };
    Ok(NakayamaPermutation {
        algebra: algebra.clone(),
        index_images,
        quotient_socles,
    })
}

fn nakayama_indices_from_socles<F: Field>(
    socles: &[PathAlgebraModule<F>],
) -> Result<Vec<usize>, EnvelopingModuleError> {
    socles
        .iter()
        .enumerate()
        .map(|(projective, socle)| {
            socle
                .dimension_vector()
                .iter()
                .position(|&dimension| dimension == 1)
                .map(|position| position + 1)
                .ok_or(
                    EnvelopingModuleError::SelfinjectiveProjectiveSocleHasNoSimplePosition {
                        projective: projective + 1,
                    },
                )
        })
        .collect()
}

/// QPA's captured Nakayama-automorphism closure.
///
/// The matrix rows are the canonical-basis coordinates of the images of the
/// canonical basis. QPA computes it as `P * (TransposedMat(P)^-1)` from its
/// Frobenius pairing matrix `P`.
#[derive(Clone, Debug)]
pub struct NakayamaAutomorphism<F: Field> {
    algebra: QuiverAlgebra<F>,
    canonical_basis: Vec<QuiverAlgebraElement<F>>,
    matrix: DenseMatrix<F>,
}

impl<F: Field> NakayamaAutomorphism<F> {
    /// Return the algebra captured by QPA's closure.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Return the canonical-basis matrix of the automorphism.
    #[must_use]
    pub const fn matrix(&self) -> &DenseMatrix<F> {
        &self.matrix
    }

    /// Apply QPA's Nakayama-automorphism closure.
    ///
    /// # Errors
    ///
    /// Preserves the source closure's explicit rejection of an element that
    /// does not belong to the captured algebra and propagates exact canonical
    /// coordinate arithmetic.
    pub fn image(
        &self,
        element: &QuiverAlgebraElement<F>,
    ) -> Result<QuiverAlgebraElement<F>, EnvelopingModuleError> {
        if !self.algebra.contains_element(element) {
            return Err(EnvelopingModuleError::NakayamaAutomorphismElementOutsideAlgebra);
        }
        let coordinates =
            canonical_basis_coefficients(&self.algebra, &self.canonical_basis, element)?;
        let image_coordinates = multiply_row_by_dense_matrix(&coordinates, &self.matrix)?;
        linear_combination_of_algebra_basis(
            &self.algebra,
            &self.canonical_basis,
            &image_coordinates,
        )
    }
}

/// QPA's captured Frobenius bilinear-form closure.
#[derive(Clone, Debug)]
pub struct FrobeniusForm<F: Field> {
    quotient: QuotientPathAlgebra<F>,
    canonical_basis: Vec<QuiverAlgebraElement<F>>,
    matrix: DenseMatrix<F>,
}

impl<F: Field> FrobeniusForm<F> {
    /// Return the quotient algebra captured by the closure.
    #[must_use]
    pub const fn algebra(&self) -> &QuotientPathAlgebra<F> {
        &self.quotient
    }

    /// Return QPA's multiplication-pairing matrix `P`.
    #[must_use]
    pub const fn matrix(&self) -> &DenseMatrix<F> {
        &self.matrix
    }

    /// Evaluate `pi(x * y)` through QPA's `x * P * y` closure.
    ///
    /// # Errors
    ///
    /// Both arguments must belong to the identical quotient algebra captured
    /// by the source closure.
    pub fn evaluate(
        &self,
        left: &QuiverAlgebraElement<F>,
        right: &QuiverAlgebraElement<F>,
    ) -> Result<F::Element, EnvelopingModuleError> {
        let algebra = QuiverAlgebra::Quotient(self.quotient.clone());
        if !algebra.contains_element(left) || !algebra.contains_element(right) {
            return Err(EnvelopingModuleError::FrobeniusFormElementOutsideAlgebra);
        }
        let left = canonical_basis_coefficients(&algebra, &self.canonical_basis, left)?;
        let right = canonical_basis_coefficients(&algebra, &self.canonical_basis, right)?;
        let paired = multiply_row_by_dense_matrix(&left, &self.matrix)?;
        Ok(dot_product(algebra.field(), &paired, &right))
    }
}

/// QPA's captured Frobenius linear-functional closure.
#[derive(Clone, Debug)]
pub struct FrobeniusLinearFunctional<F: Field> {
    quotient: QuotientPathAlgebra<F>,
    canonical_basis: Vec<QuiverAlgebraElement<F>>,
    coordinates: Vec<F::Element>,
}

impl<F: Field> FrobeniusLinearFunctional<F> {
    /// Return the quotient algebra captured by the closure.
    #[must_use]
    pub const fn algebra(&self) -> &QuotientPathAlgebra<F> {
        &self.quotient
    }

    /// Return the functional values on QPA's canonical basis.
    #[must_use]
    pub fn canonical_basis_values(&self) -> &[F::Element] {
        &self.coordinates
    }

    /// Evaluate QPA's `pi_map` closure.
    ///
    /// # Errors
    ///
    /// The argument must belong to the identical quotient algebra captured by
    /// the source closure.
    pub fn evaluate(
        &self,
        element: &QuiverAlgebraElement<F>,
    ) -> Result<F::Element, EnvelopingModuleError> {
        let algebra = QuiverAlgebra::Quotient(self.quotient.clone());
        if !algebra.contains_element(element) {
            return Err(EnvelopingModuleError::FrobeniusFunctionalElementOutsideAlgebra);
        }
        let element = canonical_basis_coefficients(&algebra, &self.canonical_basis, element)?;
        Ok(dot_product(algebra.field(), &element, &self.coordinates))
    }
}

/// The finite integer or `infinity` returned by GAP's matrix `Order` method.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum NakayamaAutomorphismOrder {
    /// The automorphism has this finite positive order.
    Finite(Integer),
    /// The automorphism has infinite order.
    Infinite,
}

/// Translate QPA's cached `NakayamaAutomorphism(A)` attribute.
///
/// A nonselfinjective algebra returns and caches boolean `false`. For a
/// quotient, the implementation follows QPA's projective-socle construction,
/// GAP's `BaseSteinitzVectors` ordering, the resulting Frobenius functional,
/// and `P * (P^T)^-1`. Successful evaluation also binds `FrobeniusForm` and
/// `FrobeniusLinearFunctional`, exactly as the source method does. A finite
/// selfinjective path algebra has no arrows and returns the identity closure.
///
/// # Errors
///
/// Propagates selfinjectivity dispatch, module/socle construction, canonical
/// coordinates, quotient multiplication, and exact matrix arithmetic.
pub fn nakayama_automorphism<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<NakayamaAutomorphism<F>>, EnvelopingModuleError> {
    match algebra.cached_nakayama_automorphism_attribute() {
        CachedAttribute::Value(matrix) => {
            return Ok(Some(rebuild_nakayama_automorphism(algebra, matrix)?));
        }
        CachedAttribute::False => return Ok(None),
        CachedAttribute::Unbound => {}
        CachedAttribute::Fail => {
            return Err(EnvelopingModuleError::InvalidNakayamaAutomorphismCacheState);
        }
    }

    match is_selfinjective_algebra(algebra)? {
        Some(false) => {
            let _ = algebra.cache_nakayama_automorphism_attribute(None);
            return Ok(None);
        }
        None => return Err(EnvelopingModuleError::NakayamaAutomorphismSelfinjectivityUnavailable),
        Some(true) => {}
    }

    match algebra {
        QuiverAlgebra::Path(_) => {
            let canonical_basis = finite_canonical_basis(algebra)?
                .ok_or(EnvelopingModuleError::MissingNakayamaCanonicalBasis)?;
            let matrix = DenseMatrix::identity(canonical_basis.len(), algebra.field().clone());
            let matrix = algebra
                .cache_nakayama_automorphism_attribute(Some(matrix))
                .ok_or(EnvelopingModuleError::NakayamaAutomorphismCacheConflict)?;
            Ok(Some(NakayamaAutomorphism {
                algebra: algebra.clone(),
                canonical_basis,
                matrix,
            }))
        }
        QuiverAlgebra::Quotient(quotient) => {
            construct_quotient_nakayama_automorphism(algebra, quotient)
        }
    }
}

/// Translate QPA's cached `FrobeniusForm(A)` quotient attribute.
///
/// Calling this first forces `NakayamaAutomorphism`, whose successful source
/// method installs the form. `Ok(None)` is the attribute's cached boolean
/// `false` for a nonselfinjective quotient.
///
/// # Errors
///
/// Propagates Nakayama-automorphism construction and cache consistency.
pub fn frobenius_form<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<Option<FrobeniusForm<F>>, EnvelopingModuleError> {
    match quotient.cached_frobenius_form_attribute() {
        CachedAttribute::Value(matrix) => {
            return rebuild_frobenius_form(quotient, matrix).map(Some);
        }
        CachedAttribute::False => return Ok(None),
        CachedAttribute::Unbound => {}
        CachedAttribute::Fail => return Err(EnvelopingModuleError::InvalidFrobeniusFormCacheState),
    }
    let algebra = QuiverAlgebra::Quotient(quotient.clone());
    if nakayama_automorphism(&algebra)?.is_none() {
        let _ = quotient.cache_frobenius_form_attribute(None);
        return Ok(None);
    }
    match quotient.cached_frobenius_form_attribute() {
        CachedAttribute::Value(matrix) => rebuild_frobenius_form(quotient, matrix).map(Some),
        CachedAttribute::False => Ok(None),
        CachedAttribute::Unbound | CachedAttribute::Fail => {
            Err(EnvelopingModuleError::MissingInstalledFrobeniusForm)
        }
    }
}

/// Translate QPA's cached `FrobeniusLinearFunctional(A)` quotient attribute.
///
/// Calling this first forces `NakayamaAutomorphism`, whose successful source
/// method installs the functional. `Ok(None)` is cached boolean `false` for a
/// nonselfinjective quotient.
///
/// # Errors
///
/// Propagates Nakayama-automorphism construction and cache consistency.
pub fn frobenius_linear_functional<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
) -> Result<Option<FrobeniusLinearFunctional<F>>, EnvelopingModuleError> {
    match quotient.cached_frobenius_linear_functional_attribute() {
        CachedAttribute::Value(coordinates) => {
            return rebuild_frobenius_functional(quotient, coordinates).map(Some);
        }
        CachedAttribute::False => return Ok(None),
        CachedAttribute::Unbound => {}
        CachedAttribute::Fail => {
            return Err(EnvelopingModuleError::InvalidFrobeniusFunctionalCacheState);
        }
    }
    let algebra = QuiverAlgebra::Quotient(quotient.clone());
    if nakayama_automorphism(&algebra)?.is_none() {
        let _ = quotient.cache_frobenius_linear_functional_attribute(None);
        return Ok(None);
    }
    match quotient.cached_frobenius_linear_functional_attribute() {
        CachedAttribute::Value(coordinates) => {
            rebuild_frobenius_functional(quotient, coordinates).map(Some)
        }
        CachedAttribute::False => Ok(None),
        CachedAttribute::Unbound | CachedAttribute::Fail => {
            Err(EnvelopingModuleError::MissingInstalledFrobeniusFunctional)
        }
    }
}

/// Translate QPA's cached `OrderOfNakayamaAutomorphism(A)` attribute.
///
/// The quotient method first constructs the canonical-basis matrix through
/// `NakayamaAutomorphism` and delegates to GAP's matrix-order semantics. Over
/// a supported prime field the exact exponent of `GL(n,p)` bounds the order;
/// over the rationals the lcm of all cyclotomic indices with totient at most
/// `n` is the exact finite-order test bound. Prime-factor removal then returns
/// the least order. A nonselfinjective path algebra returns boolean `false`;
/// QPA's quotient method instead attempts to call `false` and raises an error.
///
/// # Errors
///
/// Propagates Nakayama construction and preserves the quotient source's
/// attempted-call error on a nonselfinjective algebra.
pub fn order_of_nakayama_automorphism<F: QpaCoefficientField>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<NakayamaAutomorphismOrder>, EnvelopingModuleError> {
    match algebra.cached_order_of_nakayama_automorphism_attribute() {
        CachedAttribute::Value(value) => return Ok(Some(order_from_cache(value))),
        CachedAttribute::False => return Ok(None),
        CachedAttribute::Unbound => {}
        CachedAttribute::Fail => {
            return Err(EnvelopingModuleError::InvalidNakayamaOrderCacheState);
        }
    }

    if matches!(algebra, QuiverAlgebra::Path(_)) {
        match is_selfinjective_algebra(algebra)? {
            Some(false) => {
                let _ = algebra.cache_order_of_nakayama_automorphism_attribute(None);
                return Ok(None);
            }
            None => return Err(EnvelopingModuleError::NakayamaOrderSelfinjectivityUnavailable),
            Some(true) => {}
        }
        let value = CachedNakayamaAutomorphismOrder::Finite(Integer::one());
        let cached = algebra
            .cache_order_of_nakayama_automorphism_attribute(Some(value))
            .ok_or(EnvelopingModuleError::NakayamaOrderCacheConflict)?;
        return Ok(Some(order_from_cache(cached)));
    }

    let automorphism = nakayama_automorphism(algebra)?
        .ok_or(EnvelopingModuleError::CalledFalseNakayamaAutomorphism)?;
    let order = exact_matrix_order(automorphism.matrix(), algebra.field())?;
    let cached_value = match &order {
        NakayamaAutomorphismOrder::Finite(value) => {
            CachedNakayamaAutomorphismOrder::Finite(value.clone())
        }
        NakayamaAutomorphismOrder::Infinite => CachedNakayamaAutomorphismOrder::Infinite,
    };
    let cached = algebra
        .cache_order_of_nakayama_automorphism_attribute(Some(cached_value))
        .ok_or(EnvelopingModuleError::NakayamaOrderCacheConflict)?;
    Ok(Some(order_from_cache(cached)))
}

fn construct_quotient_nakayama_automorphism<F: Field>(
    algebra: &QuiverAlgebra<F>,
    quotient: &QuotientPathAlgebra<F>,
) -> Result<Option<NakayamaAutomorphism<F>>, EnvelopingModuleError> {
    let canonical_basis = finite_canonical_basis(algebra)?
        .ok_or(EnvelopingModuleError::MissingNakayamaCanonicalBasis)?;
    let projectives = indec_projective_modules(algebra)?;
    let mut socle_rows = Vec::with_capacity(projectives.len());
    for (position, projective) in projectives.into_iter().enumerate() {
        let inclusion = socle_of_module_inclusion(projective)?;
        let generator = minimal_generating_set_of_module(inclusion.source().clone())?
            .into_iter()
            .next()
            .ok_or(EnvelopingModuleError::EmptyProjectiveSocleGenerator {
                projective: position + 1,
            })?;
        let socle_element = inclusion.image(&generator)?;
        let parent_element = element_in_indec_projective(algebra, &socle_element, position + 1)?;
        let quotient_element = quotient.element(parent_element)?;
        socle_rows.push(canonical_basis_coefficients(
            algebra,
            &canonical_basis,
            &QuiverAlgebraElement::Quotient(quotient_element),
        )?);
    }

    let dimension = canonical_basis.len();
    let socle_matrix = DenseMatrix::from_rows(socle_rows, algebra.field().clone())?;
    let adapted_basis = base_steinitz_vectors(
        &DenseMatrix::identity(dimension, algebra.field().clone()),
        &socle_matrix,
    )?;
    let adapted_inverse = adapted_basis
        .inverse()?
        .ok_or(EnvelopingModuleError::IncompleteNakayamaSteinitzBasis)?;
    let socle_dimension = socle_matrix.nrows();
    let functional = (0..dimension)
        .map(|row| {
            (0..socle_dimension).fold(algebra.field().zero(), |mut sum, column| {
                algebra
                    .field()
                    .add_assign(&mut sum, &adapted_inverse[(row, column)]);
                sum
            })
        })
        .collect::<Vec<_>>();

    let mut pairing_entries = Vec::with_capacity(dimension.saturating_mul(dimension));
    for left in &canonical_basis {
        for right in &canonical_basis {
            let product = left.multiply(right)?;
            let coordinates = canonical_basis_coefficients(algebra, &canonical_basis, &product)?;
            pairing_entries.push(dot_product(algebra.field(), &coordinates, &functional));
        }
    }
    let pairing = DenseMatrix::from_flat(
        dimension,
        dimension,
        pairing_entries,
        algebra.field().clone(),
    )?;
    let transpose_inverse = pairing
        .transpose()
        .inverse()?
        .ok_or(EnvelopingModuleError::SingularNakayamaPairing)?;
    let matrix = pairing.multiply(&transpose_inverse)?;

    let pairing = quotient
        .cache_frobenius_form_attribute(Some(pairing))
        .ok_or(EnvelopingModuleError::FrobeniusFormCacheConflict)?;
    let functional = quotient
        .cache_frobenius_linear_functional_attribute(Some(functional))
        .ok_or(EnvelopingModuleError::FrobeniusFunctionalCacheConflict)?;
    let matrix = algebra
        .cache_nakayama_automorphism_attribute(Some(matrix))
        .ok_or(EnvelopingModuleError::NakayamaAutomorphismCacheConflict)?;
    debug_assert_eq!(pairing.nrows(), matrix.nrows());
    debug_assert_eq!(functional.len(), matrix.nrows());
    Ok(Some(NakayamaAutomorphism {
        algebra: algebra.clone(),
        canonical_basis,
        matrix,
    }))
}

fn rebuild_nakayama_automorphism<F: Field>(
    algebra: &QuiverAlgebra<F>,
    matrix: DenseMatrix<F>,
) -> Result<NakayamaAutomorphism<F>, EnvelopingModuleError> {
    let canonical_basis = finite_canonical_basis(algebra)?
        .ok_or(EnvelopingModuleError::MissingNakayamaCanonicalBasis)?;
    if matrix.nrows() != canonical_basis.len() || matrix.ncols() != canonical_basis.len() {
        return Err(
            EnvelopingModuleError::InvalidNakayamaAutomorphismMatrixShape {
                shape: (matrix.nrows(), matrix.ncols()),
                dimension: canonical_basis.len(),
            },
        );
    }
    Ok(NakayamaAutomorphism {
        algebra: algebra.clone(),
        canonical_basis,
        matrix,
    })
}

fn rebuild_frobenius_form<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
    matrix: DenseMatrix<F>,
) -> Result<FrobeniusForm<F>, EnvelopingModuleError> {
    let algebra = QuiverAlgebra::Quotient(quotient.clone());
    let canonical_basis = finite_canonical_basis(&algebra)?
        .ok_or(EnvelopingModuleError::MissingNakayamaCanonicalBasis)?;
    if matrix.nrows() != canonical_basis.len() || matrix.ncols() != canonical_basis.len() {
        return Err(EnvelopingModuleError::InvalidFrobeniusFormMatrixShape {
            shape: (matrix.nrows(), matrix.ncols()),
            dimension: canonical_basis.len(),
        });
    }
    Ok(FrobeniusForm {
        quotient: quotient.clone(),
        canonical_basis,
        matrix,
    })
}

fn rebuild_frobenius_functional<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
    coordinates: Vec<F::Element>,
) -> Result<FrobeniusLinearFunctional<F>, EnvelopingModuleError> {
    let algebra = QuiverAlgebra::Quotient(quotient.clone());
    let canonical_basis = finite_canonical_basis(&algebra)?
        .ok_or(EnvelopingModuleError::MissingNakayamaCanonicalBasis)?;
    if coordinates.len() != canonical_basis.len() {
        return Err(EnvelopingModuleError::InvalidFrobeniusFunctionalLength {
            actual: coordinates.len(),
            dimension: canonical_basis.len(),
        });
    }
    Ok(FrobeniusLinearFunctional {
        quotient: quotient.clone(),
        canonical_basis,
        coordinates,
    })
}

fn base_steinitz_vectors<F: Field>(
    basis: &DenseMatrix<F>,
    subspace_generators: &DenseMatrix<F>,
) -> Result<DenseMatrix<F>, EnvelopingModuleError> {
    if basis.nrows() == 0 {
        return DenseMatrix::from_flat(0, 0, Vec::new(), basis.field().clone()).map_err(Into::into);
    }
    if basis.ncols() != subspace_generators.ncols() || basis.field() != subspace_generators.field()
    {
        return Err(EnvelopingModuleError::InvalidNakayamaSteinitzInput);
    }
    let subspace = subspace_generators.canonical_row_space_basis();
    let missing = basis
        .nrows()
        .checked_sub(subspace.nrows())
        .ok_or(EnvelopingModuleError::InvalidNakayamaSteinitzInput)?;
    let field = basis.field();
    let dimension = basis.ncols();
    let mut remaining = (0..basis.nrows())
        .map(|row| basis.row(row).to_vec())
        .collect::<Vec<_>>();
    let mut complement = Vec::with_capacity(missing);
    let mut subspace_row = 0usize;
    let mut column = 0usize;

    while complement.len() < missing {
        let mut stop = false;
        loop {
            let mut pivot = None;
            if column < dimension
                && (subspace_row >= subspace.nrows()
                    || field.is_zero(&subspace[(subspace_row, column)]))
            {
                if let Some(position) = remaining
                    .iter()
                    .position(|row| !field.is_zero(&row[column]))
                {
                    let mut vector = remaining[position].clone();
                    let inverse = field.inv(&vector[column]);
                    for entry in &mut vector {
                        field.mul_assign(entry, &inverse);
                    }
                    complement.push(vector.clone());
                    pivot = Some(vector);
                }
            } else {
                stop = true;
                if subspace_row < subspace.nrows() {
                    let mut vector = subspace.row(subspace_row).to_vec();
                    let inverse = field.inv(&vector[column]);
                    for entry in &mut vector {
                        field.mul_assign(entry, &inverse);
                    }
                    pivot = Some(vector);
                }
            }

            if let Some(vector) = pivot {
                for row in &mut remaining {
                    if !field.is_zero(&row[column]) {
                        let scale = field.div(&row[column], &vector[column]);
                        for (entry, pivot_entry) in row.iter_mut().zip(&vector) {
                            let product = field.mul(&scale, pivot_entry);
                            field.sub_assign(entry, &product);
                        }
                    }
                }
                remaining.retain(|row| row.iter().any(|entry| !field.is_zero(entry)));
            }
            column = column
                .checked_add(1)
                .ok_or(EnvelopingModuleError::NakayamaSteinitzColumnOverflow)?;
            if stop {
                break;
            }
            if column > dimension {
                return Err(EnvelopingModuleError::IncompleteNakayamaSteinitzBasis);
            }
        }
        subspace_row = subspace_row
            .checked_add(1)
            .ok_or(EnvelopingModuleError::NakayamaSteinitzRowOverflow)?;
    }

    let mut rows = (0..subspace.nrows())
        .map(|row| subspace.row(row).to_vec())
        .collect::<Vec<_>>();
    rows.extend(complement);
    DenseMatrix::from_rows(rows, field.clone()).map_err(Into::into)
}

fn multiply_row_by_dense_matrix<F: Field>(
    row: &[F::Element],
    matrix: &DenseMatrix<F>,
) -> Result<Vec<F::Element>, EnvelopingModuleError> {
    let row = DenseMatrix::from_flat(1, row.len(), row.to_vec(), matrix.field().clone())?;
    Ok(row.multiply(matrix)?.row(0).to_vec())
}

fn dot_product<F: Field>(field: &F, left: &[F::Element], right: &[F::Element]) -> F::Element {
    debug_assert_eq!(left.len(), right.len());
    left.iter()
        .zip(right)
        .fold(field.zero(), |mut sum, (left, right)| {
            field.add_mul_assign(&mut sum, left, right);
            sum
        })
}

fn linear_combination_of_algebra_basis<F: Field>(
    algebra: &QuiverAlgebra<F>,
    basis: &[QuiverAlgebraElement<F>],
    coefficients: &[F::Element],
) -> Result<QuiverAlgebraElement<F>, EnvelopingModuleError> {
    if basis.len() != coefficients.len() {
        return Err(EnvelopingModuleError::NakayamaLinearCombinationLength {
            basis: basis.len(),
            coefficients: coefficients.len(),
        });
    }
    let mut result = algebra.zero();
    for (basis_element, coefficient) in basis.iter().zip(coefficients) {
        result = result.add(&basis_element.scale_left(coefficient)?)?;
    }
    Ok(result)
}

fn exact_matrix_order<F: QpaCoefficientField>(
    matrix: &DenseMatrix<F>,
    field: &F,
) -> Result<NakayamaAutomorphismOrder, EnvelopingModuleError> {
    if matrix.nrows() == 0 || matrix.nrows() != matrix.ncols() {
        return Err(EnvelopingModuleError::InvalidNakayamaOrderMatrixShape {
            shape: (matrix.nrows(), matrix.ncols()),
        });
    }
    let exponent = match field.qpa_characteristic() {
        QpaFieldCharacteristic::Zero => rational_finite_matrix_exponent(matrix.nrows())?,
        QpaFieldCharacteristic::Prime(prime) => {
            prime_field_general_linear_exponent(matrix.nrows(), prime)?
        }
    };
    if !matrix_is_identity(&matrix_power_integer(matrix, &exponent)?) {
        return Ok(NakayamaAutomorphismOrder::Infinite);
    }
    let mut order = exponent;
    for (prime, _) in order.factor() {
        loop {
            if (&order % &prime) != 0 {
                break;
            }
            let candidate = &order / &prime;
            if matrix_is_identity(&matrix_power_integer(matrix, &candidate)?) {
                order = candidate;
            } else {
                break;
            }
        }
    }
    Ok(NakayamaAutomorphismOrder::Finite(order))
}

fn rational_finite_matrix_exponent(dimension: usize) -> Result<Integer, EnvelopingModuleError> {
    let mut prime = 2usize;
    let mut bound = 2usize;
    let mut totient_product = 1usize;
    while totient_product < dimension {
        prime = next_prime(prime)
            .ok_or(EnvelopingModuleError::NakayamaOrderCyclotomicBoundOverflow { dimension })?;
        bound = bound
            .checked_mul(prime)
            .ok_or(EnvelopingModuleError::NakayamaOrderCyclotomicBoundOverflow { dimension })?;
        totient_product = totient_product
            .checked_mul(prime - 1)
            .ok_or(EnvelopingModuleError::NakayamaOrderCyclotomicBoundOverflow { dimension })?;
    }
    let mut exponent = Integer::one();
    for index in 1..=bound {
        if euler_totient(index) <= dimension {
            exponent = exponent.lcm(&Integer::from(index));
        }
    }
    Ok(exponent)
}

fn prime_field_general_linear_exponent(
    dimension: usize,
    prime: u64,
) -> Result<Integer, EnvelopingModuleError> {
    let prime_usize = usize::try_from(prime)
        .map_err(|_| EnvelopingModuleError::NakayamaOrderPrimeTooLarge { prime })?;
    let prime_integer = Integer::from(prime);
    let mut exponent = Integer::one();
    let mut prime_power = Integer::one();
    for _ in 1..=dimension {
        prime_power *= &prime_integer;
        exponent = exponent.lcm(&(&prime_power - 1));
    }
    let mut unipotent_exponent = 0u64;
    let mut covered = 1usize;
    while covered < dimension {
        covered = covered
            .checked_mul(prime_usize)
            .ok_or(EnvelopingModuleError::NakayamaOrderPrimePowerOverflow { dimension, prime })?;
        unipotent_exponent += 1;
    }
    Ok(exponent * prime_integer.pow(unipotent_exponent))
}

fn matrix_power_integer<F: Field>(
    matrix: &DenseMatrix<F>,
    exponent: &Integer,
) -> Result<DenseMatrix<F>, EnvelopingModuleError> {
    let mut exponent = exponent.clone();
    let mut result = DenseMatrix::identity(matrix.nrows(), matrix.field().clone());
    let mut base = matrix.clone();
    while exponent > 0 {
        if &exponent % 2 == 1 {
            result = result.multiply(&base)?;
        }
        exponent >>= 1_u8;
        if exponent > 0 {
            base = base.multiply(&base)?;
        }
    }
    Ok(result)
}

fn matrix_is_identity<F: Field>(matrix: &DenseMatrix<F>) -> bool {
    matrix == &DenseMatrix::identity(matrix.nrows(), matrix.field().clone())
}

fn next_prime(current: usize) -> Option<usize> {
    let mut candidate = current.checked_add(1)?;
    loop {
        if is_prime(candidate) {
            return Some(candidate);
        }
        candidate = candidate.checked_add(1)?;
    }
}

fn is_prime(value: usize) -> bool {
    if value < 2 {
        return false;
    }
    if value.is_multiple_of(2) {
        return value == 2;
    }
    let mut divisor = 3usize;
    while divisor <= value / divisor {
        if value.is_multiple_of(divisor) {
            return false;
        }
        divisor += 2;
    }
    true
}

fn euler_totient(mut value: usize) -> usize {
    let mut result = value;
    let mut prime = 2usize;
    while prime <= value / prime {
        if value.is_multiple_of(prime) {
            while value.is_multiple_of(prime) {
                value /= prime;
            }
            result -= result / prime;
        }
        prime += if prime == 2 { 1 } else { 2 };
    }
    if value > 1 {
        result -= result / value;
    }
    result
}

fn order_from_cache(value: CachedNakayamaAutomorphismOrder) -> NakayamaAutomorphismOrder {
    match value {
        CachedNakayamaAutomorphismOrder::Finite(value) => NakayamaAutomorphismOrder::Finite(value),
        CachedNakayamaAutomorphismOrder::Infinite => NakayamaAutomorphismOrder::Infinite,
    }
}

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
    if !element.module().is_identical(target) {
        return Err(EnvelopingModuleError::ProjectiveElementInWrongModule);
    }
    let support = element.support_vertices();
    let [vertex] = support.as_slice() else {
        return Ok(None);
    };
    let bases = basis_of_projectives(target.right_acting_algebra())?
        .ok_or(EnvelopingModuleError::ProjectiveBasisUnavailable)?;
    let projectives = indec_projective_modules(target.right_acting_algebra())?;
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
        let maps = generators
            .iter()
            .map(|generator| {
                hom_from_projective(generator, module)?
                    .ok_or(EnvelopingModuleError::MinimalGeneratorNotVertexUniform)
            })
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
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

/// Translate QPA's `LiftingMorphismFromProjective(f, g)`.
///
/// `Ok(None)` is QPA's `fail` when the projective source was not constructed
/// with `DirectSumOfQPAModules`; validation errors retain the source method's
/// order.
///
/// # Errors
///
/// Requires a projective source for `g`, a surjective `f`, and a common range,
/// then propagates exact preimage, projective-map, and map arithmetic failures.
pub fn lifting_morphism_from_projective<F: Field>(
    first: &PathAlgebraModuleHomomorphism<F>,
    second: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<PathAlgebraModuleHomomorphism<F>>, EnvelopingModuleError> {
    let projective = second.source();
    if !is_projective_module(projective)? {
        return Err(EnvelopingModuleError::LiftingSourceNotProjective);
    }
    if !first.is_surjective() {
        return Err(EnvelopingModuleError::LiftingFirstMapNotSurjective);
    }
    if first.range() != second.range() {
        return Err(EnvelopingModuleError::LiftingMapsDoNotHaveCommonRange);
    }
    if projective.dimension() == 0 {
        return Ok(Some(zero_module_mapping(
            projective.clone(),
            first.source().clone(),
        )?));
    }
    if !projective.is_direct_sum_of_modules() {
        return Ok(None);
    }
    let inclusions = qpa_core::direct_sum_inclusions(projective)?;
    let projections = direct_sum_projections(projective)?;
    let mut partial_liftings = Vec::with_capacity(inclusions.len());
    for inclusion in &inclusions {
        let generator = minimal_generating_set_of_module(inclusion.source().clone())?
            .into_iter()
            .next()
            .ok_or(EnvelopingModuleError::EmptyProjectiveSummandGenerators)?;
        let projective_element = inclusion.image(&generator)?;
        let image = second.image(&projective_element)?;
        let generated = if image.is_zero() {
            zero_module_mapping(inclusion.source().clone(), first.source().clone())?
        } else {
            let preimage = first
                .preimage_representative(&image)?
                .ok_or(EnvelopingModuleError::MissingProjectiveLiftPreimage)?;
            let map = hom_from_projective(&preimage, first.source())?
                .ok_or(EnvelopingModuleError::ProjectiveLiftNotVertexUniform)?;
            right_module_hom_over_algebra(
                inclusion.source().clone(),
                first.source().clone(),
                map.vertex_matrices().to_vec(),
            )?
        };
        partial_liftings.push(generated);
    }
    let mut composites = projections
        .iter()
        .zip(&partial_liftings)
        .map(|(projection, partial)| Ok(projection.compose(partial)?))
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?
        .into_iter();
    let mut lifting = composites
        .next()
        .ok_or(EnvelopingModuleError::EmptyProjectiveSummandGenerators)?;
    for composite in composites {
        lifting = lifting.add(&composite)?;
    }
    Ok(Some(lifting))
}

/// Translate QPA's `HomFactoringThroughProjOverAlgebra(M, N)`.
///
/// The returned basis is GAP's insertion-ordered vector-space basis of raw
/// physical map lists after composing `Hom(M, P(N))` with the projective
/// cover of `N` and removing zero rows.
///
/// # Errors
///
/// Requires a common acting algebra and propagates projective-cover, Hom,
/// exact matrix, and checked-map construction failures.
pub fn hom_factoring_through_projectives<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
) -> Result<Vec<PathAlgebraModuleHomomorphism<F>>, EnvelopingModuleError> {
    if source.right_acting_algebra() != range.right_acting_algebra() {
        return Err(EnvelopingModuleError::FactoringModulesOverDifferentAlgebras);
    }
    if source.dimension() == 0 || range.dimension() == 0 {
        return Ok(Vec::new());
    }
    let cover = projective_cover(range)?;
    let homomorphisms = hom_over_algebra(source, cover.source())?
        .ok_or(EnvelopingModuleError::InternalProjectiveHomRejected)?;
    let composites = homomorphisms
        .iter()
        .map(|homomorphism| homomorphism.compose(&cover))
        .collect::<Result<Vec<_>, _>>()?;
    let nonzero = composites
        .iter()
        .filter(|homomorphism| !homomorphism.is_zero())
        .collect::<Vec<_>>();
    if nonzero.is_empty() {
        return Ok(Vec::new());
    }
    let field = source.right_acting_algebra().field().clone();
    let row_length = nonzero[0].flat_vertex_matrix_entries().len();
    let rows = DenseMatrix::from_flat(
        nonzero.len(),
        row_length,
        nonzero
            .iter()
            .flat_map(|homomorphism| homomorphism.flat_vertex_matrix_entries())
            .collect(),
        field.clone(),
    )?;
    let basis = rows.semi_echelon_transformation().vectors;
    (0..basis.nrows())
        .map(|row| module_map_from_physical_flat(source, range, basis.row(row), &field))
        .collect()
}

/// Translate QPA's `EndModuloProjOverAlgebra(M)` natural algebra
/// homomorphism.
///
/// This constructs `EndOverAlgebra(M)`, embeds QPA's ordered basis of maps
/// factoring through projectives as block matrices, forms their two-sided
/// ideal, and returns the natural quotient map onto the stable endomorphism
/// algebra.
///
/// # Errors
///
/// QPA rejects the zero module. Otherwise projective factoring,
/// endomorphism embedding, ideal closure, and exact quotient failures
/// propagate.
pub fn end_modulo_projectives<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<ModuleEndomorphismQuotientHomomorphism<F>, EnvelopingModuleError> {
    if module.dimension() == 0 {
        return Err(EnvelopingModuleError::ZeroStableEndomorphismModule);
    }
    let endomorphisms = end_over_algebra(module)?;
    let factor_matrices = hom_factoring_through_projectives(module, module)?
        .iter()
        .map(from_hom_mm_to_end_m)
        .collect::<Result<Vec<_>, _>>()?;
    Ok(natural_endomorphism_quotient_by_ideal(
        &endomorphisms,
        &factor_matrices,
    )?)
}

fn module_map_from_physical_flat<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
    entries: &[F::Element],
    field: &F,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let expected = source
        .dimension_vector()
        .iter()
        .zip(range.dimension_vector())
        .try_fold(0_usize, |total, (&source_dimension, &range_dimension)| {
            total.checked_add(
                source_dimension
                    .max(1)
                    .checked_mul(range_dimension.max(1))?,
            )
        })
        .ok_or(EnvelopingModuleError::PhysicalMapSizeOverflow)?;
    if entries.len() != expected {
        return Err(EnvelopingModuleError::PhysicalMapLength {
            actual: entries.len(),
            expected,
        });
    }
    let mut offset = 0_usize;
    let matrices = source
        .dimension_vector()
        .iter()
        .zip(range.dimension_vector())
        .map(|(&source_dimension, &range_dimension)| {
            let rows = source_dimension.max(1);
            let columns = range_dimension.max(1);
            let length = rows
                .checked_mul(columns)
                .ok_or(EnvelopingModuleError::PhysicalMapSizeOverflow)?;
            let matrix = DenseMatrix::from_flat(
                rows,
                columns,
                entries[offset..offset + length].to_vec(),
                field.clone(),
            )?;
            offset += length;
            Ok(matrix)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    Ok(right_module_hom_over_algebra(
        source.clone(),
        range.clone(),
        matrices,
    )?)
}

/// Translate QPA's cached `DualOfModule(M)` for path-algebra matrix modules.
///
/// Arrow matrices are transposed without reordering, the dimension vector is
/// retained, and the result acts over QPA's opposite algebra. QPA's zero
/// branch returns the opposite algebra's cached [`zero_module`]. For a
/// nonzero module, QPA explicitly installs the reciprocal dual attribute on
/// the result and propagates already-true indecomposable, projective, and
/// injective properties in the corresponding dual directions; its
/// zero-module early return does not.
///
/// # Errors
///
/// Propagates opposite-algebra construction and the checked right-module
/// constructor. A generated action failing the opposite relations is reported
/// as an internal source-translation failure.
pub fn dual_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, ModuleFunctorError> {
    if let Some(cached) = module.cached_dual_module() {
        return Ok(cached);
    }
    let opposite = match module.right_acting_algebra() {
        QuiverAlgebra::Path(algebra) => QuiverAlgebra::Path(algebra.opposite()?),
        QuiverAlgebra::Quotient(algebra) => {
            QuiverAlgebra::Quotient(opposite_path_algebra(algebra)?)
        }
    };
    if module.dimension() == 0 {
        let dual = zero_module(&opposite)?;
        return Ok(module.cache_dual_module(&dual));
    }
    let dimensions = module.dimension_vector().to_vec();
    let mut specification = Vec::new();
    for (position, arrow) in opposite.quiver().arrows().enumerate() {
        let source = dimensions[arrow.source().index().get() - 1];
        let target = dimensions[arrow.target().index().get() - 1];
        if source != 0 && target != 0 {
            specification.push(NamedArrowMatrix::new(
                arrow.name(),
                module.arrow_matrices()[position].transpose(),
            ));
        }
    }
    let dual = right_module_over_path_algebra(opposite, dimensions, &specification)?
        .ok_or(ModuleFunctorError::GeneratedDualViolatesRelations)?;
    let dual = module.install_dual_module(&dual)?;
    if module.cached_is_indecomposable_module() == Some(true) {
        let _ = dual.cache_is_indecomposable_module(true);
    }
    if module.cached_is_projective_module() == Some(true) {
        let _ = dual.cache_is_injective_module(true);
    }
    if module.cached_is_injective_module() == Some(true) {
        let _ = dual.cache_is_projective_module(true);
    }
    Ok(dual)
}

/// Translate QPA's cached `TransposeOfModule(M)` construction.
///
/// This follows `functors.gi` directly: expand the projective-cover map in
/// the global canonical algebra basis, take its GAP row nullspace, quotient
/// the syzygy by the arrow-generated radical with GAP's mutable-basis
/// complement order, transpose and opposite every resulting presentation
/// entry, and form the cokernel inside the corresponding opposite
/// projective cover.
///
/// # Errors
///
/// Requires the same finite projective/canonical-basis domain as QPA and
/// propagates exact algebra, nullspace, direct-sum, subrepresentation, and
/// cokernel failures.
#[allow(clippy::too_many_lines)]
pub fn transpose_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_transpose_of_module() {
        return Ok(cached);
    }
    let algebra = module.right_acting_algebra();
    let opposite = match algebra {
        QuiverAlgebra::Path(path_algebra) => QuiverAlgebra::Path(path_algebra.opposite()?),
        QuiverAlgebra::Quotient(quotient) => {
            QuiverAlgebra::Quotient(opposite_path_algebra(quotient)?)
        }
    };
    if is_projective_module(module)? {
        let zero = zero_module(&opposite)?;
        return Ok(module.cache_transpose_of_module(&zero));
    }

    let basis =
        finite_canonical_basis(algebra)?.ok_or(EnvelopingModuleError::InfiniteTransposeAlgebra)?;
    let field = algebra.field().clone();
    let vertex_elements = algebra
        .quiver()
        .vertices()
        .map(|vertex| include_in_path_algebra(&vertex.path(), algebra))
        .collect::<Result<Vec<_>, _>>()?;
    let arrow_elements = algebra
        .quiver()
        .arrows()
        .map(|arrow| include_in_path_algebra(&arrow.path(), algebra))
        .collect::<Result<Vec<_>, _>>()?;

    // QPA's local BB uses global CanonicalBasis(A) order, filtered by left
    // vertex multiplication; it is intentionally not the target-grouped
    // BasisOfProjectives order.
    let mut projective_basis_indices = vec![Vec::new(); vertex_elements.len()];
    for (vertex, vertex_element) in vertex_elements.iter().enumerate() {
        for (basis_index, basis_element) in basis.iter().enumerate() {
            if !vertex_element.multiply(basis_element)?.is_zero() {
                projective_basis_indices[vertex].push(basis_index);
            }
        }
    }

    // CanonicalBasis(M) is assigned in the source even though the value is
    // not subsequently read.
    let _module_basis = module.basis_elements();
    let generators = minimal_generating_set_of_module(module.clone())?;
    let mut cover_vertices = Vec::with_capacity(generators.len());
    let mut cover_basis_images = Vec::new();
    for generator in &generators {
        let mut supported_vertex = None;
        for (vertex, vertex_element) in vertex_elements.iter().enumerate() {
            if !generator.act(vertex_element)?.is_zero() {
                supported_vertex = Some(vertex);
                break;
            }
        }
        let vertex =
            supported_vertex.ok_or(EnvelopingModuleError::MissingTransposeGeneratorSupport)?;
        cover_vertices.push(vertex);
        for &basis_index in &projective_basis_indices[vertex] {
            cover_basis_images.push(generator.act(&basis[basis_index])?);
        }
    }

    let cover_matrix = DenseMatrix::from_rows(
        cover_basis_images
            .iter()
            .map(PathAlgebraModuleElement::flat_external_representation)
            .collect(),
        field.clone(),
    )?;
    let solution_rows = cover_matrix.nullspace_rows();
    if solution_rows.nrows() == 0 {
        let zero = zero_module(&opposite)?;
        return Ok(module.cache_transpose_of_module(&zero));
    }
    let solution_vectors = projective_algebra_vectors_from_rows(
        algebra,
        &basis,
        &projective_basis_indices,
        &cover_vertices,
        (0..solution_rows.nrows()).map(|row| solution_rows.row(row)),
    )?;

    let mut radical_vectors = Vec::new();
    for solution in &solution_vectors {
        for arrow in &arrow_elements {
            let product = solution
                .iter()
                .map(|entry| entry.multiply(arrow))
                .collect::<Result<Vec<_>, _>>()?;
            if product.iter().any(|entry| !entry.is_zero()) {
                radical_vectors.push(product);
            }
        }
    }
    let radical_rows = radical_vectors
        .iter()
        .map(|vector| {
            let mut row = Vec::with_capacity(solution_rows.ncols());
            for (entry, &cover_vertex) in vector.iter().zip(&cover_vertices) {
                let coefficients = canonical_basis_coefficients(algebra, &basis, entry)?;
                row.extend(
                    projective_basis_indices[cover_vertex]
                        .iter()
                        .map(|&basis_index| coefficients[basis_index].clone()),
                );
            }
            Ok(row)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;

    let syzygy_basis = solution_rows.semi_echelon_transformation().vectors;
    let radical_basis = if radical_rows.is_empty() {
        DenseMatrix::from_flat(0, solution_rows.ncols(), Vec::new(), field.clone())?
    } else {
        DenseMatrix::from_rows(radical_rows, field.clone())?
            .semi_echelon_transformation()
            .vectors
    };
    let top_rows = if radical_basis.nrows() == 0 {
        // GAP's high-rank trivial-subspace method returns IdentityMapping(V).
        // The source then asks for CanonicalBasis(Range(f)) = CanonicalBasis(V),
        // which uses full triangulation rather than the ordinary Basis(V)
        // used by the nontrivial quotient method.
        let canonical = solution_rows.canonical_row_space_basis();
        (0..canonical.nrows())
            .map(|row| canonical.row(row).to_vec())
            .collect()
    } else {
        quotient_complement_in_space(&syzygy_basis, &radical_basis)?.representatives
    };
    let top_vectors = projective_algebra_vectors_from_rows(
        algebra,
        &basis,
        &projective_basis_indices,
        &cover_vertices,
        top_rows.iter().map(Vec::as_slice),
    )?;

    let projectives = indec_projective_modules(algebra)?;
    let cover_summands = cover_vertices
        .iter()
        .map(|&vertex| projectives[vertex].clone())
        .collect::<Vec<_>>();
    let _projective_cover = direct_sum_of_qpa_modules(&cover_summands)?
        .ok_or(EnvelopingModuleError::InternalTransposeDirectSumFailed)?;

    let mut syzygy_generator_vertices = Vec::with_capacity(top_vectors.len());
    for generator in &top_vectors {
        let mut support = Vec::new();
        for (vertex, vertex_element) in vertex_elements.iter().enumerate() {
            let supported = generator.iter().try_fold(false, |supported, entry| {
                Ok::<_, EnvelopingModuleError>(
                    supported || !entry.multiply(vertex_element)?.is_zero(),
                )
            })?;
            if supported {
                support.push(vertex);
            }
        }
        let [vertex] = support.as_slice() else {
            return Err(EnvelopingModuleError::NonuniformTransposeSyzygyGenerator);
        };
        syzygy_generator_vertices.push(*vertex);
    }

    let opposite_projectives = indec_projective_modules(&opposite)?;
    let opposite_cover_summands = syzygy_generator_vertices
        .iter()
        .map(|&vertex| opposite_projectives[vertex].clone())
        .collect::<Vec<_>>();
    let opposite_cover = direct_sum_of_qpa_modules(&opposite_cover_summands)?
        .ok_or(EnvelopingModuleError::InternalTransposeDirectSumFailed)?;
    let inclusions = direct_sum_inclusions(&opposite_cover)?;
    let embedded_generators = inclusions
        .iter()
        .zip(&opposite_cover_summands)
        .map(|(inclusion, summand)| {
            let generator = minimal_generating_set_of_module(summand.clone())?
                .into_iter()
                .next()
                .ok_or(EnvelopingModuleError::EmptyTransposeProjectiveGenerator)?;
            Ok(inclusion.image(&generator)?)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;

    let mut image_generators = Vec::with_capacity(generators.len());
    for cover_generator in 0..generators.len() {
        let mut image = opposite_cover.zero_element();
        for (top_generator, embedded_generator) in top_vectors.iter().zip(&embedded_generators) {
            let opposite_entry = opposite_quiver_algebra_element(&top_generator[cover_generator])?;
            image = image.add(&embedded_generator.act(&opposite_entry)?)?;
        }
        image_generators.push(image);
    }
    let inclusion = subrepresentation_inclusion(opposite_cover, &image_generators)?;
    let transpose = cokernel_of_module_homomorphism(&inclusion)?;
    if transpose.dimension() != 0 && module.cached_is_indecomposable_module() == Some(true) {
        let _ = transpose.cache_is_indecomposable_module(true);
    }
    Ok(module.cache_transpose_of_module(&transpose))
}

fn projective_algebra_vectors_from_rows<'a, F: Field + 'a>(
    algebra: &QuiverAlgebra<F>,
    basis: &[QuiverAlgebraElement<F>],
    projective_basis_indices: &[Vec<usize>],
    cover_vertices: &[usize],
    rows: impl IntoIterator<Item = &'a [F::Element]>,
) -> Result<Vec<Vec<QuiverAlgebraElement<F>>>, EnvelopingModuleError>
where
    F::Element: 'a,
{
    rows.into_iter()
        .map(|row| {
            let mut offset = 0;
            cover_vertices
                .iter()
                .map(|&vertex| {
                    let mut element = algebra.zero();
                    for &basis_index in &projective_basis_indices[vertex] {
                        let term = basis[basis_index].scale_left(&row[offset])?;
                        element = element.add(&term)?;
                        offset += 1;
                    }
                    Ok(element)
                })
                .collect::<Result<Vec<_>, EnvelopingModuleError>>()
        })
        .collect()
}

fn opposite_quiver_algebra_element<F: Field>(
    element: &QuiverAlgebraElement<F>,
) -> Result<QuiverAlgebraElement<F>, EnvelopingModuleError> {
    match element {
        QuiverAlgebraElement::Path(element) => Ok(QuiverAlgebraElement::Path(element.opposite()?)),
        QuiverAlgebraElement::Quotient(element) => Ok(QuiverAlgebraElement::Quotient(
            opposite_quotient_element(element)?,
        )),
    }
}

/// Translate QPA's cached module `DTr(M) = DualOfModule(TransposeOfModule(M))`.
///
/// # Errors
///
/// Propagates transpose and dual construction.
pub fn d_tr<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_d_tr() {
        return Ok(cached);
    }
    let transpose = transpose_of_module(module)?;
    let result = dual_of_module(&transpose)?;
    Ok(module.cache_d_tr(&result))
}

/// Translate QPA's cached module `TrD(M) = TransposeOfModule(DualOfModule(M))`.
///
/// # Errors
///
/// Propagates dual and transpose construction.
pub fn tr_d<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_tr_d() {
        return Ok(cached);
    }
    let dual = dual_of_module(module)?;
    let result = transpose_of_module(&dual)?;
    Ok(module.cache_tr_d(&result))
}

/// Translate QPA's integer-power `DTr(M, n)` dispatch.
///
/// Negative powers apply [`tr_d`], while zero returns the identical entered
/// module.
///
/// # Errors
///
/// Propagates every iteration's functor construction.
pub fn d_tr_iterated<F: Field>(
    module: &PathAlgebraModule<F>,
    exponent: i64,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    let mut result = module.clone();
    if exponent < 0 {
        for _ in 0..exponent.unsigned_abs() {
            result = tr_d(&result)?;
        }
    } else {
        for _ in 0..exponent.unsigned_abs() {
            result = d_tr(&result)?;
        }
    }
    Ok(result)
}

/// Translate QPA's integer-power `TrD(M, n)` dispatch.
///
/// Negative powers apply [`d_tr`], while zero returns the identical entered
/// module.
///
/// # Errors
///
/// Propagates every iteration's functor construction.
pub fn tr_d_iterated<F: Field>(
    module: &PathAlgebraModule<F>,
    exponent: i64,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    let mut result = module.clone();
    if exponent < 0 {
        for _ in 0..exponent.unsigned_abs() {
            result = d_tr(&result)?;
        }
    } else {
        for _ in 0..exponent.unsigned_abs() {
            result = tr_d(&result)?;
        }
    }
    Ok(result)
}

/// Translate QPA's cached `IsTauRigidModule(M)` property.
///
/// This is exactly `Length(HomOverAlgebra(M, DTr(M))) = 0`.
///
/// # Errors
///
/// Propagates `DTr` and Hom-space construction.
pub fn is_tau_rigid_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<bool, EnvelopingModuleError> {
    if let Some(cached) = module.cached_is_tau_rigid_module() {
        return Ok(cached);
    }
    let translate = d_tr(module)?;
    let hom = hom_over_algebra(module, &translate)?
        .ok_or(EnvelopingModuleError::InternalTauRigidHomRejected)?;
    Ok(module.cache_is_tau_rigid_module(hom.is_empty()))
}

/// Translate QPA's cached `DualOfModuleHomomorphism(f)` attribute.
///
/// The source transposes every vertex map in place, dualizes the source and
/// range modules, reverses them, and passes the matrices through the checked
/// module-homomorphism constructor. QPA caches only the entered map's
/// attribute; it does not install a reciprocal cache on the returned map.
///
/// # Errors
///
/// Propagates module dualization and checked map construction.
pub fn dual_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, ModuleFunctorError> {
    if let Some(cached) = homomorphism.cached_dual_homomorphism() {
        return Ok(cached);
    }
    let matrices = homomorphism
        .vertex_matrices()
        .iter()
        .map(DenseMatrix::transpose)
        .collect();
    let dual_source = dual_of_module(homomorphism.source())?;
    let dual_range = dual_of_module(homomorphism.range())?;
    let dual = right_module_hom_over_algebra(dual_range, dual_source, matrices)?;
    Ok(homomorphism.cache_dual_homomorphism(&dual))
}

/// Translate QPA's cached `TransposeOfModuleHomomorphism(f)` attribute.
///
/// This follows QPA's projective-presentation construction literally. It
/// lifts `f` to the projective covers and their first syzygy covers, applies
/// the contravariant star functor to all four maps, and takes the induced map
/// on cokernels.
///
/// # Errors
///
/// Propagates projective-cover, kernel, lifting, star, and cokernel-map
/// construction failures. The internal `fail` variants can occur only if a
/// square produced by the source construction ceases to commute.
pub fn transpose_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_transpose_homomorphism() {
        return Ok(cached);
    }

    let d0 = projective_cover(homomorphism.source())?;
    let d0_prime = projective_cover(homomorphism.range())?;
    let ker_d0 = kernel_inclusion(&d0)?;
    let ker_d0_prime = kernel_inclusion(&d0_prime)?;
    let d1_temp = projective_cover(ker_d0.source())?;
    let d1_prime_temp = projective_cover(ker_d0_prime.source())?;
    let d1 = d1_temp.compose(&ker_d0)?;
    let d1_prime = d1_prime_temp.compose(&ker_d0_prime)?;

    let f0 = lifting_morphism_from_projective(&d0_prime, &d0.compose(homomorphism)?)?
        .ok_or(EnvelopingModuleError::InternalTransposeHomomorphismLiftFailed)?;
    let f1_prime = morphism_on_kernel(&d0, &d0_prime, &f0, homomorphism)?
        .ok_or(EnvelopingModuleError::InternalTransposeHomomorphismKernelSquareFailed)?;
    let f1 = lifting_morphism_from_projective(&d1_prime_temp, &d1_temp.compose(&f1_prime)?)?
        .ok_or(EnvelopingModuleError::InternalTransposeHomomorphismLiftFailed)?;

    let d1_star = star_of_module_homomorphism(&d1)?;
    let d1_prime_star = star_of_module_homomorphism(&d1_prime)?;
    let f1_star = star_of_module_homomorphism(&f1)?;
    let f0_star = star_of_module_homomorphism(&f0)?;
    let transpose = morphism_on_cokernel(&d1_prime_star, &d1_star, &f0_star, &f1_star)?
        .ok_or(EnvelopingModuleError::InternalTransposeHomomorphismCokernelSquareFailed)?;

    Ok(homomorphism.cache_transpose_homomorphism(&transpose))
}

/// Translate QPA's cached `DTr(f)` attribute on module homomorphisms.
///
/// QPA computes the dual of the transposed map, then reconstructs it over the
/// exact cached `DTr` source and range module objects.
///
/// # Errors
///
/// Propagates transpose, dual, module `DTr`, and checked-map construction.
pub fn d_tr_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_d_tr_homomorphism() {
        return Ok(cached);
    }
    let transpose = transpose_of_module_homomorphism(homomorphism)?;
    let generated = dual_of_module_homomorphism(&transpose)?;
    let source = d_tr(homomorphism.source())?;
    let range = d_tr(homomorphism.range())?;
    let result =
        right_module_hom_over_algebra(source, range, generated.vertex_matrices().to_vec())?;
    Ok(homomorphism.cache_d_tr_homomorphism(&result))
}

/// Translate QPA's cached `TrD(f)` attribute on module homomorphisms.
///
/// QPA transposes the dual map, then reconstructs it over the exact cached
/// `TrD` source and range module objects.
///
/// # Errors
///
/// Propagates dual, transpose, module `TrD`, and checked-map construction.
pub fn tr_d_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_tr_d_homomorphism() {
        return Ok(cached);
    }
    let dual = dual_of_module_homomorphism(homomorphism)?;
    let generated = transpose_of_module_homomorphism(&dual)?;
    let source = tr_d(homomorphism.source())?;
    let range = tr_d(homomorphism.range())?;
    let result =
        right_module_hom_over_algebra(source, range, generated.vertex_matrices().to_vec())?;
    Ok(homomorphism.cache_tr_d_homomorphism(&result))
}

/// Translate QPA's integer-power `DTr(f, n)` method.
///
/// Negative powers apply [`tr_d_homomorphism`], while zero returns the
/// identical entered mapping object.
///
/// # Errors
///
/// Propagates each iterated homomorphism-functor construction.
pub fn d_tr_homomorphism_iterated<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    exponent: i64,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let mut result = homomorphism.clone();
    if exponent < 0 {
        for _ in 0..exponent.unsigned_abs() {
            result = tr_d_homomorphism(&result)?;
        }
    } else {
        for _ in 0..exponent.unsigned_abs() {
            result = d_tr_homomorphism(&result)?;
        }
    }
    Ok(result)
}

/// Translate QPA's integer-power `TrD(f, n)` method.
///
/// Negative powers apply [`d_tr_homomorphism`], while zero returns the
/// identical entered mapping object.
///
/// # Errors
///
/// Propagates each iterated homomorphism-functor construction.
pub fn tr_d_homomorphism_iterated<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
    exponent: i64,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    let mut result = homomorphism.clone();
    if exponent < 0 {
        for _ in 0..exponent.unsigned_abs() {
            result = d_tr_homomorphism(&result)?;
        }
    } else {
        for _ in 0..exponent.unsigned_abs() {
            result = tr_d_homomorphism(&result)?;
        }
    }
    Ok(result)
}

/// Translate QPA's cached `StarOfModule(M) = Hom_A(M, A)` construction.
///
/// The vertex at `i` is represented by QPA's ordered basis of
/// `Hom_A(M, e_i A)`.  For every arrow `i -> j`, the opposite-arrow action is
/// obtained by composing `Hom_A(M, e_j A)` with QPA's map `e_j A -> e_i A`
/// induced by left multiplication by that arrow, then solving in the original
/// raw physical Hom basis.
///
/// # Errors
///
/// Requires QPA's finite projective bases and propagates opposite-algebra,
/// Hom, projective-map, exact coordinate, and checked-module failures.
pub fn star_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_star_of_module() {
        return Ok(cached);
    }
    let algebra = module.right_acting_algebra();
    let opposite = match algebra {
        QuiverAlgebra::Path(path_algebra) => QuiverAlgebra::Path(path_algebra.opposite()?),
        QuiverAlgebra::Quotient(quotient) => {
            QuiverAlgebra::Quotient(opposite_path_algebra(quotient)?)
        }
    };
    let field = algebra.field().clone();
    let projectives = indec_projective_modules(algebra)?;
    let hom_spaces = projectives
        .iter()
        .map(|projective| {
            hom_over_algebra(module, projective)?
                .ok_or(EnvelopingModuleError::InternalStarHomRejected)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let hom_bases = hom_spaces
        .iter()
        .map(|basis| {
            let Some(first) = basis.first() else {
                return Ok(None);
            };
            Ok(Some(DenseMatrix::from_flat(
                basis.len(),
                first.flat_vertex_matrix_entries().len(),
                basis
                    .iter()
                    .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
                    .collect(),
                field.clone(),
            )?))
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let dimensions = hom_spaces.iter().map(Vec::len).collect::<Vec<_>>();

    // The source computes these subspace bases even though the local `Bproj`
    // value is not subsequently read.
    let _projective_bases =
        basis_of_projectives(algebra)?.ok_or(EnvelopingModuleError::ProjectiveBasisUnavailable)?;

    let arrows = algebra.quiver().arrows().collect::<Vec<_>>();
    let opposite_arrows = opposite.quiver().arrows().collect::<Vec<_>>();
    let mut left_multiplications = Vec::with_capacity(arrows.len());
    for arrow in &arrows {
        let source = arrow.source().index().get() - 1;
        let generator = minimal_generating_set_of_module(projectives[source].clone())?
            .into_iter()
            .next()
            .ok_or(EnvelopingModuleError::EmptyStarProjectiveGenerator)?;
        let arrow_element = include_in_path_algebra(&arrow.path(), algebra)?;
        let image = generator.act(&arrow_element)?;
        let left_multiplication = hom_from_projective(&image, &projectives[source])?
            .ok_or(EnvelopingModuleError::StarArrowMapNotUniform)?;
        left_multiplications.push(left_multiplication);
    }

    let mut specification = Vec::new();
    for (position, arrow) in arrows.iter().enumerate() {
        let source = arrow.source().index().get() - 1;
        let target = arrow.target().index().get() - 1;
        if dimensions[source] == 0 || dimensions[target] == 0 {
            continue;
        }
        let coordinate_basis = hom_bases[source]
            .as_ref()
            .ok_or(EnvelopingModuleError::MissingStarCoordinateBasis)?;
        let rows = hom_spaces[target]
            .iter()
            .map(|homomorphism| {
                let composite = homomorphism.compose(&left_multiplications[position])?;
                coordinate_basis
                    .solution_row(&composite.flat_vertex_matrix_entries())?
                    .ok_or(EnvelopingModuleError::StarCompositionOutsideHomBasis)
            })
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
        specification.push(NamedArrowMatrix::new(
            opposite_arrows[position].name(),
            DenseMatrix::from_rows(rows, field.clone())?,
        ));
    }
    let star = right_module_over_path_algebra(opposite, dimensions, &specification)?
        .ok_or(EnvelopingModuleError::GeneratedStarViolatesRelations)?;
    Ok(module.cache_star_of_module(&star))
}

/// Translate QPA's cached `FromIdentityToDoubleStarHomomorphism(M)` attribute.
///
/// At a basis vector `m` in vertex `i`, the source constructs the evaluation
/// map `Hom_A(M,A) -> e_i A` sending `g` to `g(m)`, converts its projective
/// values through the opposite algebra, and expresses that map in QPA's
/// ordered basis of `Hom_{A^op}(M^*, e_i A^op)`. Those coefficients are the
/// corresponding row of the natural map `M -> M^{**}`.
///
/// # Errors
///
/// Requires finite projective and canonical bases and propagates opposite
/// algebra, Hom-space, module-action, exact-coordinate, and checked-map
/// failures.
#[allow(clippy::too_many_lines)]
pub fn from_identity_to_double_star_homomorphism<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_identity_to_double_star() {
        return Ok(cached);
    }
    let algebra = module.right_acting_algebra();
    let opposite = opposite_module_algebra(algebra)?;
    let field = algebra.field().clone();
    let projectives = indec_projective_modules(algebra)?;
    let vertex_count = projectives.len();
    let hom_spaces = projectives
        .iter()
        .map(|projective| {
            hom_over_algebra(module, projective)?
                .ok_or(EnvelopingModuleError::InternalDoubleStarHomRejected)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let _hom_coordinate_bases = module_hom_flat_coordinate_bases(&hom_spaces, &field)?;

    let star = star_of_module(module)?;
    let opposite_projectives = indec_projective_modules(&opposite)?;
    let opposite_hom_spaces = opposite_projectives
        .iter()
        .map(|projective| {
            hom_over_algebra(&star, projective)?
                .ok_or(EnvelopingModuleError::InternalDoubleStarHomRejected)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let opposite_hom_coordinate_bases =
        module_hom_flat_coordinate_bases(&opposite_hom_spaces, &field)?;
    let double_star = star_of_module(&star)?;

    let opposite_projective_bases = basis_of_projectives(&opposite)?
        .ok_or(EnvelopingModuleError::ProjectiveBasisUnavailable)?;
    let opposite_canonical_basis = finite_canonical_basis(&opposite)?
        .ok_or(EnvelopingModuleError::MissingDoubleStarCanonicalBasis)?;
    let opposite_projective_basis_positions = opposite_projective_bases
        .iter()
        .enumerate()
        .map(|(projective, vertex_bases)| {
            vertex_bases
                .iter()
                .enumerate()
                .map(|(vertex, basis)| {
                    basis
                        .iter()
                        .enumerate()
                        .map(|(position, element)| {
                            opposite_canonical_basis
                                .iter()
                                .position(|candidate| candidate == element)
                                .ok_or(
                                    EnvelopingModuleError::MissingDoubleStarProjectiveBasisElement {
                                        projective: projective + 1,
                                        vertex: vertex + 1,
                                        position: position + 1,
                                    },
                                )
                        })
                        .collect::<Result<Vec<_>, EnvelopingModuleError>>()
                })
                .collect::<Result<Vec<_>, EnvelopingModuleError>>()
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;

    let module_basis = module.basis_elements();
    let module_dimensions = module.dimension_vector();
    let double_star_dimensions = double_star.dimension_vector();
    let mut matrices = Vec::with_capacity(vertex_count);
    let mut module_basis_offset = 0_usize;
    for projective_vertex in 0..vertex_count {
        let source_dimension = module_dimensions[projective_vertex];
        let range_dimension = double_star_dimensions[projective_vertex];
        if source_dimension == 0 || range_dimension == 0 {
            matrices.push(DenseMatrix::zero(
                source_dimension.max(1),
                range_dimension.max(1),
                field.clone(),
            ));
            module_basis_offset += source_dimension;
            continue;
        }
        let coordinate_basis = opposite_hom_coordinate_bases[projective_vertex]
            .as_ref()
            .ok_or(EnvelopingModuleError::MissingDoubleStarHomCoordinateBasis {
                vertex: projective_vertex + 1,
            })?;
        let mut rows = Vec::with_capacity(source_dimension);
        for module_element in
            &module_basis[module_basis_offset..module_basis_offset + source_dimension]
        {
            let mut evaluation_matrices = Vec::with_capacity(vertex_count);
            for (target_vertex, hom_basis) in hom_spaces.iter().enumerate() {
                let projective_dimension =
                    opposite_projectives[projective_vertex].dimension_vector()[target_vertex];
                if hom_basis.is_empty() || projective_dimension == 0 {
                    evaluation_matrices.push(DenseMatrix::zero(
                        hom_basis.len().max(1),
                        projective_dimension.max(1),
                        field.clone(),
                    ));
                    continue;
                }
                let value_rows = hom_basis
                    .iter()
                    .map(|homomorphism| {
                        let image = homomorphism.image(module_element)?;
                        let parent_element =
                            element_in_indec_projective(algebra, &image, target_vertex + 1)?;
                        let opposite_element =
                            opposite_projective_path_element(algebra, &parent_element)?;
                        let canonical_coordinates = canonical_basis_coefficients(
                            &opposite,
                            &opposite_canonical_basis,
                            &opposite_element,
                        )?;
                        Ok(
                            opposite_projective_basis_positions[projective_vertex][target_vertex]
                                .iter()
                                .map(|&position| canonical_coordinates[position].clone())
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
                evaluation_matrices.push(DenseMatrix::from_rows(value_rows, field.clone())?);
            }
            let evaluation = right_module_hom_over_algebra(
                star.clone(),
                opposite_projectives[projective_vertex].clone(),
                evaluation_matrices,
            )?;
            rows.push(
                coordinate_basis
                    .solution_row(&evaluation.flat_vertex_matrix_entries())?
                    .ok_or(EnvelopingModuleError::DoubleStarEvaluationOutsideHomBasis {
                        vertex: projective_vertex + 1,
                    })?,
            );
        }
        matrices.push(DenseMatrix::from_rows(rows, field.clone())?);
        module_basis_offset += source_dimension;
    }
    let evaluation = right_module_hom_over_algebra(module.clone(), double_star, matrices)?;
    Ok(module.cache_identity_to_double_star(&evaluation))
}

fn module_hom_flat_coordinate_bases<F: Field>(
    hom_spaces: &[Vec<PathAlgebraModuleHomomorphism<F>>],
    field: &F,
) -> Result<Vec<Option<DenseMatrix<F>>>, EnvelopingModuleError> {
    hom_spaces
        .iter()
        .map(|basis| {
            let Some(first) = basis.first() else {
                return Ok(None);
            };
            Ok(Some(DenseMatrix::from_flat(
                basis.len(),
                first.flat_vertex_matrix_entries().len(),
                basis
                    .iter()
                    .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
                    .collect(),
                field.clone(),
            )?))
        })
        .collect()
}

/// Translate QPA's cached contravariant `StarOfModuleHomomorphism(f)`.
///
/// # Errors
///
/// Propagates both star-module constructions, Hom composition, exact
/// coordinate solving, and the checked induced-map constructor.
pub fn star_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_star_homomorphism() {
        return Ok(cached);
    }
    let source = homomorphism.source();
    let range = homomorphism.range();
    let star_source = star_of_module(source)?;
    let star_range = star_of_module(range)?;
    let algebra = source.right_acting_algebra();
    let field = algebra.field().clone();
    let projectives = indec_projective_modules(algebra)?;

    let range_homs = projectives
        .iter()
        .map(|projective| {
            hom_over_algebra(range, projective)?
                .ok_or(EnvelopingModuleError::InternalStarHomRejected)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let source_homs = projectives
        .iter()
        .map(|projective| {
            hom_over_algebra(source, projective)?
                .ok_or(EnvelopingModuleError::InternalStarHomRejected)
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;

    let mut maps = Vec::with_capacity(projectives.len());
    for (vertex, (range_basis, source_basis)) in range_homs.iter().zip(&source_homs).enumerate() {
        if !range_basis.is_empty() && !source_basis.is_empty() {
            let coordinate_basis = DenseMatrix::from_flat(
                source_basis.len(),
                source_basis[0].flat_vertex_matrix_entries().len(),
                source_basis
                    .iter()
                    .flat_map(PathAlgebraModuleHomomorphism::flat_vertex_matrix_entries)
                    .collect(),
                field.clone(),
            )?;
            let rows = range_basis
                .iter()
                .map(|basis_map| {
                    let composite = homomorphism.compose(basis_map)?;
                    coordinate_basis
                        .solution_row(&composite.flat_vertex_matrix_entries())?
                        .ok_or(EnvelopingModuleError::StarCompositionOutsideHomBasis)
                })
                .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
            maps.push(DenseMatrix::from_rows(rows, field.clone())?);
        } else {
            maps.push(DenseMatrix::zero(
                star_range.dimension_vector()[vertex].max(1),
                star_source.dimension_vector()[vertex].max(1),
                field.clone(),
            ));
        }
    }
    let star = right_module_hom_over_algebra(star_range, star_source, maps)?;
    Ok(homomorphism.cache_star_homomorphism(&star))
}

/// Translate QPA's cached `NakayamaFunctorOfModule(M)` composition.
///
/// # Errors
///
/// Propagates [`star_of_module`] and [`dual_of_module`].
pub fn nakayama_functor_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_nakayama_functor_of_module() {
        return Ok(cached);
    }
    let star = star_of_module(module)?;
    let nakayama = dual_of_module(&star)?;
    Ok(module.cache_nakayama_functor_of_module(&nakayama))
}

/// Translate QPA's cached `OppositeNakayamaFunctorOfModule(M)` composition.
///
/// # Errors
///
/// Propagates [`dual_of_module`] and [`star_of_module`].
pub fn opposite_nakayama_functor_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = module.cached_opposite_nakayama_functor_of_module() {
        return Ok(cached);
    }
    let dual = dual_of_module(module)?;
    let nakayama = star_of_module(&dual)?;
    Ok(module.cache_opposite_nakayama_functor_of_module(&nakayama))
}

/// Translate QPA's cached Nakayama functor on a module homomorphism.
///
/// # Errors
///
/// Propagates star and dual homomorphism construction.
pub fn nakayama_functor_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_nakayama_functor_homomorphism() {
        return Ok(cached);
    }
    let star = star_of_module_homomorphism(homomorphism)?;
    let nakayama = dual_of_module_homomorphism(&star)?;
    Ok(homomorphism.cache_nakayama_functor_homomorphism(&nakayama))
}

/// Translate QPA's cached opposite Nakayama functor on a homomorphism.
///
/// # Errors
///
/// Propagates dual and star homomorphism construction.
pub fn opposite_nakayama_functor_of_module_homomorphism<F: Field>(
    homomorphism: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
    if let Some(cached) = homomorphism.cached_opposite_nakayama_functor_homomorphism() {
        return Ok(cached);
    }
    let dual = dual_of_module_homomorphism(homomorphism)?;
    let nakayama = star_of_module_homomorphism(&dual)?;
    Ok(homomorphism.cache_opposite_nakayama_functor_homomorphism(&nakayama))
}

struct NaturalTensorQuotient<F: Field> {
    relation_basis: DenseMatrix<F>,
    basis_lifts: DenseMatrix<F>,
    projection: DenseMatrix<F>,
}

fn natural_tensor_quotient<F: Field>(
    relation_rows: &DenseMatrix<F>,
) -> Result<NaturalTensorQuotient<F>, EnvelopingModuleError> {
    let field = relation_rows.field().clone();
    let ambient_dimension = relation_rows.ncols();
    let relation_basis = relation_rows.semi_echelon_transformation().vectors;
    let mut mutable_basis = relation_basis.clone();
    let mut complement = Vec::new();
    for coordinate in 0..ambient_dimension {
        let mut standard = vec![field.zero(); ambient_dimension];
        standard[coordinate] = field.one();
        let mut entries = mutable_basis.as_flat().to_vec();
        entries.extend_from_slice(&standard);
        let candidate = DenseMatrix::from_flat(
            mutable_basis.nrows() + 1,
            ambient_dimension,
            entries,
            field.clone(),
        )?
        .semi_echelon_transformation()
        .vectors;
        if candidate.nrows() > mutable_basis.nrows() {
            complement.push(standard);
            mutable_basis = candidate;
        }
    }
    if relation_basis.nrows() + complement.len() != ambient_dimension {
        return Err(EnvelopingModuleError::IncompleteTensorProductQuotientBasis);
    }

    let quotient_dimension = complement.len();
    let basis_lifts = DenseMatrix::from_flat(
        quotient_dimension,
        ambient_dimension,
        complement.iter().flatten().cloned().collect(),
        field.clone(),
    )?;
    let mut full_basis_entries = relation_basis.as_flat().to_vec();
    full_basis_entries.extend(complement.iter().flatten().cloned());
    let full_basis = DenseMatrix::from_flat(
        ambient_dimension,
        ambient_dimension,
        full_basis_entries,
        field.clone(),
    )?;
    let mut projection_entries = Vec::with_capacity(ambient_dimension * quotient_dimension);
    for coordinate in 0..ambient_dimension {
        let mut standard = vec![field.zero(); ambient_dimension];
        standard[coordinate] = field.one();
        let coefficients = full_basis
            .solution_row(&standard)?
            .ok_or(EnvelopingModuleError::IncompleteTensorProductQuotientBasis)?;
        projection_entries.extend_from_slice(
            &coefficients[ambient_dimension - quotient_dimension..ambient_dimension],
        );
    }
    let projection = DenseMatrix::from_flat(
        ambient_dimension,
        quotient_dimension,
        projection_entries,
        field,
    )?;
    Ok(NaturalTensorQuotient {
        relation_basis,
        basis_lifts,
        projection,
    })
}

fn tensor_dimension_pairing(
    left: &[usize],
    right: &[usize],
) -> Result<usize, EnvelopingModuleError> {
    left.iter().zip(right).try_fold(0_usize, |total, (&a, &b)| {
        total
            .checked_add(
                a.checked_mul(b)
                    .ok_or(EnvelopingModuleError::TensorProductDimensionOverflow)?,
            )
            .ok_or(EnvelopingModuleError::TensorProductDimensionOverflow)
    })
}

fn opposite_module_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<QuiverAlgebra<F>, EnvelopingModuleError> {
    match algebra {
        QuiverAlgebra::Path(algebra) => Ok(QuiverAlgebra::Path(algebra.opposite()?)),
        QuiverAlgebra::Quotient(algebra) => {
            Ok(QuiverAlgebra::Quotient(opposite_path_algebra(algebra)?))
        }
    }
}

fn opposite_projective_path_element<F: Field>(
    algebra: &QuiverAlgebra<F>,
    element: &PathAlgebraElement<F>,
) -> Result<QuiverAlgebraElement<F>, EnvelopingModuleError> {
    let opposite = element.opposite()?;
    match algebra {
        QuiverAlgebra::Path(_) => Ok(QuiverAlgebraElement::Path(opposite)),
        QuiverAlgebra::Quotient(algebra) => Ok(QuiverAlgebraElement::Quotient(
            opposite_path_algebra(algebra)?.element(opposite)?,
        )),
    }
}

fn tensor_projective_vertex<F: Field>(
    projective: &PathAlgebraModule<F>,
) -> Result<usize, EnvelopingModuleError> {
    let generator = minimal_generating_set_of_module(projective.clone())?
        .into_iter()
        .next()
        .ok_or(EnvelopingModuleError::EmptyTensorProjectiveGenerator)?;
    let support = generator.support_vertices();
    let [vertex] = support.as_slice() else {
        return Err(EnvelopingModuleError::NonuniformTensorProjectiveGenerator);
    };
    Ok(*vertex)
}

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

/// Failure while constructing QPA's regular bimodule over its enveloping
/// algebra.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum EnvelopingModuleError {
    /// A finite quotient must expose the canonical basis prepared by GBNP.
    #[error("a finite-dimensional quotient has no canonical basis vectors")]
    MissingFiniteCanonicalBasis,
    /// A path-algebra image contained a path absent from its canonical basis.
    #[error("an element path is absent from the finite canonical basis")]
    ElementOutsideCanonicalBasis,
    /// The algebra and element variants must remain paired.
    #[error("canonical-basis arithmetic returned another quiver-algebra variant")]
    UnexpectedElementVariant,
    /// The coordinate-swap name must resolve in the original product quiver.
    #[error("cannot find swapped product-quiver generator {name:?}")]
    MissingSwappedProductGenerator {
        /// Generated QPA product name.
        name: Box<str>,
    },
    /// QPA defines algebra Loewy length only for finite-dimensional input.
    #[error("the entered algebra is not finite dimensional")]
    InfiniteLoewyAlgebra,
    /// QPA's method requires an admissible quotient when input is not a path algebra.
    #[error("LoewyLength requires an admissible quotient of a path algebra")]
    NonadmissibleLoewyAlgebra,
    /// `SimpleModules` rejects infinite-dimensional input before constructing
    /// vertex representations.
    #[error("SimpleModules requires a finite-dimensional algebra")]
    InfiniteSimpleModuleAlgebra,
    /// The source requests projective-resolution degree `n + 1`.
    #[error("the global-dimension bound {bound} has no representable successor")]
    GlobalDimensionBoundOverflow {
        /// Entered nonnegative bound.
        bound: usize,
    },
    /// The complexity methods proceed only for admissible quotient algebras
    /// after their cached finite-global-dimension shortcut.
    #[error(
        "entered algebra is not finite dimensional or an admissible quotient of a path algebra"
    )]
    ComplexityRequiresAdmissibleQuotient,
    /// QPA's floating logarithmic expression is undefined for a nonzero
    /// projective-resolution object in these degrees.
    #[error("the complexity estimate is undefined at projective-resolution degree {degree}")]
    UndefinedComplexityEstimate {
        /// Entered arbitrary-precision complex degree.
        degree: Integer,
    },
    /// QPA's `Maximum` requires at least one vertex simple.
    #[error("ComplexityOfAlgebra found no simple modules")]
    ComplexityAlgebraHasNoSimpleModules,
    /// Both inputs to an add(M)-approximation must use the identical acting
    /// algebra.
    #[error("the two modules entered into the add(M)-approximation are not over the same algebra")]
    ApproximationModulesOverDifferentAlgebras,
    /// The common-algebra check makes QPA's internal Hom dispatch applicable.
    #[error("HomOverAlgebra unexpectedly rejected add(M)-approximation inputs")]
    ApproximationHomUnavailable,
    /// The number of radical composites must remain addressable.
    #[error("the add(M)-approximation radical-composite count overflows usize")]
    ApproximationMapCountOverflow,
    /// A nonempty quotient basis produces a nonempty direct sum in QPA.
    #[error("an internal nonempty add(M)-approximation direct sum returned fail")]
    EmptyApproximationDirectSum,
    /// Direct-sum structure maps correspond one-for-one with quotient lifts.
    #[error("add(M)-approximation structure-map and generator counts differ")]
    ApproximationMapCountMismatch,
    /// The source takes `Sum` only after a nonempty quotient-basis branch.
    #[error("an internal add(M)-approximation map sum was empty")]
    EmptyApproximationGeneratorList,
    /// A generated GAP resolution degree must fit Rust's addressable lazy
    /// cache in its direction from degree zero.
    #[error("resolution degree is not representable by usize: {degree}")]
    ResolutionDegreeNotRepresentable {
        /// Rejected arbitrary-precision degree.
        degree: Integer,
    },
    /// A panic while generating a lazy differential poisons the shared cache.
    #[error("the lazy projective-resolution cache was poisoned")]
    ProjectiveResolutionCachePoisoned,
    /// Every nonzero projective resolution starts with its degree-zero cover.
    #[error("the nonzero projective-resolution differential cache is empty")]
    EmptyProjectiveResolutionCache,
    /// A panic while generating a lazy injective differential poisons the
    /// shared cache.
    #[error("the lazy injective-resolution cache was poisoned")]
    InjectiveResolutionCachePoisoned,
    /// Every nonzero injective resolution starts with its degree-zero
    /// envelope.
    #[error("the nonzero injective-resolution differential cache is empty")]
    EmptyInjectiveResolutionCache,
    /// QPA's algebra dominant-dimension method requires a boolean
    /// selfinjectivity property after its finite-dimensionality check.
    #[error("selfinjectivity is unavailable for this finite dominant-dimension input")]
    DominantSelfinjectivityUnavailable,
    /// Every indecomposable projective of the checked finite algebra must
    /// remain in the module dominant-dimension method domain.
    #[error("a finite-algebra projective has no dominant-dimension method")]
    FiniteDominantModuleUnavailable,
    /// QPA rejects negative arguments to `NthSyzygy` before constructing the
    /// resolution.
    #[error("entered syzygy degree is negative: {degree}")]
    NegativeSyzygyDegree {
        /// Rejected arbitrary-precision degree.
        degree: Integer,
    },
    /// QPA's sum of indecomposable-projective socle dimension vectors must
    /// fit Rust's addressable module dimensions.
    #[error("the selfinjective socle multiplicity at vertex {vertex} overflows usize")]
    SelfinjectiveDimensionVectorOverflow {
        /// Zero-based quiver vertex position.
        vertex: usize,
    },
    /// QPA explicitly rejects infinite-dimensional input to its schurian
    /// property method.
    #[error("the entered algebra is not finite dimensional")]
    InfiniteSchurianAlgebra,
    /// A finite selfinjective quotient must construct its regular bimodule.
    #[error("the internal symmetric-algebra regular bimodule construction returned fail")]
    InternalSymmetricRegularModuleFailed,
    /// A finite selfinjective quotient must construct its dual regular
    /// bimodule.
    #[error("the internal symmetric-algebra dual bimodule construction returned fail")]
    InternalSymmetricDualModuleFailed,
    /// The ordered vertex-pair table used by QPA's semicommutativity loop
    /// must fit Rust's address space.
    #[error("the semicommutative vertex-pair count for {vertex_count} vertices overflows usize")]
    SemicommutativePairCountOverflow {
        /// Number of quiver vertices.
        vertex_count: usize,
    },
    /// A nonzero exact determinant must make the rational Cartan matrix
    /// invertible in QPA's following source step.
    #[error("a Cartan matrix with nonzero determinant unexpectedly had no rational inverse")]
    NonzeroCartanDeterminantNotInvertible,
    /// QPA's selfinjectivity call has no applicable method on this quotient.
    #[error("NakayamaPermutation cannot dispatch IsSelfinjectiveAlgebra on this algebra")]
    NakayamaPermutationSelfinjectivityUnavailable,
    /// The dedicated cache can contain only an unbound state, boolean false,
    /// or one list of index images.
    #[error("the NakayamaPermutation attribute cache contains an invalid fail state")]
    InvalidNakayamaPermutationCacheState,
    /// Concurrent attribute evaluation cannot bind both false and a value.
    #[error("the NakayamaPermutation attribute was concurrently cached as false")]
    NakayamaPermutationCacheConflict,
    /// The source closure accepts only modules of total dimension one.
    #[error("not a simple module was entered as an argument for the Nakayama permutation")]
    NakayamaPermutationModuleNotSimple,
    /// A foreign simple module can select beyond the captured projective list.
    #[error("Nakayama permutation module position {position} is outside 1..={vertex_count}")]
    NakayamaPermutationModulePositionOutOfRange {
        /// One-based position of the dimension-vector entry equal to one.
        position: usize,
        /// Number of captured projective socles.
        vertex_count: usize,
    },
    /// QPA's index closure accepts only the algebra's one-based simple indices.
    #[error(
        "incorrect Nakayama permutation index {index}; expected an integer in 1..={vertex_count}"
    )]
    IncorrectNakayamaPermutationIndex {
        /// Rejected arbitrary-precision integer.
        index: Integer,
        /// Number of simple-module indices.
        vertex_count: usize,
    },
    /// Selfinjectivity guarantees that every projective socle is a vertex
    /// simple and therefore has one dimension-vector entry equal to one.
    #[error("projective {projective} has no simple socle position in NakayamaPermutation")]
    SelfinjectiveProjectiveSocleHasNoSimplePosition {
        /// One-based indecomposable-projective position.
        projective: usize,
    },
    /// QPA's selfinjectivity call has no applicable method on this quotient.
    #[error("NakayamaAutomorphism cannot dispatch IsSelfinjectiveAlgebra on this algebra")]
    NakayamaAutomorphismSelfinjectivityUnavailable,
    /// The automorphism cache cannot contain GAP's `fail` sentinel.
    #[error("the NakayamaAutomorphism attribute cache contains an invalid fail state")]
    InvalidNakayamaAutomorphismCacheState,
    /// Concurrent evaluation cannot bind both boolean false and a matrix.
    #[error("the NakayamaAutomorphism attribute was concurrently cached as false")]
    NakayamaAutomorphismCacheConflict,
    /// The source closure requires an element in its captured algebra.
    #[error("the entered argument is not in the algebra")]
    NakayamaAutomorphismElementOutsideAlgebra,
    /// A successful finite selfinjective computation requires a canonical basis.
    #[error("NakayamaAutomorphism requires a finite canonical algebra basis")]
    MissingNakayamaCanonicalBasis,
    /// Every selfinjective indecomposable projective has a nonzero socle generator.
    #[error("projective {projective} has no socle generator in NakayamaAutomorphism")]
    EmptyProjectiveSocleGenerator {
        /// One-based indecomposable-projective position.
        projective: usize,
    },
    /// GAP's Steinitz primitive received inconsistent basis shapes or ranks.
    #[error("BaseSteinitzVectors received invalid Nakayama basis data")]
    InvalidNakayamaSteinitzInput,
    /// The source Steinitz extension must produce a full algebra basis.
    #[error("BaseSteinitzVectors did not complete the Nakayama socle basis")]
    IncompleteNakayamaSteinitzBasis,
    /// Addressable Steinitz column positions cannot overflow.
    #[error("the Nakayama Steinitz column counter overflowed usize")]
    NakayamaSteinitzColumnOverflow,
    /// Addressable Steinitz row positions cannot overflow.
    #[error("the Nakayama Steinitz row counter overflowed usize")]
    NakayamaSteinitzRowOverflow,
    /// A selfinjective algebra's Frobenius pairing must be nondegenerate.
    #[error("the Frobenius pairing constructed by NakayamaAutomorphism is singular")]
    SingularNakayamaPairing,
    /// The cached closure matrix must match the canonical-basis dimension.
    #[error("the cached Nakayama matrix has shape {shape:?}, expected {dimension} x {dimension}")]
    InvalidNakayamaAutomorphismMatrixShape {
        /// Cached matrix shape.
        shape: (usize, usize),
        /// Current canonical-basis dimension.
        dimension: usize,
    },
    /// Internal canonical-basis reconstruction requires equal list lengths.
    #[error(
        "Nakayama linear combination has {coefficients} coefficients for {basis} basis vectors"
    )]
    NakayamaLinearCombinationLength {
        /// Canonical-basis length.
        basis: usize,
        /// Coordinate-list length.
        coefficients: usize,
    },
    /// The Frobenius-form cache cannot contain GAP's `fail` sentinel.
    #[error("the FrobeniusForm attribute cache contains an invalid fail state")]
    InvalidFrobeniusFormCacheState,
    /// Successful Nakayama construction must install the form as a side effect.
    #[error("NakayamaAutomorphism did not install FrobeniusForm")]
    MissingInstalledFrobeniusForm,
    /// Concurrent evaluation cannot bind both false and a form matrix.
    #[error("the FrobeniusForm attribute was concurrently cached as false")]
    FrobeniusFormCacheConflict,
    /// The form closure accepts only elements of its captured quotient.
    #[error("a FrobeniusForm argument is not in the algebra")]
    FrobeniusFormElementOutsideAlgebra,
    /// The cached form matrix must match the canonical-basis dimension.
    #[error("the cached Frobenius form has shape {shape:?}, expected {dimension} x {dimension}")]
    InvalidFrobeniusFormMatrixShape {
        /// Cached matrix shape.
        shape: (usize, usize),
        /// Current canonical-basis dimension.
        dimension: usize,
    },
    /// The Frobenius-functional cache cannot contain GAP's `fail` sentinel.
    #[error("the FrobeniusLinearFunctional cache contains an invalid fail state")]
    InvalidFrobeniusFunctionalCacheState,
    /// Successful Nakayama construction must install the functional.
    #[error("NakayamaAutomorphism did not install FrobeniusLinearFunctional")]
    MissingInstalledFrobeniusFunctional,
    /// Concurrent evaluation cannot bind both false and functional coordinates.
    #[error("the FrobeniusLinearFunctional attribute was concurrently cached as false")]
    FrobeniusFunctionalCacheConflict,
    /// The functional closure accepts only elements of its captured quotient.
    #[error("the FrobeniusLinearFunctional argument is not in the algebra")]
    FrobeniusFunctionalElementOutsideAlgebra,
    /// Cached functional coordinates must match the canonical basis.
    #[error("the cached Frobenius functional has length {actual}, expected {dimension}")]
    InvalidFrobeniusFunctionalLength {
        /// Cached coordinate count.
        actual: usize,
        /// Current canonical-basis dimension.
        dimension: usize,
    },
    /// The order cache cannot contain GAP's `fail` sentinel.
    #[error("the OrderOfNakayamaAutomorphism cache contains an invalid fail state")]
    InvalidNakayamaOrderCacheState,
    /// Path-algebra order dispatch requires a selfinjectivity method.
    #[error("OrderOfNakayamaAutomorphism cannot dispatch IsSelfinjectiveAlgebra")]
    NakayamaOrderSelfinjectivityUnavailable,
    /// Concurrent order evaluation cannot bind both false and a value.
    #[error("OrderOfNakayamaAutomorphism was concurrently cached as false")]
    NakayamaOrderCacheConflict,
    /// QPA's quotient order method attempts to call a false automorphism.
    #[error("OrderOfNakayamaAutomorphism attempted to call boolean false")]
    CalledFalseNakayamaAutomorphism,
    /// GAP's matrix order operation requires a nonempty square matrix here.
    #[error("invalid Nakayama automorphism matrix shape {shape:?} for Order")]
    InvalidNakayamaOrderMatrixShape {
        /// Rejected matrix shape.
        shape: (usize, usize),
    },
    /// The proven cyclotomic-index bound must fit addressable iteration.
    #[error("the rational finite-order bound for dimension {dimension} overflows usize")]
    NakayamaOrderCyclotomicBoundOverflow {
        /// Matrix dimension.
        dimension: usize,
    },
    /// The supported prime characteristic must fit addressable arithmetic.
    #[error("prime characteristic {prime} does not fit usize")]
    NakayamaOrderPrimeTooLarge {
        /// Prime characteristic.
        prime: u64,
    },
    /// The unipotent exponent calculation must fit addressable dimensions.
    #[error("the GL({dimension},{prime}) unipotent exponent bound overflows usize")]
    NakayamaOrderPrimePowerOverflow {
        /// Matrix dimension.
        dimension: usize,
        /// Prime characteristic.
        prime: u64,
    },
    /// An annihilator basis requires QPA's finite canonical algebra basis.
    #[error("AnnihilatorOfModule requires a finite-dimensional algebra")]
    InfiniteAnnihilatorAlgebra,
    /// QPA's finite canonical path basis starts with all vertex idempotents.
    #[error("the finite canonical algebra basis does not start with the quiver vertices")]
    MissingCanonicalVertexBasis,
    /// Positive summand counts make both random projective lists nonempty.
    #[error("an internal RandomModule projective direct sum returned fail")]
    InternalRandomModuleDirectSumFailed,
    /// Both random projective sums act over the identical entered algebra.
    #[error("HomOverAlgebra rejected internal RandomModule projective sums")]
    InternalRandomModuleHomRejected,
    /// `TensorProductOfModules(M, N)` requires `N` over the opposite algebra.
    #[error("the entered modules are not over the appropriate algebras")]
    TensorProductModulesOverWrongAlgebras,
    /// Tensor-presentation dimensions must fit Rust's addressable size.
    #[error("a TensorProductOfModules dimension pairing overflows usize")]
    TensorProductDimensionOverflow,
    /// Every internal indecomposable projective has one distinguished generator.
    #[error("an internal TensorProductOfModules projective has no generator")]
    EmptyTensorProjectiveGenerator,
    /// Distinguished internal projective generators are vertex-uniform.
    #[error("an internal TensorProductOfModules projective generator is not uniform")]
    NonuniformTensorProjectiveGenerator,
    /// The source applies quotient-only nested element extraction to `P1`
    /// generators even when the acting algebra is a path algebra.
    #[error(
        "QPA's nontrivial path-algebra TensorProductOfModules branch cannot extract P1 vertices"
    )]
    TensorProductPathAlgebraP1Extraction,
    /// Both internal projective terms are nonempty in this source branch.
    #[error("an internal TensorProductOfModules projective presentation is empty")]
    EmptyTensorProjectivePresentation,
    /// A nonzero `P1 tensor N` supplies at least one physical relation row.
    #[error("the internal TensorProductOfModules relation matrix is empty")]
    EmptyTensorRelationMatrix,
    /// GAP's Steinitz extension must complete the tensor relation basis.
    #[error(
        "NaturalHomomorphismBySubspace did not produce a complete tensor-product quotient basis"
    )]
    IncompleteTensorProductQuotientBasis,
    /// QPA's `P1 tensor N = 0` closure references a never-assigned local.
    #[error("QPA TensorProductOfModules elementary closure has unassigned local size_matrix")]
    TensorProductUnboundSizeMatrix,
    /// QPA overwrites the closure's captured projection count with the last
    /// `P1` vertex and can consequently index `projectionsP0` out of range.
    #[error(
        "QPA TensorProductOfModules closure requests projection {index}, but only {count} exist"
    )]
    TensorProductClosureProjectionOutOfRange {
        /// One-based projection requested by QPA.
        index: usize,
        /// Number of projective-cover projections.
        count: usize,
    },
    /// The same captured-local bug can make QPA build too few acted summands
    /// for its subsequent `size_matrix[2]` loop.
    #[error(
        "QPA TensorProductOfModules closure requests acted summand {index}, but only {count} exist"
    )]
    TensorProductClosureActedIndexOutOfRange {
        /// One-based acted summand requested by QPA.
        index: usize,
        /// Number of acted summands constructed by the overwritten count.
        count: usize,
    },
    /// The projective cover is surjective, so every entered left element has
    /// a preimage when it belongs to the captured module.
    #[error("a tensor-product left element unexpectedly has no projective-cover preimage")]
    MissingTensorProjectiveCoverPreimage,
    /// `ExtOverAlgebra` requires one common acting algebra.
    #[error("the two modules entered are not modules over the same algebra")]
    ExtModulesOverDifferentAlgebras,
    /// A Hom call cannot reject the already checked Ext modules.
    #[error("HomOverAlgebra rejected modules in the ExtOverAlgebra construction")]
    InternalExtHomRejected,
    /// GAP's Steinitz extension must complete the image basis to the Hom basis.
    #[error("NaturalHomomorphismBySubspace did not produce a complete Ext quotient basis")]
    IncompleteExtQuotientBasis,
    /// The coefficient closure only accepts raw vectors over its original field.
    #[error("the Ext representative map uses a different coefficient field")]
    ExtCoefficientFieldMismatch,
    /// The coefficient closure accepts only maps in `Hom(Omega(M), N)`.
    #[error("the entered map is outside the Hom space used by ExtOverAlgebra")]
    ExtMapOutsideHomSpace,
    /// Internal star Hom calls use modules over the already checked algebra.
    #[error("HomOverAlgebra rejected an internal StarOfModule Hom space")]
    InternalStarHomRejected,
    /// Every indecomposable projective has its distinguished vertex generator.
    #[error("an indecomposable projective has no generator in StarOfModule")]
    EmptyStarProjectiveGenerator,
    /// Arrow multiplication preserves one-vertex support.
    #[error("left multiplication by an arrow did not define a projective map")]
    StarArrowMapNotUniform,
    /// A nonzero star vertex has a raw physical Hom basis.
    #[error("a nonzero StarOfModule vertex has no coordinate basis")]
    MissingStarCoordinateBasis,
    /// Composition with left multiplication remains in the corresponding Hom space.
    #[error("a StarOfModule composite lies outside its Hom basis")]
    StarCompositionOutsideHomBasis,
    /// Both Hom-space families used by the evaluation construction share one algebra.
    #[error("HomOverAlgebra rejected an internal double-star Hom space")]
    InternalDoubleStarHomRejected,
    /// A nonzero double-star vertex has a nonempty raw physical Hom basis.
    #[error("double-star vertex {vertex} has no Hom coordinate basis")]
    MissingDoubleStarHomCoordinateBasis {
        /// One-based vertex position.
        vertex: usize,
    },
    /// Evaluation at a module basis vector is an element of the computed double dual.
    #[error("the evaluation map at vertex {vertex} lies outside the double-star Hom basis")]
    DoubleStarEvaluationOutsideHomBasis {
        /// One-based vertex position.
        vertex: usize,
    },
    /// Finite projective bases require the opposite algebra's canonical basis.
    #[error("the opposite algebra has no canonical basis for double-star evaluation")]
    MissingDoubleStarCanonicalBasis,
    /// Every projective-basis element occurs in the opposite canonical basis.
    #[error(
        "opposite projective {projective}, vertex {vertex}, basis element {position} is absent from the canonical basis"
    )]
    MissingDoubleStarProjectiveBasisElement {
        /// One-based indecomposable-projective position.
        projective: usize,
        /// One-based representation vertex.
        vertex: usize,
        /// One-based position inside that vertex basis.
        position: usize,
    },
    /// Source-generated star actions satisfy the opposite relations.
    #[error("the generated StarOfModule does not satisfy the opposite-algebra relations")]
    GeneratedStarViolatesRelations,
    /// The explicit transpose construction requires QPA's finite canonical basis.
    #[error("TransposeOfModule requires a finite-dimensional algebra")]
    InfiniteTransposeAlgebra,
    /// Every minimal module generator is supported at a quiver vertex.
    #[error("a minimal generator has no vertex support in TransposeOfModule")]
    MissingTransposeGeneratorSupport,
    /// The source's internal projective sums are nonempty in this branch.
    #[error("an internal TransposeOfModule projective direct sum returned fail")]
    InternalTransposeDirectSumFailed,
    /// Projective covers created by the transpose-on-maps construction retain
    /// the direct-sum metadata required by QPA's lifting method.
    #[error("an internal TransposeOfModuleHomomorphism projective lifting returned fail")]
    InternalTransposeHomomorphismLiftFailed,
    /// The lifted projective-cover square commutes by construction.
    #[error("the internal TransposeOfModuleHomomorphism kernel square does not commute")]
    InternalTransposeHomomorphismKernelSquareFailed,
    /// Applying the contravariant star functor preserves the presentation
    /// square used to induce the transpose map.
    #[error("the internal TransposeOfModuleHomomorphism cokernel square does not commute")]
    InternalTransposeHomomorphismCokernelSquareFailed,
    /// Minimal syzygy generators must have exactly one right-support vertex.
    #[error("a minimal syzygy generator is not uniform in TransposeOfModule")]
    NonuniformTransposeSyzygyGenerator,
    /// Every opposite projective summand has its distinguished generator.
    #[error("an opposite projective has no generator in TransposeOfModule")]
    EmptyTransposeProjectiveGenerator,
    /// `DTr(M)` remains over the identical acting algebra of `M`.
    #[error("HomOverAlgebra rejected the internal tau-rigidity Hom space")]
    InternalTauRigidHomRejected,
    /// `SimpleModules` accepts quotients only after QPA's admissibility filter
    /// has been established.
    #[error("SimpleModules requires an admissible quotient of a path algebra")]
    NonadmissibleSimpleModuleAlgebra,
    /// A vertex-simple representation must satisfy every admissible relation.
    #[error("a generated vertex-simple module does not satisfy the algebra relations")]
    GeneratedSimpleViolatesRelations,
    /// A projective action generated from canonical basis multiplication must satisfy relations.
    #[error("a generated indecomposable projective does not satisfy the algebra relations")]
    GeneratedProjectiveViolatesRelations,
    /// Every projective right ideal contains its generating vertex idempotent.
    #[error("the canonical basis of projective {vertex} has no vertex generator")]
    MissingProjectiveVertexGenerator {
        /// One-based projective/vertex position.
        vertex: usize,
    },
    /// QPA's projectivity dimension sum exceeded Rust's addressable size.
    #[error("the projective dimension test overflows usize")]
    ProjectiveDimensionOverflow,
    /// `HomFromProjective` requires an element of the identical target module.
    #[error("the projective generator element is not in the target module")]
    ProjectiveElementInWrongModule,
    /// A projective basis cannot be extracted from infinite-dimensional input.
    #[error("BasisOfProjectives returned fail")]
    ProjectiveBasisUnavailable,
    /// Every nonzero module has at least one top/minimal generator.
    #[error("a nonzero projective cover has no summands")]
    EmptyNonzeroProjectiveCover,
    /// QPA's minimal generators are vertex-uniform.
    #[error("a minimal generator is not supported at exactly one vertex")]
    MinimalGeneratorNotVertexUniform,
    /// Projective lifting validates the source before its other inputs.
    #[error("the source of the second map is not projective")]
    LiftingSourceNotProjective,
    /// Projective lifting requires an epimorphism as its first map.
    #[error("the first map is not surjective")]
    LiftingFirstMapNotSurjective,
    /// Both maps in a projective lifting diagram must end at the same module.
    #[error("the entered maps do not have a common range")]
    LiftingMapsDoNotHaveCommonRange,
    /// Every nonzero direct projective summand has its installed generator.
    #[error("a nonzero projective summand has no minimal generator")]
    EmptyProjectiveSummandGenerators,
    /// Surjectivity guarantees a lift of every projective-generator image.
    #[error("a projective-generator image unexpectedly has no preimage")]
    MissingProjectiveLiftPreimage,
    /// A projective lift preimage must be vertex-uniform.
    #[error("a projective lifting preimage is not vertex-uniform")]
    ProjectiveLiftNotVertexUniform,
    /// Factoring-through computations require a common acting algebra.
    #[error("the factoring modules are not modules over the same algebra")]
    FactoringModulesOverDifferentAlgebras,
    /// QPA explicitly rejects the zero module for its stable endomorphism
    /// quotient operation.
    #[error("zero module given as input for EndModuleProjOverAlgebra")]
    ZeroStableEndomorphismModule,
    /// A Hom call cannot reject the checked projective-cover algebra.
    #[error("HomOverAlgebra rejected the projective cover's common algebra")]
    InternalProjectiveHomRejected,
    /// Flattened physical map data exceeded Rust's addressable size.
    #[error("the physical module-map size overflows usize")]
    PhysicalMapSizeOverflow,
    /// A raw physical map vector must fill every vertex block exactly.
    #[error("the physical module-map vector has length {actual}, expected {expected}")]
    PhysicalMapLength {
        /// Entered flat-vector length.
        actual: usize,
        /// Required total physical block length.
        expected: usize,
    },
    /// Projective positions are one-based quiver-vertex positions.
    #[error("projective vertex {vertex} is outside 1..={vertex_count}")]
    ProjectiveVertexOutOfRange {
        /// Entered vertex position.
        vertex: usize,
        /// Number of vertices/projectives.
        vertex_count: usize,
    },
    /// QPA flattens physical module coordinates against a basis with a zero
    /// placeholder at every zero-dimensional vertex space.
    #[error("projective element has {actual} physical coordinates, expected {expected}")]
    ProjectiveElementPhysicalCoordinateCount {
        /// Entered flattened physical coordinate count.
        actual: usize,
        /// Flattened projective-basis count.
        expected: usize,
    },
    /// The two source methods require an acyclic path algebra or an admissible
    /// quotient, respectively.
    #[error("the entered algebra is not finite dimensional in QPA's projective-element domain")]
    ElementInIndecProjectiveAlgebraNotFinite,
    /// The algebra-element-to-projective conversion requires the identical
    /// path algebra or quotient family.
    #[error("the projective algebra element belongs to another algebra")]
    ForeignProjectiveAlgebraElement,
    /// QPA rejects zero in `ElementIn_vA_AsElementInIndecProj`.
    #[error("the entered projective algebra element is zero")]
    ZeroProjectiveAlgebraElement,
    /// A nonzero element must have exactly one left-support vertex.
    #[error("the projective algebra element is not left uniform")]
    ProjectiveAlgebraElementNotLeftUniform,
    /// QPA's coefficient-bearing successive tips must occur in the selected
    /// projective-basis component.
    #[error(
        "a coefficient-bearing tip is absent from projective basis component ({left_vertex}, {right_vertex})"
    )]
    ProjectiveLeadingTermMissing {
        /// Generating/left-support vertex.
        left_vertex: usize,
        /// Right-support component vertex.
        right_vertex: usize,
    },
    /// A distinguished generator of an indecomposable projective must have
    /// exactly one support vertex.
    #[error("a projective summand generator is not vertex uniform")]
    ProjectiveGeneratorNotVertexUniform,
    /// GAP's `IsMatrix` method domain excludes an empty element matrix here.
    #[error("the projective element matrix must be nonempty")]
    EmptyProjectiveElementMatrix,
    /// GAP matrices are rectangular homogeneous lists.
    #[error("the projective element matrix is ragged")]
    RaggedProjectiveElementMatrix,
    /// Matrix rows/columns and the two entered support lists must agree.
    #[error(
        "projective matrix shape {rows}x{columns} disagrees with {range_vertices} range and {source_vertices} source vertices"
    )]
    ProjectiveMatrixVertexListDimensions {
        /// Matrix row count.
        rows: usize,
        /// Matrix column count.
        columns: usize,
        /// Entered range-vertex count.
        range_vertices: usize,
        /// Entered source-vertex count.
        source_vertices: usize,
    },
    /// Every matrix entry must belong to the identical entered algebra.
    #[error("a projective matrix entry belongs to another algebra")]
    ForeignProjectiveMatrixElement,
    /// Every nonzero matrix entry must be both left and right uniform.
    #[error("a projective matrix entry is not uniform")]
    ProjectiveMatrixElementNotUniform,
    /// All nonzero entries in a row must have one left support.
    #[error("projective matrix row {row} is not homogeneous")]
    ProjectiveMatrixRowNotHomogeneous {
        /// One-based row position.
        row: usize,
    },
    /// All nonzero entries in a column must have one right support.
    #[error("projective matrix column {column} is not homogeneous")]
    ProjectiveMatrixColumnNotHomogeneous {
        /// One-based column position.
        column: usize,
    },
    /// Detected supports must equal the entered lists, including their order.
    #[error("the entered vertex lists disagree with projective matrix support")]
    ProjectiveMatrixSupportDiscrepancy,
    /// Nonempty source/range vertex lists make these internal direct sums
    /// infallible in QPA's method domain.
    #[error("an internal nonempty projective direct sum returned fail")]
    InternalProjectiveDirectSumFailed,
    /// A nonzero uniform matrix entry must produce a projective homomorphism.
    #[error("HomFromProjective rejected an internally uniform matrix entry")]
    InternalProjectiveMatrixHomRejected,
    /// The configured QPA/GBNP finite-basis path failed.
    #[error(transparent)]
    Backend(#[from] PathalgGbnpError),
    /// Enveloping-algebra construction failed.
    #[error(transparent)]
    Tensor(#[from] TensorProductError),
    /// Product-quiver projection failed.
    #[error(transparent)]
    Product(#[from] QuiverProductError),
    /// Opposite-path reconstruction failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// Arithmetic in `OriginalPathAlgebra(A)` failed.
    #[error(transparent)]
    ParentPathAlgebra(#[from] PathAlgebraError),
    /// Path inclusion in the source algebra failed.
    #[error(transparent)]
    PathTensor(#[from] PathAlgebraTensorError),
    /// Source-algebra multiplication failed.
    #[error(transparent)]
    Element(#[from] QuiverAlgebraElementError),
    /// Quotient canonical-basis coordinates failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Exact matrix construction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// The checked QPA right-module constructor rejected its generated data.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// Generic QPA module dualization failed.
    #[error(transparent)]
    Functor(#[from] ModuleFunctorError),
}

/// Failure in QPA's generic module-dual functor.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ModuleFunctorError {
    /// The generated transposed actions should satisfy the opposite relations.
    #[error("the generated dual module does not satisfy the opposite-algebra relations")]
    GeneratedDualViolatesRelations,
    /// Parent path-algebra opposite construction failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Quotient opposite construction through QPA's GBNP backend failed.
    #[error(transparent)]
    Backend(#[from] PathalgGbnpError),
    /// The checked module constructor or reciprocal cache failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, Rational, RationalField},
    };
    use qpa_core::path_module::right_fac_approximation_via_all_basis_images;
    use qpa_core::{
        AlgebraGlobalDimension, ArrowSpec, CachedAttribute, DenseMatrix, ModuleDimensionSeries,
        NamedArrowMatrix, NamedArrowSpec, PathAlgebra, PathAlgebraModule, PathAlgebraModuleElement,
        PathAlgebraModuleHomomorphism, PathAlgebraModuleHomomorphismPair, PathAlgebraModuleMatrix,
        QpaCoefficientField, Quiver, QuiverAlgebra, QuiverAlgebraElement,
        StructureConstantsAlgebra, VertexIndex, dimension_vector_partial_order,
        direct_sum_of_qpa_modules, end_over_algebra, from_end_m_to_hom_mm, from_hom_mm_to_end_m,
        identity_module_mapping, include_in_path_algebra, is_indecomposable_module,
        is_semisimple_module, is_simple_qpa_module, kernel_inclusion, lifting_inclusion_morphisms,
        minimal_generating_set_of_module, radical_series, right_fac_approximation,
        right_module_hom_over_algebra, right_module_over_path_algebra, socle_of_module,
        socle_of_module_inclusion, top_of_module, top_of_module_projection, trace_of_module,
        underlying_linear_map, zero_module, zero_module_mapping,
    };
    use rand::{RngCore, SeedableRng, rngs::StdRng};
    use serde_json::{Value, json};

    use super::{
        AlmostSplitSide, DominantDimensionResult, EnvelopingModuleError, ExtOverAlgebraResult,
        GlobalDimensionOfAlgebraResult, ModuleLengthError, NakayamaAutomorphismOrder,
        PredecessorMultiplicity, PredecessorsOfModuleResult, TensorProductOfModulesResult,
        algebra_as_module_of_enveloping_algebra, algebra_as_module_over_enveloping_algebra,
        all_indec_modules_of_length_at_most, all_modules_of_length_at_most,
        all_modules_of_length_plus_one, almost_split_sequence, almost_split_sequence_from_side,
        annihilator_of_module, basis_of_projectives, canonical_basis_coefficients, cartan_matrix,
        complexity_of_algebra, complexity_of_module, coxeter_matrix, coxeter_polynomial, d_tr,
        d_tr_homomorphism, d_tr_homomorphism_iterated, d_tr_iterated,
        dominant_dimension_of_algebra, dominant_dimension_of_module,
        dual_of_algebra_as_module_over_enveloping_algebra, dual_of_module,
        dual_of_module_homomorphism, element_in_indec_projective,
        element_in_v_a_as_element_in_indec_proj, end_modulo_projectives, ext_over_algebra,
        finite_canonical_basis, finite_qpa_default_basis, first_syzygy, frobenius_form,
        frobenius_linear_functional, from_identity_to_double_star_homomorphism,
        from_matrix_to_homomorphism_of_projectives, global_dimension_from_known_properties,
        global_dimension_of_algebra, gorenstein_dimension, gorenstein_dimension_of_algebra,
        hom_factoring_through_projectives, hom_from_projective, indec_injective_modules,
        indec_projective_modules, inj_dimension, inj_dimension_of_module, injective_envelope,
        irreducible_morphisms_ending_in, irreducible_morphisms_starting_in, is_exceptional_module,
        is_finite_global_dimension_from_known_properties, is_hereditary_algebra,
        is_injective_module, is_omega_periodic, is_omega_periodic_with_progress,
        is_projective_module, is_rigid_module, is_schurian_algebra, is_selfinjective_algebra,
        is_semicommutative_algebra, is_symmetric_algebra, is_tau_periodic, is_tau_rigid_module,
        is_weakly_symmetric_algebra, left_approximation_by_add_m, lifting_morphism_from_projective,
        matrix_of_homomorphism_between_projectives, minimal_left_approximation,
        minimal_right_approximation, nakayama_automorphism, nakayama_functor_of_module,
        nakayama_functor_of_module_homomorphism, nakayama_permutation,
        opposite_nakayama_functor_of_module, opposite_nakayama_functor_of_module_homomorphism,
        order_of_nakayama_automorphism, predecessors_of_module, proj_dimension,
        proj_dimension_of_module, projective_cover, random_module, right_approximation_by_add_m,
        simple_modules, socle_series, star_of_module, star_of_module_homomorphism,
        tensor_product_of_modules, tr_d, tr_d_homomorphism, tr_d_homomorphism_iterated,
        tr_d_iterated, transpose_of_module, transpose_of_module_homomorphism,
    };
    use crate::{admissible_quotient_is_semisimple, enveloping_algebra, quotient_by_relations};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn q(value: i64) -> Rational {
        Q.nth(Integer::from(value))
    }

    fn global_dimension_json(result: Option<GlobalDimensionOfAlgebraResult>) -> Value {
        match result {
            Some(GlobalDimensionOfAlgebraResult::Finite(dimension)) => dimension.into(),
            Some(GlobalDimensionOfAlgebraResult::Infinite) => "infinity".into(),
            Some(GlobalDimensionOfAlgebraResult::ExceedsBound) => false.into(),
            None => Value::Null,
        }
    }

    fn cached_global_dimension_json(result: Option<AlgebraGlobalDimension>) -> Value {
        match result {
            Some(AlgebraGlobalDimension::Finite(dimension)) => dimension.into(),
            Some(AlgebraGlobalDimension::Infinite) => "infinity".into(),
            None => Value::Null,
        }
    }

    fn bounded_dimension_json(result: Option<usize>) -> Value {
        result.map_or(Value::Bool(false), Into::into)
    }

    fn dominant_dimension_json(result: Option<DominantDimensionResult>) -> Value {
        match result {
            Some(DominantDimensionResult::Finite(dimension)) => dimension.into(),
            Some(DominantDimensionResult::Infinite) => "infinity".into(),
            Some(DominantDimensionResult::ExceedsBound) => false.into(),
            None => Value::Null,
        }
    }

    fn fresh_a2_module<F: Field>(field: &F, dimension_vector: Vec<usize>) -> PathAlgebraModule<F> {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        right_module_over_path_algebra(algebra, dimension_vector, &[])
            .unwrap()
            .unwrap()
    }

    fn fresh_a2_zero<F: Field>(field: &F) -> PathAlgebraModule<F> {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "z")]).unwrap(),
            )
            .unwrap(),
        );
        zero_module(&algebra).unwrap()
    }

    fn fresh_dual_simple<F: Field>(field: &F) -> PathAlgebraModule<F> {
        let algebra = fresh_dual_numbers(field);
        right_module_over_path_algebra(
            algebra,
            vec![1],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(vec![vec![field.zero()]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap()
    }

    fn fresh_dual_numbers<F: Field>(field: &F) -> QuiverAlgebra<F> {
        let parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = parent.generators()[1].clone();
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.pow(2).unwrap()]).unwrap())
    }

    fn fresh_radical_square_zero_a3<F: Field>(field: &F) -> QuiverAlgebra<F> {
        let parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(
                3,
                vec![
                    ArrowSpec::named(vertex(1), vertex(2), "x"),
                    ArrowSpec::named(vertex(2), vertex(3), "y"),
                ],
            )
            .unwrap(),
        )
        .unwrap();
        let arrows = parent.quiver().arrows().collect::<Vec<_>>();
        let x = parent.element(arrows[0].path()).unwrap();
        let y = parent.element(arrows[1].path()).unwrap();
        QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x.multiply(&y).unwrap()]).unwrap())
    }

    fn isolated_path_algebra<F: Field>(field: F, vertex_count: usize) -> QuiverAlgebra<F> {
        let vertices = (1..=vertex_count)
            .map(|index| format!("v{index}").into_boxed_str())
            .collect();
        QuiverAlgebra::Path(
            PathAlgebra::new(field, Quiver::from_named(vertices, Vec::new()).unwrap()).unwrap(),
        )
    }

    struct ZeroRng;

    impl RngCore for ZeroRng {
        fn next_u32(&mut self) -> u32 {
            0
        }

        fn next_u64(&mut self) -> u64 {
            0
        }

        fn fill_bytes(&mut self, destination: &mut [u8]) {
            destination.fill(0);
        }
    }

    #[test]
    fn regular_bimodule_for_dual_numbers_has_qpa_row_actions() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(Q, quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );

        let module = algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        let again = algebra_as_module_of_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        let enveloping = enveloping_algebra(&algebra).unwrap();

        assert_eq!(module.dimension_vector(), &[2]);
        assert_eq!(module.arrow_matrices().len(), 2);
        assert_eq!(
            module.arrow_matrices()[0].as_flat(),
            &[q(0), q(1), q(0), q(0)]
        );
        assert_eq!(
            module.arrow_matrices()[1].as_flat(),
            &[q(0), q(1), q(0), q(0)]
        );
        assert!(module.is_identical(&again));
        assert!(
            module
                .retained_tensor_product()
                .unwrap()
                .is_identical(&enveloping)
        );

        let generic_dual = dual_of_module(&module).unwrap();
        let generic_back = dual_of_module(&generic_dual).unwrap();
        let bimodule_dual = dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        assert_eq!(
            generic_dual.arrow_matrices()[0].as_flat(),
            &[q(0), q(0), q(1), q(0)]
        );
        assert!(generic_back.is_identical(&module));
        assert_eq!(
            bimodule_dual.arrow_matrices()[0].as_flat(),
            &[q(0), q(0), q(1), q(0)]
        );
    }

    #[test]
    fn dual_of_module_homomorphism_transposes_reverses_and_caches_like_qpa() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra,
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![q(1)], vec![q(0)]], Q).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let map = right_module_hom_over_algebra(
            module.clone(),
            module,
            vec![
                DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
                DenseMatrix::identity(1, Q),
            ],
        )
        .unwrap();

        let dual = dual_of_module_homomorphism(&map).unwrap();
        assert_eq!(dual.source().dimension_vector(), &[2, 1]);
        assert_eq!(dual.range().dimension_vector(), &[2, 1]);
        assert_eq!(
            dual.vertex_matrices(),
            &[
                DenseMatrix::from_rows(vec![vec![q(1), q(0)], vec![q(0), q(0)]], Q).unwrap(),
                DenseMatrix::identity(1, Q),
            ]
        );
        assert!(dual.is_identical(&dual_of_module_homomorphism(&map).unwrap()));
        let double_dual = dual_of_module_homomorphism(&dual).unwrap();
        assert_eq!(double_dual, map);
        assert!(
            !double_dual.is_identical(&map),
            "the QPA source does not install a reciprocal map attribute"
        );

        let e = |value: i64| Z2.nth(Integer::from(value));
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra,
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![e(1)], vec![e(0)]], Z2).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let map = right_module_hom_over_algebra(
            module.clone(),
            module,
            vec![
                DenseMatrix::from_rows(vec![vec![e(1), e(0)], vec![e(0), e(0)]], Z2).unwrap(),
                DenseMatrix::identity(1, Z2),
            ],
        )
        .unwrap();
        assert_eq!(
            dual_of_module_homomorphism(&map).unwrap().vertex_matrices(),
            map.vertex_matrices()
        );
    }

    #[test]
    fn regular_bimodule_supports_characteristic_two_and_zero_spaces() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let module = algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();

        assert_eq!(module.dimension_vector(), &[1, 1, 0, 1]);
        assert_eq!(module.dimension(), 3);
        assert_eq!(module.arrow_matrices().len(), 4);
        assert!(module.arrow_matrices().iter().all(|matrix| {
            matrix
                .as_flat()
                .iter()
                .all(|entry| Z2.to_integer(entry) == 0 || Z2.to_integer(entry) == 1)
        }));

        let dual = dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        assert_eq!(dual.dimension_vector(), &[1, 0, 1, 1]);
    }

    #[test]
    fn cyclic_path_algebra_returns_qpa_fail_before_enveloping_construction() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());

        assert!(
            algebra_as_module_over_enveloping_algebra(&algebra)
                .unwrap()
                .is_none()
        );
        assert!(algebra.cached_enveloping_algebra().is_none());
        assert!(
            dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
                .unwrap()
                .is_none()
        );
        assert!(algebra.has_dual_enveloping_module_attribute());
    }

    fn matrix_json<F: Field>(
        matrix: &DenseMatrix<F>,
        coefficient: impl Fn(&F::Element) -> String,
    ) -> Value {
        json!({
            "shape": [matrix.nrows(), matrix.ncols()],
            "rows": (0..matrix.nrows())
                .map(|row| matrix.row(row).iter().map(&coefficient).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
        })
    }

    fn module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "dimension_vector": module.dimension_vector(),
            "dimension": module.dimension(),
            "arrows": module
                .right_acting_algebra()
                .quiver()
                .arrows()
                .map(|arrow| arrow.name().to_owned())
                .collect::<Vec<_>>(),
            "matrices": module
                .arrow_matrices()
                .iter()
                .map(|matrix| matrix_json(matrix, coefficient))
                .collect::<Vec<_>>(),
        })
    }

    fn length_module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!({
            "dimension": module.dimension_vector(),
            "matrices": module.matrices_of_path_algebra_module().iter().map(|matrix| {
                match matrix {
                    PathAlgebraModuleMatrix::ZeroScalar => json!([]),
                    PathAlgebraModuleMatrix::Matrix(matrix) => json!(
                        (0..matrix.nrows()).map(|row| matrix.row(row).iter()
                            .map(coefficient).collect::<Vec<_>>()).collect::<Vec<_>>()
                    ),
                }
            }).collect::<Vec<_>>(),
        })
    }

    fn length_module_list_json<F: Field>(
        modules: &[PathAlgebraModule<F>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!(
            modules
                .iter()
                .map(|module| length_module_json(module, coefficient))
                .collect::<Vec<_>>()
        )
    }

    #[allow(clippy::needless_pass_by_value)]
    fn modulehom_length_case<F: QpaCoefficientField>(
        field: F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let simples = simple_modules(&a2).unwrap();
        let a2_plus = all_modules_of_length_plus_one(&simples).unwrap();
        let a2_at_most_two = all_modules_of_length_at_most(&a2, &Integer::from(2)).unwrap();
        let a2_indec_two = all_indec_modules_of_length_at_most(&a2, &Integer::from(2)).unwrap();
        let a2_negative = all_modules_of_length_at_most(&a2, &Integer::from(-1)).unwrap();
        let a2_zero = all_modules_of_length_at_most(&a2, &Integer::from(0)).unwrap();

        let loop_parent = PathAlgebra::new(
            field,
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_simples = simple_modules(&dual).unwrap();
        let dual_length_two = all_modules_of_length_plus_one(&dual_simples).unwrap();
        let dual_length_three = all_modules_of_length_plus_one(&dual_length_two).unwrap();

        json!({
            "a2_plus": length_module_list_json(&a2_plus, coefficient),
            "a2_at_most_two": length_module_list_json(&a2_at_most_two, coefficient),
            "a2_indec_two": length_module_list_json(&a2_indec_two, coefficient),
            "a2_negative": length_module_list_json(&a2_negative, coefficient),
            "a2_zero": length_module_list_json(&a2_zero, coefficient),
            "dual_length_two": length_module_list_json(&dual_length_two, coefficient),
            "dual_length_three": length_module_list_json(&dual_length_three, coefficient),
        })
    }

    #[test]
    fn pinned_qpa_modulehom_length_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_length.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let rational_a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                Q,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let rational_zero = zero_module(&rational_a2).unwrap();
        let rational_projective = right_module_over_path_algebra(
            rational_a2,
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![Q.one()]], Q).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let actual = json!({
            "schema_version": 1,
            "kind": "modulehom_length",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": modulehom_length_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": modulehom_length_case(field5, |entry| json!(
                field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "rational_indecomposable": {
                "zero": is_indecomposable_module(&rational_zero).unwrap().unwrap(),
                "projective": is_indecomposable_module(&rational_projective).unwrap().unwrap(),
            },
        });

        assert_eq!(actual, expected);
    }

    #[test]
    fn module_length_enumeration_preserves_qpa_validation_order() {
        assert!(matches!(
            all_modules_of_length_plus_one::<Z2>(&[]),
            Err(ModuleLengthError::EmptyModuleList)
        ));

        let rational = QuiverAlgebra::Path(
            PathAlgebra::new(Q, Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap()).unwrap(),
        );
        assert!(matches!(
            all_modules_of_length_at_most(&rational, &Integer::from(1)),
            Err(ModuleLengthError::FiniteFieldRequired)
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
                Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
            )
            .unwrap(),
        );
        let first = simple_modules(&first_algebra).unwrap()[0].clone();
        let foreign = simple_modules(&second_algebra).unwrap()[0].clone();
        assert!(matches!(
            all_modules_of_length_plus_one(&[first.clone(), foreign]),
            Err(ModuleLengthError::ModulesOverDifferentAlgebras)
        ));
        let double = direct_sum_of_qpa_modules(&[first.clone(), first.clone()])
            .unwrap()
            .unwrap();
        assert!(matches!(
            all_modules_of_length_plus_one(&[first, double]),
            Err(ModuleLengthError::ModulesOfDifferentLengths)
        ));
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn module_complexity_case<F: QpaCoefficientField>(
        field: F,
        include_finite_only_case: bool,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let zero_module = zero_module(&a2).unwrap();
        let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let projective = right_module_over_path_algebra(
            a2.clone(),
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let semisimple_pair = right_module_over_path_algebra(
            a2,
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        let loop_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_simple = right_module_over_path_algebra(
            dual.clone(),
            vec![1],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let dual_regular = right_module_over_path_algebra(
            dual.clone(),
            vec![2],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(
                    vec![vec![zero.clone(), one.clone()], vec![zero.clone(), zero]],
                    field,
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        let common = json!({
            "exceptional_zero": is_exceptional_module(&zero_module).unwrap().unwrap(),
            "exceptional_simple_one": is_exceptional_module(&simple_one).unwrap().unwrap(),
            "exceptional_projective": is_exceptional_module(&projective).unwrap().unwrap(),
            "exceptional_dual_simple": is_exceptional_module(&dual_simple).unwrap().unwrap(),
            "complexity_a2_simple_n2": complexity_of_module(&simple_one, &Integer::from(2)).unwrap(),
            "complexity_dual_simple_n2": complexity_of_module(&dual_simple, &Integer::from(2)).unwrap(),
            "complexity_dual_simple_n3": complexity_of_module(&dual_simple, &Integer::from(3)).unwrap(),
            "complexity_dual_simple_n10": complexity_of_module(&dual_simple, &Integer::from(10)).unwrap(),
            "complexity_dual_regular_n2": complexity_of_module(&dual_regular, &Integer::from(2)).unwrap(),
            "complexity_dual_algebra_n2": complexity_of_algebra(&dual, &Integer::from(2)).unwrap(),
            "complexity_dual_algebra_n10": complexity_of_algebra(&dual, &Integer::from(10)).unwrap(),
        });
        if include_finite_only_case {
            json!({
                "common": common,
                "exceptional_semisimple_pair": is_exceptional_module(&semisimple_pair)
                    .unwrap().unwrap(),
            })
        } else {
            common
        }
    }

    #[test]
    fn pinned_qpa_module_complexity_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/module_complexity.schema-v1.json"
        ))
        .unwrap();
        let actual = json!({
            "schema_version": 1,
            "kind": "module_complexity",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": module_complexity_case(Z2, true),
            "quinary": module_complexity_case(Zp::new(5), true),
            "rational": module_complexity_case(Q, false),
        });

        assert_eq!(actual, expected);
    }

    #[test]
    fn module_complexity_preserves_qpa_dispatch_and_cache_order() {
        let loop_parent = PathAlgebra::new(
            Q,
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_simple = simple_modules(&dual).unwrap()[0].clone();
        let dual_zero = zero_module(&dual).unwrap();
        assert_eq!(
            complexity_of_module(&dual_simple, &Integer::from(-2)).unwrap(),
            0
        );
        assert_eq!(
            complexity_of_module(&dual_zero, &Integer::from(-1)).unwrap(),
            0
        );
        assert!(matches!(
            complexity_of_module(&dual_simple, &Integer::from(1)),
            Err(EnvelopingModuleError::UndefinedComplexityEstimate { .. })
        ));

        let nonadmissible_parent = PathAlgebra::new(
            Q,
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "y")]).unwrap(),
        )
        .unwrap();
        let y = nonadmissible_parent.generators()[1].clone();
        let nonadmissible =
            QuiverAlgebra::Quotient(quotient_by_relations(nonadmissible_parent, &[y]).unwrap());
        assert!(matches!(
            complexity_of_algebra(&nonadmissible, &Integer::from(2)),
            Err(EnvelopingModuleError::ComplexityRequiresAdmissibleQuotient)
        ));

        let a2_parent = PathAlgebra::new(
            Q,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap();
        let a2_quotient = QuiverAlgebra::Quotient(quotient_by_relations(a2_parent, &[]).unwrap());
        assert!(matches!(
            global_dimension_of_algebra(&a2_quotient, 2).unwrap(),
            Some(GlobalDimensionOfAlgebraResult::Finite(1))
        ));
        assert_eq!(
            complexity_of_algebra(&a2_quotient, &Integer::from(1)).unwrap(),
            0
        );

        let rational_pair = right_module_over_path_algebra(
            QuiverAlgebra::Path(
                PathAlgebra::new(
                    Q,
                    Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")])
                        .unwrap(),
                )
                .unwrap(),
            ),
            vec![0, 2],
            &[],
        )
        .unwrap()
        .unwrap();
        assert_eq!(is_exceptional_module(&rational_pair).unwrap(), None);
        let rational_simple = right_module_over_path_algebra(
            rational_pair.right_acting_algebra().clone(),
            vec![1, 0],
            &[],
        )
        .unwrap()
        .unwrap();
        assert_eq!(is_exceptional_module(&rational_simple).unwrap(), Some(true));
        assert_eq!(rational_simple.cached_is_exceptional_module(), Some(true));
        assert_eq!(is_exceptional_module(&rational_simple).unwrap(), Some(true));
    }

    fn periodicity_json(period: Option<Integer>) -> Value {
        period.map_or(Value::Bool(false), |degree| {
            serde_json::from_str(&degree.to_string()).unwrap()
        })
    }

    #[allow(clippy::needless_pass_by_value)]
    fn modulehomalg_periodicity_case<F: QpaCoefficientField>(field: F) -> Value {
        let zero = field.zero();
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let zero_module = zero_module(&a2).unwrap();
        let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(a2, vec![0, 1], &[])
            .unwrap()
            .unwrap();

        let loop_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_simple = right_module_over_path_algebra(
            dual,
            vec![1],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(vec![vec![zero]], field).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        json!({
            "omega_zero_n0": periodicity_json(is_omega_periodic(&zero_module, &Integer::from(0)).unwrap()),
            "omega_zero_negative": periodicity_json(is_omega_periodic(&zero_module, &Integer::from(-2)).unwrap()),
            "omega_zero_n3": periodicity_json(is_omega_periodic(&zero_module, &Integer::from(3)).unwrap()),
            "omega_a2_simple_one_n4": periodicity_json(is_omega_periodic(&simple_one, &Integer::from(4)).unwrap()),
            "omega_a2_simple_two_n4": periodicity_json(is_omega_periodic(&simple_two, &Integer::from(4)).unwrap()),
            "omega_dual_simple_n1": periodicity_json(is_omega_periodic(&dual_simple, &Integer::from(1)).unwrap()),
            "omega_dual_simple_n4": periodicity_json(is_omega_periodic(&dual_simple, &Integer::from(4)).unwrap()),
            "tau_zero_n0": periodicity_json(is_tau_periodic(&zero_module, &Integer::from(0)).unwrap()),
            "tau_zero_negative": periodicity_json(is_tau_periodic(&zero_module, &Integer::from(-2)).unwrap()),
            "tau_zero_n3": periodicity_json(is_tau_periodic(&zero_module, &Integer::from(3)).unwrap()),
            "tau_a2_simple_one_n4": periodicity_json(is_tau_periodic(&simple_one, &Integer::from(4)).unwrap()),
            "tau_a2_simple_two_n4": periodicity_json(is_tau_periodic(&simple_two, &Integer::from(4)).unwrap()),
            "tau_dual_simple_n1": periodicity_json(is_tau_periodic(&dual_simple, &Integer::from(1)).unwrap()),
            "tau_dual_simple_n4": periodicity_json(is_tau_periodic(&dual_simple, &Integer::from(4)).unwrap()),
        })
    }

    #[test]
    fn pinned_qpa_modulehomalg_periodicity_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehomalg_periodicity.schema-v1.json"
        ))
        .unwrap();
        let actual = json!({
            "schema_version": 1,
            "kind": "modulehomalg_periodicity",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": modulehomalg_periodicity_case(Z2),
            "quinary": modulehomalg_periodicity_case(Zp::new(5)),
            "rational": modulehomalg_periodicity_case(Q),
        });

        assert_eq!(actual, expected);
    }

    #[test]
    fn omega_periodicity_exposes_qpa_progress_without_global_output() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let simple_one = right_module_over_path_algebra(algebra, vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let mut progress = Vec::new();
        assert_eq!(
            is_omega_periodic_with_progress(&simple_one, &Integer::from(4), |degree| {
                progress.push(degree.clone());
            })
            .unwrap(),
            None
        );
        assert_eq!(progress, (1..=4).map(Integer::from).collect::<Vec<_>>());

        progress.clear();
        assert_eq!(
            is_omega_periodic_with_progress(&simple_one, &Integer::from(-1), |degree| {
                progress.push(degree.clone());
            })
            .unwrap(),
            None
        );
        assert!(progress.is_empty());
    }

    fn double_star_map_json<F: Field>(
        map: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!({
            "source_dimension": map.source().dimension_vector(),
            "range_dimension": map.range().dimension_vector(),
            "matrices": map.vertex_matrices().iter().map(|matrix| json!(
                (0..matrix.nrows()).map(|row| matrix.row(row).iter()
                    .map(coefficient).collect::<Vec<_>>()).collect::<Vec<_>>()
            )).collect::<Vec<_>>(),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehomalg_approximation_case<F: QpaCoefficientField>(
        field: F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let a2_zero = zero_module(&a2).unwrap();
        let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let projective = right_module_over_path_algebra(
            a2,
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let repeated = direct_sum_of_qpa_modules(&[simple_one.clone(), simple_one.clone()])
            .unwrap()
            .unwrap();

        let loop_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_simple = right_module_over_path_algebra(
            dual.clone(),
            vec![1],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let dual_regular = right_module_over_path_algebra(
            dual,
            vec![2],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(
                    vec![vec![zero.clone(), one.clone()], vec![zero.clone(), zero]],
                    field,
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        let evaluate = |map| double_star_map_json(&map, coefficient);
        json!({
            "right_zero": evaluate(right_approximation_by_add_m(&simple_one, &a2_zero).unwrap()),
            "right_projective_to_simple": evaluate(right_approximation_by_add_m(&projective, &simple_one).unwrap()),
            "right_simple_to_projective": evaluate(right_approximation_by_add_m(&simple_one, &projective).unwrap()),
            "right_repeated": evaluate(right_approximation_by_add_m(&repeated, &simple_one).unwrap()),
            "minimal_right_repeated": evaluate(minimal_right_approximation(&repeated, &simple_one).unwrap()),
            "left_zero": evaluate(left_approximation_by_add_m(&a2_zero, &simple_one).unwrap()),
            "left_projective_to_simple": evaluate(left_approximation_by_add_m(&projective, &simple_one).unwrap()),
            "left_simple_to_projective": evaluate(left_approximation_by_add_m(&simple_one, &projective).unwrap()),
            "left_repeated": evaluate(left_approximation_by_add_m(&simple_one, &repeated).unwrap()),
            "minimal_left_repeated": evaluate(minimal_left_approximation(&simple_one, &repeated).unwrap()),
            "right_dual_regular_to_simple": evaluate(right_approximation_by_add_m(&dual_regular, &dual_simple).unwrap()),
            "left_dual_simple_to_regular": evaluate(left_approximation_by_add_m(&dual_simple, &dual_regular).unwrap()),
        })
    }

    #[test]
    fn pinned_qpa_modulehomalg_approximation_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehomalg_approximations.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "modulehomalg_approximations",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": modulehomalg_approximation_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": modulehomalg_approximation_case(field5, |entry| json!(
                field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "rational": modulehomalg_approximation_case(Q, |entry| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });

        assert_eq!(actual, expected);
    }

    #[test]
    fn add_m_approximations_reject_foreign_algebras_before_zero_shortcuts() {
        let first = isolated_path_algebra(Z2, 1);
        let second = isolated_path_algebra(Z2, 1);
        let first_zero = zero_module(&first).unwrap();
        let second_zero = zero_module(&second).unwrap();
        assert!(matches!(
            right_approximation_by_add_m(&first_zero, &second_zero),
            Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
        ));
        assert!(matches!(
            left_approximation_by_add_m(&first_zero, &second_zero),
            Err(EnvelopingModuleError::ApproximationModulesOverDifferentAlgebras)
        ));
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_double_star_case<F: QpaCoefficientField>(
        field: F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let a2_zero = zero_module(&a2).unwrap();
        let a2_simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let a2_simple_two = right_module_over_path_algebra(a2.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let a2_projective = right_module_over_path_algebra(
            a2.clone(),
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one.clone()]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let a2_semisimple_pair = right_module_over_path_algebra(
            a2,
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        let loop_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_simple = right_module_over_path_algebra(
            dual.clone(),
            vec![1],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(vec![vec![zero.clone()]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let dual_regular = right_module_over_path_algebra(
            dual,
            vec![2],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(
                    vec![vec![zero.clone(), one.clone()], vec![zero.clone(), zero]],
                    field,
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        let evaluate = |module: &PathAlgebraModule<F>| {
            double_star_map_json(
                &from_identity_to_double_star_homomorphism(module).unwrap(),
                coefficient,
            )
        };
        json!({
            "a2_zero": evaluate(&a2_zero),
            "a2_simple_one": evaluate(&a2_simple_one),
            "a2_simple_two": evaluate(&a2_simple_two),
            "a2_projective": evaluate(&a2_projective),
            "a2_semisimple_pair": evaluate(&a2_semisimple_pair),
            "dual_simple": evaluate(&dual_simple),
            "dual_regular": evaluate(&dual_regular),
        })
    }

    #[test]
    fn pinned_qpa_modulehom_double_star_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_double_star.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "modulehom_double_star",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": modulehom_double_star_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": modulehom_double_star_case(field5, |entry| json!(
                field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "rational": modulehom_double_star_case(Q, |entry| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });

        assert_eq!(actual, expected);
    }

    #[test]
    fn identity_to_double_star_homomorphism_is_cached_with_star_modules() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Q,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let module = right_module_over_path_algebra(
            algebra,
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![Q.one()]], Q).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let first = from_identity_to_double_star_homomorphism(&module).unwrap();
        let second = from_identity_to_double_star_homomorphism(&module).unwrap();
        let double_star = star_of_module(&star_of_module(&module).unwrap()).unwrap();

        assert!(first.is_identical(&second));
        assert!(first.source().is_identical(&module));
        assert!(first.range().is_identical(&double_star));
    }

    fn almost_split_sequence_json<F: Field>(
        sequence: Option<PathAlgebraModuleHomomorphismPair<F>>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        match sequence {
            None => Value::Null,
            Some((monomorphism, epimorphism)) => json!({
                "monomorphism": double_star_map_json(&monomorphism, coefficient),
                "epimorphism": double_star_map_json(&epimorphism, coefficient),
            }),
        }
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn module_ar_almost_split_case<F: QpaCoefficientField>(
        field: F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(a2.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let projective = right_module_over_path_algebra(
            a2,
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        let loop_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_simple = right_module_over_path_algebra(
            dual,
            vec![1],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(vec![vec![zero]], field).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        json!({
            "ending_simple_one": almost_split_sequence_json(
                almost_split_sequence(&simple_one).unwrap(), coefficient),
            "ending_simple_two": almost_split_sequence_json(
                almost_split_sequence(&simple_two).unwrap(), coefficient),
            "ending_projective": almost_split_sequence_json(
                almost_split_sequence(&projective).unwrap(), coefficient),
            "ending_dual_simple": almost_split_sequence_json(
                almost_split_sequence(&dual_simple).unwrap(), coefficient),
            "starting_simple_one": almost_split_sequence_json(
                almost_split_sequence_from_side(&simple_one, AlmostSplitSide::Starting).unwrap(),
                coefficient),
            "starting_simple_two": almost_split_sequence_json(
                almost_split_sequence_from_side(&simple_two, AlmostSplitSide::Starting).unwrap(),
                coefficient),
            "starting_dual_simple": almost_split_sequence_json(
                almost_split_sequence_from_side(&dual_simple, AlmostSplitSide::Starting).unwrap(),
                coefficient),
        })
    }

    #[test]
    fn pinned_qpa_module_ar_almost_split_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/module_ar_almost_split.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "module_ar_almost_split",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": module_ar_almost_split_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": module_ar_almost_split_case(field5, |entry| json!(
                field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "rational": module_ar_almost_split_case(Q, |entry| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });

        assert_eq!(actual, expected);
    }

    #[test]
    fn almost_split_sequence_caches_values_and_projective_fail() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let nonprojective = simple_modules(&algebra).unwrap()[0].clone();
        let projective = simple_modules(&algebra).unwrap()[1].clone();
        let first = almost_split_sequence(&nonprojective).unwrap().unwrap();
        let second = almost_split_sequence(&nonprojective).unwrap().unwrap();

        assert!(first.0.is_identical(&second.0));
        assert!(first.1.is_identical(&second.1));
        assert!(almost_split_sequence(&projective).unwrap().is_none());
        assert!(matches!(
            projective.cached_almost_split_sequence(),
            CachedAttribute::Fail
        ));
        assert!(almost_split_sequence(&projective).unwrap().is_none());
    }

    fn irreducible_maps_json<F: Field>(
        maps: &[PathAlgebraModuleHomomorphism<F>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!(
            maps.iter()
                .map(|map| double_star_map_json(map, coefficient))
                .collect::<Vec<_>>()
        )
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn module_ar_irreducible_case<F: QpaCoefficientField>(
        field: F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(a2.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let projective = right_module_over_path_algebra(
            a2,
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let loop_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_simple = right_module_over_path_algebra(
            dual,
            vec![1],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(vec![vec![zero]], field).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        json!({
            "ending_simple_one": irreducible_maps_json(
                &irreducible_morphisms_ending_in(&simple_one).unwrap(), coefficient),
            "ending_projective": irreducible_maps_json(
                &irreducible_morphisms_ending_in(&projective).unwrap(), coefficient),
            "ending_dual_simple": irreducible_maps_json(
                &irreducible_morphisms_ending_in(&dual_simple).unwrap(), coefficient),
            "starting_simple_two": irreducible_maps_json(
                &irreducible_morphisms_starting_in(&simple_two).unwrap(), coefficient),
            "starting_dual_simple": irreducible_maps_json(
                &irreducible_morphisms_starting_in(&dual_simple).unwrap(), coefficient),
        })
    }

    #[test]
    fn pinned_qpa_module_ar_irreducible_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/module_ar_irreducible.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "module_ar_irreducible",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": module_ar_irreducible_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": module_ar_irreducible_case(field5, |entry| json!(
                field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
        });

        assert_eq!(actual, expected);
    }

    #[test]
    fn irreducible_morphisms_preserve_qpa_field_error_and_live_cache() {
        let binary = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let simples = simple_modules(&binary).unwrap();
        let first = irreducible_morphisms_ending_in(&simples[0]).unwrap();
        let second = irreducible_morphisms_ending_in(&simples[0]).unwrap();
        assert_eq!(first.len(), second.len());
        assert!(
            first
                .iter()
                .zip(&second)
                .all(|(left, right)| left.is_identical(right))
        );
        assert!(matches!(
            irreducible_morphisms_ending_in(&simples[1]),
            Err(super::AuslanderReitenError::ZeroIrreducibleSourceDecomposition)
        ));

        let rational = QuiverAlgebra::Path(
            PathAlgebra::new(
                Q,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let rational_simple = simple_modules(&rational).unwrap()[0].clone();
        assert!(matches!(
            irreducible_morphisms_ending_in(&rational_simple),
            Err(super::AuslanderReitenError::FiniteFieldRequired)
        ));
        assert!(matches!(
            irreducible_morphisms_starting_in(&rational_simple),
            Err(super::AuslanderReitenError::FiniteFieldRequired)
        ));
    }

    fn predecessor_module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!({
            "dimension": module.dimension_vector(),
            "matrices": module.arrow_matrices().iter().map(|matrix| json!(
                (0..matrix.nrows()).map(|row| matrix.row(row).iter()
                    .map(coefficient).collect::<Vec<_>>()).collect::<Vec<_>>()
            )).collect::<Vec<_>>(),
        })
    }

    fn predecessors_json<F: Field>(
        result: PredecessorsOfModuleResult<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        match result {
            PredecessorsOfModuleResult::Module(module) => json!({
                "kind": "module",
                "module": predecessor_module_json(&module, coefficient),
            }),
            PredecessorsOfModuleResult::Layers(result) => json!({
                "kind": "layers",
                "layers": result.layers.iter().map(|layer| layer.iter()
                    .map(|module| predecessor_module_json(module, coefficient))
                    .collect::<Vec<_>>()).collect::<Vec<_>>(),
                "valuations": result.valuations.iter().map(|layer| layer.iter()
                    .map(|valuation| {
                        let multiplicity = match valuation.multiplicity {
                            PredecessorMultiplicity::Pair { first, second } => json!([
                                first.map_or(Value::Bool(false), Value::from),
                                second.map_or(Value::Bool(false), Value::from),
                            ]),
                            PredecessorMultiplicity::Scalar(value) => json!(value),
                        };
                        json!([
                            valuation.predecessor_position,
                            valuation.successor_position,
                            multiplicity,
                        ])
                    }).collect::<Vec<_>>()).collect::<Vec<_>>(),
            }),
        }
    }

    #[allow(clippy::needless_pass_by_value)]
    fn module_ar_predecessors_case<F: QpaCoefficientField>(
        field: F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let simple_one = right_module_over_path_algebra(a2, vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let loop_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_simple = right_module_over_path_algebra(
            dual,
            vec![1],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(vec![vec![field.zero()]], field).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        json!({
            "a2_n1": predecessors_json(
                predecessors_of_module(&simple_one, &Integer::from(1)).unwrap(), coefficient),
            "a2_n2": predecessors_json(
                predecessors_of_module(&simple_one, &Integer::from(2)).unwrap(), coefficient),
            "a2_n3": predecessors_json(
                predecessors_of_module(&simple_one, &Integer::from(3)).unwrap(), coefficient),
            "a2_n4": predecessors_json(
                predecessors_of_module(&simple_one, &Integer::from(4)).unwrap(), coefficient),
            "dual_n2": predecessors_json(
                predecessors_of_module(&dual_simple, &Integer::from(2)).unwrap(), coefficient),
            "dual_n3": predecessors_json(
                predecessors_of_module(&dual_simple, &Integer::from(3)).unwrap(), coefficient),
            "dual_n4": predecessors_json(
                predecessors_of_module(&dual_simple, &Integer::from(4)).unwrap(), coefficient),
        })
    }

    #[test]
    fn pinned_qpa_module_ar_predecessors_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/module_ar_predecessors.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "module_ar_predecessors",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": module_ar_predecessors_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": module_ar_predecessors_case(field5, |entry| json!(
                field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
        });

        assert_eq!(actual, expected);
    }

    #[test]
    fn predecessors_preserve_projective_check_and_n_one_return() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let simples = simple_modules(&algebra).unwrap();
        assert!(matches!(
            predecessors_of_module(&simples[1], &Integer::from(1)),
            Err(super::AuslanderReitenError::PredecessorModuleProjective)
        ));
        assert!(matches!(
            predecessors_of_module(&simples[0], &Integer::from(0)),
            Err(super::AuslanderReitenError::PredecessorBoundBelowOne { .. })
        ));
        let PredecessorsOfModuleResult::Module(returned) =
            predecessors_of_module(&simples[0], &Integer::from(1)).unwrap()
        else {
            panic!("QPA must return the module itself at n = 1");
        };
        assert!(returned.is_identical(&simples[0]));
    }

    fn module_map_json<F: Field>(
        map: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "source_dimension": map.source().dimension_vector(),
            "range_dimension": map.range().dimension_vector(),
            "matrices": map
                .vertex_matrices()
                .iter()
                .map(|matrix| matrix_json(matrix, coefficient))
                .collect::<Vec<_>>(),
        })
    }

    fn projective_module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "dimension": module.dimension_vector(),
            "matrices": module
                .arrow_matrices()
                .iter()
                .map(|matrix| matrix_json(matrix, coefficient))
                .collect::<Vec<_>>(),
        })
    }

    fn module_elements_json<F: Field>(
        elements: &[PathAlgebraModuleElement<F>],
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!(
            elements
                .iter()
                .map(|element| element
                    .external_representation()
                    .iter()
                    .map(|row| row.iter().map(coefficient).collect::<Vec<_>>())
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>()
        )
    }

    fn path_element_json<F: Field>(
        element: &qpa_core::PathAlgebraElement<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!(
            element
                .external_representation()
                .terms()
                .iter()
                .map(|(walk, scalar)| json!({
                    "walk": walk,
                    "coefficient": coefficient(scalar),
                }))
                .collect::<Vec<_>>()
        )
    }

    fn path_element_matrix_json<F: Field>(
        matrix: &[Vec<qpa_core::PathAlgebraElement<F>>],
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!(
            matrix
                .iter()
                .map(|row| row
                    .iter()
                    .map(|element| path_element_json(element, coefficient))
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>()
        )
    }

    fn physical_module_element_json<F: Field>(
        element: &PathAlgebraModuleElement<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "dimension": element.module().dimension_vector(),
            "support": element.support_vertices(),
            "coordinates": element
                .flat_external_representation()
                .iter()
                .map(coefficient)
                .collect::<Vec<_>>(),
        })
    }

    fn entered_path_matrix_equals<F: Field>(
        returned: &[Vec<qpa_core::PathAlgebraElement<F>>],
        entered: &[Vec<QuiverAlgebraElement<F>>],
    ) -> bool {
        returned.len() == entered.len()
            && returned
                .iter()
                .zip(entered)
                .all(|(returned_row, entered_row)| {
                    returned_row.len() == entered_row.len()
                        && returned_row.iter().zip(entered_row).all(
                            |(returned_element, entered_element)| {
                                matches!(entered_element, QuiverAlgebraElement::Path(element)
                                if returned_element == element)
                            },
                        )
                })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_functor_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one.clone()], vec![zero.clone()]], field.clone())
                    .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let map = right_module_hom_over_algebra(
            module.clone(),
            module.clone(),
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![one.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::identity(1, field.clone()),
            ],
        )
        .unwrap();
        let minimal = minimal_generating_set_of_module(module.clone()).unwrap();
        let top_projection = top_of_module_projection(module.clone()).unwrap();
        let top = top_of_module(module.clone()).unwrap();
        let socle_inclusion = socle_of_module_inclusion(module.clone()).unwrap();
        let dual_map = dual_of_module_homomorphism(&map).unwrap();
        let target = right_module_over_path_algebra(
            algebra.clone(),
            vec![1, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one.clone(), one]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let zero_module = right_module_over_path_algebra(algebra, vec![0, 0], &[])
            .unwrap()
            .unwrap();
        let zero_socle = socle_of_module_inclusion(zero_module.clone()).unwrap();
        let zero_trace = trace_of_module(&simple_one, &target).unwrap();
        let zero_source_trace = trace_of_module(&zero_module, &target).unwrap();
        let zero_target_trace = trace_of_module(&simple_two, &zero_module).unwrap();
        let nonzero_trace = trace_of_module(&simple_two, &target).unwrap();
        let legacy_trace =
            right_fac_approximation_via_all_basis_images(&simple_two, &target).unwrap();
        let ordered_trace = trace_of_module(&module, &module).unwrap();
        let ordered_legacy =
            right_fac_approximation_via_all_basis_images(&module, &module).unwrap();

        json!({
            "minimal_generators": {
                "count": minimal.len(),
                "elements": module_elements_json(&minimal, coefficient),
                "zero_count": minimal_generating_set_of_module(zero_module.clone()).unwrap().len(),
                "top_projection_cached_identity": top_projection.is_identical(
                    &top_of_module_projection(module.clone()).unwrap()),
                "top_cached_identity": top.is_identical(&top_of_module(module).unwrap()),
                "top_projection_range_identity": top_projection.range().is_identical(&top),
            },
            "socle": {
                "inclusion": module_map_json(&socle_inclusion, coefficient),
                "object_dimension": socle_of_module(socle_inclusion.range().clone())
                    .unwrap()
                    .dimension_vector(),
                "zero_inclusion": module_map_json(&zero_socle, coefficient),
                "zero_source_identity": zero_socle.source().is_identical(&zero_module),
            },
            "dual_map": {
                "cached_identity": dual_map.is_identical(
                    &dual_of_module_homomorphism(&map).unwrap()),
                "source_dimension": dual_map.source().dimension_vector(),
                "range_dimension": dual_map.range().dimension_vector(),
                "matrices": dual_map.vertex_matrices().iter()
                    .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
            },
            "trace": {
                "zero": module_map_json(&zero_trace, coefficient),
                "zero_source": module_map_json(&zero_source_trace, coefficient),
                "zero_target": module_map_json(&zero_target_trace, coefficient),
                "nonzero": module_map_json(&nonzero_trace, coefficient),
                "legacy": module_map_json(&legacy_trace, coefficient),
                "legacy_equal": legacy_trace == nonzero_trace,
                "ordered": module_map_json(&ordered_trace, coefficient),
                "ordered_legacy": module_map_json(&ordered_legacy, coefficient),
                "ordered_legacy_equal": ordered_legacy == ordered_trace,
                "right_fac_equal": right_fac_approximation(&simple_two, &target).unwrap()
                    == nonzero_trace,
            },
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_projective_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let projectives = indec_projective_modules(&algebra).unwrap();
        let bases = basis_of_projectives(&algebra).unwrap().unwrap();
        let cached_projectives = indec_projective_modules(&algebra).unwrap();
        let basis_paths = bases
            .iter()
            .map(|projective| {
                projective
                    .iter()
                    .map(|vertex_basis| {
                        vertex_basis
                            .iter()
                            .map(|element| element.leading_monomial().to_string())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(
                    vec![
                        vec![one.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let generators = minimal_generating_set_of_module(module.clone()).unwrap();
        let cover = projective_cover(&module).unwrap();
        let from_generator = hom_from_projective(&generators[0], &module)
            .unwrap()
            .unwrap();
        let projective_sum = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
        let lifting = lifting_morphism_from_projective(&cover, &cover)
            .unwrap()
            .unwrap();
        let nondirect_lifting = lifting_morphism_from_projective(&cover, &from_generator).unwrap();
        let zero_module = right_module_over_path_algebra(algebra.clone(), vec![0, 0], &[])
            .unwrap()
            .unwrap();
        let zero_to_module = zero_module_mapping(zero_module.clone(), module.clone()).unwrap();
        let zero_lifting = lifting_morphism_from_projective(&cover, &zero_to_module)
            .unwrap()
            .unwrap();
        let factoring = hom_factoring_through_projectives(&module, &module).unwrap();
        let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let layer = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
            .unwrap()
            .unwrap();
        let first = right_module_hom_over_algebra(
            layer.clone(),
            module.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone())
                    .unwrap(),
                DenseMatrix::from_rows(vec![vec![one.clone(), zero.clone()]], field.clone())
                    .unwrap(),
            ],
        )
        .unwrap();
        let second = right_module_hom_over_algebra(
            layer,
            module.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone())
                    .unwrap(),
                DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone())
                    .unwrap(),
            ],
        )
        .unwrap();
        let kernel_module = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let kernel_inclusion = right_module_hom_over_algebra(
            kernel_module,
            module.clone(),
            vec![
                DenseMatrix::from_rows(vec![vec![zero.clone(), one.clone()]], field.clone())
                    .unwrap(),
                DenseMatrix::zero(1, 2, field.clone()),
            ],
        )
        .unwrap();
        let inclusion_lift = lifting_inclusion_morphisms(&first, &kernel_inclusion)
            .unwrap()
            .unwrap();
        let inclusion_fail = lifting_inclusion_morphisms(&first, &second).unwrap();

        let loop_quiver =
            Quiver::from_named(vec!["w".into()], vec![NamedArrowSpec::named("w", "w", "x")])
                .unwrap();
        let path_algebra = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let loop_element = path_algebra
            .element(path_algebra.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(
                path_algebra,
                &[loop_element.multiply(&loop_element).unwrap()],
            )
            .unwrap(),
        );
        let dual_projective = indec_projective_modules(&dual_numbers).unwrap()[0].clone();
        let dual_basis = basis_of_projectives(&dual_numbers).unwrap().unwrap();
        let dual_basis_paths = dual_basis
            .iter()
            .map(|projective| {
                projective
                    .iter()
                    .map(|vertex_basis| {
                        vertex_basis
                            .iter()
                            .map(|element| element.leading_monomial().to_string())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let dual_simple = right_module_over_path_algebra(dual_numbers, vec![1], &[])
            .unwrap()
            .unwrap();
        let dual_cover = projective_cover(&dual_simple).unwrap();

        json!({
            "a2": {
                "basis_paths": basis_paths,
                "projectives": projectives
                    .iter()
                    .map(|projective| projective_module_json(projective, coefficient))
                    .collect::<Vec<_>>(),
                "cache_identity": projectives
                    .iter()
                    .zip(&cached_projectives)
                    .map(|(projective, cached)| projective.is_identical(cached))
                    .collect::<Vec<_>>(),
                "installed_properties": [
                    projectives[0].cached_is_projective_module() == Some(true),
                    projectives[1].cached_is_projective_module() == Some(true),
                    projectives[0].cached_is_indecomposable_module() == Some(true),
                    projectives[1].cached_is_indecomposable_module() == Some(true),
                ],
                "projectivity": [
                    is_projective_module(&projective_sum).unwrap(),
                    is_projective_module(&simple_one).unwrap(),
                    is_projective_module(&simple_two).unwrap(),
                    is_projective_module(&zero_module).unwrap(),
                ],
                "minimal_generator_count": generators.len(),
                "from_generator": module_map_json(&from_generator, coefficient),
                "cover": module_map_json(&cover, coefficient),
                "cover_cache_identity": cover.is_identical(&projective_cover(&module).unwrap()),
                "lifting": module_map_json(&lifting, coefficient),
                "lifting_triangle": lifting.compose(&cover).unwrap() == cover,
                "nondirect_lifting_is_fail": nondirect_lifting.is_none(),
                "zero_lifting": module_map_json(&zero_lifting, coefficient),
                "factoring_count": factoring.len(),
                "factoring": factoring
                    .iter()
                    .map(|map| module_map_json(map, coefficient))
                    .collect::<Vec<_>>(),
                "factoring_zero_counts": [
                    hom_factoring_through_projectives(&zero_module, &module).unwrap().len(),
                    hom_factoring_through_projectives(&module, &zero_module).unwrap().len(),
                ],
                "inclusion_lift": module_map_json(&inclusion_lift, coefficient),
                "inclusion_triangle": inclusion_lift.compose(&first).unwrap() == kernel_inclusion,
                "inclusion_fail": inclusion_fail.is_none(),
            },
            "dual_numbers": {
                "basis_paths": dual_basis_paths,
                "projective": projective_module_json(&dual_projective, coefficient),
                "generator_count": minimal_generating_set_of_module(dual_projective.clone())
                    .unwrap().len(),
                "is_projective": is_projective_module(&dual_projective).unwrap(),
                "simple_is_projective": is_projective_module(&dual_simple).unwrap(),
                "cover": module_map_json(&dual_cover, coefficient),
            },
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_projective_matrix_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let path_algebra = PathAlgebra::new(field.clone(), quiver).unwrap();
        let algebra = QuiverAlgebra::Path(path_algebra);
        let vertices = algebra.quiver().vertices().collect::<Vec<_>>();
        let arrow = algebra.quiver().arrows().next().unwrap();
        let v1 = include_in_path_algebra(&vertices[0].path(), &algebra).unwrap();
        let v2 = include_in_path_algebra(&vertices[1].path(), &algebra).unwrap();
        let a = include_in_path_algebra(&arrow.path(), &algebra).unwrap();
        let entered_matrix = vec![vec![v1.clone(), a.clone()], vec![algebra.zero(), v2]];
        let map =
            from_matrix_to_homomorphism_of_projectives(&algebra, &entered_matrix, &[1, 2], &[1, 2])
                .unwrap();
        let returned_matrix = matrix_of_homomorphism_between_projectives(&map).unwrap();
        let uniform = v1.add(&a).unwrap();
        let projective_element =
            element_in_v_a_as_element_in_indec_proj(&algebra, &uniform).unwrap();
        let recovered = element_in_indec_projective(&algebra, &projective_element, 1).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let sum = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
        let zero_module = right_module_over_path_algebra(algebra.clone(), vec![0, 0], &[])
            .unwrap()
            .unwrap();
        let zero_to_sum = matrix_of_homomorphism_between_projectives(
            &zero_module_mapping(zero_module.clone(), sum.clone()).unwrap(),
        )
        .unwrap();
        let sum_to_zero = matrix_of_homomorphism_between_projectives(
            &zero_module_mapping(sum, zero_module.clone()).unwrap(),
        )
        .unwrap();
        let zero_to_zero = matrix_of_homomorphism_between_projectives(
            &zero_module_mapping(zero_module.clone(), zero_module).unwrap(),
        )
        .unwrap();
        let single =
            from_matrix_to_homomorphism_of_projectives(&algebra, &[vec![a]], &[2], &[1]).unwrap();

        let quotient_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(3), "y"),
            ],
        )
        .unwrap();
        let quotient_parent = PathAlgebra::new(field, quotient_quiver).unwrap();
        let quotient_arrows = quotient_parent.quiver().arrows().collect::<Vec<_>>();
        let x = quotient_parent.element(quotient_arrows[0].path()).unwrap();
        let y = quotient_parent.element(quotient_arrows[1].path()).unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(quotient_parent.clone(), &[x.multiply(&y).unwrap()]).unwrap(),
        );
        let quotient_vertices = quotient.quiver().vertices().collect::<Vec<_>>();
        let qv1 = include_in_path_algebra(&quotient_vertices[0].path(), &quotient).unwrap();
        let qv2 = include_in_path_algebra(&quotient_vertices[1].path(), &quotient).unwrap();
        let qx = include_in_path_algebra(&quotient_arrows[0].path(), &quotient).unwrap();
        let qy = include_in_path_algebra(&quotient_arrows[1].path(), &quotient).unwrap();
        let quotient_matrix = vec![vec![qx.clone(), quotient.zero()], vec![qv2, qy]];
        let quotient_map = from_matrix_to_homomorphism_of_projectives(
            &quotient,
            &quotient_matrix,
            &[2, 3],
            &[1, 2],
        )
        .unwrap();
        let quotient_back = matrix_of_homomorphism_between_projectives(&quotient_map).unwrap();
        let quotient_uniform = qv1.add(&qx).unwrap();
        let quotient_projective_element =
            element_in_v_a_as_element_in_indec_proj(&quotient, &quotient_uniform).unwrap();
        let quotient_recovered =
            element_in_indec_projective(&quotient, &quotient_projective_element, 1).unwrap();

        json!({
            "path": {
                "map": module_map_json(&map, coefficient),
                "matrix": path_element_matrix_json(&returned_matrix, coefficient),
                "roundtrip_equal": entered_path_matrix_equals(
                    &returned_matrix,
                    &entered_matrix,
                ),
                "underlying": matrix_json(&underlying_linear_map(&map).unwrap(), coefficient),
                "uniform_projective_element": physical_module_element_json(
                    &projective_element,
                    coefficient,
                ),
                "uniform_recovered": path_element_json(&recovered, coefficient),
                "single": module_map_json(&single, coefficient),
                "single_underlying": matrix_json(
                    &underlying_linear_map(&single).unwrap(),
                    coefficient,
                ),
                "zero_to_sum": path_element_matrix_json(&zero_to_sum, coefficient),
                "sum_to_zero": path_element_matrix_json(&sum_to_zero, coefficient),
                "zero_to_zero": path_element_matrix_json(&zero_to_zero, coefficient),
            },
            "quotient": {
                "map": module_map_json(&quotient_map, coefficient),
                "matrix_parent_representatives": path_element_matrix_json(
                    &quotient_back,
                    coefficient,
                ),
                "underlying": matrix_json(
                    &underlying_linear_map(&quotient_map).unwrap(),
                    coefficient,
                ),
                "uniform_projective_element": physical_module_element_json(
                    &quotient_projective_element,
                    coefficient,
                ),
                "uniform_recovered_parent": path_element_json(
                    &quotient_recovered,
                    coefficient,
                ),
            },
        })
    }

    fn attribute_module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "dimension": module.dimension_vector(),
            "matrices": module
                .matrices_of_path_algebra_module()
                .iter()
                .map(|entry| match entry {
                    PathAlgebraModuleMatrix::ZeroScalar => json!("0"),
                    PathAlgebraModuleMatrix::Matrix(matrix) => matrix_json(matrix, coefficient),
                })
                .collect::<Vec<_>>(),
        })
    }

    fn attribute_module_list_json<F: Field>(
        modules: &[PathAlgebraModule<F>],
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Vec<Value> {
        modules
            .iter()
            .map(|module| attribute_module_json(module, coefficient))
            .collect()
    }

    fn star_module_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module.cached_star_of_module().is_some();
        let value = star_of_module(module).unwrap();
        let has_after = module.cached_star_of_module().is_some();
        let repeated = star_of_module(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn nakayama_module_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module.cached_nakayama_functor_of_module().is_some();
        let value = nakayama_functor_of_module(module).unwrap();
        let has_after = module.cached_nakayama_functor_of_module().is_some();
        let repeated = nakayama_functor_of_module(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn opposite_nakayama_module_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module
            .cached_opposite_nakayama_functor_of_module()
            .is_some();
        let value = opposite_nakayama_functor_of_module(module).unwrap();
        let has_after = module
            .cached_opposite_nakayama_functor_of_module()
            .is_some();
        let repeated = opposite_nakayama_functor_of_module(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn module_functor_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let star = star_module_attribute_json(module, coefficient);
        let nakayama = nakayama_module_attribute_json(module, coefficient);
        let opposite_nakayama = opposite_nakayama_module_attribute_json(module, coefficient);
        json!({
            "star": star,
            "nakayama": nakayama,
            "opposite_nakayama": opposite_nakayama,
        })
    }

    fn star_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism.cached_star_homomorphism().is_some();
        let value = star_of_module_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism.cached_star_homomorphism().is_some();
        let repeated = star_of_module_homomorphism(homomorphism).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn nakayama_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism
            .cached_nakayama_functor_homomorphism()
            .is_some();
        let value = nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism
            .cached_nakayama_functor_homomorphism()
            .is_some();
        let repeated = nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn opposite_nakayama_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism
            .cached_opposite_nakayama_functor_homomorphism()
            .is_some();
        let value = opposite_nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism
            .cached_opposite_nakayama_functor_homomorphism()
            .is_some();
        let repeated = opposite_nakayama_functor_of_module_homomorphism(homomorphism).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn map_functor_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let star = star_map_attribute_json(homomorphism, coefficient);
        let nakayama = nakayama_map_attribute_json(homomorphism, coefficient);
        let opposite_nakayama = opposite_nakayama_map_attribute_json(homomorphism, coefficient);
        json!({
            "star": star,
            "nakayama": nakayama,
            "opposite_nakayama": opposite_nakayama,
        })
    }

    fn transpose_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module.cached_transpose_of_module().is_some();
        let value = transpose_of_module(module).unwrap();
        let has_after = module.cached_transpose_of_module().is_some();
        let repeated = transpose_of_module(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
            "indecomposable_bound": value.cached_is_indecomposable_module().is_some(),
            "indecomposable_value": value.cached_is_indecomposable_module(),
        })
    }

    fn d_tr_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module.cached_d_tr().is_some();
        let value = d_tr(module).unwrap();
        let has_after = module.cached_d_tr().is_some();
        let repeated = d_tr(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn tr_d_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = module.cached_tr_d().is_some();
        let value = tr_d(module).unwrap();
        let has_after = module.cached_tr_d().is_some();
        let repeated = tr_d(module).unwrap();
        json!({
            "has_before": has_before,
            "value": attribute_module_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
        })
    }

    fn transpose_functor_attribute_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let transpose = transpose_attribute_json(module, coefficient);
        let dtr = d_tr_attribute_json(module, coefficient);
        let trd = tr_d_attribute_json(module, coefficient);
        let tau_has_before = module.cached_is_tau_rigid_module().is_some();
        let tau_value = is_tau_rigid_module(module).unwrap();
        let tau_has_after = module.cached_is_tau_rigid_module().is_some();
        json!({
            "transpose": transpose,
            "dtr": dtr,
            "trd": trd,
            "tau_rigid": {
                "has_before": tau_has_before,
                "value": tau_value,
                "has_after": tau_has_after,
            },
        })
    }

    fn functor_iterations_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let d_zero = d_tr_iterated(module, 0).unwrap();
        let d_one = d_tr_iterated(module, 1).unwrap();
        let d_one_generators = minimal_generating_set_of_module(d_one.clone()).unwrap();
        let d_minus_one = d_tr_iterated(module, -1).unwrap();
        let d_two = d_tr_iterated(module, 2).unwrap();
        let t_zero = tr_d_iterated(module, 0).unwrap();
        let t_one = tr_d_iterated(module, 1).unwrap();
        let t_minus_one = tr_d_iterated(module, -1).unwrap();
        let t_two = tr_d_iterated(module, 2).unwrap();
        json!({
            "dtr_zero_identity": d_zero.is_identical(module),
            "dtr_one_identity": d_one.is_identical(&d_tr(module).unwrap()),
            "dtr_one_generators": d_one_generators.iter().map(|element|
                element.external_representation().iter().map(|row|
                    row.iter().map(coefficient).collect::<Vec<_>>()).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
            "dtr_minus_one_identity": d_minus_one.is_identical(&tr_d(module).unwrap()),
            "dtr_two": attribute_module_json(&d_two, coefficient),
            "trd_zero_identity": t_zero.is_identical(module),
            "trd_one_identity": t_one.is_identical(&tr_d(module).unwrap()),
            "trd_minus_one_identity": t_minus_one.is_identical(&d_tr(module).unwrap()),
            "trd_two": attribute_module_json(&t_two, coefficient),
        })
    }

    fn transpose_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism.cached_transpose_homomorphism().is_some();
        let value = transpose_of_module_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism.cached_transpose_homomorphism().is_some();
        let repeated = transpose_of_module_homomorphism(homomorphism).unwrap();
        let source_attribute = transpose_of_module(homomorphism.range()).unwrap();
        let range_attribute = transpose_of_module(homomorphism.source()).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
            "source_attribute_identity": value.source().is_identical(&source_attribute),
            "range_attribute_identity": value.range().is_identical(&range_attribute),
        })
    }

    fn d_tr_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism.cached_d_tr_homomorphism().is_some();
        let value = d_tr_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism.cached_d_tr_homomorphism().is_some();
        let repeated = d_tr_homomorphism(homomorphism).unwrap();
        let source_attribute = d_tr(homomorphism.source()).unwrap();
        let range_attribute = d_tr(homomorphism.range()).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
            "source_attribute_identity": value.source().is_identical(&source_attribute),
            "range_attribute_identity": value.range().is_identical(&range_attribute),
        })
    }

    fn tr_d_map_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let has_before = homomorphism.cached_tr_d_homomorphism().is_some();
        let value = tr_d_homomorphism(homomorphism).unwrap();
        let has_after = homomorphism.cached_tr_d_homomorphism().is_some();
        let repeated = tr_d_homomorphism(homomorphism).unwrap();
        let source_attribute = tr_d(homomorphism.source()).unwrap();
        let range_attribute = tr_d(homomorphism.range()).unwrap();
        json!({
            "has_before": has_before,
            "value": module_map_json(&value, coefficient),
            "has_after": has_after,
            "repeat_identity": value.is_identical(&repeated),
            "source_attribute_identity": value.source().is_identical(&source_attribute),
            "range_attribute_identity": value.range().is_identical(&range_attribute),
        })
    }

    fn ar_map_functor_attribute_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "entered": module_map_json(homomorphism, coefficient),
            "transpose": transpose_map_attribute_json(homomorphism, coefficient),
            "dtr": d_tr_map_attribute_json(homomorphism, coefficient),
            "trd": tr_d_map_attribute_json(homomorphism, coefficient),
        })
    }

    fn ar_map_iterations_json<F: Field>(
        homomorphism: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let d_zero = d_tr_homomorphism_iterated(homomorphism, 0).unwrap();
        let d_one = d_tr_homomorphism_iterated(homomorphism, 1).unwrap();
        let d_minus_one = d_tr_homomorphism_iterated(homomorphism, -1).unwrap();
        let d_two = d_tr_homomorphism_iterated(homomorphism, 2).unwrap();
        let t_zero = tr_d_homomorphism_iterated(homomorphism, 0).unwrap();
        let t_one = tr_d_homomorphism_iterated(homomorphism, 1).unwrap();
        let t_minus_one = tr_d_homomorphism_iterated(homomorphism, -1).unwrap();
        let t_two = tr_d_homomorphism_iterated(homomorphism, 2).unwrap();
        json!({
            "dtr_zero_identity": d_zero.is_identical(homomorphism),
            "dtr_one_identity": d_one.is_identical(&d_tr_homomorphism(homomorphism).unwrap()),
            "dtr_minus_one_identity": d_minus_one
                .is_identical(&tr_d_homomorphism(homomorphism).unwrap()),
            "dtr_two": module_map_json(&d_two, coefficient),
            "trd_zero_identity": t_zero.is_identical(homomorphism),
            "trd_one_identity": t_one.is_identical(&tr_d_homomorphism(homomorphism).unwrap()),
            "trd_minus_one_identity": t_minus_one
                .is_identical(&d_tr_homomorphism(homomorphism).unwrap()),
            "trd_two": module_map_json(&t_two, coefficient),
        })
    }

    fn tensor_space_json<F: Field>(
        tensor: &TensorProductOfModulesResult<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "dimension": tensor.dimension(),
            "basis": tensor
                .basis_vectors()
                .iter()
                .map(|row| row.iter().map(coefficient).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
        })
    }

    fn caught_tensor_element_json<F: Field>(
        tensor: &TensorProductOfModulesResult<F>,
        left: &PathAlgebraModuleElement<F>,
        right: &PathAlgebraModuleElement<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        match tensor.elementary_tensor(left, right) {
            Ok(value) => json!({
                "success": true,
                "value": value.iter().map(coefficient).collect::<Vec<_>>(),
            }),
            Err(_) => json!({"success": false, "value": Value::Null}),
        }
    }

    fn property_bound<F: Field>(module: &PathAlgebraModule<F>) -> Vec<bool> {
        vec![
            module.cached_is_injective_module().is_some(),
            module.cached_is_simple_qpa_module().is_some(),
            module.cached_is_semisimple_module().is_some(),
        ]
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn module_families_properties_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let injectives = indec_injective_modules(&algebra).unwrap();
        let injectives_again = indec_injective_modules(&algebra).unwrap();
        let simples = simple_modules(&algebra).unwrap();
        let simples_again = simple_modules(&algebra).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let zero = zero_module(&algebra).unwrap();
        let semisimple_sum = direct_sum_of_qpa_modules(&simples).unwrap().unwrap();
        let injective_installed = injectives
            .iter()
            .map(|module| {
                vec![
                    module.cached_is_injective_module() == Some(true),
                    module.cached_is_indecomposable_module() == Some(true),
                    module.cached_is_projective_module().is_some(),
                ]
            })
            .collect::<Vec<_>>();
        let simple_installed = simples
            .iter()
            .map(|module| {
                vec![
                    module.cached_is_simple_qpa_module() == Some(true),
                    module.cached_is_indecomposable_module() == Some(true),
                    module.cached_is_semisimple_module() == Some(true),
                ]
            })
            .collect::<Vec<_>>();
        let injective_values = injectives
            .iter()
            .chain(&simples)
            .chain(&projectives)
            .chain(std::iter::once(&zero))
            .map(|module| is_injective_module(module).unwrap())
            .collect::<Vec<_>>();
        let simple_values = simples
            .iter()
            .chain(&projectives)
            .chain([&zero, &semisimple_sum])
            .map(is_simple_qpa_module)
            .collect::<Vec<_>>();
        let semisimple_values = simples
            .iter()
            .chain(&projectives)
            .chain([&zero, &semisimple_sum])
            .map(|module| is_semisimple_module(module).unwrap())
            .collect::<Vec<_>>();
        let properties_bound = injectives
            .iter()
            .chain(&simples)
            .chain(&projectives)
            .chain([&zero, &semisimple_sum])
            .map(property_bound)
            .collect::<Vec<_>>();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field, loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_injectives = indec_injective_modules(&quotient).unwrap();
        let quotient_injectives_again = indec_injective_modules(&quotient).unwrap();
        let quotient_simples = simple_modules(&quotient).unwrap();
        let quotient_simples_again = simple_modules(&quotient).unwrap();
        let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
        let quotient_zero = zero_module(&quotient).unwrap();
        let quotient_injective_installed = quotient_injectives
            .iter()
            .map(|module| {
                vec![
                    module.cached_is_injective_module() == Some(true),
                    module.cached_is_indecomposable_module() == Some(true),
                    module.cached_is_projective_module().is_some(),
                ]
            })
            .collect::<Vec<_>>();
        let quotient_simple_installed = quotient_simples
            .iter()
            .map(|module| {
                vec![
                    module.cached_is_simple_qpa_module() == Some(true),
                    module.cached_is_indecomposable_module() == Some(true),
                    module.cached_is_semisimple_module() == Some(true),
                ]
            })
            .collect::<Vec<_>>();
        let quotient_injective_values = quotient_injectives
            .iter()
            .chain(&quotient_simples)
            .chain([&quotient_projective, &quotient_zero])
            .map(|module| is_injective_module(module).unwrap())
            .collect::<Vec<_>>();
        let quotient_simple_values = quotient_simples
            .iter()
            .chain([&quotient_projective, &quotient_zero])
            .map(is_simple_qpa_module)
            .collect::<Vec<_>>();
        let quotient_semisimple_values = quotient_simples
            .iter()
            .chain([&quotient_projective, &quotient_zero])
            .map(|module| is_semisimple_module(module).unwrap())
            .collect::<Vec<_>>();
        let quotient_properties_bound = quotient_injectives
            .iter()
            .chain(&quotient_simples)
            .chain([&quotient_projective, &quotient_zero])
            .map(property_bound)
            .collect::<Vec<_>>();

        json!({
            "path": {
                "injectives": attribute_module_list_json(&injectives, coefficient),
                "injective_cache_identity": injectives.iter().zip(&injectives_again)
                    .map(|(first, second)| first.is_identical(second)).collect::<Vec<_>>(),
                "injective_installed": injective_installed,
                "simples": attribute_module_list_json(&simples, coefficient),
                "simple_cache_identity": simples.iter().zip(&simples_again)
                    .map(|(first, second)| first.is_identical(second)).collect::<Vec<_>>(),
                "simple_installed": simple_installed,
                "injective_values": injective_values,
                "simple_values": simple_values,
                "semisimple_values": semisimple_values,
                "property_bound": properties_bound,
            },
            "quotient": {
                "injectives": attribute_module_list_json(&quotient_injectives, coefficient),
                "injective_cache_identity": quotient_injectives.iter()
                    .zip(&quotient_injectives_again).map(|(first, second)|
                        first.is_identical(second)).collect::<Vec<_>>(),
                "injective_installed": quotient_injective_installed,
                "simples": attribute_module_list_json(&quotient_simples, coefficient),
                "simple_cache_identity": quotient_simples.iter().zip(&quotient_simples_again)
                    .map(|(first, second)| first.is_identical(second)).collect::<Vec<_>>(),
                "simple_installed": quotient_simple_installed,
                "injective_values": quotient_injective_values,
                "simple_values": quotient_simple_values,
                "semisimple_values": quotient_semisimple_values,
                "property_bound": quotient_properties_bound,
            },
        })
    }

    fn series_json(series: ModuleDimensionSeries) -> Value {
        match series {
            ModuleDimensionSeries::ZeroDimensionVector(vector) => json!(vector),
            ModuleDimensionSeries::Layers(layers) => json!(layers),
        }
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn module_series_case<F: Field>(field: F) -> Value {
        let quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let projectives = indec_projective_modules(&algebra).unwrap();
        let path_simple = simple_modules(&algebra).unwrap()[0].clone();
        let zero = zero_module(&algebra).unwrap();
        let semisimple = right_module_over_path_algebra(algebra, vec![1, 1, 1], &[])
            .unwrap()
            .unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field, loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let x_squared = x.multiply(&x).unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x_squared.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
        let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
        let quotient_zero = zero_module(&quotient).unwrap();

        json!({
            "path": {
                "dimensions": projectives.iter().map(PathAlgebraModule::dimension)
                    .chain([zero.dimension(), semisimple.dimension()]).collect::<Vec<_>>(),
                "partial_order": [
                    dimension_vector_partial_order(&projectives[2], &projectives[0]),
                    dimension_vector_partial_order(&projectives[0], &projectives[1]),
                    dimension_vector_partial_order(&projectives[0], &projectives[0]),
                    dimension_vector_partial_order(&projectives[1], &projectives[0]),
                ],
                "radical_series": projectives.iter().chain(std::iter::once(&semisimple))
                    .map(|module| series_json(radical_series(module).unwrap()))
                    .collect::<Vec<_>>(),
                "socle_series": projectives.iter().chain(std::iter::once(&semisimple))
                    .map(|module| series_json(socle_series(module).unwrap()))
                    .collect::<Vec<_>>(),
                "zero_radical_series": series_json(radical_series(&zero).unwrap()),
                "zero_socle_series": series_json(socle_series(&zero).unwrap()),
            },
            "quotient": {
                "dimensions": [quotient_projective.dimension(), quotient_simple.dimension(),
                    quotient_zero.dimension()],
                "partial_order": [
                    dimension_vector_partial_order(&quotient_simple, &quotient_projective),
                    dimension_vector_partial_order(&quotient_projective, &quotient_simple),
                ],
                "cross_length_partial_order": [
                    dimension_vector_partial_order(&quotient_simple, &path_simple),
                    dimension_vector_partial_order(&path_simple, &quotient_simple),
                    dimension_vector_partial_order(&quotient_simple, &projectives[0]),
                    dimension_vector_partial_order(&projectives[0], &quotient_simple),
                    dimension_vector_partial_order(&quotient_simple, &projectives[1]),
                    dimension_vector_partial_order(&projectives[1], &quotient_simple),
                ],
                "radical_projective": series_json(radical_series(&quotient_projective).unwrap()),
                "socle_projective": series_json(socle_series(&quotient_projective).unwrap()),
                "radical_simple": series_json(radical_series(&quotient_simple).unwrap()),
                "socle_simple": series_json(socle_series(&quotient_simple).unwrap()),
                "zero_radical_series": series_json(radical_series(&quotient_zero).unwrap()),
                "zero_socle_series": series_json(socle_series(&quotient_zero).unwrap()),
            },
        })
    }

    fn annihilator_json<F: Field>(
        algebra: &QuiverAlgebra<F>,
        elements: &[QuiverAlgebraElement<F>],
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let canonical_basis = finite_canonical_basis(algebra).unwrap().unwrap();
        let default_basis = finite_qpa_default_basis(algebra, &canonical_basis).unwrap();
        let basis_change = DenseMatrix::from_rows(
            default_basis
                .iter()
                .map(|element| {
                    canonical_basis_coefficients(algebra, &canonical_basis, element).unwrap()
                })
                .collect(),
            algebra.field().clone(),
        )
        .unwrap();
        json!(
            elements
                .iter()
                .map(|element| {
                    let canonical_coordinates =
                        canonical_basis_coefficients(algebra, &canonical_basis, element).unwrap();
                    basis_change
                        .solution_row(&canonical_coordinates)
                        .unwrap()
                        .unwrap()
                        .iter()
                        .map(coefficient)
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        )
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn module_annihilator_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let projectives = indec_projective_modules(&algebra).unwrap();
        let faithful = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
        let simples = simple_modules(&algebra).unwrap();
        let zero = zero_module(&algebra).unwrap();
        let arrow_zero = right_module_over_path_algebra(algebra.clone(), vec![1, 1], &[])
            .unwrap()
            .unwrap();
        let zero_annihilator = annihilator_of_module(&zero).unwrap();
        let simple1_annihilator = annihilator_of_module(&simples[0]).unwrap();
        let simple2_annihilator = annihilator_of_module(&simples[1]).unwrap();
        let arrow_zero_annihilator = annihilator_of_module(&arrow_zero).unwrap();
        let faithful_annihilator = annihilator_of_module(&faithful).unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
        let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
        let quotient_zero = zero_module(&quotient).unwrap();
        let quotient_two_zero = right_module_over_path_algebra(quotient.clone(), vec![2], &[])
            .unwrap()
            .unwrap();
        let quotient_projective_annihilator = annihilator_of_module(&quotient_projective).unwrap();
        let quotient_simple_annihilator = annihilator_of_module(&quotient_simple).unwrap();
        let quotient_zero_annihilator = annihilator_of_module(&quotient_zero).unwrap();
        let quotient_two_zero_annihilator = annihilator_of_module(&quotient_two_zero).unwrap();

        let bound_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "c"),
                ArrowSpec::named(vertex(2), vertex(3), "d"),
            ],
        )
        .unwrap();
        let bound_parent = PathAlgebra::new(field, bound_quiver).unwrap();
        let bound_arrows = bound_parent.quiver().arrows().collect::<Vec<_>>();
        let c = bound_parent.element(bound_arrows[0].path()).unwrap();
        let d = bound_parent.element(bound_arrows[1].path()).unwrap();
        let bound_quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(bound_parent, &[c.multiply(&d).unwrap()]).unwrap(),
        );
        let bound_zero = zero_module(&bound_quotient).unwrap();
        let bound_simple = simple_modules(&bound_quotient).unwrap()[2].clone();
        let bound_zero_annihilator = annihilator_of_module(&bound_zero).unwrap();
        let bound_simple_annihilator = annihilator_of_module(&bound_simple).unwrap();

        json!({
            "path": {
                "basis_dimension": finite_canonical_basis(&algebra).unwrap().unwrap().len(),
                "zero": annihilator_json(&algebra, &zero_annihilator, coefficient),
                "simple1": annihilator_json(&algebra, &simple1_annihilator, coefficient),
                "simple2": annihilator_json(&algebra, &simple2_annihilator, coefficient),
                "arrow_zero": annihilator_json(&algebra, &arrow_zero_annihilator, coefficient),
                "faithful": annihilator_json(&algebra, &faithful_annihilator, coefficient),
                "cache_bound": [
                    zero.cached_annihilator_of_module().is_some(),
                    simples[0].cached_annihilator_of_module().is_some(),
                    simples[1].cached_annihilator_of_module().is_some(),
                    arrow_zero.cached_annihilator_of_module().is_some(),
                    faithful.cached_annihilator_of_module().is_some(),
                ],
            },
            "quotient": {
                "basis_dimension": finite_canonical_basis(&quotient).unwrap().unwrap().len(),
                "projective": annihilator_json(
                    &quotient,
                    &quotient_projective_annihilator,
                    coefficient,
                ),
                "simple": annihilator_json(
                    &quotient,
                    &quotient_simple_annihilator,
                    coefficient,
                ),
                "zero": annihilator_json(&quotient, &quotient_zero_annihilator, coefficient),
                "two_zero": annihilator_json(
                    &quotient,
                    &quotient_two_zero_annihilator,
                    coefficient,
                ),
                "cache_bound": [
                    quotient_projective.cached_annihilator_of_module().is_some(),
                    quotient_simple.cached_annihilator_of_module().is_some(),
                    quotient_zero.cached_annihilator_of_module().is_some(),
                    quotient_two_zero.cached_annihilator_of_module().is_some(),
                ],
            },
            "bound_quotient": {
                "basis_dimension": finite_canonical_basis(&bound_quotient)
                    .unwrap().unwrap().len(),
                "zero": annihilator_json(
                    &bound_quotient,
                    &bound_zero_annihilator,
                    coefficient,
                ),
                "simple3": annihilator_json(
                    &bound_quotient,
                    &bound_simple_annihilator,
                    coefficient,
                ),
                "cache_bound": [
                    bound_zero.cached_annihilator_of_module().is_some(),
                    bound_simple.cached_annihilator_of_module().is_some(),
                ],
            },
        })
    }

    fn ext_json<F: Field>(
        ext: &ExtOverAlgebraResult<F>,
        field: &F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let basis_coordinates = if ext.basis().is_empty() {
            Value::Null
        } else {
            json!(
                ext.basis()
                    .iter()
                    .map(|homomorphism| ext
                        .coefficients(homomorphism)
                        .unwrap()
                        .unwrap()
                        .iter()
                        .map(coefficient)
                        .collect::<Vec<_>>())
                    .collect::<Vec<_>>()
            )
        };
        let (entered_coefficients, combo_coordinates) = if ext.basis().is_empty() {
            (Value::Null, Value::Null)
        } else {
            let entered = (0..ext.dimension())
                .map(|index| field.nth(Integer::from(i64::try_from(index + 2).unwrap())))
                .collect::<Vec<_>>();
            let mut terms = ext
                .basis()
                .iter()
                .zip(&entered)
                .map(|(homomorphism, scalar)| homomorphism.scale(scalar).unwrap());
            let mut combo = terms.next().unwrap();
            for term in terms {
                combo = combo.add(&term).unwrap();
            }
            let coordinates = ext.coefficients(&combo).unwrap().unwrap();
            (
                json!(entered.iter().map(coefficient).collect::<Vec<_>>()),
                json!(coordinates.iter().map(coefficient).collect::<Vec<_>>()),
            )
        };
        json!({
            "syzygy_inclusion": module_map_json(ext.syzygy_inclusion(), coefficient),
            "dimension": ext.dimension(),
            "basis": ext.basis().iter()
                .map(|homomorphism| module_map_json(homomorphism, coefficient))
                .collect::<Vec<_>>(),
            "basis_coordinates": basis_coordinates,
            "entered_coefficients": entered_coefficients,
            "combo_coordinates": combo_coordinates,
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn functors_star_nakayama_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let simples = simple_modules(&algebra).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let zero = zero_module(&algebra).unwrap();
        let cover = projective_cover(&simples[0]).unwrap();
        let inclusion = kernel_inclusion(&cover).unwrap();

        let mixed = right_module_over_path_algebra(
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
        let mixed_map = right_module_hom_over_algebra(
            mixed.clone(),
            mixed.clone(),
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![field.one(), field.zero()],
                        vec![field.zero(), field.zero()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::from_rows(vec![vec![field.one()]], field.clone()).unwrap(),
            ],
        )
        .unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
        let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
        let quotient_x = include_in_path_algebra(
            &quotient.quiver().arrows().next().unwrap().path(),
            &quotient,
        )
        .unwrap();
        let quotient_generator =
            minimal_generating_set_of_module(quotient_projective.clone()).unwrap()[0].clone();
        let radical_image = quotient_generator.act(&quotient_x).unwrap();
        let quotient_radical_map = hom_from_projective(&radical_image, &quotient_projective)
            .unwrap()
            .unwrap();

        json!({
            "modules": {
                "a2_s2": module_functor_attribute_json(&simples[1], coefficient),
                "a2_p1": module_functor_attribute_json(&projectives[0], coefficient),
                "a2_zero": module_functor_attribute_json(&zero, coefficient),
                "a2_mixed": module_functor_attribute_json(&mixed, coefficient),
                "dual_simple": module_functor_attribute_json(&quotient_simple, coefficient),
                "dual_projective": module_functor_attribute_json(
                    &quotient_projective,
                    coefficient,
                ),
            },
            "maps": {
                "a2_inclusion": map_functor_attribute_json(&inclusion, coefficient),
                "a2_endomorphism": map_functor_attribute_json(&mixed_map, coefficient),
                "dual_radical": map_functor_attribute_json(&quotient_radical_map, coefficient),
            },
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn functors_transpose_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let simples = simple_modules(&algebra).unwrap();
        let mixed = right_module_over_path_algebra(
            algebra,
            vec![2, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![field.one()], vec![field.zero()]], field.clone())
                    .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        let kronecker_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "p"),
                ArrowSpec::named(vertex(1), vertex(2), "q"),
            ],
        )
        .unwrap();
        let kronecker =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
        let kronecker_simple = simple_modules(&kronecker).unwrap()[0].clone();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();

        let bound_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let bound_parent = PathAlgebra::new(field, bound_quiver).unwrap();
        let bound_arrows = bound_parent.quiver().arrows().collect::<Vec<_>>();
        let a = bound_parent.element(bound_arrows[0].path()).unwrap();
        let b = bound_parent.element(bound_arrows[1].path()).unwrap();
        let bound = QuiverAlgebra::Quotient(
            quotient_by_relations(bound_parent, &[a.multiply(&b).unwrap()]).unwrap(),
        );
        let bound_simple = simple_modules(&bound).unwrap()[0].clone();

        json!({
            "a2_s1": transpose_functor_attribute_json(&simples[0], coefficient),
            "a2_s2_projective": transpose_functor_attribute_json(&simples[1], coefficient),
            "a2_mixed": transpose_functor_attribute_json(&mixed, coefficient),
            "kronecker_s1": transpose_functor_attribute_json(&kronecker_simple, coefficient),
            "dual_simple": transpose_functor_attribute_json(&dual_simple, coefficient),
            "bound_a3_s1": transpose_functor_attribute_json(&bound_simple, coefficient),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn functors_ar_maps_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let a2_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());
        let a2_module = right_module_over_path_algebra(a2, vec![2, 0], &[])
            .unwrap()
            .unwrap();
        let a2_map = right_module_hom_over_algebra(
            a2_module.clone(),
            a2_module,
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![field.one(), field.one()],
                        vec![field.zero(), field.one()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::zero(1, 1, field.clone()),
            ],
        )
        .unwrap();

        let kronecker_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "p"),
                ArrowSpec::named(vertex(1), vertex(2), "q"),
            ],
        )
        .unwrap();
        let kronecker =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
        let kronecker_module = right_module_over_path_algebra(kronecker, vec![2, 0], &[])
            .unwrap()
            .unwrap();
        let kronecker_map = right_module_hom_over_algebra(
            kronecker_module.clone(),
            kronecker_module,
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![field.one(), field.one()],
                        vec![field.zero(), field.one()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::zero(1, 1, field.clone()),
            ],
        )
        .unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();
        let dual_identity = identity_module_mapping(dual_simple).unwrap();

        let bound_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let bound_parent = PathAlgebra::new(field, bound_quiver).unwrap();
        let bound_arrows = bound_parent.quiver().arrows().collect::<Vec<_>>();
        let a = bound_parent.element(bound_arrows[0].path()).unwrap();
        let b = bound_parent.element(bound_arrows[1].path()).unwrap();
        let bound = QuiverAlgebra::Quotient(
            quotient_by_relations(bound_parent, &[a.multiply(&b).unwrap()]).unwrap(),
        );
        let bound_simple = simple_modules(&bound).unwrap()[0].clone();
        let bound_identity = identity_module_mapping(bound_simple).unwrap();

        json!({
            "a2_shear": ar_map_functor_attribute_json(&a2_map, coefficient),
            "kronecker_shear": ar_map_functor_attribute_json(&kronecker_map, coefficient),
            "dual_simple_identity": ar_map_functor_attribute_json(&dual_identity, coefficient),
            "bound_a3_identity": ar_map_functor_attribute_json(&bound_identity, coefficient),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn functors_ar_map_iterations_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let a2_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());
        let a2_module = right_module_over_path_algebra(a2, vec![2, 0], &[])
            .unwrap()
            .unwrap();
        let a2_map = right_module_hom_over_algebra(
            a2_module.clone(),
            a2_module,
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![field.one(), field.one()],
                        vec![field.zero(), field.one()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::zero(1, 1, field.clone()),
            ],
        )
        .unwrap();

        let kronecker_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "p"),
                ArrowSpec::named(vertex(1), vertex(2), "q"),
            ],
        )
        .unwrap();
        let kronecker =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
        let kronecker_module = right_module_over_path_algebra(kronecker, vec![2, 0], &[])
            .unwrap()
            .unwrap();
        let kronecker_map = right_module_hom_over_algebra(
            kronecker_module.clone(),
            kronecker_module,
            vec![
                DenseMatrix::from_rows(
                    vec![
                        vec![field.one(), field.one()],
                        vec![field.zero(), field.one()],
                    ],
                    field.clone(),
                )
                .unwrap(),
                DenseMatrix::zero(1, 1, field.clone()),
            ],
        )
        .unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_parent = PathAlgebra::new(field, loop_quiver).unwrap();
        let x = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();
        let dual_identity = identity_module_mapping(dual_simple).unwrap();

        json!({
            "a2_shear": ar_map_iterations_json(&a2_map, coefficient),
            "kronecker_shear": ar_map_iterations_json(&kronecker_map, coefficient),
            "dual_simple_identity": ar_map_iterations_json(&dual_identity, coefficient),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn functors_tensor_modules_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let a2_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());
        let simples = simple_modules(&a2).unwrap();
        let projectives = indec_projective_modules(&a2).unwrap();

        let zero_right = dual_of_module(&simples[1]).unwrap();
        let zero_tensor = tensor_product_of_modules(&simples[0], &zero_right).unwrap();
        let zero_element = zero_tensor
            .elementary_tensor(
                &simples[0].basis_elements()[0],
                &zero_right.basis_elements()[0],
            )
            .unwrap();

        let projective_right = dual_of_module(&simples[0]).unwrap();
        let projective_tensor =
            tensor_product_of_modules(&projectives[0], &projective_right).unwrap();
        let projective_element = caught_tensor_element_json(
            &projective_tensor,
            &projectives[0].basis_elements()[0],
            &projective_right.basis_elements()[0],
            coefficient,
        );
        assert_eq!(
            projective_tensor
                .elementary_tensor(
                    &projectives[0].basis_elements()[0],
                    &projective_right.basis_elements()[0],
                )
                .unwrap_err(),
            super::EnvelopingModuleError::TensorProductUnboundSizeMatrix
        );

        let path_main =
            tensor_product_of_modules(&simples[0], &dual_of_module(&projectives[0]).unwrap());
        assert_eq!(
            path_main.as_ref().unwrap_err(),
            &super::EnvelopingModuleError::TensorProductPathAlgebraP1Extraction
        );
        let wrong_algebras = tensor_product_of_modules(&simples[0], &simples[0]);
        assert_eq!(
            wrong_algebras.as_ref().unwrap_err(),
            &super::EnvelopingModuleError::TensorProductModulesOverWrongAlgebras
        );

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let cubic = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(3).unwrap()]).unwrap(),
        );
        let cubic_module = right_module_over_path_algebra(
            cubic.clone(),
            vec![2],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(
                    vec![
                        vec![field.zero(), field.one()],
                        vec![field.zero(), field.zero()],
                    ],
                    field.clone(),
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let cubic_projective = indec_projective_modules(&cubic).unwrap()[0].clone();
        let cubic_right = dual_of_module(&cubic_projective).unwrap();
        let cubic_tensor = tensor_product_of_modules(&cubic_module, &cubic_right).unwrap();
        let cubic_elements = cubic_module
            .basis_elements()
            .iter()
            .map(|left| {
                cubic_right
                    .basis_elements()
                    .iter()
                    .map(|right| {
                        caught_tensor_element_json(&cubic_tensor, left, right, coefficient)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        let bound_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let bound_parent = PathAlgebra::new(field, bound_quiver).unwrap();
        let bound_arrows = bound_parent.quiver().arrows().collect::<Vec<_>>();
        let a = bound_parent.element(bound_arrows[0].path()).unwrap();
        let b = bound_parent.element(bound_arrows[1].path()).unwrap();
        let bound = QuiverAlgebra::Quotient(
            quotient_by_relations(bound_parent, &[a.multiply(&b).unwrap()]).unwrap(),
        );
        let bound_simple = simple_modules(&bound).unwrap()[0].clone();
        let bound_projective = indec_projective_modules(&bound).unwrap()[0].clone();
        let bound_right = dual_of_module(&bound_projective).unwrap();
        let bound_tensor = tensor_product_of_modules(&bound_simple, &bound_right).unwrap();
        let bound_element = caught_tensor_element_json(
            &bound_tensor,
            &bound_simple.basis_elements()[0],
            &bound_right.basis_elements()[0],
            coefficient,
        );
        assert!(matches!(
            bound_tensor.elementary_tensor(
                &bound_simple.basis_elements()[0],
                &bound_right.basis_elements()[0]
            ),
            Err(
                super::EnvelopingModuleError::TensorProductClosureProjectionOutOfRange {
                    index: 2,
                    count: 1
                }
            )
        ));

        json!({
            "easy_zero": {
                "space": tensor_space_json(&zero_tensor, coefficient),
                "element": zero_element.iter().map(coefficient).collect::<Vec<_>>(),
            },
            "projective_branch": {
                "space": tensor_space_json(&projective_tensor, coefficient),
                "element": projective_element,
            },
            "path_main_constructor_success": path_main.is_ok(),
            "wrong_algebra_constructor_success": wrong_algebras.is_ok(),
            "cubic_quotient": {
                "space": tensor_space_json(&cubic_tensor, coefficient),
                "elements": cubic_elements,
            },
            "bound_capture_bug": {
                "space": tensor_space_json(&bound_tensor, coefficient),
                "element": bound_element,
            },
        })
    }

    #[allow(clippy::needless_pass_by_value)]
    fn functors_iterations_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let a2_simple = simple_modules(&algebra).unwrap()[0].clone();

        let kronecker_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "p"),
                ArrowSpec::named(vertex(1), vertex(2), "q"),
            ],
        )
        .unwrap();
        let kronecker =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
        let kronecker_simple = simple_modules(&kronecker).unwrap()[0].clone();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field, loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let dual_simple = simple_modules(&dual_numbers).unwrap()[0].clone();

        json!({
            "a2_s1": functor_iterations_json(&a2_simple, coefficient),
            "kronecker_s1": functor_iterations_json(&kronecker_simple, coefficient),
            "dual_simple": functor_iterations_json(&dual_simple, coefficient),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_ext_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let simples = simple_modules(&algebra).unwrap();
        let zero = zero_module(&algebra).unwrap();
        let a2_12 = ext_over_algebra(&simples[0], &simples[1]).unwrap();
        let a2_11 = ext_over_algebra(&simples[0], &simples[0]).unwrap();
        let a2_22 = ext_over_algebra(&simples[1], &simples[1]).unwrap();
        let a2_02 = ext_over_algebra(&zero, &simples[1]).unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_simple = simple_modules(&quotient).unwrap()[0].clone();
        let quotient_projective = indec_projective_modules(&quotient).unwrap()[0].clone();
        let dual_simple_ext = ext_over_algebra(&quotient_simple, &quotient_simple).unwrap();
        let dual_projective_ext = ext_over_algebra(&quotient_simple, &quotient_projective).unwrap();
        let quotient_sum =
            direct_sum_of_qpa_modules(&[quotient_simple.clone(), quotient_projective.clone()])
                .unwrap()
                .unwrap();
        let dual_sum = ext_over_algebra(&quotient_simple, &quotient_sum).unwrap();

        let kronecker_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "p"),
                ArrowSpec::named(vertex(1), vertex(2), "q"),
            ],
        )
        .unwrap();
        let kronecker =
            QuiverAlgebra::Path(PathAlgebra::new(field.clone(), kronecker_quiver).unwrap());
        let kronecker_simples = simple_modules(&kronecker).unwrap();
        let kronecker_ext = ext_over_algebra(&kronecker_simples[0], &kronecker_simples[1]).unwrap();

        let rigid_json = |module: &PathAlgebraModule<F>| {
            let has_before = module.cached_is_rigid_module().is_some();
            let value = is_rigid_module(module).unwrap();
            let has_after = module.cached_is_rigid_module().is_some();
            json!({
                "has_before": has_before,
                "value": value,
                "has_after": has_after,
            })
        };

        json!({
            "a2_s1_s2": ext_json(&a2_12, &field, coefficient),
            "a2_s1_s1": ext_json(&a2_11, &field, coefficient),
            "a2_s2_s2": ext_json(&a2_22, &field, coefficient),
            "a2_zero_s2": ext_json(&a2_02, &field, coefficient),
            "dual_simple_simple": ext_json(&dual_simple_ext, &field, coefficient),
            "dual_simple_projective": ext_json(&dual_projective_ext, &field, coefficient),
            "dual_simple_sum": ext_json(&dual_sum, &field, coefficient),
            "kronecker_s1_s2": ext_json(&kronecker_ext, &field, coefficient),
            "rigid": {
                "a2_s1": rigid_json(&simples[0]),
                "a2_s2": rigid_json(&simples[1]),
                "a2_zero": rigid_json(&zero),
                "dual_simple": rigid_json(&quotient_simple),
                "dual_projective": rigid_json(&quotient_projective),
                "dual_sum": rigid_json(&quotient_sum),
                "kronecker_s1": rigid_json(&kronecker_simples[0]),
                "kronecker_s2": rigid_json(&kronecker_simples[1]),
            },
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn modulehom_syzygy_injective_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let zero = zero_module(&algebra).unwrap();
        let zero_again = zero_module(&algebra).unwrap();
        let simple1 = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple2 = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let syzygy1 = first_syzygy(&simple1).unwrap();
        let syzygy2 = first_syzygy(&simple2).unwrap();
        let syzygy0 = first_syzygy(&zero).unwrap();
        let envelope1 = injective_envelope(&simple1).unwrap();
        let envelope2 = injective_envelope(&simple2).unwrap();
        let envelope0 = injective_envelope(&zero).unwrap();

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(field, loop_quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let quotient_zero = zero_module(&quotient).unwrap();
        let quotient_simple = right_module_over_path_algebra(quotient.clone(), vec![1], &[])
            .unwrap()
            .unwrap();
        let quotient_syzygy = first_syzygy(&quotient_simple).unwrap();
        let quotient_envelope = injective_envelope(&quotient_simple).unwrap();

        json!({
            "path": {
                "zero": attribute_module_json(&zero, coefficient),
                "zero_cache_identity": zero.is_identical(&zero_again),
                "syzygy_simple1": attribute_module_json(&syzygy1, coefficient),
                "syzygy_simple2": attribute_module_json(&syzygy2, coefficient),
                "syzygy_zero": attribute_module_json(&syzygy0, coefficient),
                "syzygy_cache_identity": [
                    syzygy1.is_identical(&first_syzygy(&simple1).unwrap()),
                    syzygy2.is_identical(&first_syzygy(&simple2).unwrap()),
                    syzygy0.is_identical(&zero),
                ],
                "injective_simple1": module_map_json(&envelope1, coefficient),
                "injective_simple2": module_map_json(&envelope2, coefficient),
                "injective_zero": module_map_json(&envelope0, coefficient),
                "injective_cache_identity": [
                    envelope1.is_identical(&injective_envelope(&simple1).unwrap()),
                    envelope2.is_identical(&injective_envelope(&simple2).unwrap()),
                    envelope0.is_identical(&injective_envelope(&zero).unwrap()),
                ],
                "injective_source_identity": [
                    envelope1.source().is_identical(&simple1),
                    envelope2.source().is_identical(&simple2),
                    envelope0.source().is_identical(&zero),
                ],
            },
            "quotient": {
                "zero": attribute_module_json(&quotient_zero, coefficient),
                "syzygy_simple": attribute_module_json(&quotient_syzygy, coefficient),
                "injective_simple": module_map_json(&quotient_envelope, coefficient),
                "cache_identity": [
                    quotient_zero.is_identical(&zero_module(&quotient).unwrap()),
                    quotient_syzygy.is_identical(&first_syzygy(&quotient_simple).unwrap()),
                    quotient_envelope
                        .is_identical(&injective_envelope(&quotient_simple).unwrap()),
                ],
            },
        })
    }

    fn stable_module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        print_zero_vertex_map: bool,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let natural = end_modulo_projectives(module).unwrap();
        let factoring = hom_factoring_through_projectives(module, module).unwrap();
        let factor_matrices = factoring
            .iter()
            .map(|map| from_hom_mm_to_end_m(map).unwrap())
            .collect::<Vec<_>>();
        let quotient_basis = natural.range().basis_elements();
        let images = natural
            .source()
            .basis_matrices()
            .iter()
            .map(|matrix| {
                natural
                    .image(matrix)
                    .unwrap()
                    .coordinates()
                    .iter()
                    .map(coefficient)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let one = natural.range().one().unwrap().map(|element| {
            element
                .coordinates()
                .iter()
                .map(coefficient)
                .collect::<Vec<_>>()
        });
        let multiplication = quotient_basis
            .iter()
            .map(|left| {
                quotient_basis
                    .iter()
                    .map(|right| {
                        natural
                            .range()
                            .multiply(left, right)
                            .unwrap()
                            .coordinates()
                            .iter()
                            .map(coefficient)
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let roundtrip_matrices = natural
            .source()
            .basis_matrices()
            .iter()
            .map(|matrix| {
                from_hom_mm_to_end_m(&from_end_m_to_hom_mm(module, matrix).unwrap()).unwrap()
            })
            .collect::<Vec<_>>();
        let mut value = json!({
            "end_basis": natural.source().basis_matrices().iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
            "factor_matrices": factor_matrices.iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
            "ideal_basis": natural.range().ideal_basis_matrices().iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
            "kernel_basis": natural.kernel_basis_matrices().iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
            "source_cache_identity": natural.source().is_identical(&end_over_algebra(module).unwrap()),
            "dimension": natural.range().dimension(),
            "images": images,
            "one": one,
            "multiplication": multiplication,
            "roundtrip_matrices": roundtrip_matrices.iter()
                .map(|matrix| matrix_json(matrix, coefficient)).collect::<Vec<_>>(),
        });
        if print_zero_vertex_map {
            value["zero_vertex_map"] = module_map_json(
                &from_end_m_to_hom_mm(module, &natural.source().basis_matrices()[0]).unwrap(),
                coefficient,
            );
        }
        value
    }

    #[allow(clippy::needless_pass_by_value)]
    fn modulehom_stable_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), quiver).unwrap());
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 2],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(
                    vec![
                        vec![one.clone(), zero.clone()],
                        vec![zero.clone(), zero.clone()],
                    ],
                    field.clone(),
                )
                .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let simple_one = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(algebra, vec![0, 1], &[])
            .unwrap()
            .unwrap();

        let loop_quiver =
            Quiver::from_named(vec!["w".into()], vec![NamedArrowSpec::named("w", "w", "x")])
                .unwrap();
        let path_algebra = PathAlgebra::new(field, loop_quiver).unwrap();
        let loop_element = path_algebra
            .element(path_algebra.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(
                path_algebra,
                &[loop_element.multiply(&loop_element).unwrap()],
            )
            .unwrap(),
        );
        let dual_simple = right_module_over_path_algebra(dual_numbers.clone(), vec![1], &[])
            .unwrap()
            .unwrap();
        let dual_projective = indec_projective_modules(&dual_numbers).unwrap()[0].clone();

        json!({
            "a2_module": stable_module_json(&module, false, coefficient),
            "a2_simple_one": stable_module_json(&simple_one, true, coefficient),
            "a2_simple_two": stable_module_json(&simple_two, true, coefficient),
            "dual_simple": stable_module_json(&dual_simple, false, coefficient),
            "dual_projective": stable_module_json(&dual_projective, false, coefficient),
        })
    }

    #[test]
    fn pinned_qpa_modulehom_projective_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_projectives.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_projective_case(Q, ToString::to_string);
        let binary = modulehom_projective_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_modulehom_projective_matrix_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_projective_matrices.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_projective_matrix_case(Q, ToString::to_string);
        let binary = modulehom_projective_matrix_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_modulehom_syzygy_injective_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_syzygy_injective.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_syzygy_injective_case(Q, ToString::to_string);
        let binary = modulehom_syzygy_injective_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_module_families_properties_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/module_families_properties.schema-v1.json"
        ))
        .unwrap();
        let rational = module_families_properties_case(Q, ToString::to_string);
        let binary = module_families_properties_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_module_series_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/module_series.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(module_series_case(Q), expected["rational"]);
        assert_eq!(module_series_case(Z2), expected["binary"]);
    }

    #[test]
    fn pinned_qpa_module_annihilator_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/module_annihilator.schema-v1.json"
        ))
        .unwrap();
        let rational = module_annihilator_case(Q, ToString::to_string);
        let binary = module_annihilator_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_modulehom_ext_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_ext.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_ext_case(Q, ToString::to_string);
        let binary = modulehom_ext_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_star_and_nakayama_functor_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/functors_star_nakayama.schema-v1.json"
        ))
        .unwrap();
        let rational = functors_star_nakayama_case(Q, ToString::to_string);
        let binary = functors_star_nakayama_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_transpose_and_ar_functor_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/functors_transpose.schema-v1.json"
        ))
        .unwrap();
        let rational = functors_transpose_case(Q, ToString::to_string);
        let binary = functors_transpose_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_ar_functors_on_homomorphisms_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/functors_ar_maps.schema-v1.json"
        ))
        .unwrap();
        let rational = functors_ar_maps_case(Q, ToString::to_string);
        let binary = functors_ar_maps_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_ar_functor_homomorphism_iterations_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/functors_ar_map_iterations.schema-v1.json"
        ))
        .unwrap();
        let rational = functors_ar_map_iterations_case(Q, ToString::to_string);
        let binary = functors_ar_map_iterations_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_tensor_product_of_modules_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/functors_tensor_modules.schema-v1.json"
        ))
        .unwrap();
        let rational = functors_tensor_modules_case(Q, ToString::to_string);
        let binary = functors_tensor_modules_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_ar_functor_iterations_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/functors_iterations.schema-v1.json"
        ))
        .unwrap();
        let rational = functors_iterations_case(Q, ToString::to_string);
        let binary = functors_iterations_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn annihilator_of_module_reports_its_finite_basis_domain() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let module = right_module_over_path_algebra(algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        assert_eq!(
            annihilator_of_module(&module).unwrap_err(),
            super::EnvelopingModuleError::InfiniteAnnihilatorAlgebra
        );
    }

    #[test]
    fn ext_over_algebra_preserves_domain_and_coefficient_boundaries() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let simples = simple_modules(&algebra).unwrap();
        let ext = ext_over_algebra(&simples[0], &simples[1]).unwrap();
        let outside = identity_module_mapping(simples[0].clone()).unwrap();
        assert_eq!(
            ext.coefficients(&outside).unwrap_err(),
            super::EnvelopingModuleError::ExtMapOutsideHomSpace
        );
        let zero_ext = ext_over_algebra(&simples[0], &simples[0]).unwrap();
        assert_eq!(zero_ext.coefficients(&outside).unwrap(), None);

        let other_quiver = Quiver::from_indexed(1, Vec::new()).unwrap();
        let other_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, other_quiver).unwrap());
        let other_simple = simple_modules(&other_algebra).unwrap()[0].clone();
        assert_eq!(
            ext_over_algebra(&simples[0], &other_simple).unwrap_err(),
            super::EnvelopingModuleError::ExtModulesOverDifferentAlgebras
        );
    }

    #[test]
    fn simple_modules_preserves_qpa_source_domain_errors() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(Q, quiver).unwrap();
        assert_eq!(
            simple_modules(&QuiverAlgebra::Path(parent.clone())).unwrap_err(),
            super::EnvelopingModuleError::InfiniteSimpleModuleAlgebra
        );
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let quotient = QuiverAlgebra::Quotient(quotient_by_relations(parent, &[x]).unwrap());
        assert_eq!(
            simple_modules(&quotient).unwrap_err(),
            super::EnvelopingModuleError::NonadmissibleSimpleModuleAlgebra
        );
    }

    #[test]
    fn projective_matrix_conversions_report_qpa_domain_errors() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        assert_eq!(
            element_in_v_a_as_element_in_indec_proj(&algebra, &algebra.zero()),
            Err(super::EnvelopingModuleError::ZeroProjectiveAlgebraElement)
        );
        let v1 = include_in_path_algebra(
            &algebra.quiver().vertices().next().unwrap().path(),
            &algebra,
        )
        .unwrap();
        assert_eq!(
            from_matrix_to_homomorphism_of_projectives(&algebra, &[vec![v1]], &[2], &[1],),
            Err(super::EnvelopingModuleError::ProjectiveMatrixSupportDiscrepancy)
        );
        let scaled_arrow =
            include_in_path_algebra(&algebra.quiver().arrows().next().unwrap().path(), &algebra)
                .unwrap()
                .scale_left(&q(2))
                .unwrap();
        assert_eq!(
            element_in_v_a_as_element_in_indec_proj(&algebra, &scaled_arrow),
            Err(super::EnvelopingModuleError::ProjectiveLeadingTermMissing {
                left_vertex: 1,
                right_vertex: 2,
            })
        );
        assert_eq!(
            from_matrix_to_homomorphism_of_projectives(&algebra, &[], &[], &[]),
            Err(super::EnvelopingModuleError::EmptyProjectiveElementMatrix)
        );
    }

    #[test]
    fn pinned_qpa_modulehom_stable_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_stable.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_stable_case(Q, ToString::to_string);
        let binary = modulehom_stable_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn stable_endomorphism_quotient_rejects_the_zero_module_like_qpa() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let zero = right_module_over_path_algebra(algebra, vec![0, 0], &[])
            .unwrap()
            .unwrap();

        assert_eq!(
            end_modulo_projectives(&zero).unwrap_err(),
            super::EnvelopingModuleError::ZeroStableEndomorphismModule
        );
    }

    #[test]
    fn pinned_qpa_modulehom_functor_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_functors.schema-v1.json"
        ))
        .unwrap();
        let rational = modulehom_functor_case(Q, ToString::to_string);
        let binary = modulehom_functor_case(Z2, |entry| {
            if Z2.is_zero(entry) {
                "0*Z(2)".to_owned()
            } else {
                "Z(2)^0".to_owned()
            }
        });
        assert_eq!(rational, expected["rational"]);
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_regular_enveloping_module_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalgtensor_regular_module.schema-v1.json"
        ))
        .unwrap();

        let loop_quiver =
            Quiver::from_named(vec!["v".into()], vec![NamedArrowSpec::named("v", "v", "x")])
                .unwrap();
        let parent = PathAlgebra::new(Q, loop_quiver.clone()).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let rational_algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let rational_cached_before = rational_algebra.has_regular_enveloping_module_attribute();
        let rational_module = algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let rational_again = algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let rational_synonym = algebra_as_module_of_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let rational_enveloping = enveloping_algebra(&rational_algebra).unwrap();
        let rational_json = json!({
            "attribute_identity": [
                rational_cached_before,
                rational_algebra.has_regular_enveloping_module_attribute(),
                rational_module.is_identical(&rational_again),
                rational_module.is_identical(&rational_synonym),
                rational_module.right_acting_algebra() == rational_enveloping.algebra(),
            ],
            "module": module_json(&rational_module, ToString::to_string),
        });

        let path_quiver = Quiver::from_named(
            vec!["u1".into(), "u2".into()],
            vec![NamedArrowSpec::named("u1", "u2", "a")],
        )
        .unwrap();
        let binary_algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, path_quiver).unwrap());
        let binary_module = algebra_as_module_over_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_again = algebra_as_module_of_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_enveloping = enveloping_algebra(&binary_algebra).unwrap();
        let binary_json = json!({
            "attribute_identity": [
                binary_module.is_identical(&binary_again),
                binary_module.right_acting_algebra() == binary_enveloping.algebra(),
            ],
            "module": module_json(&binary_module, |entry| {
                if Z2.is_zero(entry) {
                    "0*Z(2)".to_owned()
                } else {
                    "Z(2)^0".to_owned()
                }
            }),
        });

        let infinite_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, loop_quiver).unwrap());
        let infinite_result = algebra_as_module_over_enveloping_algebra(&infinite_algebra).unwrap();
        let infinite_json = json!({
            "is_fail": infinite_result.is_none(),
            "attribute_bound": infinite_algebra.has_regular_enveloping_module_attribute(),
            "enveloping_bound": infinite_algebra.cached_enveloping_algebra().is_some(),
        });

        assert_eq!(rational_json, expected["rational_quotient"]);
        assert_eq!(binary_json, expected["binary_path"]);
        assert_eq!(infinite_json, expected["infinite_path"]);
    }

    #[test]
    fn pinned_qpa_dual_enveloping_module_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalgtensor_dual_module.schema-v1.json"
        ))
        .unwrap();

        let loop_quiver =
            Quiver::from_named(vec!["v".into()], vec![NamedArrowSpec::named("v", "v", "x")])
                .unwrap();
        let parent = PathAlgebra::new(Q, loop_quiver.clone()).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let rational_algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let cached_before = rational_algebra.has_dual_enveloping_module_attribute();
        let bimodule_dual = dual_of_algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let bimodule_again = dual_of_algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let regular = algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let generic_dual = dual_of_module(&regular).unwrap();
        let generic_again = dual_of_module(&regular).unwrap();
        let generic_back = dual_of_module(&generic_dual).unwrap();
        let rational_enveloping = enveloping_algebra(&rational_algebra).unwrap();
        let rational_json = json!({
            "attribute_identity": [
                cached_before,
                rational_algebra.has_dual_enveloping_module_attribute(),
                bimodule_dual.is_identical(&bimodule_again),
                bimodule_dual.right_acting_algebra() == rational_enveloping.algebra(),
            ],
            "generic_identity": [
                generic_dual.is_identical(&generic_again),
                generic_back.is_identical(&regular),
            ],
            "generic_dual": module_json(&generic_dual, ToString::to_string),
            "bimodule_dual": module_json(&bimodule_dual, ToString::to_string),
        });

        let path_quiver = Quiver::from_named(
            vec!["u1".into(), "u2".into()],
            vec![NamedArrowSpec::named("u1", "u2", "a")],
        )
        .unwrap();
        let binary_algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, path_quiver).unwrap());
        let binary_dual = dual_of_algebra_as_module_over_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_again = dual_of_algebra_as_module_over_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_enveloping = enveloping_algebra(&binary_algebra).unwrap();
        let binary_json = json!({
            "attribute_identity": [
                binary_dual.is_identical(&binary_again),
                binary_dual.right_acting_algebra() == binary_enveloping.algebra(),
            ],
            "bimodule_dual": module_json(&binary_dual, |entry| {
                if Z2.is_zero(entry) {
                    "0*Z(2)".to_owned()
                } else {
                    "Z(2)^0".to_owned()
                }
            }),
        });

        let infinite_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, loop_quiver).unwrap());
        let infinite_result =
            dual_of_algebra_as_module_over_enveloping_algebra(&infinite_algebra).unwrap();
        let infinite_json = json!({
            "is_fail": infinite_result.is_none(),
            "attribute_bound": infinite_algebra.has_dual_enveloping_module_attribute(),
            "regular_bound": infinite_algebra.has_regular_enveloping_module_attribute(),
            "enveloping_bound": infinite_algebra.cached_enveloping_algebra().is_some(),
        });

        assert_eq!(rational_json, expected["rational_quotient"]);
        assert_eq!(binary_json, expected["binary_path"]);
        assert_eq!(infinite_json, expected["infinite_path"]);
    }

    #[test]
    fn pinned_qpa_random_module_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/module_random.schema-v1.json"
        ))
        .unwrap();
        let one = std::num::NonZeroUsize::new(1).unwrap();

        let mut binary_rng = ZeroRng;
        let binary =
            random_module(&isolated_path_algebra(Z2, 1), one, one, &mut binary_rng).unwrap();
        let binary_json = json!({
            "success": true,
            "is_fail": binary.is_none(),
            "dimension_vector": binary.as_ref().map(PathAlgebraModule::dimension_vector),
        });

        let mut ternary_rng = StdRng::seed_from_u64(1);
        let ternary = random_module(
            &isolated_path_algebra(Zp::new(3), 1),
            one,
            one,
            &mut ternary_rng,
        )
        .unwrap();
        let ternary_json = json!({
            "success": true,
            "is_fail": ternary.is_none(),
            "dimension_vector": ternary.as_ref().map(PathAlgebraModule::dimension_vector),
        });

        let mut rational_rng = StdRng::seed_from_u64(1);
        let rational =
            random_module(&isolated_path_algebra(Q, 1), one, one, &mut rational_rng).unwrap();
        let rational_json = json!({
            "success": true,
            "is_fail": rational.is_none(),
            "dimension_vector": rational.as_ref().map(PathAlgebraModule::dimension_vector),
        });

        let mut isolated_rng = StdRng::seed_from_u64(1);
        let isolated =
            random_module(&isolated_path_algebra(Z2, 100), one, one, &mut isolated_rng).unwrap();
        let isolated_json = json!({
            "success": true,
            "is_fail": isolated.is_none(),
            "dimension_vector": isolated.as_ref().map(PathAlgebraModule::dimension_vector),
        });

        assert_eq!(binary_json, expected["binary"]);
        assert_eq!(ternary_json, expected["ternary"]);
        assert_eq!(rational_json, expected["rational"]);
        assert_eq!(isolated_json, expected["isolated"]);
    }

    fn modulehomalg_gorenstein_dimension_case<F: QpaCoefficientField>(field: &F) -> Value {
        let discrete = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
            )
            .unwrap(),
        );
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let radical_square = fresh_radical_square_zero_a3(field);
        let known_finite = fresh_radical_square_zero_a3(field);
        let dual = fresh_dual_numbers(field);
        let known_infinite = fresh_dual_numbers(field);
        let negative_one = Integer::from(-1);
        let zero_bound = Integer::from(0);
        let two_bound = Integer::from(2);

        let discrete_bound_0 = gorenstein_dimension_of_algebra(&discrete, &zero_bound).unwrap();
        let a2_bound_0 = gorenstein_dimension_of_algebra(&a2, &zero_bound).unwrap();
        let a2_attribute = gorenstein_dimension(&a2).unwrap();
        let a2_cached_negative = gorenstein_dimension_of_algebra(&a2, &negative_one).unwrap();

        let radical_initially_bound = gorenstein_dimension(&radical_square).is_some();
        let radical_bound_0 =
            gorenstein_dimension_of_algebra(&radical_square, &zero_bound).unwrap();
        let radical_bound_after_0 = gorenstein_dimension(&radical_square).is_some();
        let radical_bound_2 = gorenstein_dimension_of_algebra(&radical_square, &two_bound).unwrap();
        let radical_attribute = gorenstein_dimension(&radical_square).unwrap();
        let radical_cached_bound_0 =
            gorenstein_dimension_of_algebra(&radical_square, &zero_bound).unwrap();

        let known_finite_setup = global_dimension_of_algebra(&known_finite, 2).unwrap();
        let known_finite_negative =
            gorenstein_dimension_of_algebra(&known_finite, &negative_one).unwrap();

        let dual_bound_0 = gorenstein_dimension_of_algebra(&dual, &zero_bound).unwrap();
        let dual_attribute = gorenstein_dimension(&dual).unwrap();
        let known_infinite_setup = global_dimension_of_algebra(&known_infinite, 0).unwrap();
        let known_infinite_bound_0 =
            gorenstein_dimension_of_algebra(&known_infinite, &zero_bound).unwrap();

        json!({
            "discrete_bound_0": bounded_dimension_json(discrete_bound_0),
            "a2_bound_0": bounded_dimension_json(a2_bound_0),
            "a2_attribute": a2_attribute,
            "a2_cached_negative": bounded_dimension_json(a2_cached_negative),
            "radical_initially_bound": radical_initially_bound,
            "radical_bound_0": bounded_dimension_json(radical_bound_0),
            "radical_bound_after_0": radical_bound_after_0,
            "radical_bound_2": bounded_dimension_json(radical_bound_2),
            "radical_attribute": radical_attribute,
            "radical_cached_bound_0": bounded_dimension_json(radical_cached_bound_0),
            "known_finite_setup": global_dimension_json(known_finite_setup),
            "known_finite_negative": bounded_dimension_json(known_finite_negative),
            "dual_bound_0": bounded_dimension_json(dual_bound_0),
            "dual_attribute": dual_attribute,
            "known_infinite_setup": global_dimension_json(known_infinite_setup),
            "known_infinite_bound_0": bounded_dimension_json(known_infinite_bound_0),
        })
    }

    #[test]
    fn pinned_qpa_modulehomalg_gorenstein_dimension_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehomalg_gorenstein_dimension.schema-v1.json"
        ))
        .unwrap();
        let actual = json!({
            "schema_version": 1,
            "kind": "modulehomalg_gorenstein_dimension",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": modulehomalg_gorenstein_dimension_case(&Z2),
            "quinary": modulehomalg_gorenstein_dimension_case(&Zp::new(5)),
            "rational": modulehomalg_gorenstein_dimension_case(&Q),
        });
        assert_eq!(actual, expected);
    }

    fn modulehomalg_dominant_dimension_case<F: QpaCoefficientField>(field: &F) -> Value {
        let zero_scalar = field.zero();
        let one_scalar = field.one();
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let simple_two = right_module_over_path_algebra(a2.clone(), vec![0, 1], &[])
            .unwrap()
            .unwrap();
        let projective = right_module_over_path_algebra(
            a2.clone(),
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one_scalar]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let zero = zero_module(&a2).unwrap();

        let loop_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_simple = right_module_over_path_algebra(
            dual.clone(),
            vec![1],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(vec![vec![zero_scalar]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let discrete = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
            )
            .unwrap(),
        );
        let negative_one = Integer::from(-1);
        let zero_bound = Integer::from(0);
        let one_bound = Integer::from(1);
        let two_bound = Integer::from(2);

        json!({
            "simple_one_bound_0": dominant_dimension_json(
                dominant_dimension_of_module(&simple_one, &zero_bound).unwrap()),
            "simple_one_bound_2": dominant_dimension_json(
                dominant_dimension_of_module(&simple_one, &two_bound).unwrap()),
            "simple_one_bound_negative": dominant_dimension_json(
                dominant_dimension_of_module(&simple_one, &negative_one).unwrap()),
            "simple_two_bound_0": dominant_dimension_json(
                dominant_dimension_of_module(&simple_two, &zero_bound).unwrap()),
            "simple_two_bound_1": dominant_dimension_json(
                dominant_dimension_of_module(&simple_two, &one_bound).unwrap()),
            "simple_two_bound_2": dominant_dimension_json(
                dominant_dimension_of_module(&simple_two, &two_bound).unwrap()),
            "simple_two_bound_negative": dominant_dimension_json(
                dominant_dimension_of_module(&simple_two, &negative_one).unwrap()),
            "projective_injective": dominant_dimension_json(
                dominant_dimension_of_module(&projective, &zero_bound).unwrap()),
            "zero": dominant_dimension_json(
                dominant_dimension_of_module(&zero, &zero_bound).unwrap()),
            "dual_simple_bound_2": dominant_dimension_json(
                dominant_dimension_of_module(&dual_simple, &two_bound).unwrap()),
            "a2_bound_0": dominant_dimension_json(
                dominant_dimension_of_algebra(&a2, &zero_bound).unwrap()),
            "a2_bound_1": dominant_dimension_json(
                dominant_dimension_of_algebra(&a2, &one_bound).unwrap()),
            "a2_bound_2": dominant_dimension_json(
                dominant_dimension_of_algebra(&a2, &two_bound).unwrap()),
            "dual_algebra": dominant_dimension_json(
                dominant_dimension_of_algebra(&dual, &two_bound).unwrap()),
            "discrete_algebra": dominant_dimension_json(
                dominant_dimension_of_algebra(&discrete, &zero_bound).unwrap()),
        })
    }

    #[test]
    fn pinned_qpa_modulehomalg_dominant_dimension_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehomalg_dominant_dimension.schema-v1.json"
        ))
        .unwrap();
        let actual = json!({
            "schema_version": 1,
            "kind": "modulehomalg_dominant_dimension",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": modulehomalg_dominant_dimension_case(&Z2),
            "quinary": modulehomalg_dominant_dimension_case(&Zp::new(5)),
            "rational": modulehomalg_dominant_dimension_case(&Q),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn dominant_dimension_preserves_infinite_dimensional_try_next_method() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
            )
            .unwrap(),
        );
        let module = right_module_over_path_algebra(
            algebra.clone(),
            vec![1],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(vec![vec![Z2.zero()]], Z2).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        assert_eq!(
            dominant_dimension_of_module(&module, &Integer::from(0)).unwrap(),
            None
        );
        assert_eq!(
            dominant_dimension_of_algebra(&algebra, &Integer::from(0)).unwrap(),
            None
        );
    }

    fn modulehomalg_dimensions_case<F: QpaCoefficientField>(field: &F) -> Value {
        let negative_one = Integer::from(-1);
        let zero_bound = Integer::from(0);
        let one_bound = Integer::from(1);
        let two_bound = Integer::from(2);

        let pd_one = fresh_a2_module(field, vec![1, 0]);
        let pd_one_initially_bound = proj_dimension(&pd_one).is_some();
        let pd_one_bound_0 = proj_dimension_of_module(&pd_one, &zero_bound).unwrap();
        let pd_one_bound_after_0 = proj_dimension(&pd_one).is_some();
        let pd_one_bound_1 = proj_dimension_of_module(&pd_one, &one_bound).unwrap();
        let pd_one_attribute = proj_dimension(&pd_one).unwrap();
        let pd_one_cached_bound_0 = proj_dimension_of_module(&pd_one, &zero_bound).unwrap();
        let pd_one_cached_bound_2 = proj_dimension_of_module(&pd_one, &two_bound).unwrap();

        let pd_two = fresh_a2_module(field, vec![0, 1]);
        let pd_two_bound_0 = proj_dimension_of_module(&pd_two, &zero_bound).unwrap();
        let pd_dual = fresh_dual_simple(field);
        let pd_dual_bound_2 = proj_dimension_of_module(&pd_dual, &two_bound).unwrap();
        let pd_negative = fresh_a2_module(field, vec![0, 1]);
        let pd_negative_bound = proj_dimension_of_module(&pd_negative, &negative_one).unwrap();
        let pd_negative_bound_after = proj_dimension(&pd_negative).is_some();
        let pd_zero = fresh_a2_zero(field);
        let pd_zero_bound_0 = proj_dimension_of_module(&pd_zero, &zero_bound).unwrap();

        let id_two = fresh_a2_module(field, vec![0, 1]);
        let id_two_initially_bound = inj_dimension(&id_two).is_some();
        let id_two_bound_0 = inj_dimension_of_module(&id_two, &zero_bound).unwrap();
        let id_two_bound_after_0 = inj_dimension(&id_two).is_some();
        let id_two_bound_1 = inj_dimension_of_module(&id_two, &one_bound).unwrap();
        let id_two_attribute = inj_dimension(&id_two).unwrap();
        let id_two_cached_bound_0 = inj_dimension_of_module(&id_two, &zero_bound).unwrap();
        let id_two_cached_bound_2 = inj_dimension_of_module(&id_two, &two_bound).unwrap();

        let id_one = fresh_a2_module(field, vec![1, 0]);
        let id_one_bound_0 = inj_dimension_of_module(&id_one, &zero_bound).unwrap();
        let id_dual = fresh_dual_simple(field);
        let id_dual_bound_2 = inj_dimension_of_module(&id_dual, &two_bound).unwrap();
        let id_negative = fresh_a2_module(field, vec![1, 0]);
        let id_negative_bound = inj_dimension_of_module(&id_negative, &negative_one).unwrap();
        let id_negative_bound_after = inj_dimension(&id_negative).is_some();
        let id_zero = fresh_a2_zero(field);
        let id_zero_bound_0 = inj_dimension_of_module(&id_zero, &zero_bound).unwrap();

        let shortcut = fresh_a2_module(field, vec![0, 1]);
        let shortcut_dual = dual_of_module(&shortcut).unwrap();
        let id_dual_cache_setup = proj_dimension_of_module(&shortcut_dual, &one_bound).unwrap();
        let id_dual_cache_shortcut_bound_0 =
            inj_dimension_of_module(&shortcut, &zero_bound).unwrap();
        let id_dual_cache_shortcut_bound_after = inj_dimension(&shortcut).is_some();

        json!({
            "pd_one_initially_bound": pd_one_initially_bound,
            "pd_one_bound_0": bounded_dimension_json(pd_one_bound_0),
            "pd_one_bound_after_0": pd_one_bound_after_0,
            "pd_one_bound_1": bounded_dimension_json(pd_one_bound_1),
            "pd_one_attribute": pd_one_attribute,
            "pd_one_cached_bound_0": bounded_dimension_json(pd_one_cached_bound_0),
            "pd_one_cached_bound_2": bounded_dimension_json(pd_one_cached_bound_2),
            "pd_two_bound_0": bounded_dimension_json(pd_two_bound_0),
            "pd_dual_bound_2": bounded_dimension_json(pd_dual_bound_2),
            "pd_negative_bound": bounded_dimension_json(pd_negative_bound),
            "pd_negative_bound_after": pd_negative_bound_after,
            "pd_zero_bound_0": bounded_dimension_json(pd_zero_bound_0),
            "id_two_initially_bound": id_two_initially_bound,
            "id_two_bound_0": bounded_dimension_json(id_two_bound_0),
            "id_two_bound_after_0": id_two_bound_after_0,
            "id_two_bound_1": bounded_dimension_json(id_two_bound_1),
            "id_two_attribute": id_two_attribute,
            "id_two_cached_bound_0": bounded_dimension_json(id_two_cached_bound_0),
            "id_two_cached_bound_2": bounded_dimension_json(id_two_cached_bound_2),
            "id_one_bound_0": bounded_dimension_json(id_one_bound_0),
            "id_dual_bound_2": bounded_dimension_json(id_dual_bound_2),
            "id_negative_bound": bounded_dimension_json(id_negative_bound),
            "id_negative_bound_after": id_negative_bound_after,
            "id_zero_bound_0": bounded_dimension_json(id_zero_bound_0),
            "id_dual_cache_setup": bounded_dimension_json(id_dual_cache_setup),
            "id_dual_cache_shortcut_bound_0": bounded_dimension_json(
                id_dual_cache_shortcut_bound_0),
            "id_dual_cache_shortcut_bound_after": id_dual_cache_shortcut_bound_after,
        })
    }

    #[test]
    fn pinned_qpa_modulehomalg_dimensions_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehomalg_dimensions.schema-v1.json"
        ))
        .unwrap();
        let actual = json!({
            "schema_version": 1,
            "kind": "modulehomalg_dimensions",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": modulehomalg_dimensions_case(&Z2),
            "quinary": modulehomalg_dimensions_case(&Zp::new(5)),
            "rational": modulehomalg_dimensions_case(&Q),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn pinned_qpa_global_dimension_of_algebra_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/modulehom_global_dimension.schema-v1.json"
        ))
        .unwrap();

        let discrete = isolated_path_algebra(Z2, 1);
        let a2_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let a2 = QuiverAlgebra::Path(PathAlgebra::new(Z2, a2_quiver).unwrap());

        let a3_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(3), "y"),
            ],
        )
        .unwrap();
        let a3_parent = PathAlgebra::new(Z2, a3_quiver).unwrap();
        let arrows = a3_parent.quiver().arrows().collect::<Vec<_>>();
        let x = a3_parent.element(arrows[0].path()).unwrap();
        let y = a3_parent.element(arrows[1].path()).unwrap();
        let radical_square_zero = QuiverAlgebra::Quotient(
            quotient_by_relations(a3_parent, &[x.multiply(&y).unwrap()]).unwrap(),
        );

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
        let loop_parent = PathAlgebra::new(Z2, loop_quiver).unwrap();
        let loop_arrow = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[loop_arrow.multiply(&loop_arrow).unwrap()])
                .unwrap(),
        );

        let result = json!({
            "discrete": [global_dimension_json(global_dimension_of_algebra(&discrete, 0).unwrap())],
            "path_a2": [
                global_dimension_json(global_dimension_of_algebra(&a2, 0).unwrap()),
                global_dimension_json(global_dimension_of_algebra(&a2, 3).unwrap())
            ],
            "radical_square_zero_a3": [
                global_dimension_json(global_dimension_of_algebra(&radical_square_zero, 0).unwrap()),
                global_dimension_json(global_dimension_of_algebra(&radical_square_zero, 1).unwrap()),
                global_dimension_json(global_dimension_of_algebra(&radical_square_zero, 2).unwrap()),
                global_dimension_json(global_dimension_of_algebra(&radical_square_zero, 0).unwrap())
            ],
            "dual_numbers": [global_dimension_json(global_dimension_of_algebra(&dual_numbers, 0).unwrap())]
        });
        assert_eq!(result, expected["f2"]);
    }

    #[allow(clippy::too_many_lines)]
    fn algebra_homological_property_cases<F: QpaCoefficientField>(field: F) -> Value {
        let zero = field.zero();
        let one = field.one();
        let semisimple_structure = StructureConstantsAlgebra::new(
            field.clone(),
            2,
            vec![
                one.clone(),
                zero.clone(),
                zero.clone(),
                zero.clone(),
                zero.clone(),
                zero.clone(),
                zero.clone(),
                one.clone(),
            ],
        )
        .unwrap();
        let dual_structure = StructureConstantsAlgebra::new(
            field.clone(),
            2,
            vec![
                one.clone(),
                zero.clone(),
                zero.clone(),
                one.clone(),
                zero.clone(),
                one.clone(),
                zero.clone(),
                zero.clone(),
            ],
        )
        .unwrap();

        let structure_semisimple_known_before = semisimple_structure
            .cached_is_semisimple_algebra()
            .is_some();
        let structure_semisimple_value = semisimple_structure.is_semisimple_algebra().unwrap();
        let structure_semisimple_known_after = semisimple_structure
            .cached_is_semisimple_algebra()
            .is_some();
        let structure_global_known_before =
            semisimple_structure.cached_global_dimension().is_some();
        let structure_global_dimension = cached_global_dimension_json(
            semisimple_structure.global_dimension_from_known_properties(),
        );
        let structure_global_known_after = semisimple_structure.cached_global_dimension().is_some();

        let structure_dual_known_before = dual_structure.cached_is_semisimple_algebra().is_some();
        let structure_dual_value = dual_structure.is_semisimple_algebra().unwrap();
        let structure_dual_known_after = dual_structure.cached_is_semisimple_algebra().is_some();

        let discrete = isolated_path_algebra(field.clone(), 1);
        let a2_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let a2_parent = PathAlgebra::new(field.clone(), a2_quiver).unwrap();
        let a2 = QuiverAlgebra::Path(a2_parent.clone());
        let a2_quotient =
            QuiverAlgebra::Quotient(quotient_by_relations(a2_parent.clone(), &[]).unwrap());

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
        let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let loop_path = QuiverAlgebra::Path(loop_parent.clone());
        let t = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let t_squared = t.multiply(&t).unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent.clone(), std::slice::from_ref(&t_squared)).unwrap(),
        );
        let dual_numbers_reverse = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, std::slice::from_ref(&t_squared)).unwrap(),
        );

        let a3_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(3), "y"),
            ],
        )
        .unwrap();
        let a3_parent = PathAlgebra::new(field, a3_quiver).unwrap();
        let arrows = a3_parent.quiver().arrows().collect::<Vec<_>>();
        let x = a3_parent.element(arrows[0].path()).unwrap();
        let y = a3_parent.element(arrows[1].path()).unwrap();
        let radical_square_zero_a3 = QuiverAlgebra::Quotient(
            quotient_by_relations(a3_parent, &[x.multiply(&y).unwrap()]).unwrap(),
        );

        let QuiverAlgebra::Quotient(a2_quotient_object) = &a2_quotient else {
            unreachable!()
        };
        let a2_quotient_semisimple = admissible_quotient_is_semisimple(a2_quotient_object).unwrap();
        let a2_quotient_hereditary = is_hereditary_algebra(&a2_quotient).unwrap().unwrap();

        let QuiverAlgebra::Quotient(dual_numbers_object) = &dual_numbers else {
            unreachable!()
        };
        let dual_global_initially = dual_numbers_object.cached_global_dimension().is_some();
        let dual_selfinjective = is_selfinjective_algebra(&dual_numbers).unwrap().unwrap();
        let dual_global_after_selfinjective =
            dual_numbers_object.cached_global_dimension().is_some();
        let dual_semisimple = admissible_quotient_is_semisimple(dual_numbers_object).unwrap();
        let dual_global_after_semisimple = dual_numbers_object.cached_global_dimension().is_some();
        let dual_global_dimension =
            cached_global_dimension_json(global_dimension_from_known_properties(&dual_numbers));
        let dual_global_after_global = dual_numbers_object.cached_global_dimension().is_some();
        let dual_hereditary = is_hereditary_algebra(&dual_numbers).unwrap().unwrap();

        let QuiverAlgebra::Quotient(dual_reverse_object) = &dual_numbers_reverse else {
            unreachable!()
        };
        let dual_reverse_semisimple =
            admissible_quotient_is_semisimple(dual_reverse_object).unwrap();
        let dual_reverse_global_after_semisimple =
            dual_reverse_object.cached_global_dimension().is_some();
        let dual_reverse_selfinjective = is_selfinjective_algebra(&dual_numbers_reverse)
            .unwrap()
            .unwrap();
        let dual_reverse_global_after_selfinjective =
            dual_reverse_object.cached_global_dimension().is_some();

        let QuiverAlgebra::Quotient(a3_object) = &radical_square_zero_a3 else {
            unreachable!()
        };
        let a3_selfinjective = is_selfinjective_algebra(&radical_square_zero_a3)
            .unwrap()
            .unwrap();
        let a3_semisimple = admissible_quotient_is_semisimple(a3_object).unwrap();
        let a3_hereditary = is_hereditary_algebra(&radical_square_zero_a3)
            .unwrap()
            .unwrap();
        let a3_global_after_hereditary = a3_object.cached_global_dimension().is_some();
        let a3_bounded =
            global_dimension_json(global_dimension_of_algebra(&radical_square_zero_a3, 2).unwrap());
        let a3_global_after_bounded = a3_object.cached_global_dimension().is_some();

        json!({
            "structure_semisimple": {
                "known_before": structure_semisimple_known_before,
                "value": structure_semisimple_value,
                "known_after": structure_semisimple_known_after,
                "global_known_before": structure_global_known_before,
                "global_dimension": structure_global_dimension,
                "global_known_after": structure_global_known_after,
                "finite_global_dimension": semisimple_structure
                    .is_finite_global_dimension_from_known_properties().unwrap()
            },
            "structure_dual_numbers": {
                "known_before": structure_dual_known_before,
                "value": structure_dual_value,
                "known_after": structure_dual_known_after,
                "global_known": dual_structure.cached_global_dimension().is_some()
            },
            "path_discrete": {
                "selfinjective": is_selfinjective_algebra(&discrete).unwrap().unwrap(),
                "semisimple": discrete.quiver().number_of_arrows() == 0,
                "hereditary": is_hereditary_algebra(&discrete).unwrap().unwrap(),
                "finite_global_dimension":
                    is_finite_global_dimension_from_known_properties(&discrete).unwrap(),
                "global_dimension": cached_global_dimension_json(
                    global_dimension_from_known_properties(&discrete))
            },
            "path_a2": {
                "selfinjective": is_selfinjective_algebra(&a2).unwrap().unwrap(),
                "semisimple": a2.quiver().number_of_arrows() == 0,
                "hereditary": is_hereditary_algebra(&a2).unwrap().unwrap(),
                "finite_global_dimension":
                    is_finite_global_dimension_from_known_properties(&a2).unwrap(),
                "global_dimension": cached_global_dimension_json(
                    global_dimension_from_known_properties(&a2))
            },
            "path_loop": {
                "selfinjective": is_selfinjective_algebra(&loop_path).unwrap().unwrap(),
                "hereditary": is_hereditary_algebra(&loop_path).unwrap().unwrap(),
                "finite_global_dimension":
                    is_finite_global_dimension_from_known_properties(&loop_path).unwrap(),
                "global_dimension": cached_global_dimension_json(
                    global_dimension_from_known_properties(&loop_path))
            },
            "quotient_a2_zero_ideal": {
                "semisimple": a2_quotient_semisimple,
                "hereditary": a2_quotient_hereditary,
                "finite_global_dimension":
                    is_finite_global_dimension_from_known_properties(&a2_quotient).unwrap(),
                "global_dimension": cached_global_dimension_json(
                    global_dimension_from_known_properties(&a2_quotient))
            },
            "quotient_dual_numbers": {
                "global_known_initially": dual_global_initially,
                "selfinjective": dual_selfinjective,
                "global_known_after_selfinjective": dual_global_after_selfinjective,
                "semisimple": dual_semisimple,
                "global_known_after_semisimple": dual_global_after_semisimple,
                "global_dimension": dual_global_dimension,
                "global_known_after_global": dual_global_after_global,
                "hereditary": dual_hereditary
            },
            "quotient_dual_numbers_semisimple_first": {
                "semisimple": dual_reverse_semisimple,
                "global_known_after_semisimple": dual_reverse_global_after_semisimple,
                "selfinjective": dual_reverse_selfinjective,
                "global_known_after_selfinjective": dual_reverse_global_after_selfinjective,
                "global_dimension": cached_global_dimension_json(
                    global_dimension_from_known_properties(&dual_numbers_reverse))
            },
            "quotient_radical_square_zero_a3": {
                "selfinjective": a3_selfinjective,
                "semisimple": a3_semisimple,
                "hereditary": a3_hereditary,
                "global_known_after_hereditary": a3_global_after_hereditary,
                "bounded_global_dimension": a3_bounded,
                "global_known_after_bounded": a3_global_after_bounded
            }
        })
    }

    #[test]
    fn pinned_qpa_algebra_homological_properties_match() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/algebra_homological_properties.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(algebra_homological_property_cases(Z2), expected["f2"]);
        assert_eq!(algebra_homological_property_cases(Q), expected["rational"]);

        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let parent = PathAlgebra::new(Z2, quiver).unwrap();
        let arrow = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let nonadmissible = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, std::slice::from_ref(&arrow)).unwrap(),
        );
        assert_eq!(is_selfinjective_algebra(&nonadmissible).unwrap(), None);
        assert_eq!(is_hereditary_algebra(&nonadmissible).unwrap(), None);
    }

    fn symmetry_property_json<F: Field>(algebra: &QuiverAlgebra<F>) -> Value {
        json!({
            "symmetric": is_symmetric_algebra(algebra).unwrap().unwrap(),
            "weakly_symmetric": is_weakly_symmetric_algebra(algebra).unwrap().unwrap(),
            "schurian": is_schurian_algebra(algebra).unwrap().unwrap(),
            "semicommutative": is_semicommutative_algebra(algebra).unwrap().unwrap()
        })
    }

    #[allow(clippy::too_many_lines)]
    fn pathalg_symmetry_property_cases<F: Field>(field: F) -> Value {
        let discrete = isolated_path_algebra(field.clone(), 1);
        let a2_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
        let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let loop_path = QuiverAlgebra::Path(loop_parent.clone());
        let loop_arrow = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[loop_arrow.multiply(&loop_arrow).unwrap()])
                .unwrap(),
        );

        let diamond_quiver = Quiver::from_indexed(
            4,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(1), vertex(3), "b"),
                ArrowSpec::named(vertex(2), vertex(4), "c"),
                ArrowSpec::named(vertex(3), vertex(4), "d"),
            ],
        )
        .unwrap();
        let diamond_parent = PathAlgebra::new(field.clone(), diamond_quiver).unwrap();
        let diamond_arrows = diamond_parent.quiver().arrows().collect::<Vec<_>>();
        let diamond_a = diamond_parent.element(diamond_arrows[0].path()).unwrap();
        let diamond_b = diamond_parent.element(diamond_arrows[1].path()).unwrap();
        let diamond_c = diamond_parent.element(diamond_arrows[2].path()).unwrap();
        let diamond_d = diamond_parent.element(diamond_arrows[3].path()).unwrap();
        let ac = diamond_a.multiply(&diamond_c).unwrap();
        let bd = diamond_b.multiply(&diamond_d).unwrap();
        let diamond_commutative = QuiverAlgebra::Quotient(
            quotient_by_relations(diamond_parent.clone(), &[ac.sub(&bd).unwrap()]).unwrap(),
        );
        let diamond_one_path_zero = QuiverAlgebra::Quotient(
            quotient_by_relations(diamond_parent.clone(), std::slice::from_ref(&ac)).unwrap(),
        );
        let diamond_path = QuiverAlgebra::Path(diamond_parent);

        let cycle_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(3), "y"),
                ArrowSpec::named(vertex(3), vertex(1), "z"),
            ],
        )
        .unwrap();
        let cycle_parent = PathAlgebra::new(field, cycle_quiver).unwrap();
        let cycle_arrows = cycle_parent.quiver().arrows().collect::<Vec<_>>();
        let cycle_x = cycle_parent.element(cycle_arrows[0].path()).unwrap();
        let cycle_y = cycle_parent.element(cycle_arrows[1].path()).unwrap();
        let cycle_z = cycle_parent.element(cycle_arrows[2].path()).unwrap();
        let radical_square_zero_cycle = QuiverAlgebra::Quotient(
            quotient_by_relations(
                cycle_parent,
                &[
                    cycle_x.multiply(&cycle_y).unwrap(),
                    cycle_y.multiply(&cycle_z).unwrap(),
                    cycle_z.multiply(&cycle_x).unwrap(),
                ],
            )
            .unwrap(),
        );

        json!({
            "path_discrete": symmetry_property_json(&discrete),
            "path_a2": symmetry_property_json(&a2),
            "path_loop": {
                "symmetric": is_symmetric_algebra(&loop_path).unwrap().unwrap(),
                "weakly_symmetric": is_weakly_symmetric_algebra(&loop_path).unwrap().unwrap()
            },
            "path_diamond": symmetry_property_json(&diamond_path),
            "quotient_dual_numbers": symmetry_property_json(&dual_numbers),
            "quotient_radical_square_zero_cycle":
                symmetry_property_json(&radical_square_zero_cycle),
            "quotient_diamond_commutative": symmetry_property_json(&diamond_commutative),
            "quotient_diamond_one_path_zero": symmetry_property_json(&diamond_one_path_zero)
        })
    }

    #[test]
    fn pinned_qpa_pathalg_symmetry_properties_match() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_symmetry_properties.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(pathalg_symmetry_property_cases(Z2), expected["f2"]);
        assert_eq!(pathalg_symmetry_property_cases(Q), expected["rational"]);

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
        let loop_path = QuiverAlgebra::Path(PathAlgebra::new(Z2, loop_quiver).unwrap());
        assert_eq!(
            is_schurian_algebra(&loop_path),
            Err(EnvelopingModuleError::InfiniteSchurianAlgebra)
        );
        assert_eq!(
            is_semicommutative_algebra(&loop_path),
            Err(EnvelopingModuleError::InfiniteSchurianAlgebra)
        );

        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let parent = PathAlgebra::new(Z2, quiver).unwrap();
        let arrow = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let nonadmissible = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, std::slice::from_ref(&arrow)).unwrap(),
        );
        assert_eq!(is_symmetric_algebra(&nonadmissible).unwrap(), None);
        assert_eq!(is_weakly_symmetric_algebra(&nonadmissible).unwrap(), None);
        assert_eq!(is_schurian_algebra(&nonadmissible).unwrap(), None);
        assert_eq!(is_semicommutative_algebra(&nonadmissible).unwrap(), None);
    }

    fn rational_matrix_rows_json(matrix: &DenseMatrix<RationalField>) -> Value {
        Value::Array(
            (0..matrix.nrows())
                .map(|row| {
                    Value::Array(
                        matrix
                            .row(row)
                            .iter()
                            .map(|entry| Value::String(entry.to_string()))
                            .collect(),
                    )
                })
                .collect(),
        )
    }

    fn cartan_coxeter_case_json<F: Field>(algebra: &QuiverAlgebra<F>) -> Value {
        let cartan_initially_known = algebra.cached_cartan_matrix_attribute().is_bound();
        let matrix_initially_known = algebra.cached_coxeter_matrix_attribute().is_bound();
        let polynomial_initially_known = algebra.cached_coxeter_polynomial_attribute().is_bound();
        let cartan = cartan_matrix(algebra).unwrap().unwrap();
        let matrix = coxeter_matrix(algebra).unwrap().unwrap();
        let polynomial = coxeter_polynomial(algebra).unwrap().unwrap();
        json!({
            "cartan_initially_known": cartan_initially_known,
            "matrix_initially_known": matrix_initially_known,
            "polynomial_initially_known": polynomial_initially_known,
            "cartan": cartan,
            "coxeter": rational_matrix_rows_json(&matrix),
            "polynomial_coefficients": polynomial
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>(),
            "cartan_known_after": algebra.cached_cartan_matrix_attribute().is_bound(),
            "matrix_known_after": algebra.cached_coxeter_matrix_attribute().is_bound(),
            "polynomial_known_after": algebra.cached_coxeter_polynomial_attribute().is_bound()
        })
    }

    #[allow(clippy::too_many_lines)]
    fn pathalg_cartan_coxeter_cases<F: Field>(field: F) -> Value {
        let discrete = isolated_path_algebra(field.clone(), 2);
        let a3_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(3), "y"),
            ],
        )
        .unwrap();
        let a3 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a3_quiver).unwrap());

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
        let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let t = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[t.multiply(&t).unwrap()]).unwrap(),
        );

        let diamond_quiver = Quiver::from_indexed(
            4,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(1), vertex(3), "b"),
                ArrowSpec::named(vertex(2), vertex(4), "c"),
                ArrowSpec::named(vertex(3), vertex(4), "d"),
            ],
        )
        .unwrap();
        let diamond_parent = PathAlgebra::new(field, diamond_quiver).unwrap();
        let arrows = diamond_parent.quiver().arrows().collect::<Vec<_>>();
        let diamond_a = diamond_parent.element(arrows[0].path()).unwrap();
        let diamond_b = diamond_parent.element(arrows[1].path()).unwrap();
        let diamond_c = diamond_parent.element(arrows[2].path()).unwrap();
        let diamond_d = diamond_parent.element(arrows[3].path()).unwrap();
        let relation = diamond_a
            .multiply(&diamond_c)
            .unwrap()
            .sub(&diamond_b.multiply(&diamond_d).unwrap())
            .unwrap();
        let commutative_diamond =
            QuiverAlgebra::Quotient(quotient_by_relations(diamond_parent, &[relation]).unwrap());

        json!({
            "path_discrete": cartan_coxeter_case_json(&discrete),
            "path_a3": cartan_coxeter_case_json(&a3),
            "quotient_dual_numbers": cartan_coxeter_case_json(&dual_numbers),
            "quotient_commutative_diamond":
                cartan_coxeter_case_json(&commutative_diamond)
        })
    }

    #[test]
    fn pinned_qpa_pathalg_cartan_coxeter_match() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_cartan_coxeter.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(pathalg_cartan_coxeter_cases(Z2), expected["f2"]);
        assert_eq!(pathalg_cartan_coxeter_cases(Zp::new(5)), expected["f5"]);
        assert_eq!(pathalg_cartan_coxeter_cases(Q), expected["rational"]);

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
        let loop_path = QuiverAlgebra::Path(PathAlgebra::new(Z2, loop_quiver).unwrap());
        assert_eq!(
            loop_path.cached_cartan_matrix_attribute(),
            CachedAttribute::Unbound
        );
        assert_eq!(cartan_matrix(&loop_path).unwrap(), None);
        assert_eq!(
            loop_path.cached_cartan_matrix_attribute(),
            CachedAttribute::Fail
        );
        assert_eq!(coxeter_matrix(&loop_path).unwrap(), None);
        assert_eq!(
            loop_path.cached_coxeter_matrix_attribute(),
            CachedAttribute::Fail
        );
        assert_eq!(coxeter_polynomial(&loop_path).unwrap(), None);
        assert_eq!(
            loop_path.cached_coxeter_polynomial_attribute(),
            CachedAttribute::Fail
        );

        let arrow_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let nonadmissible_parent = PathAlgebra::new(Z2, arrow_quiver).unwrap();
        let arrow = nonadmissible_parent
            .element(
                nonadmissible_parent
                    .quiver()
                    .arrows()
                    .next()
                    .unwrap()
                    .path(),
            )
            .unwrap();
        let nonadmissible =
            QuiverAlgebra::Quotient(quotient_by_relations(nonadmissible_parent, &[arrow]).unwrap());
        assert_eq!(cartan_matrix(&nonadmissible).unwrap(), None);
        assert_eq!(coxeter_matrix(&nonadmissible).unwrap(), None);
        assert_eq!(coxeter_polynomial(&nonadmissible).unwrap(), None);

        let cycle_quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(1), "y"),
            ],
        )
        .unwrap();
        let cycle_parent = PathAlgebra::new(Z2, cycle_quiver).unwrap();
        let cycle_arrows = cycle_parent.quiver().arrows().collect::<Vec<_>>();
        let x = cycle_parent.element(cycle_arrows[0].path()).unwrap();
        let y = cycle_parent.element(cycle_arrows[1].path()).unwrap();
        let singular = QuiverAlgebra::Quotient(
            quotient_by_relations(
                cycle_parent,
                &[x.multiply(&y).unwrap(), y.multiply(&x).unwrap()],
            )
            .unwrap(),
        );
        assert_eq!(
            cartan_matrix(&singular).unwrap(),
            Some(vec![vec![1, 1], vec![1, 1]])
        );
        assert_eq!(coxeter_matrix(&singular).unwrap(), None);
        assert_eq!(
            singular.cached_coxeter_matrix_attribute(),
            CachedAttribute::Fail
        );
        assert_eq!(coxeter_polynomial(&singular).unwrap(), None);
        assert_eq!(
            singular.cached_coxeter_polynomial_attribute(),
            CachedAttribute::Fail
        );
    }

    fn nakayama_permutation_case_json<F: Field>(algebra: &QuiverAlgebra<F>) -> Value {
        let initially_known = algebra.cached_nakayama_permutation_attribute().is_bound();
        let Some(permutation) = nakayama_permutation(algebra).unwrap() else {
            return json!({
                "initially_known": initially_known,
                "value": false,
                "known_after": algebra
                    .cached_nakayama_permutation_attribute()
                    .is_bound()
            });
        };
        let simples = simple_modules(algebra).unwrap();
        let module_images = simples
            .iter()
            .map(|simple| {
                permutation
                    .image_of_module(simple)
                    .unwrap()
                    .dimension_vector()
                    .to_vec()
            })
            .collect::<Vec<_>>();
        let module_identical = simples
            .iter()
            .map(|simple| {
                permutation
                    .image_of_module(simple)
                    .unwrap()
                    .is_identical(simple)
            })
            .collect::<Vec<_>>();
        let index_images = (1..=simples.len())
            .map(|index| permutation.image_of_index(&Integer::from(index)).unwrap())
            .collect::<Vec<_>>();
        json!({
            "initially_known": initially_known,
            "index_images": index_images,
            "module_images": module_images,
            "module_identical": module_identical,
            "known_after": algebra
                .cached_nakayama_permutation_attribute()
                .is_bound()
        })
    }

    #[allow(clippy::too_many_lines)]
    fn pathalg_nakayama_permutation_cases<F: Field>(field: F) -> Value {
        let discrete = isolated_path_algebra(field.clone(), 3);
        let a2_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
        let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let loop_arrow = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[loop_arrow.multiply(&loop_arrow).unwrap()])
                .unwrap(),
        );

        let a3_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "u"),
                ArrowSpec::named(vertex(2), vertex(3), "v"),
            ],
        )
        .unwrap();
        let a3_parent = PathAlgebra::new(field.clone(), a3_quiver).unwrap();
        let a3_arrows = a3_parent.quiver().arrows().collect::<Vec<_>>();
        let first_arrow = a3_parent.element(a3_arrows[0].path()).unwrap();
        let second_arrow = a3_parent.element(a3_arrows[1].path()).unwrap();
        let nonselfinjective = QuiverAlgebra::Quotient(
            quotient_by_relations(a3_parent, &[first_arrow.multiply(&second_arrow).unwrap()])
                .unwrap(),
        );

        let cycle_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(3), "y"),
                ArrowSpec::named(vertex(3), vertex(1), "z"),
            ],
        )
        .unwrap();
        let cycle_parent = PathAlgebra::new(field, cycle_quiver).unwrap();
        let arrows = cycle_parent.quiver().arrows().collect::<Vec<_>>();
        let cycle_x = cycle_parent.element(arrows[0].path()).unwrap();
        let cycle_y = cycle_parent.element(arrows[1].path()).unwrap();
        let cycle_z = cycle_parent.element(arrows[2].path()).unwrap();
        let cycle = QuiverAlgebra::Quotient(
            quotient_by_relations(
                cycle_parent,
                &[
                    cycle_x.multiply(&cycle_y).unwrap(),
                    cycle_y.multiply(&cycle_z).unwrap(),
                    cycle_z.multiply(&cycle_x).unwrap(),
                ],
            )
            .unwrap(),
        );

        json!({
            "path_discrete": nakayama_permutation_case_json(&discrete),
            "path_a2": nakayama_permutation_case_json(&a2),
            "quotient_dual_numbers": nakayama_permutation_case_json(&dual_numbers),
            "quotient_nonselfinjective": nakayama_permutation_case_json(&nonselfinjective),
            "quotient_radical_square_zero_cycle": nakayama_permutation_case_json(&cycle)
        })
    }

    #[test]
    fn pinned_qpa_pathalg_nakayama_permutation_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_nakayama_permutation.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(pathalg_nakayama_permutation_cases(Z2), expected["f2"]);
        assert_eq!(
            pathalg_nakayama_permutation_cases(Zp::new(5)),
            expected["f5"]
        );
        assert_eq!(pathalg_nakayama_permutation_cases(Q), expected["rational"]);

        let cycle_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(3), "y"),
                ArrowSpec::named(vertex(3), vertex(1), "z"),
            ],
        )
        .unwrap();
        let cycle_parent = PathAlgebra::new(Z2, cycle_quiver).unwrap();
        let arrows = cycle_parent.quiver().arrows().collect::<Vec<_>>();
        let x = cycle_parent.element(arrows[0].path()).unwrap();
        let y = cycle_parent.element(arrows[1].path()).unwrap();
        let z = cycle_parent.element(arrows[2].path()).unwrap();
        let cycle = QuiverAlgebra::Quotient(
            quotient_by_relations(
                cycle_parent,
                &[
                    x.multiply(&y).unwrap(),
                    y.multiply(&z).unwrap(),
                    z.multiply(&x).unwrap(),
                ],
            )
            .unwrap(),
        );
        let permutation = nakayama_permutation(&cycle).unwrap().unwrap();
        let cycle_simples = simple_modules(&cycle).unwrap();
        let nonsimple = direct_sum_of_qpa_modules(&cycle_simples[..2])
            .unwrap()
            .unwrap();
        assert_eq!(
            permutation.image_of_module(&nonsimple),
            Err(EnvelopingModuleError::NakayamaPermutationModuleNotSimple)
        );
        for index in [0, -1, 4] {
            assert_eq!(
                permutation.image_of_index(&Integer::from(index)),
                Err(EnvelopingModuleError::IncorrectNakayamaPermutationIndex {
                    index: Integer::from(index),
                    vertex_count: 3
                })
            );
        }

        let foreign = isolated_path_algebra(Z2, 4);
        let foreign_simple = simple_modules(&foreign).unwrap()[3].clone();
        assert_eq!(
            permutation.image_of_module(&foreign_simple),
            Err(
                EnvelopingModuleError::NakayamaPermutationModulePositionOutOfRange {
                    position: 4,
                    vertex_count: 3
                }
            )
        );

        let arrow_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let parent = PathAlgebra::new(Z2, arrow_quiver).unwrap();
        let arrow = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let nonadmissible = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, std::slice::from_ref(&arrow)).unwrap(),
        );
        assert_eq!(
            nakayama_permutation(&nonadmissible).unwrap_err(),
            EnvelopingModuleError::NakayamaPermutationSelfinjectivityUnavailable
        );
    }

    fn radical_square_zero_cycle<F: Field>(field: F) -> QuiverAlgebra<F> {
        let quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "x"),
                ArrowSpec::named(vertex(2), vertex(3), "y"),
                ArrowSpec::named(vertex(3), vertex(1), "z"),
            ],
        )
        .unwrap();
        let parent = PathAlgebra::new(field, quiver).unwrap();
        let arrows = parent.quiver().arrows().collect::<Vec<_>>();
        let x = parent.element(arrows[0].path()).unwrap();
        let y = parent.element(arrows[1].path()).unwrap();
        let z = parent.element(arrows[2].path()).unwrap();
        QuiverAlgebra::Quotient(
            quotient_by_relations(
                parent,
                &[
                    x.multiply(&y).unwrap(),
                    y.multiply(&z).unwrap(),
                    z.multiply(&x).unwrap(),
                ],
            )
            .unwrap(),
        )
    }

    fn quantum_complete_intersection<F: Field>(field: F, scalar: &F::Element) -> QuiverAlgebra<F> {
        let quiver = Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(1), "y"),
            ],
        )
        .unwrap();
        let parent = PathAlgebra::new(field, quiver).unwrap();
        let arrows = parent.quiver().arrows().collect::<Vec<_>>();
        let x = parent.element(arrows[0].path()).unwrap();
        let y = parent.element(arrows[1].path()).unwrap();
        let xy = x.multiply(&y).unwrap();
        let yx = y.multiply(&x).unwrap();
        QuiverAlgebra::Quotient(
            quotient_by_relations(
                parent,
                &[
                    x.multiply(&x).unwrap(),
                    y.multiply(&y).unwrap(),
                    xy.sub(&yx.scale_left(scalar).unwrap()).unwrap(),
                ],
            )
            .unwrap(),
        )
    }

    fn gap_z2_scalar(value: <Z2 as numerica::domains::Set>::Element) -> String {
        if Z2.is_zero(&value) {
            "0*Z(2)".to_owned()
        } else {
            "Z(2)^0".to_owned()
        }
    }

    fn gap_zp_scalar(field: &Zp, value: <Zp as numerica::domains::Set>::Element) -> String {
        if field.is_zero(&value) {
            return format!("0*Z({})", field.get_prime());
        }
        let elements = field.qpa_elements().unwrap().unwrap();
        let position = elements
            .iter()
            .position(|candidate| candidate == &value)
            .unwrap();
        let exponent = position - 1;
        match exponent {
            1 => format!("Z({})", field.get_prime()),
            _ => format!("Z({})^{exponent}", field.get_prime()),
        }
    }

    fn scalar_matrix_json<F: Field>(
        matrix: &DenseMatrix<F>,
        scalar: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!(
            (0..matrix.nrows())
                .map(|row| matrix.row(row).iter().map(scalar).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        )
    }

    fn nakayama_order_json(value: &NakayamaAutomorphismOrder) -> Value {
        match value {
            NakayamaAutomorphismOrder::Finite(order) => {
                serde_json::from_str(&order.to_string()).unwrap()
            }
            NakayamaAutomorphismOrder::Infinite => json!("infinity"),
        }
    }

    fn nakayama_path_case_json<F: QpaCoefficientField>(
        algebra: &QuiverAlgebra<F>,
        scalar: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let auto_initially_known = algebra.cached_nakayama_automorphism_attribute().is_bound();
        let order_initially_known = algebra
            .cached_order_of_nakayama_automorphism_attribute()
            .is_bound();
        let automorphism = nakayama_automorphism(algebra).unwrap();
        let order = order_of_nakayama_automorphism(algebra).unwrap();
        json!({
            "auto_initially_known": auto_initially_known,
            "order_initially_known": order_initially_known,
            "automorphism": automorphism.as_ref().map_or(Value::Bool(false),
                |value| scalar_matrix_json(value.matrix(), scalar)),
            "order": order.as_ref().map_or(Value::Bool(false), nakayama_order_json),
            "auto_known_after": algebra.cached_nakayama_automorphism_attribute().is_bound(),
            "order_known_after": algebra
                .cached_order_of_nakayama_automorphism_attribute().is_bound()
        })
    }

    #[derive(Clone, Copy)]
    enum FrobeniusTrigger {
        Automorphism,
        Form,
        Functional,
    }

    fn nakayama_quotient_case_json<F: QpaCoefficientField>(
        algebra: &QuiverAlgebra<F>,
        trigger: FrobeniusTrigger,
        scalar: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let QuiverAlgebra::Quotient(quotient) = algebra else {
            unreachable!()
        };
        let auto_initially_known = algebra.cached_nakayama_automorphism_attribute().is_bound();
        let form_initially_known = quotient.cached_frobenius_form_attribute().is_bound();
        let functional_initially_known = quotient
            .cached_frobenius_linear_functional_attribute()
            .is_bound();
        let order_initially_known = algebra
            .cached_order_of_nakayama_automorphism_attribute()
            .is_bound();

        let (automorphism, mut form, mut functional) = match trigger {
            FrobeniusTrigger::Automorphism => (nakayama_automorphism(algebra).unwrap(), None, None),
            FrobeniusTrigger::Form => {
                let form = frobenius_form(quotient).unwrap();
                let automorphism = nakayama_automorphism(algebra).unwrap();
                (automorphism, form, None)
            }
            FrobeniusTrigger::Functional => {
                let functional = frobenius_linear_functional(quotient).unwrap();
                let automorphism = nakayama_automorphism(algebra).unwrap();
                (automorphism, None, functional)
            }
        };
        let auto_known_after_trigger = algebra.cached_nakayama_automorphism_attribute().is_bound();
        let form_known_after_trigger = quotient.cached_frobenius_form_attribute().is_bound();
        let functional_known_after_trigger = quotient
            .cached_frobenius_linear_functional_attribute()
            .is_bound();

        if form.is_none() {
            form = frobenius_form(quotient).unwrap();
        }
        if functional.is_none() {
            functional = frobenius_linear_functional(quotient).unwrap();
        }
        let mut result = json!({
            "auto_initially_known": auto_initially_known,
            "form_initially_known": form_initially_known,
            "functional_initially_known": functional_initially_known,
            "order_initially_known": order_initially_known,
            "auto_known_after_trigger": auto_known_after_trigger,
            "form_known_after_trigger": form_known_after_trigger,
            "functional_known_after_trigger": functional_known_after_trigger,
            "automorphism": automorphism.as_ref().map_or(Value::Bool(false),
                |value| scalar_matrix_json(value.matrix(), scalar)),
            "form": form.as_ref().map_or(Value::Bool(false),
                |value| scalar_matrix_json(value.matrix(), scalar)),
            "functional": functional.as_ref().map_or(Value::Bool(false), |value| {
                json!(value.canonical_basis_values().iter().map(scalar).collect::<Vec<_>>())
            }),
        });
        if automorphism.is_some() {
            result["order"] =
                nakayama_order_json(&order_of_nakayama_automorphism(algebra).unwrap().unwrap());
        }
        result["auto_known_after"] =
            json!(algebra.cached_nakayama_automorphism_attribute().is_bound());
        result["form_known_after"] = json!(quotient.cached_frobenius_form_attribute().is_bound());
        result["functional_known_after"] = json!(
            quotient
                .cached_frobenius_linear_functional_attribute()
                .is_bound()
        );
        result["order_known_after"] = json!(
            algebra
                .cached_order_of_nakayama_automorphism_attribute()
                .is_bound()
        );
        result
    }

    #[allow(clippy::too_many_lines)]
    fn pathalg_nakayama_automorphism_cases<F: QpaCoefficientField>(
        field: F,
        quantum_scalar: &F::Element,
        scalar: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        let discrete = isolated_path_algebra(field.clone(), 2);
        let a2_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let a2 = QuiverAlgebra::Path(PathAlgebra::new(field.clone(), a2_quiver).unwrap());

        let loop_quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "t")]).unwrap();
        let loop_parent = PathAlgebra::new(field.clone(), loop_quiver).unwrap();
        let loop_arrow = loop_parent
            .element(loop_parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let dual_numbers = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[loop_arrow.multiply(&loop_arrow).unwrap()])
                .unwrap(),
        );

        let a3_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "u"),
                ArrowSpec::named(vertex(2), vertex(3), "v"),
            ],
        )
        .unwrap();
        let a3_parent = PathAlgebra::new(field.clone(), a3_quiver).unwrap();
        let arrows = a3_parent.quiver().arrows().collect::<Vec<_>>();
        let u = a3_parent.element(arrows[0].path()).unwrap();
        let v = a3_parent.element(arrows[1].path()).unwrap();
        let nonselfinjective = QuiverAlgebra::Quotient(
            quotient_by_relations(a3_parent, &[u.multiply(&v).unwrap()]).unwrap(),
        );

        let cycle = radical_square_zero_cycle(field.clone());
        let quantum = quantum_complete_intersection(field, quantum_scalar);

        json!({
            "path_discrete": nakayama_path_case_json(&discrete, scalar),
            "path_a2": nakayama_path_case_json(&a2, scalar),
            "quotient_dual_numbers": nakayama_quotient_case_json(
                &dual_numbers, FrobeniusTrigger::Form, scalar),
            "quotient_nonselfinjective": nakayama_quotient_case_json(
                &nonselfinjective, FrobeniusTrigger::Automorphism, scalar),
            "quotient_radical_square_zero_cycle": nakayama_quotient_case_json(
                &cycle, FrobeniusTrigger::Automorphism, scalar),
            "quotient_quantum_complete_intersection": nakayama_quotient_case_json(
                &quantum, FrobeniusTrigger::Functional, scalar)
        })
    }

    #[test]
    fn pinned_qpa_pathalg_nakayama_automorphism_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_nakayama_automorphism.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            pathalg_nakayama_automorphism_cases(Z2, &Z2.one(), |value| gap_z2_scalar(*value)),
            expected["f2"]
        );
        let field5 = Zp::new(5);
        assert_eq!(
            pathalg_nakayama_automorphism_cases(field5.clone(), &field5.to_element(2), |value| {
                gap_zp_scalar(&field5, *value)
            },),
            expected["f5"]
        );
        assert_eq!(
            pathalg_nakayama_automorphism_cases(Q, &q(2), ToString::to_string),
            expected["rational"]
        );
    }

    #[test]
    fn pinned_qpa_pathalg_nakayama_automorphism_dispatch_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_nakayama_automorphism_dispatch.schema-v1.json"
        ))
        .unwrap();
        let cycle = radical_square_zero_cycle(Z2);
        let QuiverAlgebra::Quotient(cycle_quotient) = &cycle else {
            unreachable!()
        };
        let automorphism = nakayama_automorphism(&cycle).unwrap().unwrap();
        let form = frobenius_form(cycle_quotient).unwrap().unwrap();
        let functional = frobenius_linear_functional(cycle_quotient)
            .unwrap()
            .unwrap();
        let foreign = isolated_path_algebra(Z2, 1).one().unwrap();
        let one = cycle.one().unwrap();

        let a3_quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "u"),
                ArrowSpec::named(vertex(2), vertex(3), "v"),
            ],
        )
        .unwrap();
        let a3_parent = PathAlgebra::new(Z2, a3_quiver).unwrap();
        let arrows = a3_parent.quiver().arrows().collect::<Vec<_>>();
        let u = a3_parent.element(arrows[0].path()).unwrap();
        let v = a3_parent.element(arrows[1].path()).unwrap();
        let nonselfinjective = QuiverAlgebra::Quotient(
            quotient_by_relations(a3_parent, &[u.multiply(&v).unwrap()]).unwrap(),
        );

        let arrow_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let parent = PathAlgebra::new(Z2, arrow_quiver).unwrap();
        let arrow = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let nonadmissible = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, std::slice::from_ref(&arrow)).unwrap(),
        );
        let QuiverAlgebra::Quotient(nonself_quotient) = &nonselfinjective else {
            unreachable!()
        };
        let QuiverAlgebra::Quotient(nonadmissible_quotient) = &nonadmissible else {
            unreachable!()
        };

        let actual = json!({
            "schema_version": 1,
            "source": {
                "gap": "4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa": "1.37@8d5f847258d4853700ce3a5d947db2a1441109f1"
            },
            "foreign_auto_success": automorphism.image(&foreign).is_ok(),
            "foreign_form_left_success": form.evaluate(&foreign, &one).is_ok(),
            "foreign_form_right_success": form.evaluate(&one, &foreign).is_ok(),
            "foreign_functional_success": functional.evaluate(&foreign).is_ok(),
            "nonself_auto_is_false": nakayama_automorphism(&nonselfinjective)
                .unwrap().is_none(),
            "nonself_form_is_false": frobenius_form(nonself_quotient).unwrap().is_none(),
            "nonself_functional_is_false": frobenius_linear_functional(nonself_quotient)
                .unwrap().is_none(),
            "nonself_order_success": order_of_nakayama_automorphism(&nonselfinjective).is_ok(),
            "nonadmissible_auto_success": nakayama_automorphism(&nonadmissible).is_ok(),
            "nonadmissible_form_success": frobenius_form(nonadmissible_quotient).is_ok(),
            "nonadmissible_functional_success":
                frobenius_linear_functional(nonadmissible_quotient).is_ok(),
            "nonadmissible_order_success":
                order_of_nakayama_automorphism(&nonadmissible).is_ok()
        });
        assert_eq!(actual, expected);
        assert!(matches!(
            automorphism.image(&foreign),
            Err(EnvelopingModuleError::NakayamaAutomorphismElementOutsideAlgebra)
        ));
        assert_eq!(
            order_of_nakayama_automorphism(&nonselfinjective),
            Err(EnvelopingModuleError::CalledFalseNakayamaAutomorphism)
        );
        assert!(matches!(
            nakayama_automorphism(&nonadmissible),
            Err(EnvelopingModuleError::NakayamaAutomorphismSelfinjectivityUnavailable)
        ));
    }

    #[test]
    fn nakayama_automorphism_frobenius_data_and_order_follow_qpa() {
        let cycle = radical_square_zero_cycle(Q);
        let QuiverAlgebra::Quotient(cycle_quotient) = &cycle else {
            unreachable!()
        };
        assert!(!cycle.cached_nakayama_automorphism_attribute().is_bound());
        assert!(!cycle_quotient.cached_frobenius_form_attribute().is_bound());
        let form = frobenius_form(cycle_quotient).unwrap().unwrap();
        assert!(cycle.cached_nakayama_automorphism_attribute().is_bound());
        assert!(
            cycle_quotient
                .cached_frobenius_linear_functional_attribute()
                .is_bound()
        );
        let automorphism = nakayama_automorphism(&cycle).unwrap().unwrap();
        let expected_automorphism = DenseMatrix::from_rows(
            vec![
                vec![q(0), q(1), q(0), q(0), q(0), q(0)],
                vec![q(0), q(0), q(1), q(0), q(0), q(0)],
                vec![q(1), q(0), q(0), q(0), q(0), q(0)],
                vec![q(0), q(0), q(0), q(0), q(1), q(0)],
                vec![q(0), q(0), q(0), q(0), q(0), q(1)],
                vec![q(0), q(0), q(0), q(1), q(0), q(0)],
            ],
            Q,
        )
        .unwrap();
        assert_eq!(automorphism.matrix(), &expected_automorphism);
        assert_eq!(
            order_of_nakayama_automorphism(&cycle).unwrap(),
            Some(NakayamaAutomorphismOrder::Finite(Integer::from(3)))
        );
        let functional = frobenius_linear_functional(cycle_quotient)
            .unwrap()
            .unwrap();
        assert_eq!(
            functional.canonical_basis_values(),
            &[q(0), q(0), q(0), q(1), q(1), q(1)]
        );
        let basis = finite_canonical_basis(&cycle).unwrap().unwrap();
        for (row, left) in basis.iter().enumerate() {
            for (column, right) in basis.iter().enumerate() {
                assert_eq!(
                    form.evaluate(left, right).unwrap(),
                    form.matrix()[(row, column)]
                );
            }
            assert_eq!(
                functional.evaluate(left).unwrap(),
                functional.canonical_basis_values()[row]
            );
            assert_eq!(
                canonical_basis_coefficients(&cycle, &basis, &automorphism.image(left).unwrap())
                    .unwrap(),
                expected_automorphism.row(row)
            );
        }

        let quantum_q = quantum_complete_intersection(Q, &q(2));
        let quantum_q_auto = nakayama_automorphism(&quantum_q).unwrap().unwrap();
        assert_eq!(
            quantum_q_auto.matrix(),
            &DenseMatrix::from_rows(
                vec![
                    vec![q(1), q(0), q(0), q(0)],
                    vec![q(0), q(2), q(0), q(0)],
                    vec![q(0), q(0), Q.to_element(1.into(), 2.into(), true), q(0)],
                    vec![q(0), q(0), q(0), q(1)],
                ],
                Q,
            )
            .unwrap()
        );
        assert_eq!(
            order_of_nakayama_automorphism(&quantum_q).unwrap(),
            Some(NakayamaAutomorphismOrder::Infinite)
        );

        let field5 = Zp::new(5);
        let quantum_f5 = quantum_complete_intersection(field5.clone(), &field5.to_element(2));
        assert_eq!(
            order_of_nakayama_automorphism(&quantum_f5).unwrap(),
            Some(NakayamaAutomorphismOrder::Finite(Integer::from(4)))
        );
        let quantum_f2 = quantum_complete_intersection(Z2, &Z2.one());
        assert_eq!(
            order_of_nakayama_automorphism(&quantum_f2).unwrap(),
            Some(NakayamaAutomorphismOrder::Finite(Integer::one()))
        );

        let discrete = isolated_path_algebra(Q, 2);
        assert_eq!(
            order_of_nakayama_automorphism(&discrete).unwrap(),
            Some(NakayamaAutomorphismOrder::Finite(Integer::one()))
        );
        let a2_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let a2 = QuiverAlgebra::Path(PathAlgebra::new(Q, a2_quiver).unwrap());
        assert!(nakayama_automorphism(&a2).unwrap().is_none());
        assert!(order_of_nakayama_automorphism(&a2).unwrap().is_none());
    }
}
