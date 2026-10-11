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

type ProjectiveCoordinates = std::collections::HashMap<Vec<qpa_core::ArrowIndex>, (usize, usize)>;

struct ProjectiveBasisData<F: Field> {
    basis: Vec<QuiverAlgebraElement<F>>,
    intervals: Vec<Vec<Vec<usize>>>,
    // Lazy per-projective indexes: basis-only requests and live cache hits
    // do not build them. The empty word is the starting vertex's idempotent.
    coordinates: Vec<std::sync::OnceLock<ProjectiveCoordinates>>,
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
    let projectives = (0..algebra.quiver().number_of_vertices())
        .map(|vertex| projective_from_basis_data(algebra, &data, vertex))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(algebra.merge_indec_projective_modules(&projectives))
}

// Both full-family requests and HomFromProjective enter here. Recover live
// members individually: an expired sibling must not rebuild a live module.
fn projective_from_basis_data<F: Field>(
    algebra: &QuiverAlgebra<F>,
    data: &ProjectiveBasisData<F>,
    projective_vertex: usize,
) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_indec_projective_module(projective_vertex) {
        return Ok(cached);
    }
    let basis = &data.basis;
    let intervals = &data.intervals[projective_vertex];
    let coordinates = data.coordinates[projective_vertex].get_or_init(|| {
        intervals.iter().enumerate().flat_map(|(target, indices)| {
            indices.iter().enumerate().map(move |(coordinate, &index)| {
                let path = data.basis[index].leading_monomial();
                (path.walk().expect("a nonzero canonical basis path").to_vec(), (target, coordinate))
            })
        }).collect()
    });
    let dimensions = intervals.iter().map(Vec::len).collect::<Vec<_>>();
    let mut specification = Vec::new();
    for arrow in algebra.quiver().arrows() {
        let source = arrow.source().index().get() - 1;
        let target = arrow.target().index().get() - 1;
        if dimensions[source] == 0 || dimensions[target] == 0 {
            continue;
        }
        let arrow_element = include_in_path_algebra(&arrow.path(), algebra)?;
        let rows = intervals[source]
            .iter()
            .map(|&basis_index| {
                let image = basis[basis_index].multiply(&arrow_element)?;
                let representative = match &image {
                    QuiverAlgebraElement::Path(element) => element,
                    QuiverAlgebraElement::Quotient(element) => {
                        if !element.is_normal_form() {
                            return Err(QuotientError::ElementsNotInNormalForm.into());
                        }
                        element.representative()
                    }
                };
                let mut row = vec![algebra.field().zero(); dimensions[target]];
                for term in representative.terms() {
                    let PathEndpoint::Vertex(left) = term.path().source() else {
                        return Err(EnvelopingModuleError::ElementOutsideCanonicalBasis);
                    };
                    if left.index().get() - 1 != projective_vertex {
                        return Err(EnvelopingModuleError::ElementOutsideCanonicalBasis);
                    }
                    let &(right, coordinate) = coordinates
                        .get(term.path().walk().expect("a nonzero product term"))
                        .ok_or(EnvelopingModuleError::ElementOutsideCanonicalBasis)?;
                    if right != target {
                        return Err(EnvelopingModuleError::ElementOutsideCanonicalBasis);
                    }
                    row[coordinate] = term.coefficient().clone();
                }
                Ok(row)
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
    algebra.merge_indec_projective_module(projective_vertex, &projective)
        .ok_or(EnvelopingModuleError::ProjectiveVertexOutOfRange {
            vertex: projective_vertex + 1,
            vertex_count: algebra.quiver().number_of_vertices(),
        })
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

impl QpaRandomScalar for qpa_core::GaloisField {
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

/// Translate QPA's `N_RigidModule(M, n)` operation.
///
/// Starting with `N = M`, the source tests `ExtOverAlgebra(N, M)[2]` and then
/// replaces `N` by `NthSyzygy(N, 1)`. Thus it tests the self-Ext groups in
/// degrees one through `n`. Zero and negative bounds return `true` without
/// computing Ext. On every vanishing Ext group QPA computes the next syzygy,
/// including after the final requested degree; this translation preserves
/// that evaluation order.
///
/// # Errors
///
/// Propagates every Ext, projective-resolution, and syzygy construction
/// failure.
pub fn n_rigid_module<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<bool, EnvelopingModuleError> {
    let one = Integer::one();
    let mut degree = one.clone();
    let mut current = module.clone();
    while &degree <= bound {
        if ext_dimension_over_algebra(&current, module)? > 0 {
            return Ok(false);
        }
        degree += &one;
        current = nth_syzygy(&current, &one)?;
    }
    Ok(true)
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
        QpaFieldCharacteristic::Prime(prime) => finite_field_general_linear_exponent(
            matrix.nrows(),
            prime,
            &field.qpa_order().expect("positive-characteristic field is finite"),
        )?,
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

fn finite_field_general_linear_exponent(
    dimension: usize,
    prime: u64,
    field_order: &Integer,
) -> Result<Integer, EnvelopingModuleError> {
    let prime_usize = usize::try_from(prime)
        .map_err(|_| EnvelopingModuleError::NakayamaOrderPrimeTooLarge { prime })?;
    let prime_integer = Integer::from(prime);
    let mut exponent = Integer::one();
    let mut field_order_power = Integer::one();
    for _ in 1..=dimension {
        field_order_power *= field_order;
        exponent = exponent.lcm(&(&field_order_power - 1));
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
