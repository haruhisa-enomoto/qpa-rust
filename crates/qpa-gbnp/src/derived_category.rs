//! Projective/injective complex descriptions and helper maps translated from
//! QPA's `derivedcat.gd` and `derivedcat.gi`.

use std::{fmt::Write as _, sync::Arc};

use numerica::domains::{Field, integer::Integer};
use qpa_core::{
    DenseMatrix, ExtendedInteger, InfiniteList, PathAlgebraError, PathAlgebraModule,
    PathAlgebraModuleError, PathAlgebraModuleHomomorphism, PathAlgebraTensorError,
    QpaCoefficientField, QuiverAlgebra, QuiverAlgebraElement, QuiverAlgebraElementError,
    direct_sum_inclusions, direct_sum_of_qpa_modules, direct_sum_projections,
    include_in_path_algebra, inverse_of_isomorphism, kernel_inclusion,
    minimal_generating_set_of_module, pull_back, right_module_hom_over_algebra, socle_of_module,
    top_of_module, zero_module, zero_module_mapping,
};
use thiserror::Error;

use crate::{
    complex::{
        ComplexError, ComplexHalf, ModuleComplex, complex, finite_complex,
        projective_resolution_complex, stalk_complex, yoneda_product,
    },
    module::{
        EnvelopingModuleError, basis_of_projectives, hom_from_projective, indec_injective_modules,
        indec_projective_modules, is_injective_module, is_projective_module,
        opposite_module_algebra, projective_cover,
    },
    opposite_quotient_element,
    pathalg::PathalgGbnpError,
};

/// Failures exposed by QPA's projective/injective complex helper layer.
#[derive(Debug, Error)]
pub enum DerivedCategoryError {
    /// A complex operation failed.
    #[error(transparent)]
    Complex(#[from] ComplexError),
    /// A module or module-map operation failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// Parent path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Arithmetic in the path-algebra/quotient union failed.
    #[error(transparent)]
    AlgebraElement(#[from] QuiverAlgebraElementError),
    /// Inclusion of a path into a quotient/path-algebra union failed.
    #[error(transparent)]
    Tensor(#[from] PathAlgebraTensorError),
    /// A GBNP-backed opposite quotient construction failed.
    #[error(transparent)]
    Pathalg(#[from] PathalgGbnpError),
    /// A projective/injective module construction failed.
    #[error(transparent)]
    Enveloping(#[from] EnvelopingModuleError),
    /// The source description method requires finite bounds.
    #[error("complex entered is not finite")]
    NonFiniteComplex,
    /// QPA's projective/injective printer requires a finite lower bound.
    #[error("projective/injective complex printer requires a finite lower bound")]
    PrinterNeedsFiniteLowerBound,
    /// QPA's infinite printer requires a finite highest known degree.
    #[error("projective/injective complex printer requires a finite highest known degree")]
    PrinterNeedsFiniteHighestKnownDegree,
    /// `DirectSumMinusSummands` was called on a module without direct-sum data.
    #[error("DirectSumMinusSummands requires a direct-sum module")]
    NotDirectSumModule,
    /// A one-based removed-summand position lies outside the source list.
    #[error("summand position {position} lies outside 1..={summand_count}")]
    SummandPositionOutOfRange {
        /// Entered one-based position.
        position: usize,
        /// Number of direct summands.
        summand_count: usize,
    },
    /// QPA's sequential zero-component insertion reached an invalid position.
    #[error(
        "summand position {position} cannot be inserted into a component list of length {length}"
    )]
    ComponentInsertionOutOfRange {
        /// Entered one-based position.
        position: usize,
        /// Current component-list length.
        length: usize,
    },
    /// Component lists must match when QPA multiplies a row of maps by a
    /// column of maps.
    #[error("map component counts differ: {left} and {right}")]
    MapComponentCountMismatch {
        /// Number of row components.
        left: usize,
        /// Number of column components.
        right: usize,
    },
    /// The source helper cannot add an empty list of component composites.
    #[error("map component list is empty")]
    EmptyMapComponentList,
    /// `DescriptionOfProjectiveModule` returned GAP's `false` sentinel where
    /// the reduction method requires a projective vertex number.
    #[error("a differential summand is not described by an indecomposable projective vertex")]
    NonProjectiveSummandDescription,
    /// A one-based cancellation position is outside the component grid.
    #[error(
        "cancellation component ({source_position},{range_position}) lies outside a {source_count} by {range_count} map grid"
    )]
    CancellationPositionOutOfRange {
        /// Entered one-based source-summand position.
        source_position: usize,
        /// Entered one-based range-summand position.
        range_position: usize,
        /// Number of source summands.
        source_count: usize,
        /// Number of range summands.
        range_count: usize,
    },
    /// A component grid does not have the expected rectangular dimensions.
    #[error("map component grid is not rectangular")]
    NonRectangularMapComponentGrid,
    /// `ProjectiveResolutionOfComplex` performs integer arithmetic on the
    /// lower bound without guarding GAP's infinity objects.
    #[error("ProjectiveResolutionOfComplex requires a finite lower bound")]
    ProjectiveResolutionNeedsFiniteLowerBound,
    /// The nonfinite-upper branch requires a finite highest known degree.
    #[error("ProjectiveResolutionOfComplex requires a finite highest known degree")]
    ProjectiveResolutionNeedsFiniteHighestKnownDegree,
    /// A nonnegative GAP loop bound did not fit a Rust address-sized index.
    #[error("ProjectiveResolutionOfComplex loop bound is not addressable: {value}")]
    ProjectiveResolutionLoopBoundNotAddressable {
        /// GAP integer that could not be used as a list/range index.
        value: Integer,
    },
    /// The two maps passed to QPA's `PullBack` unexpectedly had different
    /// ranges inside the resolution construction.
    #[error("PullBack returned fail inside ProjectiveResolutionOfComplex")]
    ProjectiveResolutionPullBackReturnedFail,
    /// The conversion method performs integer arithmetic on the lower bound.
    #[error("ProjectiveToInjectiveComplex requires a finite lower bound")]
    ProjectiveToInjectiveNeedsFiniteLowerBound,
    /// The explicitly finite conversion method requires a finite upper bound.
    #[error("ProjectiveToInjectiveFiniteComplex requires a finite upper bound")]
    ProjectiveToInjectiveFiniteNeedsFiniteUpperBound,
    /// QPA's range loop cannot use a nonfinite highest-known degree.
    #[error("ProjectiveToInjectiveComplex requires a finite highest known degree")]
    ProjectiveToInjectiveNeedsFiniteHighestKnownDegree,
    /// QPA indexed the first entry of an empty projective description.
    #[error("{operation} indexed an empty projective description")]
    EmptyProjectiveDescription {
        /// Source operation that indexed the empty list.
        operation: &'static str,
    },
    /// A dimension-vector entry was reused as an invalid one-based vertex.
    #[error("projective description entry {vertex} lies outside 1..={vertex_count}")]
    ProjectiveDescriptionVertexOutOfRange {
        /// Entry interpreted by QPA as a vertex number.
        vertex: usize,
        /// Number of available indecomposable modules.
        vertex_count: usize,
    },
    /// A star-helper description did not match its component map list.
    #[error("star component position {position} lies outside 1..={component_count}")]
    StarComponentPositionOutOfRange {
        /// One-based requested component.
        position: usize,
        /// Available component maps.
        component_count: usize,
    },
    /// The selected projective basis and module coordinates did not agree.
    #[error(
        "star projective basis has {basis} elements but the module component has {coordinates} coordinates"
    )]
    StarProjectiveBasisCoordinateMismatch {
        /// Number of basis elements.
        basis: usize,
        /// Number of coordinates supplied by the module element.
        coordinates: usize,
    },
    /// QPA's filtered opposite-projective generator list was empty.
    #[error("StarOfMapBetweenIndecProjectives found no opposite generator at vertex {vertex}")]
    StarOppositeGeneratorMissing {
        /// Expected one-based supporting vertex.
        vertex: usize,
    },
    /// QPA's internal `HomFromProjective` unexpectedly returned `fail`.
    #[error("HomFromProjective returned fail inside StarOfMapBetweenIndecProjectives")]
    StarHomFromProjectiveReturnedFail,
}

/// QPA's four-element result from `CheckForIsomorphism`.
#[derive(Clone, Debug)]
pub struct ProjectiveIsomorphismComponent<F: Field> {
    map: PathAlgebraModuleHomomorphism<F>,
    vertex: usize,
    source_position: usize,
    range_position: usize,
}

impl<F: Field> ProjectiveIsomorphismComponent<F> {
    /// The isomorphism between the selected indecomposable summands.
    #[must_use]
    pub const fn map(&self) -> &PathAlgebraModuleHomomorphism<F> {
        &self.map
    }

    /// One-based vertex describing the common projective summand.
    #[must_use]
    pub const fn vertex(&self) -> usize {
        self.vertex
    }

    /// One-based position in the differential source.
    #[must_use]
    pub const fn source_position(&self) -> usize {
        self.source_position
    }

    /// One-based position in the differential range.
    #[must_use]
    pub const fn range_position(&self) -> usize {
        self.range_position
    }
}

fn finite_bound(value: Option<ExtendedInteger>) -> Option<Integer> {
    match value {
        Some(ExtendedInteger::Finite(value)) => Some(value),
        _ => None,
    }
}

fn integer_interval(first: &Integer, last: &Integer) -> Vec<Integer> {
    let mut result = Vec::new();
    let mut degree = first.clone();
    while degree <= *last {
        result.push(degree.clone());
        degree += Integer::from(1);
    }
    result
}

/// Translate QPA's cached `IsProjectiveComplex` property.
///
/// The source returns `false` without inspecting objects unless both bounds
/// are finite.  Otherwise it tests every object in the closed degree range.
///
/// # Errors
///
/// Propagates complex bounds/object lookup and module-projectivity failures.
pub fn is_projective_complex<F>(complex: &ModuleComplex<F>) -> Result<bool, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    if let Some(cached) = complex.cached_is_projective_complex() {
        return Ok(cached);
    }
    let Some(lower) = finite_bound(complex.lower_bound()?) else {
        return Ok(complex.cache_is_projective_complex(false));
    };
    let Some(upper) = finite_bound(complex.upper_bound()?) else {
        return Ok(complex.cache_is_projective_complex(false));
    };
    for degree in integer_interval(&lower, &upper) {
        if !is_projective_module(&complex.object(&degree)?)? {
            return Ok(complex.cache_is_projective_complex(false));
        }
    }
    Ok(complex.cache_is_projective_complex(true))
}

/// Translate QPA's cached `IsInjectiveComplex` property.
///
/// # Errors
///
/// Propagates complex bounds/object lookup and module-injectivity failures.
pub fn is_injective_complex<F>(complex: &ModuleComplex<F>) -> Result<bool, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    if let Some(cached) = complex.cached_is_injective_complex() {
        return Ok(cached);
    }
    let Some(lower) = finite_bound(complex.lower_bound()?) else {
        return Ok(complex.cache_is_injective_complex(false));
    };
    let Some(upper) = finite_bound(complex.upper_bound()?) else {
        return Ok(complex.cache_is_injective_complex(false));
    };
    for degree in integer_interval(&lower, &upper) {
        if !is_injective_module(&complex.object(&degree)?)? {
            return Ok(complex.cache_is_injective_complex(false));
        }
    }
    Ok(complex.cache_is_injective_complex(true))
}

/// Translate `CompareWithIndecProjective`.
///
/// `None` is GAP's boolean `false`.  A successful value is the one-based
/// vertex position of the unique simple in the module top.
///
/// # Errors
///
/// Propagates projectivity and top construction.
pub fn compare_with_indec_projective<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<Option<usize>, DerivedCategoryError> {
    if !is_projective_module(module)? {
        return Ok(None);
    }
    let top = top_of_module(module.clone())?;
    if top.dimension_vector().iter().sum::<usize>() > 1 {
        return Ok(None);
    }
    Ok(top
        .dimension_vector()
        .iter()
        .position(|&multiplicity| multiplicity == 1)
        .map(|position| position + 1))
}

/// Translate `CompareWithIndecInjective`.
///
/// `None` is GAP's boolean `false`; a successful value is one-based.
///
/// # Errors
///
/// Propagates injectivity and socle construction.
pub fn compare_with_indec_injective<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<Option<usize>, DerivedCategoryError> {
    if !is_injective_module(module)? {
        return Ok(None);
    }
    let socle = socle_of_module(module.clone())?;
    if socle.dimension_vector().iter().sum::<usize>() > 1 {
        return Ok(None);
    }
    Ok(socle
        .dimension_vector()
        .iter()
        .position(|&multiplicity| multiplicity == 1)
        .map(|position| position + 1))
}

/// Translate `DescriptionOfProjectiveModule`.
///
/// Entries retain QPA's possibility of a boolean `false` from
/// `CompareWithIndecProjective` as `None`.
///
/// # Errors
///
/// Propagates projectivity and top construction for each summand.
pub fn description_of_projective_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<Vec<Option<usize>>, DerivedCategoryError> {
    if module.dimension() == 0 {
        return Ok(Vec::new());
    }
    if let Some(summands) = module.direct_sum_summands() {
        return summands.iter().map(compare_with_indec_projective).collect();
    }
    Ok(vec![compare_with_indec_projective(module)?])
}

/// Translate `DescriptionOfProjOrInjComplexInDegree`.
///
/// QPA returns the top dimension vector for projectives and the socle
/// dimension vector for injectives.  A zero object is represented by `[]`.
///
/// # Errors
///
/// Propagates lazy object lookup and top/socle construction.
pub fn description_of_proj_or_inj_complex_in_degree<F>(
    complex: &ModuleComplex<F>,
    degree: &Integer,
    projective: bool,
) -> Result<Vec<usize>, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    let object = complex.object(degree)?;
    if object.dimension() == 0 {
        return Ok(Vec::new());
    }
    let semisimple = if projective {
        top_of_module(object)?
    } else {
        socle_of_module(object)?
    };
    Ok(semisimple.dimension_vector().to_vec())
}

/// Translate `DescriptionOfProjComplexInDegree`.
///
/// # Errors
///
/// Propagates [`description_of_proj_or_inj_complex_in_degree`].
pub fn description_of_proj_complex_in_degree<F>(
    complex: &ModuleComplex<F>,
    degree: &Integer,
) -> Result<Vec<usize>, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    description_of_proj_or_inj_complex_in_degree(complex, degree, true)
}

/// Translate `DescriptionOfInjComplexInDegree`.
///
/// # Errors
///
/// Propagates [`description_of_proj_or_inj_complex_in_degree`].
pub fn description_of_inj_complex_in_degree<F>(
    complex: &ModuleComplex<F>,
    degree: &Integer,
) -> Result<Vec<usize>, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    description_of_proj_or_inj_complex_in_degree(complex, degree, false)
}

/// Translate `DescriptionOfFiniteProjOrInjComplex`.
///
/// # Errors
///
/// Returns [`DerivedCategoryError::NonFiniteComplex`] unless both bounds are
/// finite, then propagates every per-degree description failure.
pub fn description_of_finite_proj_or_inj_complex<F>(
    complex: &ModuleComplex<F>,
    projective: bool,
) -> Result<Vec<Vec<usize>>, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    let lower =
        finite_bound(complex.lower_bound()?).ok_or(DerivedCategoryError::NonFiniteComplex)?;
    let upper =
        finite_bound(complex.upper_bound()?).ok_or(DerivedCategoryError::NonFiniteComplex)?;
    integer_interval(&lower, &upper)
        .iter()
        .map(|degree| description_of_proj_or_inj_complex_in_degree(complex, degree, projective))
        .collect()
}

/// Translate `DescriptionOfFiniteProjComplex`.
///
/// # Errors
///
/// Propagates [`description_of_finite_proj_or_inj_complex`].
pub fn description_of_finite_proj_complex<F>(
    complex: &ModuleComplex<F>,
) -> Result<Vec<Vec<usize>>, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    description_of_finite_proj_or_inj_complex(complex, true)
}

/// Translate `DescriptionOfFiniteInjComplex`.
///
/// # Errors
///
/// Propagates [`description_of_finite_proj_or_inj_complex`].
pub fn description_of_finite_inj_complex<F>(
    complex: &ModuleComplex<F>,
) -> Result<Vec<Vec<usize>>, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    description_of_finite_proj_or_inj_complex(complex, false)
}

fn qpa_proj_or_inj_complex_string<F>(
    complex: &ModuleComplex<F>,
    projective: bool,
) -> Result<String, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    let upper = complex.upper_bound()?;
    let finite = matches!(upper, Some(ExtendedInteger::Finite(_)));
    let mut degree = if finite {
        finite_bound(upper).ok_or(DerivedCategoryError::NonFiniteComplex)?
    } else {
        finite_bound(complex.highest_known_degree()?)
            .ok_or(DerivedCategoryError::PrinterNeedsFiniteHighestKnownDegree)?
    };
    let lower = finite_bound(complex.lower_bound()?)
        .ok_or(DerivedCategoryError::PrinterNeedsFiniteLowerBound)?;
    let descriptions = if finite {
        let mut values = description_of_finite_proj_or_inj_complex(complex, projective)?;
        values.reverse();
        values
    } else {
        let mut values = integer_interval(&lower, &degree)
            .iter()
            .map(|position| {
                description_of_proj_or_inj_complex_in_degree(complex, position, projective)
            })
            .collect::<Result<Vec<_>, _>>()?;
        values.reverse();
        values
    };
    let symbol = if projective { 'P' } else { 'I' };
    let mut output = if finite {
        "0 ->".to_owned()
    } else {
        "--- ->".to_owned()
    };
    for description in descriptions {
        if description.is_empty() {
            continue;
        }
        write!(output, " {degree}: ").expect("writing to a String cannot fail");
        let mut first = true;
        for (vertex, multiplicity) in description.iter().copied().enumerate() {
            if multiplicity == 0 {
                continue;
            }
            if !first {
                output.push_str(" + ");
            }
            if multiplicity != 1 {
                write!(output, "{multiplicity}").expect("writing to a String cannot fail");
            }
            write!(output, "{symbol}{}", vertex + 1).expect("writing to a String cannot fail");
            first = false;
        }
        output.push_str(" ->");
        degree -= Integer::from(1);
    }
    output.push_str(" 0");
    Ok(output)
}

/// Return QPA's specialized `PrintObj` rendering for a projective complex.
///
/// # Errors
///
/// Propagates bounds and projective-description failures.
pub fn qpa_projective_complex_string<F>(
    complex: &ModuleComplex<F>,
) -> Result<String, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    qpa_proj_or_inj_complex_string(complex, true)
}

/// Return QPA's specialized `PrintObj` rendering for an injective complex.
///
/// # Errors
///
/// Propagates bounds and injective-description failures.
pub fn qpa_injective_complex_string<F>(
    complex: &ModuleComplex<F>,
) -> Result<String, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    qpa_proj_or_inj_complex_string(complex, false)
}

/// Translate `FindAllMapComponents` in source-summand-major order.
///
/// # Errors
///
/// Requires direct-sum metadata on both modules and propagates map
/// composition failures.
pub fn find_all_map_components<F: Field>(
    map: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Vec<Vec<PathAlgebraModuleHomomorphism<F>>>, DerivedCategoryError> {
    let inclusions = direct_sum_inclusions(map.source())?;
    let projections = direct_sum_projections(map.range())?;
    inclusions
        .iter()
        .map(|inclusion| {
            let first_part = inclusion.compose(map)?;
            projections
                .iter()
                .map(|projection| Ok(first_part.compose(projection)?))
                .collect()
        })
        .collect()
}

/// Translate `CheckForIsomorphism(d)`.
///
/// QPA intersects the projective-vertex descriptions, then scans common
/// vertices, source positions, and range positions in that order.  The first
/// component that is an isomorphism is returned; `None` is its empty list.
///
/// # Errors
///
/// Propagates projective descriptions and component extraction, and rejects a
/// `false` summand description in the method's intended projective domain.
pub fn check_for_isomorphism<F: Field>(
    differential: &PathAlgebraModuleHomomorphism<F>,
) -> Result<Option<ProjectiveIsomorphismComponent<F>>, DerivedCategoryError> {
    let source_description = description_of_projective_module(differential.source())?
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or(DerivedCategoryError::NonProjectiveSummandDescription)?;
    let range_description = description_of_projective_module(differential.range())?
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or(DerivedCategoryError::NonProjectiveSummandDescription)?;
    let mut mutual = source_description
        .iter()
        .copied()
        .filter(|vertex| range_description.contains(vertex))
        .collect::<Vec<_>>();
    mutual.sort_unstable();
    mutual.dedup();
    if mutual.is_empty() {
        return Ok(None);
    }
    let components = find_all_map_components(differential)?;
    for vertex in mutual {
        let source_positions = source_description
            .iter()
            .enumerate()
            .filter_map(|(position, &entry)| (entry == vertex).then_some(position));
        let range_positions = range_description
            .iter()
            .enumerate()
            .filter_map(|(position, &entry)| (entry == vertex).then_some(position))
            .collect::<Vec<_>>();
        for source_position in source_positions {
            for &range_position in &range_positions {
                let map = components
                    .get(source_position)
                    .and_then(|row| row.get(range_position))
                    .ok_or(DerivedCategoryError::NonRectangularMapComponentGrid)?;
                if map.is_isomorphism() {
                    return Ok(Some(ProjectiveIsomorphismComponent {
                        map: map.clone(),
                        vertex,
                        source_position: source_position + 1,
                        range_position: range_position + 1,
                    }));
                }
            }
        }
    }
    Ok(None)
}

fn cancellation_indices<F: Field>(
    components: &[Vec<PathAlgebraModuleHomomorphism<F>>],
    cancellation: &ProjectiveIsomorphismComponent<F>,
) -> Result<(usize, usize), DerivedCategoryError> {
    let source_count = components.len();
    let range_count = components.first().map_or(0, Vec::len);
    if components.iter().any(|row| row.len() != range_count) {
        return Err(DerivedCategoryError::NonRectangularMapComponentGrid);
    }
    let source_position = cancellation.source_position.checked_sub(1);
    let range_position = cancellation.range_position.checked_sub(1);
    let (Some(source_position), Some(range_position)) = (source_position, range_position) else {
        return Err(DerivedCategoryError::CancellationPositionOutOfRange {
            source_position: cancellation.source_position,
            range_position: cancellation.range_position,
            source_count,
            range_count,
        });
    };
    if source_position >= source_count || range_position >= range_count {
        return Err(DerivedCategoryError::CancellationPositionOutOfRange {
            source_position: cancellation.source_position,
            range_position: cancellation.range_position,
            source_count,
            range_count,
        });
    }
    Ok((source_position, range_position))
}

fn map_from_source_components<F: Field>(
    source: &PathAlgebraModule<F>,
    range: &PathAlgebraModule<F>,
    components: &[Vec<PathAlgebraModuleHomomorphism<F>>],
) -> Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError> {
    let projections = direct_sum_projections(source)?;
    let inclusions = direct_sum_inclusions(range)?;
    if components.len() != projections.len()
        || components.iter().any(|row| row.len() != inclusions.len())
    {
        return Err(DerivedCategoryError::NonRectangularMapComponentGrid);
    }
    let mut terms = components
        .iter()
        .enumerate()
        .flat_map(|(source_position, row)| {
            let projections = &projections;
            let inclusions = &inclusions;
            row.iter()
                .enumerate()
                .map(move |(range_position, component)| {
                    projections[source_position]
                        .compose(component)?
                        .compose(&inclusions[range_position])
                })
        });
    let Some(first) = terms.next() else {
        return Err(DerivedCategoryError::EmptyMapComponentList);
    };
    let mut sum = first?;
    for term in terms {
        sum = sum.add(&term?)?;
    }
    Ok(sum)
}

/// Translate `ReduceDifferential(d,list,reducedModules)`.
///
/// After deleting the selected source row and range column, QPA subtracts
/// `left * inverse(selected) * right` from every surviving component.  This
/// is the literal map-valued Schur complement in the source method.
///
/// # Errors
///
/// Propagates component extraction, exact inverse, composition, addition, and
/// checked direct-sum map construction.
pub fn reduce_differential<F: Field>(
    differential: &PathAlgebraModuleHomomorphism<F>,
    cancellation: &ProjectiveIsomorphismComponent<F>,
    reduced_source: &PathAlgebraModule<F>,
    reduced_range: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError> {
    if reduced_source.dimension() == 0 || reduced_range.dimension() == 0 {
        return Ok(zero_module_mapping(
            reduced_source.clone(),
            reduced_range.clone(),
        )?);
    }
    let components = find_all_map_components(differential)?;
    let (selected_source, selected_range) = cancellation_indices(&components, cancellation)?;
    let inverse = inverse_of_isomorphism(&components[selected_source][selected_range])?;
    let source_positions = (0..components.len())
        .filter(|&position| position != selected_source)
        .collect::<Vec<_>>();
    let range_positions = (0..components[0].len())
        .filter(|&position| position != selected_range)
        .collect::<Vec<_>>();
    let mut reduced_components = Vec::with_capacity(source_positions.len());
    for &source_position in &source_positions {
        let mut row = Vec::with_capacity(range_positions.len());
        for &range_position in &range_positions {
            let correction = components[source_position][selected_range]
                .compose(&inverse)?
                .compose(&components[selected_source][range_position])?;
            row.push(components[source_position][range_position].add(&correction.negate()?)?);
        }
        reduced_components.push(row);
    }
    map_from_source_components(reduced_source, reduced_range, &reduced_components)
}

/// Translate `ReducePreviousDifferential(d,inclusion)`.
///
/// # Errors
///
/// Propagates checked map composition.
pub fn reduce_previous_differential<F: Field>(
    differential: &PathAlgebraModuleHomomorphism<F>,
    inclusion: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError> {
    Ok(inclusion.compose(differential)?)
}

/// Translate `ReduceTMap(t,d,list,incl)`.
///
/// This applies the same cancellation row operation to the quasi-isomorphism
/// component `t`: restrict to the reduced source, then subtract the route
/// through the cancelled range summand, its inverse, and the removed source
/// component of `t`.
///
/// # Errors
///
/// Propagates direct-sum component extraction, inverse, composition, addition,
/// and cancellation-position failures.
pub fn reduce_t_map<F: Field>(
    map: &PathAlgebraModuleHomomorphism<F>,
    differential: &PathAlgebraModuleHomomorphism<F>,
    cancellation: &ProjectiveIsomorphismComponent<F>,
    inclusion: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError> {
    let restricted = inclusion.compose(map)?;
    let original_components = find_all_map_components(differential)?;
    let (selected_source, selected_range) =
        cancellation_indices(&original_components, cancellation)?;
    let inverse = inverse_of_isomorphism(&original_components[selected_source][selected_range])?;
    let source_inclusions = direct_sum_inclusions(differential.source())?;
    let removed_map = source_inclusions[selected_source].compose(map)?;
    let reduced_differential = inclusion.compose(differential)?;
    let reduced_components = find_all_map_components(&reduced_differential)?;
    let reduced_projections = direct_sum_projections(inclusion.source())?;
    if reduced_components.len() != reduced_projections.len()
        || reduced_components
            .iter()
            .any(|row| selected_range >= row.len())
    {
        return Err(DerivedCategoryError::NonRectangularMapComponentGrid);
    }
    let mut terms = reduced_projections
        .iter()
        .zip(&reduced_components)
        .map(|(projection, row)| {
            projection
                .compose(&row[selected_range])?
                .compose(&inverse)?
                .compose(&removed_map)
        });
    let Some(first) = terms.next() else {
        return Ok(restricted);
    };
    let mut correction = first?;
    for term in terms {
        correction = correction.add(&term?)?;
    }
    Ok(restricted.add(&correction.negate()?)?)
}

/// Translate `ReducePreviousTMap(t,inclusion)`.
///
/// # Errors
///
/// Propagates checked map composition.
pub fn reduce_previous_t_map<F: Field>(
    map: &PathAlgebraModuleHomomorphism<F>,
    inclusion: &PathAlgebraModuleHomomorphism<F>,
) -> Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError> {
    Ok(inclusion.compose(map)?)
}

fn cache_projective_complex<F>(complex: ModuleComplex<F>) -> ModuleComplex<F>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    let _ = complex.cache_is_projective_complex(true);
    complex
}

fn cache_injective_complex<F>(complex: ModuleComplex<F>) -> ModuleComplex<F>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    let _ = complex.cache_is_injective_complex(true);
    complex
}

/// Translate QPA's `ProjectiveResolutionOfComplex`.
///
/// This is the source pullback/projective-cover construction, including its
/// repeated cancellation of isomorphic projective summands, stalk shortcut,
/// finite early exits, and the final syzygy/projective-resolution Yoneda
/// splice.  The source's documented one-cancellation-only overlap behavior is
/// retained.
///
/// # Errors
///
/// Propagates all complex, projective-cover, pullback, kernel, direct-sum,
/// cancellation, and infinite-list failures.  GAP infinity arithmetic and an
/// unaddressable loop range are represented by typed errors.
#[allow(clippy::too_many_lines)]
pub fn projective_resolution_of_complex<F>(
    complex: &ModuleComplex<F>,
) -> Result<ModuleComplex<F>, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    if is_projective_complex(complex)? {
        return Ok(complex.clone());
    }

    let category = complex.category().clone();
    let lower = finite_bound(complex.lower_bound()?)
        .ok_or(DerivedCategoryError::ProjectiveResolutionNeedsFiniteLowerBound)?;
    let terminal = if let Some(upper) = finite_bound(complex.upper_bound()?) {
        upper
    } else {
        finite_bound(complex.highest_known_degree()?)
            .ok_or(DerivedCategoryError::ProjectiveResolutionNeedsFiniteHighestKnownDegree)?
    };
    let steps_integer = &terminal - &lower + Integer::from(1);
    let steps = steps_integer.to_string().parse::<usize>().map_err(|_| {
        DerivedCategoryError::ProjectiveResolutionLoopBoundNotAddressable {
            value: steps_integer.clone(),
        }
    })?;

    if matches!(
        complex.length()?,
        Some(ExtendedInteger::Finite(ref length)) if *length == 1
    ) {
        let resolution = projective_resolution_complex(&complex.object(&lower)?)?
            .brutal_truncation_below(&Integer::from(0))?
            .shift_unsigned(&(-lower))?;
        return Ok(cache_projective_complex(resolution));
    }

    let first_comparison = projective_cover(&complex.object(&lower)?)?;
    let first_pullback = pull_back(
        &complex.differential(&(&lower + Integer::from(1)))?,
        &first_comparison,
    )?
    .ok_or(DerivedCategoryError::ProjectiveResolutionPullBackReturnedFail)?;
    let first_cover = projective_cover(first_pullback.0.source())?;
    let mut differentials = vec![first_cover.compose(&first_pullback.0)?];
    let mut comparisons = vec![first_comparison, first_cover.compose(&first_pullback.1)?];

    while let Some(cancellation) = check_for_isomorphism(&differentials[0])? {
        let source_reduction = direct_sum_minus_summands(
            differentials[0].source(),
            &[cancellation.source_position()],
        )?;
        let range_reduction =
            direct_sum_minus_summands(differentials[0].range(), &[cancellation.range_position()])?;
        let reduced = reduce_differential(
            &differentials[0],
            &cancellation,
            source_reduction.0.source(),
            range_reduction.0.source(),
        )?;
        comparisons[1] = reduce_t_map(
            &comparisons[1],
            &differentials[0],
            &cancellation,
            &source_reduction.0,
        )?;
        comparisons[0] = reduce_previous_t_map(&comparisons[0], &range_reduction.0)?;
        differentials[0] = reduced;
    }

    for position in 0..steps {
        let kernel = kernel_inclusion(&differentials[position])?;
        let comparison = kernel.compose(&comparisons[position + 1])?;
        let object_degree = &lower + Integer::from(position + 2);
        if kernel.source().dimension() == 0 && complex.object(&object_degree)?.dimension() == 0 {
            let finite =
                finite_complex(category.clone(), &lower + Integer::from(1), differentials)?;
            return Ok(cache_projective_complex(finite));
        }

        let pullback = pull_back(&complex.differential(&object_degree)?, &comparison)?
            .ok_or(DerivedCategoryError::ProjectiveResolutionPullBackReturnedFail)?;
        let cover = projective_cover(pullback.0.source())?;
        if cover.source().dimension() == 0 {
            let finite =
                finite_complex(category.clone(), &lower + Integer::from(1), differentials)?;
            return Ok(cache_projective_complex(finite));
        }

        differentials.push(cover.compose(&pullback.0)?.compose(&kernel)?);
        comparisons.push(cover.compose(&pullback.1)?);
        let differential_position = position + 1;
        while let Some(cancellation) = check_for_isomorphism(&differentials[differential_position])?
        {
            let source_reduction = direct_sum_minus_summands(
                differentials[differential_position].source(),
                &[cancellation.source_position()],
            )?;
            let range_reduction = direct_sum_minus_summands(
                differentials[differential_position].range(),
                &[cancellation.range_position()],
            )?;
            let reduced = reduce_differential(
                &differentials[differential_position],
                &cancellation,
                source_reduction.0.source(),
                range_reduction.0.source(),
            )?;
            comparisons[differential_position + 1] = reduce_t_map(
                &comparisons[differential_position + 1],
                &differentials[differential_position],
                &cancellation,
                &source_reduction.0,
            )?;
            comparisons[differential_position] =
                reduce_previous_t_map(&comparisons[differential_position], &range_reduction.0)?;
            differentials[differential_position] = reduced;
            differentials[differential_position - 1] = reduce_previous_differential(
                &differentials[differential_position - 1],
                &range_reduction.0,
            )?;
        }
    }

    let steps_degree = Integer::from(steps);
    let temporary = finite_complex(category.clone(), &lower + Integer::from(1), differentials)?
        .syzygy_truncation(&(&lower + &steps_degree + Integer::from(1)))?;
    let tail = projective_resolution_complex(
        &temporary.object(&(&lower + &steps_degree + Integer::from(2)))?,
    )?;
    let mut result = yoneda_product(&tail, &temporary)?;

    let overlap_degree = &lower + &steps_degree + Integer::from(2);
    if let Some(cancellation) = check_for_isomorphism(&result.differential(&overlap_degree)?)? {
        let differential = result.differential(&overlap_degree)?;
        let source_reduction =
            direct_sum_minus_summands(differential.source(), &[cancellation.source_position()])?;
        let range_reduction =
            direct_sum_minus_summands(differential.range(), &[cancellation.range_position()])?;
        let reduced = reduce_differential(
            &differential,
            &cancellation,
            source_reduction.0.source(),
            range_reduction.0.source(),
        )?;
        let previous_degree = &lower + &steps_degree + Integer::from(1);
        let previous = reduce_previous_differential(
            &result.differential(&previous_degree)?,
            &range_reduction.0,
        )?;
        let differential_list = result.differentials();
        let positive =
            differential_list.positive_part_from(&(&lower + &steps_degree + Integer::from(3)))?;
        let replacement = InfiniteList::finite(previous_degree, vec![previous, reduced]);
        let negative = differential_list.negative_part_from(&(&lower + &steps_degree))?;
        result = result.with_differential_list(
            category,
            InfiniteList::concatenate(&[positive, replacement, negative])?,
        );
    }
    Ok(cache_projective_complex(result))
}

/// Translate `MatricesOfDualMap` by transposing every stored vertex matrix.
#[must_use]
pub fn matrices_of_dual_map<F: Field>(
    map: &PathAlgebraModuleHomomorphism<F>,
) -> Vec<DenseMatrix<F>> {
    map.vertex_matrices()
        .iter()
        .map(DenseMatrix::transpose)
        .collect()
}

fn modules_from_projective_description<F: Field>(
    modules: &[PathAlgebraModule<F>],
    description: &[usize],
    operation: &'static str,
) -> Result<PathAlgebraModule<F>, DerivedCategoryError> {
    if description.is_empty() {
        return Err(DerivedCategoryError::EmptyProjectiveDescription { operation });
    }
    let summands = description
        .iter()
        .map(|&vertex| {
            modules
                .get(vertex.checked_sub(1).ok_or(
                    DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange {
                        vertex,
                        vertex_count: modules.len(),
                    },
                )?)
                .cloned()
                .ok_or(
                    DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange {
                        vertex,
                        vertex_count: modules.len(),
                    },
                )
        })
        .collect::<Result<Vec<_>, _>>()?;
    direct_sum_of_qpa_modules(&summands)?
        .ok_or(DerivedCategoryError::EmptyProjectiveDescription { operation })
}

fn sum_nonempty_maps<F: Field>(
    maps: impl IntoIterator<Item = Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError>>,
) -> Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError> {
    let mut maps = maps.into_iter();
    let Some(first) = maps.next() else {
        return Err(DerivedCategoryError::EmptyMapComponentList);
    };
    let mut sum = first?;
    for map in maps {
        sum = sum.add(&map?)?;
    }
    Ok(sum)
}

fn opposite_quiver_algebra_element<F: Field>(
    element: &QuiverAlgebraElement<F>,
) -> Result<QuiverAlgebraElement<F>, DerivedCategoryError> {
    match element {
        QuiverAlgebraElement::Path(element) => Ok(QuiverAlgebraElement::Path(element.opposite()?)),
        QuiverAlgebraElement::Quotient(element) => Ok(QuiverAlgebraElement::Quotient(
            opposite_quotient_element(element)?,
        )),
    }
}

fn linear_combination_of_projective_basis<F: Field>(
    algebra: &QuiverAlgebra<F>,
    basis: &[QuiverAlgebraElement<F>],
    coefficients: &[F::Element],
) -> Result<QuiverAlgebraElement<F>, DerivedCategoryError> {
    if basis.len() != coefficients.len() {
        return Err(
            DerivedCategoryError::StarProjectiveBasisCoordinateMismatch {
                basis: basis.len(),
                coordinates: coefficients.len(),
            },
        );
    }
    let mut result = algebra.zero();
    for (element, coefficient) in basis.iter().zip(coefficients) {
        result = result.add(&element.scale_left(coefficient)?)?;
    }
    Ok(result)
}

/// Translate `StarOfMapBetweenIndecProjectives`.
///
/// The method reconstructs the image of the distinguished generator inside
/// `BasisOfProjectives`, reverses that algebra element, acts on the matching
/// opposite-projective generator, and combines the resulting component maps
/// in QPA's entered summand order.
///
/// # Errors
///
/// Propagates projective bases, opposite-algebra arithmetic, module actions,
/// direct sums, component maps, and the source method's unchecked indexing.
#[allow(clippy::too_many_lines)]
pub fn star_of_map_between_indec_projectives<F: Field>(
    map: &PathAlgebraModuleHomomorphism<F>,
    source_vertex: usize,
    range_description: &[usize],
) -> Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError> {
    let algebra = map.source().right_acting_algebra();
    let projective_bases =
        basis_of_projectives(algebra)?.ok_or(EnvelopingModuleError::ProjectiveBasisUnavailable)?;
    let opposite = opposite_module_algebra(algebra)?;
    let opposite_projectives = indec_projective_modules(&opposite)?;
    let opposite_range = opposite_projectives
        .get(source_vertex.checked_sub(1).ok_or(
            DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange {
                vertex: source_vertex,
                vertex_count: opposite_projectives.len(),
            },
        )?)
        .ok_or(
            DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange {
                vertex: source_vertex,
                vertex_count: opposite_projectives.len(),
            },
        )?;

    let summands = if map.range().is_direct_sum_of_modules() {
        direct_sum_projections(map.range())?
            .iter()
            .map(|projection| Ok(map.compose(projection)?))
            .collect::<Result<Vec<_>, DerivedCategoryError>>()?
    } else {
        vec![map.clone()]
    };
    let mut opposite_maps = Vec::with_capacity(range_description.len());
    for (position, &range_vertex) in range_description.iter().enumerate() {
        let summand = summands.get(position).ok_or(
            DerivedCategoryError::StarComponentPositionOutOfRange {
                position: position + 1,
                component_count: summands.len(),
            },
        )?;
        let opposite_source = opposite_projectives
            .get(range_vertex.checked_sub(1).ok_or(
                DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange {
                    vertex: range_vertex,
                    vertex_count: opposite_projectives.len(),
                },
            )?)
            .ok_or(
                DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange {
                    vertex: range_vertex,
                    vertex_count: opposite_projectives.len(),
                },
            )?;
        if summand.is_zero() {
            opposite_maps.push(zero_module_mapping(
                opposite_source.clone(),
                opposite_range.clone(),
            )?);
            continue;
        }

        let generator = minimal_generating_set_of_module(summand.source().clone())?
            .into_iter()
            .next()
            .ok_or(DerivedCategoryError::StarOppositeGeneratorMissing {
                vertex: source_vertex,
            })?;
        let image = summand.image(&generator)?;
        let basis = projective_bases
            .get(range_vertex.checked_sub(1).ok_or(
                DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange {
                    vertex: range_vertex,
                    vertex_count: projective_bases.len(),
                },
            )?)
            .and_then(|components| components.get(source_vertex.wrapping_sub(1)))
            .ok_or(
                DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange {
                    vertex: source_vertex,
                    vertex_count: projective_bases.len(),
                },
            )?;
        let coordinates = image
            .external_representation()
            .get(source_vertex.wrapping_sub(1))
            .ok_or(
                DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange {
                    vertex: source_vertex,
                    vertex_count: image.external_representation().len(),
                },
            )?;
        let algebra_element = linear_combination_of_projective_basis(algebra, basis, coordinates)?;
        let opposite_element = opposite_quiver_algebra_element(&algebra_element)?;
        let opposite_generator = minimal_generating_set_of_module(opposite_range.clone())?
            .into_iter()
            .find(|element| {
                let vertex = source_vertex.wrapping_sub(1);
                let dimension = opposite_range
                    .dimension_vector()
                    .get(vertex)
                    .copied()
                    .unwrap_or(0);
                element
                    .external_representation()
                    .get(vertex)
                    .is_some_and(|row| {
                        row.iter()
                            .take(dimension)
                            .any(|entry| !opposite.field().is_zero(entry))
                    })
            })
            .ok_or(DerivedCategoryError::StarOppositeGeneratorMissing {
                vertex: source_vertex,
            })?;
        let acted = opposite_generator.act(&opposite_element)?;
        let support_vertex = opposite
            .quiver()
            .vertices()
            .nth(range_vertex.wrapping_sub(1))
            .ok_or(
                DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange {
                    vertex: range_vertex,
                    vertex_count: opposite.quiver().number_of_vertices(),
                },
            )?;
        let idempotent = include_in_path_algebra(&support_vertex.path(), &opposite)?;
        let support = acted.act(&idempotent)?;
        let component = if support.is_zero() {
            zero_module_mapping(opposite_source.clone(), opposite_range.clone())?
        } else {
            hom_from_projective(&support, opposite_range)?
                .ok_or(DerivedCategoryError::StarHomFromProjectiveReturnedFail)?
        };
        opposite_maps.push(component);
    }

    let sources = opposite_maps
        .iter()
        .map(|component| component.source().clone())
        .collect::<Vec<_>>();
    let source = direct_sum_of_qpa_modules(&sources)?.ok_or(
        DerivedCategoryError::EmptyProjectiveDescription {
            operation: "StarOfMapBetweenIndecProjectives",
        },
    )?;
    let projections = direct_sum_projections(&source)?;
    sum_nonempty_maps(
        projections
            .iter()
            .zip(&opposite_maps)
            .map(|(projection, component)| Ok(projection.compose(component)?)),
    )
}

/// Translate `StarOfMapBetweenDecompProjectives`.
///
/// # Errors
///
/// Propagates summand extraction, each indecomposable star construction,
/// QPA's source-identification bridge, and direct-sum map assembly.
pub fn star_of_map_between_decomp_projectives<F: Field>(
    map: &PathAlgebraModuleHomomorphism<F>,
    source_description: &[usize],
    range_description: &[usize],
) -> Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError> {
    let source_inclusions = direct_sum_inclusions(map.source())?;
    let mut maps = Vec::with_capacity(source_description.len());
    for (position, &vertex) in source_description.iter().enumerate() {
        let inclusion = source_inclusions.get(position).ok_or(
            DerivedCategoryError::StarComponentPositionOutOfRange {
                position: position + 1,
                component_count: source_inclusions.len(),
            },
        )?;
        maps.push(star_of_map_between_indec_projectives(
            &inclusion.compose(map)?,
            vertex,
            range_description,
        )?);
    }
    if let Some(first) = maps.first().cloned() {
        for map in maps.iter_mut().skip(1) {
            let projections = direct_sum_projections(first.source())?;
            let inclusions = direct_sum_inclusions(map.source())?;
            let bridge = sum_nonempty_maps(
                projections
                    .iter()
                    .zip(&inclusions)
                    .map(|(projection, inclusion)| Ok(projection.compose(inclusion)?)),
            )?;
            *map = bridge.compose(map)?;
        }
    }
    let ranges = maps
        .iter()
        .map(|component| component.range().clone())
        .collect::<Vec<_>>();
    let range = direct_sum_of_qpa_modules(&ranges)?.ok_or(
        DerivedCategoryError::EmptyProjectiveDescription {
            operation: "StarOfMapBetweenDecompProjectives",
        },
    )?;
    let inclusions = direct_sum_inclusions(&range)?;
    sum_nonempty_maps(
        maps.iter()
            .zip(&inclusions)
            .map(|(component, inclusion)| Ok(component.compose(inclusion)?)),
    )
}

/// Translate `StarOfMapBetweenProjectives` and its length-based dispatch.
///
/// # Errors
///
/// Propagates the selected indecomposable/decomposable source method. An empty
/// source description retains QPA's unchecked first-entry failure.
pub fn star_of_map_between_projectives<F: Field>(
    map: &PathAlgebraModuleHomomorphism<F>,
    source_description: &[usize],
    range_description: &[usize],
) -> Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError> {
    if source_description.len() <= 1 {
        let &vertex =
            source_description
                .first()
                .ok_or(DerivedCategoryError::EmptyProjectiveDescription {
                    operation: "StarOfMapBetweenProjectives",
                })?;
        star_of_map_between_indec_projectives(map, vertex, range_description)
    } else {
        star_of_map_between_decomp_projectives(map, source_description, range_description)
    }
}

/// Translate `ProjectiveToInjectiveFiniteComplex`.
///
/// QPA's dimension-vector-as-vertex-list behavior is retained literally.
///
/// # Errors
///
/// Requires finite bounds and propagates projective descriptions, star maps,
/// transposition, direct sums, and checked finite-complex construction.
pub fn projective_to_injective_finite_complex<F>(
    projective_complex: &ModuleComplex<F>,
) -> Result<ModuleComplex<F>, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    let start = finite_bound(projective_complex.lower_bound()?)
        .ok_or(DerivedCategoryError::ProjectiveToInjectiveNeedsFiniteLowerBound)?;
    let stop = finite_bound(projective_complex.upper_bound()?)
        .ok_or(DerivedCategoryError::ProjectiveToInjectiveFiniteNeedsFiniteUpperBound)?;
    let injectives =
        indec_injective_modules(projective_complex.object(&start)?.right_acting_algebra())?;
    let descriptions = description_of_finite_proj_complex(projective_complex)?;
    if matches!(
        projective_complex.length()?,
        Some(ExtendedInteger::Finite(ref length)) if *length == 1
    ) {
        let description =
            descriptions
                .first()
                .ok_or(DerivedCategoryError::EmptyProjectiveDescription {
                    operation: "ProjectiveToInjectiveFiniteComplex",
                })?;
        let &vertex =
            description
                .first()
                .ok_or(DerivedCategoryError::EmptyProjectiveDescription {
                    operation: "ProjectiveToInjectiveFiniteComplex",
                })?;
        let injective = injectives.get(vertex.wrapping_sub(1)).ok_or(
            DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange {
                vertex,
                vertex_count: injectives.len(),
            },
        )?;
        let result = stalk_complex(projective_complex.category().clone(), injective, start)?;
        return Ok(cache_injective_complex(result));
    }

    let first_description =
        descriptions
            .first()
            .ok_or(DerivedCategoryError::EmptyProjectiveDescription {
                operation: "ProjectiveToInjectiveFiniteComplex",
            })?;
    let mut previous_injective = modules_from_projective_description(
        &injectives,
        first_description,
        "ProjectiveToInjectiveFiniteComplex",
    )?;
    let mut maps = Vec::new();
    for (position, degree) in integer_interval(&(&start + Integer::from(1)), &stop)
        .into_iter()
        .enumerate()
    {
        let range_description = descriptions.get(position).ok_or(
            DerivedCategoryError::StarComponentPositionOutOfRange {
                position: position + 1,
                component_count: descriptions.len(),
            },
        )?;
        let source_description = descriptions.get(position + 1).ok_or(
            DerivedCategoryError::StarComponentPositionOutOfRange {
                position: position + 2,
                component_count: descriptions.len(),
            },
        )?;
        let differential = projective_complex.differential(&degree)?;
        let star =
            star_of_map_between_projectives(&differential, source_description, range_description)?;
        let source_injective = modules_from_projective_description(
            &injectives,
            source_description,
            "ProjectiveToInjectiveFiniteComplex",
        )?;
        maps.push(right_module_hom_over_algebra(
            source_injective.clone(),
            previous_injective,
            matrices_of_dual_map(&star),
        )?);
        previous_injective = source_injective;
    }
    let result = finite_complex(
        projective_complex.category().clone(),
        &start + Integer::from(1),
        maps,
    )?;
    Ok(cache_injective_complex(result))
}

fn derived_lazy_error(error: &DerivedCategoryError) -> ComplexError {
    ComplexError::DerivedCategoryLazyDifferential(error.to_string())
}

/// Translate `ProjectiveToInjectiveComplex`.
///
/// Finite input dispatches to
/// [`projective_to_injective_finite_complex`]. Otherwise the known prefix and
/// QPA's stored degree-function tail are constructed literally.
///
/// # Errors
///
/// Requires finite lower/highest-known degrees and propagates all finite
/// prefix construction errors. Failures during later lazy generation are
/// retained as [`ComplexError::DerivedCategoryLazyDifferential`].
pub fn projective_to_injective_complex<F>(
    projective_complex: &ModuleComplex<F>,
) -> Result<ModuleComplex<F>, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    if finite_bound(projective_complex.upper_bound()?).is_some() {
        return projective_to_injective_finite_complex(projective_complex);
    }
    let start = finite_bound(projective_complex.lower_bound()?)
        .ok_or(DerivedCategoryError::ProjectiveToInjectiveNeedsFiniteLowerBound)?;
    let mut stop = finite_bound(projective_complex.highest_known_degree()?)
        .ok_or(DerivedCategoryError::ProjectiveToInjectiveNeedsFiniteHighestKnownDegree)?;
    if start == stop {
        stop += Integer::from(1);
    }
    let algebra = projective_complex
        .object(&start)?
        .right_acting_algebra()
        .clone();
    let injectives = indec_injective_modules(&algebra)?;
    let mut range_description = description_of_proj_complex_in_degree(projective_complex, &start)?;
    let mut range_injective = modules_from_projective_description(
        &injectives,
        &range_description,
        "ProjectiveToInjectiveComplex",
    )?;
    let mut middle = Vec::new();
    for degree in integer_interval(&(&start + Integer::from(1)), &stop) {
        let source_description =
            description_of_proj_complex_in_degree(projective_complex, &degree)?;
        if source_description.is_empty() {
            break;
        }
        let differential = projective_complex.differential(&degree)?;
        let star = star_of_map_between_projectives(
            &differential,
            &source_description,
            &range_description,
        )?;
        let source_injective = modules_from_projective_description(
            &injectives,
            &source_description,
            "ProjectiveToInjectiveComplex",
        )?;
        middle.push(right_module_hom_over_algebra(
            source_injective.clone(),
            range_injective,
            matrices_of_dual_map(&star),
        )?);
        range_injective = source_injective;
        range_description = source_description;
    }

    let input = projective_complex.clone();
    let lazy_injectives = Arc::new(injectives);
    let lazy_algebra = algebra;
    let positive = ComplexHalf::Position {
        function: Arc::new(move |output, degree| {
            (|| -> Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError> {
                let range_object = output.object(&(degree - Integer::from(1)))?;
                let source_description = description_of_proj_complex_in_degree(&input, degree)?;
                let range_description =
                    description_of_proj_complex_in_degree(&input, &(degree - Integer::from(1)))?;
                let differential = input.differential(degree)?;
                if differential.is_zero() {
                    return Ok(zero_module_mapping(
                        zero_module(&lazy_algebra)?,
                        range_object,
                    )?);
                }
                let star = star_of_map_between_projectives(
                    &differential,
                    &source_description,
                    &range_description,
                )?;
                let source = modules_from_projective_description(
                    &lazy_injectives,
                    &source_description,
                    "ProjectiveToInjectiveComplex",
                )?;
                Ok(right_module_hom_over_algebra(
                    source,
                    range_object,
                    matrices_of_dual_map(&star),
                )?)
            })()
            .map_err(|error| derived_lazy_error(&error))
        }),
        store_values: true,
    };
    let result = complex(
        projective_complex.category().clone(),
        &start + Integer::from(1),
        middle,
        positive,
        ComplexHalf::Zero,
    )?;
    Ok(cache_injective_complex(result))
}

/// Translate `TauOfComplex`.
///
/// # Errors
///
/// Propagates projective resolution, signed shift, projective-to-injective
/// conversion, and both `CutComplexAbove` scans.
pub fn tau_of_complex<F>(
    complex: &ModuleComplex<F>,
) -> Result<ModuleComplex<F>, DerivedCategoryError>
where
    F: QpaCoefficientField + Send + Sync + 'static,
    F::Element: Send + Sync,
{
    let projective = projective_resolution_of_complex(complex)?.shift(&Integer::from(1))?;
    let injective = projective_to_injective_complex(&projective)?.cut_complex_above()?;
    let tau = projective_resolution_of_complex(&injective)?.cut_complex_above()?;
    Ok(cache_projective_complex(tau))
}

fn sum_component_composites<F: Field>(
    left: &[PathAlgebraModuleHomomorphism<F>],
    right: &[PathAlgebraModuleHomomorphism<F>],
) -> Result<PathAlgebraModuleHomomorphism<F>, DerivedCategoryError> {
    if left.len() != right.len() {
        return Err(DerivedCategoryError::MapComponentCountMismatch {
            left: left.len(),
            right: right.len(),
        });
    }
    let mut terms = left
        .iter()
        .zip(right)
        .map(|(first, second)| first.compose(second));
    let Some(first) = terms.next() else {
        return Err(DerivedCategoryError::EmptyMapComponentList);
    };
    let mut sum = first?;
    for term in terms {
        sum = sum.add(&term?)?;
    }
    Ok(sum)
}

/// Translate `DirectSumMinusSummands(M,list)`.
///
/// Positions in `removed` are one-based and are inserted sequentially exactly
/// as in the GAP source.  The result is the canonical inclusion of the
/// complement into `M` and projection from `M` onto the complement.
///
/// # Errors
///
/// Requires direct-sum metadata, validates source list positions, and
/// propagates direct-sum, zero-map, composition, and addition failures.
pub fn direct_sum_minus_summands<F: Field>(
    module: &PathAlgebraModule<F>,
    removed: &[usize],
) -> Result<
    (
        PathAlgebraModuleHomomorphism<F>,
        PathAlgebraModuleHomomorphism<F>,
    ),
    DerivedCategoryError,
> {
    let summands = module
        .direct_sum_summands()
        .ok_or(DerivedCategoryError::NotDirectSumModule)?;
    for &position in removed {
        if position == 0 || position > summands.len() {
            return Err(DerivedCategoryError::SummandPositionOutOfRange {
                position,
                summand_count: summands.len(),
            });
        }
    }
    let redundant = removed
        .iter()
        .map(|&position| summands[position - 1].clone())
        .collect::<Vec<_>>();
    let complement = summands
        .iter()
        .enumerate()
        .filter(|(position, _)| !removed.contains(&(position + 1)))
        .map(|(_, summand)| summand.clone())
        .collect::<Vec<_>>();
    if complement.is_empty() {
        let zero = zero_module(module.right_acting_algebra())?;
        return Ok((
            zero_module_mapping(zero.clone(), module.clone())?,
            zero_module_mapping(module.clone(), zero)?,
        ));
    }
    let complement_module = direct_sum_of_qpa_modules(&complement)?
        .ok_or(DerivedCategoryError::EmptyMapComponentList)?;

    let mut source_components = direct_sum_projections(&complement_module)?;
    for (removed_index, (&position, summand)) in removed.iter().zip(&redundant).enumerate() {
        if position > source_components.len() + 1 {
            return Err(DerivedCategoryError::ComponentInsertionOutOfRange {
                position,
                length: source_components.len(),
            });
        }
        source_components.insert(
            position - 1,
            zero_module_mapping(complement_module.clone(), summand.clone())?,
        );
        debug_assert_eq!(
            removed_index + 1,
            source_components.len() - complement.len()
        );
    }
    let inclusion = sum_component_composites(&source_components, &direct_sum_inclusions(module)?)?;

    let mut range_components = direct_sum_inclusions(&complement_module)?;
    for (&position, summand) in removed.iter().zip(&redundant) {
        if position > range_components.len() + 1 {
            return Err(DerivedCategoryError::ComponentInsertionOutOfRange {
                position,
                length: range_components.len(),
            });
        }
        range_components.insert(
            position - 1,
            zero_module_mapping(summand.clone(), complement_module.clone())?,
        );
    }
    let projection = sum_component_composites(&direct_sum_projections(module)?, &range_components)?;
    Ok((inclusion, projection))
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field,
        finite_field::{FiniteFieldCore, Z2},
        rational::{Q, Rational},
    };
    use qpa_core::{
        ArrowSpec, DenseMatrix, PathAlgebra, PathAlgebraModuleHomomorphism, Quiver, QuiverAlgebra,
        VertexIndex, direct_sum_of_qpa_modules, identity_module_mapping, zero_module_mapping,
    };
    use serde_json::{Value, json};

    use super::*;
    use crate::{
        ComplexHalf, cat_of_right_algebra_modules, complex, indec_injective_modules,
        indec_projective_modules, quotient_by_relations, simple_modules, stalk_complex,
        zero_complex,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn matrix_json<F: Field>(
        matrix: &DenseMatrix<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!(
            (0..matrix.nrows())
                .map(|row| matrix.row(row).iter().map(coefficient).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        )
    }

    fn matrices_json<F: Field>(
        matrices: &[DenseMatrix<F>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!(
            matrices
                .iter()
                .map(|matrix| matrix_json(matrix, coefficient))
                .collect::<Vec<_>>()
        )
    }

    fn map_json<F: Field>(
        map: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!([
            map.source().dimension_vector(),
            map.range().dimension_vector(),
            matrices_json(map.vertex_matrices(), coefficient),
        ])
    }

    fn maybe_description_json(description: Vec<Option<usize>>) -> Value {
        Value::Array(
            description
                .into_iter()
                .map(|entry| entry.map_or_else(|| json!("false"), |entry| json!(entry)))
                .collect(),
        )
    }

    fn isomorphism_json<F: Field>(
        value: Option<ProjectiveIsomorphismComponent<F>>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        value.map_or_else(
            || json!([]),
            |value| {
                json!([
                    map_json(value.map(), coefficient),
                    value.vertex(),
                    value.source_position(),
                    value.range_position(),
                ])
            },
        )
    }

    #[allow(clippy::needless_pass_by_value)]
    fn derivedcat_reduction_case<F>(
        field: F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value
    where
        F: QpaCoefficientField + Send + Sync + 'static,
        F::Element: Send + Sync,
    {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(field, Quiver::from_indexed(1, Vec::new()).unwrap()).unwrap(),
        );
        let projective = indec_projective_modules(&algebra).unwrap()[0].clone();
        let source = direct_sum_of_qpa_modules(&[projective.clone(), projective.clone()])
            .unwrap()
            .unwrap();
        let range = direct_sum_of_qpa_modules(&[projective.clone(), projective.clone()])
            .unwrap()
            .unwrap();
        let one = identity_module_mapping(projective.clone()).unwrap();
        let zero = zero_module_mapping(projective.clone(), projective).unwrap();
        let differential = map_from_source_components(
            &source,
            &range,
            &[vec![one.clone(), one.clone()], vec![one.clone(), zero]],
        )
        .unwrap();
        let no_iso =
            check_for_isomorphism(&zero_module_mapping(source.clone(), range.clone()).unwrap())
                .unwrap();
        let cancellation = check_for_isomorphism(&differential).unwrap().unwrap();
        let source_reduced =
            direct_sum_minus_summands(&source, &[cancellation.source_position()]).unwrap();
        let range_reduced =
            direct_sum_minus_summands(&range, &[cancellation.range_position()]).unwrap();
        let reduced = reduce_differential(
            &differential,
            &cancellation,
            source_reduced.0.source(),
            range_reduced.0.source(),
        )
        .unwrap();
        let reduced_cancellation = check_for_isomorphism(&reduced).unwrap().unwrap();
        let source_zero =
            direct_sum_minus_summands(reduced.source(), &[reduced_cancellation.source_position()])
                .unwrap();
        let range_zero =
            direct_sum_minus_summands(reduced.range(), &[reduced_cancellation.range_position()])
                .unwrap();
        let reduced_zero = reduce_differential(
            &reduced,
            &reduced_cancellation,
            source_zero.0.source(),
            range_zero.0.source(),
        )
        .unwrap();
        let t_map = direct_sum_projections(&source).unwrap()[0].clone();
        let reduced_t =
            reduce_t_map(&t_map, &differential, &cancellation, &source_reduced.0).unwrap();
        let range_projections = direct_sum_projections(&range).unwrap();
        let previous = range_projections[0].add(&range_projections[1]).unwrap();
        let reduced_previous = reduce_previous_differential(&previous, &range_reduced.0).unwrap();
        let source_projections = direct_sum_projections(&source).unwrap();
        let previous_t = source_projections[0].add(&source_projections[1]).unwrap();
        let reduced_previous_t = reduce_previous_t_map(&previous_t, &source_reduced.0).unwrap();

        json!({
            "differential": map_json(&differential, coefficient),
            "no_iso": isomorphism_json(no_iso, coefficient),
            "iso": isomorphism_json(Some(cancellation), coefficient),
            "reduced": map_json(&reduced, coefficient),
            "reduced_iso": isomorphism_json(Some(reduced_cancellation), coefficient),
            "reduced_zero": map_json(&reduced_zero, coefficient),
            "reduced_t": map_json(&reduced_t, coefficient),
            "reduced_previous": map_json(&reduced_previous, coefficient),
            "reduced_previous_t": map_json(&reduced_previous_t, coefficient),
        })
    }

    fn extended_string(value: &ExtendedInteger) -> String {
        match value {
            ExtendedInteger::NegativeInfinity => "-infinity".to_owned(),
            ExtendedInteger::Finite(value) => value.to_string(),
            ExtendedInteger::PositiveInfinity => "infinity".to_owned(),
        }
    }

    fn optional_extended_string(value: Option<ExtendedInteger>) -> String {
        value.map_or_else(|| "fail".to_owned(), |value| extended_string(&value))
    }

    fn derivedcat_complex_json<F>(
        complex: &ModuleComplex<F>,
        degrees: &[i64],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value
    where
        F: QpaCoefficientField + Send + Sync + 'static,
        F::Element: Send + Sync,
    {
        json!({
            "upper": optional_extended_string(complex.upper_bound().unwrap()),
            "lower": optional_extended_string(complex.lower_bound().unwrap()),
            "length": optional_extended_string(complex.length().unwrap()),
            "projective": is_projective_complex(complex).unwrap(),
            "objects": degrees.iter().map(|degree| complex
                .object(&Integer::from(*degree)).unwrap().dimension_vector().to_vec())
                .collect::<Vec<_>>(),
            "differentials": degrees.iter().map(|degree| map_json(
                &complex.differential(&Integer::from(*degree)).unwrap(), coefficient))
                .collect::<Vec<_>>(),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn derivedcat_projective_resolution_case<F>(
        field: F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value
    where
        F: QpaCoefficientField + Send + Sync + 'static,
        F::Element: Send + Sync,
    {
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let category = cat_of_right_algebra_modules(&a2).unwrap();
        let simple = simple_modules(&a2).unwrap()[0].clone();
        let projective = indec_projective_modules(&a2).unwrap()[0].clone();
        let projective_complex =
            stalk_complex(category.clone(), &projective, Integer::from(0)).unwrap();
        let already = projective_resolution_of_complex(&projective_complex).unwrap();
        let stalk = stalk_complex(category.clone(), &simple, Integer::from(2)).unwrap();
        let stalk_resolution = projective_resolution_of_complex(&stalk).unwrap();
        let zero_complex = finite_complex(
            category.clone(),
            Integer::from(1),
            vec![zero_module_mapping(simple.clone(), simple.clone()).unwrap()],
        )
        .unwrap();
        let zero_resolution = projective_resolution_of_complex(&zero_complex).unwrap();

        let loop_parent = PathAlgebra::new(
            field,
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_category = cat_of_right_algebra_modules(&dual).unwrap();
        let dual_simple = simple_modules(&dual).unwrap()[0].clone();
        let dual_zero_complex = finite_complex(
            dual_category,
            Integer::from(1),
            vec![zero_module_mapping(dual_simple.clone(), dual_simple).unwrap()],
        )
        .unwrap();
        let dual_zero_resolution = projective_resolution_of_complex(&dual_zero_complex).unwrap();
        let degrees = [-2, -1, 0, 1, 2, 3, 4, 5, 6];

        json!({
            "already_identical": already == projective_complex,
            "already": derivedcat_complex_json(&already, &degrees, coefficient),
            "stalk": derivedcat_complex_json(&stalk_resolution, &degrees, coefficient),
            "zero_map": derivedcat_complex_json(&zero_resolution, &degrees, coefficient),
            "dual_zero_map": derivedcat_complex_json(
                &dual_zero_resolution, &degrees, coefficient),
        })
    }

    fn derivedcat_property_complex_json<F>(
        complex: &ModuleComplex<F>,
        degrees: &[i64],
        property: bool,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value
    where
        F: QpaCoefficientField + Send + Sync + 'static,
        F::Element: Send + Sync,
    {
        json!({
            "upper": optional_extended_string(complex.upper_bound().unwrap()),
            "lower": optional_extended_string(complex.lower_bound().unwrap()),
            "length": optional_extended_string(complex.length().unwrap()),
            "property": property,
            "objects": degrees.iter().map(|degree| complex
                .object(&Integer::from(*degree)).unwrap().dimension_vector().to_vec())
                .collect::<Vec<_>>(),
            "differentials": degrees.iter().map(|degree| map_json(
                &complex.differential(&Integer::from(*degree)).unwrap(), coefficient))
                .collect::<Vec<_>>(),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn derivedcat_tau_case<F>(field: F, coefficient: impl Fn(&F::Element) -> Value + Copy) -> Value
    where
        F: QpaCoefficientField + Send + Sync + 'static,
        F::Element: Send + Sync,
    {
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let category = cat_of_right_algebra_modules(&a2).unwrap();
        let projectives = indec_projective_modules(&a2).unwrap();
        let sum = direct_sum_of_qpa_modules(&projectives).unwrap().unwrap();
        let identity = identity_module_mapping(sum).unwrap();
        let star_decomposable =
            star_of_map_between_projectives(&identity, &[1, 2], &[1, 2]).unwrap();
        let finite = finite_complex(category.clone(), Integer::from(1), vec![identity]).unwrap();
        let finite_injective = projective_to_injective_complex(&finite).unwrap();
        let stalk = stalk_complex(category, &projectives[0], Integer::from(2)).unwrap();
        let stalk_injective = projective_to_injective_finite_complex(&stalk).unwrap();

        let loop_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_category = cat_of_right_algebra_modules(&dual).unwrap();
        let dual_projective = indec_projective_modules(&dual).unwrap()[0].clone();
        let generator =
            minimal_generating_set_of_module(dual_projective.clone()).unwrap()[0].clone();
        let quotient_x =
            include_in_path_algebra(&dual.quiver().arrows().next().unwrap().path(), &dual).unwrap();
        let x_image = generator.act(&quotient_x).unwrap();
        let x_map = hom_from_projective(&x_image, &dual_projective)
            .unwrap()
            .unwrap();
        let star_indecomposable = star_of_map_between_projectives(&x_map, &[1], &[1]).unwrap();
        let lazy_x = x_map.clone();
        let infinite = complex(
            dual_category,
            Integer::from(0),
            vec![x_map],
            ComplexHalf::Position {
                function: Arc::new(move |_, _| Ok(lazy_x.clone())),
                store_values: true,
            },
            ComplexHalf::Zero,
        )
        .unwrap();
        let infinite_injective = projective_to_injective_complex(&infinite).unwrap();

        let semisimple = QuiverAlgebra::Path(
            PathAlgebra::new(field, Quiver::from_indexed(1, Vec::new()).unwrap()).unwrap(),
        );
        let semisimple_category = cat_of_right_algebra_modules(&semisimple).unwrap();
        let semisimple_projective = indec_projective_modules(&semisimple).unwrap()[0].clone();
        let tau = tau_of_complex(
            &stalk_complex(
                semisimple_category,
                &semisimple_projective,
                Integer::from(0),
            )
            .unwrap(),
        )
        .unwrap();
        let degrees = [-2, -1, 0, 1, 2, 3, 4, 5];

        json!({
            "star_decomposable": map_json(&star_decomposable, coefficient),
            "star_indecomposable": map_json(&star_indecomposable, coefficient),
            "finite": derivedcat_property_complex_json(
                &finite_injective,
                &degrees,
                is_injective_complex(&finite_injective).unwrap(),
                coefficient,
            ),
            "stalk": derivedcat_property_complex_json(
                &stalk_injective,
                &degrees,
                is_injective_complex(&stalk_injective).unwrap(),
                coefficient,
            ),
            "infinite": derivedcat_property_complex_json(
                &infinite_injective,
                &degrees,
                is_injective_complex(&infinite_injective).unwrap(),
                coefficient,
            ),
            "tau": derivedcat_property_complex_json(
                &tau,
                &degrees,
                is_projective_complex(&tau).unwrap(),
                coefficient,
            ),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn derivedcat_case<F>(field: F, coefficient: impl Fn(&F::Element) -> Value + Copy) -> Value
    where
        F: QpaCoefficientField + Send + Sync + 'static,
        F::Element: Send + Sync,
    {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                field,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let category = cat_of_right_algebra_modules(&algebra).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let injectives = indec_injective_modules(&algebra).unwrap();
        let simples = simple_modules(&algebra).unwrap();
        let nonprojective = simples
            .iter()
            .find(|module| !is_projective_module(module).unwrap())
            .unwrap()
            .clone();
        let noninjective = simples
            .iter()
            .find(|module| !is_injective_module(module).unwrap())
            .unwrap()
            .clone();
        let projective_sum = direct_sum_of_qpa_modules(&[
            projectives[0].clone(),
            projectives[1].clone(),
            projectives[0].clone(),
        ])
        .unwrap()
        .unwrap();
        let injective_sum = direct_sum_of_qpa_modules(&[
            injectives[1].clone(),
            injectives[0].clone(),
            injectives[1].clone(),
        ])
        .unwrap()
        .unwrap();
        let projective_complex =
            stalk_complex(category.clone(), &projective_sum, Integer::from(1)).unwrap();
        let injective_complex =
            stalk_complex(category.clone(), &injective_sum, Integer::from(-1)).unwrap();
        let nonprojective_complex =
            stalk_complex(category.clone(), &nonprojective, Integer::from(0)).unwrap();
        let noninjective_complex =
            stalk_complex(category, &noninjective, Integer::from(0)).unwrap();
        let projective_cached_before = projective_complex.cached_is_projective_complex().is_some();
        let injective_cached_before = injective_complex.cached_is_injective_complex().is_some();
        let identity = identity_module_mapping(projective_sum.clone()).unwrap();
        let components = find_all_map_components(&identity).unwrap();
        let complement = direct_sum_minus_summands(&projective_sum, &[2]).unwrap();
        let all_removed = direct_sum_minus_summands(&projective_sum, &[1, 2, 3]).unwrap();

        json!({
            "projective_complex": is_projective_complex(&projective_complex).unwrap(),
            "injective_complex": is_injective_complex(&injective_complex).unwrap(),
            "nonprojective_complex": is_projective_complex(&nonprojective_complex).unwrap(),
            "noninjective_complex": is_injective_complex(&noninjective_complex).unwrap(),
            "cache": [
                projective_cached_before,
                projective_complex.cached_is_projective_complex().is_some(),
                injective_cached_before,
                injective_complex.cached_is_injective_complex().is_some(),
            ],
            "compare_projective": [
                compare_with_indec_projective(&projectives[0]).unwrap().map_or_else(
                    || json!("false"), |value| json!(value)),
                compare_with_indec_projective(&projective_sum).unwrap().map_or_else(
                    || json!("false"), |value| json!(value)),
                compare_with_indec_projective(&nonprojective).unwrap().map_or_else(
                    || json!("false"), |value| json!(value)),
            ],
            "compare_injective": [
                compare_with_indec_injective(&injectives[1]).unwrap().map_or_else(
                    || json!("false"), |value| json!(value)),
                compare_with_indec_injective(&injective_sum).unwrap().map_or_else(
                    || json!("false"), |value| json!(value)),
                compare_with_indec_injective(&noninjective).unwrap().map_or_else(
                    || json!("false"), |value| json!(value)),
            ],
            "projective_descriptions": [
                maybe_description_json(description_of_projective_module(&projectives[0]).unwrap()),
                maybe_description_json(description_of_projective_module(&projective_sum).unwrap()),
                maybe_description_json(description_of_projective_module(
                    &zero_module(&algebra).unwrap()).unwrap()),
            ],
            "degree_descriptions": [
                description_of_proj_complex_in_degree(
                    &projective_complex, &Integer::from(0)).unwrap(),
                description_of_proj_complex_in_degree(
                    &projective_complex, &Integer::from(1)).unwrap(),
                description_of_proj_complex_in_degree(
                    &projective_complex, &Integer::from(2)).unwrap(),
                description_of_inj_complex_in_degree(
                    &injective_complex, &Integer::from(-1)).unwrap(),
            ],
            "finite_descriptions": [
                description_of_finite_proj_complex(&projective_complex).unwrap(),
                description_of_finite_inj_complex(&injective_complex).unwrap(),
            ],
            "components": components.iter().map(|row| row.iter()
                .map(|map| map_json(map, coefficient)).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
            "dual_matrices": matrices_json(
                &matrices_of_dual_map(&complement.0), coefficient),
            "complement": [
                map_json(&complement.0, coefficient),
                map_json(&complement.1, coefficient),
            ],
            "all_removed": [
                map_json(&all_removed.0, coefficient),
                map_json(&all_removed.1, coefficient),
            ],
        })
    }

    #[test]
    fn pinned_qpa_derivedcat_foundation_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/derivedcat_foundation.schema-v1.json"
        ))
        .unwrap();
        let actual = json!({
            "schema_version": 1,
            "kind": "derivedcat_foundation",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": derivedcat_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "rational": derivedcat_case(Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn pinned_qpa_derivedcat_reduction_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/derivedcat_reduction.schema-v1.json"
        ))
        .unwrap();
        let actual = json!({
            "schema_version": 1,
            "kind": "derivedcat_reduction",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": derivedcat_reduction_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "rational": derivedcat_reduction_case(Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn pinned_qpa_derivedcat_projective_resolution_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/derivedcat_projective_resolution.schema-v1.json"
        ))
        .unwrap();
        let actual = json!({
            "schema_version": 1,
            "kind": "derivedcat_projective_resolution",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": derivedcat_projective_resolution_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "rational": derivedcat_projective_resolution_case(Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn pinned_qpa_derivedcat_tau_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/derivedcat_tau.schema-v1.json"
        ))
        .unwrap();
        let actual = json!({
            "schema_version": 1,
            "kind": "derivedcat_tau",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": derivedcat_tau_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "rational": derivedcat_tau_case(Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn derivedcat_tau_source_failures_are_typed() {
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                Q,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let category = cat_of_right_algebra_modules(&a2).unwrap();
        let projective = indec_projective_modules(&a2).unwrap()[0].clone();
        let identity = identity_module_mapping(projective.clone()).unwrap();
        let finite = finite_complex(category, Integer::from(1), vec![identity.clone()]).unwrap();
        assert!(matches!(
            projective_to_injective_finite_complex(&finite),
            Err(DerivedCategoryError::ProjectiveDescriptionVertexOutOfRange { vertex: 0, .. })
        ));
        assert!(matches!(
            star_of_map_between_projectives(&identity, &[], &[1]),
            Err(DerivedCategoryError::EmptyProjectiveDescription { .. })
        ));

        let loop_parent = PathAlgebra::new(
            Q,
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_category = cat_of_right_algebra_modules(&dual).unwrap();
        let dual_projective = indec_projective_modules(&dual).unwrap()[0].clone();
        let generator =
            minimal_generating_set_of_module(dual_projective.clone()).unwrap()[0].clone();
        let quotient_x =
            include_in_path_algebra(&dual.quiver().arrows().next().unwrap().path(), &dual).unwrap();
        let x_map = hom_from_projective(&generator.act(&quotient_x).unwrap(), &dual_projective)
            .unwrap()
            .unwrap();
        let repeating = complex(
            dual_category,
            Integer::from(0),
            vec![x_map.clone()],
            ComplexHalf::Repeat(vec![x_map]),
            ComplexHalf::Zero,
        )
        .unwrap();
        assert!(matches!(
            projective_to_injective_complex(&repeating),
            Err(DerivedCategoryError::ProjectiveToInjectiveNeedsFiniteHighestKnownDegree)
        ));

        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/derivedcat_tau_errors.schema-v1.json"
        ))
        .unwrap();
        assert!(
            fixture["success"]
                .as_object()
                .unwrap()
                .values()
                .all(|value| value == &json!(false))
        );
    }

    #[test]
    fn derivedcat_resolution_source_failures_are_typed() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Q,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let category = cat_of_right_algebra_modules(&algebra).unwrap();
        let simple = simple_modules(&algebra).unwrap()[0].clone();
        let identity_complex = finite_complex(
            category.clone(),
            Integer::from(1),
            vec![identity_module_mapping(simple.clone()).unwrap()],
        )
        .unwrap();
        assert!(projective_resolution_of_complex(&identity_complex).is_err());

        let zero = zero_module_mapping(simple.clone(), simple).unwrap();
        let periodic = complex(
            category.clone(),
            Integer::from(0),
            vec![zero.clone()],
            ComplexHalf::Repeat(vec![zero.clone()]),
            ComplexHalf::Repeat(vec![zero]),
        )
        .unwrap();
        assert!(matches!(
            projective_resolution_of_complex(&periodic),
            Err(DerivedCategoryError::ProjectiveResolutionNeedsFiniteLowerBound)
        ));
        assert!(matches!(
            projective_resolution_of_complex(&zero_complex(category).unwrap()),
            Err(DerivedCategoryError::ProjectiveResolutionNeedsFiniteLowerBound)
        ));

        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/derivedcat_resolution_errors.schema-v1.json"
        ))
        .unwrap();
        assert!(
            fixture["success"]
                .as_object()
                .unwrap()
                .values()
                .all(|value| value == &json!(false))
        );
    }

    #[test]
    fn qpa_derivedcat_specialized_printing_matches() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Q,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let category = cat_of_right_algebra_modules(&algebra).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let injectives = indec_injective_modules(&algebra).unwrap();
        let projective_sum = direct_sum_of_qpa_modules(&[
            projectives[0].clone(),
            projectives[1].clone(),
            projectives[0].clone(),
        ])
        .unwrap()
        .unwrap();
        let injective_sum = direct_sum_of_qpa_modules(&[
            injectives[1].clone(),
            injectives[0].clone(),
            injectives[1].clone(),
        ])
        .unwrap()
        .unwrap();
        let projective_complex =
            stalk_complex(category.clone(), &projective_sum, Integer::from(1)).unwrap();
        let injective_complex = stalk_complex(category, &injective_sum, Integer::from(-1)).unwrap();
        assert!(is_projective_complex(&projective_complex).unwrap());
        assert!(is_injective_complex(&injective_complex).unwrap());
        let actual = format!(
            "{}\n{}\n",
            qpa_projective_complex_string(&projective_complex).unwrap(),
            qpa_injective_complex_string(&injective_complex).unwrap(),
        );
        assert_eq!(
            actual,
            include_str!("../../../tests/fixtures/derivedcat_foundation_print.txt")
        );
    }

    #[test]
    fn derivedcat_foundation_source_failures_are_typed() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Q,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let category = cat_of_right_algebra_modules(&algebra).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let sum = direct_sum_of_qpa_modules(&[projectives[0].clone(), projectives[1].clone()])
            .unwrap()
            .unwrap();
        let zero = zero_module_mapping(projectives[0].clone(), projectives[0].clone()).unwrap();
        let periodic = complex(
            category,
            Integer::from(0),
            vec![zero.clone()],
            ComplexHalf::Repeat(vec![zero.clone()]),
            ComplexHalf::Repeat(vec![zero]),
        )
        .unwrap();
        assert!(matches!(
            description_of_finite_proj_complex(&periodic),
            Err(DerivedCategoryError::NonFiniteComplex)
        ));
        assert!(matches!(
            direct_sum_minus_summands(&projectives[0], &[1]),
            Err(DerivedCategoryError::NotDirectSumModule)
        ));
        assert!(matches!(
            direct_sum_minus_summands(&sum, &[3]),
            Err(DerivedCategoryError::SummandPositionOutOfRange { .. })
        ));
        assert!(
            find_all_map_components(&identity_module_mapping(projectives[0].clone()).unwrap())
                .is_err()
        );
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/derivedcat_foundation_errors.schema-v1.json"
        ))
        .unwrap();
        assert!(
            fixture["success"]
                .as_object()
                .unwrap()
                .values()
                .all(|value| value == &json!(false))
        );
    }
}
