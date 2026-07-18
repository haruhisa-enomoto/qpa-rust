//! Projective/injective complex descriptions and helper maps translated from
//! QPA's `derivedcat.gd` and `derivedcat.gi`.

use std::fmt::Write as _;

use numerica::domains::{Field, integer::Integer};
use qpa_core::{
    DenseMatrix, ExtendedInteger, PathAlgebraModule, PathAlgebraModuleError,
    PathAlgebraModuleHomomorphism, QpaCoefficientField, direct_sum_inclusions,
    direct_sum_of_qpa_modules, direct_sum_projections, socle_of_module, top_of_module, zero_module,
    zero_module_mapping,
};
use thiserror::Error;

use crate::{
    complex::{ComplexError, ModuleComplex},
    module::{EnvelopingModuleError, is_injective_module, is_projective_module},
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
        indec_projective_modules, simple_modules, stalk_complex,
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
