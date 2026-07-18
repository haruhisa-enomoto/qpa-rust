//! QPA string-algebra word primitives, bounded strings, and bands.
//!
//! The algorithms follow QPA 1.37 `lib/stringalgebra.gi`, including its
//! one-character direct/inverse alphabet, heap-ordered string enumeration,
//! cyclic-band representative order, and KMP prefix-array convention.

use numerica::domains::Field;
use qpa_core::{Arrow, ArrowSpec, Quiver, QuiverError, Vertex, VertexIndex};
use thiserror::Error;

use crate::{
    QpaStringAlgebra, QpaStringExtension, ValidityError, is_valid_string, qpa_string_direct_left,
    qpa_string_inverse_left, qpa_string_inverse_right, qpa_string_sigma_epsilon,
    validity::ensure_acceptable_quiver,
};

/// Failures made implicit by GAP indexing in `stringalgebra.gi`.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum StringAlgebraError {
    /// `SourceOfArrow` or `TargetOfArrow` rejects an alphabet position.
    #[error("The input string is invalid")]
    InvalidArrow,
    /// `IsABand` explicitly rejects an invalid string.
    #[error("The input string is invalid")]
    InvalidString,
    /// `MaximumLengthOfDirectString` indexes the first entry of an empty list.
    #[error("MaximumLengthOfDirectString requires at least one arrow")]
    EmptyArrowSet,
    /// QPA's dense enumeration loop indexed beyond its appended tree.
    #[error("the QPA string-enumeration tree ended before its level boundary")]
    SparseEnumerationTree,
    /// The binary heap boundary `2^(level+1)-1` is not representable.
    #[error("the requested string-enumeration level is too large")]
    EnumerationLevelOverflow,
    /// A hook/cohook helper returned QPA's accidental arrow-list value where a
    /// string is required.
    #[error("a string extension returned an arrow list instead of a string")]
    NonStringExtension,
    /// `BridgeQuiver` explicitly rejects a non-domestic string algebra.
    #[error("The String Algebra represented by the input is not a Domestic String Algebra.")]
    NonDomestic,
    /// Quiver construction failed at QPA's domesticity or bridge graph stage.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// String validity or its underlying algebra backend failed.
    #[error(transparent)]
    Validity(#[from] ValidityError),
}

/// Translate QPA's `computeLPSarray` zero-based prefix lengths.
#[must_use]
pub fn compute_lps_array(pattern: &str) -> Vec<usize> {
    let pattern = pattern.chars().collect::<Vec<_>>();
    let mut lps = vec![0; pattern.len()];
    let mut length = 0_usize;
    let mut index = 1_usize;
    while index < pattern.len() {
        if pattern[index] == pattern[length] {
            length += 1;
            lps[index] = length;
            index += 1;
        } else if length != 0 {
            length = lps[length - 1];
        } else {
            index += 1;
        }
    }
    lps
}

/// Translate QPA's orientation-aware `SourceOfArrow(Q,ch)`.
///
/// # Errors
///
/// Returns QPA's invalid-input diagnostic for a character outside the arrow
/// list or absent from its first-character name search.
pub fn source_of_string_arrow(
    quiver: &Quiver,
    character: char,
) -> Result<Vertex, StringAlgebraError> {
    let (arrow, inverse) = named_arrow(quiver, character)?;
    Ok(if inverse {
        arrow.target()
    } else {
        arrow.source()
    })
}

/// Translate QPA's orientation-aware `TargetOfArrow(Q,ch)`.
///
/// # Errors
///
/// Returns QPA's invalid-input diagnostic for a character outside the arrow
/// list or absent from its first-character name search.
pub fn target_of_string_arrow(
    quiver: &Quiver,
    character: char,
) -> Result<Vertex, StringAlgebraError> {
    let (arrow, inverse) = named_arrow(quiver, character)?;
    Ok(if inverse {
        arrow.source()
    } else {
        arrow.target()
    })
}

/// Translate QPA's left rotations from `CyclicPermutationOfABand`.
#[must_use]
pub fn cyclic_permutations_of_band(band: &str) -> Vec<String> {
    let mut characters = band.chars().collect::<Vec<_>>();
    let mut rotations = Vec::with_capacity(characters.len());
    for _ in 0..characters.len() {
        characters.rotate_left(1);
        rotations.push(characters.iter().collect());
    }
    rotations
}

/// Translate `MaximumLengthOfDirectString` by repeated direct-left extension.
///
/// Like QPA, this has no artificial termination bound when a direct extension
/// cycle exists.
///
/// # Errors
///
/// Propagates sigma/epsilon and validity failures and reports QPA's empty-arrow
/// indexing failure.
pub fn maximum_length_of_direct_string<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
) -> Result<usize, StringAlgebraError> {
    let quiver = algebra.qpa_underlying_path_algebra().quiver();
    let signs = qpa_string_sigma_epsilon(algebra)?;
    let mut lengths = Vec::with_capacity(quiver.number_of_arrows());
    for position in 1..=quiver.number_of_arrows() {
        let mut word = direct_character(position)?.to_string();
        loop {
            match qpa_string_direct_left(algebra, &word, &signs)? {
                QpaStringExtension::CannotPerform => {
                    lengths.push(word.chars().count());
                    break;
                }
                QpaStringExtension::Word(next) => word = next,
                QpaStringExtension::ArrowList(_) => {
                    return Err(StringAlgebraError::NonStringExtension);
                }
            }
        }
    }
    lengths
        .into_iter()
        .max()
        .ok_or(StringAlgebraError::EmptyArrowSet)
}

/// Translate `NumberOfJoints` as the number of one-arrow inverse-right moves.
///
/// # Errors
///
/// Propagates sigma/epsilon, alphabet, and validity failures.
pub fn number_of_joints<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
) -> Result<usize, StringAlgebraError> {
    let quiver = algebra.qpa_underlying_path_algebra().quiver();
    let signs = qpa_string_sigma_epsilon(algebra)?;
    let mut count = 0_usize;
    for position in 1..=quiver.number_of_arrows() {
        let word = direct_character(position)?.to_string();
        if !matches!(
            qpa_string_inverse_right(algebra, &word, &signs)?,
            QpaStringExtension::CannotPerform
        ) {
            count += 1;
        }
    }
    Ok(count)
}

/// Translate QPA's `IsABand` primitive/cyclic validity test.
///
/// # Errors
///
/// Propagates string-algebra validity and arrow-endpoint failures. An invalid
/// string produces QPA's explicit error rather than Boolean `false`.
pub fn is_a_band<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
) -> Result<bool, StringAlgebraError> {
    if !is_valid_string(algebra, input)? {
        return Err(StringAlgebraError::InvalidString);
    }
    if input.starts_with('(') || input.chars().count() == 1 {
        return Ok(false);
    }
    let quiver = algebra.qpa_underlying_path_algebra().quiver();
    let characters = input.chars().collect::<Vec<_>>();
    let Some(&last) = characters.last() else {
        return Err(StringAlgebraError::InvalidString);
    };
    let source = source_of_string_arrow(quiver, last)?;
    let target = target_of_string_arrow(quiver, characters[0])?;
    if source != target {
        return Ok(false);
    }
    let has_inverse = characters.iter().any(char::is_ascii_uppercase);
    let has_direct = characters.iter().any(char::is_ascii_lowercase);
    if !has_inverse || !has_direct {
        return Ok(false);
    }
    let square = format!("{input}{input}");
    if !is_valid_string(algebra, &square)? {
        return Ok(false);
    }
    let lps = compute_lps_array(input);
    let Some(&border) = lps.last() else {
        return Err(StringAlgebraError::InvalidString);
    };
    Ok(border == 0 || characters.len() % (characters.len() - border) != 0)
}

/// Translate QPA's heap-ordered `StringsLessThan(A,level)` enumeration.
///
/// The GAP method dispatches only on positive integers, but its installed body
/// contains a level-zero branch; Rust exposes that source branch at `level=0`.
///
/// # Errors
///
/// Propagates validity and extension failures and reports QPA's unchecked tree
/// and integer-boundary failures.
pub fn strings_less_than<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    level: usize,
) -> Result<Vec<String>, StringAlgebraError> {
    let quiver = algebra.qpa_underlying_path_algebra().quiver();
    if level == 0 {
        let mut strings = Vec::with_capacity(2 * quiver.number_of_vertices());
        for position in 1..=quiver.number_of_vertices() {
            strings.push(format!("({position},1)"));
            strings.push(format!("({position},-1)"));
        }
        return Ok(strings);
    }
    let boundary = 2_usize
        .checked_pow(
            u32::try_from(level + 1).map_err(|_| StringAlgebraError::EnumerationLevelOverflow)?,
        )
        .and_then(|value| value.checked_sub(1))
        .ok_or(StringAlgebraError::EnumerationLevelOverflow)?;
    let signs = qpa_string_sigma_epsilon(algebra)?;
    let mut strings = Vec::new();
    for position in 1..=quiver.number_of_vertices() {
        for sign in [1_i8, -1_i8] {
            let root = format!("({position},{sign})");
            let tree = make_string_tree(algebra, &root, boundary, &signs)?;
            strings.extend(tree.into_iter().map(|(_, word)| word));
        }
    }
    Ok(strings)
}

/// Translate `BandsLessThan` by filtering QPA's ordered string list.
///
/// # Errors
///
/// Propagates bounded enumeration and band-validity failures.
pub fn bands_less_than<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    level: usize,
) -> Result<Vec<String>, StringAlgebraError> {
    let strings = strings_less_than(algebra, level)?;
    let mut bands = Vec::new();
    for string in strings {
        if is_a_band(algebra, &string)? {
            bands.push(string);
        }
    }
    Ok(bands)
}

/// Translate `BandRepresentativesLessThan` in QPA's first-seen rotation order.
///
/// # Errors
///
/// Propagates bounded band enumeration failures.
pub fn band_representatives_less_than<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    level: usize,
) -> Result<Vec<String>, StringAlgebraError> {
    let bands = bands_less_than(algebra, level)?;
    let mut seen_rotations = Vec::<String>::new();
    let mut representatives = Vec::new();
    for band in bands {
        let Some(first) = band.chars().next() else {
            return Err(StringAlgebraError::InvalidString);
        };
        let Some(last) = band.chars().last() else {
            return Err(StringAlgebraError::InvalidString);
        };
        if !seen_rotations.contains(&band) && (first as u32) >= 97 && (last as u32) < 97 {
            representatives.push(band.clone());
            seen_rotations.extend(cyclic_permutations_of_band(&band));
        }
    }
    Ok(representatives)
}

/// Translate QPA's finite band-concatenation test
/// `IsDomesticStringAlgebra`.
///
/// QPA constructs this graph on every enumerated cyclic permutation of every
/// band, not only on band representatives; Rust preserves that distinction.
///
/// # Errors
///
/// Propagates string-algebra validation, bounded enumeration, and QPA's quiver
/// construction failures (including its zero-band `Quiver(0,[])` branch).
pub fn is_domestic_string_algebra<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
) -> Result<bool, StringAlgebraError> {
    ensure_string_algebra_input(algebra)?;
    let level = string_algebra_search_level(algebra)?;
    let bands = bands_less_than(algebra, level)?;
    Ok(domesticity_quiver(algebra, &bands)?.is_acyclic())
}

/// Translate QPA's complete `BridgeQuiver` construction.
///
/// The implementation follows the source's band-free substring search,
/// weak-bridge order, zero-string concatenation behavior, KMP removal stage,
/// and final arrow filtering.
///
/// # Errors
///
/// Propagates validation, enumeration, endpoint, and quiver-construction
/// failures. A cyclic domesticity graph returns QPA's explicit non-domestic
/// error.
pub fn bridge_quiver<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
) -> Result<Quiver, StringAlgebraError> {
    ensure_string_algebra_input(algebra)?;
    let level = string_algebra_search_level(algebra)?;
    let strings = strings_less_than(algebra, level)?;
    let mut bands = Vec::new();
    for string in &strings {
        if is_a_band(algebra, string)? {
            bands.push(string.clone());
        }
    }
    let representatives = band_representatives_from_bands(&bands)?;
    if !domesticity_quiver(algebra, &bands)?.is_acyclic() {
        return Err(StringAlgebraError::NonDomestic);
    }
    let (weak_bridges, indexed_weak_bridges) =
        weak_bridge_data(algebra, &representatives, &strings)?;
    reduce_weak_bridges(
        algebra,
        &representatives,
        &weak_bridges,
        &indexed_weak_bridges,
    )
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct WeakBridge {
    left_band: String,
    right_band: String,
    word: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct IndexedWeakBridge {
    source: usize,
    target: usize,
    word: String,
}

fn ensure_string_algebra_input<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
) -> Result<(), StringAlgebraError> {
    if !algebra
        .qpa_is_string_algebra()
        .map_err(ValidityError::from)?
    {
        return Err(ValidityError::NotStringAlgebra.into());
    }
    ensure_acceptable_quiver(algebra.qpa_underlying_path_algebra().quiver())?;
    Ok(())
}

fn string_algebra_search_level<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
) -> Result<usize, StringAlgebraError> {
    let joints = number_of_joints(algebra)?;
    let maximum = maximum_length_of_direct_string(algebra)?;
    maximum
        .checked_mul(joints)
        .and_then(|value| value.checked_mul(2))
        .ok_or(StringAlgebraError::EnumerationLevelOverflow)
}

fn domesticity_quiver<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    bands: &[String],
) -> Result<Quiver, StringAlgebraError> {
    let mut arrows = Vec::new();
    let mut count = 0_usize;
    for (left_index, left) in bands.iter().enumerate() {
        for (right_index, right) in bands.iter().enumerate() {
            if left_index != right_index {
                let candidate = format!("{left}{right}");
                if is_valid_string(algebra, &candidate)? && !candidate.starts_with('(') {
                    count += 1;
                    arrows.push(ArrowSpec::named(
                        vertex_index(right_index + 1)?,
                        vertex_index(left_index + 1)?,
                        format!("a{count}"),
                    ));
                }
            }
        }
    }
    Ok(Quiver::from_indexed(bands.len(), arrows)?)
}

fn band_representatives_from_bands(bands: &[String]) -> Result<Vec<String>, StringAlgebraError> {
    let mut rotations = Vec::<String>::new();
    let mut representatives = Vec::new();
    for band in bands {
        let Some(first) = band.chars().next() else {
            return Err(StringAlgebraError::InvalidString);
        };
        let Some(last) = band.chars().last() else {
            return Err(StringAlgebraError::InvalidString);
        };
        if !rotations.contains(band) && (first as u32) >= 97 && (last as u32) < 97 {
            representatives.push(band.clone());
            rotations.extend(cyclic_permutations_of_band(band));
        }
    }
    Ok(representatives)
}

fn is_band_free_string<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
) -> Result<bool, StringAlgebraError> {
    if input.starts_with('(') {
        return Ok(true);
    }
    let quiver = algebra.qpa_underlying_path_algebra().quiver();
    let characters = input.chars().collect::<Vec<_>>();
    for left in 0..characters.len() {
        for right in left + 1..characters.len() {
            if target_of_string_arrow(quiver, characters[left])?
                == source_of_string_arrow(quiver, characters[right])?
            {
                let candidate = characters[left..=right].iter().collect::<String>();
                if is_a_band(algebra, &candidate)? {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

fn band_free_strings<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    source: usize,
    target: usize,
    strings: &[String],
) -> Result<Vec<String>, StringAlgebraError> {
    let quiver = algebra.qpa_underlying_path_algebra().quiver();
    if source > quiver.number_of_vertices() || target > quiver.number_of_vertices() {
        return Ok(Vec::new());
    }
    let source_digit = source.to_string().chars().next();
    let mut result = Vec::new();
    for string in strings {
        let characters = string.chars().collect::<Vec<_>>();
        if characters.first() == Some(&'(') {
            if source == target && characters.get(1).copied() == source_digit {
                result.push(string.clone());
            }
        } else if let (Some(&first), Some(&last)) = (characters.first(), characters.last())
            && target_of_string_arrow(quiver, first)?.index().get() == target
            && source_of_string_arrow(quiver, last)?.index().get() == source
            && is_band_free_string(algebra, string)?
        {
            result.push(string.clone());
        }
    }
    Ok(result)
}

fn weak_bridge_data<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    representatives: &[String],
    strings: &[String],
) -> Result<(Vec<WeakBridge>, Vec<IndexedWeakBridge>), StringAlgebraError> {
    let quiver = algebra.qpa_underlying_path_algebra().quiver();
    let mut weak = Vec::new();
    let mut indexed = Vec::new();
    for (left_index, left) in representatives.iter().enumerate() {
        for (right_index, right) in representatives.iter().enumerate() {
            let right_first = right
                .chars()
                .next()
                .ok_or(StringAlgebraError::InvalidString)?;
            let left_last = left
                .chars()
                .last()
                .ok_or(StringAlgebraError::InvalidString)?;
            let source = target_of_string_arrow(quiver, right_first)?.index().get();
            let target = source_of_string_arrow(quiver, left_last)?.index().get();
            for word in band_free_strings(algebra, source, target, strings)? {
                let candidate = format!("{left}{word}{right}");
                if is_valid_string(algebra, &candidate)? && !candidate.starts_with('(') {
                    weak.push(WeakBridge {
                        left_band: right.clone(),
                        right_band: left.clone(),
                        word: word.clone(),
                    });
                    indexed.push(IndexedWeakBridge {
                        source: right_index + 1,
                        target: left_index + 1,
                        word,
                    });
                }
            }
        }
    }
    let specs = indexed
        .iter()
        .map(indexed_bridge_spec)
        .collect::<Result<Vec<_>, _>>()?;
    let _source_weak_quiver = Quiver::from_indexed(representatives.len(), specs)?;
    Ok((weak, indexed))
}

fn reduce_weak_bridges<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    representatives: &[String],
    weak: &[WeakBridge],
    indexed: &[IndexedWeakBridge],
) -> Result<Quiver, StringAlgebraError> {
    let mut retained = vec![true; weak.len()];
    for left in weak {
        for right in weak {
            if left.left_band == right.right_band {
                let joined = format!("{}{}", left.word, right.word);
                if is_valid_string(algebra, &joined)? && joined.starts_with('(') {
                    if is_band_free_string(algebra, &joined)? {
                        if let Some(position) = weak.iter().position(|bridge| bridge.word == joined)
                        {
                            retained[position] = false;
                        }
                    } else {
                        let mut remainder = String::new();
                        for rotation in cyclic_permutations_of_band(&left.left_band) {
                            if let Some(position) = kmp_search(&rotation, &joined) {
                                let joined_chars = joined.chars().collect::<Vec<_>>();
                                remainder = joined_chars[..position]
                                    .iter()
                                    .chain(
                                        joined_chars[position + rotation.chars().count()..].iter(),
                                    )
                                    .collect();
                                break;
                            }
                        }
                        if let Some(position) =
                            weak.iter().position(|bridge| bridge.word == remainder)
                        {
                            retained[position] = false;
                        }
                    }
                }
            }
        }
    }
    let specs = indexed
        .iter()
        .zip(retained)
        .filter(|(_, retain)| *retain)
        .map(|(bridge, _)| indexed_bridge_spec(bridge))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Quiver::from_indexed(representatives.len(), specs)?)
}

fn kmp_search(pattern: &str, text: &str) -> Option<usize> {
    let pattern = pattern.chars().collect::<Vec<_>>();
    let text = text.chars().collect::<Vec<_>>();
    if pattern.is_empty() {
        return Some(0);
    }
    text.windows(pattern.len())
        .position(|window| window == pattern)
}

fn indexed_bridge_spec(bridge: &IndexedWeakBridge) -> Result<ArrowSpec, StringAlgebraError> {
    Ok(ArrowSpec::named(
        vertex_index(bridge.source)?,
        vertex_index(bridge.target)?,
        bridge.word.clone(),
    ))
}

fn vertex_index(position: usize) -> Result<VertexIndex, StringAlgebraError> {
    VertexIndex::new(position).ok_or(StringAlgebraError::SparseEnumerationTree)
}

fn make_string_tree<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    root: &str,
    boundary: usize,
    signs: &crate::StringSigmaEpsilon,
) -> Result<Vec<(usize, String)>, StringAlgebraError> {
    let mut tree = vec![(1_usize, root.to_owned())];
    let direct = extension_word(qpa_string_direct_left(algebra, root, signs)?)?;
    let inverse = extension_word(qpa_string_inverse_left(algebra, root, signs)?)?;
    if let Some(word) = direct {
        tree.push((2, word));
    }
    if let Some(word) = inverse {
        tree.push((3, word));
    }
    if tree.len() == 1 {
        return Ok(Vec::new());
    }

    let mut cursor = 1_usize;
    loop {
        let (heap_position, word) = tree
            .get(cursor)
            .cloned()
            .ok_or(StringAlgebraError::SparseEnumerationTree)?;
        let left_position = heap_position
            .checked_mul(2)
            .ok_or(StringAlgebraError::EnumerationLevelOverflow)?;
        if left_position > boundary {
            break;
        }
        if let Some(next) = extension_word(qpa_string_direct_left(algebra, &word, signs)?)? {
            tree.push((left_position, next));
        }
        if let Some(next) = extension_word(qpa_string_inverse_left(algebra, &word, signs)?)? {
            tree.push((left_position + 1, next));
        }
        cursor += 1;
    }
    Ok(tree)
}

fn extension_word(extension: QpaStringExtension) -> Result<Option<String>, StringAlgebraError> {
    match extension {
        QpaStringExtension::CannotPerform => Ok(None),
        QpaStringExtension::Word(word) => Ok(Some(word)),
        QpaStringExtension::ArrowList(_) => Err(StringAlgebraError::NonStringExtension),
    }
}

fn named_arrow(quiver: &Quiver, character: char) -> Result<(Arrow, bool), StringAlgebraError> {
    let inverse = character.is_ascii_uppercase();
    let normalized = if inverse {
        char::from_u32(character as u32 + 32).ok_or(StringAlgebraError::InvalidArrow)?
    } else {
        character
    };
    let one_based = (normalized as u32)
        .checked_sub('a' as u32)
        .and_then(|value| value.checked_add(1))
        .ok_or(StringAlgebraError::InvalidArrow)?;
    if usize::try_from(one_based).map_err(|_| StringAlgebraError::InvalidArrow)?
        > quiver.number_of_arrows()
    {
        return Err(StringAlgebraError::InvalidArrow);
    }
    let arrow = quiver
        .arrows()
        .find(|arrow| arrow.name().starts_with(normalized))
        .ok_or(StringAlgebraError::InvalidArrow)?;
    Ok((arrow, inverse))
}

fn direct_character(position: usize) -> Result<char, StringAlgebraError> {
    let position = u32::try_from(position).map_err(|_| StringAlgebraError::InvalidArrow)?;
    char::from_u32(position + 96).ok_or(StringAlgebraError::InvalidArrow)
}

#[cfg(test)]
mod tests {
    use numerica::domains::{Field, finite_field::Z2, rational::Q};
    use qpa_core::{ArrowSpec, PathAlgebra, Quiver, VertexIndex};
    use serde_json::{Value, json};

    use crate::quotient_by_relations;

    use super::{
        band_representatives_less_than, bands_less_than, bridge_quiver, compute_lps_array,
        cyclic_permutations_of_band, is_a_band, is_domestic_string_algebra,
        maximum_length_of_direct_string, number_of_joints, source_of_string_arrow,
        strings_less_than, target_of_string_arrow,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn fixture_case<F: Field>(field: F) -> Value {
        let quiver = Quiver::from_indexed(
            4,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(1), vertex(2), "b"),
                ArrowSpec::named(vertex(2), vertex(3), "c"),
                ArrowSpec::named(vertex(3), vertex(4), "d"),
                ArrowSpec::named(vertex(3), vertex(4), "e"),
            ],
        )
        .unwrap();
        let parent = PathAlgebra::new(field, quiver.clone()).unwrap();
        let generators = parent.generators();
        let relations = [
            generators[5].multiply(&generators[6]).unwrap(),
            generators[6].multiply(&generators[7]).unwrap(),
        ];
        let algebra = quotient_by_relations(parent, &relations).unwrap();
        let endpoints = ['a', 'A', 'e', 'E'].map(|character| {
            [
                source_of_string_arrow(&quiver, character)
                    .unwrap()
                    .to_string(),
                target_of_string_arrow(&quiver, character)
                    .unwrap()
                    .to_string(),
            ]
        });
        let band_inputs = ["eca", "Ab", "aB", "De", "(1,1)", "A", "AbAb"];
        let bridge = bridge_quiver(&algebra).unwrap();
        let bridge_arrows = bridge
            .arrows()
            .map(|arrow| {
                [
                    arrow.source().to_string(),
                    arrow.target().to_string(),
                    arrow.to_string(),
                ]
            })
            .collect::<Vec<_>>();
        json!({
            "lps": [compute_lps_array("ababaca"), compute_lps_array("aaaa"), compute_lps_array("")],
            "endpoints": endpoints,
            "cyclic": cyclic_permutations_of_band("Abc"),
            "maximum_direct": maximum_length_of_direct_string(&algebra).unwrap(),
            "joints": number_of_joints(&algebra).unwrap(),
            "band_inputs": band_inputs,
            "is_band": band_inputs.map(|input| is_a_band(&algebra, input).unwrap()),
            "strings_1": strings_less_than(&algebra, 1).unwrap(),
            "strings_2": strings_less_than(&algebra, 2).unwrap(),
            "bands_3": bands_less_than(&algebra, 3).unwrap(),
            "representatives_3": band_representatives_less_than(&algebra, 3).unwrap(),
            "domestic": is_domestic_string_algebra(&algebra).unwrap(),
            "bridge_vertices": bridge.vertices().map(|vertex| vertex.to_string()).collect::<Vec<_>>(),
            "bridge_arrows": bridge_arrows
        })
    }

    #[test]
    fn pinned_qpa_stringalgebra_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/stringalgebra.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(fixture_case(Z2), expected["binary"]);
        assert_eq!(fixture_case(Q), expected["rational"]);
    }
}
