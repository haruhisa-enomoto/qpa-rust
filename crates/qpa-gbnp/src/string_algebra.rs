//! QPA string-algebra word primitives, bounded strings, and bands.
//!
//! The algorithms follow QPA 1.37 `lib/stringalgebra.gi`, including its
//! one-character direct/inverse alphabet, heap-ordered string enumeration,
//! cyclic-band representative order, and KMP prefix-array convention.

use numerica::domains::Field;
use qpa_core::{Arrow, Quiver, Vertex};
use thiserror::Error;

use crate::{
    QpaStringAlgebra, QpaStringExtension, ValidityError, is_valid_string, qpa_string_direct_left,
    qpa_string_inverse_left, qpa_string_inverse_right, qpa_string_sigma_epsilon,
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
        band_representatives_less_than, bands_less_than, compute_lps_array,
        cyclic_permutations_of_band, is_a_band, maximum_length_of_direct_string, number_of_joints,
        source_of_string_arrow, strings_less_than, target_of_string_arrow,
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
            "representatives_3": band_representatives_less_than(&algebra, 3).unwrap()
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
