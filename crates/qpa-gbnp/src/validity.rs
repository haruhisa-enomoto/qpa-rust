//! QPA string validity, hook/cohook moves, and local AR arrows.
//!
//! This is a direct translation of QPA 1.37 `lib/validity.gi`.  In
//! particular, it retains QPA's one-character arrow alphabet, its
//! `(vertex,sign)` notation for lazy strings, its concrete sigma/epsilon
//! choice, and the source order of hook/cohook searches.

use std::collections::HashSet;

use numerica::domains::Field;
use qpa_core::{
    Arrow, ArrowIndex, Path, PathAlgebra, PathAlgebraElement, PathError, Quiver, QuiverError,
    QuotientPathAlgebra, Vertex,
};
use thiserror::Error;

use crate::{
    PathalgGbnpError,
    pathalg::{path_algebra_is_string, quotient_is_string},
};

const CANNOT_PERFORM: &str = "Cannot Perform The Operation";

/// The algebra surface accepted by QPA's `[ IsAlgebra, IsString ]` methods in
/// `validity.gi`.
pub trait QpaStringAlgebra<F: Field> {
    /// Return QPA's underlying path algebra.
    fn qpa_underlying_path_algebra(&self) -> &PathAlgebra<F>;

    /// Return QPA's `RelationsOfAlgebra` list.
    fn qpa_defining_relations(&self) -> &[PathAlgebraElement<F>];

    /// Dispatch QPA's `IsStringAlgebra` property.
    ///
    /// # Errors
    ///
    /// Propagates the quotient's default GBNP predicate failures.
    fn qpa_is_string_algebra(&self) -> Result<bool, PathalgGbnpError>;
}

impl<F: Field> QpaStringAlgebra<F> for PathAlgebra<F> {
    fn qpa_underlying_path_algebra(&self) -> &PathAlgebra<F> {
        self
    }

    fn qpa_defining_relations(&self) -> &[PathAlgebraElement<F>] {
        &[]
    }

    fn qpa_is_string_algebra(&self) -> Result<bool, PathalgGbnpError> {
        Ok(path_algebra_is_string(self))
    }
}

impl<F: Field> QpaStringAlgebra<F> for QuotientPathAlgebra<F> {
    fn qpa_underlying_path_algebra(&self) -> &PathAlgebra<F> {
        self.original_path_algebra()
    }

    fn qpa_defining_relations(&self) -> &[PathAlgebraElement<F>] {
        self.relations()
    }

    fn qpa_is_string_algebra(&self) -> Result<bool, PathalgGbnpError> {
        quotient_is_string(self)
    }
}

/// QPA's two sign functions on the direct-arrow alphabet.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StringSigmaEpsilon {
    sigma: Vec<i8>,
    epsilon: Vec<i8>,
}

impl StringSigmaEpsilon {
    /// Return sigma in QPA arrow order.
    #[must_use]
    pub fn sigma(&self) -> &[i8] {
        &self.sigma
    }

    /// Return epsilon in QPA arrow order.
    #[must_use]
    pub fn epsilon(&self) -> &[i8] {
        &self.epsilon
    }
}

/// The heterogeneous string result of QPA's four internal extension
/// functions.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum QpaStringExtension {
    /// QPA's literal `"Cannot Perform The Operation"` result.
    CannotPerform,
    /// A positive-length string word.
    Word(String),
    /// The source bug in `QPAStringDirectRight`: one zero-string branch
    /// returns the one-element incoming-arrow list instead of its name.
    ArrowList(Vec<Arrow>),
}

impl QpaStringExtension {
    /// Return QPA's sentinel text or string word when this is a textual result.
    #[must_use]
    pub fn as_qpa_text(&self) -> Option<&str> {
        match self {
            Self::CannotPerform => Some(CANNOT_PERFORM),
            Self::Word(word) => Some(word),
            Self::ArrowList(_) => None,
        }
    }
}

/// Failures made implicit by GAP list indexing or method dispatch in
/// `validity.gi`.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ValidityError {
    /// QPA explicitly rejects a non-string algebra.
    #[error("The input algebra is not a string algebra")]
    NotStringAlgebra,
    /// QPA explicitly requires vertices `v1`, ... and arrows `a`, ... .
    #[error("The Quiver of A is not acceptable")]
    UnacceptableQuiver,
    /// A defining relation has no first magma element for QPA to inspect.
    #[error("a defining relation has no first monomial path")]
    EmptyRelation,
    /// `QPAStringQuiverRho` found no arrow for a relation character.
    #[error("relation character {character:?} does not name an arrow")]
    UnknownRelationCharacter {
        /// Character that failed QPA's nested arrow-name search.
        character: char,
    },
    /// A string helper received an empty positive-length word.
    #[error("the empty input has no first or last string letter")]
    EmptyInput,
    /// QPA's concrete sigma/epsilon construction indexed a missing incidence.
    #[error("QPAStringSigmaEps requires two outgoing arrows at an incidence-free vertex")]
    MissingSourceArrows,
    /// A string or arrow name is outside QPA's generated one-character alphabet.
    #[error("character {character:?} is outside the QPA arrow alphabet")]
    InvalidArrowCharacter {
        /// Rejected character.
        character: char,
    },
    /// A malformed lazy-string prefix triggered QPA's unchecked list access.
    #[error("the zero-length string notation is malformed")]
    MalformedZeroString,
    /// QPA explicitly rejects an invalid string in `LocalARQuiver`.
    #[error("The input string is invalid")]
    InvalidString,
    /// QPA's fixed-width lazy-string deletion attempted an impossible unbind.
    #[error("the lazy-string markers cannot be removed from this input")]
    InvalidLazyStringMarkers,
    /// A local-AR iteration encountered QPA's accidental arrow-list result.
    #[error("a string extension returned an arrow list instead of a string")]
    NonStringExtension,
    /// Path multiplication failed while reconstructing relation paths.
    #[error(transparent)]
    Path(#[from] PathError),
    /// A source-translated quiver incidence lookup failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// String-algebra property dispatch or its GBNP backend failed.
    #[error(transparent)]
    Pathalg(#[from] PathalgGbnpError),
}

/// Translate `QPAStringQuiverRho(Q,rho)` from relation words to quiver paths.
///
/// # Errors
///
/// Returns the checked counterpart of QPA's empty match/product and path
/// multiplication failures.
pub fn qpa_string_quiver_rho(
    quiver: &Quiver,
    relations: &[String],
) -> Result<Vec<Path>, ValidityError> {
    let arrows = quiver.arrows().collect::<Vec<_>>();
    relations
        .iter()
        .map(|relation| {
            let mut factors = Vec::new();
            for character in relation.chars() {
                let before = factors.len();
                for arrow in &arrows {
                    if first_character(arrow.name()) == Some(character) {
                        factors.push(arrow.clone());
                    }
                }
                if factors.len() == before {
                    return Err(ValidityError::UnknownRelationCharacter { character });
                }
            }
            if factors.is_empty() {
                return Err(ValidityError::EmptyRelation);
            }
            Path::from_arrow_product(quiver, &factors)?.ok_or(ValidityError::EmptyRelation)
        })
        .collect()
}

/// Translate QPA's concrete `QPAStringSigmaEps` choice.
///
/// # Errors
///
/// Returns checked errors for QPA's unchecked relation and incidence list
/// indexing.
pub fn qpa_string_sigma_epsilon<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
) -> Result<StringSigmaEpsilon, ValidityError> {
    let parent = algebra.qpa_underlying_path_algebra();
    let quiver = parent.quiver();
    let relation_paths = relation_paths(algebra)?;
    let mut sigma = vec![0_i8; quiver.number_of_arrows()];
    let mut epsilon = vec![0_i8; quiver.number_of_arrows()];

    let composable = |outgoing: &Arrow, incoming: &Arrow| -> Result<bool, ValidityError> {
        if outgoing.source() != incoming.target() {
            return Ok(false);
        }
        let product = incoming.path().multiply(&outgoing.path())?;
        Ok(!relation_paths.iter().any(|relation| relation == &product))
    };

    for vertex in quiver.vertices() {
        let incoming = quiver.incoming_arrows(&vertex)?;
        let outgoing = quiver.outgoing_arrows(&vertex)?;
        match incoming.as_slice() {
            [first, second] => {
                epsilon[arrow_position(first, quiver)?] = 1;
                epsilon[arrow_position(second, quiver)?] = -1;
                match outgoing.as_slice() {
                    [only] => {
                        if composable(only, first)? {
                            sigma[arrow_position(only, quiver)?] = -1;
                        } else if composable(only, second)? {
                            sigma[arrow_position(only, quiver)?] = 1;
                        }
                    }
                    [left, right] => {
                        if composable(left, first)? || composable(right, second)? {
                            sigma[arrow_position(left, quiver)?] = -1;
                            sigma[arrow_position(right, quiver)?] = 1;
                        } else {
                            sigma[arrow_position(left, quiver)?] = 1;
                            sigma[arrow_position(right, quiver)?] = -1;
                        }
                    }
                    _ => {}
                }
            }
            [only_incoming] => match outgoing.as_slice() {
                [only_outgoing] => {
                    epsilon[arrow_position(only_incoming, quiver)?] = 1;
                    sigma[arrow_position(only_outgoing, quiver)?] = -1;
                }
                [left, right] => {
                    sigma[arrow_position(left, quiver)?] = 1;
                    sigma[arrow_position(right, quiver)?] = -1;
                    epsilon[arrow_position(only_incoming, quiver)?] =
                        if composable(left, only_incoming)? {
                            -1
                        } else {
                            1
                        };
                }
                _ => {}
            },
            [] => {
                let [left, right, ..] = outgoing.as_slice() else {
                    return Err(ValidityError::MissingSourceArrows);
                };
                sigma[arrow_position(left, quiver)?] = 1;
                sigma[arrow_position(right, quiver)?] = -1;
            }
            _ => {}
        }
    }
    Ok(StringSigmaEpsilon { sigma, epsilon })
}

/// Translate QPA's `IsValidString(A,input_str)` Boolean result.
///
/// # Errors
///
/// Preserves QPA's explicit algebra/quiver errors and reports its unchecked
/// list-index failures as typed errors.
pub fn is_valid_string<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
) -> Result<bool, ValidityError> {
    if !algebra.qpa_is_string_algebra()? {
        return Err(ValidityError::NotStringAlgebra);
    }
    let parent = algebra.qpa_underlying_path_algebra();
    let quiver = parent.quiver();
    ensure_acceptable_quiver(quiver)?;
    let relations = relation_paths(algebra)?;
    let signs = qpa_string_sigma_epsilon(algebra)?;
    is_valid_string_after_checks(quiver, &relations, &signs, input)
}

fn is_valid_string_after_checks(
    quiver: &Quiver,
    relations: &[Path],
    signs: &StringSigmaEpsilon,
    input: &str,
) -> Result<bool, ValidityError> {
    let mut letters = input.chars().collect::<Vec<_>>();
    if letters.is_empty() {
        return Err(ValidityError::EmptyInput);
    }
    let alphabet = string_alphabet(quiver)?;
    let mut initial_alphabet = alphabet.clone();
    initial_alphabet.extend("0123456789(-,)".chars());
    if letters
        .iter()
        .any(|character| !initial_alphabet.contains(character))
    {
        return Ok(false);
    }

    if letters[0] == '(' && parse_zero_string(quiver, &letters)?.is_some() {
        return Ok(true);
    }

    remove_lazy_markers(&mut letters, quiver, signs)?;
    if letters
        .iter()
        .any(|character| !alphabet.contains(character))
    {
        return Ok(false);
    }
    if letters.is_empty() {
        return Err(ValidityError::EmptyInput);
    }
    if !is_reduced_walk(&letters) || !is_valid_walk(quiver, &letters)? {
        return Ok(false);
    }

    for relation in relations {
        let Some(walk) = relation.walk() else {
            return Err(ValidityError::EmptyRelation);
        };
        if walk.is_empty() {
            return Err(ValidityError::EmptyRelation);
        }
        let mut reversed = Vec::with_capacity(walk.len());
        let mut inverse = Vec::with_capacity(walk.len());
        for index in walk.iter().rev() {
            reversed.push(arrow_character(quiver, *index)?);
        }
        for index in walk {
            inverse.push(inverse_character(arrow_character(quiver, *index)?)?);
        }
        if contains_char_slice(&letters, &reversed) || contains_char_slice(&letters, &inverse) {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Translate `QPAStringDirectLeft`.
///
/// # Errors
///
/// Propagates QPA validity, alphabet, and incidence failures.
pub fn qpa_string_direct_left<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
    signs: &StringSigmaEpsilon,
) -> Result<QpaStringExtension, ValidityError> {
    let valid = is_valid_string(algebra, input)?;
    let quiver = algebra.qpa_underlying_path_algebra().quiver();
    let letters = input.chars().collect::<Vec<_>>();
    if valid && letters.first() == Some(&'(') {
        let Some((vertex, sign)) = parse_zero_string(quiver, &letters)? else {
            return Err(ValidityError::MalformedZeroString);
        };
        let outgoing = quiver.outgoing_arrows(&vertex)?;
        return choose_signed_arrow(&outgoing, &signs.sigma, -sign, false, quiver);
    }
    let first = *letters.first().ok_or(ValidityError::EmptyInput)?;
    let (underlying, direct) = arrow_for_letter(quiver, first)?;
    let outgoing = if direct {
        quiver.outgoing_arrows(&underlying.target())?
    } else {
        quiver.outgoing_arrows(&underlying.source())?
    };
    if direct {
        try_prepend_direct(algebra, input, &outgoing)
    } else if let [left, right] = outgoing.as_slice() {
        let selected = if inverse_character(arrow_character(quiver, left.index())?)? == first {
            right
        } else {
            left
        };
        validated_prepend(algebra, input, arrow_character(quiver, selected.index())?)
    } else {
        Ok(QpaStringExtension::CannotPerform)
    }
}

/// Translate `QPAStringInverseLeft`.
///
/// # Errors
///
/// Propagates QPA validity, alphabet, and incidence failures.
pub fn qpa_string_inverse_left<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
    signs: &StringSigmaEpsilon,
) -> Result<QpaStringExtension, ValidityError> {
    let valid = is_valid_string(algebra, input)?;
    let quiver = algebra.qpa_underlying_path_algebra().quiver();
    let letters = input.chars().collect::<Vec<_>>();
    if valid && letters.first() == Some(&'(') {
        let Some((vertex, sign)) = parse_zero_string(quiver, &letters)? else {
            return Err(ValidityError::MalformedZeroString);
        };
        let incoming = quiver.incoming_arrows(&vertex)?;
        return choose_signed_arrow(&incoming, &signs.epsilon, -sign, true, quiver);
    }
    let first = *letters.first().ok_or(ValidityError::EmptyInput)?;
    let (underlying, direct) = arrow_for_letter(quiver, first)?;
    let incoming = if direct {
        quiver.incoming_arrows(&underlying.target())?
    } else {
        quiver.incoming_arrows(&underlying.source())?
    };
    if direct {
        if let [left, right] = incoming.as_slice() {
            let selected = if arrow_character(quiver, left.index())? == first {
                right
            } else {
                left
            };
            validated_prepend(
                algebra,
                input,
                inverse_character(arrow_character(quiver, selected.index())?)?,
            )
        } else {
            Ok(QpaStringExtension::CannotPerform)
        }
    } else {
        try_prepend_inverse(algebra, input, &incoming, quiver)
    }
}

/// Translate `QPAStringDirectRight`.
///
/// # Errors
///
/// Propagates QPA validity, alphabet, and incidence failures.
pub fn qpa_string_direct_right<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
    signs: &StringSigmaEpsilon,
) -> Result<QpaStringExtension, ValidityError> {
    let valid = is_valid_string(algebra, input)?;
    let quiver = algebra.qpa_underlying_path_algebra().quiver();
    let letters = input.chars().collect::<Vec<_>>();
    if valid && letters.first() == Some(&'(') {
        let Some((vertex, sign)) = parse_zero_string(quiver, &letters)? else {
            return Err(ValidityError::MalformedZeroString);
        };
        let incoming = quiver.incoming_arrows(&vertex)?;
        if sign == 1 && incoming.len() == 1 {
            let position = arrow_position(&incoming[0], quiver)?;
            if signs.epsilon[position] == 1 {
                return Ok(QpaStringExtension::ArrowList(incoming));
            }
        }
        return choose_signed_arrow(&incoming, &signs.epsilon, sign, false, quiver);
    }
    let last = *letters.last().ok_or(ValidityError::EmptyInput)?;
    let (underlying, direct) = arrow_for_letter(quiver, last)?;
    let incoming = if direct {
        quiver.incoming_arrows(&underlying.source())?
    } else {
        quiver.incoming_arrows(&underlying.target())?
    };
    if direct {
        try_append_direct(algebra, input, &incoming)
    } else if let [left, right] = incoming.as_slice() {
        let selected = if inverse_character(arrow_character(quiver, left.index())?)? == last {
            right
        } else {
            left
        };
        validated_append(algebra, input, arrow_character(quiver, selected.index())?)
    } else {
        Ok(QpaStringExtension::CannotPerform)
    }
}

/// Translate `QPAStringInverseRight`.
///
/// # Errors
///
/// Propagates QPA validity, alphabet, and incidence failures.
pub fn qpa_string_inverse_right<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
    signs: &StringSigmaEpsilon,
) -> Result<QpaStringExtension, ValidityError> {
    let valid = is_valid_string(algebra, input)?;
    let quiver = algebra.qpa_underlying_path_algebra().quiver();
    let letters = input.chars().collect::<Vec<_>>();
    if valid && letters.first() == Some(&'(') {
        let Some((vertex, sign)) = parse_zero_string(quiver, &letters)? else {
            return Err(ValidityError::MalformedZeroString);
        };
        let outgoing = quiver.outgoing_arrows(&vertex)?;
        return choose_signed_arrow(&outgoing, &signs.sigma, sign, true, quiver);
    }
    let last = *letters.last().ok_or(ValidityError::EmptyInput)?;
    let (underlying, direct) = arrow_for_letter(quiver, last)?;
    let outgoing = if direct {
        quiver.outgoing_arrows(&underlying.source())?
    } else {
        quiver.outgoing_arrows(&underlying.target())?
    };
    if direct {
        if let [left, right] = outgoing.as_slice() {
            let selected = if arrow_character(quiver, left.index())? == last {
                right
            } else {
                left
            };
            validated_append(
                algebra,
                input,
                inverse_character(arrow_character(quiver, selected.index())?)?,
            )
        } else {
            Ok(QpaStringExtension::CannotPerform)
        }
    } else {
        try_append_inverse(algebra, input, &outgoing, quiver)
    }
}

/// Translate QPA's `LocalARQuiver(A,input_str)` list of irreducible-map
/// descriptions.
///
/// As in QPA, the four hook/cohook loops have no artificial iteration bound.
///
/// # Errors
///
/// Preserves QPA's algebra, quiver, and invalid-input diagnostics and
/// propagates extension failures.
pub fn local_ar_quiver<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
) -> Result<Vec<String>, ValidityError> {
    if !is_valid_string(algebra, input)? {
        return Err(ValidityError::InvalidString);
    }
    let signs = qpa_string_sigma_epsilon(algebra)?;
    let mut arrows = Vec::new();
    if let Some(word) =
        extend_until_blocked(qpa_string_inverse_left(algebra, input, &signs)?, |word| {
            qpa_string_direct_left(algebra, word, &signs)
        })?
    {
        arrows.push(format!(
            "The canonical embedding StringModule({input}) to StringModule({word})"
        ));
    }
    if let Some(word) =
        extend_until_blocked(qpa_string_direct_left(algebra, input, &signs)?, |word| {
            qpa_string_inverse_left(algebra, word, &signs)
        })?
    {
        arrows.push(format!(
            "The canonical projection StringModule({word}) to StringModule({input})"
        ));
    }
    if let Some(word) =
        extend_until_blocked(qpa_string_inverse_right(algebra, input, &signs)?, |word| {
            qpa_string_direct_right(algebra, word, &signs)
        })?
    {
        arrows.push(format!(
            "The canonical projection StringModule({word}) to StringModule({input})"
        ));
    }
    if let Some(word) =
        extend_until_blocked(qpa_string_direct_right(algebra, input, &signs)?, |word| {
            qpa_string_inverse_right(algebra, word, &signs)
        })?
    {
        arrows.push(format!(
            "The canonical embedding StringModule({input}) to StringModule({word})"
        ));
    }
    Ok(arrows)
}

fn relation_paths<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
) -> Result<Vec<Path>, ValidityError> {
    algebra
        .qpa_defining_relations()
        .iter()
        .map(|relation| {
            relation
                .terms()
                .first()
                .map(|term| term.path().clone())
                .ok_or(ValidityError::EmptyRelation)
        })
        .collect()
}

fn ensure_acceptable_quiver(quiver: &Quiver) -> Result<(), ValidityError> {
    for vertex in quiver.vertices() {
        let decimal = vertex.index().get().to_string();
        let expected = decimal.chars().next().expect("positive vertex index");
        let mut name = vertex.name().chars();
        if name.next() != Some('v') || name.next() != Some(expected) || name.next().is_some() {
            return Err(ValidityError::UnacceptableQuiver);
        }
    }
    for arrow in quiver.arrows() {
        let position = u32::try_from(arrow.index().get())
            .map_err(|_| ValidityError::InvalidArrowCharacter { character: '\0' })?;
        let expected = char::from_u32(position + 96)
            .ok_or(ValidityError::InvalidArrowCharacter { character: '\0' })?;
        let mut name = arrow.name().chars();
        if name.next() != Some(expected) || name.next().is_some() {
            return Err(ValidityError::UnacceptableQuiver);
        }
    }
    Ok(())
}

fn string_alphabet(quiver: &Quiver) -> Result<HashSet<char>, ValidityError> {
    let mut alphabet = HashSet::with_capacity(2 * quiver.number_of_arrows());
    for arrow in quiver.arrows() {
        let direct = arrow_character(quiver, arrow.index())?;
        alphabet.insert(direct);
        alphabet.insert(inverse_character(direct)?);
    }
    Ok(alphabet)
}

fn remove_lazy_markers(
    letters: &mut Vec<char>,
    quiver: &Quiver,
    signs: &StringSigmaEpsilon,
) -> Result<(), ValidityError> {
    let length = letters.len();
    let mut removed = 0_usize;
    if length > 2 {
        for index in 0..length - 2 {
            if letters[index] == '-'
                && letters[index + 1] == '1'
                && letter_position(quiver, letters[index + 2]).is_ok()
            {
                let (position, direct) = letter_position(quiver, letters[index + 2])?;
                if (direct && signs.epsilon[position] == 1)
                    || (!direct && signs.sigma[position] == 1)
                {
                    mark_invalid(letters);
                    return Ok(());
                }
                for shift in index..length - 2 {
                    letters[shift] = letters[shift + 2];
                }
                removed += 2;
                if letters[index + 1] == '1' {
                    letters[index + 1] = '2';
                }
            }
        }
    }
    if length > 1 {
        for index in 0..length - 1 {
            if letters[index] == '1' && letter_position(quiver, letters[index + 1]).is_ok() {
                let (position, direct) = letter_position(quiver, letters[index + 1])?;
                if (direct && signs.epsilon[position] == -1)
                    || (!direct && signs.sigma[position] == -1)
                {
                    mark_invalid(letters);
                    return Ok(());
                }
                for shift in index..length - 1 {
                    letters[shift] = letters[shift + 1];
                }
                removed += 1;
            }
        }
    }
    if length > 1 {
        if letters[length - 2] == '-' && letters[length - 1] == '1' {
            if length < 3 {
                return Err(ValidityError::InvalidLazyStringMarkers);
            }
            if let Ok((position, direct)) = letter_position(quiver, letters[length - 3]) {
                if (direct && signs.sigma[position] == -1)
                    || (!direct && signs.epsilon[position] == -1)
                {
                    mark_invalid(letters);
                    return Ok(());
                }
                removed += 2;
            }
        } else if letters[length - 1] == '1' && letter_position(quiver, letters[length - 2]).is_ok()
        {
            let (position, direct) = letter_position(quiver, letters[length - 2])?;
            if (direct && signs.sigma[position] == 1) || (!direct && signs.epsilon[position] == 1) {
                mark_invalid(letters);
                return Ok(());
            }
            removed += 1;
        }
    }
    let retained = length
        .checked_sub(removed)
        .ok_or(ValidityError::InvalidLazyStringMarkers)?;
    letters.truncate(retained);
    Ok(())
}

// QPA returns `false` immediately for these sign conflicts.  A private marker
// outside the arrow alphabet lets the caller reach that same Boolean branch
// without turning the source result into a Rust error.
fn mark_invalid(letters: &mut Vec<char>) {
    letters.clear();
    letters.push('\0');
}

fn is_reduced_walk(letters: &[char]) -> bool {
    letters
        .windows(2)
        .all(|pair| (pair[0] as i64 - pair[1] as i64).unsigned_abs() != 32)
}

fn is_valid_walk(quiver: &Quiver, letters: &[char]) -> Result<bool, ValidityError> {
    for pair in letters.windows(2) {
        let (left, left_direct) = arrow_for_letter(quiver, pair[0])?;
        let (right, right_direct) = arrow_for_letter(quiver, pair[1])?;
        let composable = match (left_direct, right_direct) {
            (true, true) => left.source() == right.target(),
            (true, false) => left.source() == right.source(),
            (false, true) => left.target() == right.target(),
            (false, false) => left.target() == right.source(),
        };
        if !composable {
            return Ok(false);
        }
    }
    Ok(true)
}

fn parse_zero_string(
    quiver: &Quiver,
    letters: &[char],
) -> Result<Option<(Vertex, i8)>, ValidityError> {
    if letters.first() != Some(&'(') {
        return Ok(None);
    }
    let digits = quiver.number_of_vertices().to_string().len();
    let digit_end = 1 + digits;
    if letters.len() <= digit_end {
        return Err(ValidityError::MalformedZeroString);
    }
    let mut number = 0_usize;
    for character in &letters[1..digit_end] {
        let Some(digit) = character.to_digit(10) else {
            return Ok(None);
        };
        number = number * 10 + digit as usize;
    }
    if letters.get(digit_end) != Some(&',') {
        return Ok(None);
    }
    let sign = match (
        letters.get(digit_end + 1),
        letters.get(digit_end + 2),
        letters.get(digit_end + 3),
    ) {
        (Some('1'), Some(')'), _) => 1,
        (Some('-'), Some('1'), Some(')')) => -1,
        _ => return Ok(None),
    };
    let Some(vertex) = qpa_core::VertexIndex::new(number).and_then(|index| quiver.vertex(index))
    else {
        return Ok(None);
    };
    Ok(Some((vertex, sign)))
}

fn choose_signed_arrow(
    arrows: &[Arrow],
    values: &[i8],
    desired: i8,
    inverse: bool,
    quiver: &Quiver,
) -> Result<QpaStringExtension, ValidityError> {
    let selected = match arrows {
        [only] if values[arrow_position(only, quiver)?] == desired => Some(only),
        [left, right] => {
            if values[arrow_position(left, quiver)?] == desired {
                Some(left)
            } else {
                Some(right)
            }
        }
        _ => None,
    };
    let Some(selected) = selected else {
        return Ok(QpaStringExtension::CannotPerform);
    };
    let character = arrow_character(quiver, selected.index())?;
    let character = if inverse {
        inverse_character(character)?
    } else {
        character
    };
    Ok(QpaStringExtension::Word(character.to_string()))
}

fn try_prepend_direct<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
    arrows: &[Arrow],
) -> Result<QpaStringExtension, ValidityError> {
    for arrow in arrows.iter().take(2) {
        let result = validated_prepend(
            algebra,
            input,
            arrow_character(
                algebra.qpa_underlying_path_algebra().quiver(),
                arrow.index(),
            )?,
        )?;
        if !matches!(result, QpaStringExtension::CannotPerform) {
            return Ok(result);
        }
    }
    Ok(QpaStringExtension::CannotPerform)
}

fn try_prepend_inverse<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
    arrows: &[Arrow],
    quiver: &Quiver,
) -> Result<QpaStringExtension, ValidityError> {
    for arrow in arrows.iter().take(2) {
        let result = validated_prepend(
            algebra,
            input,
            inverse_character(arrow_character(quiver, arrow.index())?)?,
        )?;
        if !matches!(result, QpaStringExtension::CannotPerform) {
            return Ok(result);
        }
    }
    Ok(QpaStringExtension::CannotPerform)
}

fn try_append_direct<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
    arrows: &[Arrow],
) -> Result<QpaStringExtension, ValidityError> {
    for arrow in arrows.iter().take(2) {
        let result = validated_append(
            algebra,
            input,
            arrow_character(
                algebra.qpa_underlying_path_algebra().quiver(),
                arrow.index(),
            )?,
        )?;
        if !matches!(result, QpaStringExtension::CannotPerform) {
            return Ok(result);
        }
    }
    Ok(QpaStringExtension::CannotPerform)
}

fn try_append_inverse<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
    arrows: &[Arrow],
    quiver: &Quiver,
) -> Result<QpaStringExtension, ValidityError> {
    for arrow in arrows.iter().take(2) {
        let result = validated_append(
            algebra,
            input,
            inverse_character(arrow_character(quiver, arrow.index())?)?,
        )?;
        if !matches!(result, QpaStringExtension::CannotPerform) {
            return Ok(result);
        }
    }
    Ok(QpaStringExtension::CannotPerform)
}

fn validated_prepend<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
    character: char,
) -> Result<QpaStringExtension, ValidityError> {
    let candidate = format!("{character}{input}");
    Ok(if is_valid_string(algebra, &candidate)? {
        QpaStringExtension::Word(candidate)
    } else {
        QpaStringExtension::CannotPerform
    })
}

fn validated_append<F: Field, A: QpaStringAlgebra<F>>(
    algebra: &A,
    input: &str,
    character: char,
) -> Result<QpaStringExtension, ValidityError> {
    let candidate = format!("{input}{character}");
    Ok(if is_valid_string(algebra, &candidate)? {
        QpaStringExtension::Word(candidate)
    } else {
        QpaStringExtension::CannotPerform
    })
}

fn extend_until_blocked(
    initial: QpaStringExtension,
    mut extend: impl FnMut(&str) -> Result<QpaStringExtension, ValidityError>,
) -> Result<Option<String>, ValidityError> {
    let mut word = match initial {
        QpaStringExtension::CannotPerform => return Ok(None),
        QpaStringExtension::Word(word) => word,
        QpaStringExtension::ArrowList(_) => return Err(ValidityError::NonStringExtension),
    };
    loop {
        match extend(&word)? {
            QpaStringExtension::CannotPerform => return Ok(Some(word)),
            QpaStringExtension::Word(next) => word = next,
            QpaStringExtension::ArrowList(_) => return Err(ValidityError::NonStringExtension),
        }
    }
}

fn arrow_for_letter(quiver: &Quiver, character: char) -> Result<(Arrow, bool), ValidityError> {
    let (position, direct) = letter_position(quiver, character)?;
    let index =
        ArrowIndex::new(position + 1).ok_or(ValidityError::InvalidArrowCharacter { character })?;
    let arrow = quiver
        .arrow(index)
        .ok_or(ValidityError::InvalidArrowCharacter { character })?;
    Ok((arrow, direct))
}

fn letter_position(quiver: &Quiver, character: char) -> Result<(usize, bool), ValidityError> {
    let code = character as u32;
    let direct = code > 96;
    let one_based = if direct {
        code.checked_sub(96)
    } else {
        code.checked_sub(64)
    }
    .ok_or(ValidityError::InvalidArrowCharacter { character })? as usize;
    if one_based == 0 || one_based > quiver.number_of_arrows() {
        return Err(ValidityError::InvalidArrowCharacter { character });
    }
    Ok((one_based - 1, direct))
}

fn arrow_position(arrow: &Arrow, quiver: &Quiver) -> Result<usize, ValidityError> {
    let character = arrow_character(quiver, arrow.index())?;
    Ok(letter_position(quiver, character)?.0)
}

fn arrow_character(quiver: &Quiver, index: ArrowIndex) -> Result<char, ValidityError> {
    let arrow = quiver
        .arrow(index)
        .ok_or(ValidityError::InvalidArrowCharacter { character: '\0' })?;
    first_character(arrow.name()).ok_or(ValidityError::InvalidArrowCharacter { character: '\0' })
}

fn inverse_character(character: char) -> Result<char, ValidityError> {
    char::from_u32(
        (character as u32)
            .checked_sub(32)
            .ok_or(ValidityError::InvalidArrowCharacter { character })?,
    )
    .ok_or(ValidityError::InvalidArrowCharacter { character })
}

fn first_character(name: &str) -> Option<char> {
    name.chars().next()
}

fn contains_char_slice(haystack: &[char], needle: &[char]) -> bool {
    !needle.is_empty()
        && haystack
            .windows(needle.len())
            .any(|window| window == needle)
}

#[cfg(test)]
mod tests {
    use numerica::domains::{Field, finite_field::Z2, rational::Q};
    use qpa_core::{ArrowSpec, PathAlgebra, Quiver, VertexIndex};
    use serde_json::{Value, json};

    use crate::quotient_by_relations;

    use super::{
        QpaStringExtension, is_valid_string, local_ar_quiver, qpa_string_direct_left,
        qpa_string_direct_right, qpa_string_inverse_left, qpa_string_inverse_right,
        qpa_string_quiver_rho, qpa_string_sigma_epsilon,
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
        let signs = qpa_string_sigma_epsilon(&algebra).unwrap();
        let inputs = [
            "eca", "Ab", "bc", "CB", "cb", "BC", "dc", "CD", "(1,1)", "(4,-1)", "(5,1)", "aa",
            "aA", "ca", "1a", "-1a", "a1", "a-1",
        ];
        let extension_text = |extension: QpaStringExtension| {
            extension
                .as_qpa_text()
                .expect("fixture branches return text")
                .to_owned()
        };
        let extensions = [
            qpa_string_direct_left(&algebra, "eca", &signs).unwrap(),
            qpa_string_inverse_left(&algebra, "eca", &signs).unwrap(),
            qpa_string_direct_right(&algebra, "eca", &signs).unwrap(),
            qpa_string_inverse_right(&algebra, "eca", &signs).unwrap(),
        ]
        .map(extension_text);
        let zero_extensions = [
            qpa_string_direct_left(&algebra, "(1,1)", &signs).unwrap(),
            qpa_string_inverse_left(&algebra, "(2,-1)", &signs).unwrap(),
            qpa_string_direct_right(&algebra, "(4,1)", &signs).unwrap(),
            qpa_string_inverse_right(&algebra, "(3,-1)", &signs).unwrap(),
        ]
        .map(extension_text);
        let rho = qpa_string_quiver_rho(&quiver, &["bc".to_owned(), "cd".to_owned()])
            .unwrap()
            .iter()
            .map(|path| {
                path.walk()
                    .unwrap()
                    .iter()
                    .map(|index| quiver.arrow(*index).unwrap().name().to_owned())
                    .collect::<String>()
            })
            .collect::<Vec<_>>();
        json!({
            "rho": rho,
            "sigma": signs.sigma(),
            "epsilon": signs.epsilon(),
            "inputs": inputs,
            "valid": inputs.map(|input| is_valid_string(&algebra, input).unwrap()),
            "extensions": extensions,
            "zero_extensions": zero_extensions,
            "local_ar": local_ar_quiver(&algebra, "eca").unwrap()
        })
    }

    #[test]
    fn pinned_qpa_validity_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/validity.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(fixture_case(Z2), expected["binary"]);
        assert_eq!(fixture_case(Q), expected["rational"]);
    }
}
