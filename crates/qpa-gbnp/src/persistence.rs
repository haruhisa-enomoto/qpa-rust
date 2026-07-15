//! QPA-compatible `SaveAlgebra` and safe `ReadAlgebra` text persistence.
//!
//! The writer intentionally reproduces QPA 1.37's line-oriented format. The
//! reader parses that format directly; unlike GAP it never evaluates input as
//! code. The initial inferred-field boundary is the exact rational field and
//! Numerica's 32-bit prime fields, including its specialized binary field.

use std::{collections::HashMap, fs, path::Path, str::FromStr};

use numerica::domains::{
    Field, Ring,
    finite_field::{FiniteFieldCore, PrimitiveRootIterator, Z2, Zp},
    integer::Integer,
    rational::{Q, RationalField},
};
use qpa_core::{
    AlgebraDimension, NamedArrowSpec, PathAlgebra, PathAlgebraError, Quiver, QuiverAlgebra,
    QuiverError,
};
use thiserror::Error;

use crate::{PathalgGbnpError, quotient_by_relations};

/// Exact coefficient fields that QPA's text format can write in this port.
///
/// QPA writes rational numbers literally and nonzero finite-field elements as
/// powers of `PrimitiveRoot(K)`. Implementations preserve those encodings.
pub trait QpaTextField: Field {
    /// Return QPA's field line (`Rationals` or `GF(p)`).
    ///
    /// # Errors
    ///
    /// Returns an error when the field lies outside the admitted QPA text
    /// contract.
    fn qpa_field_text(&self) -> Result<String, QpaPersistenceError>;

    /// Format one nonzero sparse coefficient exactly as QPA `SaveAlgebra`.
    ///
    /// # Errors
    ///
    /// Returns an error for zero sparse coefficients or a field whose QPA
    /// primitive-root encoding cannot be constructed.
    fn qpa_coefficient_text(
        &self,
        coefficient: &Self::Element,
    ) -> Result<String, QpaPersistenceError>;

    /// Parse one coefficient emitted for this field by QPA `SaveAlgebra`.
    ///
    /// # Errors
    ///
    /// Returns an error when the token is not a valid QPA coefficient for this
    /// field.
    fn parse_qpa_coefficient(
        &self,
        coefficient: &str,
    ) -> Result<Self::Element, QpaPersistenceError>;
}

impl QpaTextField for RationalField {
    fn qpa_field_text(&self) -> Result<String, QpaPersistenceError> {
        Ok("Rationals".to_owned())
    }

    fn qpa_coefficient_text(
        &self,
        coefficient: &Self::Element,
    ) -> Result<String, QpaPersistenceError> {
        if coefficient.denominator_ref() == &Integer::from(1) {
            Ok(coefficient.numerator_ref().to_string())
        } else {
            Ok(format!(
                "{}/{}",
                coefficient.numerator_ref(),
                coefficient.denominator_ref()
            ))
        }
    }

    fn parse_qpa_coefficient(
        &self,
        coefficient: &str,
    ) -> Result<Self::Element, QpaPersistenceError> {
        let (numerator, denominator) = coefficient
            .split_once('/')
            .map_or((coefficient, "1"), |parts| parts);
        let numerator =
            Integer::from_str(numerator).map_err(|_| QpaPersistenceError::InvalidCoefficient {
                coefficient: coefficient.into(),
            })?;
        let denominator = Integer::from_str(denominator).map_err(|_| {
            QpaPersistenceError::InvalidCoefficient {
                coefficient: coefficient.into(),
            }
        })?;
        if denominator == 0 {
            return Err(QpaPersistenceError::ZeroRationalDenominator);
        }
        Ok(self.to_element(numerator, denominator, true))
    }
}

impl QpaTextField for Z2 {
    fn qpa_field_text(&self) -> Result<String, QpaPersistenceError> {
        Ok("GF(2)".to_owned())
    }

    fn qpa_coefficient_text(
        &self,
        coefficient: &Self::Element,
    ) -> Result<String, QpaPersistenceError> {
        if self.is_zero(coefficient) {
            Err(QpaPersistenceError::ZeroSparseCoefficient)
        } else {
            Ok("z^0".to_owned())
        }
    }

    fn parse_qpa_coefficient(
        &self,
        coefficient: &str,
    ) -> Result<Self::Element, QpaPersistenceError> {
        let exponent = parse_primitive_power(coefficient)?;
        if exponent == 0 {
            Ok(self.one())
        } else {
            Ok(self.pow(&self.one(), exponent))
        }
    }
}

impl QpaTextField for Zp {
    fn qpa_field_text(&self) -> Result<String, QpaPersistenceError> {
        ensure_supported_odd_prime(self)?;
        Ok(format!("GF({})", self.get_prime()))
    }

    fn qpa_coefficient_text(
        &self,
        coefficient: &Self::Element,
    ) -> Result<String, QpaPersistenceError> {
        ensure_supported_odd_prime(self)?;
        if self.is_zero(coefficient) {
            return Err(QpaPersistenceError::ZeroSparseCoefficient);
        }
        let characteristic = self.get_prime();
        let prime = u64::from(characteristic);
        let primitive_root = u64::from(primitive_root(characteristic)?);
        let residue = u64::from(self.from_element(coefficient));
        let exponent = discrete_log_prime_field(primitive_root, residue, prime).ok_or(
            QpaPersistenceError::FiniteFieldDiscreteLog {
                residue,
                characteristic: prime,
            },
        )?;
        Ok(format!("z^{exponent}"))
    }

    fn parse_qpa_coefficient(
        &self,
        coefficient: &str,
    ) -> Result<Self::Element, QpaPersistenceError> {
        ensure_supported_odd_prime(self)?;
        let exponent = parse_primitive_power(coefficient)?;
        let root = primitive_root(self.get_prime())?;
        Ok(self.pow(&self.to_element(root), exponent))
    }
}

/// An algebra inferred by QPA `ReadAlgebra` over the initially supported
/// Numerica coefficient domains.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum LoadedQpaAlgebra {
    /// An algebra over the exact rational field.
    Rational(QuiverAlgebra<RationalField>),
    /// An algebra over Numerica's specialized `GF(2)` field.
    Binary(QuiverAlgebra<Z2>),
    /// An algebra over an odd 32-bit prime field.
    Prime(QuiverAlgebra<Zp>),
}

impl LoadedQpaAlgebra {
    /// Return whether the saved type reconstructed an actual path algebra.
    #[must_use]
    pub const fn is_path_algebra(&self) -> bool {
        match self {
            Self::Rational(algebra) => algebra.is_path_algebra(),
            Self::Binary(algebra) => algebra.is_path_algebra(),
            Self::Prime(algebra) => algebra.is_path_algebra(),
        }
    }

    /// Return the finite dimension when the loaded value is a quotient with a
    /// complete finite-nontip basis.
    ///
    /// # Errors
    ///
    /// Propagates quotient dimension failures; path algebras return `None`.
    pub fn quotient_dimension(&self) -> Result<Option<AlgebraDimension>, QpaPersistenceError> {
        match self {
            Self::Rational(algebra) => algebra
                .as_quotient()
                .map(qpa_core::QuotientPathAlgebra::dimension)
                .transpose(),
            Self::Binary(algebra) => algebra
                .as_quotient()
                .map(qpa_core::QuotientPathAlgebra::dimension)
                .transpose(),
            Self::Prime(algebra) => algebra
                .as_quotient()
                .map(qpa_core::QuotientPathAlgebra::dimension)
                .transpose(),
        }
        .map_err(Into::into)
    }
}

/// Render the exact line-oriented text produced by QPA `SaveAlgebra` after
/// its overwrite/file-opening branch.
///
/// # Errors
///
/// Returns an error for a coefficient field outside the supported text
/// contract or an invalid sparse zero coefficient.
pub fn qpa_algebra_text<F: QpaTextField>(
    algebra: &QuiverAlgebra<F>,
) -> Result<String, QpaPersistenceError> {
    let parent = algebra.original_path_algebra();
    let quiver = parent.quiver();
    let relations = algebra.relations();
    let mut output = String::new();
    output.push_str(if relations.is_empty() {
        "IsPathAlgebra\n"
    } else {
        "IsQuotientOfPathAlgebra\n"
    });
    output.push_str("Vertices:\n");
    for (position, vertex) in quiver.vertices().enumerate() {
        if position > 0 {
            if position % 10 == 0 {
                output.push_str(",\n");
            } else {
                output.push_str(", ");
            }
        }
        output.push_str(vertex.name());
    }
    output.push('\n');
    output.push_str("Arrows:\n");
    for (position, arrow) in quiver.arrows().enumerate() {
        if position > 0 {
            if position % 6 == 0 {
                output.push_str(",\n");
            } else {
                output.push_str(", ");
            }
        }
        output.push_str(arrow.name());
        output.push(':');
        output.push_str(arrow.source().name());
        output.push_str("->");
        output.push_str(arrow.target().name());
    }
    if quiver.number_of_arrows() > 0 {
        output.push('\n');
    }
    output.push_str("Field:\n");
    output.push_str(&parent.field().qpa_field_text()?);
    output.push('\n');

    if !relations.is_empty() {
        output.push_str("Relations:\n");
        for (relation_position, relation) in relations.iter().enumerate() {
            for (term_position, term) in relation.terms().iter().enumerate() {
                if term_position > 0 {
                    output.push_str(" + ");
                }
                output.push('(');
                output.push_str(&parent.field().qpa_coefficient_text(term.coefficient())?);
                output.push_str(")*");
                if let Some(walk) = term.path().walk() {
                    for (arrow_position, index) in walk.iter().enumerate() {
                        if arrow_position > 0 {
                            output.push('*');
                        }
                        let arrow = quiver.arrow(*index).ok_or(
                            QpaPersistenceError::InvalidRelationArrowIndex { index: index.get() },
                        )?;
                        output.push_str(arrow.name());
                    }
                }
            }
            if relation_position + 1 < relations.len() {
                output.push(',');
            }
            output.push('\n');
        }
    }
    Ok(output)
}

/// Translate QPA `SaveAlgebra(A, file, overwrite)` with an explicit callback
/// for its interactive `query` branch.
///
/// The callback is invoked only when the normalized option starts with
/// `query` and the target exists. `Ok(false)` declines replacement. QPA's
/// source bug for `query` on a new file is returned as a typed error.
///
/// # Errors
///
/// Propagates format, callback, and filesystem failures.
pub fn save_qpa_algebra_file<F, Query>(
    algebra: &QuiverAlgebra<F>,
    file: impl AsRef<Path>,
    overwrite: &str,
    query_existing: Query,
) -> Result<bool, QpaPersistenceError>
where
    F: QpaTextField,
    Query: FnOnce() -> Result<bool, QpaPersistenceError>,
{
    let text = qpa_algebra_text(algebra)?;
    let option = normalize_whitespace(overwrite);
    if overwrite.chars().count() < 4 {
        return Err(QpaPersistenceError::InvalidOverwriteOption {
            option: overwrite.into(),
        });
    }
    let file = file.as_ref();
    let exists = file.try_exists()?;
    if option.starts_with("keep") {
        if exists {
            return Ok(false);
        }
    } else if option.starts_with("query") {
        if !exists {
            return Err(QpaPersistenceError::QueryNewFileUnboundOutput);
        }
        if !query_existing()? {
            return Ok(false);
        }
    } else if !option.starts_with("delete") {
        return Err(QpaPersistenceError::InvalidOverwriteOption {
            option: overwrite.into(),
        });
    }
    fs::write(file, text)?;
    Ok(true)
}

/// Safely parse QPA `SaveAlgebra` text and reconstruct its algebra through the
/// configured default GBNP quotient backend.
///
/// # Errors
///
/// Returns stable format/field errors and propagates quiver, path-algebra, and
/// default GBNP failures.
pub fn parse_qpa_algebra_text(text: &str) -> Result<LoadedQpaAlgebra, QpaPersistenceError> {
    let presentation = ParsedPresentation::parse(text)?;
    match presentation.field {
        ParsedField::Rational => Ok(LoadedQpaAlgebra::Rational(presentation.instantiate(Q)?)),
        ParsedField::Prime(2) => Ok(LoadedQpaAlgebra::Binary(presentation.instantiate(Z2)?)),
        ParsedField::Prime(characteristic) => {
            if characteristic < 3 || characteristic.is_multiple_of(2) {
                return Err(QpaPersistenceError::UnsupportedPrimeField { characteristic });
            }
            let field = Zp::new(characteristic);
            ensure_supported_odd_prime(&field)?;
            Ok(LoadedQpaAlgebra::Prime(presentation.instantiate(field)?))
        }
    }
}

/// Translate QPA `ReadAlgebra(file)` using the safe text parser.
///
/// # Errors
///
/// Propagates file reading and [`parse_qpa_algebra_text`] failures.
pub fn read_qpa_algebra_file(
    file: impl AsRef<Path>,
) -> Result<LoadedQpaAlgebra, QpaPersistenceError> {
    parse_qpa_algebra_text(&fs::read_to_string(file)?)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SavedAlgebraType {
    Path,
    Quotient,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ParsedField {
    Rational,
    Prime(u32),
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ParsedArrow {
    source: Box<str>,
    target: Box<str>,
    name: Box<str>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ParsedTerm {
    coefficient: Box<str>,
    arrows: Vec<Box<str>>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ParsedPresentation {
    algebra_type: SavedAlgebraType,
    vertices: Vec<Box<str>>,
    arrows: Vec<ParsedArrow>,
    field: ParsedField,
    relations: Vec<Vec<ParsedTerm>>,
}

impl ParsedPresentation {
    fn parse(text: &str) -> Result<Self, QpaPersistenceError> {
        let lines = text.lines().collect::<Vec<_>>();
        let mut position = 0_usize;
        let algebra_type = match normalized_line(&lines, &mut position, "algebra type")?.as_str() {
            "IsPathAlgebra" => SavedAlgebraType::Path,
            "IsQuotientOfPathAlgebra" => SavedAlgebraType::Quotient,
            value => {
                return Err(QpaPersistenceError::UnknownAlgebraType {
                    value: value.into(),
                });
            }
        };
        require_header_prefix(&lines, &mut position, "Vertices")?;
        let vertices = parse_vertices(&lines, &mut position)?;
        let arrows = parse_arrows(&lines, &mut position)?;
        let field = parse_field(&lines, &mut position)?;
        let relations = if algebra_type == SavedAlgebraType::Path {
            Vec::new()
        } else {
            require_header_prefix(&lines, &mut position, "Relations")?;
            parse_relations(&lines[position..])?
        };

        Ok(Self {
            algebra_type,
            vertices,
            arrows,
            field,
            relations,
        })
    }

    fn instantiate<F: QpaTextField>(
        &self,
        field: F,
    ) -> Result<QuiverAlgebra<F>, QpaPersistenceError> {
        let quiver = Quiver::from_named(
            self.vertices.clone(),
            self.arrows
                .iter()
                .map(|arrow| {
                    NamedArrowSpec::named(
                        arrow.source.clone(),
                        arrow.target.clone(),
                        arrow.name.clone(),
                    )
                })
                .collect(),
        )?;
        let parent = PathAlgebra::new(field, quiver.clone())?;
        if self.algebra_type == SavedAlgebraType::Path || self.relations.is_empty() {
            return Ok(QuiverAlgebra::Path(parent));
        }
        let arrows_by_name = quiver
            .arrows()
            .map(|arrow| (Box::<str>::from(arrow.name()), arrow.path()))
            .collect::<HashMap<_, _>>();
        let mut relations = Vec::with_capacity(self.relations.len());
        for relation in &self.relations {
            let mut value = parent.zero();
            for term in relation {
                let coefficient = parent.field().parse_qpa_coefficient(&term.coefficient)?;
                let mut monomial = parent.one();
                for arrow_name in &term.arrows {
                    let path = arrows_by_name.get(arrow_name).ok_or_else(|| {
                        QpaPersistenceError::UnknownArrowName {
                            name: arrow_name.clone(),
                        }
                    })?;
                    monomial = monomial.multiply(&parent.element(path.clone())?)?;
                }
                value = value.add(&monomial.scale_left(&coefficient)?)?;
            }
            relations.push(value);
        }
        Ok(QuiverAlgebra::Quotient(quotient_by_relations(
            parent, &relations,
        )?))
    }
}

fn parse_vertices(
    lines: &[&str],
    position: &mut usize,
) -> Result<Vec<Box<str>>, QpaPersistenceError> {
    let mut vertices = Vec::new();
    loop {
        let line = lines
            .get(*position)
            .ok_or(QpaPersistenceError::UnexpectedEnd {
                section: "vertices",
            })?;
        *position += 1;
        if normalize_whitespace(line).starts_with("Arrows") {
            break;
        }
        vertices.extend(
            line.split(',')
                .map(str::trim)
                .filter(|name| !name.is_empty())
                .map(Box::<str>::from),
        );
    }
    if vertices.is_empty() {
        Err(QpaPersistenceError::NoVertices)
    } else {
        Ok(vertices)
    }
}

fn parse_arrows(
    lines: &[&str],
    position: &mut usize,
) -> Result<Vec<ParsedArrow>, QpaPersistenceError> {
    let mut arrows = Vec::new();
    loop {
        let line = lines
            .get(*position)
            .ok_or(QpaPersistenceError::UnexpectedEnd { section: "arrows" })?;
        *position += 1;
        if normalize_whitespace(line).starts_with("Field") {
            return Ok(arrows);
        }
        for arrow in line
            .split(',')
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            let (name, endpoints) =
                arrow
                    .split_once(':')
                    .ok_or_else(|| QpaPersistenceError::InvalidArrow {
                        arrow: arrow.into(),
                    })?;
            let (source, target) =
                endpoints
                    .split_once("->")
                    .ok_or_else(|| QpaPersistenceError::InvalidArrow {
                        arrow: arrow.into(),
                    })?;
            arrows.push(ParsedArrow {
                source: source.trim().into(),
                target: target.trim().into(),
                name: name.trim().into(),
            });
        }
    }
}

fn parse_field(lines: &[&str], position: &mut usize) -> Result<ParsedField, QpaPersistenceError> {
    let field = normalized_line(lines, position, "field")?;
    if field == "Rationals" {
        return Ok(ParsedField::Rational);
    }
    if let Some(characteristic) = field
        .strip_prefix("GF(")
        .and_then(|value| value.strip_suffix(')'))
    {
        return Ok(ParsedField::Prime(characteristic.parse().map_err(
            |_| QpaPersistenceError::UnsupportedField {
                field: field.clone().into(),
            },
        )?));
    }
    Err(QpaPersistenceError::UnsupportedField {
        field: field.into(),
    })
}

fn parse_relations(lines: &[&str]) -> Result<Vec<Vec<ParsedTerm>>, QpaPersistenceError> {
    let mut relations = Vec::new();
    for line in lines {
        let compact = line
            .chars()
            .filter(|character| !character.is_whitespace() && !matches!(character, '\'' | '\\'))
            .collect::<String>();
        for relation in compact.split(',').filter(|relation| !relation.is_empty()) {
            relations.push(parse_relation(relation)?);
        }
    }
    Ok(relations)
}

fn parse_relation(relation: &str) -> Result<Vec<ParsedTerm>, QpaPersistenceError> {
    let mut terms = Vec::new();
    for term in relation.split('+').filter(|term| !term.is_empty()) {
        let mut factors = term.split('*');
        let coefficient = factors
            .next()
            .and_then(|value| value.strip_prefix('('))
            .and_then(|value| value.strip_suffix(')'))
            .ok_or_else(|| QpaPersistenceError::InvalidRelation {
                relation: relation.into(),
            })?;
        terms.push(ParsedTerm {
            coefficient: coefficient.into(),
            arrows: factors
                .filter(|factor| !factor.is_empty())
                .map(Box::<str>::from)
                .collect(),
        });
    }
    if terms.is_empty() {
        Err(QpaPersistenceError::InvalidRelation {
            relation: relation.into(),
        })
    } else {
        Ok(terms)
    }
}

fn require_header_prefix(
    lines: &[&str],
    position: &mut usize,
    expected: &'static str,
) -> Result<(), QpaPersistenceError> {
    let actual = normalized_line(lines, position, expected)?;
    if actual.starts_with(expected) {
        Ok(())
    } else {
        Err(QpaPersistenceError::MissingHeader {
            expected,
            actual: actual.into(),
        })
    }
}

fn normalized_line(
    lines: &[&str],
    position: &mut usize,
    section: &'static str,
) -> Result<String, QpaPersistenceError> {
    let line = lines
        .get(*position)
        .ok_or(QpaPersistenceError::UnexpectedEnd { section })?;
    *position += 1;
    Ok(normalize_whitespace(line))
}

fn normalize_whitespace(value: &str) -> String {
    value.split_whitespace().collect::<Vec<_>>().join(" ")
}

fn ensure_supported_odd_prime(field: &Zp) -> Result<(), QpaPersistenceError> {
    let characteristic = field.get_prime();
    if characteristic >= 3 && !characteristic.is_multiple_of(2) && field.is_prime_field(0) {
        Ok(())
    } else {
        Err(QpaPersistenceError::UnsupportedPrimeField { characteristic })
    }
}

fn primitive_root(characteristic: u32) -> Result<u32, QpaPersistenceError> {
    let characteristic64 = u64::from(characteristic);
    let root = PrimitiveRootIterator::new(characteristic64).next().ok_or(
        QpaPersistenceError::MissingPrimitiveRoot {
            characteristic: characteristic64,
        },
    )?;
    u32::try_from(root).map_err(|_| QpaPersistenceError::MissingPrimitiveRoot {
        characteristic: characteristic64,
    })
}

fn parse_primitive_power(value: &str) -> Result<u64, QpaPersistenceError> {
    value
        .strip_prefix("z^")
        .ok_or_else(|| QpaPersistenceError::InvalidCoefficient {
            coefficient: value.into(),
        })?
        .parse()
        .map_err(|_| QpaPersistenceError::InvalidCoefficient {
            coefficient: value.into(),
        })
}

fn discrete_log_prime_field(generator: u64, value: u64, prime: u64) -> Option<u64> {
    let order = prime.checked_sub(1)?;
    let step = integer_sqrt_ceil(order);
    let mut baby_steps = HashMap::with_capacity(usize::try_from(step).ok()?);
    let mut current = 1_u64;
    for exponent in 0..step {
        baby_steps.entry(current).or_insert(exponent);
        current = mul_mod(current, generator, prime);
    }
    let inverse_step = pow_mod(generator, order - (step % order), prime);
    let mut giant = value;
    for giant_exponent in 0..=step {
        if let Some(&baby_exponent) = baby_steps.get(&giant) {
            let exponent = giant_exponent
                .checked_mul(step)?
                .checked_add(baby_exponent)?;
            if exponent < order {
                return Some(exponent);
            }
        }
        giant = mul_mod(giant, inverse_step, prime);
    }
    None
}

fn integer_sqrt_ceil(value: u64) -> u64 {
    if value <= 1 {
        return value;
    }
    let mut low = 1_u64;
    let mut high = value.min(1_u64 << 32);
    while low < high {
        let middle = low + (high - low) / 2;
        if middle >= value.div_ceil(middle) {
            high = middle;
        } else {
            low = middle + 1;
        }
    }
    low
}

fn pow_mod(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1_u64;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mul_mod(result, base, modulus);
        }
        exponent >>= 1;
        if exponent > 0 {
            base = mul_mod(base, base, modulus);
        }
    }
    result
}

fn mul_mod(left: u64, right: u64, modulus: u64) -> u64 {
    u64::try_from((u128::from(left) * u128::from(right)) % u128::from(modulus))
        .expect("a residue modulo a u64 modulus always fits in u64")
}

/// Failure in QPA-compatible algebra text persistence.
#[derive(Debug, Error)]
pub enum QpaPersistenceError {
    /// Filesystem access failed.
    #[error(transparent)]
    Io(#[from] std::io::Error),
    /// QPA rejects overwrite strings shorter than four characters.
    #[error("invalid SaveAlgebra overwrite option: {option}")]
    InvalidOverwriteOption {
        /// Rejected source string.
        option: Box<str>,
    },
    /// QPA's `query` branch fails to bind its output stream for a new file.
    #[error("QPA query overwrite on a new file leaves output unbound")]
    QueryNewFileUnboundOutput,
    /// The first line is not one of QPA's two saved algebra types.
    #[error("unknown saved QPA algebra type: {value}")]
    UnknownAlgebraType {
        /// Rejected type line.
        value: Box<str>,
    },
    /// A required line was absent.
    #[error("unexpected end of QPA algebra data while reading {section}")]
    UnexpectedEnd {
        /// Section being parsed.
        section: &'static str,
    },
    /// A required QPA section header was not present.
    #[error("expected QPA header {expected}, found {actual}")]
    MissingHeader {
        /// Required header.
        expected: &'static str,
        /// Actual normalized line.
        actual: Box<str>,
    },
    /// QPA quivers require at least one saved vertex.
    #[error("saved QPA algebra has no vertices")]
    NoVertices,
    /// An arrow record lacks QPA's `name:source->target` form.
    #[error("invalid saved QPA arrow: {arrow}")]
    InvalidArrow {
        /// Rejected arrow record.
        arrow: Box<str>,
    },
    /// The field line is outside the initial persistence contract.
    #[error("unsupported saved QPA field: {field}")]
    UnsupportedField {
        /// Rejected field line.
        field: Box<str>,
    },
    /// `GF(n)` is not a supported 32-bit prime field.
    #[error("GF({characteristic}) is not a supported prime field")]
    UnsupportedPrimeField {
        /// Rejected characteristic.
        characteristic: u32,
    },
    /// A relation does not follow QPA's coefficient/walk grammar.
    #[error("invalid saved QPA relation: {relation}")]
    InvalidRelation {
        /// Rejected relation text.
        relation: Box<str>,
    },
    /// A coefficient token could not be decoded.
    #[error("invalid saved QPA coefficient: {coefficient}")]
    InvalidCoefficient {
        /// Rejected coefficient token.
        coefficient: Box<str>,
    },
    /// A rational coefficient had denominator zero.
    #[error("saved QPA rational coefficient has denominator zero")]
    ZeroRationalDenominator,
    /// Sparse element invariants exclude zero coefficients from saved terms.
    #[error("a saved sparse relation unexpectedly contains a zero coefficient")]
    ZeroSparseCoefficient,
    /// A relation names an arrow absent from the reconstructed quiver.
    #[error("saved QPA relation names unknown arrow {name}")]
    UnknownArrowName {
        /// Missing name.
        name: Box<str>,
    },
    /// A stored path referenced an invalid local arrow position.
    #[error("saved QPA relation has invalid arrow index {index}")]
    InvalidRelationArrowIndex {
        /// One-based arrow index.
        index: usize,
    },
    /// Numerica could not supply a primitive root for the prime field.
    #[error("no primitive root found for GF({characteristic})")]
    MissingPrimitiveRoot {
        /// Prime characteristic.
        characteristic: u64,
    },
    /// A nonzero field residue had no logarithm to the selected primitive root.
    #[error("no primitive-root logarithm for {residue} in GF({characteristic})")]
    FiniteFieldDiscreteLog {
        /// Nonzero standard residue.
        residue: u64,
        /// Prime characteristic.
        characteristic: u64,
    },
    /// Named-quiver reconstruction failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// Path-algebra construction or arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Default GBNP quotient reconstruction failed.
    #[error(transparent)]
    Gbnp(#[from] PathalgGbnpError),
    /// Quotient dimension access failed.
    #[error(transparent)]
    Quotient(#[from] qpa_core::QuotientError),
}

#[cfg(test)]
mod tests {
    use std::{fs, path::PathBuf};

    use numerica::domains::{
        Ring,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::Q,
    };
    use qpa_core::{AlgebraDimension, NamedArrowSpec, PathAlgebra, Quiver, QuiverAlgebra};
    use serde_json::Value;

    use super::{
        LoadedQpaAlgebra, QpaPersistenceError, QpaTextField, parse_qpa_algebra_text,
        qpa_algebra_text, save_qpa_algebra_file,
    };
    use crate::quotient_by_relations;

    fn chain_quiver() -> Quiver {
        Quiver::from_named(
            vec!["v1".into(), "v2".into(), "v3".into()],
            vec![
                NamedArrowSpec::named("v1", "v2", "a"),
                NamedArrowSpec::named("v2", "v3", "b"),
                NamedArrowSpec::named("v1", "v3", "c"),
            ],
        )
        .unwrap()
    }

    fn path_fixture() -> QuiverAlgebra<numerica::domains::rational::RationalField> {
        let vertices = (1..=11)
            .map(|index| index.to_string().into_boxed_str())
            .collect::<Vec<_>>();
        let arrows = (1..=10)
            .map(|index| {
                NamedArrowSpec::named(
                    index.to_string(),
                    (index + 1).to_string(),
                    format!("a_{index}"),
                )
            })
            .collect();
        QuiverAlgebra::Path(
            PathAlgebra::new(Q, Quiver::from_named(vertices, arrows).unwrap()).unwrap(),
        )
    }

    #[test]
    fn pinned_qpa_persistence_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_persistence.schema-v1.json"
        ))
        .unwrap();
        let path = path_fixture();
        assert_eq!(qpa_algebra_text(&path).unwrap(), expected["path"]["text"]);
        let loaded_path =
            parse_qpa_algebra_text(expected["path"]["text"].as_str().unwrap()).unwrap();
        let LoadedQpaAlgebra::Rational(loaded_path) = loaded_path else {
            panic!("QPA path fixture must infer the rational field");
        };
        assert!(loaded_path.is_path_algebra());
        assert_eq!(
            loaded_path
                .original_path_algebra()
                .quiver()
                .number_of_vertices(),
            usize::try_from(expected["path"]["read_vertices"].as_u64().unwrap()).unwrap()
        );
        assert_eq!(
            loaded_path
                .original_path_algebra()
                .quiver()
                .number_of_arrows(),
            usize::try_from(expected["path"]["read_arrows"].as_u64().unwrap()).unwrap()
        );

        let quiver = chain_quiver();
        let parent = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let generators = parent.generators();
        let half = Q.to_element(Integer::from(1), Integer::from(2), true);
        let relation = generators[3]
            .multiply(&generators[4])
            .unwrap()
            .sub(&generators[5].scale_left(&half).unwrap())
            .unwrap();
        let rational = QuiverAlgebra::Quotient(quotient_by_relations(parent, &[relation]).unwrap());
        assert_eq!(
            qpa_algebra_text(&rational).unwrap(),
            expected["rational"]["text"]
        );

        let parent2 = PathAlgebra::new(Z2, quiver.clone()).unwrap();
        let generators2 = parent2.generators();
        let relation2 = generators2[3]
            .multiply(&generators2[4])
            .unwrap()
            .add(&generators2[5])
            .unwrap();
        let binary = QuiverAlgebra::Quotient(quotient_by_relations(parent2, &[relation2]).unwrap());
        assert_eq!(qpa_algebra_text(&binary).unwrap(), expected["f2"]["text"]);

        let field5 = Zp::new(5);
        let parent5 = PathAlgebra::new(field5.clone(), quiver).unwrap();
        let generators5 = parent5.generators();
        let relation5 = generators5[3]
            .multiply(&generators5[4])
            .unwrap()
            .add(&generators5[5].scale_left(&field5.to_element(2)).unwrap())
            .unwrap();
        let prime = QuiverAlgebra::Quotient(quotient_by_relations(parent5, &[relation5]).unwrap());
        assert_eq!(qpa_algebra_text(&prime).unwrap(), expected["f5"]["text"]);

        for (fixture, expected_variant) in [
            (&expected["rational"], "rational"),
            (&expected["f2"], "binary"),
            (&expected["f5"], "prime"),
        ] {
            let loaded = parse_qpa_algebra_text(fixture["text"].as_str().unwrap()).unwrap();
            assert_eq!(
                loaded.quotient_dimension().unwrap(),
                Some(AlgebraDimension::Finite(Integer::from(
                    fixture["read_dimension"].as_i64().unwrap()
                )))
            );
            assert!(matches!(
                (&loaded, expected_variant),
                (LoadedQpaAlgebra::Rational(_), "rational")
                    | (LoadedQpaAlgebra::Binary(_), "binary")
                    | (LoadedQpaAlgebra::Prime(_), "prime")
            ));
        }
    }

    #[test]
    fn finite_field_text_uses_gap_primitive_root_powers() {
        let field = Zp::new(5);
        assert_eq!(
            field.qpa_coefficient_text(&field.to_element(1)).unwrap(),
            "z^0"
        );
        assert_eq!(
            field.qpa_coefficient_text(&field.to_element(2)).unwrap(),
            "z^1"
        );
        assert_eq!(
            field.qpa_coefficient_text(&field.to_element(4)).unwrap(),
            "z^2"
        );
        assert_eq!(
            field.qpa_coefficient_text(&field.to_element(3)).unwrap(),
            "z^3"
        );
        for exponent in 0..4 {
            let encoded = format!("z^{exponent}");
            let decoded = field.parse_qpa_coefficient(&encoded).unwrap();
            assert_eq!(field.qpa_coefficient_text(&decoded).unwrap(), encoded);
        }
    }

    #[test]
    fn prime_field_logs_match_pinned_qpa_primitive_roots() {
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_prime_field_logs.schema-v1.json"
        ))
        .unwrap();
        for field_logs in fixture["logs"].as_array().unwrap() {
            let characteristic = u32::try_from(field_logs[0].as_u64().unwrap()).unwrap();
            let logarithms = field_logs[1].as_array().unwrap();
            if characteristic == 2 {
                assert_eq!(
                    Z2.qpa_coefficient_text(&Z2.one()).unwrap(),
                    format!("z^{}", logarithms[0].as_u64().unwrap())
                );
                continue;
            }
            let field = Zp::new(characteristic);
            for (offset, logarithm) in logarithms.iter().enumerate() {
                assert_eq!(
                    field
                        .qpa_coefficient_text(
                            &field.to_element(u32::try_from(offset + 1).unwrap()),
                        )
                        .unwrap(),
                    format!("z^{}", logarithm.as_u64().unwrap()),
                    "GF({characteristic}) residue {}",
                    offset + 1
                );
            }
        }
    }

    #[test]
    fn save_overwrite_branches_preserve_qpa_results_and_query_bug() {
        let algebra = path_fixture();
        let path = temporary_path();
        let _ = fs::remove_file(&path);
        assert!(save_qpa_algebra_file(&algebra, &path, "delete old", || Ok(true)).unwrap());
        assert!(!save_qpa_algebra_file(&algebra, &path, "keep existing", || Ok(true)).unwrap());
        assert!(!save_qpa_algebra_file(&algebra, &path, "query", || Ok(false)).unwrap());
        assert!(save_qpa_algebra_file(&algebra, &path, "query", || Ok(true)).unwrap());
        fs::remove_file(&path).unwrap();

        assert!(matches!(
            save_qpa_algebra_file(&algebra, &path, "query", || Ok(true)),
            Err(QpaPersistenceError::QueryNewFileUnboundOutput)
        ));
        assert!(matches!(
            save_qpa_algebra_file(&algebra, &path, "bad", || Ok(true)),
            Err(QpaPersistenceError::InvalidOverwriteOption { .. })
        ));
    }

    #[test]
    fn safe_reader_rejects_unsupported_and_malformed_input() {
        let extension_field = "IsPathAlgebra\nVertices:\nv1\nArrows:\nField:\nGF(4)\n";
        assert!(matches!(
            parse_qpa_algebra_text(extension_field),
            Err(QpaPersistenceError::UnsupportedPrimeField { characteristic: 4 })
        ));
        assert!(matches!(
            parse_qpa_algebra_text("not GAP code\n"),
            Err(QpaPersistenceError::UnknownAlgebraType { .. })
        ));
        let unknown_arrow = "IsQuotientOfPathAlgebra\nVertices:\nv1\nArrows:\nField:\nRationals\nRelations:\n(1)*x\n";
        assert!(matches!(
            parse_qpa_algebra_text(unknown_arrow),
            Err(QpaPersistenceError::UnknownArrowName { .. })
        ));
        let zero_denominator = "IsQuotientOfPathAlgebra\nVertices:\nv1\nArrows:\nField:\nRationals\nRelations:\n(1/0)*\n";
        assert!(matches!(
            parse_qpa_algebra_text(zero_denominator),
            Err(QpaPersistenceError::ZeroRationalDenominator)
        ));
    }

    #[test]
    fn reader_preserves_qpa_prefix_headers_and_empty_quotient_dispatch() {
        let text = concat!(
            "IsQuotientOfPathAlgebra\n",
            "Vertices from QPA\n",
            "v1\n",
            "Arrows from QPA\n",
            "Field from QPA\n",
            "Rationals\n",
            "Relations from QPA\n",
        );
        let LoadedQpaAlgebra::Rational(algebra) = parse_qpa_algebra_text(text).unwrap() else {
            panic!("QPA prefix-header fixture must infer the rational field");
        };
        assert!(algebra.is_path_algebra());
        assert!(algebra.relations().is_empty());
        assert!(algebra.as_quotient().is_none());
    }

    fn temporary_path() -> PathBuf {
        std::env::temp_dir().join(format!(
            "qpa-rust-persistence-test-{}.txt",
            std::process::id()
        ))
    }
}
