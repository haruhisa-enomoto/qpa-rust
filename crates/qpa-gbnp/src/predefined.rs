//! QPA's basic predefined, poset, and Brauer-configuration algebras.
//!
//! These are source-section translations of QPA 1.37
//! `lib/pathalgpredef.gi`. Basic nonempty presentations retain their original
//! relation ideal before attaching the configured GBNP basis; QPA constructors
//! that divide by a relation list retain that distinct quotient path. QPA's
//! installed constructor filters and attributes are represented explicitly
//! rather than inferred from the resulting presentation.

use std::fmt::Display;

use numerica::domains::{Field, integer::Integer};
use qpa_core::{
    ArrowSpec, NamedArrowSpec, PathAlgebra, PathAlgebraElement, PathAlgebraError, PathAlgebraIdeal,
    PathEndpoint, Poset, PosetError, Quiver, QuiverAlgebra, QuiverError, QuotientError,
    QuotientPathAlgebra, VertexIndex,
};
use thiserror::Error;

use crate::{QpaGbnpError, gbnp_groebner_basis};

/// The diagnostic printed by QPA before it returns an invalid admissible
/// sequence unchanged.
pub const INVALID_ADMISSIBLE_SEQUENCE_DIAGNOSTIC: &str =
    "The admissible sequence is not a valid one.\n";

/// The constructor filter installed on a predefined QPA algebra.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum PredefinedAlgebraKind {
    /// `IsNakayamaAlgebra`.
    Nakayama,
    /// `IsCanonicalAlgebra`.
    Canonical,
    /// `IsKroneckerAlgebra`.
    Kronecker,
}

/// A path algebra or admissible quotient carrying QPA's constructor filter.
///
/// QPA does not recognize these three classes from arbitrary presentations in
/// this source block.  It sets a filter on the object returned by the
/// constructor.  This wrapper makes that provenance observable without adding
/// an independent recognition algorithm.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PredefinedAlgebra<F: Field> {
    kind: PredefinedAlgebraKind,
    algebra: QuiverAlgebra<F>,
}

impl<F: Field> PredefinedAlgebra<F> {
    fn new(kind: PredefinedAlgebraKind, algebra: QuiverAlgebra<F>) -> Self {
        Self { kind, algebra }
    }

    /// Return the constructor filter set by QPA.
    #[must_use]
    pub const fn kind(&self) -> PredefinedAlgebraKind {
        self.kind
    }

    /// Borrow the constructed path algebra or admissible quotient.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Consume the constructor tag and return the underlying quiver algebra.
    #[must_use]
    pub fn into_algebra(self) -> QuiverAlgebra<F> {
        self.algebra
    }

    /// Return whether the `IsNakayamaAlgebra` filter is installed.
    #[must_use]
    pub const fn has_nakayama_algebra_filter(&self) -> bool {
        matches!(self.kind, PredefinedAlgebraKind::Nakayama)
    }

    /// Return QPA's constructor-known `IsNakayamaAlgebra` value.
    ///
    /// `None` means the constructor did not install the filter; QPA has no
    /// recognition method here from which to infer `false`.
    #[must_use]
    pub const fn is_nakayama_algebra(&self) -> Option<bool> {
        if self.has_nakayama_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }

    /// Return whether the `IsCanonicalAlgebra` filter is installed.
    #[must_use]
    pub const fn has_canonical_algebra_filter(&self) -> bool {
        matches!(self.kind, PredefinedAlgebraKind::Canonical)
    }

    /// Return QPA's constructor-known `IsCanonicalAlgebra` value.
    ///
    /// `None` means the constructor did not install the filter; QPA has no
    /// recognition method here from which to infer `false`.
    #[must_use]
    pub const fn is_canonical_algebra(&self) -> Option<bool> {
        if self.has_canonical_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }

    /// Return whether the `IsKroneckerAlgebra` filter is installed.
    #[must_use]
    pub const fn has_kronecker_algebra_filter(&self) -> bool {
        matches!(self.kind, PredefinedAlgebraKind::Kronecker)
    }

    /// Return QPA's constructor-known `IsKroneckerAlgebra` value.
    ///
    /// `None` means the constructor did not install the filter; QPA has no
    /// recognition method here from which to infer `false`.
    #[must_use]
    pub const fn is_kronecker_algebra(&self) -> Option<bool> {
        if self.has_kronecker_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }

    /// Return QPA's constructor-implied `IsFiniteTypeAlgebra` value.
    ///
    /// `None` means this source block installs no implication for the
    /// constructor; it does not assert mathematical falsity.
    #[must_use]
    pub const fn is_finite_type_algebra(&self) -> Option<bool> {
        if self.has_nakayama_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }

    /// Return QPA's constructor-implied finite-global-dimension value.
    ///
    /// `None` means this source block installs no implication for the
    /// constructor; it does not assert mathematical falsity.
    #[must_use]
    pub const fn is_finite_global_dimension_algebra(&self) -> Option<bool> {
        if self.has_canonical_algebra_filter() || self.has_kronecker_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }

    /// Return QPA's constructor-implied `IsHereditaryAlgebra` value.
    ///
    /// `None` means this source block installs no implication for the
    /// constructor; it does not assert mathematical falsity.
    #[must_use]
    pub const fn is_hereditary_algebra(&self) -> Option<bool> {
        if self.has_kronecker_algebra_filter() {
            Some(true)
        } else {
            None
        }
    }
}

/// QPA's two possible return types from `NakayamaAlgebra`.
///
/// QPA returns an invalid sequence unchanged after printing
/// [`INVALID_ADMISSIBLE_SEQUENCE_DIAGNOSTIC`].  That branch is not a GAP
/// error, so Rust retains it as a successful typed result.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum NakayamaAlgebraResult<F: Field> {
    /// A path algebra or admissible quotient tagged as Nakayama.
    Algebra(PredefinedAlgebra<F>),
    /// The original invalid sequence returned unchanged by QPA.
    InvalidSequence(Vec<Integer>),
}

impl<F: Field> NakayamaAlgebraResult<F> {
    /// Borrow the constructed algebra, if the sequence was admissible.
    #[must_use]
    pub const fn as_algebra(&self) -> Option<&PredefinedAlgebra<F>> {
        match self {
            Self::Algebra(algebra) => Some(algebra),
            Self::InvalidSequence(_) => None,
        }
    }

    /// Return the exact QPA diagnostic for the invalid-sequence branch.
    #[must_use]
    pub const fn diagnostic(&self) -> Option<&'static str> {
        match self {
            Self::Algebra(_) => None,
            Self::InvalidSequence(_) => Some(INVALID_ADMISSIBLE_SEQUENCE_DIAGNOSTIC),
        }
    }
}

/// Translate `NakayamaAlgebra(K, admiss_seq)` using QPA's default GBNP
/// backend.
///
/// The input uses Numerica's arbitrary-precision [`Integer`] so the intended
/// GAP integer domain is not artificially narrowed.  An invalid nonempty
/// sequence is returned in [`NakayamaAlgebraResult::InvalidSequence`], exactly
/// as in QPA.
///
/// # Errors
///
/// The empty list reproduces QPA's out-of-range list access as a typed error.
/// Valid entries that cannot be represented by the platform's allocation
/// index type are also rejected.  Quiver, path-algebra, GBNP, and quotient
/// failures are propagated.
pub fn nakayama_algebra<F: Field>(
    field: F,
    admissible_sequence: &[Integer],
) -> Result<NakayamaAlgebraResult<F>, PredefinedAlgebraError> {
    if admissible_sequence.is_empty() {
        return Err(PredefinedAlgebraError::EmptyAdmissibleSequence);
    }
    if !is_admissible_sequence(admissible_sequence) {
        return Ok(NakayamaAlgebraResult::InvalidSequence(
            admissible_sequence.to_vec(),
        ));
    }

    // For a singleton list QPA skips the loop containing its positivity test.
    // Consequently `[0]` and negative singletons reach the cyclic branch,
    // generate `One(KQ)` as their sole relation, and fail in the configured
    // Gröbner backend. Preserve that call path instead of underflowing a Rust
    // unsigned length conversion.
    if admissible_sequence.len() == 1 && admissible_sequence[0] < 1 {
        let only_vertex = vertex(1)?;
        let quiver = Quiver::from_indexed(1, vec![ArrowSpec::unnamed(only_vertex, only_vertex)])?;
        let path_algebra = PathAlgebra::new(field, quiver)?;
        let relation = path_algebra.one();
        let algebra = factor_predefined_presentation(path_algebra, &[relation])?;
        return Ok(NakayamaAlgebraResult::Algebra(PredefinedAlgebra::new(
            PredefinedAlgebraKind::Nakayama,
            algebra,
        )));
    }

    let sequence = admissible_sequence
        .iter()
        .enumerate()
        .map(|(position, value)| {
            usize::try_from(value.clone()).map_err(|_| {
                PredefinedAlgebraError::SizeNotRepresentable {
                    parameter: "admissible-sequence entry",
                    position: position + 1,
                    value: value.clone(),
                }
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let number_of_vertices = sequence.len();

    let (path_algebra, relations) = if sequence[number_of_vertices - 1] == 1 {
        linear_nakayama_presentation(field, &sequence)?
    } else {
        cyclic_nakayama_presentation(field, &sequence)?
    };
    let algebra = factor_predefined_presentation(path_algebra, &relations)?;
    Ok(NakayamaAlgebraResult::Algebra(PredefinedAlgebra::new(
        PredefinedAlgebraKind::Nakayama,
        algebra,
    )))
}

/// Translate QPA's sequence-first `NakayamaAlgebra(admiss_seq, K)` overload.
///
/// # Errors
///
/// Delegates to [`nakayama_algebra`].
pub fn nakayama_algebra_sequence_first<F: Field>(
    admissible_sequence: &[Integer],
    field: F,
) -> Result<NakayamaAlgebraResult<F>, PredefinedAlgebraError> {
    nakayama_algebra(field, admissible_sequence)
}

/// Translate `CanonicalAlgebra(K, weights, relcoeff)` using QPA's default
/// GBNP backend.
///
/// Coefficients have type `F::Element`, so QPA's membership check is enforced
/// by Rust's type system.  The number, order, and names of all vertices,
/// arrows, arm products, and relations follow the pinned source.
///
/// # Errors
///
/// Returns QPA-compatible typed errors for mismatched input lengths, weights
/// below two, and zero relation coefficients.  Structural and GBNP failures
/// are propagated.
pub fn canonical_algebra<F: Field>(
    field: F,
    weights: &[Integer],
    relation_coefficients: &[F::Element],
) -> Result<PredefinedAlgebra<F>, PredefinedAlgebraError> {
    if weights.len() != relation_coefficients.len() + 2 {
        return Err(PredefinedAlgebraError::CanonicalArmCoefficientMismatch {
            arms: weights.len(),
            coefficients: relation_coefficients.len(),
        });
    }
    if let Some((position, weight)) = weights.iter().enumerate().find(|(_, weight)| *weight < &2) {
        return Err(PredefinedAlgebraError::CanonicalWeightBelowTwo {
            position: position + 1,
            weight: weight.clone(),
        });
    }
    if let Some(position) = relation_coefficients
        .iter()
        .position(|coefficient| field.is_zero(coefficient))
    {
        return Err(PredefinedAlgebraError::ZeroCanonicalRelationCoefficient {
            position: position + 1,
        });
    }

    let weights = weights
        .iter()
        .enumerate()
        .map(|(position, value)| {
            usize::try_from(value.clone()).map_err(|_| {
                PredefinedAlgebraError::SizeNotRepresentable {
                    parameter: "canonical weight",
                    position: position + 1,
                    value: value.clone(),
                }
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut vertices = vec![Box::<str>::from("v")];
    for (arm, &weight) in weights.iter().enumerate() {
        for position in 2..=weight {
            vertices.push(format!("v{}{position}", arm + 1).into_boxed_str());
        }
    }
    vertices.push(Box::<str>::from("w"));

    let mut arrows = Vec::new();
    for (arm, &weight) in weights.iter().enumerate() {
        let arm = arm + 1;
        arrows.push(NamedArrowSpec::named(
            "v",
            format!("v{arm}2"),
            format!("a{arm}"),
        ));
        for position in 2..weight {
            arrows.push(NamedArrowSpec::named(
                format!("v{arm}{position}"),
                format!("v{arm}{}", position + 1),
                format!("a{arm}{position}"),
            ));
        }
        arrows.push(NamedArrowSpec::named(
            format!("v{arm}{weight}"),
            "w",
            format!("b{arm}"),
        ));
    }

    let quiver = Quiver::from_named(vertices, arrows)?;
    let path_algebra = PathAlgebra::new(field, quiver)?;
    if weights.len() == 2 {
        return Ok(PredefinedAlgebra::new(
            PredefinedAlgebraKind::Canonical,
            QuiverAlgebra::Path(path_algebra),
        ));
    }

    let generators = path_algebra.generators();
    let arrow_generators = &generators[path_algebra.quiver().number_of_vertices()..];
    let mut arms = Vec::with_capacity(weights.len());
    let mut start = 0_usize;
    for &weight in &weights {
        let end = start
            .checked_add(weight)
            .ok_or(PredefinedAlgebraError::GeneratedIndexOverflow)?;
        arms.push(product_of_arrows(
            &path_algebra,
            arrow_generators,
            start..end,
        )?);
        start = end;
    }

    let mut relations = Vec::with_capacity(weights.len() - 2);
    for arm in 2..weights.len() {
        let relation = arms[arm]
            .sub(&arms[1])?
            .add(&arms[0].scale_left(&relation_coefficients[arm - 2])?)?;
        relations.push(relation);
    }
    let algebra = factor_predefined_presentation(path_algebra, &relations)?;
    Ok(PredefinedAlgebra::new(
        PredefinedAlgebraKind::Canonical,
        algebra,
    ))
}

/// Translate QPA's two-argument `CanonicalAlgebra(K, weights)` overload.
///
/// # Errors
///
/// QPA accepts this overload only for exactly two weights; construction then
/// delegates to [`canonical_algebra`] with an empty coefficient list.
pub fn canonical_algebra_two_weights<F: Field>(
    field: F,
    weights: &[Integer],
) -> Result<PredefinedAlgebra<F>, PredefinedAlgebraError> {
    if weights.len() != 2 {
        return Err(PredefinedAlgebraError::CanonicalTwoArgumentWeightCount {
            weights: weights.len(),
        });
    }
    canonical_algebra(field, weights, &[])
}

/// Translate `KroneckerAlgebra(K, n)`.
///
/// # Errors
///
/// QPA rejects fewer than two arrows.  Counts outside the platform allocation
/// index range and quiver/path-algebra construction failures are also
/// reported.
pub fn kronecker_algebra<F: Field>(
    field: F,
    number_of_arrows: &Integer,
) -> Result<PredefinedAlgebra<F>, PredefinedAlgebraError> {
    if number_of_arrows < &2 {
        return Err(PredefinedAlgebraError::KroneckerNeedsTwoArrows {
            number_of_arrows: number_of_arrows.clone(),
        });
    }
    let count = usize::try_from(number_of_arrows.clone()).map_err(|_| {
        PredefinedAlgebraError::SizeNotRepresentable {
            parameter: "Kronecker arrow count",
            position: 1,
            value: number_of_arrows.clone(),
        }
    })?;
    let source = vertex(1)?;
    let target = vertex(2)?;
    let arrows = (1..=count)
        .map(|position| ArrowSpec::named(source, target, format!("a{position}")))
        .collect();
    let quiver = Quiver::from_indexed(2, arrows)?;
    let path_algebra = PathAlgebra::new(field, quiver)?;
    Ok(PredefinedAlgebra::new(
        PredefinedAlgebraKind::Kronecker,
        QuiverAlgebra::Path(path_algebra),
    ))
}

/// A quiver algebra constructed by QPA's `PosetAlgebra`, together with its
/// identity-preserving `PosetOfPosetAlgebra` attribute.
#[derive(Clone, Debug)]
pub struct PosetAlgebra<F: Field, Element> {
    algebra: QuiverAlgebra<F>,
    poset: Poset<Element>,
}

/// A named vertex and positive multiplicity in QPA's Brauer-configuration
/// input.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BrauerConfigurationVertex {
    /// QPA vertex name.
    pub name: String,
    /// QPA vertex multiplicity.
    pub multiplicity: Integer,
}

/// A named Brauer polygon followed by its vertex occurrences.
///
/// Occurrences are retained rather than deduplicated because QPA's type-two
/// relation loop treats repeated vertices in a polygon separately.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BrauerConfigurationPolygon {
    /// Polygon name, later used as a quiver vertex name.
    pub name: String,
    /// Ordered vertex occurrences in the polygon.
    pub vertices: Vec<String>,
}

/// Typed form of QPA's `[vertices, polygons, orientations]` input.
///
/// [`brauer_configuration_algebra`] sorts orientation rows in place at the
/// same point as the pinned QPA method, including before returning a final
/// containment-validation failure.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BrauerConfiguration {
    /// Named vertices with multiplicities.
    pub vertices: Vec<BrauerConfigurationVertex>,
    /// Named polygons with vertex occurrences.
    pub polygons: Vec<BrauerConfigurationPolygon>,
    /// One polygon-name orientation row for each vertex.
    pub orientations: Vec<Vec<String>>,
}

/// A path algebra or relation-list quotient returned by QPA's
/// `BrauerConfigurationAlgebra` and tagged symmetric by the constructor.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BrauerConfigurationAlgebra<F: Field> {
    algebra: QuiverAlgebra<F>,
}

impl<F: Field> BrauerConfigurationAlgebra<F> {
    /// Borrow the constructed path algebra or quotient.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Return the `IsSymmetricAlgebra` value set by the constructor.
    #[must_use]
    pub const fn is_symmetric_algebra(&self) -> Option<bool> {
        Some(true)
    }

    /// Consume the symmetric tag and return the underlying quiver algebra.
    #[must_use]
    pub fn into_algebra(self) -> QuiverAlgebra<F> {
        self.algebra
    }
}

/// QPA's algebra-or-`fail` result from `BrauerConfigurationAlgebra`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BrauerConfigurationAlgebraResult<F: Field> {
    /// A valid configuration's symmetric quiver algebra.
    Algebra(BrauerConfigurationAlgebra<F>),
    /// QPA's `fail`, paired with the validation branch that printed it.
    InvalidConfiguration(BrauerConfigurationValidationError),
}

impl<F: Field> BrauerConfigurationAlgebraResult<F> {
    /// Borrow the constructed algebra, if validation succeeded.
    #[must_use]
    pub const fn as_algebra(&self) -> Option<&BrauerConfigurationAlgebra<F>> {
        match self {
            Self::Algebra(algebra) => Some(algebra),
            Self::InvalidConfiguration(_) => None,
        }
    }

    /// Return the exact diagnostic printed by QPA before `fail`.
    #[must_use]
    pub const fn diagnostic(&self) -> Option<&'static str> {
        match self {
            Self::Algebra(_) => None,
            Self::InvalidConfiguration(error) => Some(error.diagnostic()),
        }
    }
}

/// A representable validation failure in QPA's Brauer-configuration method.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BrauerConfigurationValidationError {
    /// The vertex list is empty.
    NoVertices,
    /// A vertex multiplicity is less than one.
    NonPositiveMultiplicity {
        /// One-based vertex position.
        position: usize,
    },
    /// Vertex names are not duplicate-free.
    DuplicateVertexName {
        /// Repeated name.
        name: String,
    },
    /// A polygon has no vertex occurrence after its name.
    PolygonWithoutVertex {
        /// One-based polygon position.
        position: usize,
    },
    /// A polygon name is also a Brauer vertex name.
    PolygonNameMatchesVertex {
        /// Conflicting name.
        name: String,
    },
    /// A polygon contains a name absent from the vertex list.
    PolygonContainsUnknownVertex {
        /// Polygon name.
        polygon: String,
        /// Unknown vertex name.
        vertex: String,
    },
    /// Polygon names are not duplicate-free.
    DuplicatePolygonName {
        /// Repeated name.
        name: String,
    },
    /// The number of orientations differs from the number of vertices.
    OrientationCountMismatch {
        /// Required number of rows.
        expected: usize,
        /// Supplied number of rows.
        actual: usize,
    },
    /// An orientation row is empty.
    EmptyOrientation {
        /// One-based vertex/orientation position.
        position: usize,
    },
    /// An orientation names an unknown polygon.
    OrientationContainsUnknownPolygon {
        /// One-based vertex/orientation position.
        position: usize,
        /// Unknown polygon name.
        polygon: String,
    },
    /// A sorted orientation is not exactly the polygon-containment list.
    OrientationContainmentMismatch {
        /// One-based vertex/orientation position.
        position: usize,
    },
}

impl BrauerConfigurationValidationError {
    /// Return QPA's exact printed diagnostic for this `fail` branch.
    #[must_use]
    pub const fn diagnostic(&self) -> &'static str {
        match self {
            Self::NoVertices => "There must be at least 1 vertex\n",
            Self::NonPositiveMultiplicity { .. } => "Multiplicities must be positive integers.\n",
            Self::DuplicateVertexName { .. } => "Each vertex must have a unique name.\n",
            Self::PolygonWithoutVertex { .. } => {
                "Edges or Polygons should be lists of length 2 or greater.\n"
            }
            Self::PolygonNameMatchesVertex { .. } => {
                "An edge or polygon cannot have the same name as a vertex.\n"
            }
            Self::PolygonContainsUnknownVertex { .. } => {
                "An edge or polygon cannot contain a vertex which was not listed.\n"
            }
            Self::DuplicatePolygonName { .. } => "Edges or polygons must have unique names.\n",
            Self::OrientationCountMismatch { .. } => {
                "There must be an orientation corresponding to each vertex.\n"
            }
            Self::EmptyOrientation { .. } => "Orientations must be lists of length 1 or greater.\n",
            Self::OrientationContainsUnknownPolygon { .. } => {
                "Orientations may only contain edges or polygons which were listed.\n"
            }
            Self::OrientationContainmentMismatch { .. } => {
                "Orientations must contain exactly the edges or polygons which contain the corresponding vertex.\n"
            }
        }
    }
}

/// Translate `BrauerConfigurationAlgebra(K, configuration)` through QPA's
/// three relation families and default GBNP backend.
///
/// Rust's typed input enforces the source's outer-list, row-list, string, and
/// integer shape checks.  Every remaining validation branch returns
/// [`BrauerConfigurationAlgebraResult::InvalidConfiguration`], corresponding
/// to QPA's printed diagnostic followed by `fail`.  QPA sorts orientation rows
/// in place during its final validation loop; this function deliberately does
/// the same.
///
/// # Errors
///
/// Propagates named-quiver construction, path arithmetic, generated lookup,
/// QPA-to-GBNP conversion, and quotient construction failures.
pub fn brauer_configuration_algebra<F: Field>(
    field: F,
    configuration: &mut BrauerConfiguration,
) -> Result<BrauerConfigurationAlgebraResult<F>, PredefinedAlgebraError> {
    if let Err(error) = validate_brauer_configuration(configuration) {
        return Ok(BrauerConfigurationAlgebraResult::InvalidConfiguration(
            error,
        ));
    }

    let BrauerConfigurationPresentation {
        path_algebra,
        relation_families:
            [
                mut type_one_relations,
                type_two_relations,
                type_three_relations,
            ],
    } = brauer_configuration_algebra_presentation(field, configuration)?;
    type_one_relations.extend(type_two_relations);
    type_one_relations.extend(type_three_relations);
    let algebra = factor_relation_list_presentation(path_algebra, &type_one_relations)?;
    Ok(BrauerConfigurationAlgebraResult::Algebra(
        BrauerConfigurationAlgebra { algebra },
    ))
}

fn validate_brauer_configuration(
    configuration: &mut BrauerConfiguration,
) -> Result<(), BrauerConfigurationValidationError> {
    validate_brauer_vertices(&configuration.vertices)?;
    let vertex_names = configuration
        .vertices
        .iter()
        .map(|vertex| vertex.name.clone())
        .collect::<Vec<_>>();
    validate_brauer_polygons(&configuration.polygons, &vertex_names)?;
    let polygon_names = configuration
        .polygons
        .iter()
        .map(|polygon| polygon.name.clone())
        .collect::<Vec<_>>();
    validate_brauer_orientations(
        &configuration.vertices,
        &configuration.polygons,
        &polygon_names,
        &mut configuration.orientations,
    )
}

fn validate_brauer_vertices(
    vertices: &[BrauerConfigurationVertex],
) -> Result<(), BrauerConfigurationValidationError> {
    if vertices.is_empty() {
        return Err(BrauerConfigurationValidationError::NoVertices);
    }
    if let Some((position, _)) = vertices
        .iter()
        .enumerate()
        .find(|(_, vertex)| vertex.multiplicity < 1)
    {
        return Err(
            BrauerConfigurationValidationError::NonPositiveMultiplicity {
                position: position + 1,
            },
        );
    }
    if let Some(name) = first_duplicate(vertices.iter().map(|vertex| &vertex.name)) {
        return Err(BrauerConfigurationValidationError::DuplicateVertexName { name: name.clone() });
    }
    Ok(())
}

fn validate_brauer_polygons(
    polygons: &[BrauerConfigurationPolygon],
    vertex_names: &[String],
) -> Result<(), BrauerConfigurationValidationError> {
    for (position, polygon) in polygons.iter().enumerate() {
        if polygon.vertices.is_empty() {
            return Err(BrauerConfigurationValidationError::PolygonWithoutVertex {
                position: position + 1,
            });
        }
        if vertex_names.contains(&polygon.name) {
            return Err(
                BrauerConfigurationValidationError::PolygonNameMatchesVertex {
                    name: polygon.name.clone(),
                },
            );
        }
        if let Some(vertex) = polygon
            .vertices
            .iter()
            .find(|vertex| !vertex_names.contains(vertex))
        {
            return Err(
                BrauerConfigurationValidationError::PolygonContainsUnknownVertex {
                    polygon: polygon.name.clone(),
                    vertex: vertex.to_string(),
                },
            );
        }
    }
    if let Some(name) = first_duplicate(polygons.iter().map(|polygon| &polygon.name)) {
        return Err(BrauerConfigurationValidationError::DuplicatePolygonName {
            name: name.clone(),
        });
    }
    Ok(())
}

fn validate_brauer_orientations(
    vertices: &[BrauerConfigurationVertex],
    polygons: &[BrauerConfigurationPolygon],
    polygon_names: &[String],
    orientations: &mut [Vec<String>],
) -> Result<(), BrauerConfigurationValidationError> {
    if orientations.len() != vertices.len() {
        return Err(
            BrauerConfigurationValidationError::OrientationCountMismatch {
                expected: vertices.len(),
                actual: orientations.len(),
            },
        );
    }
    for (position, orientation) in orientations.iter().enumerate() {
        if orientation.is_empty() {
            return Err(BrauerConfigurationValidationError::EmptyOrientation {
                position: position + 1,
            });
        }
        if let Some(polygon) = orientation
            .iter()
            .find(|polygon| !polygon_names.contains(polygon))
        {
            return Err(
                BrauerConfigurationValidationError::OrientationContainsUnknownPolygon {
                    position: position + 1,
                    polygon: polygon.to_string(),
                },
            );
        }
    }
    for (position, (orientation, vertex)) in orientations.iter_mut().zip(vertices).enumerate() {
        let mut containing = polygons
            .iter()
            .filter(|polygon| polygon.vertices.contains(&vertex.name))
            .map(|polygon| polygon.name.clone())
            .collect::<Vec<_>>();
        orientation.sort();
        containing.sort();
        if orientation != &containing {
            return Err(
                BrauerConfigurationValidationError::OrientationContainmentMismatch {
                    position: position + 1,
                },
            );
        }
    }
    Ok(())
}

struct BrauerConfigurationPresentation<F: Field> {
    path_algebra: PathAlgebra<F>,
    relation_families: [Vec<PathAlgebraElement<F>>; 3],
}

fn brauer_configuration_algebra_presentation<F: Field>(
    field: F,
    configuration: &BrauerConfiguration,
) -> Result<BrauerConfigurationPresentation<F>, PredefinedAlgebraError> {
    let (quiver, special_cycles) = brauer_quiver_and_special_cycles(configuration)?;
    let path_algebra = PathAlgebra::new(field, quiver)?;
    let generators = path_algebra.generators();
    let arrow_elements = &generators[path_algebra.quiver().number_of_vertices()..];
    let type_one_relations = brauer_type_one_relations(
        arrow_elements,
        &special_cycles,
        configuration.vertices.len(),
    )?;
    let type_two_relations = brauer_type_two_relations(
        &path_algebra,
        arrow_elements,
        &special_cycles,
        configuration,
    )?;
    let type_three_relations = brauer_type_three_relations(
        &path_algebra,
        arrow_elements,
        &special_cycles,
        &configuration.vertices,
    )?;
    Ok(BrauerConfigurationPresentation {
        path_algebra,
        relation_families: [type_one_relations, type_two_relations, type_three_relations],
    })
}

fn brauer_quiver_and_special_cycles(
    configuration: &BrauerConfiguration,
) -> Result<(Quiver, Vec<Vec<usize>>), PredefinedAlgebraError> {
    let mut arrows = Vec::new();
    let mut special_cycles = Vec::with_capacity(configuration.vertices.len());
    for (vertex_index, vertex) in configuration.vertices.iter().enumerate() {
        let orientation = &configuration.orientations[vertex_index];
        if orientation.len() > 1 || vertex.multiplicity > 1 {
            let mut cycle = Vec::with_capacity(orientation.len());
            for position in 0..orientation.len() {
                let next = (position + 1) % orientation.len();
                arrows.push(NamedArrowSpec::named(
                    orientation[position].clone(),
                    orientation[next].clone(),
                    format!(
                        "{}_{}_from_{}_to_{}",
                        vertex.name,
                        position + 1,
                        orientation[position],
                        orientation[next]
                    ),
                ));
                cycle.push(arrows.len() - 1);
            }
            special_cycles.push(cycle);
        } else {
            special_cycles.push(Vec::new());
        }
    }
    let quiver_vertices = configuration
        .polygons
        .iter()
        .map(|polygon| polygon.name.clone().into_boxed_str())
        .collect();
    let quiver = Quiver::from_named(quiver_vertices, arrows)?;
    Ok((quiver, special_cycles))
}

fn brauer_type_one_relations<F: Field>(
    arrow_elements: &[PathAlgebraElement<F>],
    special_cycles: &[Vec<usize>],
    number_of_vertices: usize,
) -> Result<Vec<PathAlgebraElement<F>>, PredefinedAlgebraError> {
    let mut type_one_relations = Vec::new();
    for first in 0..number_of_vertices.saturating_sub(1) {
        for second in first + 1..number_of_vertices {
            for &left in &special_cycles[first] {
                for &right in &special_cycles[second] {
                    type_one_relations.push(
                        generated_arrow(arrow_elements, left)?
                            .multiply(generated_arrow(arrow_elements, right)?)?,
                    );
                    type_one_relations.push(
                        generated_arrow(arrow_elements, right)?
                            .multiply(generated_arrow(arrow_elements, left)?)?,
                    );
                }
            }
        }
    }
    Ok(type_one_relations)
}

fn brauer_type_two_relations<F: Field>(
    path_algebra: &PathAlgebra<F>,
    arrow_elements: &[PathAlgebraElement<F>],
    special_cycles: &[Vec<usize>],
    configuration: &BrauerConfiguration,
) -> Result<Vec<PathAlgebraElement<F>>, PredefinedAlgebraError> {
    let vertex_names = configuration
        .vertices
        .iter()
        .map(|vertex| vertex.name.clone())
        .collect::<Vec<_>>();
    let mut type_two_relations = Vec::new();
    for (polygon_index, polygon) in configuration.polygons.iter().enumerate() {
        let polygon_vertex = path_algebra
            .quiver()
            .vertex(
                VertexIndex::new(polygon_index + 1)
                    .ok_or(PredefinedAlgebraError::GeneratedIndexOverflow)?,
            )
            .ok_or(PredefinedAlgebraError::GeneratedBrauerLookupMissing)?;
        if polygon_vertex.in_degree() <= 1 {
            continue;
        }
        for first_occurrence in 0..polygon.vertices.len().saturating_sub(1) {
            let first_vertex = vertex_names
                .iter()
                .position(|name| name == &polygon.vertices[first_occurrence])
                .ok_or(PredefinedAlgebraError::GeneratedBrauerLookupMissing)?;
            if special_cycles[first_vertex].is_empty() {
                continue;
            }
            for second_occurrence in first_occurrence + 1..polygon.vertices.len() {
                let second_vertex = vertex_names
                    .iter()
                    .position(|name| name == &polygon.vertices[second_occurrence])
                    .ok_or(PredefinedAlgebraError::GeneratedBrauerLookupMissing)?;
                if special_cycles[second_vertex].is_empty() {
                    continue;
                }
                let first_start = configuration.orientations[first_vertex]
                    .iter()
                    .position(|name| name == &polygon.name)
                    .ok_or(PredefinedAlgebraError::GeneratedBrauerLookupMissing)?;
                let second_start = configuration.orientations[second_vertex]
                    .iter()
                    .position(|name| name == &polygon.name)
                    .ok_or(PredefinedAlgebraError::GeneratedBrauerLookupMissing)?;
                let first_cycle = brauer_special_cycle(
                    path_algebra,
                    arrow_elements,
                    &special_cycles[first_vertex],
                    first_start,
                )?;
                let second_cycle = brauer_special_cycle(
                    path_algebra,
                    arrow_elements,
                    &special_cycles[second_vertex],
                    second_start,
                )?;
                type_two_relations.push(
                    pow_path_element_integer(
                        &first_cycle,
                        &configuration.vertices[first_vertex].multiplicity,
                    )?
                    .sub(&pow_path_element_integer(
                        &second_cycle,
                        &configuration.vertices[second_vertex].multiplicity,
                    )?)?,
                );
            }
        }
    }
    Ok(type_two_relations)
}

fn brauer_type_three_relations<F: Field>(
    path_algebra: &PathAlgebra<F>,
    arrow_elements: &[PathAlgebraElement<F>],
    special_cycles: &[Vec<usize>],
    vertices: &[BrauerConfigurationVertex],
) -> Result<Vec<PathAlgebraElement<F>>, PredefinedAlgebraError> {
    let mut type_three_relations = Vec::new();
    for (vertex_index, cycle) in special_cycles.iter().enumerate() {
        for start in 0..cycle.len() {
            let special_cycle = brauer_special_cycle(path_algebra, arrow_elements, cycle, start)?;
            type_three_relations.push(
                pow_path_element_integer(&special_cycle, &vertices[vertex_index].multiplicity)?
                    .multiply(generated_arrow(arrow_elements, cycle[start])?)?,
            );
        }
    }
    Ok(type_three_relations)
}

impl<F: Field, Element> PosetAlgebra<F, Element> {
    /// Borrow the constructed path algebra or relation-list quotient.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Return QPA's identity-preserving `PosetOfPosetAlgebra` attribute.
    #[must_use]
    pub const fn poset_of_poset_algebra(&self) -> &Poset<Element> {
        &self.poset
    }

    /// Return QPA's constructor-known `IsPosetAlgebra` value.
    #[must_use]
    pub const fn is_poset_algebra(&self) -> Option<bool> {
        Some(true)
    }

    /// Return the finite-global-dimension filter set by `PosetAlgebra`.
    #[must_use]
    pub const fn is_finite_global_dimension_algebra(&self) -> Option<bool> {
        Some(true)
    }

    /// Consume the attribute wrapper and return the underlying quiver algebra.
    #[must_use]
    pub fn into_algebra(self) -> QuiverAlgebra<F> {
        self.algebra
    }
}

/// Translate `PosetAlgebra(K, P)` using Rust display strings for QPA's
/// `String(p)` calls.
///
/// The quiver is the Hasse quiver.  Relations follow QPA's computation of a
/// basis of `J^2`, its endpoint subspaces, and the repeated `B[1] - B[k]`
/// loop before the relation list is passed to the default GBNP backend.
/// Cloning the supplied [`Poset`] handle preserves the identity observed by
/// `PosetOfPosetAlgebra`.
///
/// Use [`poset_algebra_with_string_function`] when an element's Rust
/// [`Display`] output is not its intended GAP `String` representation.
///
/// # Errors
///
/// Propagates poset queries, named-quiver construction, path arithmetic,
/// QPA-to-GBNP conversion, and quotient construction failures.
pub fn poset_algebra<F, Element>(
    field: F,
    poset: &Poset<Element>,
) -> Result<PosetAlgebra<F, Element>, PredefinedAlgebraError>
where
    F: Field,
    Element: Clone + Display + Ord,
{
    poset_algebra_with_string_function(field, poset, ToString::to_string)
}

/// Translate `PosetAlgebra(K, P)` with an explicit implementation of QPA's
/// `String(p)` conversion.
///
/// This is the same source algorithm as [`poset_algebra`]; only the boundary
/// used to name vertices and arrows is supplied by the caller.
///
/// # Errors
///
/// Propagates poset queries, named-quiver construction, path arithmetic,
/// QPA-to-GBNP conversion, and quotient construction failures.
pub fn poset_algebra_with_string_function<F, Element, StringFunction>(
    field: F,
    poset: &Poset<Element>,
    mut string_of_element: StringFunction,
) -> Result<PosetAlgebra<F, Element>, PredefinedAlgebraError>
where
    F: Field,
    Element: Clone + Ord,
    StringFunction: FnMut(&Element) -> String,
{
    let (path_algebra, relations) =
        poset_algebra_presentation(poset, &mut string_of_element, field)?;
    let algebra = factor_relation_list_presentation(path_algebra, &relations)?;
    Ok(PosetAlgebra {
        algebra,
        poset: poset.clone(),
    })
}

fn poset_algebra_presentation<F, Element, StringFunction>(
    poset: &Poset<Element>,
    string_of_element: &mut StringFunction,
    field: F,
) -> Result<(PathAlgebra<F>, Vec<PathAlgebraElement<F>>), PredefinedAlgebraError>
where
    F: Field,
    Element: Clone + Ord,
    StringFunction: FnMut(&Element) -> String,
{
    let elements = poset.underlying_set();
    let names = elements.iter().map(string_of_element).collect::<Vec<_>>();
    let vertices = names
        .iter()
        .map(|name| name.clone().into_boxed_str())
        .collect();
    let mut arrows = Vec::new();
    for source in 0..elements.len() {
        for target in 0..elements.len() {
            if poset
                .minimal_relations()
                .iter()
                .any(|relation| relation == &(elements[source].clone(), elements[target].clone()))
            {
                arrows.push(NamedArrowSpec::named(
                    names[source].clone(),
                    names[target].clone(),
                    format!("{}_{}", names[source], names[target]),
                ));
            }
        }
    }
    let quiver = Quiver::from_named(vertices, arrows)?;
    let path_algebra = PathAlgebra::new(field, quiver)?;

    // `Basis(ProductSpace(I, I))` for the arrow ideal I consists of every
    // path of length at least two in canonical path order.
    let radical_square_basis = path_algebra
        .quiver()
        .paths()
        .filter(|path| path.length().is_some_and(|length| length >= 2))
        .map(|path| {
            let PathEndpoint::Vertex(source) = path.source() else {
                return Err(PredefinedAlgebraError::GeneratedPathEndpointMissing);
            };
            let PathEndpoint::Vertex(target) = path.target() else {
                return Err(PredefinedAlgebraError::GeneratedPathEndpointMissing);
            };
            Ok((
                source.index().get() - 1,
                target.index().get() - 1,
                path_algebra.element(path)?,
            ))
        })
        .collect::<Result<Vec<_>, PredefinedAlgebraError>>()?;

    let mut endpoint_spaces = Vec::new();
    for source in 0..elements.len() {
        for target in 0..elements.len() {
            if source != target && poset.less_than_or_equal(&elements[source], &elements[target])? {
                let basis = radical_square_basis
                    .iter()
                    .filter(|(left, right, _)| *left == source && *right == target)
                    .map(|(_, _, element)| element.clone())
                    .collect::<Vec<_>>();
                if basis.len() > 1 {
                    endpoint_spaces.push((source, target, basis));
                }
            }
        }
    }

    let mut relations = Vec::new();
    for (source, _target, basis) in endpoint_spaces {
        // The pinned method shadows its saved target with this loop variable.
        // Hence each basis-difference list is appended n-1 times.
        for repeated_target in 0..elements.len() {
            if source != repeated_target && basis.len() > 1 {
                for other in &basis[1..] {
                    relations.push(basis[0].sub(other)?);
                }
            }
        }
    }
    Ok((path_algebra, relations))
}

fn is_admissible_sequence(sequence: &[Integer]) -> bool {
    let one = Integer::from(1);
    for pair in sequence.windows(2) {
        let previous_minus_one = &pair[0] - &one;
        if pair[1] < previous_minus_one || previous_minus_one < one {
            return false;
        }
    }
    sequence[0] >= &sequence[sequence.len() - 1] - &one
}

fn linear_nakayama_presentation<F: Field>(
    field: F,
    sequence: &[usize],
) -> Result<(PathAlgebra<F>, Vec<PathAlgebraElement<F>>), PredefinedAlgebraError> {
    let number_of_vertices = sequence.len();
    let mut arrows = Vec::with_capacity(number_of_vertices.saturating_sub(1));
    for position in 1..number_of_vertices {
        arrows.push(ArrowSpec::unnamed(vertex(position)?, vertex(position + 1)?));
    }
    let quiver = Quiver::from_indexed(number_of_vertices, arrows)?;
    let path_algebra = PathAlgebra::new(field, quiver)?;
    let generators = path_algebra.generators();
    let arrow_generators = &generators[number_of_vertices..];
    let mut relations = Vec::new();
    for position in 0..number_of_vertices.saturating_sub(1) {
        if sequence[position] < number_of_vertices - position
            && sequence[position] - 1 != sequence[position + 1]
        {
            let end = position
                .checked_add(sequence[position])
                .ok_or(PredefinedAlgebraError::GeneratedIndexOverflow)?;
            relations.push(product_of_arrows(
                &path_algebra,
                arrow_generators,
                position..end,
            )?);
        }
    }
    Ok((path_algebra, relations))
}

fn cyclic_nakayama_presentation<F: Field>(
    field: F,
    sequence: &[usize],
) -> Result<(PathAlgebra<F>, Vec<PathAlgebraElement<F>>), PredefinedAlgebraError> {
    let number_of_vertices = sequence.len();
    let mut arrows = Vec::with_capacity(number_of_vertices);
    for position in 1..number_of_vertices {
        arrows.push(ArrowSpec::unnamed(vertex(position)?, vertex(position + 1)?));
    }
    arrows.push(ArrowSpec::unnamed(vertex(number_of_vertices)?, vertex(1)?));
    let quiver = Quiver::from_indexed(number_of_vertices, arrows)?;
    let path_algebra = PathAlgebra::new(field, quiver)?;
    let generators = path_algebra.generators();
    let arrow_generators = &generators[number_of_vertices..];
    let mut relations = Vec::new();

    for position in 0..number_of_vertices {
        let next = (position + 1) % number_of_vertices;
        if sequence[position] - 1 == sequence[next] {
            continue;
        }
        let quotient = sequence[position] / number_of_vertices;
        let remainder = sequence[position] % number_of_vertices;
        let cycle = product_of_arrows(
            &path_algebra,
            arrow_generators,
            (position..number_of_vertices).chain(0..position),
        )?;
        let mut relation = path_algebra.one();
        for _ in 0..quotient {
            relation = relation.multiply(&cycle)?;
        }
        relation = append_cyclic_remainder(relation, arrow_generators, position, remainder)?;
        relations.push(relation);
    }
    Ok((path_algebra, relations))
}

fn append_cyclic_remainder<F: Field>(
    mut product: PathAlgebraElement<F>,
    arrows: &[PathAlgebraElement<F>],
    start: usize,
    length: usize,
) -> Result<PathAlgebraElement<F>, PredefinedAlgebraError> {
    for offset in 0..length {
        let index = (start + offset) % arrows.len();
        let factor = arrows
            .get(index)
            .ok_or(PredefinedAlgebraError::GeneratedArrowIndexOutOfRange { index })?;
        product = product.multiply(factor)?;
    }
    Ok(product)
}

fn product_of_arrows<F: Field>(
    path_algebra: &PathAlgebra<F>,
    arrows: &[PathAlgebraElement<F>],
    indices: impl IntoIterator<Item = usize>,
) -> Result<PathAlgebraElement<F>, PredefinedAlgebraError> {
    let mut product = path_algebra.one();
    for index in indices {
        let factor = arrows
            .get(index)
            .ok_or(PredefinedAlgebraError::GeneratedArrowIndexOutOfRange { index })?;
        product = product.multiply(factor)?;
    }
    Ok(product)
}

fn first_duplicate<'a>(values: impl IntoIterator<Item = &'a String>) -> Option<&'a String> {
    let mut seen = Vec::new();
    for value in values {
        if seen.contains(&value) {
            return Some(value);
        }
        seen.push(value);
    }
    None
}

fn generated_arrow<F: Field>(
    arrows: &[PathAlgebraElement<F>],
    index: usize,
) -> Result<&PathAlgebraElement<F>, PredefinedAlgebraError> {
    arrows
        .get(index)
        .ok_or(PredefinedAlgebraError::GeneratedArrowIndexOutOfRange { index })
}

fn brauer_special_cycle<F: Field>(
    path_algebra: &PathAlgebra<F>,
    arrows: &[PathAlgebraElement<F>],
    cycle: &[usize],
    start: usize,
) -> Result<PathAlgebraElement<F>, PredefinedAlgebraError> {
    if start >= cycle.len() {
        return Err(PredefinedAlgebraError::GeneratedBrauerLookupMissing);
    }
    product_of_arrows(
        path_algebra,
        arrows,
        cycle[start..].iter().chain(&cycle[..start]).copied(),
    )
}

fn pow_path_element_integer<F: Field>(
    element: &PathAlgebraElement<F>,
    exponent: &Integer,
) -> Result<PathAlgebraElement<F>, PredefinedAlgebraError> {
    let mut exponent = exponent.clone();
    let mut result = element.algebra().one();
    let mut base = element.clone();
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

fn factor_relation_list_presentation<F: Field>(
    path_algebra: PathAlgebra<F>,
    relations: &[PathAlgebraElement<F>],
) -> Result<QuiverAlgebra<F>, PredefinedAlgebraError> {
    if relations.is_empty() {
        return Ok(QuiverAlgebra::Path(path_algebra));
    }
    let complete_relations = gbnp_groebner_basis(relations, &path_algebra)?;
    Ok(QuiverAlgebra::Quotient(
        QuotientPathAlgebra::by_complete_groebner_relations(path_algebra, &complete_relations)?,
    ))
}

fn factor_predefined_presentation<F: Field>(
    path_algebra: PathAlgebra<F>,
    relations: &[PathAlgebraElement<F>],
) -> Result<QuiverAlgebra<F>, PredefinedAlgebraError> {
    if relations.is_empty() {
        return Ok(QuiverAlgebra::Path(path_algebra));
    }
    let complete_relations = gbnp_groebner_basis(relations, &path_algebra)?;
    let ideal = PathAlgebraIdeal::new(path_algebra, relations.to_vec())?;
    ideal.with_complete_groebner_basis(&complete_relations)?;
    let quotient = QuotientPathAlgebra::by_ideal(&ideal)?;
    let _ = quotient.is_admissible()?;
    Ok(QuiverAlgebra::Quotient(quotient))
}

fn vertex(index: usize) -> Result<VertexIndex, PredefinedAlgebraError> {
    VertexIndex::new(index).ok_or(PredefinedAlgebraError::GeneratedIndexOverflow)
}

/// Failure in QPA's predefined algebra constructors.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PredefinedAlgebraError {
    /// QPA indexes both the first and last entry of an empty list.
    #[error("the admissible sequence must not be empty")]
    EmptyAdmissibleSequence,
    /// Canonical algebras have two more arms than relation coefficients.
    #[error(
        "number of canonical arms ({arms}) does not match relation coefficients ({coefficients})"
    )]
    CanonicalArmCoefficientMismatch {
        /// Number of supplied weights.
        arms: usize,
        /// Number of supplied relation coefficients.
        coefficients: usize,
    },
    /// Every canonical weight must be at least two.
    #[error("canonical weight {position} is below two: {weight}")]
    CanonicalWeightBelowTwo {
        /// One-based list position.
        position: usize,
        /// Rejected weight.
        weight: Integer,
    },
    /// Every canonical relation coefficient must be nonzero.
    #[error("canonical relation coefficient {position} is zero")]
    ZeroCanonicalRelationCoefficient {
        /// One-based list position.
        position: usize,
    },
    /// QPA's two-argument canonical overload accepts exactly two weights.
    #[error("the two-argument canonical constructor needs two weights, received {weights}")]
    CanonicalTwoArgumentWeightCount {
        /// Number of supplied weights.
        weights: usize,
    },
    /// A Kronecker quiver must have at least two parallel arrows.
    #[error("the Kronecker arrow count must be at least two: {number_of_arrows}")]
    KroneckerNeedsTwoArrows {
        /// Rejected arrow count.
        number_of_arrows: Integer,
    },
    /// An arbitrary-precision source count cannot index a Rust allocation.
    #[error("{parameter} at position {position} cannot be represented as a platform size: {value}")]
    SizeNotRepresentable {
        /// Kind of count being converted.
        parameter: &'static str,
        /// One-based position, or one for a scalar count.
        position: usize,
        /// Rejected arbitrary-precision value.
        value: Integer,
    },
    /// A generated endpoint or range overflowed the platform index type.
    #[error("a generated predefined-algebra index overflowed")]
    GeneratedIndexOverflow,
    /// A translated source range did not address a generated arrow.
    #[error("generated arrow index {index} is outside the predefined quiver")]
    GeneratedArrowIndexOutOfRange {
        /// Zero-based internal arrow index.
        index: usize,
    },
    /// A generated nonzero path did not expose vertex endpoints.
    #[error("a generated poset-algebra path has no vertex endpoint")]
    GeneratedPathEndpointMissing,
    /// A validated Brauer configuration did not resolve a generated object.
    #[error("a generated Brauer-configuration lookup unexpectedly failed")]
    GeneratedBrauerLookupMissing,
    /// A source-translated poset query failed.
    #[error(transparent)]
    Poset(#[from] PosetError),
    /// Quiver construction failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// Path-algebra construction or arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// QPA-to-GBNP conversion or `SGrobner` failed.
    #[error(transparent)]
    Gbnp(#[from] QpaGbnpError),
    /// Ideal or quotient construction failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring,
        finite_field::Z2,
        integer::Integer,
        rational::{Q, RationalField},
    };
    use qpa_core::{AlgebraDimension, PathAlgebraElement, Poset, QuiverAlgebra};
    use serde_json::{Value, json};

    use super::{
        BrauerConfiguration, BrauerConfigurationAlgebra, BrauerConfigurationAlgebraResult,
        BrauerConfigurationPolygon, BrauerConfigurationValidationError, BrauerConfigurationVertex,
        INVALID_ADMISSIBLE_SEQUENCE_DIAGNOSTIC, NakayamaAlgebraResult, PredefinedAlgebra,
        PredefinedAlgebraError, brauer_configuration_algebra,
        brauer_configuration_algebra_presentation, canonical_algebra,
        canonical_algebra_two_weights, kronecker_algebra, nakayama_algebra,
        nakayama_algebra_sequence_first, poset_algebra, poset_algebra_presentation,
        poset_algebra_with_string_function,
    };

    fn integers(values: &[i64]) -> Vec<Integer> {
        values.iter().copied().map(Integer::from).collect()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the predefined-algebra fixture uses only small integers")
    }

    fn rational_coefficient(value: &<RationalField as numerica::domains::Set>::Element) -> Value {
        json!([
            small_integer(value.numerator_ref()),
            small_integer(value.denominator_ref())
        ])
    }

    fn element_json<F: Field>(
        element: &PathAlgebraElement<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().external_representation(),
                        "coefficient": coefficient_json(term.coefficient())
                    })
                })
                .collect(),
        )
    }

    fn algebra_json<F: Field>(
        algebra: &PredefinedAlgebra<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        let mut value = quiver_algebra_json(algebra.algebra(), coefficient_json);
        let admissible_quotient = match algebra.algebra() {
            QuiverAlgebra::Path(_) => false,
            QuiverAlgebra::Quotient(quotient) => quotient.is_admissible().unwrap(),
        };
        value
            .as_object_mut()
            .expect("algebra fixture is a JSON object")
            .insert(
                "admissible_quotient".to_owned(),
                Value::Bool(admissible_quotient),
            );
        value
    }

    fn quiver_algebra_json<F: Field>(
        algebra: &QuiverAlgebra<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        let original = algebra.original_path_algebra();
        let dimension = match algebra {
            QuiverAlgebra::Path(path_algebra) => path_algebra
                .canonical_basis()
                .expect("oracle predefined path algebra is finite")
                .len(),
            QuiverAlgebra::Quotient(quotient) => {
                let AlgebraDimension::Finite(dimension) = quotient.dimension().unwrap() else {
                    panic!("oracle predefined quotient is finite");
                };
                usize::try_from(dimension).unwrap()
            }
        };
        json!({
            "path_algebra": algebra.is_path_algebra(),
            "dimension": dimension,
            "vertices": original
                .quiver()
                .vertices()
                .map(|vertex| vertex.name().to_owned())
                .collect::<Vec<_>>(),
            "arrows": original
                .quiver()
                .arrows()
                .map(|arrow| json!({
                    "name": arrow.name(),
                    "source": arrow.source().index().get(),
                    "target": arrow.target().index().get()
                }))
                .collect::<Vec<_>>(),
            "relations": algebra
                .relations()
                .iter()
                .map(|relation| element_json(relation, coefficient_json))
                .collect::<Vec<_>>()
        })
    }

    fn constructed<F: Field>(result: NakayamaAlgebraResult<F>) -> PredefinedAlgebra<F> {
        let NakayamaAlgebraResult::Algebra(algebra) = result else {
            panic!("fixture sequence must be admissible");
        };
        algebra
    }

    fn brauer_vertex(name: &str, multiplicity: i64) -> BrauerConfigurationVertex {
        BrauerConfigurationVertex {
            name: name.to_owned(),
            multiplicity: Integer::from(multiplicity),
        }
    }

    fn brauer_polygon(name: &str, vertices: &[&str]) -> BrauerConfigurationPolygon {
        BrauerConfigurationPolygon {
            name: name.to_owned(),
            vertices: vertices.iter().map(ToString::to_string).collect(),
        }
    }

    fn brauer_orientations(rows: &[&[&str]]) -> Vec<Vec<String>> {
        rows.iter()
            .map(|row| row.iter().map(ToString::to_string).collect())
            .collect()
    }

    fn trivial_brauer_configuration() -> BrauerConfiguration {
        BrauerConfiguration {
            vertices: vec![brauer_vertex("v", 1)],
            polygons: vec![brauer_polygon("e", &["v"])],
            orientations: brauer_orientations(&[&["e"]]),
        }
    }

    fn documented_brauer_configuration() -> BrauerConfiguration {
        BrauerConfiguration {
            vertices: vec![
                brauer_vertex("v1", 1),
                brauer_vertex("v2", 1),
                brauer_vertex("v3", 2),
            ],
            polygons: vec![
                brauer_polygon("e1", &["v1", "v2"]),
                brauer_polygon("e2", &["v2", "v3"]),
            ],
            orientations: brauer_orientations(&[&["e1"], &["e2", "e1"], &["e2"]]),
        }
    }

    fn repeated_brauer_configuration() -> BrauerConfiguration {
        BrauerConfiguration {
            vertices: vec![
                brauer_vertex("a", 1),
                brauer_vertex("b", 1),
                brauer_vertex("c", 1),
            ],
            polygons: vec![
                brauer_polygon("1", &["a", "a"]),
                brauer_polygon("2", &["a", "b"]),
                brauer_polygon("3", &["b", "a"]),
                brauer_polygon("4", &["b", "c"]),
            ],
            orientations: brauer_orientations(&[&["1", "2", "3"], &["2", "4", "3"], &["4"]]),
        }
    }

    fn constructed_brauer<F: Field>(
        result: BrauerConfigurationAlgebraResult<F>,
    ) -> BrauerConfigurationAlgebra<F> {
        let BrauerConfigurationAlgebraResult::Algebra(algebra) = result else {
            panic!("fixture configuration must be valid");
        };
        algebra
    }

    fn invalid_brauer(
        configuration: &mut BrauerConfiguration,
    ) -> BrauerConfigurationValidationError {
        let result = brauer_configuration_algebra(Q, configuration).unwrap();
        let BrauerConfigurationAlgebraResult::InvalidConfiguration(error) = result else {
            panic!("invalid fixture unexpectedly constructed an algebra");
        };
        error
    }

    #[test]
    fn invalid_and_error_branches_follow_qpa_dispatch() {
        assert_eq!(
            nakayama_algebra(Q, &[]),
            Err(PredefinedAlgebraError::EmptyAdmissibleSequence)
        );

        let invalid_sequence = integers(&[1, 1]);
        let invalid = nakayama_algebra(Q, &invalid_sequence).unwrap();
        assert_eq!(
            invalid,
            NakayamaAlgebraResult::InvalidSequence(invalid_sequence.clone())
        );
        assert_eq!(
            invalid.diagnostic(),
            Some(INVALID_ADMISSIBLE_SEQUENCE_DIAGNOSTIC)
        );
        assert!(invalid.as_algebra().is_none());
        assert!(matches!(
            nakayama_algebra_sequence_first(&invalid_sequence, Q).unwrap(),
            NakayamaAlgebraResult::InvalidSequence(sequence) if sequence == invalid_sequence
        ));
        assert!(matches!(
            nakayama_algebra(Q, &integers(&[0])),
            Err(PredefinedAlgebraError::Gbnp(_))
        ));
        assert!(matches!(
            nakayama_algebra(Q, &integers(&[-2])),
            Err(PredefinedAlgebraError::Gbnp(_))
        ));

        let two = Q.nth(Integer::from(2));
        assert_eq!(
            canonical_algebra(Q, &integers(&[2, 2]), std::slice::from_ref(&two)),
            Err(PredefinedAlgebraError::CanonicalArmCoefficientMismatch {
                arms: 2,
                coefficients: 1
            })
        );
        assert_eq!(
            canonical_algebra(Q, &integers(&[2, 1]), &[]),
            Err(PredefinedAlgebraError::CanonicalWeightBelowTwo {
                position: 2,
                weight: Integer::from(1)
            })
        );
        assert_eq!(
            canonical_algebra(Q, &integers(&[2, 2, 2]), &[Q.zero()]),
            Err(PredefinedAlgebraError::ZeroCanonicalRelationCoefficient { position: 1 })
        );
        assert_eq!(
            canonical_algebra_two_weights(Q, &integers(&[2, 2, 2])),
            Err(PredefinedAlgebraError::CanonicalTwoArgumentWeightCount { weights: 3 })
        );
        assert_eq!(
            kronecker_algebra(Q, &Integer::from(1)),
            Err(PredefinedAlgebraError::KroneckerNeedsTwoArrows {
                number_of_arrows: Integer::from(1)
            })
        );

        let nakayama = constructed(nakayama_algebra(Q, &integers(&[1])).unwrap());
        assert_eq!(nakayama.is_nakayama_algebra(), Some(true));
        assert_eq!(nakayama.is_canonical_algebra(), None);
        assert_eq!(nakayama.is_kronecker_algebra(), None);
        assert_eq!(nakayama.is_finite_global_dimension_algebra(), None);
        assert_eq!(nakayama.is_hereditary_algebra(), None);
    }

    #[test]
    fn pinned_qpa_predefined_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_predefined.schema-v1.json"
        ))
        .unwrap();

        let nakayama_singleton = constructed(nakayama_algebra(Q, &integers(&[1])).unwrap());
        let nakayama_singleton_quotient =
            constructed(nakayama_algebra(Q, &integers(&[2])).unwrap());
        let nakayama_linear_path = constructed(nakayama_algebra(Q, &integers(&[3, 2, 1])).unwrap());
        let nakayama_linear_quotient =
            constructed(nakayama_algebra(Q, &integers(&[3, 3, 2, 1])).unwrap());
        let nakayama_cyclic = constructed(nakayama_algebra(Q, &integers(&[3, 3, 3])).unwrap());
        let nakayama_long_cyclic =
            constructed(nakayama_algebra(Q, &integers(&[7, 6, 7, 7])).unwrap());
        let canonical_two = canonical_algebra_two_weights(Q, &integers(&[2, 3])).unwrap();
        let canonical_four = canonical_algebra(
            Q,
            &integers(&[2, 2, 3, 4]),
            &[Q.nth(Integer::from(2)), Q.nth(Integer::from(-3))],
        )
        .unwrap();

        let rational = json!({
            "nakayama_singleton": algebra_json(&nakayama_singleton, &rational_coefficient),
            "nakayama_singleton_quotient": algebra_json(
                &nakayama_singleton_quotient,
                &rational_coefficient
            ),
            "nakayama_linear_path": algebra_json(&nakayama_linear_path, &rational_coefficient),
            "nakayama_linear_quotient": algebra_json(
                &nakayama_linear_quotient,
                &rational_coefficient
            ),
            "nakayama_cyclic": algebra_json(&nakayama_cyclic, &rational_coefficient),
            "nakayama_long_cyclic": algebra_json(&nakayama_long_cyclic, &rational_coefficient),
            "canonical_two": algebra_json(&canonical_two, &rational_coefficient),
            "canonical_four": algebra_json(&canonical_four, &rational_coefficient),
            "properties": {
                "nakayama": [
                    nakayama_cyclic.is_nakayama_algebra().unwrap(),
                    nakayama_cyclic.is_finite_type_algebra().unwrap()
                ],
                "canonical": [
                    canonical_four.is_canonical_algebra().unwrap(),
                    canonical_four
                        .is_finite_global_dimension_algebra()
                        .unwrap()
                ]
            }
        });

        let nakayama_cyclic_f2 = constructed(nakayama_algebra(Z2, &integers(&[3, 3, 3])).unwrap());
        let canonical_f2 =
            canonical_algebra(Z2, &integers(&[2, 3, 4]), std::slice::from_ref(&Z2.one())).unwrap();
        let kronecker_f2 = kronecker_algebra(Z2, &Integer::from(3)).unwrap();
        let z2_coefficient = |coefficient: &u8| Value::from(*coefficient);
        let f2 = json!({
            "nakayama_cyclic": algebra_json(&nakayama_cyclic_f2, &z2_coefficient),
            "canonical": algebra_json(&canonical_f2, &z2_coefficient),
            "kronecker": algebra_json(&kronecker_f2, &z2_coefficient),
            "properties": {
                "nakayama": [
                    nakayama_cyclic_f2.is_nakayama_algebra().unwrap(),
                    nakayama_cyclic_f2.is_finite_type_algebra().unwrap()
                ],
                "canonical": [
                    canonical_f2.is_canonical_algebra().unwrap(),
                    canonical_f2
                        .is_finite_global_dimension_algebra()
                        .unwrap()
                ],
                "kronecker": [
                    kronecker_f2.is_kronecker_algebra().unwrap(),
                    kronecker_f2
                        .is_finite_global_dimension_algebra()
                        .unwrap(),
                    kronecker_f2.is_hereditary_algebra().unwrap()
                ]
            }
        });

        assert_eq!(rational, expected["rational"]);
        assert_eq!(f2, expected["f2"]);
    }

    #[test]
    fn poset_presentation_preserves_qpa_shadowed_relation_loop() {
        let diamond =
            Poset::new(vec![4, 2, 1, 3], &[vec![1, 2, 3], vec![2, 4], vec![3, 4]]).unwrap();
        let mut integer_string = |value: &i32| value.to_string();
        let (_, diamond_relations) =
            poset_algebra_presentation(&diamond, &mut integer_string, Q).unwrap();
        assert_eq!(diamond_relations.len(), 3);
        assert!(
            diamond_relations
                .iter()
                .all(|relation| relation == &diamond_relations[0])
        );

        let parallel = Poset::new(
            vec![1, 2, 3, 4, 5],
            &[vec![1, 2, 3, 4], vec![2, 5], vec![3, 5], vec![4, 5]],
        )
        .unwrap();
        let (_, parallel_relations) =
            poset_algebra_presentation(&parallel, &mut integer_string, Q).unwrap();
        assert_eq!(parallel_relations.len(), 8);
        assert_eq!(parallel_relations[0], parallel_relations[2]);
        assert_eq!(parallel_relations[1], parallel_relations[3]);

        let named =
            poset_algebra_with_string_function(Q, &diamond, |value| format!("p{value}")).unwrap();
        assert!(named.poset_of_poset_algebra().is_identical(&diamond));
        assert_eq!(named.is_poset_algebra(), Some(true));
        assert_eq!(named.is_finite_global_dimension_algebra(), Some(true));
        assert_eq!(
            named
                .algebra()
                .original_path_algebra()
                .quiver()
                .vertices()
                .map(|vertex| vertex.name().to_owned())
                .collect::<Vec<_>>(),
            ["p1", "p2", "p3", "p4"]
        );
    }

    #[test]
    fn pinned_qpa_poset_algebra_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_poset_algebra.schema-v1.json"
        ))
        .unwrap();

        let singleton_poset = Poset::new(vec!["only"], &[]).unwrap();
        let singleton = poset_algebra(Q, &singleton_poset).unwrap();
        let chain_poset =
            Poset::new(vec!["c", "a", "b"], &[vec!["a", "b"], vec!["b", "c"]]).unwrap();
        let chain = poset_algebra(Q, &chain_poset).unwrap();
        let diamond_poset =
            Poset::new(vec![4, 2, 1, 3], &[vec![1, 2, 3], vec![2, 4], vec![3, 4]]).unwrap();
        let diamond = poset_algebra(Q, &diamond_poset).unwrap();
        let rational = json!({
            "singleton": quiver_algebra_json(singleton.algebra(), &rational_coefficient),
            "chain": quiver_algebra_json(chain.algebra(), &rational_coefficient),
            "diamond": quiver_algebra_json(diamond.algebra(), &rational_coefficient),
            "attribute_identity": [
                singleton.poset_of_poset_algebra().is_identical(&singleton_poset),
                chain.poset_of_poset_algebra().is_identical(&chain_poset),
                diamond.poset_of_poset_algebra().is_identical(&diamond_poset)
            ],
            "properties": [
                diamond.is_poset_algebra().unwrap(),
                diamond.is_finite_global_dimension_algebra().unwrap()
            ]
        });

        let parallel_poset = Poset::new(
            vec![1, 2, 3, 4, 5],
            &[vec![1, 2, 3, 4], vec![2, 5], vec![3, 5], vec![4, 5]],
        )
        .unwrap();
        let parallel = poset_algebra(Z2, &parallel_poset).unwrap();
        let non_graded_poset = Poset::new(
            vec![1, 2, 3, 4, 5],
            &[vec![1, 2, 3], vec![2, 5], vec![3, 4], vec![4, 5]],
        )
        .unwrap();
        let non_graded = poset_algebra(Z2, &non_graded_poset).unwrap();
        let z2_coefficient = |coefficient: &u8| Value::from(*coefficient);
        let f2 = json!({
            "parallel": quiver_algebra_json(parallel.algebra(), &z2_coefficient),
            "non_graded": quiver_algebra_json(non_graded.algebra(), &z2_coefficient),
            "attribute_identity": [
                parallel.poset_of_poset_algebra().is_identical(&parallel_poset),
                non_graded
                    .poset_of_poset_algebra()
                    .is_identical(&non_graded_poset)
            ],
            "properties": [
                non_graded.is_poset_algebra().unwrap(),
                non_graded.is_finite_global_dimension_algebra().unwrap()
            ]
        });

        assert_eq!(rational, expected["rational"]);
        assert_eq!(f2, expected["f2"]);
    }

    #[test]
    fn brauer_configuration_vertex_and_polygon_validation_follow_qpa() {
        let mut no_vertices = BrauerConfiguration {
            vertices: Vec::new(),
            polygons: Vec::new(),
            orientations: Vec::new(),
        };
        assert_eq!(
            invalid_brauer(&mut no_vertices),
            BrauerConfigurationValidationError::NoVertices
        );

        let mut nonpositive = trivial_brauer_configuration();
        nonpositive.vertices[0].multiplicity = Integer::from(0);
        assert_eq!(
            invalid_brauer(&mut nonpositive),
            BrauerConfigurationValidationError::NonPositiveMultiplicity { position: 1 }
        );

        let mut duplicate_vertex = trivial_brauer_configuration();
        duplicate_vertex.vertices.push(brauer_vertex("v", 1));
        assert_eq!(
            invalid_brauer(&mut duplicate_vertex),
            BrauerConfigurationValidationError::DuplicateVertexName { name: "v".into() }
        );

        let mut empty_polygon = trivial_brauer_configuration();
        empty_polygon.polygons[0].vertices.clear();
        assert_eq!(
            invalid_brauer(&mut empty_polygon),
            BrauerConfigurationValidationError::PolygonWithoutVertex { position: 1 }
        );

        let mut polygon_vertex_name = trivial_brauer_configuration();
        polygon_vertex_name.polygons[0].name = "v".into();
        polygon_vertex_name.orientations[0][0] = "v".into();
        assert_eq!(
            invalid_brauer(&mut polygon_vertex_name),
            BrauerConfigurationValidationError::PolygonNameMatchesVertex { name: "v".into() }
        );

        let mut unknown_vertex = trivial_brauer_configuration();
        unknown_vertex.polygons[0].vertices[0] = "x".into();
        assert_eq!(
            invalid_brauer(&mut unknown_vertex),
            BrauerConfigurationValidationError::PolygonContainsUnknownVertex {
                polygon: "e".into(),
                vertex: "x".into()
            }
        );

        let mut duplicate_polygon = trivial_brauer_configuration();
        duplicate_polygon.polygons.push(brauer_polygon("e", &["v"]));
        assert_eq!(
            invalid_brauer(&mut duplicate_polygon),
            BrauerConfigurationValidationError::DuplicatePolygonName { name: "e".into() }
        );
    }

    #[test]
    fn brauer_configuration_orientation_validation_and_mutation_follow_qpa() {
        let mut orientation_count = trivial_brauer_configuration();
        orientation_count.orientations.clear();
        assert_eq!(
            invalid_brauer(&mut orientation_count),
            BrauerConfigurationValidationError::OrientationCountMismatch {
                expected: 1,
                actual: 0
            }
        );

        let mut empty_orientation = trivial_brauer_configuration();
        empty_orientation.orientations[0].clear();
        assert_eq!(
            invalid_brauer(&mut empty_orientation),
            BrauerConfigurationValidationError::EmptyOrientation { position: 1 }
        );

        let mut unknown_polygon = trivial_brauer_configuration();
        unknown_polygon.orientations[0][0] = "x".into();
        assert_eq!(
            invalid_brauer(&mut unknown_polygon),
            BrauerConfigurationValidationError::OrientationContainsUnknownPolygon {
                position: 1,
                polygon: "x".into()
            }
        );

        let mut containment = BrauerConfiguration {
            vertices: vec![brauer_vertex("a", 1), brauer_vertex("b", 1)],
            polygons: vec![brauer_polygon("z", &["a"]), brauer_polygon("x", &["b"])],
            orientations: brauer_orientations(&[&["z"], &["z", "x"]]),
        };
        let containment_error = invalid_brauer(&mut containment);
        assert_eq!(
            containment_error,
            BrauerConfigurationValidationError::OrientationContainmentMismatch { position: 2 }
        );
        assert_eq!(
            containment.orientations,
            brauer_orientations(&[&["z"], &["x", "z"]])
        );
    }

    #[test]
    fn brauer_configuration_diagnostics_follow_qpa() {
        let expected_diagnostics = [
            (
                BrauerConfigurationValidationError::NoVertices,
                "There must be at least 1 vertex\n",
            ),
            (
                BrauerConfigurationValidationError::NonPositiveMultiplicity { position: 1 },
                "Multiplicities must be positive integers.\n",
            ),
            (
                BrauerConfigurationValidationError::DuplicateVertexName { name: "v".into() },
                "Each vertex must have a unique name.\n",
            ),
            (
                BrauerConfigurationValidationError::PolygonWithoutVertex { position: 1 },
                "Edges or Polygons should be lists of length 2 or greater.\n",
            ),
            (
                BrauerConfigurationValidationError::PolygonNameMatchesVertex { name: "v".into() },
                "An edge or polygon cannot have the same name as a vertex.\n",
            ),
            (
                BrauerConfigurationValidationError::PolygonContainsUnknownVertex {
                    polygon: "e".into(),
                    vertex: "x".into(),
                },
                "An edge or polygon cannot contain a vertex which was not listed.\n",
            ),
            (
                BrauerConfigurationValidationError::DuplicatePolygonName { name: "e".into() },
                "Edges or polygons must have unique names.\n",
            ),
            (
                BrauerConfigurationValidationError::OrientationCountMismatch {
                    expected: 1,
                    actual: 0,
                },
                "There must be an orientation corresponding to each vertex.\n",
            ),
            (
                BrauerConfigurationValidationError::EmptyOrientation { position: 1 },
                "Orientations must be lists of length 1 or greater.\n",
            ),
            (
                BrauerConfigurationValidationError::OrientationContainsUnknownPolygon {
                    position: 1,
                    polygon: "x".into(),
                },
                "Orientations may only contain edges or polygons which were listed.\n",
            ),
            (
                BrauerConfigurationValidationError::OrientationContainmentMismatch { position: 2 },
                "Orientations must contain exactly the edges or polygons which contain the corresponding vertex.\n",
            ),
        ];
        for (error, diagnostic) in expected_diagnostics {
            assert_eq!(error.diagnostic(), diagnostic);
        }
    }

    #[test]
    fn brauer_configuration_raw_relation_families_follow_qpa_loops() {
        let mut trivial = trivial_brauer_configuration();
        constructed_brauer(brauer_configuration_algebra(Q, &mut trivial).unwrap());
        let trivial_presentation = brauer_configuration_algebra_presentation(Q, &trivial).unwrap();
        assert_eq!(
            trivial_presentation
                .relation_families
                .each_ref()
                .map(Vec::len),
            [0, 0, 0]
        );

        let mut documented = documented_brauer_configuration();
        constructed_brauer(brauer_configuration_algebra(Q, &mut documented).unwrap());
        let documented_presentation =
            brauer_configuration_algebra_presentation(Q, &documented).unwrap();
        assert_eq!(
            documented_presentation
                .relation_families
                .each_ref()
                .map(Vec::len),
            [4, 1, 3]
        );
        assert_eq!(
            documented_presentation.relation_families[0]
                .iter()
                .filter(|relation| relation.is_zero())
                .count(),
            2
        );

        let mut repeated = repeated_brauer_configuration();
        constructed_brauer(brauer_configuration_algebra(Z2, &mut repeated).unwrap());
        let repeated_presentation =
            brauer_configuration_algebra_presentation(Z2, &repeated).unwrap();
        assert_eq!(
            repeated_presentation
                .relation_families
                .each_ref()
                .map(Vec::len),
            [18, 2, 6]
        );
    }

    #[test]
    fn pinned_qpa_brauer_configuration_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_brauer_configuration.schema-v1.json"
        ))
        .unwrap();

        let mut trivial_configuration = trivial_brauer_configuration();
        let trivial = constructed_brauer(
            brauer_configuration_algebra(Q, &mut trivial_configuration).unwrap(),
        );
        let mut documented_configuration = documented_brauer_configuration();
        let documented = constructed_brauer(
            brauer_configuration_algebra(Q, &mut documented_configuration).unwrap(),
        );
        let rational = json!({
            "trivial": quiver_algebra_json(trivial.algebra(), &rational_coefficient),
            "documented": quiver_algebra_json(documented.algebra(), &rational_coefficient),
            "orientations_after": documented_configuration.orientations,
            "properties": [
                trivial.is_symmetric_algebra().unwrap(),
                documented.is_symmetric_algebra().unwrap()
            ]
        });

        let mut repeated_configuration = repeated_brauer_configuration();
        let repeated = constructed_brauer(
            brauer_configuration_algebra(Z2, &mut repeated_configuration).unwrap(),
        );
        let z2_coefficient = |coefficient: &u8| Value::from(*coefficient);
        let f2 = json!({
            "repeated": quiver_algebra_json(repeated.algebra(), &z2_coefficient),
            "orientations_after": repeated_configuration.orientations,
            "properties": [repeated.is_symmetric_algebra().unwrap()]
        });

        assert_eq!(rational, expected["rational"]);
        assert_eq!(f2, expected["f2"]);
    }
}
