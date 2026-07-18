//! Quivers and paths with QPA's identity and multiplication conventions.
//!
//! This module translates the foundational methods in QPA 1.37
//! `lib/quiver.gi`: the three `Quiver` constructors, path external
//! representations, path multiplication, the default length-left-lexicographic
//! ordering, `NextPath`, and the opposite/separated/double constructions.

use std::{
    cmp::Ordering,
    collections::{HashSet, VecDeque},
    fmt,
    hash::{Hash, Hasher},
    num::NonZeroUsize,
    sync::{
        Arc, OnceLock, RwLock, Weak,
        atomic::{AtomicU64, Ordering as AtomicOrdering},
    },
};

use serde::{Deserialize, Serialize};
use thiserror::Error;

use numerica::domains::integer::Integer;

use crate::{OrderingError, QuiverOrdering};

static NEXT_QUIVER_ID: AtomicU64 = AtomicU64::new(1);

/// Rust logging target corresponding to QPA's declared `InfoQuiver` class.
pub const INFO_QUIVER: &str = "qpa::quiver";

/// Translate QPA's public `UniqueQuiverName` counter.
///
/// Quiver construction consumes values from the same counter, just as the
/// source uses this function to name every newly allocated path family.
#[must_use]
pub fn unique_quiver_name() -> String {
    format!("Quiver_{}", next_quiver_id())
}

fn next_quiver_id() -> u64 {
    NEXT_QUIVER_ID.fetch_add(1, AtomicOrdering::Relaxed)
}

/// A one-based vertex index, matching QPA's `VertexIndex` result.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(transparent)]
pub struct VertexIndex(NonZeroUsize);

impl VertexIndex {
    /// Construct a one-based vertex index.
    #[must_use]
    pub const fn new(index: usize) -> Option<Self> {
        match NonZeroUsize::new(index) {
            Some(index) => Some(Self(index)),
            None => None,
        }
    }

    /// Return the QPA-compatible one-based index.
    #[must_use]
    pub const fn get(self) -> usize {
        self.0.get()
    }

    pub(crate) const fn from_zero_based(index: usize) -> Self {
        let Some(index) = NonZeroUsize::new(index + 1) else {
            unreachable!();
        };
        Self(index)
    }

    pub(crate) const fn zero_based(self) -> usize {
        self.get() - 1
    }
}

/// A one-based arrow index, matching QPA's `ArrowIndex` result.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(transparent)]
pub struct ArrowIndex(NonZeroUsize);

impl ArrowIndex {
    /// Construct a one-based arrow index.
    #[must_use]
    pub const fn new(index: usize) -> Option<Self> {
        match NonZeroUsize::new(index) {
            Some(index) => Some(Self(index)),
            None => None,
        }
    }

    /// Return the QPA-compatible one-based arrow index.
    #[must_use]
    pub const fn get(self) -> usize {
        self.0.get()
    }

    pub(crate) const fn from_zero_based(index: usize) -> Self {
        let Some(index) = NonZeroUsize::new(index + 1) else {
            unreachable!();
        };
        Self(index)
    }

    pub(crate) const fn zero_based(self) -> usize {
        self.get() - 1
    }
}

/// An arrow specification using QPA-compatible one-based vertex indices.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct ArrowSpec {
    /// Source vertex.
    pub source: VertexIndex,
    /// Target vertex.
    pub target: VertexIndex,
    /// Optional explicit name. If omitted, QPA's `a1`, `a2`, ... names apply.
    pub name: Option<Box<str>>,
}

impl ArrowSpec {
    /// Construct an explicitly named arrow specification.
    #[must_use]
    pub fn named(source: VertexIndex, target: VertexIndex, name: impl Into<Box<str>>) -> Self {
        Self {
            source,
            target,
            name: Some(name.into()),
        }
    }

    /// Construct an arrow to which QPA assigns its default name.
    #[must_use]
    pub const fn unnamed(source: VertexIndex, target: VertexIndex) -> Self {
        Self {
            source,
            target,
            name: None,
        }
    }
}

/// An arrow specification whose endpoints are vertex names.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct NamedArrowSpec {
    /// Source vertex name.
    pub source: Box<str>,
    /// Target vertex name.
    pub target: Box<str>,
    /// Optional explicit name. If omitted, QPA's `a1`, `a2`, ... names apply.
    pub name: Option<Box<str>>,
}

/// An orientation symbol used by QPA's Dynkin-quiver constructors.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum DynkinOrientation {
    /// QPA's `"l"` orientation.
    Left,
    /// QPA's `"r"` orientation.
    Right,
    /// QPA's `"u"` orientation for the exceptional arm of `E_n`.
    Up,
    /// QPA's `"d"` orientation for the exceptional arm of `E_n`.
    Down,
}

/// A Dynkin diagram family recognized by QPA.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum DynkinFamily {
    /// Type A.
    A,
    /// Type D.
    D,
    /// Type E.
    E,
}

impl DynkinFamily {
    /// Parse the three strings dispatched by QPA's generic `DynkinQuiver`.
    /// An unknown string reaches the end of the GAP method and returns `fail`.
    #[must_use]
    pub fn from_qpa_name(name: &str) -> Option<Self> {
        match name {
            "A" => Some(Self::A),
            "D" => Some(Self::D),
            "E" => Some(Self::E),
            _ => None,
        }
    }
}

impl fmt::Display for DynkinFamily {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::A => "A",
            Self::D => "D",
            Self::E => "E",
        })
    }
}

/// The Dynkin type printed by QPA's `IsDynkinQuiver` method.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct DynkinType {
    /// Dynkin diagram family.
    pub family: DynkinFamily,
    /// Number of vertices.
    pub rank: usize,
}

impl fmt::Display for DynkinType {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}_{}", self.family, self.rank)
    }
}

/// QPA's diagnostic when `IsDynkinQuiver` receives a disconnected quiver.
pub const DYNKIN_NOT_CONNECTED_MESSAGE: &str = "Quiver is not connected.\n";

/// QPA's diagnostic when `IsDynkinQuiver` rejects the quiver by arrow count.
pub const DYNKIN_UNORIENTED_CYCLE_MESSAGE: &str = "Quiver contains an (un)oriented cycle.\n";

/// QPA's Boolean Dynkin decision together with its ordered printed output.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DynkinDecision {
    dynkin_type: Option<DynkinType>,
    messages: Vec<Box<str>>,
}

impl DynkinDecision {
    /// Return the recognized Dynkin type, or `None` for QPA's `false`.
    #[must_use]
    pub const fn dynkin_type(&self) -> Option<DynkinType> {
        self.dynkin_type
    }

    /// Return QPA's Boolean result.
    #[must_use]
    pub const fn is_dynkin(&self) -> bool {
        self.dynkin_type.is_some()
    }

    /// Return the lines printed by QPA in source order.
    #[must_use]
    pub fn messages(&self) -> &[Box<str>] {
        &self.messages
    }
}

/// The finite/infinite result of QPA's `Size` method for a quiver.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum QuiverCardinality {
    /// An arbitrary-precision finite size, including QPA's zero path.
    Finite(Integer),
    /// QPA's `infinity`, occurring exactly when the quiver has a cycle.
    Infinite,
}

/// QPA's second argument to `NextPath` after the initial integer zero.
///
/// `Start` represents the source call `NextPath(Q, 0)`. `Queue` is the
/// shallow-copied list returned in the second component of a successful call.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum NextPathPosition {
    /// QPA's initial integer-zero position.
    Start,
    /// The remaining breadth-first path queue.
    Queue(Vec<Path>),
}

impl NamedArrowSpec {
    /// Construct an explicitly named arrow specification.
    #[must_use]
    pub fn named(
        source: impl Into<Box<str>>,
        target: impl Into<Box<str>>,
        name: impl Into<Box<str>>,
    ) -> Self {
        Self {
            source: source.into(),
            target: target.into(),
            name: Some(name.into()),
        }
    }

    /// Construct an arrow to which QPA assigns its default name.
    #[must_use]
    pub fn unnamed(source: impl Into<Box<str>>, target: impl Into<Box<str>>) -> Self {
        Self {
            source: source.into(),
            target: target.into(),
            name: None,
        }
    }
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct ArrowData {
    source: VertexIndex,
    target: VertexIndex,
    name: Box<str>,
}

#[derive(Clone, Copy, Debug)]
enum OppositeNameMap {
    AppendSuffix,
    RemoveSuffix,
}

struct QuiverData {
    id: u64,
    vertex_names: Box<[Box<str>]>,
    arrows: Box<[ArrowData]>,
    adjacency: Box<[Box<[usize]>]>,
    incoming: Box<[Box<[ArrowIndex]>]>,
    outgoing: Box<[Box<[ArrowIndex]>]>,
    neighbors: Box<[Box<[VertexIndex]>]>,
    special_arrow_names: bool,
    default_ordering_identity: RwLock<Option<Weak<()>>>,
    opposite_weak: RwLock<Option<Weak<QuiverData>>>,
    opposite_strong: OnceLock<Quiver>,
    opposite_name_map: OnceLock<OppositeNameMap>,
}

/// A QPA quiver object.
///
/// Cloning this handle preserves QPA object identity. Constructing another
/// quiver from the same specifications creates a distinct path family.
#[derive(Clone)]
pub struct Quiver(Arc<QuiverData>);

impl Quiver {
    /// Translate `Quiver(N, arrow_spec)` from QPA.
    ///
    /// Vertices are named `v1`, `v2`, ... . All arrow specifications must
    /// either have explicit names or omit them, mirroring the uniform GAP list
    /// row length required by QPA.
    ///
    /// # Errors
    ///
    /// Returns an error for zero vertices, mixed named/unnamed arrows, or an
    /// endpoint outside `1..=number_of_vertices`.
    pub fn from_indexed(
        number_of_vertices: usize,
        arrows: Vec<ArrowSpec>,
    ) -> Result<Self, QuiverError> {
        if number_of_vertices == 0 {
            return Err(QuiverError::IndexedQuiverNeedsVertex);
        }
        ensure_uniform_naming(arrows.iter().map(|arrow| arrow.name.is_some()))?;

        let vertex_names = (1..=number_of_vertices)
            .map(|index| format!("v{index}").into_boxed_str())
            .collect();
        let resolved = arrows
            .into_iter()
            .enumerate()
            .map(|(position, arrow)| {
                validate_vertex_index(arrow.source, number_of_vertices)?;
                validate_vertex_index(arrow.target, number_of_vertices)?;
                Ok(ArrowData {
                    source: arrow.source,
                    target: arrow.target,
                    name: arrow
                        .name
                        .unwrap_or_else(|| format!("a{}", position + 1).into_boxed_str()),
                })
            })
            .collect::<Result<Vec<_>, QuiverError>>()?;

        let special_arrow_names = resolved.iter().enumerate().all(|(position, arrow)| {
            u8::try_from(position)
                .ok()
                .and_then(|position| b'a'.checked_add(position))
                .is_some_and(|name| arrow.name.as_bytes() == [name])
        });
        Self::from_resolved(vertex_names, resolved, special_arrow_names)
    }

    /// Translate `Quiver(vertex_names, arrow_spec)` from QPA.
    ///
    /// QPA resolves a repeated vertex name to its first position; this
    /// constructor intentionally does the same.
    ///
    /// # Errors
    ///
    /// Returns an error for mixed named/unnamed arrows or an unknown endpoint.
    pub fn from_named(
        vertex_names: Vec<Box<str>>,
        arrows: Vec<NamedArrowSpec>,
    ) -> Result<Self, QuiverError> {
        ensure_uniform_naming(arrows.iter().map(|arrow| arrow.name.is_some()))?;

        let resolved = arrows
            .into_iter()
            .enumerate()
            .map(|(position, arrow)| {
                let source = position_of_name(&vertex_names, &arrow.source).ok_or_else(|| {
                    QuiverError::UnknownVertexName {
                        name: arrow.source.clone(),
                    }
                })?;
                let target = position_of_name(&vertex_names, &arrow.target).ok_or_else(|| {
                    QuiverError::UnknownVertexName {
                        name: arrow.target.clone(),
                    }
                })?;
                Ok(ArrowData {
                    source,
                    target,
                    name: arrow
                        .name
                        .unwrap_or_else(|| format!("a{}", position + 1).into_boxed_str()),
                })
            })
            .collect::<Result<Vec<_>, QuiverError>>()?;

        Self::from_resolved(vertex_names, resolved, false)
    }

    /// Translate `Quiver(adjacency_matrix)` from QPA.
    ///
    /// Arrows are generated in row-major order and named `a1`, `a2`, ... .
    ///
    /// # Errors
    ///
    /// Returns an error if the matrix is not square or the number of arrows
    /// cannot be represented in memory.
    pub fn from_adjacency(adjacency: Vec<Vec<usize>>) -> Result<Self, QuiverError> {
        let size = adjacency.len();
        if let Some((row, actual)) = adjacency
            .iter()
            .enumerate()
            .find_map(|(row, entries)| (entries.len() != size).then_some((row, entries.len())))
        {
            return Err(QuiverError::NonSquareAdjacency {
                row,
                actual,
                expected: size,
            });
        }

        let vertex_names = (1..=size)
            .map(|index| format!("v{index}").into_boxed_str())
            .collect::<Vec<_>>();
        let mut arrows = Vec::new();
        for (source, row) in adjacency.iter().enumerate() {
            for (target, &multiplicity) in row.iter().enumerate() {
                for _ in 0..multiplicity {
                    let position = arrows
                        .len()
                        .checked_add(1)
                        .ok_or(QuiverError::ArrowCountOverflow)?;
                    arrows.push(ArrowData {
                        source: VertexIndex::from_zero_based(source),
                        target: VertexIndex::from_zero_based(target),
                        name: format!("a{position}").into_boxed_str(),
                    });
                }
            }
        }

        Self::from_resolved_with_adjacency(vertex_names, arrows, false, adjacency)
    }

    fn from_resolved(
        vertex_names: Vec<Box<str>>,
        arrows: Vec<ArrowData>,
        special_arrow_names: bool,
    ) -> Result<Self, QuiverError> {
        let mut adjacency = vec![vec![0_usize; vertex_names.len()]; vertex_names.len()];
        for arrow in &arrows {
            let entry = &mut adjacency[arrow.source.zero_based()][arrow.target.zero_based()];
            *entry = entry
                .checked_add(1)
                .ok_or(QuiverError::ArrowCountOverflow)?;
        }
        Self::from_resolved_with_adjacency(vertex_names, arrows, special_arrow_names, adjacency)
    }

    fn from_resolved_with_adjacency(
        vertex_names: Vec<Box<str>>,
        arrows: Vec<ArrowData>,
        special_arrow_names: bool,
        adjacency: Vec<Vec<usize>>,
    ) -> Result<Self, QuiverError> {
        let mut incoming = vec![Vec::new(); vertex_names.len()];
        let mut outgoing = vec![Vec::new(); vertex_names.len()];
        let mut neighbors = vec![Vec::new(); vertex_names.len()];

        for (position, arrow) in arrows.iter().enumerate() {
            validate_vertex_index(arrow.source, vertex_names.len())?;
            validate_vertex_index(arrow.target, vertex_names.len())?;
            let arrow_index = ArrowIndex::from_zero_based(position);
            outgoing[arrow.source.zero_based()].push(arrow_index);
            incoming[arrow.target.zero_based()].push(arrow_index);
            let source_neighbors = &mut neighbors[arrow.source.zero_based()];
            if !source_neighbors.contains(&arrow.target) {
                source_neighbors.push(arrow.target);
            }
        }

        Ok(Self(Arc::new(QuiverData {
            id: next_quiver_id(),
            vertex_names: vertex_names.into_boxed_slice(),
            arrows: arrows.into_boxed_slice(),
            adjacency: adjacency
                .into_iter()
                .map(Vec::into_boxed_slice)
                .collect::<Vec<_>>()
                .into_boxed_slice(),
            incoming: incoming
                .into_iter()
                .map(Vec::into_boxed_slice)
                .collect::<Vec<_>>()
                .into_boxed_slice(),
            outgoing: outgoing
                .into_iter()
                .map(Vec::into_boxed_slice)
                .collect::<Vec<_>>()
                .into_boxed_slice(),
            neighbors: neighbors
                .into_iter()
                .map(Vec::into_boxed_slice)
                .collect::<Vec<_>>()
                .into_boxed_slice(),
            special_arrow_names,
            default_ordering_identity: RwLock::new(None),
            opposite_weak: RwLock::new(None),
            opposite_strong: OnceLock::new(),
            opposite_name_map: OnceLock::new(),
        })))
    }

    /// Return the number of vertices.
    #[must_use]
    pub fn number_of_vertices(&self) -> usize {
        self.0.vertex_names.len()
    }

    /// Return the number of arrows.
    #[must_use]
    pub fn number_of_arrows(&self) -> usize {
        self.0.arrows.len()
    }

    /// Return whether two handles denote the identical QPA quiver object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Return whether QPA selected its special `IsQuiverSA` category.
    #[must_use]
    pub fn has_special_single_letter_arrow_names(&self) -> bool {
        self.0.special_arrow_names
    }

    /// Return a vertex by its one-based QPA index.
    #[must_use]
    pub fn vertex(&self, index: VertexIndex) -> Option<Vertex> {
        (index.get() <= self.number_of_vertices()).then(|| Vertex {
            quiver: self.clone(),
            index,
        })
    }

    /// Return an arrow by its one-based QPA arrow index.
    #[must_use]
    pub fn arrow(&self, index: ArrowIndex) -> Option<Arrow> {
        (index.get() <= self.number_of_arrows()).then(|| Arrow {
            quiver: self.clone(),
            index,
        })
    }

    /// Iterate over vertices in QPA generator order.
    #[must_use]
    pub fn vertices(&self) -> impl ExactSizeIterator<Item = Vertex> + '_ {
        (0..self.number_of_vertices()).map(|index| Vertex {
            quiver: self.clone(),
            index: VertexIndex::from_zero_based(index),
        })
    }

    /// Iterate over arrows in QPA generator order.
    #[must_use]
    pub fn arrows(&self) -> impl ExactSizeIterator<Item = Arrow> + '_ {
        (0..self.number_of_arrows()).map(|index| Arrow {
            quiver: self.clone(),
            index: ArrowIndex::from_zero_based(index),
        })
    }

    /// Iterate `GeneratorsOfQuiver`: vertices first, then arrows.
    pub fn generators(&self) -> impl Iterator<Item = Generator> + '_ {
        self.vertices()
            .map(Generator::Vertex)
            .chain(self.arrows().map(Generator::Arrow))
    }

    /// Return QPA's adjacency matrix with rows as sources and columns as targets.
    #[must_use]
    pub fn adjacency_matrix(&self) -> &[Box<[usize]>] {
        &self.0.adjacency
    }

    /// Return the distinguished zero path.
    #[must_use]
    pub fn zero(&self) -> Path {
        Path {
            quiver: self.clone(),
            representation: PathRepresentation::Zero,
        }
    }

    /// Look up QPA's record component by name.
    ///
    /// QPA assigns vertices first and arrows second, so the last generator
    /// having a repeated name wins.
    #[must_use]
    pub fn generator_by_name(&self, name: &str) -> Option<Generator> {
        if let Some((position, _)) = self
            .0
            .arrows
            .iter()
            .enumerate()
            .rev()
            .find(|(_, arrow)| arrow.name.as_ref() == name)
        {
            return self
                .arrow(ArrowIndex::from_zero_based(position))
                .map(Generator::Arrow);
        }
        self.0
            .vertex_names
            .iter()
            .enumerate()
            .rev()
            .find(|(_, vertex_name)| vertex_name.as_ref() == name)
            .and_then(|(position, _)| self.vertex(VertexIndex::from_zero_based(position)))
            .map(Generator::Vertex)
    }

    /// Return incoming arrows in their QPA construction order.
    ///
    /// # Errors
    ///
    /// Returns [`QuiverError::DifferentQuivers`] for a foreign vertex.
    pub fn incoming_arrows(&self, vertex: &Vertex) -> Result<Vec<Arrow>, QuiverError> {
        self.ensure_same_quiver(&vertex.quiver)?;
        Ok(self.0.incoming[vertex.index.zero_based()]
            .iter()
            .filter_map(|&index| self.arrow(index))
            .collect())
    }

    /// Return outgoing arrows in their QPA construction order.
    ///
    /// # Errors
    ///
    /// Returns [`QuiverError::DifferentQuivers`] for a foreign vertex.
    pub fn outgoing_arrows(&self, vertex: &Vertex) -> Result<Vec<Arrow>, QuiverError> {
        self.ensure_same_quiver(&vertex.quiver)?;
        Ok(self.0.outgoing[vertex.index.zero_based()]
            .iter()
            .filter_map(|&index| self.arrow(index))
            .collect())
    }

    /// Return distinct targets of outgoing arrows in first-occurrence order.
    ///
    /// # Errors
    ///
    /// Returns [`QuiverError::DifferentQuivers`] for a foreign vertex.
    pub fn neighbors(&self, vertex: &Vertex) -> Result<Vec<Vertex>, QuiverError> {
        self.ensure_same_quiver(&vertex.quiver)?;
        Ok(self.0.neighbors[vertex.index.zero_based()]
            .iter()
            .filter_map(|&index| self.vertex(index))
            .collect())
    }

    /// Compare paths using QPA's default length-left-lexicographic ordering.
    ///
    /// # Errors
    ///
    /// Returns an error if either path belongs to another quiver.
    pub fn compare_paths(&self, left: &Path, right: &Path) -> Result<Ordering, PathError> {
        left.ensure_quiver(self)?;
        right.ensure_quiver(self)?;

        let by_length = left.ordering_length().cmp(&right.ordering_length());
        if by_length != Ordering::Equal {
            return Ok(by_length);
        }
        if left.is_zero() {
            return Ok(if right.is_zero() {
                Ordering::Equal
            } else {
                Ordering::Less
            });
        }
        if right.is_zero() {
            return Ok(Ordering::Greater);
        }
        Ok(left
            .external_representation()
            .cmp(&right.external_representation()))
    }

    /// Return QPA's default `OrderingOfQuiver` value.
    ///
    /// # Errors
    ///
    /// Construction is infallible for a valid quiver, but the checked result
    /// composes with the other ordering constructors.
    pub fn ordering(&self) -> Result<QuiverOrdering, OrderingError> {
        let mut slot = self
            .0
            .default_ordering_identity
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let identity = slot.as_ref().and_then(Weak::upgrade).unwrap_or_else(|| {
            let identity = Arc::new(());
            *slot = Some(Arc::downgrade(&identity));
            identity
        });
        QuiverOrdering::default_for_with_identity(self, identity)
    }

    /// Translate QPA's `OrderedBy(Q, O)` family-copy operation.
    ///
    /// The returned wrapper retains the identical ordering object `O`, while
    /// all vertices, arrows, and paths belong to a fresh quiver family. QPA
    /// deliberately does not require `O` to have been constructed from `Q`.
    ///
    /// # Errors
    ///
    /// Propagates only an internal reconstruction error. The copied endpoint
    /// indices and adjacency data have already been validated by `Q`.
    pub fn ordered_by(&self, ordering: QuiverOrdering) -> Result<OrderedQuiver, QuiverError> {
        let copied = Self::from_resolved_with_adjacency(
            self.0.vertex_names.to_vec(),
            self.0.arrows.to_vec(),
            false,
            self.0.adjacency.iter().map(|row| row.to_vec()).collect(),
        )?;
        Ok(OrderedQuiver {
            quiver: copied,
            ordering,
        })
    }

    /// Iterate paths in QPA's `NextPath` breadth-first order.
    ///
    /// The zero path is not yielded, matching QPA's iterator. Cyclic quivers
    /// produce an infinite iterator.
    #[must_use]
    pub fn paths(&self) -> PathIterator {
        PathIterator {
            quiver: self.clone(),
            position: NextPathPosition::Start,
        }
    }

    /// Translate one call to QPA's public `NextPath(Q, position)` helper.
    ///
    /// The input position is not mutated: QPA begins with `ShallowCopy` and
    /// returns the new queue beside the selected path.
    ///
    /// # Errors
    ///
    /// A queue beginning with the distinguished zero path reaches a QPA
    /// missing-method error when it requests outgoing arrows of that endpoint;
    /// Rust reports that branch explicitly.
    pub fn next_path(
        &self,
        position: &NextPathPosition,
    ) -> Result<Option<(Path, NextPathPosition)>, PathError> {
        let mut queue = match position {
            NextPathPosition::Start => self
                .vertices()
                .map(|vertex| vertex.path())
                .collect::<VecDeque<_>>(),
            NextPathPosition::Queue(queue) => queue.iter().cloned().collect(),
        };
        let Some(result) = queue.pop_front() else {
            return Ok(None);
        };
        let PathEndpoint::Vertex(target) = result.target() else {
            return Err(PathError::NextPathZeroEndpoint);
        };
        for arrow_index in &target.quiver.0.outgoing[target.index.zero_based()] {
            let arrow =
                target
                    .quiver
                    .arrow(*arrow_index)
                    .ok_or(PathError::GeneratorIndexOutOfRange {
                        index: target.quiver.number_of_vertices() + arrow_index.get(),
                        number_of_generators: target.quiver.number_of_vertices()
                            + target.quiver.number_of_arrows(),
                    })?;
            queue.push_back(result.multiply(&arrow.path())?);
        }
        Ok(Some((
            result,
            NextPathPosition::Queue(queue.into_iter().collect()),
        )))
    }

    /// Translate one-based indexing of QPA's quiver enumerator.
    #[must_use]
    pub fn path_at(&self, number: usize) -> Option<Path> {
        if number == 0 {
            return None;
        }
        let mut position = NextPathPosition::Start;
        let mut result = None;
        for _ in 0..number {
            let (path, next) = self.next_path(&position).ok()??;
            result = Some(path);
            position = next;
        }
        result
    }

    /// Return QPA's concise `ViewObj` text.
    #[must_use]
    pub fn view_string(&self) -> String {
        format!(
            "<quiver with {} vertices and {} arrows>",
            self.number_of_vertices(),
            self.number_of_arrows()
        )
    }

    /// Test acyclicity using the depth-first search in
    /// `quivergraphalgorithms.gi`.
    #[must_use]
    pub fn is_acyclic(&self) -> bool {
        self.topological_order().is_some()
    }

    /// Translate QPA's immediate `IsFinite` method for an acyclic quiver.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.is_acyclic()
    }

    /// Return QPA's DFS topological order, or `None` when a cycle is found.
    #[must_use]
    pub fn topological_order(&self) -> Option<Vec<Vertex>> {
        let mut colors = vec![VisitColor::White; self.number_of_vertices()];
        let mut sorted = Vec::with_capacity(self.number_of_vertices());
        for vertex in 0..self.number_of_vertices() {
            if colors[vertex] == VisitColor::White
                && !visit_acyclic(self, vertex, &mut colors, &mut sorted)
            {
                return None;
            }
        }
        sorted.reverse();
        Some(
            sorted
                .into_iter()
                .filter_map(|index| self.vertex(VertexIndex::from_zero_based(index)))
                .collect(),
        )
    }

    /// Translate QPA's `Size` method for quivers.
    ///
    /// A finite result counts every path and also QPA's distinguished zero
    /// element. A cyclic quiver has cardinality `infinity`.
    #[must_use]
    pub fn cardinality(&self) -> QuiverCardinality {
        let Some(sorted) = self.topological_order() else {
            return QuiverCardinality::Infinite;
        };
        let mut ending_at = vec![Integer::from(0); self.number_of_vertices()];
        let mut total = Integer::from(1);

        for vertex in sorted {
            let source = vertex.index.zero_based();
            let contribution = &ending_at[source] + 1_u8;
            for arrow in &self.0.outgoing[source] {
                let target = self.0.arrows[arrow.zero_based()].target.zero_based();
                ending_at[target] += contribution.clone();
            }
            total += contribution;
        }
        QuiverCardinality::Finite(total)
    }

    /// Test connectedness of the underlying unoriented graph using QPA's DFS.
    #[must_use]
    pub fn is_connected(&self) -> bool {
        if self.number_of_vertices() == 0 {
            return true;
        }
        let mut visited = vec![false; self.number_of_vertices()];
        visit_connected(self, 0, &mut visited, &mut Vec::new());
        visited.into_iter().all(std::convert::identity)
    }

    /// Test QPA's `IsUAcyclicQuiver` property.
    ///
    /// Parallel arrows and loops are unoriented cycles of length two and one,
    /// respectively, exactly as in the source method.
    #[must_use]
    pub fn is_unoriented_acyclic(&self) -> bool {
        let mut colors = vec![VisitColor::White; self.number_of_vertices()];
        for vertex in 0..self.number_of_vertices() {
            if colors[vertex] == VisitColor::White
                && !visit_unoriented_acyclic(self, vertex, None, &mut colors)
            {
                return false;
            }
        }
        true
    }

    /// Test QPA's `IsTreeQuiver` condition.
    #[must_use]
    pub fn is_tree(&self) -> bool {
        self.number_of_vertices() > 0
            && self.number_of_arrows() == self.number_of_vertices() - 1
            && self.is_connected()
    }

    /// Return the type recognized by QPA's `IsDynkinQuiver` method.
    ///
    /// QPA prints this type as a side effect and returns a Boolean. Rust
    /// returns the same classification as structured data; `None` corresponds
    /// to QPA's `false`.
    #[must_use]
    pub fn dynkin_type(&self) -> Option<DynkinType> {
        self.dynkin_decision().dynkin_type()
    }

    /// Translate `IsDynkinQuiver`, including its printed side effects.
    #[must_use]
    pub fn dynkin_decision(&self) -> DynkinDecision {
        if !self.is_connected() {
            return DynkinDecision {
                dynkin_type: None,
                messages: vec![DYNKIN_NOT_CONNECTED_MESSAGE.into()],
            };
        }
        if self.number_of_vertices().checked_sub(1) != Some(self.number_of_arrows()) {
            return DynkinDecision {
                dynkin_type: None,
                messages: vec![DYNKIN_UNORIENTED_CYCLE_MESSAGE.into()],
            };
        }

        let dynkin_type = self.classify_dynkin_tree();
        let messages = dynkin_type
            .map(|value| vec![format!("{value}\n").into_boxed_str()])
            .unwrap_or_default();
        DynkinDecision {
            dynkin_type,
            messages,
        }
    }

    fn classify_dynkin_tree(&self) -> Option<DynkinType> {
        let high_degree = self
            .vertices()
            .filter(|vertex| vertex.in_degree() + vertex.out_degree() > 2)
            .collect::<Vec<_>>();
        if high_degree.is_empty() {
            return Some(DynkinType {
                family: DynkinFamily::A,
                rank: self.number_of_vertices(),
            });
        }
        if high_degree.len() != 1 || high_degree[0].in_degree() + high_degree[0].out_degree() != 3 {
            return None;
        }

        let fork = &high_degree[0];
        let mut incident = self.0.incoming[fork.index.zero_based()].to_vec();
        incident.extend_from_slice(&self.0.outgoing[fork.index.zero_based()]);
        let mut arms = incident
            .into_iter()
            .map(|arrow| self.dynkin_arm_length(fork.index, arrow))
            .collect::<Vec<_>>();
        arms.sort_unstable();

        let family = if arms.get(0..2) == Some([1, 1].as_slice()) {
            DynkinFamily::D
        } else if matches!(arms.as_slice(), [1, 2, 2..=4]) {
            DynkinFamily::E
        } else {
            return None;
        };
        Some(DynkinType {
            family,
            rank: self.number_of_vertices(),
        })
    }

    /// Translate QPA's `FullSubquiver(Q, vertices)`.
    ///
    /// Vertex order and duplicate names are preserved. Every arrow whose two
    /// endpoints occur in the selected list is retained in original order.
    ///
    /// # Errors
    ///
    /// Returns an error if a selected vertex belongs to a different quiver or
    /// if the QPA-compatible named constructor fails.
    pub fn full_subquiver(&self, vertices: &[Vertex]) -> Result<Self, QuiverError> {
        for vertex in vertices {
            self.ensure_same_quiver(vertex.quiver())?;
        }
        let selected = vertices.iter().map(Vertex::index).collect::<HashSet<_>>();
        let names = vertices
            .iter()
            .map(|vertex| vertex.name().into())
            .collect::<Vec<Box<str>>>();
        let arrows = self
            .0
            .arrows
            .iter()
            .filter(|arrow| selected.contains(&arrow.source) && selected.contains(&arrow.target))
            .map(|arrow| {
                NamedArrowSpec::named(
                    self.0.vertex_names[arrow.source.zero_based()].as_ref(),
                    self.0.vertex_names[arrow.target.zero_based()].as_ref(),
                    arrow.name.clone(),
                )
            })
            .collect();
        Self::from_named(names, arrows)
    }

    /// Translate QPA's `ConnectedComponentsOfQuiver` method.
    ///
    /// Component vertices appear in the DFS postorder used by QPA.
    ///
    /// # Errors
    ///
    /// Propagates construction failures from [`Self::full_subquiver`].
    pub fn connected_components(&self) -> Result<Vec<Self>, QuiverError> {
        let mut visited = vec![false; self.number_of_vertices()];
        let mut result = Vec::new();
        for start in 0..self.number_of_vertices() {
            if !visited[start] {
                let mut component = Vec::new();
                visit_connected(self, start, &mut visited, &mut component);
                let vertices = component
                    .into_iter()
                    .filter_map(|index| self.vertex(VertexIndex::from_zero_based(index)))
                    .collect::<Vec<_>>();
                result.push(self.full_subquiver(&vertices)?);
            }
        }
        Ok(result)
    }

    /// Translate QPA's generic `DynkinQuiver` dispatcher.
    ///
    /// # Errors
    ///
    /// Returns a stable validation error for an invalid rank, orientation
    /// count, or orientation symbol for the requested family.
    pub fn dynkin(
        family: DynkinFamily,
        rank: usize,
        orientation: &[DynkinOrientation],
    ) -> Result<Self, QuiverError> {
        match family {
            DynkinFamily::A => Self::dynkin_a(rank, orientation),
            DynkinFamily::D => Self::dynkin_d(rank, orientation),
            DynkinFamily::E => Self::dynkin_e(rank, orientation),
        }
    }

    /// Translate QPA's string-dispatched `DynkinQuiver` including its
    /// implicit `fail` for a family other than `"A"`, `"D"`, or `"E"`.
    ///
    /// # Errors
    ///
    /// Propagates the selected Dynkin constructor's validation errors.
    pub fn dynkin_from_name(
        family: &str,
        rank: usize,
        orientation: &[DynkinOrientation],
    ) -> Result<Option<Self>, QuiverError> {
        DynkinFamily::from_qpa_name(family)
            .map(|family| Self::dynkin(family, rank, orientation))
            .transpose()
    }

    /// Translate QPA's `DynkinQuiverAn`.
    ///
    /// # Errors
    ///
    /// Requires `rank >= 1`, exactly `rank - 1` left/right symbols.
    pub fn dynkin_a(rank: usize, orientation: &[DynkinOrientation]) -> Result<Self, QuiverError> {
        if rank == 0 {
            return Err(QuiverError::InvalidDynkinRank {
                family: DynkinFamily::A,
                rank,
            });
        }
        validate_orientation_count(DynkinFamily::A, rank - 1, orientation.len())?;
        validate_horizontal_orientations(DynkinFamily::A, orientation)?;

        let names = dynkin_vertex_names(rank);
        let arrows = (0..rank - 1)
            .map(|position| {
                dynkin_line_arrow(&names, position, orientation[position], position + 1)
            })
            .collect();
        Self::from_named(names, arrows)
    }

    /// Translate QPA's `DynkinQuiverDn`.
    ///
    /// # Errors
    ///
    /// Requires `rank >= 4`, exactly `rank - 1` left/right symbols.
    pub fn dynkin_d(rank: usize, orientation: &[DynkinOrientation]) -> Result<Self, QuiverError> {
        if rank < 4 {
            return Err(QuiverError::InvalidDynkinRank {
                family: DynkinFamily::D,
                rank,
            });
        }
        validate_orientation_count(DynkinFamily::D, rank - 1, orientation.len())?;
        validate_horizontal_orientations(DynkinFamily::D, orientation)?;

        let names = dynkin_vertex_names(rank);
        let mut arrows = Vec::with_capacity(rank - 1);
        for (position, &symbol) in orientation.iter().enumerate().take(2) {
            let (source, target) = if symbol == DynkinOrientation::Left {
                (2, position)
            } else {
                (position, 2)
            };
            arrows.push(NamedArrowSpec::named(
                names[source].as_ref(),
                names[target].as_ref(),
                format!("a_{}", position + 1),
            ));
        }
        for (position, &symbol) in orientation.iter().enumerate().take(rank - 1).skip(2) {
            arrows.push(dynkin_line_arrow(&names, position, symbol, position + 1));
        }
        Self::from_named(names, arrows)
    }

    /// Translate QPA's `DynkinQuiverEn`.
    ///
    /// # Errors
    ///
    /// Requires rank 6, 7, or 8; the first `rank - 2` symbols must be
    /// left/right and the final symbol up/down.
    pub fn dynkin_e(rank: usize, orientation: &[DynkinOrientation]) -> Result<Self, QuiverError> {
        if !matches!(rank, 6..=8) {
            return Err(QuiverError::InvalidDynkinRank {
                family: DynkinFamily::E,
                rank,
            });
        }
        validate_orientation_count(DynkinFamily::E, rank - 1, orientation.len())?;
        validate_horizontal_orientations(DynkinFamily::E, &orientation[..rank - 2])?;
        if !matches!(
            orientation[rank - 2],
            DynkinOrientation::Up | DynkinOrientation::Down
        ) {
            return Err(QuiverError::InvalidDynkinOrientation {
                family: DynkinFamily::E,
                position: rank - 1,
            });
        }

        let names = dynkin_vertex_names(rank);
        let mut arrows = (0..rank - 2)
            .map(|position| {
                dynkin_line_arrow(&names, position, orientation[position], position + 1)
            })
            .collect::<Vec<_>>();
        let (source, target) = if orientation[rank - 2] == DynkinOrientation::Down {
            (rank - 1, 2)
        } else {
            (2, rank - 1)
        };
        arrows.push(NamedArrowSpec::named(
            names[source].as_ref(),
            names[target].as_ref(),
            format!("a_{}", rank - 1),
        ));
        Self::from_named(names, arrows)
    }

    fn dynkin_arm_length(&self, fork: VertexIndex, first_arrow: ArrowIndex) -> usize {
        let mut vertex = fork;
        let mut arrow = first_arrow;
        let mut length = 0;
        loop {
            let data = &self.0.arrows[arrow.zero_based()];
            vertex = if data.target == vertex {
                data.source
            } else {
                data.target
            };
            length += 1;

            let mut remaining = self.0.incoming[vertex.zero_based()].to_vec();
            remaining.extend_from_slice(&self.0.outgoing[vertex.zero_based()]);
            remaining.retain(|candidate| *candidate != arrow);
            if remaining.len() != 1 {
                return length;
            }
            arrow = remaining[0];
        }
    }

    /// Construct and cache QPA's opposite quiver.
    ///
    /// The cache preserves `OppositeQuiver(OppositeQuiver(Q))` object identity
    /// without creating a strong reference cycle.
    ///
    /// # Errors
    ///
    /// Returns a construction error only if internal QPA-derived specifications
    /// cannot be represented.
    pub fn opposite(&self) -> Result<Self, QuiverError> {
        if let Some(opposite) = self.0.opposite_strong.get() {
            return Ok(opposite.clone());
        }
        if let Some(opposite) = self
            .0
            .opposite_weak
            .read()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .as_ref()
            .and_then(Weak::upgrade)
        {
            return Ok(Self(opposite));
        }

        let vertices = self
            .0
            .vertex_names
            .iter()
            .map(|name| format!("{name}_op").into_boxed_str())
            .collect::<Vec<_>>();
        let arrows = self
            .0
            .arrows
            .iter()
            .map(|arrow| {
                NamedArrowSpec::named(
                    &*vertices[arrow.target.zero_based()],
                    &*vertices[arrow.source.zero_based()],
                    format!("{}_op", arrow.name),
                )
            })
            .collect();
        let opposite = Self::from_named(vertices, arrows)?;
        if opposite.0.opposite_strong.set(self.clone()).is_err() {
            return Err(QuiverError::OppositeCacheRace);
        }
        if opposite
            .0
            .opposite_name_map
            .set(OppositeNameMap::RemoveSuffix)
            .is_err()
        {
            return Err(QuiverError::OppositeCacheRace);
        }

        let _ = self.0.opposite_name_map.set(OppositeNameMap::AppendSuffix);
        let mut slot = self
            .0
            .opposite_weak
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if let Some(cached) = slot.as_ref().and_then(Weak::upgrade) {
            return Ok(Self(cached));
        }
        *slot = Some(Arc::downgrade(&opposite.0));
        Ok(opposite)
    }

    /// Apply QPA's `OppositeQuiverNameMap(Q)` to one generator name.
    ///
    /// Constructing the opposite first installs the same append/remove map
    /// attributes as QPA's `OppositeQuiver` method.
    ///
    /// # Errors
    ///
    /// Propagates opposite-quiver construction. On an opposite quiver QPA's
    /// `deop_name` closure removes the final three characters unconditionally;
    /// it does not validate that they spell `_op`.
    pub fn opposite_quiver_name(&self, name: &str) -> Result<Box<str>, PathError> {
        let _ = self.opposite().map_err(PathError::QuiverConstruction)?;
        let name_map = *self
            .0
            .opposite_name_map
            .get()
            .ok_or(PathError::MissingOppositeNameMap)?;
        Ok(match name_map {
            OppositeNameMap::AppendSuffix => format!("{name}_op").into_boxed_str(),
            OppositeNameMap::RemoveSuffix => {
                let mut boundaries = name
                    .char_indices()
                    .map(|(index, _)| index)
                    .collect::<Vec<_>>();
                boundaries.push(name.len());
                let keep = boundaries
                    .len()
                    .checked_sub(4)
                    .ok_or_else(|| PathError::OppositeNameTooShort { name: name.into() })?;
                name[..boundaries[keep]].into()
            }
        })
    }

    /// Translate QPA's `OppositePath` method.
    ///
    /// QPA performs this operation through generator names, including its
    /// record-component collision behavior, rather than directly through
    /// indices. This implementation follows that route.
    ///
    /// # Errors
    ///
    /// Returns an error for a foreign or zero path, or if name collisions make
    /// the mapped generator product invalid.
    pub fn opposite_path(&self, path: &Path) -> Result<Path, PathError> {
        path.ensure_quiver(self)?;
        if path.is_zero() {
            return Err(PathError::ZeroHasNoWalk);
        }
        let opposite = self.opposite().map_err(PathError::QuiverConstruction)?;
        let names = match &path.representation {
            PathRepresentation::Vertex(index) => {
                vec![self.0.vertex_names[index.zero_based()].as_ref()]
            }
            PathRepresentation::Walk(walk) => walk
                .iter()
                .rev()
                .map(|index| self.0.arrows[index.zero_based()].name.as_ref())
                .collect(),
            PathRepresentation::FormalGeneratorWalk(walk) => walk
                .iter()
                .rev()
                .map(|position| {
                    if *position <= self.number_of_vertices() {
                        self.0.vertex_names[*position - 1].as_ref()
                    } else {
                        self.0.arrows[*position - self.number_of_vertices() - 1]
                            .name
                            .as_ref()
                    }
                })
                .collect(),
            PathRepresentation::Zero => unreachable!(),
        };
        let mut mapped = Vec::with_capacity(names.len());
        for name in names {
            let mapped_name = self.opposite_quiver_name(name)?;
            mapped.push(
                opposite
                    .generator_by_name(&mapped_name)
                    .ok_or(PathError::MissingOppositeGenerator { name: mapped_name })?,
            );
        }
        Path::product_of_generators(&opposite, &mapped)
    }

    /// Translate QPA's `SeparatedQuiver` construction.
    ///
    /// # Errors
    ///
    /// Propagates construction errors caused by QPA-compatible name lookup.
    pub fn separated(&self) -> Result<Self, QuiverError> {
        let mut vertices = self.0.vertex_names.to_vec();
        vertices.extend(
            self.0
                .vertex_names
                .iter()
                .map(|name| format!("{name}'").into_boxed_str()),
        );
        let arrows = self
            .0
            .arrows
            .iter()
            .map(|arrow| {
                NamedArrowSpec::named(
                    &*self.0.vertex_names[arrow.source.zero_based()],
                    format!("{}'", self.0.vertex_names[arrow.target.zero_based()]),
                    arrow.name.clone(),
                )
            })
            .collect();
        Self::from_named(vertices, arrows)
    }

    /// Translate QPA's `DoubleQuiver` construction.
    ///
    /// # Errors
    ///
    /// Propagates construction errors caused by QPA-compatible name lookup.
    pub fn double(&self) -> Result<Self, QuiverError> {
        let mut arrows = Vec::with_capacity(self.number_of_arrows().saturating_mul(2));
        for arrow in &self.0.arrows {
            let source = &self.0.vertex_names[arrow.source.zero_based()];
            let target = &self.0.vertex_names[arrow.target.zero_based()];
            arrows.push(NamedArrowSpec::named(
                source.as_ref(),
                target.as_ref(),
                arrow.name.clone(),
            ));
            arrows.push(NamedArrowSpec::named(
                target.as_ref(),
                source.as_ref(),
                format!("{}bar", arrow.name),
            ));
        }
        Self::from_named(self.0.vertex_names.to_vec(), arrows)
    }

    fn ensure_same_quiver(&self, other: &Self) -> Result<(), QuiverError> {
        if self.is_identical(other) {
            Ok(())
        } else {
            Err(QuiverError::DifferentQuivers)
        }
    }
}

impl fmt::Debug for Quiver {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Quiver")
            .field("id", &self.0.id)
            .field("vertices", &self.0.vertex_names)
            .field("arrows", &self.0.arrows)
            .finish_non_exhaustive()
    }
}

impl fmt::Display for Quiver {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("Quiver( [")?;
        for (position, vertex) in self.vertices().enumerate() {
            if position > 0 {
                formatter.write_str(",")?;
            }
            write!(formatter, "\"{}\"", vertex.name())?;
        }
        formatter.write_str("], [")?;
        for (position, arrow) in self.arrows().enumerate() {
            if position > 0 {
                formatter.write_str(",")?;
            }
            write!(
                formatter,
                "[\"{}\",\"{}\",\"{}\"]",
                arrow.source().name(),
                arrow.target().name(),
                arrow.name()
            )?;
        }
        formatter.write_str("] )")
    }
}

/// The fresh quiver family and retained ordering produced by QPA's
/// `OrderedBy(Q, O)` method.
#[derive(Clone, Debug)]
pub struct OrderedQuiver {
    quiver: Quiver,
    ordering: QuiverOrdering,
}

impl OrderedQuiver {
    /// Return the newly allocated quiver family.
    #[must_use]
    pub const fn quiver(&self) -> &Quiver {
        &self.quiver
    }

    /// Return the identical ordering handle supplied to `OrderedBy`.
    #[must_use]
    pub const fn ordering(&self) -> &QuiverOrdering {
        &self.ordering
    }

    /// Compare paths through the retained source ordering.
    ///
    /// # Errors
    ///
    /// Both paths must belong to the new family, and its generator count must
    /// match that of the quiver for which the ordering was constructed.
    pub fn compare_paths(&self, left: &Path, right: &Path) -> Result<Ordering, OrderingError> {
        self.ordering
            .compare_on_equally_generated_quiver(&self.quiver, left, right)
    }

    /// Consume the wrapper into a quiver and an ordering rebound to it.
    ///
    /// This is the typed Rust route corresponding to passing QPA's ordered
    /// quiver into `PathAlgebra`.
    ///
    /// # Errors
    ///
    /// Returns an error if the supplied source ordering has a different
    /// generator count.
    pub fn into_parts(self) -> Result<(Quiver, QuiverOrdering), OrderingError> {
        let ordering = self.ordering.rebase_to(&self.quiver)?;
        Ok((self.quiver, ordering))
    }
}

impl std::ops::Deref for OrderedQuiver {
    type Target = Quiver;

    fn deref(&self) -> &Self::Target {
        &self.quiver
    }
}

impl fmt::Display for OrderedQuiver {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.quiver.fmt(formatter)
    }
}

impl PartialEq for Quiver {
    fn eq(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}

impl Eq for Quiver {}

impl Hash for Quiver {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.id.hash(state);
    }
}

/// A quiver vertex handle.
#[derive(Clone, Debug)]
pub struct Vertex {
    quiver: Quiver,
    index: VertexIndex,
}

impl Vertex {
    /// Return the containing quiver.
    #[must_use]
    pub const fn quiver(&self) -> &Quiver {
        &self.quiver
    }

    /// Return QPA's one-based vertex and generator index.
    #[must_use]
    pub const fn index(&self) -> VertexIndex {
        self.index
    }

    /// Return the vertex name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.quiver.0.vertex_names[self.index.zero_based()]
    }

    /// Regard this vertex as a length-zero nonzero path.
    #[must_use]
    pub fn path(&self) -> Path {
        Path {
            quiver: self.quiver.clone(),
            representation: PathRepresentation::Vertex(self.index),
        }
    }

    /// Return the in-degree, counting parallel arrows.
    #[must_use]
    pub fn in_degree(&self) -> usize {
        self.quiver.0.incoming[self.index.zero_based()].len()
    }

    /// Return the out-degree, counting parallel arrows.
    #[must_use]
    pub fn out_degree(&self) -> usize {
        self.quiver.0.outgoing[self.index.zero_based()].len()
    }
}

impl fmt::Display for Vertex {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

impl PartialEq for Vertex {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index && self.quiver == other.quiver
    }
}

impl Eq for Vertex {}

impl Hash for Vertex {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.quiver.hash(state);
        self.index.hash(state);
    }
}

/// A quiver arrow handle.
#[derive(Clone, Debug)]
pub struct Arrow {
    quiver: Quiver,
    index: ArrowIndex,
}

impl Arrow {
    /// Return the containing quiver.
    #[must_use]
    pub const fn quiver(&self) -> &Quiver {
        &self.quiver
    }

    /// Return QPA's one-based index among arrows.
    #[must_use]
    pub const fn index(&self) -> ArrowIndex {
        self.index
    }

    /// Return QPA's one-based generator index (vertices first, then arrows).
    #[must_use]
    pub fn generator_index(&self) -> usize {
        self.quiver.number_of_vertices() + self.index.get()
    }

    /// Return the arrow name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.quiver.0.arrows[self.index.zero_based()].name
    }

    /// Return the source vertex.
    #[must_use]
    pub fn source(&self) -> Vertex {
        let index = self.quiver.0.arrows[self.index.zero_based()].source;
        Vertex {
            quiver: self.quiver.clone(),
            index,
        }
    }

    /// Return the target vertex.
    #[must_use]
    pub fn target(&self) -> Vertex {
        let index = self.quiver.0.arrows[self.index.zero_based()].target;
        Vertex {
            quiver: self.quiver.clone(),
            index,
        }
    }

    /// Regard this arrow as a length-one path.
    #[must_use]
    pub fn path(&self) -> Path {
        Path {
            quiver: self.quiver.clone(),
            representation: PathRepresentation::Walk(vec![self.index].into_boxed_slice()),
        }
    }
}

impl fmt::Display for Arrow {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

impl PartialEq for Arrow {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index && self.quiver == other.quiver
    }
}

impl Eq for Arrow {}

impl Hash for Arrow {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.quiver.hash(state);
        self.index.hash(state);
    }
}

/// A vertex or arrow in QPA generator order.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum Generator {
    /// A length-zero vertex path.
    Vertex(Vertex),
    /// A length-one arrow path.
    Arrow(Arrow),
}

impl Generator {
    /// Return the generator name.
    #[must_use]
    pub fn name(&self) -> &str {
        match self {
            Self::Vertex(vertex) => vertex.name(),
            Self::Arrow(arrow) => arrow.name(),
        }
    }

    /// Convert this generator to a path.
    #[must_use]
    pub fn path(&self) -> Path {
        match self {
            Self::Vertex(vertex) => vertex.path(),
            Self::Arrow(arrow) => arrow.path(),
        }
    }
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
enum PathRepresentation {
    Zero,
    Vertex(VertexIndex),
    Walk(Box<[ArrowIndex]>),
    /// QPA's `PathFromArrowList` also accepts a homogeneous list containing
    /// vertices.  Such an object is a formal path: its walk is the supplied
    /// generator list and need not be composable.
    FormalGeneratorWalk(Box<[usize]>),
}

/// A zero path, vertex, arrow, or composite path within one quiver family.
#[derive(Clone, Debug)]
pub struct Path {
    quiver: Quiver,
    representation: PathRepresentation,
}

impl Path {
    /// Return the containing quiver.
    #[must_use]
    pub const fn quiver(&self) -> &Quiver {
        &self.quiver
    }

    /// Return whether this is QPA's distinguished zero path.
    #[must_use]
    pub const fn is_zero(&self) -> bool {
        matches!(self.representation, PathRepresentation::Zero)
    }

    /// Return whether this path is a quiver vertex.
    #[must_use]
    pub const fn is_vertex(&self) -> bool {
        matches!(self.representation, PathRepresentation::Vertex(_))
    }

    /// Return whether this path is an arrow generator.
    #[must_use]
    pub fn is_arrow(&self) -> bool {
        matches!(&self.representation, PathRepresentation::Walk(walk) if walk.len() == 1)
    }

    /// Return `LengthOfPath`, which is unbound for QPA's zero path.
    #[must_use]
    pub fn length(&self) -> Option<usize> {
        match &self.representation {
            PathRepresentation::Zero => None,
            PathRepresentation::Vertex(_) => Some(0),
            PathRepresentation::Walk(walk) => Some(walk.len()),
            PathRepresentation::FormalGeneratorWalk(walk) => Some(walk.len()),
        }
    }

    /// Return arrow indices in `WalkOfPath` order.
    ///
    /// QPA leaves `WalkOfPath` unbound on its zero path and returns an empty
    /// walk for vertices. A formal vertex-containing `PathFromArrowList` walk
    /// is available through [`Self::generator_walk`] instead.
    #[must_use]
    pub fn walk(&self) -> Option<&[ArrowIndex]> {
        match &self.representation {
            PathRepresentation::Zero | PathRepresentation::FormalGeneratorWalk(_) => None,
            PathRepresentation::Vertex(_) => Some(&[]),
            PathRepresentation::Walk(walk) => Some(walk),
        }
    }

    /// Return QPA's complete `WalkOfPath` as quiver generators.
    ///
    /// Unlike [`Self::walk`], this represents the formal vertex-containing
    /// walks constructed by `PathFromArrowList`.  The zero path has no bound
    /// walk and vertices have the empty walk, exactly as in QPA.
    #[must_use]
    pub fn generator_walk(&self) -> Option<Vec<Generator>> {
        match &self.representation {
            PathRepresentation::Zero => None,
            PathRepresentation::Vertex(_) => Some(Vec::new()),
            PathRepresentation::Walk(walk) => Some(
                walk.iter()
                    .filter_map(|index| self.quiver.arrow(*index))
                    .map(Generator::Arrow)
                    .collect(),
            ),
            PathRepresentation::FormalGeneratorWalk(walk) => Some(
                walk.iter()
                    .map(|&position| self.generator_at_position_unchecked(position))
                    .collect(),
            ),
        }
    }

    /// Return the source endpoint. For QPA's zero path this is the zero path
    /// itself rather than a vertex.
    #[must_use]
    pub fn source(&self) -> PathEndpoint {
        match &self.representation {
            PathRepresentation::Zero => PathEndpoint::Zero(self.quiver.zero()),
            PathRepresentation::Vertex(index) => PathEndpoint::Vertex(Vertex {
                quiver: self.quiver.clone(),
                index: *index,
            }),
            PathRepresentation::Walk(walk) => {
                let Some(arrow) = walk
                    .first()
                    .and_then(|index| self.quiver.0.arrows.get(index.zero_based()))
                else {
                    return PathEndpoint::Zero(self.quiver.zero());
                };
                PathEndpoint::Vertex(Vertex {
                    quiver: self.quiver.clone(),
                    index: arrow.source,
                })
            }
            PathRepresentation::FormalGeneratorWalk(walk) => PathEndpoint::Vertex(Vertex {
                quiver: self.quiver.clone(),
                index: self.generator_source_vertex_index(walk[0]),
            }),
        }
    }

    /// Return the target endpoint. For QPA's zero path this is the zero path
    /// itself rather than a vertex.
    #[must_use]
    pub fn target(&self) -> PathEndpoint {
        match &self.representation {
            PathRepresentation::Zero => PathEndpoint::Zero(self.quiver.zero()),
            PathRepresentation::Vertex(index) => PathEndpoint::Vertex(Vertex {
                quiver: self.quiver.clone(),
                index: *index,
            }),
            PathRepresentation::Walk(walk) => {
                let Some(arrow) = walk
                    .last()
                    .and_then(|index| self.quiver.0.arrows.get(index.zero_based()))
                else {
                    return PathEndpoint::Zero(self.quiver.zero());
                };
                PathEndpoint::Vertex(Vertex {
                    quiver: self.quiver.clone(),
                    index: arrow.target,
                })
            }
            PathRepresentation::FormalGeneratorWalk(walk) => PathEndpoint::Vertex(Vertex {
                quiver: self.quiver.clone(),
                index: self.generator_target_vertex_index(walk[walk.len() - 1]),
            }),
        }
    }

    /// Return QPA's `ExtRepOfObj` list.
    #[must_use]
    pub fn external_representation(&self) -> Vec<usize> {
        match &self.representation {
            PathRepresentation::Zero => vec![0],
            PathRepresentation::Vertex(index) => vec![index.get()],
            PathRepresentation::Walk(walk) => walk
                .iter()
                .map(|index| self.quiver.number_of_vertices() + index.get())
                .collect(),
            PathRepresentation::FormalGeneratorWalk(walk) => walk.to_vec(),
        }
    }

    /// Translate QPA's `ObjByExtRep` method.
    ///
    /// # Errors
    ///
    /// QPA's generic empty product returns `fail`, represented by `None`.
    /// Invalid generator positions or a product between different generator
    /// categories are returned as explicit errors.
    pub fn from_external_representation(
        quiver: &Quiver,
        representation: &[usize],
    ) -> Result<Option<Self>, PathError> {
        if representation.is_empty() {
            return Ok(None);
        }
        let mut generators = Vec::with_capacity(representation.len());
        let mut has_zero = false;
        for &position in representation {
            if position == 0 {
                has_zero = true;
                continue;
            }
            if position <= quiver.number_of_vertices() {
                let index =
                    VertexIndex::new(position).ok_or(PathError::GeneratorIndexOutOfRange {
                        index: position,
                        number_of_generators: quiver.number_of_vertices()
                            + quiver.number_of_arrows(),
                    })?;
                generators.push(Generator::Vertex(Vertex {
                    quiver: quiver.clone(),
                    index,
                }));
            } else {
                let arrow_position = position - quiver.number_of_vertices();
                let Some(index) = ArrowIndex::new(arrow_position) else {
                    unreachable!();
                };
                generators.push(Generator::Arrow(quiver.arrow(index).ok_or(
                    PathError::GeneratorIndexOutOfRange {
                        index: position,
                        number_of_generators: quiver.number_of_vertices()
                            + quiver.number_of_arrows(),
                    },
                )?));
            }
        }
        if has_zero {
            return Ok(Some(quiver.zero()));
        }
        Self::product_of_generators(quiver, &generators).map(Some)
    }

    /// Translate the list-of-arrows branch of QPA's public `Path` function.
    ///
    /// This branch returns GAP's generic `Product(arrow_list)`: multiplication
    /// checks composability and an empty list has no quiver identity.
    ///
    /// # Errors
    ///
    /// Returns an error for an arrow from another quiver. QPA's `fail` for an
    /// empty generic product is represented by `None`.
    pub fn from_arrow_product(
        quiver: &Quiver,
        arrows: &[Arrow],
    ) -> Result<Option<Self>, PathError> {
        if arrows.is_empty() {
            return Ok(None);
        }
        let generators = arrows
            .iter()
            .cloned()
            .map(Generator::Arrow)
            .collect::<Vec<_>>();
        Self::product_of_generators(quiver, &generators).map(Some)
    }

    /// Translate QPA's `PathFromArrowList` helper exactly.
    ///
    /// In particular, QPA returns zero for an empty list and does not validate
    /// composability for a list of two or more arrows.
    ///
    /// # Errors
    ///
    /// Returns an error if an arrow belongs to another quiver.
    pub fn from_arrow_list(quiver: &Quiver, arrows: &[Arrow]) -> Result<Self, PathError> {
        for arrow in arrows {
            if !arrow.quiver.is_identical(quiver) {
                return Err(PathError::DifferentQuivers);
            }
        }
        Ok(match arrows {
            [] => quiver.zero(),
            [arrow] => arrow.path(),
            _ => Self {
                quiver: quiver.clone(),
                representation: PathRepresentation::Walk(
                    arrows
                        .iter()
                        .map(Arrow::index)
                        .collect::<Vec<_>>()
                        .into_boxed_slice(),
                ),
            },
        })
    }

    /// Translate QPA's full `PathFromArrowList` surface.
    ///
    /// The name in QPA is historical: a singleton vertex is returned as that
    /// vertex, and a longer homogeneous generator list becomes a formal path
    /// without a composability check.  Empty input returns the zero path.
    ///
    /// # Errors
    ///
    /// Rust requires every supplied generator to belong to `quiver`; QPA can
    /// construct malformed cross-family formal paths that fail in later path
    /// operations.
    pub fn from_generator_list(
        quiver: &Quiver,
        generators: &[Generator],
    ) -> Result<Self, PathError> {
        for generator in generators {
            if !generator.path().quiver().is_identical(quiver) {
                return Err(PathError::DifferentQuivers);
            }
        }
        Ok(match generators {
            [] => quiver.zero(),
            [generator] => generator.path(),
            _ if generators
                .iter()
                .all(|generator| matches!(generator, Generator::Arrow(_))) =>
            {
                Self {
                    quiver: quiver.clone(),
                    representation: PathRepresentation::Walk(
                        generators
                            .iter()
                            .map(|generator| match generator {
                                Generator::Arrow(arrow) => arrow.index(),
                                Generator::Vertex(_) => unreachable!(),
                            })
                            .collect::<Vec<_>>()
                            .into_boxed_slice(),
                    ),
                }
            }
            _ => Self {
                quiver: quiver.clone(),
                representation: PathRepresentation::FormalGeneratorWalk(
                    generators
                        .iter()
                        .map(|generator| match generator {
                            Generator::Vertex(vertex) => vertex.index().get(),
                            Generator::Arrow(arrow) => {
                                quiver.number_of_vertices() + arrow.index().get()
                            }
                        })
                        .collect::<Vec<_>>()
                        .into_boxed_slice(),
                ),
            },
        })
    }

    /// Multiply paths in QPA's left-to-right convention: `a * b` means first
    /// `a`, then `b`.
    ///
    /// # Errors
    ///
    /// QPA installs multiplication only for paths in an identical family;
    /// foreign-quiver operands therefore return an explicit error.
    pub fn multiply(&self, right: &Self) -> Result<Self, PathError> {
        if !self.quiver.is_identical(&right.quiver) {
            return Err(PathError::DifferentQuivers);
        }
        if self.is_zero() || right.is_zero() {
            return Ok(self.quiver.zero());
        }

        match (&self.representation, &right.representation) {
            (PathRepresentation::Vertex(left), PathRepresentation::Vertex(right)) => {
                Ok(if left == right {
                    self.clone()
                } else {
                    self.quiver.zero()
                })
            }
            (PathRepresentation::Vertex(left), _) => {
                Ok(if right.source_vertex_index() == Some(*left) {
                    right.clone()
                } else {
                    self.quiver.zero()
                })
            }
            (_, PathRepresentation::Vertex(right)) => {
                Ok(if self.target_vertex_index() == Some(*right) {
                    self.clone()
                } else {
                    self.quiver.zero()
                })
            }
            (
                PathRepresentation::Walk(_) | PathRepresentation::FormalGeneratorWalk(_),
                PathRepresentation::Walk(_) | PathRepresentation::FormalGeneratorWalk(_),
            ) => {
                if self.target_vertex_index() != right.source_vertex_index() {
                    return Ok(self.quiver.zero());
                }
                let left = self.generator_positions();
                let right = right.generator_positions();
                let mut walk = Vec::with_capacity(left.len().saturating_add(right.len()));
                walk.extend(left);
                walk.extend(right);
                Ok(Self::from_generator_positions_unchecked(&self.quiver, walk))
            }
            (PathRepresentation::Zero, _) | (_, PathRepresentation::Zero) => unreachable!(),
        }
    }

    fn product_of_generators(quiver: &Quiver, generators: &[Generator]) -> Result<Self, PathError> {
        let Some((first, rest)) = generators.split_first() else {
            return Err(PathError::EmptyProduct);
        };
        let mut product = first.path();
        product.ensure_quiver(quiver)?;
        for generator in rest {
            product = product.multiply(&generator.path())?;
        }
        Ok(product)
    }

    fn ordering_length(&self) -> usize {
        self.length().unwrap_or(0)
    }

    fn source_vertex_index(&self) -> Option<VertexIndex> {
        match &self.representation {
            PathRepresentation::Zero => None,
            PathRepresentation::Vertex(index) => Some(*index),
            PathRepresentation::Walk(walk) => self
                .quiver
                .arrow(walk[0])
                .map(|arrow| arrow.source().index()),
            PathRepresentation::FormalGeneratorWalk(walk) => {
                Some(self.generator_source_vertex_index(walk[0]))
            }
        }
    }

    fn target_vertex_index(&self) -> Option<VertexIndex> {
        match &self.representation {
            PathRepresentation::Zero => None,
            PathRepresentation::Vertex(index) => Some(*index),
            PathRepresentation::Walk(walk) => self
                .quiver
                .arrow(*walk.last().expect("walk paths are nonempty"))
                .map(|arrow| arrow.target().index()),
            PathRepresentation::FormalGeneratorWalk(walk) => {
                Some(self.generator_target_vertex_index(walk[walk.len() - 1]))
            }
        }
    }

    fn generator_positions(&self) -> Vec<usize> {
        match &self.representation {
            PathRepresentation::Walk(walk) => walk
                .iter()
                .map(|index| self.quiver.number_of_vertices() + index.get())
                .collect(),
            PathRepresentation::FormalGeneratorWalk(walk) => walk.to_vec(),
            PathRepresentation::Zero | PathRepresentation::Vertex(_) => Vec::new(),
        }
    }

    fn from_generator_positions_unchecked(quiver: &Quiver, positions: Vec<usize>) -> Self {
        if positions
            .iter()
            .all(|position| *position > quiver.number_of_vertices())
        {
            Self {
                quiver: quiver.clone(),
                representation: PathRepresentation::Walk(
                    positions
                        .into_iter()
                        .map(|position| {
                            ArrowIndex::new(position - quiver.number_of_vertices())
                                .expect("arrow generator positions are positive")
                        })
                        .collect::<Vec<_>>()
                        .into_boxed_slice(),
                ),
            }
        } else {
            Self {
                quiver: quiver.clone(),
                representation: PathRepresentation::FormalGeneratorWalk(
                    positions.into_boxed_slice(),
                ),
            }
        }
    }

    fn generator_at_position_unchecked(&self, position: usize) -> Generator {
        if position <= self.quiver.number_of_vertices() {
            VertexIndex::new(position)
                .and_then(|index| self.quiver.vertex(index))
                .map(Generator::Vertex)
                .expect("formal vertex positions belong to the path quiver")
        } else {
            ArrowIndex::new(position - self.quiver.number_of_vertices())
                .and_then(|index| self.quiver.arrow(index))
                .map(Generator::Arrow)
                .expect("formal arrow positions belong to the path quiver")
        }
    }

    fn generator_source_vertex_index(&self, position: usize) -> VertexIndex {
        match self.generator_at_position_unchecked(position) {
            Generator::Vertex(vertex) => vertex.index(),
            Generator::Arrow(arrow) => arrow.source().index(),
        }
    }

    fn generator_target_vertex_index(&self, position: usize) -> VertexIndex {
        match self.generator_at_position_unchecked(position) {
            Generator::Vertex(vertex) => vertex.index(),
            Generator::Arrow(arrow) => arrow.target().index(),
        }
    }

    fn ensure_quiver(&self, quiver: &Quiver) -> Result<(), PathError> {
        if self.quiver.is_identical(quiver) {
            Ok(())
        } else {
            Err(PathError::DifferentQuivers)
        }
    }
}

impl fmt::Display for Path {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.representation {
            PathRepresentation::Zero => formatter.write_str("0"),
            PathRepresentation::Vertex(index) => {
                formatter.write_str(&self.quiver.0.vertex_names[index.zero_based()])
            }
            PathRepresentation::Walk(walk) => {
                let mut groups = walk.iter().peekable();
                let mut first = true;
                while let Some(index) = groups.next() {
                    if !first {
                        formatter.write_str("*")?;
                    }
                    first = false;
                    formatter.write_str(&self.quiver.0.arrows[index.zero_based()].name)?;
                    let mut exponent = 1_usize;
                    while groups.peek().is_some_and(|next| *next == index) {
                        groups.next();
                        exponent += 1;
                    }
                    if exponent > 1 {
                        write!(formatter, "^{exponent}")?;
                    }
                }
                Ok(())
            }
            PathRepresentation::FormalGeneratorWalk(walk) => {
                let mut groups = walk.iter().peekable();
                let mut first = true;
                while let Some(&position) = groups.next() {
                    if !first {
                        formatter.write_str("*")?;
                    }
                    first = false;
                    match self.generator_at_position_unchecked(position) {
                        Generator::Vertex(vertex) => formatter.write_str(vertex.name())?,
                        Generator::Arrow(arrow) => formatter.write_str(arrow.name())?,
                    }
                    let mut exponent = 1_usize;
                    while groups.peek().is_some_and(|next| **next == position) {
                        groups.next();
                        exponent += 1;
                    }
                    if exponent > 1 {
                        write!(formatter, "^{exponent}")?;
                    }
                }
                Ok(())
            }
        }
    }
}

impl PartialEq for Path {
    fn eq(&self, other: &Self) -> bool {
        self.quiver == other.quiver && self.representation == other.representation
    }
}

impl Eq for Path {}

impl Hash for Path {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.quiver.hash(state);
        self.representation.hash(state);
    }
}

/// The endpoint value returned by QPA's `SourceOfPath` and `TargetOfPath`.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum PathEndpoint {
    /// A nonzero path endpoint.
    Vertex(Vertex),
    /// The zero path is its own source and target in QPA.
    Zero(Path),
}

/// QPA's breadth-first path iterator state.
#[derive(Clone, Debug)]
pub struct PathIterator {
    quiver: Quiver,
    position: NextPathPosition,
}

impl PathIterator {
    /// Return QPA's `IsDoneIterator` result.
    ///
    /// In particular, the untouched iterator of an empty quiver is not done:
    /// the pinned QPA `NextIterator` returns `fail` without replacing its
    /// initial integer-zero position.
    #[must_use]
    pub fn is_done(&self) -> bool {
        matches!(&self.position, NextPathPosition::Queue(queue) if queue.is_empty())
    }

    /// Return the current QPA `NextPath` position.
    #[must_use]
    pub const fn position(&self) -> &NextPathPosition {
        &self.position
    }
}

impl Iterator for PathIterator {
    type Item = Path;

    fn next(&mut self) -> Option<Self::Item> {
        let (result, position) = self
            .quiver
            .next_path(&self.position)
            .expect("an internally produced NextPath queue is valid")?;
        self.position = position;
        Some(result)
    }
}

/// Failure to construct a QPA quiver.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum QuiverError {
    /// QPA's indexed constructor requires a positive integer.
    #[error("Quiver(N, arrows) requires at least one vertex")]
    IndexedQuiverNeedsVertex,
    /// GAP arrow specification rows must all have the same length.
    #[error("arrow specifications cannot mix explicit and default names")]
    MixedArrowNaming,
    /// A one-based endpoint lies outside the quiver.
    #[error("vertex index {index} is outside 1..={number_of_vertices}")]
    VertexIndexOutOfRange {
        /// Rejected one-based index.
        index: usize,
        /// Number of vertices.
        number_of_vertices: usize,
    },
    /// A named endpoint was not found.
    #[error("cannot find vertex named {name}")]
    UnknownVertexName {
        /// Rejected name.
        name: Box<str>,
    },
    /// An adjacency matrix must be square.
    #[error("adjacency row {row} has length {actual}, expected {expected}")]
    NonSquareAdjacency {
        /// Zero-based row position.
        row: usize,
        /// Actual row length.
        actual: usize,
        /// Expected row length.
        expected: usize,
    },
    /// Arrow multiplicities overflow addressable storage.
    #[error("the requested number of arrows overflows addressable storage")]
    ArrowCountOverflow,
    /// An operation received a vertex from another path family.
    #[error("the objects belong to different quivers")]
    DifferentQuivers,
    /// Concurrent opposite construction ended without a live cache entry.
    #[error("concurrent opposite-quiver cache initialization failed")]
    OppositeCacheRace,
    /// A Dynkin constructor received a rank outside QPA's accepted range.
    #[error("{family}_{rank} is not accepted by QPA's Dynkin constructor")]
    InvalidDynkinRank {
        /// Requested family.
        family: DynkinFamily,
        /// Requested rank.
        rank: usize,
    },
    /// A Dynkin constructor received the wrong number of orientation symbols.
    #[error("{family} orientation has length {actual}, expected {expected}")]
    DynkinOrientationCount {
        /// Requested family.
        family: DynkinFamily,
        /// Actual number of symbols.
        actual: usize,
        /// Required number of symbols.
        expected: usize,
    },
    /// A symbol is not accepted at this Dynkin orientation position.
    #[error("invalid {family} orientation at one-based position {position}")]
    InvalidDynkinOrientation {
        /// Requested family.
        family: DynkinFamily,
        /// One-based orientation position.
        position: usize,
    },
}

/// Failure to construct or combine QPA paths.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PathError {
    /// QPA only installs these operations for an identical path family.
    #[error("the paths belong to different quivers")]
    DifferentQuivers,
    /// A quiver has no multiplicative identity for an empty generic product.
    #[error("an empty generator product is undefined")]
    EmptyProduct,
    /// A generator position lies outside the vertex-plus-arrow list.
    #[error("generator index {index} is outside 1..={number_of_generators}")]
    GeneratorIndexOutOfRange {
        /// Rejected one-based position.
        index: usize,
        /// Number of available generators.
        number_of_generators: usize,
    },
    /// QPA does not bind `WalkOfPath` on the zero path.
    #[error("the zero path has no WalkOfPath value")]
    ZeroHasNoWalk,
    /// QPA's remove-last-three-characters closure received a shorter name.
    #[error("opposite generator name {name} has fewer than three characters")]
    OppositeNameTooShort {
        /// Rejected generator name.
        name: Box<str>,
    },
    /// QPA's component lookup did not find a mapped generator.
    #[error("cannot find opposite generator named {name}")]
    MissingOppositeGenerator {
        /// Missing mapped name.
        name: Box<str>,
    },
    /// An opposite was constructed without QPA's cached name map.
    #[error("the opposite-quiver name map is unavailable")]
    MissingOppositeNameMap,
    /// `NextPath` was explicitly called with a zero path at the queue head.
    #[error("QPA cannot request outgoing arrows of the zero-path endpoint")]
    NextPathZeroEndpoint,
    /// Opposite-quiver construction failed.
    #[error(transparent)]
    QuiverConstruction(#[from] QuiverError),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum VisitColor {
    White,
    Gray,
    Black,
}

fn visit_acyclic(
    quiver: &Quiver,
    vertex: usize,
    colors: &mut [VisitColor],
    sorted: &mut Vec<usize>,
) -> bool {
    colors[vertex] = VisitColor::Gray;
    for neighbor in &quiver.0.neighbors[vertex] {
        let neighbor = neighbor.zero_based();
        if colors[neighbor] == VisitColor::Gray {
            return false;
        }
        if colors[neighbor] == VisitColor::White && !visit_acyclic(quiver, neighbor, colors, sorted)
        {
            return false;
        }
    }
    sorted.push(vertex);
    colors[vertex] = VisitColor::Black;
    true
}

fn visit_connected(
    quiver: &Quiver,
    vertex: usize,
    visited: &mut [bool],
    postorder: &mut Vec<usize>,
) {
    visited[vertex] = true;
    let mut adjacent = quiver.0.neighbors[vertex]
        .iter()
        .map(|index| index.zero_based())
        .collect::<Vec<_>>();
    for arrow in &quiver.0.incoming[vertex] {
        let source = quiver.0.arrows[arrow.zero_based()].source.zero_based();
        if !adjacent.contains(&source) {
            adjacent.push(source);
        }
    }
    for neighbor in adjacent {
        if !visited[neighbor] {
            visit_connected(quiver, neighbor, visited, postorder);
        }
    }
    postorder.push(vertex);
}

fn visit_unoriented_acyclic(
    quiver: &Quiver,
    vertex: usize,
    predecessor: Option<usize>,
    colors: &mut [VisitColor],
) -> bool {
    colors[vertex] = VisitColor::Gray;
    let mut adjacent = Vec::new();
    for arrow in &quiver.0.incoming[vertex] {
        let source = quiver.0.arrows[arrow.zero_based()].source.zero_based();
        if adjacent.contains(&source) {
            return false;
        }
        adjacent.push(source);
    }
    for arrow in &quiver.0.outgoing[vertex] {
        let target = quiver.0.arrows[arrow.zero_based()].target.zero_based();
        if adjacent.contains(&target) {
            return false;
        }
        adjacent.push(target);
    }

    for neighbor in adjacent {
        if colors[neighbor] == VisitColor::Gray && predecessor != Some(neighbor) {
            return false;
        }
        if colors[neighbor] == VisitColor::White
            && !visit_unoriented_acyclic(quiver, neighbor, Some(vertex), colors)
        {
            return false;
        }
    }
    colors[vertex] = VisitColor::Black;
    true
}

fn dynkin_vertex_names(rank: usize) -> Vec<Box<str>> {
    (1..=rank)
        .map(|index| index.to_string().into_boxed_str())
        .collect()
}

fn dynkin_line_arrow(
    names: &[Box<str>],
    position: usize,
    orientation: DynkinOrientation,
    arrow_number: usize,
) -> NamedArrowSpec {
    let (source, target) = if orientation == DynkinOrientation::Left {
        (position + 1, position)
    } else {
        (position, position + 1)
    };
    NamedArrowSpec::named(
        names[source].as_ref(),
        names[target].as_ref(),
        format!("a_{arrow_number}"),
    )
}

fn validate_orientation_count(
    family: DynkinFamily,
    expected: usize,
    actual: usize,
) -> Result<(), QuiverError> {
    if actual == expected {
        Ok(())
    } else {
        Err(QuiverError::DynkinOrientationCount {
            family,
            actual,
            expected,
        })
    }
}

fn validate_horizontal_orientations(
    family: DynkinFamily,
    orientation: &[DynkinOrientation],
) -> Result<(), QuiverError> {
    if let Some(position) = orientation
        .iter()
        .position(|symbol| !matches!(symbol, DynkinOrientation::Left | DynkinOrientation::Right))
    {
        Err(QuiverError::InvalidDynkinOrientation {
            family,
            position: position + 1,
        })
    } else {
        Ok(())
    }
}

fn ensure_uniform_naming(mut names: impl Iterator<Item = bool>) -> Result<(), QuiverError> {
    let Some(first) = names.next() else {
        return Ok(());
    };
    if names.any(|named| named != first) {
        Err(QuiverError::MixedArrowNaming)
    } else {
        Ok(())
    }
}

fn validate_vertex_index(index: VertexIndex, count: usize) -> Result<(), QuiverError> {
    if index.get() <= count {
        Ok(())
    } else {
        Err(QuiverError::VertexIndexOutOfRange {
            index: index.get(),
            number_of_vertices: count,
        })
    }
}

fn position_of_name(names: &[Box<str>], name: &str) -> Option<VertexIndex> {
    names
        .iter()
        .position(|candidate| candidate.as_ref() == name)
        .map(VertexIndex::from_zero_based)
}

#[cfg(test)]
mod tests {
    use std::fmt::Write as _;

    use numerica::domains::integer::Integer;
    use serde_json::{Value, json};

    use super::{
        ArrowIndex, ArrowSpec, DynkinFamily, DynkinOrientation, DynkinType, NamedArrowSpec,
        NextPathPosition, Path, PathEndpoint, Quiver, QuiverCardinality, VertexIndex,
        unique_quiver_name,
    };
    use crate::{Generator, QuiverOrdering};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn arrow(index: usize) -> ArrowIndex {
        ArrowIndex::new(index).unwrap()
    }

    fn line_with_return() -> Quiver {
        Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
                ArrowSpec::named(vertex(3), vertex(2), "c"),
            ],
        )
        .unwrap()
    }

    fn quiver_from_edges(number_of_vertices: usize, edges: &[(usize, usize, &str)]) -> Quiver {
        Quiver::from_indexed(
            number_of_vertices,
            edges
                .iter()
                .map(|(source, target, name)| {
                    ArrowSpec::named(vertex(*source), vertex(*target), *name)
                })
                .collect(),
        )
        .unwrap()
    }

    fn quiver_record(quiver: &Quiver) -> Value {
        json!({
            "vertices": quiver
                .vertices()
                .map(|vertex| vertex.name().to_owned())
                .collect::<Vec<_>>(),
            "arrows": quiver
                .arrows()
                .map(|arrow| json!({
                    "name": arrow.name(),
                    "source": arrow.source().name(),
                    "target": arrow.target().name()
                }))
                .collect::<Vec<_>>()
        })
    }

    fn cardinality_json(cardinality: QuiverCardinality) -> Value {
        match cardinality {
            QuiverCardinality::Finite(size) => json!(size.to_string().parse::<u64>().unwrap()),
            QuiverCardinality::Infinite => json!("infinity"),
        }
    }

    #[test]
    fn indexed_constructor_preserves_qpa_generator_numbering() {
        let quiver = line_with_return();
        assert_eq!(quiver.number_of_vertices(), 3);
        assert_eq!(quiver.number_of_arrows(), 3);
        assert!(quiver.has_special_single_letter_arrow_names());
        assert_eq!(quiver.vertex(vertex(2)).unwrap().name(), "v2");
        assert_eq!(quiver.arrow(arrow(2)).unwrap().generator_index(), 5);
        assert_eq!(
            quiver.adjacency_matrix(),
            &[
                Box::from([0, 1, 0]),
                Box::from([0, 0, 1]),
                Box::from([0, 1, 0])
            ]
        );
    }

    #[test]
    fn path_multiplication_is_left_to_right_and_noncomposable_is_zero() {
        let quiver = line_with_return();
        let a = quiver.arrow(arrow(1)).unwrap().path();
        let b = quiver.arrow(arrow(2)).unwrap().path();
        let c = quiver.arrow(arrow(3)).unwrap().path();

        let ab = a.multiply(&b).unwrap();
        assert_eq!(ab.to_string(), "a*b");
        assert_eq!(ab.length(), Some(2));
        assert_eq!(ab.external_representation(), vec![4, 5]);
        assert_eq!(b.multiply(&a).unwrap(), quiver.zero());
        assert_eq!(ab.multiply(&c).unwrap().to_string(), "a*b*c");
    }

    #[test]
    fn vertices_are_local_units_and_zero_is_its_own_endpoint() {
        let quiver = line_with_return();
        let v1 = quiver.vertex(vertex(1)).unwrap().path();
        let v2 = quiver.vertex(vertex(2)).unwrap().path();
        let a = quiver.arrow(arrow(1)).unwrap().path();

        assert_eq!(v1.multiply(&a).unwrap(), a);
        assert_eq!(v2.multiply(&a).unwrap(), quiver.zero());
        assert_eq!(v1.multiply(&v1).unwrap(), v1);
        assert!(matches!(quiver.zero().source(), PathEndpoint::Zero(_)));
        assert_eq!(quiver.zero().length(), None);
        assert_eq!(quiver.zero().walk(), None);
    }

    #[test]
    fn external_representations_round_trip_qpa_generators_and_walks() {
        let quiver = line_with_return();
        for representation in [vec![0], vec![2], vec![4], vec![4, 5, 6]] {
            let path = Path::from_external_representation(&quiver, &representation)
                .unwrap()
                .unwrap();
            assert_eq!(path.external_representation(), representation);
        }
        assert_eq!(
            Path::from_external_representation(&quiver, &[]).unwrap(),
            None
        );
    }

    #[test]
    fn default_ordering_is_length_then_left_lexicographic() {
        let quiver = line_with_return();
        let first_ordering = quiver.ordering().unwrap();
        let second_ordering = quiver.ordering().unwrap();
        assert!(first_ordering.is_identical(&second_ordering));
        let zero = quiver.zero();
        let v1 = quiver.vertex(vertex(1)).unwrap().path();
        let v2 = quiver.vertex(vertex(2)).unwrap().path();
        let a = quiver.arrow(arrow(1)).unwrap().path();
        let ab = a.multiply(&quiver.arrow(arrow(2)).unwrap().path()).unwrap();

        assert!(quiver.compare_paths(&zero, &v1).unwrap().is_lt());
        assert!(quiver.compare_paths(&v1, &v2).unwrap().is_lt());
        assert!(quiver.compare_paths(&v2, &a).unwrap().is_lt());
        assert!(quiver.compare_paths(&a, &ab).unwrap().is_lt());
    }

    #[test]
    fn next_path_order_is_qpa_breadth_first_order() {
        let quiver = Quiver::from_named(
            vec!["x".into(), "y".into()],
            vec![
                NamedArrowSpec::named("x", "y", "a"),
                NamedArrowSpec::named("x", "x", "b"),
            ],
        )
        .unwrap();
        let first_six = quiver
            .paths()
            .take(6)
            .map(|path| path.to_string())
            .collect::<Vec<_>>();
        assert_eq!(first_six, ["x", "y", "a", "b", "b*a", "b^2"]);
    }

    #[test]
    fn opposite_is_cached_and_reverses_walks() {
        let quiver = line_with_return();
        let opposite = quiver.opposite().unwrap();
        assert!(opposite.opposite().unwrap().is_identical(&quiver));
        let ab = quiver
            .arrow(arrow(1))
            .unwrap()
            .path()
            .multiply(&quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let opposite_path = quiver.opposite_path(&ab).unwrap();
        assert_eq!(opposite_path.to_string(), "b_op*a_op");
        assert_eq!(opposite_path.external_representation(), vec![5, 4]);
    }

    #[test]
    fn acyclicity_uses_qpa_outgoing_neighbor_graph() {
        let cyclic = line_with_return();
        assert!(!cyclic.is_acyclic());

        let acyclic = Quiver::from_adjacency(vec![vec![0, 1], vec![0, 0]]).unwrap();
        assert!(acyclic.is_acyclic());
        assert_eq!(
            acyclic
                .topological_order()
                .unwrap()
                .iter()
                .map(super::Vertex::name)
                .collect::<Vec<_>>(),
            ["v1", "v2"]
        );
    }

    #[test]
    fn graph_properties_and_cardinality_follow_qpa_examples() {
        let line =
            Quiver::from_adjacency(vec![vec![0, 1, 0], vec![0, 0, 1], vec![0, 0, 0]]).unwrap();
        assert_eq!(
            line.cardinality(),
            QuiverCardinality::Finite(Integer::from(7))
        );
        assert!(line.is_connected());
        assert!(line.is_unoriented_acyclic());
        assert!(line.is_tree());
        assert_eq!(
            line.dynkin_type(),
            Some(DynkinType {
                family: DynkinFamily::A,
                rank: 3
            })
        );

        let parallel =
            Quiver::from_adjacency(vec![vec![0, 2, 0], vec![0, 0, 1], vec![0, 0, 0]]).unwrap();
        assert_eq!(
            parallel.cardinality(),
            QuiverCardinality::Finite(Integer::from(9))
        );
        assert!(parallel.is_connected());
        assert!(!parallel.is_unoriented_acyclic());
        assert!(!parallel.is_tree());
        assert_eq!(
            line_with_return().cardinality(),
            QuiverCardinality::Infinite
        );
    }

    #[test]
    fn full_subquiver_and_components_preserve_qpa_dfs_order() {
        let quiver = line_with_return();
        let subquiver = quiver
            .full_subquiver(&[
                quiver.vertex(vertex(3)).unwrap(),
                quiver.vertex(vertex(2)).unwrap(),
            ])
            .unwrap();
        assert_eq!(
            subquiver
                .vertices()
                .map(|vertex| vertex.name().to_owned())
                .collect::<Vec<_>>(),
            ["v3", "v2"]
        );
        assert_eq!(
            subquiver
                .arrows()
                .map(|arrow| arrow.name().to_owned())
                .collect::<Vec<_>>(),
            ["b", "c"]
        );

        let disconnected = Quiver::from_named(
            vec!["p".into(), "q".into(), "r".into()],
            vec![NamedArrowSpec::named("p", "q", "x")],
        )
        .unwrap();
        let components = disconnected.connected_components().unwrap();
        assert_eq!(components.len(), 2);
        assert_eq!(
            components[0]
                .vertices()
                .map(|vertex| vertex.name().to_owned())
                .collect::<Vec<_>>(),
            ["q", "p"]
        );
        assert_eq!(
            components[1]
                .vertices()
                .map(|vertex| vertex.name().to_owned())
                .collect::<Vec<_>>(),
            ["r"]
        );

        let foreign = Quiver::from_indexed(1, Vec::new()).unwrap();
        assert_eq!(
            disconnected.full_subquiver(&[foreign.vertex(vertex(1)).unwrap()]),
            Err(super::QuiverError::DifferentQuivers)
        );
    }

    #[test]
    fn dynkin_constructors_translate_qpa_orientations() {
        use DynkinOrientation::{Down, Left, Right};

        let a = Quiver::dynkin_a(4, &[Right, Left, Right]).unwrap();
        let d = Quiver::dynkin_d(5, &[Right, Left, Right, Left]).unwrap();
        let e = Quiver::dynkin_e(6, &[Right, Left, Right, Left, Down]).unwrap();

        let arrow_specs = |quiver: &Quiver| {
            quiver
                .arrows()
                .map(|arrow| {
                    (
                        arrow.source().name().to_owned(),
                        arrow.target().name().to_owned(),
                        arrow.name().to_owned(),
                    )
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(
            arrow_specs(&a),
            [
                ("1".into(), "2".into(), "a_1".into()),
                ("3".into(), "2".into(), "a_2".into()),
                ("3".into(), "4".into(), "a_3".into())
            ]
        );
        assert_eq!(
            arrow_specs(&d),
            [
                ("1".into(), "3".into(), "a_1".into()),
                ("3".into(), "2".into(), "a_2".into()),
                ("3".into(), "4".into(), "a_3".into()),
                ("5".into(), "4".into(), "a_4".into())
            ]
        );
        assert_eq!(
            arrow_specs(&e),
            [
                ("1".into(), "2".into(), "a_1".into()),
                ("3".into(), "2".into(), "a_2".into()),
                ("3".into(), "4".into(), "a_3".into()),
                ("5".into(), "4".into(), "a_4".into()),
                ("6".into(), "3".into(), "a_5".into())
            ]
        );
        assert_eq!(
            d.dynkin_type(),
            Some(DynkinType {
                family: DynkinFamily::D,
                rank: 5
            })
        );
        assert_eq!(
            e.dynkin_type(),
            Some(DynkinType {
                family: DynkinFamily::E,
                rank: 6
            })
        );
    }

    #[test]
    fn dynkin_decision_matches_pinned_qpa_output() {
        let cases = [
            (
                "A4",
                quiver_from_edges(4, &[(1, 2, "a"), (2, 3, "b"), (3, 4, "c")]),
            ),
            (
                "D5",
                quiver_from_edges(5, &[(1, 3, "a"), (2, 3, "b"), (3, 4, "c"), (4, 5, "d")]),
            ),
            (
                "E8",
                quiver_from_edges(
                    8,
                    &[
                        (1, 2, "a1"),
                        (1, 3, "b1"),
                        (3, 4, "b2"),
                        (1, 5, "c1"),
                        (5, 6, "c2"),
                        (6, 7, "c3"),
                        (7, 8, "c4"),
                    ],
                ),
            ),
            ("DISCONNECTED", quiver_from_edges(2, &[])),
            ("CYCLE", quiver_from_edges(2, &[(1, 2, "a"), (2, 1, "b")])),
            (
                "FOUR_ARMS",
                quiver_from_edges(5, &[(1, 2, "a"), (1, 3, "b"), (1, 4, "c"), (1, 5, "d")]),
            ),
            (
                "LONG_STAR",
                quiver_from_edges(
                    9,
                    &[
                        (1, 2, "a1"),
                        (1, 3, "b1"),
                        (3, 4, "b2"),
                        (1, 5, "c1"),
                        (5, 6, "c2"),
                        (6, 7, "c3"),
                        (7, 8, "c4"),
                        (8, 9, "c5"),
                    ],
                ),
            ),
        ];

        let mut actual = String::new();
        for (label, quiver) in cases {
            let decision = quiver.dynkin_decision();
            writeln!(actual, "BEGIN_{label}").unwrap();
            for message in decision.messages() {
                actual.push_str(message);
            }
            writeln!(actual, "RESULT_{label}={}", decision.is_dynkin()).unwrap();
        }
        assert_eq!(
            actual,
            include_str!("../../../tests/fixtures/quiver_dynkin_output.txt")
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_quiver_graph_fixture_matches() {
        use DynkinOrientation::{Down, Left, Right};

        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/quiver_graph.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1"
            })
        );

        let line =
            Quiver::from_adjacency(vec![vec![0, 1, 0], vec![0, 0, 1], vec![0, 0, 0]]).unwrap();
        let parallel =
            Quiver::from_adjacency(vec![vec![0, 2, 0], vec![0, 0, 1], vec![0, 0, 0]]).unwrap();
        let cyclic = line_with_return();
        let properties = |quiver: &Quiver| {
            json!({
                "size": cardinality_json(quiver.cardinality()),
                "finite": quiver.is_finite(),
                "connected": quiver.is_connected(),
                "unoriented_acyclic": quiver.is_unoriented_acyclic(),
                "tree": quiver.is_tree()
            })
        };
        assert_eq!(
            json!({
                "line": properties(&line),
                "parallel": properties(&parallel),
                "cyclic": properties(&cyclic)
            }),
            expected["properties"]
        );

        let disconnected = Quiver::from_named(
            vec!["p".into(), "q".into(), "r".into()],
            vec![NamedArrowSpec::named("p", "q", "x")],
        )
        .unwrap();
        assert_eq!(
            Value::Array(
                disconnected
                    .connected_components()
                    .unwrap()
                    .iter()
                    .map(quiver_record)
                    .collect()
            ),
            expected["components"]
        );
        assert_eq!(
            quiver_record(
                &cyclic
                    .full_subquiver(&[
                        cyclic.vertex(vertex(3)).unwrap(),
                        cyclic.vertex(vertex(2)).unwrap(),
                    ])
                    .unwrap()
            ),
            expected["full_subquiver"]
        );

        let a4 = Quiver::dynkin_a(4, &[Right, Left, Right]).unwrap();
        let d5 = Quiver::dynkin_d(5, &[Right, Left, Right, Left]).unwrap();
        let e6 = Quiver::dynkin_e(6, &[Right, Left, Right, Left, Down]).unwrap();
        assert_eq!(
            json!({
                "a4": quiver_record(&a4),
                "d5": quiver_record(&d5),
                "e6": quiver_record(&e6)
            }),
            expected["dynkin"]
        );
        assert_eq!(
            quiver_record(&cyclic.separated().unwrap()),
            expected["separated"]
        );
        assert_eq!(quiver_record(&cyclic.double().unwrap()), expected["double"]);

        let universe = [
            ArrowSpec::named(vertex(1), vertex(1), "a"),
            ArrowSpec::named(vertex(1), vertex(2), "b"),
            ArrowSpec::named(vertex(2), vertex(1), "c"),
            ArrowSpec::named(vertex(2), vertex(2), "d"),
        ];
        let mut exhaustive = Vec::new();
        for quiver_mask in 0_u8..16 {
            let specifications = universe
                .iter()
                .enumerate()
                .filter(|(index, _)| quiver_mask & (1 << index) != 0)
                .map(|(_, specification)| specification.clone())
                .collect();
            let quiver = Quiver::from_indexed(2, specifications).unwrap();
            let topological = quiver
                .topological_order()
                .map(|vertices| {
                    vertices
                        .into_iter()
                        .map(|vertex| vertex.index().get())
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();
            exhaustive.push(json!([
                quiver_mask,
                quiver.is_acyclic(),
                quiver.is_finite(),
                quiver.is_connected(),
                quiver.is_unoriented_acyclic(),
                quiver.is_tree(),
                topological,
            ]));
        }
        assert_eq!(json!(exhaustive), expected["exhaustive"]);
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_quiver_path_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/quiver_paths.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            })
        );

        let quiver = line_with_return();
        let a = quiver.arrow(arrow(1)).unwrap().path();
        let b = quiver.arrow(arrow(2)).unwrap().path();
        let c = quiver.arrow(arrow(3)).unwrap().path();
        let ab = a.multiply(&b).unwrap();
        let abc = ab.multiply(&c).unwrap();

        let endpoint = |endpoint: PathEndpoint| match endpoint {
            PathEndpoint::Vertex(vertex) => vertex.to_string(),
            PathEndpoint::Zero(path) => path.to_string(),
        };
        let path_record = |label: &str, path: &Path| {
            json!({
                "label": label,
                "display": path.to_string(),
                "external_representation": path.external_representation(),
                "length": path.length(),
                "walk": path.walk().map(|walk| {
                    walk.iter().map(|index| index.get()).collect::<Vec<_>>()
                }),
                "source": endpoint(path.source()),
                "target": endpoint(path.target())
            })
        };
        let actual_quiver = json!({
            "vertices": quiver.vertices().map(|vertex| vertex.name().to_owned()).collect::<Vec<_>>(),
            "arrows": quiver.arrows().map(|arrow| json!({
                "name": arrow.name(),
                "index": arrow.index().get(),
                "generator_index": arrow.generator_index(),
                "source": arrow.source().index().get(),
                "target": arrow.target().index().get()
            })).collect::<Vec<_>>(),
            "adjacency": quiver.adjacency_matrix().iter().map(|row| row.to_vec()).collect::<Vec<_>>(),
            "special_arrow_names": quiver.has_special_single_letter_arrow_names(),
            "acyclic": quiver.is_acyclic()
        });
        assert_eq!(actual_quiver, expected["quiver"]);

        assert_eq!(
            json!([
                path_record("zero", &quiver.zero()),
                path_record("v2", &quiver.vertex(vertex(2)).unwrap().path()),
                path_record("a", &a),
                path_record("ab", &ab),
                path_record("abc", &abc)
            ]),
            expected["paths"]
        );
        assert_eq!(
            json!({
                "a_b": a.multiply(&b).unwrap().external_representation(),
                "b_a": b.multiply(&a).unwrap().external_representation(),
                "v1_a": quiver.vertex(vertex(1)).unwrap().path().multiply(&a).unwrap().external_representation(),
                "v2_a": quiver.vertex(vertex(2)).unwrap().path().multiply(&a).unwrap().external_representation(),
                "a_v2": a.multiply(&quiver.vertex(vertex(2)).unwrap().path()).unwrap().external_representation(),
                "a_v1": a.multiply(&quiver.vertex(vertex(1)).unwrap().path()).unwrap().external_representation()
            }),
            expected["products"]
        );
        assert_eq!(
            json!({
                "zero_lt_v1": quiver.compare_paths(&quiver.zero(), &quiver.vertex(vertex(1)).unwrap().path()).unwrap().is_lt(),
                "v1_lt_v2": quiver.compare_paths(&quiver.vertex(vertex(1)).unwrap().path(), &quiver.vertex(vertex(2)).unwrap().path()).unwrap().is_lt(),
                "v2_lt_a": quiver.compare_paths(&quiver.vertex(vertex(2)).unwrap().path(), &a).unwrap().is_lt(),
                "a_lt_ab": quiver.compare_paths(&a, &ab).unwrap().is_lt()
            }),
            expected["ordering"]
        );

        let iterator_quiver = Quiver::from_named(
            vec!["x".into(), "y".into()],
            vec![
                NamedArrowSpec::named("x", "y", "a"),
                NamedArrowSpec::named("x", "x", "b"),
            ],
        )
        .unwrap();
        assert_eq!(
            json!(
                iterator_quiver
                    .paths()
                    .take(8)
                    .map(|path| path.to_string())
                    .collect::<Vec<_>>()
            ),
            expected["iterator_first"]
        );

        let opposite = quiver.opposite().unwrap();
        let opposite_ab = quiver.opposite_path(&ab).unwrap();
        assert_eq!(
            json!({
                "vertices": opposite.vertices().map(|vertex| vertex.name().to_owned()).collect::<Vec<_>>(),
                "arrows": opposite.arrows().map(|arrow| arrow.name().to_owned()).collect::<Vec<_>>(),
                "ab_display": opposite_ab.to_string(),
                "ab_external_representation": opposite_ab.external_representation(),
                "double_opposite_is_identical": opposite.opposite().unwrap().is_identical(&quiver)
            }),
            expected["opposite"]
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_quiver_surface_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/quiver_surface.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1"
            })
        );

        let quiver = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(1), vertex(1), "b"),
            ],
        )
        .unwrap();
        let indexed = Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::unnamed(vertex(1), vertex(2)),
                ArrowSpec::unnamed(vertex(2), vertex(1)),
            ],
        )
        .unwrap();
        let named_quiver = Quiver::from_named(
            vec!["p".into(), "q".into()],
            vec![
                NamedArrowSpec::unnamed("p", "q"),
                NamedArrowSpec::unnamed("q", "p"),
            ],
        )
        .unwrap();
        let adjacency = Quiver::from_adjacency(vec![vec![0, 2], vec![0, 0]]).unwrap();
        let first_default_ordering = quiver.ordering().unwrap();
        let second_default_ordering = quiver.ordering().unwrap();
        let names = |quiver: &Quiver| {
            quiver
                .vertices()
                .map(|value| value.name().to_owned())
                .collect::<Vec<_>>()
        };
        let arrow_names = |quiver: &Quiver| {
            quiver
                .arrows()
                .map(|value| value.name().to_owned())
                .collect::<Vec<_>>()
        };
        assert_eq!(
            json!({
                "indexed_vertices": names(&indexed),
                "indexed_arrows": arrow_names(&indexed),
                "named_vertices": names(&named_quiver),
                "named_arrows": arrow_names(&named_quiver),
                "adjacency_vertices": names(&adjacency),
                "adjacency_arrows": arrow_names(&adjacency),
                "adjacency_matrix": adjacency.adjacency_matrix().iter().map(|row| row.to_vec()).collect::<Vec<_>>(),
                "special_named_arrows": quiver.has_special_single_letter_arrow_names(),
                "default_ordering_cached": first_default_ordering.is_identical(&second_default_ordering)
            }),
            expected["constructors"]
        );

        let v1 = quiver.vertex(vertex(1)).unwrap();
        let v2 = quiver.vertex(vertex(2)).unwrap();
        let a = quiver.arrow(arrow(1)).unwrap();
        let b = quiver.arrow(arrow(2)).unwrap();
        assert_eq!(
            json!({
                "incoming_y": quiver.incoming_arrows(&v2).unwrap().iter().map(ToString::to_string).collect::<Vec<_>>(),
                "outgoing_x": quiver.outgoing_arrows(&v1).unwrap().iter().map(ToString::to_string).collect::<Vec<_>>(),
                "in_degree_y": v2.in_degree(),
                "out_degree_x": v1.out_degree(),
                "neighbors_x": quiver.neighbors(&v1).unwrap().iter().map(ToString::to_string).collect::<Vec<_>>(),
                "piece_lookup_is_identical": matches!(quiver.generator_by_name("a"), Some(Generator::Arrow(found)) if found == a)
            }),
            expected["attributes"]
        );
        assert_eq!(
            json!({
                "quiver": quiver.to_string(),
                "view": quiver.view_string(),
                "power_path": b.path().multiply(&b.path()).unwrap().multiply(&a.path()).unwrap().to_string()
            }),
            expected["printing"]
        );

        let mut preferred = quiver.generators().collect::<Vec<_>>();
        preferred.reverse();
        let ordering = QuiverOrdering::right_lexicographic(&quiver, &preferred).unwrap();
        let ordered = quiver.ordered_by(ordering.clone()).unwrap();
        let ordered_v1 = ordered.vertex(vertex(1)).unwrap().path();
        let ordered_v2 = ordered.vertex(vertex(2)).unwrap().path();
        let ordered_a = ordered.arrow(arrow(1)).unwrap().path();
        let ordered_b = ordered.arrow(arrow(2)).unwrap().path();
        assert_eq!(
            json!({
                "ordering_is_identical": ordered.ordering().is_identical(&ordering),
                "new_family": !ordered_v1.quiver().is_identical(&quiver),
                "vertices": names(&ordered),
                "arrows": arrow_names(&ordered),
                "v1_lt_v2": ordered.compare_paths(&ordered_v1, &ordered_v2).unwrap().is_lt(),
                "a_lt_b": ordered.compare_paths(&ordered_a, &ordered_b).unwrap().is_lt(),
                "ba": ordered_b.multiply(&ordered_a).unwrap().to_string()
            }),
            expected["ordered"]
        );
        let opposite = quiver.opposite().unwrap();
        assert_eq!(
            json!({
                "append": quiver.opposite_quiver_name("plain").unwrap(),
                "remove_without_suffix_check": opposite.opposite_quiver_name("abcd").unwrap()
            }),
            expected["opposite_name_map"]
        );

        let first = quiver.next_path(&NextPathPosition::Start).unwrap().unwrap();
        let second = quiver.next_path(&first.1).unwrap().unwrap();
        let empty = Quiver::from_named(Vec::new(), Vec::new()).unwrap();
        let mut empty_iterator = empty.paths();
        let empty_done_before = empty_iterator.is_done();
        let empty_next_is_fail = empty_iterator.next().is_none();
        let empty_done_after = empty_iterator.is_done();
        let queue_strings = |position: &NextPathPosition| match position {
            NextPathPosition::Start => Vec::new(),
            NextPathPosition::Queue(queue) => {
                queue.iter().map(ToString::to_string).collect::<Vec<_>>()
            }
        };
        assert_eq!(
            json!({
                "first": first.0.to_string(),
                "first_queue": queue_strings(&first.1),
                "second": second.0.to_string(),
                "second_queue": queue_strings(&second.1),
                "enumerator_5": quiver.path_at(5).unwrap().to_string(),
                "empty_next_is_fail": empty_next_is_fail,
                "empty_done_before": empty_done_before,
                "empty_done_after": empty_done_after
            }),
            expected["next_path"]
        );

        assert_eq!(
            json!({
                "vertex": v2.index().get(),
                "arrow": b.index().get(),
                "generator_vertex": v2.index().get(),
                "generator_arrow": b.generator_index(),
                "quiver_containing_path": b.path().multiply(&a.path()).unwrap().quiver().is_identical(&quiver)
            }),
            expected["indices"]
        );
        let unique1 = unique_quiver_name();
        let unique2 = unique_quiver_name();
        assert_eq!(
            json!({
                "different": unique1 != unique2,
                "prefixes": unique1.starts_with("Quiver_") && unique2.starts_with("Quiver_")
            }),
            expected["unique_names"]
        );
    }

    #[test]
    fn pinned_qpa_quiver_dispatch_fixture_matches() {
        use DynkinOrientation::{Right, Up};

        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/quiver_dispatch.schema-v1.json"
        ))
        .unwrap();
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let foreign = Quiver::from_indexed(1, Vec::new()).unwrap();

        assert_eq!(
            json!([
                Quiver::from_indexed(0, Vec::new()).is_ok(),
                Quiver::from_indexed(
                    2,
                    vec![
                        ArrowSpec::unnamed(vertex(1), vertex(2)),
                        ArrowSpec::named(vertex(2), vertex(1), "b")
                    ]
                )
                .is_ok(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(3), "a")]).is_ok(),
                Quiver::from_named(
                    vec!["x".into(), "y".into()],
                    vec![NamedArrowSpec::named("x", "z", "a")]
                )
                .is_ok(),
                Quiver::from_adjacency(vec![vec![0, 1], vec![0]]).is_ok()
            ]),
            expected["constructor_success"]
        );

        let empty_product = Path::from_arrow_product(&quiver, &[]);
        let foreign_product = quiver
            .arrow(arrow(1))
            .unwrap()
            .path()
            .multiply(&foreign.vertex(vertex(1)).unwrap().path());
        let empty_extrep = Path::from_external_representation(&quiver, &[]);
        let next_zero = quiver.next_path(&NextPathPosition::Queue(vec![quiver.zero()]));
        let opposite_zero = quiver.opposite_path(&quiver.zero());
        assert_eq!(
            json!([
                empty_product.is_ok(),
                foreign_product.is_ok(),
                empty_extrep.is_ok(),
                next_zero.is_ok(),
                opposite_zero.is_ok()
            ]),
            expected["path_success"]
        );
        assert_eq!(
            empty_product.unwrap().is_some(),
            expected["path_empty_result_is_one"]
        );
        let extrep_mixed_zero = Path::from_external_representation(&quiver, &[0, 1, 3]);
        assert_eq!(
            json!({
                "success": extrep_mixed_zero.is_ok(),
                "is_zero": extrep_mixed_zero.unwrap().unwrap().is_zero()
            }),
            expected["extrep_mixed_zero"]
        );

        let foreign_ordering = QuiverOrdering::right_lexicographic(
            &foreign,
            &foreign.generators().collect::<Vec<_>>(),
        )
        .unwrap();
        let foreign_ordered = quiver.ordered_by(foreign_ordering).unwrap();
        let foreign_order_compare = foreign_ordered.compare_paths(
            &foreign_ordered.vertex(vertex(1)).unwrap().path(),
            &foreign_ordered.vertex(vertex(2)).unwrap().path(),
        );
        assert_eq!(
            json!(true),
            expected["foreign_ordered_construction_success"]
        );
        assert_eq!(
            json!(foreign_order_compare.is_ok()),
            expected["foreign_order_compare_success"]
        );

        assert_eq!(
            json!([
                Quiver::dynkin_a(0, &[]).is_ok(),
                Quiver::dynkin_a(3, &[Right]).is_ok(),
                Quiver::dynkin_a(2, &[Right, DynkinOrientation::Left]).is_ok(),
                Quiver::dynkin_a(2, &[Up]).is_ok(),
                Quiver::dynkin_d(3, &[Right, Right]).is_ok(),
                Quiver::dynkin_e(5, &[]).is_ok(),
                Quiver::dynkin_e(6, &[Right, Right, Right, Right, Right]).is_ok(),
                Quiver::dynkin_from_name("X", 1, &[]).is_ok()
            ]),
            expected["dynkin_success"]
        );
        assert_eq!(Quiver::dynkin_from_name("X", 1, &[]).unwrap(), None);
    }
}
