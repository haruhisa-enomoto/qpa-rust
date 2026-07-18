//! QPA's dynamic suffix-tree dictionary.
//!
//! This module translates `lib/groebnerdyndict.gi`: its online suffix-tree
//! phases, suffix links, terminal leaves, marked nodes, dynamic-forest updates,
//! deletion/compression rules, and the two pattern-occurrence queries. The
//! data lives in typed Rust arenas, but labels retain QPA's one-based inclusive
//! descriptors so each mutation remains directly traceable to the source.

use std::{
    fmt,
    sync::atomic::{AtomicU64, Ordering as AtomicOrdering},
};

use thiserror::Error;

use crate::{DynamicForest, DynamicTree, DynamicTreeError};

/// Rust marker for QPA's `InfoDynamicDictionary` information class.
pub const INFO_DYNAMIC_DICTIONARY: &str = "InfoDynamicDictionary";

static NEXT_SUFFIX_TREE_IDENTITY: AtomicU64 = AtomicU64::new(1);

/// Stable identity of a suffix-tree node.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct SuffixTreeNodeId {
    owner: u64,
    index: usize,
}

/// Stable identity of a suffix-tree edge.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct SuffixTreeEdgeId {
    owner: u64,
    index: usize,
}

/// A label token, including QPA's reserved empty-list termination symbol.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SuffixTreeToken<Symbol> {
    /// One symbol from an inserted pattern.
    Symbol(Symbol),
    /// QPA's appended `[]` token.
    Termination,
}

/// Whether a suffix-tree node has QPA's internal or external representation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SuffixTreeNodeKind {
    /// A branching node with edges, a suffix link, and a marked flag.
    Internal,
    /// A leaf containing one or more pattern/start descriptors.
    External,
}

/// One QPA one-based `[pattern_number, position]` occurrence.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SuffixTreeOccurrence {
    /// One-based insertion number of the stored pattern.
    pub pattern_number: usize,
    /// One-based start position in a pattern or queried sequence.
    pub position: usize,
}

#[derive(Clone, Copy, Debug)]
struct LabelSpan {
    pattern_number: usize,
    start: usize,
    end: Option<usize>,
}

#[derive(Clone, Debug)]
struct StoredPattern<Symbol> {
    symbols: Vec<Symbol>,
    tokens: Vec<SuffixTreeToken<Symbol>>,
    active: bool,
}

#[derive(Debug)]
struct SuffixTreeNode {
    parent_edge: Option<SuffixTreeEdgeId>,
    dynamic_tree: Option<DynamicTree<SuffixTreeNodeId>>,
    representation: SuffixTreeNodeRepresentation,
}

#[derive(Debug)]
enum SuffixTreeNodeRepresentation {
    Internal {
        label: Option<LabelSpan>,
        edges: Vec<SuffixTreeEdgeId>,
        suffix_link: Option<SuffixTreeNodeId>,
        marked: bool,
    },
    External {
        occurrences: Vec<SuffixTreeOccurrence>,
    },
}

#[derive(Clone, Copy, Debug)]
struct SuffixTreeEdge {
    label: LabelSpan,
    parent: Option<SuffixTreeNodeId>,
    child: Option<SuffixTreeNodeId>,
}

#[derive(Clone, Copy, Debug)]
enum ExtensionLocus {
    Node(SuffixTreeNodeId),
    Edge(SuffixTreeEdgeId, usize),
}

#[derive(Clone, Copy, Debug)]
struct ExtensionResult {
    internal: bool,
    current_edge: Option<SuffixTreeEdgeId>,
    location: usize,
    next_node: Option<SuffixTreeNodeId>,
    current_node: SuffixTreeNodeId,
    new_node: Option<SuffixTreeNodeId>,
}

/// QPA's mutable generalized suffix tree used by `TipReduceGroebnerBasis`.
pub struct SuffixTree<Symbol> {
    identity: u64,
    patterns: Vec<StoredPattern<Symbol>>,
    root: SuffixTreeNodeId,
    dead_length: usize,
    nodes: Vec<SuffixTreeNode>,
    edges: Vec<SuffixTreeEdge>,
    dynamic_forest: DynamicForest<SuffixTreeNodeId>,
}

impl<Symbol> SuffixTree<Symbol>
where
    Symbol: Clone + Eq,
{
    /// Translate QPA's `CreateSuffixTree`.
    ///
    /// # Panics
    ///
    /// Panics only if one process exhausts all nonzero 64-bit suffix-tree
    /// identities.
    #[must_use]
    pub fn new() -> Self {
        let identity = NEXT_SUFFIX_TREE_IDENTITY.fetch_add(1, AtomicOrdering::Relaxed);
        assert_ne!(identity, 0, "exhausted suffix-tree identities");
        let dynamic_forest = DynamicForest::new();
        let root = SuffixTreeNodeId {
            owner: identity,
            index: 0,
        };
        let root_dynamic_tree = dynamic_forest.create(Some(root));
        Self {
            identity,
            patterns: Vec::new(),
            root,
            dead_length: 0,
            nodes: vec![SuffixTreeNode {
                parent_edge: None,
                dynamic_tree: Some(root_dynamic_tree),
                representation: SuffixTreeNodeRepresentation::Internal {
                    label: None,
                    edges: Vec::new(),
                    suffix_link: None,
                    marked: false,
                },
            }],
            edges: Vec::new(),
            dynamic_forest,
        }
    }

    /// Identity of QPA's root node.
    #[must_use]
    pub const fn root_node(&self) -> SuffixTreeNodeId {
        self.root
    }

    /// Number of patterns ever inserted; QPA retains slots after deletion.
    #[must_use]
    pub fn number_of_patterns(&self) -> usize {
        self.patterns.len()
    }

    /// QPA's currently unused `deadLength` component.
    #[must_use]
    pub const fn dead_length(&self) -> usize {
        self.dead_length
    }

    /// Translate suffix-tree indexing `tree[pattern_number]`.
    ///
    /// QPA removes the appended termination token and retains deleted pattern
    /// slots, so this accessor does the same.
    ///
    /// # Errors
    ///
    /// Returns an error when the one-based pattern number does not exist.
    pub fn pattern(&self, pattern_number: usize) -> Result<&[Symbol], SuffixTreeError> {
        Ok(&self.stored_pattern(pattern_number)?.symbols)
    }

    /// Whether a retained QPA pattern slot has not been deleted.
    ///
    /// # Errors
    ///
    /// Returns an error when the one-based pattern number does not exist.
    pub fn is_pattern_active(&self, pattern_number: usize) -> Result<bool, SuffixTreeError> {
        Ok(self.stored_pattern(pattern_number)?.active)
    }

    /// Return the representation category of a node.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown node identity.
    pub fn node_kind(&self, node: SuffixTreeNodeId) -> Result<SuffixTreeNodeKind, SuffixTreeError> {
        Ok(match &self.node(node)?.representation {
            SuffixTreeNodeRepresentation::Internal { .. } => SuffixTreeNodeKind::Internal,
            SuffixTreeNodeRepresentation::External { .. } => SuffixTreeNodeKind::External,
        })
    }

    /// Translate QPA's `ParentNode` (`None` represents QPA's `0`).
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown node or a detached corrupt edge.
    pub fn parent_node(
        &self,
        node: SuffixTreeNodeId,
    ) -> Result<Option<SuffixTreeNodeId>, SuffixTreeError> {
        let Some(parent_edge) = self.node(node)?.parent_edge else {
            return Ok(None);
        };
        Ok(self.edge(parent_edge)?.parent)
    }

    /// Translate both methods of `LabelOfSuffixTreeNode`.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid identity, label, or empty external leaf.
    pub fn label_of_node(
        &self,
        node: SuffixTreeNodeId,
    ) -> Result<Vec<SuffixTreeToken<Symbol>>, SuffixTreeError> {
        match &self.node(node)?.representation {
            SuffixTreeNodeRepresentation::Internal { label: None, .. } => Ok(Vec::new()),
            SuffixTreeNodeRepresentation::Internal {
                label: Some(label), ..
            } => self.internal_label(*label),
            SuffixTreeNodeRepresentation::External { occurrences } => {
                let occurrence = occurrences
                    .first()
                    .ok_or(SuffixTreeError::EmptyExternalNode)?;
                let pattern = self.stored_pattern(occurrence.pattern_number)?;
                if occurrence.position == 0 || occurrence.position > pattern.tokens.len() {
                    return Err(SuffixTreeError::InvalidLabelRange);
                }
                Ok(pattern.tokens[occurrence.position - 1..].to_vec())
            }
        }
    }

    /// Translate `LabelOfSuffixTreeEdge`.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown edge or invalid label range.
    pub fn label_of_edge(
        &self,
        edge: SuffixTreeEdgeId,
    ) -> Result<Vec<SuffixTreeToken<Symbol>>, SuffixTreeError> {
        self.edge_label(self.edge(edge)?.label)
    }

    /// Ordered outgoing edges of an internal node.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown or external node.
    pub fn node_edges(
        &self,
        node: SuffixTreeNodeId,
    ) -> Result<&[SuffixTreeEdgeId], SuffixTreeError> {
        match &self.node(node)?.representation {
            SuffixTreeNodeRepresentation::Internal { edges, .. } => Ok(edges),
            SuffixTreeNodeRepresentation::External { .. } => {
                Err(SuffixTreeError::ExpectedInternalNode)
            }
        }
    }

    /// Suffix link of an internal node (`None` represents QPA's `0`).
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown or external node.
    pub fn node_suffix_link(
        &self,
        node: SuffixTreeNodeId,
    ) -> Result<Option<SuffixTreeNodeId>, SuffixTreeError> {
        match &self.node(node)?.representation {
            SuffixTreeNodeRepresentation::Internal { suffix_link, .. } => Ok(*suffix_link),
            SuffixTreeNodeRepresentation::External { .. } => {
                Err(SuffixTreeError::ExpectedInternalNode)
            }
        }
    }

    /// QPA's marked flag for an internal node.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown or external node.
    pub fn node_is_marked(&self, node: SuffixTreeNodeId) -> Result<bool, SuffixTreeError> {
        match &self.node(node)?.representation {
            SuffixTreeNodeRepresentation::Internal { marked, .. } => Ok(*marked),
            SuffixTreeNodeRepresentation::External { .. } => {
                Err(SuffixTreeError::ExpectedInternalNode)
            }
        }
    }

    /// Pattern/start descriptors stored in an external node.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown or internal node.
    pub fn external_occurrences(
        &self,
        node: SuffixTreeNodeId,
    ) -> Result<&[SuffixTreeOccurrence], SuffixTreeError> {
        match &self.node(node)?.representation {
            SuffixTreeNodeRepresentation::External { occurrences } => Ok(occurrences),
            SuffixTreeNodeRepresentation::Internal { .. } => {
                Err(SuffixTreeError::ExpectedExternalNode)
            }
        }
    }

    /// Parent endpoint of an edge, or `None` after QPA detaches it.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown edge.
    pub fn edge_parent(
        &self,
        edge: SuffixTreeEdgeId,
    ) -> Result<Option<SuffixTreeNodeId>, SuffixTreeError> {
        Ok(self.edge(edge)?.parent)
    }

    /// Child endpoint of an edge, or `None` after QPA detaches it.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown edge.
    pub fn edge_child(
        &self,
        edge: SuffixTreeEdgeId,
    ) -> Result<Option<SuffixTreeNodeId>, SuffixTreeError> {
        Ok(self.edge(edge)?.child)
    }

    /// Translate `SuffixTreeEdgeStartingAt`, preserving edge insertion order.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown/external node or corrupt label.
    pub fn edge_starting_at(
        &self,
        node: SuffixTreeNodeId,
        token: &SuffixTreeToken<Symbol>,
    ) -> Result<Option<SuffixTreeEdgeId>, SuffixTreeError> {
        for edge in self.node_edges(node)? {
            let label = self.label_of_edge(*edge)?;
            if label.first() == Some(token) {
                return Ok(Some(*edge));
            }
        }
        Ok(None)
    }

    /// Translate `CreateInternalSuffixTreeNode` for an existing pattern slot.
    ///
    /// `start` and `end` use QPA's one-based inclusive convention.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown pattern or invalid label range.
    pub fn create_internal_node(
        &mut self,
        pattern_number: usize,
        start: usize,
        end: usize,
    ) -> Result<SuffixTreeNodeId, SuffixTreeError> {
        self.validate_internal_span(pattern_number, start, end)?;
        Ok(self.push_internal_node(Some(LabelSpan {
            pattern_number,
            start,
            end: Some(end),
        })))
    }

    /// Translate `CreateExternalSuffixTreeNode` for an existing pattern slot.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown pattern or invalid one-based start.
    pub fn create_external_node(
        &mut self,
        pattern_number: usize,
        start: usize,
    ) -> Result<SuffixTreeNodeId, SuffixTreeError> {
        self.validate_external_start(pattern_number, start)?;
        Ok(self.push_external_node(pattern_number, start))
    }

    /// Translate `CreateSuffixTreeEdge` and update both endpoint records.
    ///
    /// `end == None` represents QPA's `-1` terminal edge end.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid identities, endpoint categories, or label.
    pub fn create_edge(
        &mut self,
        parent: SuffixTreeNodeId,
        child: SuffixTreeNodeId,
        pattern_number: usize,
        start: usize,
        end: Option<usize>,
    ) -> Result<SuffixTreeEdgeId, SuffixTreeError> {
        self.node(parent)?;
        self.node(child)?;
        if self.node_kind(parent)? != SuffixTreeNodeKind::Internal {
            return Err(SuffixTreeError::ExpectedInternalNode);
        }
        self.validate_edge_span(pattern_number, start, end)?;
        Ok(self.push_edge(
            parent,
            child,
            LabelSpan {
                pattern_number,
                start,
                end,
            },
        ))
    }

    /// Translate `AddInternalSuffixTreeNode` by splitting an edge.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid identities, location, or node category.
    pub fn add_internal_node(
        &mut self,
        edge: SuffixTreeEdgeId,
        location: usize,
        node: SuffixTreeNodeId,
    ) -> Result<(), SuffixTreeError> {
        if location == 0 || self.node_kind(node)? != SuffixTreeNodeKind::Internal {
            return Err(SuffixTreeError::InvalidSplitLocation);
        }
        let old_edge = *self.edge(edge)?;
        let old_child = old_edge.child.ok_or(SuffixTreeError::DetachedEdge)?;
        let old_label_length = self.label_of_edge(edge)?.len();
        if location >= old_label_length {
            return Err(SuffixTreeError::InvalidSplitLocation);
        }

        self.edges[edge.index].child = Some(node);
        self.nodes[node.index].parent_edge = Some(edge);
        self.edges[edge.index].label.end = Some(old_edge.label.start + location - 1);
        self.push_edge(
            node,
            old_child,
            LabelSpan {
                pattern_number: old_edge.label.pattern_number,
                start: old_edge.label.start + location,
                end: old_edge.label.end,
            },
        );
        Ok(())
    }

    /// Translate `AddExternalSuffixTreeNode`.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid identities/categories or an invalid label.
    pub fn add_external_node(
        &mut self,
        parent: SuffixTreeNodeId,
        child: SuffixTreeNodeId,
    ) -> Result<SuffixTreeEdgeId, SuffixTreeError> {
        if self.node_kind(parent)? != SuffixTreeNodeKind::Internal {
            return Err(SuffixTreeError::ExpectedInternalNode);
        }
        let occurrence = self
            .external_occurrences(child)?
            .first()
            .copied()
            .ok_or(SuffixTreeError::EmptyExternalNode)?;
        let start = self.label_of_node(parent)?.len() + occurrence.position;
        self.validate_edge_span(occurrence.pattern_number, start, None)?;
        Ok(self.push_edge(
            parent,
            child,
            LabelSpan {
                pattern_number: occurrence.pattern_number,
                start,
                end: None,
            },
        ))
    }

    /// Translate `InsertPatternIntoSuffixTree` phase-for-phase.
    ///
    /// The returned pattern number and all internal positions are one-based,
    /// matching QPA.
    ///
    /// # Errors
    ///
    /// Returns a typed source failure for QPA's unbound-variable empty-pattern
    /// case or if a suffix-tree/dynamic-tree invariant is violated.
    #[allow(
        clippy::comparison_chain,
        clippy::if_not_else,
        clippy::many_single_char_names,
        clippy::too_many_lines
    )]
    pub fn insert_pattern(&mut self, pattern: Vec<Symbol>) -> Result<usize, SuffixTreeError> {
        if pattern.is_empty() {
            return Err(SuffixTreeError::PinnedEmptyPatternInsertion);
        }
        let pattern_number = self.patterns.len() + 1;
        let mut tokens = pattern
            .iter()
            .cloned()
            .map(SuffixTreeToken::Symbol)
            .collect::<Vec<_>>();
        tokens.push(SuffixTreeToken::Termination);
        self.patterns.push(StoredPattern {
            symbols: pattern,
            tokens,
            active: true,
        });

        let pattern_tokens = self.stored_pattern(pattern_number)?.tokens.clone();
        let m = pattern_tokens.len();
        let mut i = 0;
        let mut current_node = self.root;
        let mut current_edge = None;
        let mut location = 0;

        while i < m - 1 {
            location = 0;
            current_edge = self.edge_starting_at(current_node, &pattern_tokens[i])?;
            let Some(edge) = current_edge else {
                break;
            };
            let label = self.label_of_edge(edge)?;
            let next_node = self.edge_child_required(edge)?;
            i += 1;
            location = 1;

            while location < label.len() && i < m - 1 && pattern_tokens[i] == label[location] {
                i += 1;
                location += 1;
            }

            if location == label.len() {
                current_node = next_node;
                location = 0;
            } else {
                break;
            }
        }

        let mut j_i = 0;
        while i < m {
            let mut phase_done = false;
            let mut j = j_i + 1;
            let extension = if location != 0 {
                ExtensionLocus::Edge(
                    current_edge.ok_or(SuffixTreeError::InvalidSuffixTreeState)?,
                    location,
                )
            } else {
                ExtensionLocus::Node(current_node)
            };
            let extension = self.apply_extension(extension, pattern_number, i, j)?;
            let mut w = extension.internal.then_some(extension.current_node);
            let mut v = extension.current_node;
            if extension.new_node.is_none() {
                j_i = j - 1;
                current_node = extension
                    .next_node
                    .ok_or(SuffixTreeError::InvalidSuffixTreeState)?;
                current_edge = extension.current_edge;
                location = extension.location;
                phase_done = true;
            } else {
                j_i += 1;
                current_node = self.root;
                location = 0;
            }
            j += 1;

            while j <= i + 1 && !phase_done {
                let mut label;
                if self.node_suffix_link(v)?.is_none() && v != self.root {
                    let parent_edge = self.parent_edge_required(v)?;
                    label = self.label_of_edge(parent_edge)?;
                    v = self.edge_parent_required(parent_edge)?;
                } else {
                    label = Vec::new();
                }

                if v != self.root {
                    current_node = self
                        .node_suffix_link(v)?
                        .ok_or(SuffixTreeError::MissingSuffixLink)?;
                } else {
                    label = qpa_pattern_slice(&pattern_tokens, j, i);
                    current_node = v;
                }

                let mut remaining = label.len();
                let mut label_position = 0;
                location = 0;
                while remaining > 0 {
                    let edge = self
                        .edge_starting_at(current_node, &label[label_position])?
                        .ok_or(SuffixTreeError::MissingGuaranteedEdge)?;
                    current_edge = Some(edge);
                    let edge_length = self.label_of_edge(edge)?.len();
                    if edge_length < remaining {
                        label_position += edge_length;
                        current_node = self.edge_child_required(edge)?;
                    } else if edge_length == remaining {
                        current_node = self.edge_child_required(edge)?;
                        location = 0;
                    } else {
                        location = remaining;
                    }
                    remaining = remaining.saturating_sub(edge_length);
                }

                let extension = if location != 0 {
                    ExtensionLocus::Edge(
                        current_edge.ok_or(SuffixTreeError::InvalidSuffixTreeState)?,
                        location,
                    )
                } else {
                    ExtensionLocus::Node(current_node)
                };
                let extension = self.apply_extension(extension, pattern_number, i, j)?;
                let x = extension.internal.then_some(extension.current_node);
                v = extension.current_node;
                if extension.new_node.is_none() {
                    j_i = j - 1;
                    current_node = extension
                        .next_node
                        .ok_or(SuffixTreeError::InvalidSuffixTreeState)?;
                    current_edge = extension.current_edge;
                    location = extension.location;
                    phase_done = true;
                } else {
                    j_i += 1;
                    current_node = self.root;
                    location = 0;
                }
                if let Some(previous_internal) = w {
                    let actual_label = self.label_of_node(v)?;
                    let expected_label = qpa_pattern_slice(&pattern_tokens, j, i);
                    if actual_label != expected_label {
                        return Err(SuffixTreeError::IncorrectSuffixLink {
                            phase: i,
                            extension: j,
                            actual_length: actual_label.len(),
                            expected_length: expected_label.len(),
                        });
                    }
                    self.set_suffix_link(previous_internal, v)?;
                }
                w = x;
                j += 1;
            }
            i += 1;
        }

        Ok(pattern_number)
    }

    /// Translate `DeletePatternFromSuffixTree`, including QPA's dynamic-forest
    /// updates and unary-node edge compression.
    ///
    /// # Errors
    ///
    /// Returns an error for an unknown/already-deleted slot or a violated
    /// source invariant.
    #[allow(clippy::if_not_else, clippy::too_many_lines)]
    pub fn delete_pattern(&mut self, pattern_number: usize) -> Result<(), SuffixTreeError> {
        let stored = self.stored_pattern(pattern_number)?;
        if !stored.active {
            return Err(SuffixTreeError::PatternAlreadyDeleted);
        }
        let pattern = stored.tokens.clone();
        let m = pattern.len();
        let mut current_node = self.root;
        let mut label = pattern;

        for i in 1..=m {
            let mut remaining = label.len();
            let mut label_position = 0;
            while remaining > 0 {
                let current_edge = self
                    .edge_starting_at(current_node, &label[label_position])?
                    .ok_or(SuffixTreeError::MissingGuaranteedEdge)?;
                let edge_length = self.label_of_edge(current_edge)?.len();
                if edge_length < remaining {
                    label_position += edge_length;
                    current_node = self.edge_child_required(current_edge)?;
                } else {
                    current_node = self.edge_child_required(current_edge)?;
                }
                remaining = remaining.saturating_sub(edge_length);
            }

            let u = current_node;
            let v = self
                .parent_node(u)?
                .ok_or(SuffixTreeError::MissingParentNode)?;
            let children = self.qpa_child_count(v)?;
            self.remove_external_occurrence(u, v, pattern_number, i)?;

            if v != self.root && children <= 2 {
                let edges = self.node_edges(v)?.to_vec();
                if edges.len() != 1 {
                    return Err(SuffixTreeError::ExpectedUnaryInternalNode);
                }
                let w = self.edge_child_required(edges[0])?;
                self.delete_from_dynamic_forest(v, true, Some(w))?;
                let v_parent_edge = self.parent_edge_required(v)?;
                let w_parent_edge = self.parent_edge_required(w)?;
                let prefix_length = self.label_of_edge(v_parent_edge)?.len();
                let mut merged_label = self.edge(w_parent_edge)?.label;
                merged_label.start = merged_label
                    .start
                    .checked_sub(prefix_length)
                    .ok_or(SuffixTreeError::InvalidLabelRange)?;
                self.edges[v_parent_edge.index].label = merged_label;
                self.edges[v_parent_edge.index].child = Some(w);
                self.nodes[w.index].parent_edge = Some(v_parent_edge);
                self.nodes[v.index].parent_edge = None;
            }

            if i < m {
                if v == self.root {
                    current_node = v;
                    label = label.get(1..).unwrap_or_default().to_vec();
                } else {
                    current_node = self
                        .node_suffix_link(v)?
                        .ok_or(SuffixTreeError::MissingSuffixLink)?;
                    label = self.label_of_edge(self.parent_edge_required(u)?)?;
                }
            }
        }
        self.patterns[pattern_number - 1].active = false;
        Ok(())
    }

    /// Translate `AllPatternsInSequence` with QPA's clocus/beta traversal and
    /// marked-node dynamic forest.
    ///
    /// # Errors
    ///
    /// Returns an error if the maintained suffix/dynamic-tree invariants are
    /// violated.
    #[allow(clippy::if_not_else, clippy::too_many_lines)]
    pub fn all_patterns_in_sequence(
        &self,
        sequence: &[Symbol],
    ) -> Result<Vec<SuffixTreeOccurrence>, SuffixTreeError> {
        let sequence = sequence
            .iter()
            .cloned()
            .map(SuffixTreeToken::Symbol)
            .collect::<Vec<_>>();
        let n = sequence.len();
        let mut occurrences = Vec::new();
        let mut clocus = self.root;
        let mut beta = Vec::new();

        for j in 1..=n {
            if clocus != self.root {
                let mut x = self
                    .node_suffix_link(clocus)?
                    .ok_or(SuffixTreeError::MissingSuffixLink)?;
                let mut beta_head = beta.clone();
                let mut alpha = Vec::new();
                let mut current_edge = None;
                while alpha.len() < beta_head.len() {
                    beta_head = beta_head[alpha.len()..].to_vec();
                    let edge = self
                        .edge_starting_at(x, &beta_head[0])?
                        .ok_or(SuffixTreeError::MissingGuaranteedEdge)?;
                    current_edge = Some(edge);
                    alpha = self.label_of_edge(edge)?;
                    if alpha.len() < beta_head.len() {
                        x = self.edge_child_required(edge)?;
                    }
                }
                if alpha.len() > beta_head.len() {
                    clocus = x;
                    beta = beta_head;
                } else {
                    let y = if alpha.is_empty() {
                        x
                    } else {
                        self.edge_child_required(
                            current_edge.ok_or(SuffixTreeError::InvalidSuffixTreeState)?,
                        )?
                    };
                    let walked = self.walk_sequence_from(y, j, &sequence)?;
                    clocus = walked.0;
                    beta = self.beta_slice(j, clocus, walked.1, &sequence)?;
                }
            } else {
                let walked = self.walk_sequence_from(self.root, j, &sequence)?;
                clocus = walked.0;
                beta = self.beta_slice(j, clocus, walked.1, &sequence)?;
            }

            if !beta.is_empty()
                && let Some(edge) = self.edge_starting_at(clocus, &beta[0])?
            {
                let child = self.edge_child_required(edge)?;
                if self.node_kind(child)? == SuffixTreeNodeKind::External
                    && self.label_of_edge(edge)?.len() == beta.len() + 1
                {
                    self.append_pattern_starts_at(child, j, &mut occurrences)?;
                }
            }

            let mut u = clocus;
            while u != self.root {
                let dynamic_root = self
                    .dynamic_handle(u)?
                    .root()?
                    .item()
                    .ok_or(SuffixTreeError::InvalidSuffixTreeState)?;
                if dynamic_root != self.root {
                    let terminal_edge = self
                        .edge_starting_at(dynamic_root, &SuffixTreeToken::Termination)?
                        .ok_or(SuffixTreeError::MissingTerminalEdge)?;
                    let external = self.edge_child_required(terminal_edge)?;
                    self.append_pattern_starts_at(external, j, &mut occurrences)?;
                    u = self
                        .parent_node(dynamic_root)?
                        .ok_or(SuffixTreeError::MissingParentNode)?;
                } else {
                    u = self.root;
                }
            }
        }
        Ok(occurrences)
    }

    /// Translate `SequenceInPatterns`, preserving QPA's edge and leaf order.
    ///
    /// # Errors
    ///
    /// Returns an error if the suffix-tree arena is inconsistent.
    pub fn sequence_in_patterns(
        &self,
        sequence: &[Symbol],
    ) -> Result<Vec<SuffixTreeOccurrence>, SuffixTreeError> {
        let sequence = sequence
            .iter()
            .cloned()
            .map(SuffixTreeToken::Symbol)
            .collect::<Vec<_>>();
        let mut current_node = self.root;
        let mut current_edge = None;
        let mut current_label = Vec::new();
        let mut edge_length = 0;
        let mut index = 0;
        let mut location = 0;

        while index < sequence.len() {
            if location == 0 {
                let Some(edge) = self.edge_starting_at(current_node, &sequence[index])? else {
                    break;
                };
                current_edge = Some(edge);
                location = 1;
                current_label = self.label_of_edge(edge)?;
                edge_length = current_label.len();
            } else if sequence[index] == current_label[location] {
                location += 1;
            } else {
                break;
            }
            if location >= edge_length {
                current_node = self.edge_child_required(
                    current_edge.ok_or(SuffixTreeError::InvalidSuffixTreeState)?,
                )?;
                location = 0;
            }
            index += 1;
        }

        let mut occurrences = Vec::new();
        if index == sequence.len() {
            if location != 0 {
                current_node = self.edge_child_required(
                    current_edge.ok_or(SuffixTreeError::InvalidSuffixTreeState)?,
                )?;
            }
            self.collect_external_occurrences(current_node, &mut occurrences)?;
        }
        Ok(occurrences)
    }

    #[allow(clippy::if_not_else)]
    fn apply_extension(
        &mut self,
        locus: ExtensionLocus,
        pattern_number: usize,
        i: usize,
        j: usize,
    ) -> Result<ExtensionResult, SuffixTreeError> {
        let next_token = self.stored_pattern(pattern_number)?.tokens[i].clone();
        let mut result = ExtensionResult {
            internal: false,
            current_edge: None,
            location: 0,
            next_node: None,
            current_node: self.root,
            new_node: None,
        };
        let applied_edge;

        match locus {
            ExtensionLocus::Edge(current_edge, location) => {
                applied_edge = Some(current_edge);
                let label = self.label_of_edge(current_edge)?;
                if label
                    .get(location)
                    .ok_or(SuffixTreeError::InvalidSplitLocation)?
                    != &next_token
                {
                    let current_node = self.create_internal_node(pattern_number, j, i)?;
                    self.add_internal_node(current_edge, location, current_node)?;
                    let new_node = self.create_external_node(pattern_number, j)?;
                    self.add_external_node(current_node, new_node)?;
                    result.internal = true;
                    result.current_node = current_node;
                    result.new_node = Some(new_node);
                    if next_token == SuffixTreeToken::Termination {
                        let terminal_edge = self
                            .edge_starting_at(current_node, &SuffixTreeToken::Termination)?
                            .ok_or(SuffixTreeError::MissingTerminalEdge)?;
                        let terminal_child = self.edge_child_required(terminal_edge)?;
                        if self
                            .external_occurrences(terminal_child)?
                            .iter()
                            .any(|occurrence| occurrence.position == 1)
                        {
                            self.set_marked(current_node, true)?;
                        }
                    }
                    self.insert_into_dynamic_forest(current_node)?;
                    self.insert_into_dynamic_forest(new_node)?;
                } else {
                    let current_node = self.edge_parent_required(current_edge)?;
                    result.current_node = current_node;
                    if location + 1 < label.len() {
                        result.current_edge = Some(current_edge);
                        result.location = location + 1;
                        result.next_node = Some(current_node);
                    } else {
                        result.next_node = Some(self.edge_child_required(current_edge)?);
                    }
                }
            }
            ExtensionLocus::Node(current_node) => {
                result.current_node = current_node;
                let current_edge = self.edge_starting_at(current_node, &next_token)?;
                applied_edge = current_edge;
                if let Some(current_edge) = current_edge {
                    if self.label_of_edge(current_edge)?.len() > 1 {
                        result.current_edge = Some(current_edge);
                        result.location = 1;
                        result.next_node = Some(current_node);
                    } else {
                        result.next_node = Some(self.edge_child_required(current_edge)?);
                    }
                } else {
                    if self.node_kind(current_node)? == SuffixTreeNodeKind::External {
                        return Err(SuffixTreeError::ExpectedInternalNode);
                    }
                    let new_node = self.create_external_node(pattern_number, j)?;
                    self.add_external_node(current_node, new_node)?;
                    result.next_node = Some(current_node);
                    result.new_node = Some(new_node);
                    self.insert_into_dynamic_forest(new_node)?;
                }
            }
        }

        if next_token == SuffixTreeToken::Termination && result.new_node.is_none() {
            let current_edge = applied_edge.ok_or(SuffixTreeError::InvalidSuffixTreeState)?;
            let child = self.edge_child_required(current_edge)?;
            if self.node_kind(child)? != SuffixTreeNodeKind::External {
                return Err(SuffixTreeError::ExpectedExternalNode);
            }
            self.external_occurrences_mut(child)?
                .push(SuffixTreeOccurrence {
                    pattern_number,
                    position: j,
                });
            let parent = self.edge_parent_required(current_edge)?;
            if j == 1 && self.label_of_node(parent)?.len() == i {
                self.set_marked(parent, true)?;
                self.insert_into_dynamic_forest(parent)?;
            }
            result.new_node = Some(child);
        }

        Ok(result)
    }

    fn insert_into_dynamic_forest(
        &mut self,
        node: SuffixTreeNodeId,
    ) -> Result<(), SuffixTreeError> {
        let parent = self
            .parent_node(node)?
            .ok_or(SuffixTreeError::MissingParentNode)?;
        let parent_dynamic = self.dynamic_handle(parent)?;
        match self.node_kind(node)? {
            SuffixTreeNodeKind::External => {
                let dynamic = self.dynamic_forest.create(Some(node));
                self.nodes[node.index].dynamic_tree = Some(dynamic.clone());
                dynamic.link(&parent_dynamic)?;
            }
            SuffixTreeNodeKind::Internal => {
                let marked = self.node_is_marked(node)?;
                if marked && self.nodes[node.index].dynamic_tree.is_some() {
                    self.dynamic_handle(node)?.cut()?;
                    return Ok(());
                }
                if self.nodes[node.index].dynamic_tree.is_some() {
                    return Err(SuffixTreeError::InvalidSuffixTreeState);
                }
                let dynamic = self.dynamic_forest.create(Some(node));
                self.nodes[node.index].dynamic_tree = Some(dynamic.clone());
                let first_edge = *self
                    .node_edges(node)?
                    .first()
                    .ok_or(SuffixTreeError::ExpectedInternalChild)?;
                let old_child = self.edge_child_required(first_edge)?;
                let child_dynamic = self.dynamic_handle(old_child)?;
                let same_root = parent_dynamic.root()? == child_dynamic.root()?;
                if !marked {
                    if same_root {
                        child_dynamic.cut()?;
                        dynamic.link(&parent_dynamic)?;
                        child_dynamic.link(&dynamic)?;
                    } else {
                        dynamic.link(&parent_dynamic)?;
                    }
                } else if same_root {
                    child_dynamic.cut()?;
                    child_dynamic.link(&dynamic)?;
                }
            }
        }
        Ok(())
    }

    fn delete_from_dynamic_forest(
        &self,
        node: SuffixTreeNodeId,
        deleted: bool,
        replacement: Option<SuffixTreeNodeId>,
    ) -> Result<(), SuffixTreeError> {
        let parent = self
            .parent_node(node)?
            .ok_or(SuffixTreeError::MissingParentNode)?;
        if self.node_kind(node)? == SuffixTreeNodeKind::External {
            self.dynamic_handle(node)?.cut()?;
        } else if deleted {
            let replacement = replacement.ok_or(SuffixTreeError::ExpectedInternalChild)?;
            self.dynamic_handle(node)?.cut()?;
            self.dynamic_handle(replacement)?.cut()?;
            self.dynamic_handle(replacement)?
                .link(&self.dynamic_handle(parent)?)?;
        } else {
            self.dynamic_handle(node)?
                .link(&self.dynamic_handle(parent)?)?;
        }
        Ok(())
    }

    fn remove_external_occurrence(
        &mut self,
        node: SuffixTreeNodeId,
        parent: SuffixTreeNodeId,
        pattern_number: usize,
        phase: usize,
    ) -> Result<(), SuffixTreeError> {
        if self.node_kind(node)? != SuffixTreeNodeKind::External {
            return Err(SuffixTreeError::ExpectedExternalNode);
        }
        let parent_edge = self.parent_edge_required(node)?;
        let terminal_edge = self.label_of_edge(parent_edge)? == vec![SuffixTreeToken::Termination];
        if self.external_occurrences(node)?.len() > 1 {
            self.external_occurrences_mut(node)?
                .retain(|occurrence| occurrence.pattern_number != pattern_number);
        } else {
            self.delete_from_dynamic_forest(node, true, None)?;
            self.internal_edges_mut(parent)?
                .retain(|candidate| *candidate != parent_edge);
            self.edges[parent_edge.index].parent = None;
            self.edges[parent_edge.index].child = None;
        }
        if phase == 1 && terminal_edge && self.node_is_marked(parent)? {
            self.set_marked(parent, false)?;
            self.delete_from_dynamic_forest(parent, false, None)?;
        }
        Ok(())
    }

    fn walk_sequence_from(
        &self,
        start: SuffixTreeNodeId,
        query_start: usize,
        sequence: &[SuffixTreeToken<Symbol>],
    ) -> Result<(SuffixTreeNodeId, usize), SuffixTreeError> {
        let mut node = start;
        let mut last_node = node;
        let mut k = query_start + self.label_of_node(node)?.len();
        let mut current_edge = if k <= sequence.len() {
            self.edge_starting_at(node, &sequence[k - 1])?
        } else {
            None
        };

        while let Some(edge) = current_edge {
            node = self.edge_child_required(edge)?;
            let mut location = 2;
            k += 1;
            let label = self.label_of_edge(edge)?;
            while location <= label.len() && k <= sequence.len() {
                if label[location - 1] == sequence[k - 1] {
                    k += 1;
                    location += 1;
                } else {
                    break;
                }
            }
            if location > label.len() {
                last_node = node;
                current_edge = if k <= sequence.len() {
                    self.edge_starting_at(node, &sequence[k - 1])?
                } else {
                    None
                };
            } else {
                current_edge = None;
            }
        }
        Ok((last_node, k))
    }

    fn beta_slice(
        &self,
        query_start: usize,
        clocus: SuffixTreeNodeId,
        k: usize,
        sequence: &[SuffixTreeToken<Symbol>],
    ) -> Result<Vec<SuffixTreeToken<Symbol>>, SuffixTreeError> {
        let start = query_start + self.label_of_node(clocus)?.len();
        if start > sequence.len() {
            return Ok(Vec::new());
        }
        let start_index = start - 1;
        let end_index = k.saturating_sub(1).min(sequence.len());
        if start_index > end_index {
            Ok(Vec::new())
        } else {
            Ok(sequence[start_index..end_index].to_vec())
        }
    }

    fn append_pattern_starts_at(
        &self,
        external: SuffixTreeNodeId,
        query_position: usize,
        output: &mut Vec<SuffixTreeOccurrence>,
    ) -> Result<(), SuffixTreeError> {
        for occurrence in self.external_occurrences(external)? {
            if occurrence.position == 1 {
                output.push(SuffixTreeOccurrence {
                    pattern_number: occurrence.pattern_number,
                    position: query_position,
                });
            }
        }
        Ok(())
    }

    fn collect_external_occurrences(
        &self,
        node: SuffixTreeNodeId,
        output: &mut Vec<SuffixTreeOccurrence>,
    ) -> Result<(), SuffixTreeError> {
        match &self.node(node)?.representation {
            SuffixTreeNodeRepresentation::External { occurrences } => {
                output.extend(occurrences);
            }
            SuffixTreeNodeRepresentation::Internal { edges, .. } => {
                for edge in edges.clone() {
                    self.collect_external_occurrences(self.edge_child_required(edge)?, output)?;
                }
            }
        }
        Ok(())
    }

    fn qpa_child_count(&self, node: SuffixTreeNodeId) -> Result<usize, SuffixTreeError> {
        let mut count = 0;
        for edge in self.node_edges(node)? {
            let child = self.edge_child_required(*edge)?;
            count += match &self.node(child)?.representation {
                SuffixTreeNodeRepresentation::External { occurrences } => occurrences.len(),
                SuffixTreeNodeRepresentation::Internal { .. } => 1,
            };
        }
        Ok(count)
    }

    fn set_suffix_link(
        &mut self,
        node: SuffixTreeNodeId,
        suffix_link: SuffixTreeNodeId,
    ) -> Result<(), SuffixTreeError> {
        self.node(suffix_link)?;
        match &mut self.node_mut(node)?.representation {
            SuffixTreeNodeRepresentation::Internal {
                suffix_link: slot, ..
            } => {
                *slot = Some(suffix_link);
                Ok(())
            }
            SuffixTreeNodeRepresentation::External { .. } => {
                Err(SuffixTreeError::ExpectedInternalNode)
            }
        }
    }

    fn set_marked(&mut self, node: SuffixTreeNodeId, value: bool) -> Result<(), SuffixTreeError> {
        match &mut self.node_mut(node)?.representation {
            SuffixTreeNodeRepresentation::Internal { marked, .. } => {
                *marked = value;
                Ok(())
            }
            SuffixTreeNodeRepresentation::External { .. } => {
                Err(SuffixTreeError::ExpectedInternalNode)
            }
        }
    }

    fn internal_edges_mut(
        &mut self,
        node: SuffixTreeNodeId,
    ) -> Result<&mut Vec<SuffixTreeEdgeId>, SuffixTreeError> {
        match &mut self.node_mut(node)?.representation {
            SuffixTreeNodeRepresentation::Internal { edges, .. } => Ok(edges),
            SuffixTreeNodeRepresentation::External { .. } => {
                Err(SuffixTreeError::ExpectedInternalNode)
            }
        }
    }

    fn external_occurrences_mut(
        &mut self,
        node: SuffixTreeNodeId,
    ) -> Result<&mut Vec<SuffixTreeOccurrence>, SuffixTreeError> {
        match &mut self.node_mut(node)?.representation {
            SuffixTreeNodeRepresentation::External { occurrences } => Ok(occurrences),
            SuffixTreeNodeRepresentation::Internal { .. } => {
                Err(SuffixTreeError::ExpectedExternalNode)
            }
        }
    }

    fn push_internal_node(&mut self, label: Option<LabelSpan>) -> SuffixTreeNodeId {
        let node = SuffixTreeNodeId {
            owner: self.identity,
            index: self.nodes.len(),
        };
        self.nodes.push(SuffixTreeNode {
            parent_edge: None,
            dynamic_tree: None,
            representation: SuffixTreeNodeRepresentation::Internal {
                label,
                edges: Vec::new(),
                suffix_link: None,
                marked: false,
            },
        });
        node
    }

    fn push_external_node(&mut self, pattern_number: usize, start: usize) -> SuffixTreeNodeId {
        let node = SuffixTreeNodeId {
            owner: self.identity,
            index: self.nodes.len(),
        };
        self.nodes.push(SuffixTreeNode {
            parent_edge: None,
            dynamic_tree: None,
            representation: SuffixTreeNodeRepresentation::External {
                occurrences: vec![SuffixTreeOccurrence {
                    pattern_number,
                    position: start,
                }],
            },
        });
        node
    }

    fn push_edge(
        &mut self,
        parent: SuffixTreeNodeId,
        child: SuffixTreeNodeId,
        label: LabelSpan,
    ) -> SuffixTreeEdgeId {
        let edge = SuffixTreeEdgeId {
            owner: self.identity,
            index: self.edges.len(),
        };
        self.edges.push(SuffixTreeEdge {
            label,
            parent: Some(parent),
            child: Some(child),
        });
        match &mut self.nodes[parent.index].representation {
            SuffixTreeNodeRepresentation::Internal { edges, .. } => edges.push(edge),
            SuffixTreeNodeRepresentation::External { .. } => {
                unreachable!("validated suffix-tree edge parent")
            }
        }
        self.nodes[child.index].parent_edge = Some(edge);
        edge
    }

    fn dynamic_handle(
        &self,
        node: SuffixTreeNodeId,
    ) -> Result<DynamicTree<SuffixTreeNodeId>, SuffixTreeError> {
        self.node(node)?
            .dynamic_tree
            .clone()
            .ok_or(SuffixTreeError::MissingDynamicTree)
    }

    fn parent_edge_required(
        &self,
        node: SuffixTreeNodeId,
    ) -> Result<SuffixTreeEdgeId, SuffixTreeError> {
        self.node(node)?
            .parent_edge
            .ok_or(SuffixTreeError::MissingParentEdge)
    }

    fn edge_parent_required(
        &self,
        edge: SuffixTreeEdgeId,
    ) -> Result<SuffixTreeNodeId, SuffixTreeError> {
        self.edge(edge)?.parent.ok_or(SuffixTreeError::DetachedEdge)
    }

    fn edge_child_required(
        &self,
        edge: SuffixTreeEdgeId,
    ) -> Result<SuffixTreeNodeId, SuffixTreeError> {
        self.edge(edge)?.child.ok_or(SuffixTreeError::DetachedEdge)
    }

    fn stored_pattern(
        &self,
        pattern_number: usize,
    ) -> Result<&StoredPattern<Symbol>, SuffixTreeError> {
        pattern_number
            .checked_sub(1)
            .and_then(|index| self.patterns.get(index))
            .ok_or(SuffixTreeError::UnknownPattern)
    }

    fn node(&self, node: SuffixTreeNodeId) -> Result<&SuffixTreeNode, SuffixTreeError> {
        if node.owner != self.identity {
            return Err(SuffixTreeError::UnknownNode);
        }
        self.nodes
            .get(node.index)
            .ok_or(SuffixTreeError::UnknownNode)
    }

    fn node_mut(&mut self, node: SuffixTreeNodeId) -> Result<&mut SuffixTreeNode, SuffixTreeError> {
        if node.owner != self.identity {
            return Err(SuffixTreeError::UnknownNode);
        }
        self.nodes
            .get_mut(node.index)
            .ok_or(SuffixTreeError::UnknownNode)
    }

    fn edge(&self, edge: SuffixTreeEdgeId) -> Result<&SuffixTreeEdge, SuffixTreeError> {
        if edge.owner != self.identity {
            return Err(SuffixTreeError::UnknownEdge);
        }
        self.edges
            .get(edge.index)
            .ok_or(SuffixTreeError::UnknownEdge)
    }

    fn internal_label(
        &self,
        label: LabelSpan,
    ) -> Result<Vec<SuffixTreeToken<Symbol>>, SuffixTreeError> {
        let pattern = self.stored_pattern(label.pattern_number)?;
        let end = label.end.ok_or(SuffixTreeError::InvalidLabelRange)?;
        if label.start == 0 || end > pattern.tokens.len() {
            return Err(SuffixTreeError::InvalidLabelRange);
        }
        if end < label.start {
            return Ok(Vec::new());
        }
        Ok(pattern.tokens[label.start - 1..end].to_vec())
    }

    fn edge_label(
        &self,
        label: LabelSpan,
    ) -> Result<Vec<SuffixTreeToken<Symbol>>, SuffixTreeError> {
        let pattern = self.stored_pattern(label.pattern_number)?;
        let end = match label.end {
            Some(end) if end >= label.start => end,
            Some(_) | None => pattern.tokens.len(),
        };
        if label.start == 0 || label.start > pattern.tokens.len() || end > pattern.tokens.len() {
            return Err(SuffixTreeError::InvalidLabelRange);
        }
        Ok(pattern.tokens[label.start - 1..end].to_vec())
    }

    fn validate_internal_span(
        &self,
        pattern_number: usize,
        start: usize,
        end: usize,
    ) -> Result<(), SuffixTreeError> {
        let pattern = self.stored_pattern(pattern_number)?;
        if start == 0 || start > pattern.tokens.len() + 1 || end > pattern.tokens.len() {
            Err(SuffixTreeError::InvalidLabelRange)
        } else {
            Ok(())
        }
    }

    fn validate_external_start(
        &self,
        pattern_number: usize,
        start: usize,
    ) -> Result<(), SuffixTreeError> {
        let pattern = self.stored_pattern(pattern_number)?;
        if start == 0 || start > pattern.tokens.len() {
            Err(SuffixTreeError::InvalidLabelRange)
        } else {
            Ok(())
        }
    }

    fn validate_edge_span(
        &self,
        pattern_number: usize,
        start: usize,
        end: Option<usize>,
    ) -> Result<(), SuffixTreeError> {
        let pattern = self.stored_pattern(pattern_number)?;
        if start == 0
            || start > pattern.tokens.len()
            || end.is_some_and(|end| end > pattern.tokens.len())
        {
            Err(SuffixTreeError::InvalidLabelRange)
        } else {
            Ok(())
        }
    }
}

impl<Symbol> Default for SuffixTree<Symbol>
where
    Symbol: Clone + Eq,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<Symbol> fmt::Debug for SuffixTree<Symbol> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SuffixTree")
            .field("patterns", &self.patterns.len())
            .field("nodes", &self.nodes.len())
            .field("edges", &self.edges.len())
            .finish_non_exhaustive()
    }
}

fn qpa_pattern_slice<Symbol: Clone>(
    pattern: &[SuffixTreeToken<Symbol>],
    start: usize,
    end: usize,
) -> Vec<SuffixTreeToken<Symbol>> {
    if start == 0 || start > end || end > pattern.len() {
        Vec::new()
    } else {
        pattern[start - 1..end].to_vec()
    }
}

/// Typed failures replacing QPA's `fail`, assertions, and record errors.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum SuffixTreeError {
    /// No retained pattern has this one-based insertion number.
    #[error("the suffix-tree pattern number does not exist")]
    UnknownPattern,
    /// The pattern slot was already removed from the active tree.
    #[error("the suffix-tree pattern has already been deleted")]
    PatternAlreadyDeleted,
    /// No arena node has this identity.
    #[error("the suffix-tree node does not exist")]
    UnknownNode,
    /// No arena edge has this identity.
    #[error("the suffix-tree edge does not exist")]
    UnknownEdge,
    /// An operation requiring an internal node received an external leaf.
    #[error("the suffix-tree node must be internal")]
    ExpectedInternalNode,
    /// An operation requiring an external leaf received an internal node.
    #[error("the suffix-tree node must be external")]
    ExpectedExternalNode,
    /// QPA attempted to read a missing `parentEdge` component.
    #[error("the suffix-tree node has no parent edge")]
    MissingParentEdge,
    /// QPA attempted to use a node without a parent.
    #[error("the suffix-tree node has no parent node")]
    MissingParentNode,
    /// A detached edge was unexpectedly traversed.
    #[error("the suffix-tree edge is detached")]
    DetachedEdge,
    /// A QPA one-based label descriptor is malformed.
    #[error("the suffix-tree label range is invalid")]
    InvalidLabelRange,
    /// Edge splitting requires a positive proper-prefix location.
    #[error("the suffix-tree split location is invalid")]
    InvalidSplitLocation,
    /// An external leaf unexpectedly has no pattern/start descriptor.
    #[error("the external suffix-tree node has no occurrences")]
    EmptyExternalNode,
    /// An internal insertion required the preserved first child edge.
    #[error("the internal suffix-tree node has no preserved child")]
    ExpectedInternalChild,
    /// QPA's guaranteed suffix-tree path could not be followed.
    #[error("a guaranteed suffix-tree edge is missing")]
    MissingGuaranteedEdge,
    /// An internal node required a source-computed suffix link.
    #[error("a required suffix link is missing")]
    MissingSuffixLink,
    /// A marked pattern node required its terminating edge.
    #[error("a required terminating suffix-tree edge is missing")]
    MissingTerminalEdge,
    /// A suffix-tree node was not inserted into the dynamic forest.
    #[error("the suffix-tree node has no dynamic-tree handle")]
    MissingDynamicTree,
    /// QPA's suffix-link assertion failed.
    #[error(
        "QPA computed an incorrect suffix link in phase {phase}, extension {extension}: label lengths {actual_length} and {expected_length} differ"
    )]
    IncorrectSuffixLink {
        /// Zero-based source phase index.
        phase: usize,
        /// One-based source extension index.
        extension: usize,
        /// Length of the reached node label.
        actual_length: usize,
        /// Length of QPA's expected pattern slice.
        expected_length: usize,
    },
    /// QPA's deletion compression expected exactly one remaining edge.
    #[error("QPA expected a unary internal suffix-tree node")]
    ExpectedUnaryInternalNode,
    /// The pinned QPA method reads an unbound `currentEdge` for this input.
    #[error("pinned QPA cannot insert an empty suffix-tree pattern")]
    PinnedEmptyPatternInsertion,
    /// A translated source invariant was violated.
    #[error("the suffix tree is in an invalid source state")]
    InvalidSuffixTreeState,
    /// QPA's translated dynamic-tree operation failed.
    #[error(transparent)]
    DynamicTree(#[from] DynamicTreeError),
}

#[cfg(test)]
mod tests {
    use serde_json::{Value, json};

    use super::{
        INFO_DYNAMIC_DICTIONARY, SuffixTree, SuffixTreeError, SuffixTreeNodeId, SuffixTreeNodeKind,
        SuffixTreeOccurrence, SuffixTreeToken,
    };

    fn token_json(token: &SuffixTreeToken<i32>) -> Value {
        match token {
            SuffixTreeToken::Symbol(symbol) => json!(symbol),
            SuffixTreeToken::Termination => json!([]),
        }
    }

    fn token_list_json(tokens: &[SuffixTreeToken<i32>]) -> Value {
        Value::Array(tokens.iter().map(token_json).collect())
    }

    fn occurrences_json(occurrences: &[SuffixTreeOccurrence]) -> Value {
        json!(
            occurrences
                .iter()
                .map(|occurrence| [occurrence.pattern_number, occurrence.position])
                .collect::<Vec<_>>()
        )
    }

    fn node_json(tree: &SuffixTree<i32>, node: SuffixTreeNodeId) -> Value {
        let label = token_list_json(&tree.label_of_node(node).unwrap());
        match tree.node_kind(node).unwrap() {
            SuffixTreeNodeKind::External => json!({
                "kind": "external",
                "label": label,
                "occurrences": occurrences_json(tree.external_occurrences(node).unwrap()),
            }),
            SuffixTreeNodeKind::Internal => {
                let suffix = tree
                    .node_suffix_link(node)
                    .unwrap()
                    .map_or(Value::Null, |suffix| {
                        token_list_json(&tree.label_of_node(suffix).unwrap())
                    });
                let edges = tree
                    .node_edges(node)
                    .unwrap()
                    .iter()
                    .map(|edge| {
                        json!({
                            "label": token_list_json(&tree.label_of_edge(*edge).unwrap()),
                            "child": node_json(tree, tree.edge_child(*edge).unwrap().unwrap()),
                        })
                    })
                    .collect::<Vec<_>>();
                json!({
                    "kind": "internal",
                    "label": label,
                    "marked": tree.node_is_marked(node).unwrap(),
                    "suffix": suffix,
                    "edges": edges,
                })
            }
        }
    }

    fn queries_json(tree: &SuffixTree<i32>) -> Value {
        json!({
            "all_in": occurrences_json(
                &tree
                    .all_patterns_in_sequence(&[0, 1, 2, 1, 2, 1, 3])
                    .unwrap()
            ),
            "containing_12": occurrences_json(&tree.sequence_in_patterns(&[1, 2]).unwrap()),
            "containing_21": occurrences_json(&tree.sequence_in_patterns(&[2, 1]).unwrap()),
            "containing_empty": occurrences_json(&tree.sequence_in_patterns(&[]).unwrap()),
        })
    }

    #[test]
    fn pinned_qpa_suffix_tree_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/suffix_tree.schema-v1.json"
        ))
        .unwrap();
        let mut tree = SuffixTree::new();
        let patterns = [vec![1, 2, 1], vec![2, 1], vec![1, 2], vec![3]];
        let indices = patterns
            .iter()
            .cloned()
            .map(|pattern| tree.insert_pattern(pattern).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(json!(indices), expected["indices"]);
        assert_eq!(
            json!(
                (1..=tree.number_of_patterns())
                    .map(|pattern| tree.pattern(pattern).unwrap())
                    .collect::<Vec<_>>()
            ),
            expected["patterns"]
        );
        assert_eq!(queries_json(&tree), expected["before_delete"]);
        assert_eq!(
            node_json(&tree, tree.root_node()),
            expected["tree_before_delete"]
        );

        tree.delete_pattern(2).unwrap();
        assert_eq!(queries_json(&tree), expected["after_delete_2"]);
        tree.delete_pattern(1).unwrap();
        assert_eq!(queries_json(&tree), expected["after_delete_1"]);
        tree.delete_pattern(3).unwrap();
        assert_eq!(queries_json(&tree), expected["after_delete_3"]);
        assert_eq!(
            node_json(&tree, tree.root_node()),
            expected["tree_after_delete"]
        );
        assert!(!tree.is_pattern_active(1).unwrap());
        assert!(tree.is_pattern_active(4).unwrap());
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_suffix_tree_surface_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/suffix_tree_surface.schema-v1.json"
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

        let mut tree = SuffixTree::new();
        assert_eq!(tree.insert_pattern(vec![1, 2, 3]), Ok(1));
        let root = tree.root_node();
        let internal = tree.create_internal_node(1, 1, 1).unwrap();
        let external = tree.create_external_node(1, 2).unwrap();
        let edge = tree.create_edge(internal, external, 1, 2, None).unwrap();

        let split_parent = tree.create_internal_node(1, 1, 1).unwrap();
        let split_child = tree.create_external_node(1, 1).unwrap();
        let split_edge = tree
            .create_edge(split_parent, split_child, 1, 1, Some(2))
            .unwrap();
        let split_node = tree.create_internal_node(1, 1, 1).unwrap();
        tree.add_internal_node(split_edge, 1, split_node).unwrap();
        let split_remainder = tree.node_edges(split_node).unwrap()[0];

        let added_external = tree.create_external_node(1, 1).unwrap();
        let added_edge = tree.add_external_node(root, added_external).unwrap();

        let mut other_tree = SuffixTree::new();
        assert_eq!(other_tree.insert_pattern(vec![1, 2, 3]), Ok(1));
        let other_root = other_tree.root_node();
        let other_edge = other_tree.node_edges(other_root).unwrap()[0];

        assert_eq!(
            json!(INFO_DYNAMIC_DICTIONARY == "InfoDynamicDictionary"),
            expected["info_class"]
        );
        assert_eq!(json!([true, true]), expected["tree_filters"]);
        assert_eq!(
            json!([
                true,
                true,
                tree.node_kind(root).unwrap() == SuffixTreeNodeKind::Internal,
                tree.node_kind(root).unwrap() == SuffixTreeNodeKind::External,
                true,
                tree.node_kind(root).unwrap() == SuffixTreeNodeKind::Internal
            ]),
            expected["root_filters"]
        );
        assert_eq!(
            json!([
                true,
                true,
                tree.node_kind(external).unwrap() == SuffixTreeNodeKind::Internal,
                tree.node_kind(external).unwrap() == SuffixTreeNodeKind::External,
                true
            ]),
            expected["external_filters"]
        );
        assert_eq!(json!([true, true, true]), expected["edge_filters"]);
        assert_eq!(
            json!({
                "dead_length": tree.dead_length(),
                "pattern_count": tree.number_of_patterns(),
                "pattern_1": tree.pattern(1).unwrap()
            }),
            expected["components"]
        );
        assert_eq!(
            json!({
                "internal_label": token_list_json(&tree.label_of_node(internal).unwrap()),
                "external_label": token_list_json(&tree.label_of_node(external).unwrap()),
                "edge_label": token_list_json(&tree.label_of_edge(edge).unwrap()),
                "root_parent_is_zero": tree.parent_node(root).unwrap().is_none(),
                "child_parent_identity": tree.parent_node(external).unwrap() == Some(internal),
                "edge_lookup_identity": tree.edge_starting_at(internal, &SuffixTreeToken::Symbol(2)).unwrap() == Some(edge),
                "missing_lookup_is_fail": tree.edge_starting_at(internal, &SuffixTreeToken::Symbol(9)).unwrap().is_none(),
                "self_equal": edge == edge,
                "distinct_equal": edge == split_edge,
                "other_tree_root_equal": root == other_root,
                "other_tree_edge_equal": edge == other_edge
            }),
            expected["manual_edge"]
        );
        assert_eq!(
            json!({
                "prefix": token_list_json(&tree.label_of_edge(split_edge).unwrap()),
                "remainder": token_list_json(&tree.label_of_edge(split_remainder).unwrap()),
                "old_child_parent_identity": tree.parent_node(split_child).unwrap() == Some(split_node)
            }),
            expected["split"]
        );
        assert_eq!(
            json!({
                "edge_label": token_list_json(&tree.label_of_edge(added_edge).unwrap()),
                "parent_identity": tree.parent_node(added_external).unwrap() == Some(root)
            }),
            expected["add_external"]
        );

        assert_eq!(
            tree.parent_node(other_root),
            Err(SuffixTreeError::UnknownNode)
        );
        assert_eq!(
            tree.label_of_edge(other_edge),
            Err(SuffixTreeError::UnknownEdge)
        );
    }

    #[test]
    fn qpa_source_failures_and_retained_slots_are_explicit() {
        let mut tree = SuffixTree::<i32>::new();
        assert_eq!(tree.dead_length(), 0);
        assert_eq!(
            tree.insert_pattern(Vec::new()),
            Err(SuffixTreeError::PinnedEmptyPatternInsertion)
        );
        assert_eq!(tree.pattern(1), Err(SuffixTreeError::UnknownPattern));

        let pattern = tree.insert_pattern(vec![1]).unwrap();
        assert_eq!(tree.pattern(pattern).unwrap(), &[1]);
        tree.delete_pattern(pattern).unwrap();
        assert_eq!(tree.pattern(pattern).unwrap(), &[1]);
        assert_eq!(
            tree.delete_pattern(pattern),
            Err(SuffixTreeError::PatternAlreadyDeleted)
        );
    }

    #[test]
    fn disjoint_quadratic_patterns_support_qpa_gentle_relations() {
        let mut tree = SuffixTree::new();
        for pattern in [
            vec![0, 2],
            vec![2, 3],
            vec![3, 0],
            vec![4, 5],
            vec![5, 6],
            vec![6, 7],
            vec![7, 8],
        ] {
            tree.insert_pattern(pattern).unwrap();
        }
    }

    #[test]
    fn cubic_radical_products_support_qpa_gentle_relations() {
        let mut tree = SuffixTree::new();
        for pattern in [
            vec![1, 3, 4],
            vec![1, 3, 7],
            vec![2, 1, 3],
            vec![3, 4, 1],
            vec![3, 7, 8],
            vec![4, 1, 3],
            vec![5, 6, 4],
            vec![5, 6, 7],
            vec![6, 4, 1],
            vec![6, 7, 8],
            vec![7, 8, 6],
            vec![7, 8, 9],
            vec![8, 6, 7],
            vec![8, 9, 8],
            vec![9, 8, 9],
        ] {
            if let Err(error) = tree.insert_pattern(pattern.clone()) {
                panic!("failed to insert {pattern:?}: {error:?}");
            }
        }
    }

    #[test]
    fn pinned_qpa_duplicate_deletion_failure_is_preserved() {
        let mut tree = SuffixTree::new();
        for pattern in [
            vec![1, 2, 1],
            vec![2, 1],
            vec![1, 2],
            vec![3],
            vec![1, 2, 1],
        ] {
            tree.insert_pattern(pattern).unwrap();
        }
        tree.delete_pattern(2).unwrap();
        tree.delete_pattern(1).unwrap();
        assert_eq!(
            tree.all_patterns_in_sequence(&[0, 1, 2, 1, 2, 1, 3]),
            Err(SuffixTreeError::MissingGuaranteedEdge)
        );
    }
}
