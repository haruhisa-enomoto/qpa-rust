//! QPA's generic and quiver-specific static dictionaries.
//!
//! This is a direct Rust model of `lib/groebnerstadict.gi`. QPA builds an
//! Aho--Corasick-style trie from Groebner tips, completes it with quiver
//! transitions, detects cycles in the nonaccepting automaton, and uses the
//! resulting DAG to count or enumerate nontip paths.

use std::collections::VecDeque;

use numerica::domains::integer::Integer;
use thiserror::Error;

use crate::{ArrowIndex, Path, PathEndpoint, Quiver, QuiverError, VertexIndex};

/// QPA's first level-one `InfoStaticDictionary` line.
pub const STATIC_DICTIONARY_CREATING_AUTOMATA_INFO: &str = "#I  Creating finite automata...\n";

/// QPA's second level-one `InfoStaticDictionary` line.
pub const STATIC_DICTIONARY_CHECKING_CYCLE_INFO: &str = "#I  Checking for cycle...\n";

#[derive(Clone, Debug)]
struct DictionaryState<Symbol> {
    labels: Vec<Symbol>,
    children: Vec<usize>,
    parent: Option<usize>,
    label: Option<Symbol>,
    failure: usize,
    output: Option<usize>,
    patterns: Vec<usize>,
}

impl<Symbol> DictionaryState<Symbol> {
    fn root() -> Self {
        Self {
            labels: Vec::new(),
            children: Vec::new(),
            parent: None,
            label: None,
            failure: 0,
            output: None,
            patterns: Vec::new(),
        }
    }

    fn child(parent: usize, label: Symbol) -> Self {
        Self {
            labels: Vec::new(),
            children: Vec::new(),
            parent: Some(parent),
            label: Some(label),
            failure: 0,
            output: None,
            patterns: Vec::new(),
        }
    }
}

/// QPA's Aho--Corasick static dictionary for arbitrary dense sequences.
#[derive(Clone, Debug)]
pub struct StaticDictionary<Symbol> {
    patterns: Vec<Vec<Symbol>>,
    states: Vec<DictionaryState<Symbol>>,
}

impl<Symbol: Clone + Eq> StaticDictionary<Symbol> {
    /// Translate QPA's `StaticDictionary(patterns)` constructor.
    #[must_use]
    pub fn new(patterns: &[Vec<Symbol>]) -> Self {
        let mut states = vec![DictionaryState::root()];
        for (pattern_index, pattern) in patterns.iter().enumerate() {
            let mut current = 0;
            for symbol in pattern {
                let next = states[current]
                    .labels
                    .iter()
                    .position(|label| label == symbol)
                    .map(|position| states[current].children[position]);
                current = next.unwrap_or_else(|| {
                    let next = states.len();
                    states.push(DictionaryState::child(current, symbol.clone()));
                    states[current].labels.push(symbol.clone());
                    states[current].children.push(next);
                    next
                });
            }
            states[current].patterns.push(pattern_index + 1);
        }

        let mut queue = VecDeque::new();
        let root_children = states[0].children.clone();
        for child in root_children {
            states[child].failure = 0;
            queue.extend(states[child].children.iter().copied());
        }
        while let Some(current) = queue.pop_front() {
            let Some(parent) = states[current].parent else {
                continue;
            };
            let Some(symbol) = states[current].label.clone() else {
                continue;
            };
            let mut fallback = states[parent].failure;
            while fallback != 0 && !states[fallback].labels.iter().any(|label| label == &symbol) {
                fallback = states[fallback].failure;
            }
            states[current].failure = states[fallback]
                .labels
                .iter()
                .position(|label| label == &symbol)
                .map_or(0, |position| states[fallback].children[position]);
            let failure = states[current].failure;
            states[current].output = if states[failure].patterns.is_empty() {
                states[failure].output
            } else {
                Some(failure)
            };
            let children = states[current].children.clone();
            queue.extend(children);
        }

        Self {
            patterns: patterns.to_vec(),
            states,
        }
    }

    /// Return QPA's immutable `Patterns` attribute.
    #[must_use]
    pub fn patterns(&self) -> &[Vec<Symbol>] {
        &self.patterns
    }

    /// Translate QPA's dense-list `Search` method.
    ///
    /// # Errors
    ///
    /// Returns [`StaticDictionaryError::EmptyText`] for the source method's
    /// initial access to the first element of an empty dense list.
    pub fn search(&self, text: &[Symbol]) -> Result<Vec<DictionaryMatch>, StaticDictionaryError> {
        if text.is_empty() {
            return Err(StaticDictionaryError::EmptyText);
        }
        let mut position = 0;
        let mut state = 0;
        let mut result = Vec::new();
        loop {
            let mut transition = self.states[state]
                .labels
                .iter()
                .position(|label| label == &text[position]);
            while let Some(index) = transition {
                let next = self.states[state].children[index];
                let mut patterns = self.states[next].patterns.clone();
                let mut output = self.states[next].output;
                while let Some(output_state) = output {
                    patterns.extend(self.states[output_state].patterns.iter().copied());
                    output = self.states[output_state].output;
                }
                if !patterns.is_empty() {
                    result.push(DictionaryMatch {
                        end_position: position + 1,
                        pattern_indices: patterns,
                    });
                }
                state = next;
                position += 1;
                if position < text.len() {
                    transition = self.states[state]
                        .labels
                        .iter()
                        .position(|label| label == &text[position]);
                } else {
                    transition = None;
                }
            }
            if state == 0 {
                position += 1;
            }
            state = self.states[state].failure;
            if position >= text.len() {
                return Ok(result);
            }
        }
    }
}

impl<Symbol> std::fmt::Display for StaticDictionary<Symbol> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "<static dictionary containing {} pattern{}>",
            self.patterns.len(),
            if self.patterns.len() > 1 { "s" } else { "" }
        )
    }
}

#[derive(Clone, Debug)]
struct AutomatonState {
    target: VertexIndex,
    trie_transitions: Vec<Option<usize>>,
    transitions: Vec<Option<usize>>,
    failure: Option<usize>,
    own_patterns: Vec<usize>,
    outputs: Vec<usize>,
}

impl AutomatonState {
    fn new(target: VertexIndex, number_of_arrows: usize) -> Self {
        Self {
            target,
            trie_transitions: vec![None; number_of_arrows],
            transitions: vec![None; number_of_arrows],
            failure: None,
            own_patterns: Vec::new(),
            outputs: Vec::new(),
        }
    }

    fn is_tip(&self) -> bool {
        !self.own_patterns.is_empty()
    }
}

/// One QPA static-dictionary search match.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DictionaryMatch {
    /// One-based arrow position where the match ends; vertex-prefix matches use
    /// zero, matching QPA's `PrefixSearch`.
    pub end_position: usize,
    /// One-based indices into [`QuiverStaticDictionary::patterns`].
    pub pattern_indices: Vec<usize>,
}

/// A QPA quiver static dictionary built from forbidden path patterns.
#[derive(Clone, Debug)]
pub struct QuiverStaticDictionary {
    quiver: Quiver,
    patterns: Vec<Path>,
    states: Vec<AutomatonState>,
    vertex_states: Vec<usize>,
}

impl QuiverStaticDictionary {
    /// Translate QPA's `QuiverStaticDictionary(quiver, patterns)`.
    ///
    /// # Errors
    ///
    /// Returns an error for a pattern from another quiver or QPA's zero path.
    #[allow(clippy::too_many_lines)]
    pub fn new(quiver: &Quiver, patterns: &[Path]) -> Result<Self, StaticDictionaryError> {
        let number_of_arrows = quiver.number_of_arrows();
        let mut states = Vec::with_capacity(quiver.number_of_vertices() + patterns.len());
        let mut vertex_states = Vec::with_capacity(quiver.number_of_vertices());
        for vertex in quiver.vertices() {
            vertex_states.push(states.len());
            states.push(AutomatonState::new(vertex.index(), number_of_arrows));
        }

        for (pattern_index, pattern) in patterns.iter().enumerate() {
            if !pattern.quiver().is_identical(quiver) {
                return Err(StaticDictionaryError::DifferentQuivers);
            }
            let source = endpoint_index(pattern.source())?;
            let walk = pattern.walk().ok_or(StaticDictionaryError::ZeroPattern)?;
            let mut current = vertex_states[source.zero_based()];
            for arrow_index in walk {
                let position = arrow_index.zero_based();
                let next = if let Some(next) = states[current].trie_transitions[position] {
                    next
                } else {
                    let arrow = quiver
                        .arrow(*arrow_index)
                        .ok_or(StaticDictionaryError::InvalidPattern)?;
                    let next = states.len();
                    states.push(AutomatonState::new(
                        arrow.target().index(),
                        number_of_arrows,
                    ));
                    states[current].trie_transitions[position] = Some(next);
                    next
                };
                current = next;
            }
            states[current].own_patterns.push(pattern_index + 1);
        }

        let mut queue = VecDeque::new();
        for &vertex_state in &vertex_states {
            let children = states[vertex_state]
                .trie_transitions
                .iter()
                .enumerate()
                .filter_map(|(arrow, child)| child.map(|child| (arrow, child)))
                .collect::<Vec<_>>();
            for (arrow, child) in children {
                let target = quiver
                    .arrow(ArrowIndex::from_zero_based(arrow))
                    .ok_or(StaticDictionaryError::InvalidPattern)?
                    .target()
                    .index();
                states[child].failure = Some(vertex_states[target.zero_based()]);
                queue.push_back(child);
            }
        }

        while let Some(current) = queue.pop_front() {
            let failure_outputs = states[current]
                .failure
                .map(|failure| states[failure].outputs.clone())
                .unwrap_or_default();
            let own_patterns = states[current].own_patterns.clone();
            states[current].outputs.clone_from(&own_patterns);
            states[current].outputs.extend(failure_outputs);

            let children = states[current]
                .trie_transitions
                .iter()
                .enumerate()
                .filter_map(|(arrow, child)| child.map(|child| (arrow, child)))
                .collect::<Vec<_>>();
            for (arrow, child) in children {
                let failure = failure_transition(
                    &states,
                    &vertex_states,
                    quiver,
                    states[current].failure,
                    ArrowIndex::from_zero_based(arrow),
                )?;
                states[child].failure = Some(failure);
                queue.push_back(child);
            }
        }

        for &vertex_state in &vertex_states {
            let own_patterns = states[vertex_state].own_patterns.clone();
            states[vertex_state].outputs.clone_from(&own_patterns);
        }
        complete_quiver_transitions(quiver, &vertex_states, &mut states)?;

        Ok(Self {
            quiver: quiver.clone(),
            patterns: patterns.to_vec(),
            states,
            vertex_states,
        })
    }

    /// Return QPA's immutable pattern list.
    #[must_use]
    pub fn patterns(&self) -> &[Path] {
        &self.patterns
    }

    /// Return QPA's exact `ViewObj` text.
    #[must_use]
    pub fn view_string(&self) -> String {
        format!(
            "<static dictionary for quivers containing {} path{}>",
            self.patterns.len(),
            if self.patterns.len() > 1 { "s" } else { "" }
        )
    }

    /// Translate QPA's `Search(dict, path)`.
    ///
    /// # Errors
    ///
    /// Returns an error for a path from another quiver or the zero path.
    pub fn search(&self, path: &Path) -> Result<Vec<DictionaryMatch>, StaticDictionaryError> {
        self.validate_path(path)?;
        let source = endpoint_index(path.source())?;
        let walk = path.walk().ok_or(StaticDictionaryError::ZeroPattern)?;
        let mut current = self.vertex_states[source.zero_based()];
        let mut matches = Vec::new();
        let mut position = 0;
        while position < walk.len() {
            let arrow = walk[position];
            if let Some(next) = self.states[current].transitions[arrow.zero_based()] {
                current = next;
                if !self.states[current].outputs.is_empty() {
                    matches.push(DictionaryMatch {
                        end_position: position + 1,
                        pattern_indices: self.states[current].outputs.clone(),
                    });
                }
                position += 1;
            } else if let Some(failure) = self.states[current].failure {
                current = failure;
            } else {
                position += 1;
                if let Some(next_arrow) = walk.get(position) {
                    let source = self
                        .quiver
                        .arrow(*next_arrow)
                        .ok_or(StaticDictionaryError::InvalidPattern)?
                        .source()
                        .index();
                    current = self.vertex_states[source.zero_based()];
                }
            }
        }
        Ok(matches)
    }

    /// Translate QPA's `PrefixSearch(dict, path)`.
    ///
    /// # Errors
    ///
    /// Returns an error for a path from another quiver or the zero path.
    pub fn prefix_search(
        &self,
        path: &Path,
    ) -> Result<Vec<DictionaryMatch>, StaticDictionaryError> {
        self.validate_path(path)?;
        let source = endpoint_index(path.source())?;
        let walk = path.walk().ok_or(StaticDictionaryError::ZeroPattern)?;
        let mut current = self.vertex_states[source.zero_based()];
        let mut matches = Vec::new();
        if !self.states[current].own_patterns.is_empty() {
            matches.push(DictionaryMatch {
                end_position: 0,
                pattern_indices: self.states[current].own_patterns.clone(),
            });
        }
        for (position, arrow) in walk.iter().enumerate() {
            let Some(next) = self.states[current].trie_transitions[arrow.zero_based()] else {
                break;
            };
            current = next;
            if !self.states[current].own_patterns.is_empty() {
                matches.push(DictionaryMatch {
                    end_position: position + 1,
                    pattern_indices: self.states[current].own_patterns.clone(),
                });
            }
        }
        Ok(matches)
    }

    /// Translate QPA's `IsFiniteDifference` cycle test.
    #[must_use]
    pub fn is_finite_difference(&self) -> bool {
        self.topological_order().is_some()
    }

    /// Translate QPA's `DifferenceSize`; `None` represents `infinity`.
    #[must_use]
    pub fn difference_size(&self) -> Option<Integer> {
        let order = self.topological_order()?;
        let mut counts = vec![Integer::from(0); self.states.len()];
        for &vertex in &self.vertex_states {
            counts[vertex] += 1_u8;
        }
        let mut total = Integer::from(0);
        for state in order {
            let contribution = counts[state].clone();
            for target in self.transitions(state) {
                counts[target] += &contribution;
            }
            if !self.states[state].is_tip() {
                total += &contribution;
            }
        }
        Some(total)
    }

    /// Translate QPA's `DifferenceWords`; `None` represents `fail` for an
    /// infinite difference language.
    #[must_use]
    pub fn difference_words(&self) -> Option<Vec<Path>> {
        let order = self.topological_order()?;
        let mut state_paths = vec![Vec::new(); self.states.len()];
        for (position, &state) in self.vertex_states.iter().enumerate() {
            if let Some(vertex) = self.quiver.vertex(VertexIndex::from_zero_based(position)) {
                state_paths[state].push(vertex.path());
            }
        }
        let mut paths = Vec::new();
        for state in order {
            let current_paths = state_paths[state].clone();
            for (arrow_index, target) in self.indexed_transitions(state) {
                let Some(arrow) = self.quiver.arrow(arrow_index) else {
                    continue;
                };
                for path in &current_paths {
                    if let Ok(extension) = path.multiply(&arrow.path()) {
                        state_paths[target].push(extension);
                    }
                }
            }
            if !self.states[state].is_tip() {
                paths.extend(state_paths[state].iter().cloned());
            }
        }
        paths.sort_by(|left, right| {
            self.quiver
                .compare_paths(left, right)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        Some(paths)
    }

    fn validate_path(&self, path: &Path) -> Result<(), StaticDictionaryError> {
        if path.quiver().is_identical(&self.quiver) {
            Ok(())
        } else {
            Err(StaticDictionaryError::DifferentQuivers)
        }
    }

    fn transitions(&self, state: usize) -> impl Iterator<Item = usize> + '_ {
        self.indexed_transitions(state).map(|(_, target)| target)
    }

    fn indexed_transitions(&self, state: usize) -> impl Iterator<Item = (ArrowIndex, usize)> + '_ {
        self.states[state]
            .transitions
            .iter()
            .enumerate()
            .filter_map(|(arrow, target)| {
                target.map(|target| (ArrowIndex::from_zero_based(arrow), target))
            })
    }

    fn topological_order(&self) -> Option<Vec<usize>> {
        let mut indegree = vec![0_usize; self.states.len()];
        for state in 0..self.states.len() {
            for target in self.transitions(state) {
                indegree[target] += 1;
            }
        }
        let mut queue = indegree
            .iter()
            .enumerate()
            .filter_map(|(state, &degree)| (degree == 0).then_some(state))
            .collect::<VecDeque<_>>();
        let mut order = Vec::with_capacity(self.states.len());
        while let Some(state) = queue.pop_front() {
            order.push(state);
            for target in self.transitions(state) {
                indegree[target] -= 1;
                if indegree[target] == 0 {
                    queue.push_back(target);
                }
            }
        }
        (order.len() == self.states.len()).then_some(order)
    }
}

impl std::fmt::Display for QuiverStaticDictionary {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.view_string())
    }
}

fn complete_quiver_transitions(
    quiver: &Quiver,
    vertex_states: &[usize],
    states: &mut [AutomatonState],
) -> Result<(), StaticDictionaryError> {
    for state in states.iter_mut() {
        state.transitions.clone_from(&state.trie_transitions);
    }
    let mut queue = VecDeque::new();
    for &vertex_state in vertex_states {
        let vertex = quiver
            .vertex(states[vertex_state].target)
            .ok_or(StaticDictionaryError::InvalidPattern)?;
        for arrow in quiver.outgoing_arrows(&vertex)? {
            let position = arrow.index().zero_based();
            states[vertex_state].transitions[position]
                .get_or_insert(vertex_states[arrow.target().index().zero_based()]);
        }
        queue.extend(
            states[vertex_state]
                .trie_transitions
                .iter()
                .filter_map(|child| *child)
                .filter(|child| states[*child].trie_transitions.iter().any(Option::is_some)),
        );
    }
    while let Some(current) = queue.pop_front() {
        let children = states[current]
            .trie_transitions
            .iter()
            .filter_map(|child| *child)
            .filter(|child| states[*child].own_patterns.is_empty())
            .collect::<Vec<_>>();
        queue.extend(children);
        let vertex = quiver
            .vertex(states[current].target)
            .ok_or(StaticDictionaryError::InvalidPattern)?;
        for arrow in quiver.outgoing_arrows(&vertex)? {
            let position = arrow.index().zero_based();
            if states[current].transitions[position].is_none() {
                let fallback = states[current]
                    .failure
                    .and_then(|failure| states[failure].transitions[position])
                    .ok_or(StaticDictionaryError::InvalidPattern)?;
                states[current].transitions[position] = Some(fallback);
            }
        }
    }
    Ok(())
}

fn failure_transition(
    states: &[AutomatonState],
    vertex_states: &[usize],
    quiver: &Quiver,
    mut failure: Option<usize>,
    arrow_index: ArrowIndex,
) -> Result<usize, StaticDictionaryError> {
    while let Some(state) = failure {
        if let Some(target) = states[state].trie_transitions[arrow_index.zero_based()] {
            return Ok(target);
        }
        failure = states[state].failure;
    }
    let target = quiver
        .arrow(arrow_index)
        .ok_or(StaticDictionaryError::InvalidPattern)?
        .target()
        .index();
    Ok(vertex_states[target.zero_based()])
}

fn endpoint_index(endpoint: PathEndpoint) -> Result<VertexIndex, StaticDictionaryError> {
    match endpoint {
        PathEndpoint::Vertex(vertex) => Ok(vertex.index()),
        PathEndpoint::Zero(_) => Err(StaticDictionaryError::ZeroPattern),
    }
}

/// Failure to construct or query a QPA static dictionary.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum StaticDictionaryError {
    /// QPA's generic dense-list search reads position one before checking length.
    #[error("an empty dense list has no first search symbol")]
    EmptyText,
    /// Patterns and searched paths must belong to the dictionary's quiver.
    #[error("a dictionary path belongs to a different quiver")]
    DifferentQuivers,
    /// QPA's zero path has no source or walk for dictionary construction.
    #[error("the zero path cannot be a static-dictionary pattern or query")]
    ZeroPattern,
    /// A pattern contains an invalid quiver transition.
    #[error("a static-dictionary pattern is not a valid quiver path")]
    InvalidPattern,
    /// A QPA quiver adjacency query failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
}

#[cfg(test)]
mod tests {
    use serde_json::{Value, json};

    use crate::{ArrowIndex, ArrowSpec, Quiver, VertexIndex};

    use super::{
        DictionaryMatch, QuiverStaticDictionary, STATIC_DICTIONARY_CHECKING_CYCLE_INFO,
        STATIC_DICTIONARY_CREATING_AUTOMATA_INFO, StaticDictionary, StaticDictionaryError,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn arrow(index: usize) -> ArrowIndex {
        ArrowIndex::new(index).unwrap()
    }

    fn two_loops() -> Quiver {
        Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "x"),
                ArrowSpec::named(vertex(1), vertex(1), "y"),
            ],
        )
        .unwrap()
    }

    fn matches_json(matches: &[DictionaryMatch]) -> Value {
        json!(
            matches
                .iter()
                .map(|matched| json!([matched.end_position, matched.pattern_indices]))
                .collect::<Vec<_>>()
        )
    }

    #[test]
    fn pinned_qpa_static_dictionary_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/static_dictionary.schema-v1.json"
        ))
        .unwrap();

        let string_patterns =
            ["he", "she", "his", "hers"].map(|pattern| pattern.chars().collect::<Vec<_>>());
        let generic = StaticDictionary::new(&string_patterns);
        assert_eq!(
            json!({
                "patterns": generic
                    .patterns()
                    .iter()
                    .map(|pattern| pattern.iter().collect::<String>())
                    .collect::<Vec<_>>(),
                "ushers": matches_json(&generic.search(&"ushers".chars().collect::<Vec<_>>()).unwrap()),
                "overlap": matches_json(
                    &StaticDictionary::new(&[vec![1, 2], vec![2], vec![1, 2]])
                        .search(&[1, 2, 1, 2])
                        .unwrap()
                )
            }),
            expected["generic"]
        );

        let quiver = two_loops();
        let x = quiver.arrow(arrow(1)).unwrap().path();
        let y = quiver.arrow(arrow(2)).unwrap().path();
        let xx = x.multiply(&x).unwrap();
        let yx = y.multiply(&x).unwrap();
        let yy = y.multiply(&y).unwrap();
        let finite = QuiverStaticDictionary::new(&quiver, &[xx.clone(), yx, yy]).unwrap();
        let infinite = QuiverStaticDictionary::new(&quiver, std::slice::from_ref(&xx)).unwrap();
        let vertex_prefix = QuiverStaticDictionary::new(
            &quiver,
            &[quiver.vertex(vertex(1)).unwrap().path(), x.clone()],
        )
        .unwrap();
        let overlap =
            QuiverStaticDictionary::new(&quiver, &[x.clone(), x.multiply(&y).unwrap()]).unwrap();
        assert_eq!(
            json!({
                "patterns": finite.patterns().iter().map(ToString::to_string).collect::<Vec<_>>(),
                "search": matches_json(&finite.search(&y.multiply(&xx).unwrap()).unwrap()),
                "prefix": matches_json(&finite.prefix_search(&xx).unwrap()),
                "finite": finite.is_finite_difference(),
                "size": finite.difference_size().unwrap().to_string().parse::<u64>().unwrap(),
                "words": finite.difference_words().unwrap().iter().map(ToString::to_string).collect::<Vec<_>>(),
                "infinite": infinite.is_finite_difference(),
                "infinite_size": "infinity",
                "infinite_words_fail": infinite.difference_words().is_none(),
                "vertex_prefix": matches_json(&vertex_prefix.prefix_search(&x).unwrap()),
                "overlap_finite": overlap.is_finite_difference(),
                "overlap_size": "infinity"
            }),
            expected["quiver"]
        );
    }

    #[test]
    fn pinned_qpa_views_and_failures_match() {
        assert_eq!(
            format!(
                "{STATIC_DICTIONARY_CREATING_AUTOMATA_INFO}{STATIC_DICTIONARY_CHECKING_CYCLE_INFO}"
            ),
            include_str!("../../../tests/fixtures/static_dictionary_info.txt")
        );

        let output = format!(
            "GENERIC_0={}\nGENERIC_1={}\nGENERIC_2={}\nQUIVER_0={}\nQUIVER_1={}\nQUIVER_2={}\n",
            StaticDictionary::<i32>::new(&[]),
            StaticDictionary::new(&[vec![1]]),
            StaticDictionary::new(&[vec![1], vec![2]]),
            QuiverStaticDictionary::new(&two_loops(), &[]).unwrap(),
            {
                let quiver = two_loops();
                QuiverStaticDictionary::new(&quiver, &[quiver.arrow(arrow(1)).unwrap().path()])
                    .unwrap()
            },
            {
                let quiver = two_loops();
                QuiverStaticDictionary::new(
                    &quiver,
                    &[
                        quiver.vertex(vertex(1)).unwrap().path(),
                        quiver.arrow(arrow(1)).unwrap().path(),
                    ],
                )
                .unwrap()
            }
        );
        assert_eq!(
            output,
            include_str!("../../../tests/fixtures/static_dictionary_output.txt")
        );

        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/static_dictionary_errors.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["failures"],
            json!({
                "sparse_patterns_success": false,
                "empty_search_success": false,
                "foreign_pattern_success": false,
                "zero_pattern_success": false,
                "foreign_search_success": false,
                "zero_search_success": false
            })
        );
        assert_eq!(
            StaticDictionary::new(&[vec![1]]).search(&[]),
            Err(StaticDictionaryError::EmptyText)
        );
        let quiver = two_loops();
        let foreign = two_loops();
        let foreign_path = foreign.arrow(arrow(1)).unwrap().path();
        assert_eq!(
            QuiverStaticDictionary::new(&quiver, std::slice::from_ref(&foreign_path)).unwrap_err(),
            StaticDictionaryError::DifferentQuivers
        );
        assert_eq!(
            QuiverStaticDictionary::new(&quiver, &[quiver.zero()]).unwrap_err(),
            StaticDictionaryError::ZeroPattern
        );
        let dictionary =
            QuiverStaticDictionary::new(&quiver, &[quiver.arrow(arrow(1)).unwrap().path()])
                .unwrap();
        assert_eq!(
            dictionary.search(&foreign_path),
            Err(StaticDictionaryError::DifferentQuivers)
        );
        assert_eq!(
            dictionary.search(&quiver.zero()),
            Err(StaticDictionaryError::ZeroPattern)
        );
    }

    #[test]
    fn finite_difference_enumerates_qpa_nontips() {
        let quiver = two_loops();
        let x = quiver.arrow(arrow(1)).unwrap().path();
        let y = quiver.arrow(arrow(2)).unwrap().path();
        let patterns = [
            x.multiply(&x).unwrap(),
            y.multiply(&x).unwrap(),
            y.multiply(&y).unwrap(),
        ];
        let dictionary = QuiverStaticDictionary::new(&quiver, &patterns).unwrap();
        assert!(dictionary.is_finite_difference());
        assert_eq!(dictionary.difference_size().unwrap(), 4);
        assert_eq!(
            dictionary
                .difference_words()
                .unwrap()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>(),
            ["v1", "x", "y", "x*y"]
        );
    }

    #[test]
    fn unrestricted_loop_keeps_difference_infinite() {
        let quiver = two_loops();
        let x = quiver.arrow(arrow(1)).unwrap().path();
        let dictionary = QuiverStaticDictionary::new(&quiver, &[x.multiply(&x).unwrap()]).unwrap();
        assert!(!dictionary.is_finite_difference());
        assert!(dictionary.difference_size().is_none());
        assert!(dictionary.difference_words().is_none());
    }

    #[test]
    fn search_reports_one_based_pattern_and_end_positions() {
        let quiver = two_loops();
        let x = quiver.arrow(arrow(1)).unwrap().path();
        let y = quiver.arrow(arrow(2)).unwrap().path();
        let xx = x.multiply(&x).unwrap();
        let yxx = y.multiply(&xx).unwrap();
        let dictionary = QuiverStaticDictionary::new(&quiver, &[xx]).unwrap();
        assert_eq!(
            dictionary.search(&yxx).unwrap(),
            vec![DictionaryMatch {
                end_position: 3,
                pattern_indices: vec![1]
            }]
        );
    }
}
