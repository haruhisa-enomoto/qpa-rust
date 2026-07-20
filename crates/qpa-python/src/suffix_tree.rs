//! Python wrappers for QPA's dynamic suffix-tree dictionary.

use pyo3::{
    Bound, Py, PyAny, PyResult, Python,
    prelude::{pyclass, pyfunction, pymethods},
    types::PyList,
};
use qpa_core::{
    SuffixTree, SuffixTreeEdgeId, SuffixTreeNodeId, SuffixTreeNodeKind, SuffixTreeOccurrence,
    SuffixTreeToken,
};

use crate::{
    invalid,
    util::{intern_symbol, sequence_symbols, symbol_index},
};

/// QPA's mutable generalized suffix tree.
#[pyclass(name = "SuffixTree", unsendable, skip_from_py_object)]
pub(crate) struct PySuffixTree {
    inner: SuffixTree<usize>,
    patterns: Vec<Py<PyAny>>,
    representatives: Vec<Py<PyAny>>,
}

/// Stable identity of one suffix-tree node.
#[pyclass(name = "SuffixTreeNode", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PySuffixTreeNode {
    inner: SuffixTreeNodeId,
    kind: SuffixTreeNodeKind,
}

/// Stable identity of one suffix-tree edge.
#[pyclass(name = "SuffixTreeEdge", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PySuffixTreeEdge {
    inner: SuffixTreeEdgeId,
}

fn occurrences_to_python(values: &[SuffixTreeOccurrence]) -> Vec<(usize, usize)> {
    values
        .iter()
        .map(|value| (value.pattern_number, value.position))
        .collect()
}

#[pymethods]
impl PySuffixTreeNode {
    #[getter]
    #[allow(clippy::unused_self)]
    fn is_suffix_tree_element(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_suffix_tree_node(&self) -> bool {
        true
    }

    #[getter]
    fn is_internal(&self) -> bool {
        self.kind == SuffixTreeNodeKind::Internal
    }

    #[getter]
    fn is_external(&self) -> bool {
        self.kind == SuffixTreeNodeKind::External
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn is_default_internal_rep(&self) -> bool {
        self.is_internal()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

#[pymethods]
impl PySuffixTreeEdge {
    #[getter]
    #[allow(clippy::unused_self)]
    fn is_suffix_tree_element(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_suffix_tree_edge(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl PySuffixTree {
    fn wrap_node(&self, inner: SuffixTreeNodeId) -> PyResult<PySuffixTreeNode> {
        Ok(PySuffixTreeNode {
            kind: self.inner.node_kind(inner).map_err(invalid)?,
            inner,
        })
    }

    fn encode_pattern(
        &mut self,
        py: Python<'_>,
        pattern: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<usize>> {
        sequence_symbols(py, pattern)?
            .iter()
            .map(|symbol| intern_symbol(py, &mut self.representatives, symbol))
            .collect()
    }

    fn encode_sequence(&self, py: Python<'_>, sequence: &Bound<'_, PyAny>) -> PyResult<Vec<usize>> {
        let unknown = self.representatives.len();
        sequence_symbols(py, sequence)?
            .iter()
            .map(|symbol| Ok(symbol_index(py, &self.representatives, symbol)?.unwrap_or(unknown)))
            .collect()
    }

    fn encode_token(
        &self,
        py: Python<'_>,
        token: &Bound<'_, PyAny>,
    ) -> PyResult<Option<SuffixTreeToken<usize>>> {
        if token
            .cast::<PyList>()
            .is_ok_and(pyo3::types::PyListMethods::is_empty)
        {
            return Ok(Some(SuffixTreeToken::Termination));
        }
        let token = token.clone().unbind();
        Ok(symbol_index(py, &self.representatives, &token)?.map(SuffixTreeToken::Symbol))
    }

    fn tokens_to_python(
        &self,
        py: Python<'_>,
        tokens: Vec<SuffixTreeToken<usize>>,
    ) -> PyResult<Vec<Py<PyAny>>> {
        tokens
            .into_iter()
            .map(|token| match token {
                SuffixTreeToken::Symbol(index) => self
                    .representatives
                    .get(index)
                    .map(|symbol| symbol.clone_ref(py))
                    .ok_or_else(|| invalid("suffix-tree symbol token does not exist")),
                SuffixTreeToken::Termination => Ok(PyList::empty(py).into_any().unbind()),
            })
            .collect()
    }
}

#[pymethods]
impl PySuffixTree {
    #[new]
    fn new() -> Self {
        Self {
            inner: SuffixTree::new(),
            patterns: Vec::new(),
            representatives: Vec::new(),
        }
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_suffix_tree(&self) -> bool {
        true
    }

    #[getter]
    #[allow(clippy::unused_self)]
    fn is_default_rep(&self) -> bool {
        true
    }

    #[getter]
    fn root(&self) -> PyResult<PySuffixTreeNode> {
        self.wrap_node(self.inner.root_node())
    }

    #[getter]
    fn dead_length(&self) -> usize {
        self.inner.dead_length()
    }

    #[getter]
    fn number_of_patterns(&self) -> usize {
        self.inner.number_of_patterns()
    }

    #[getter]
    fn patterns(&self, py: Python<'_>) -> Vec<Py<PyAny>> {
        self.patterns
            .iter()
            .map(|pattern| pattern.clone_ref(py))
            .collect()
    }

    fn pattern(&self, py: Python<'_>, one_based_number: usize) -> PyResult<Py<PyAny>> {
        self.inner.pattern(one_based_number).map_err(invalid)?;
        self.patterns
            .get(one_based_number.saturating_sub(1))
            .map(|pattern| pattern.clone_ref(py))
            .ok_or_else(|| invalid("the suffix-tree pattern number does not exist"))
    }

    fn __getitem__(&self, py: Python<'_>, index: usize) -> PyResult<Py<PyAny>> {
        self.pattern(py, index + 1)
    }

    fn is_pattern_active(&self, one_based_number: usize) -> PyResult<bool> {
        self.inner
            .is_pattern_active(one_based_number)
            .map_err(invalid)
    }

    fn insert_pattern(&mut self, py: Python<'_>, pattern: &Bound<'_, PyAny>) -> PyResult<usize> {
        let encoded = self.encode_pattern(py, pattern)?;
        let number = self.inner.insert_pattern(encoded).map_err(invalid)?;
        self.patterns.push(pattern.clone().unbind());
        Ok(number)
    }

    fn delete_pattern(&mut self, one_based_number: usize) -> PyResult<()> {
        self.inner.delete_pattern(one_based_number).map_err(invalid)
    }

    fn all_patterns_in_sequence(
        &self,
        py: Python<'_>,
        sequence: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<(usize, usize)>> {
        let sequence = self.encode_sequence(py, sequence)?;
        Ok(occurrences_to_python(
            &self
                .inner
                .all_patterns_in_sequence(&sequence)
                .map_err(invalid)?,
        ))
    }

    fn sequence_in_patterns(
        &self,
        py: Python<'_>,
        sequence: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<(usize, usize)>> {
        let sequence = self.encode_sequence(py, sequence)?;
        Ok(occurrences_to_python(
            &self
                .inner
                .sequence_in_patterns(&sequence)
                .map_err(invalid)?,
        ))
    }

    fn parent_node(&self, node: &PySuffixTreeNode) -> PyResult<Option<PySuffixTreeNode>> {
        self.inner
            .parent_node(node.inner)
            .map_err(invalid)?
            .map(|node| self.wrap_node(node))
            .transpose()
    }

    fn label_of_node(&self, py: Python<'_>, node: &PySuffixTreeNode) -> PyResult<Vec<Py<PyAny>>> {
        self.tokens_to_python(py, self.inner.label_of_node(node.inner).map_err(invalid)?)
    }

    fn label_of_edge(&self, py: Python<'_>, edge: &PySuffixTreeEdge) -> PyResult<Vec<Py<PyAny>>> {
        self.tokens_to_python(py, self.inner.label_of_edge(edge.inner).map_err(invalid)?)
    }

    fn node_edges(&self, node: &PySuffixTreeNode) -> PyResult<Vec<PySuffixTreeEdge>> {
        Ok(self
            .inner
            .node_edges(node.inner)
            .map_err(invalid)?
            .iter()
            .map(|&inner| PySuffixTreeEdge { inner })
            .collect())
    }

    fn node_suffix_link(&self, node: &PySuffixTreeNode) -> PyResult<Option<PySuffixTreeNode>> {
        self.inner
            .node_suffix_link(node.inner)
            .map_err(invalid)?
            .map(|node| self.wrap_node(node))
            .transpose()
    }

    fn node_is_marked(&self, node: &PySuffixTreeNode) -> PyResult<bool> {
        self.inner.node_is_marked(node.inner).map_err(invalid)
    }

    fn external_occurrences(&self, node: &PySuffixTreeNode) -> PyResult<Vec<(usize, usize)>> {
        Ok(occurrences_to_python(
            self.inner
                .external_occurrences(node.inner)
                .map_err(invalid)?,
        ))
    }

    fn edge_parent(&self, edge: &PySuffixTreeEdge) -> PyResult<Option<PySuffixTreeNode>> {
        self.inner
            .edge_parent(edge.inner)
            .map_err(invalid)?
            .map(|node| self.wrap_node(node))
            .transpose()
    }

    fn edge_child(&self, edge: &PySuffixTreeEdge) -> PyResult<Option<PySuffixTreeNode>> {
        self.inner
            .edge_child(edge.inner)
            .map_err(invalid)?
            .map(|node| self.wrap_node(node))
            .transpose()
    }

    fn edge_starting_at(
        &self,
        py: Python<'_>,
        node: &PySuffixTreeNode,
        token: &Bound<'_, PyAny>,
    ) -> PyResult<Option<PySuffixTreeEdge>> {
        let Some(token) = self.encode_token(py, token)? else {
            return Ok(None);
        };
        Ok(self
            .inner
            .edge_starting_at(node.inner, &token)
            .map_err(invalid)?
            .map(|inner| PySuffixTreeEdge { inner }))
    }

    fn create_internal_node(
        &mut self,
        pattern_number: usize,
        start: usize,
        end: usize,
    ) -> PyResult<PySuffixTreeNode> {
        let inner = self
            .inner
            .create_internal_node(pattern_number, start, end)
            .map_err(invalid)?;
        self.wrap_node(inner)
    }

    fn create_external_node(
        &mut self,
        pattern_number: usize,
        start: usize,
    ) -> PyResult<PySuffixTreeNode> {
        let inner = self
            .inner
            .create_external_node(pattern_number, start)
            .map_err(invalid)?;
        self.wrap_node(inner)
    }

    fn create_edge(
        &mut self,
        parent: &PySuffixTreeNode,
        child: &PySuffixTreeNode,
        pattern_number: usize,
        start: usize,
        end: Option<usize>,
    ) -> PyResult<PySuffixTreeEdge> {
        Ok(PySuffixTreeEdge {
            inner: self
                .inner
                .create_edge(parent.inner, child.inner, pattern_number, start, end)
                .map_err(invalid)?,
        })
    }

    fn add_internal_node(
        &mut self,
        edge: &PySuffixTreeEdge,
        location: usize,
        node: &PySuffixTreeNode,
    ) -> PyResult<()> {
        self.inner
            .add_internal_node(edge.inner, location, node.inner)
            .map_err(invalid)
    }

    fn add_external_node(
        &mut self,
        parent: &PySuffixTreeNode,
        child: &PySuffixTreeNode,
    ) -> PyResult<PySuffixTreeEdge> {
        Ok(PySuffixTreeEdge {
            inner: self
                .inner
                .add_external_node(parent.inner, child.inner)
                .map_err(invalid)?,
        })
    }
}

/// Name of QPA's `InfoDynamicDictionary` information class.
#[pyfunction(name = "info_dynamic_dictionary_name")]
pub(crate) fn py_info_dynamic_dictionary_name() -> &'static str {
    qpa_core::INFO_DYNAMIC_DICTIONARY
}
