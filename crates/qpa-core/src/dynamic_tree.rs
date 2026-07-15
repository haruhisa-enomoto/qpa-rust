//! QPA's Sleator--Tarjan-style dynamic trees.
//!
//! This is a direct translation of `lib/groebnerdyntree.gi`. QPA stores the
//! represented-tree parent and the auxiliary splay-tree children in mutable
//! object fields. Rust stores the same four fields in a shared arena so those
//! cyclic relationships neither require `unsafe` nor leak reference cycles.

use std::{cell::RefCell, fmt, rc::Rc};

use thiserror::Error;

#[derive(Debug)]
struct DynamicTreeNode<Item> {
    item: Option<Item>,
    parent: Option<usize>,
    left: Option<usize>,
    right: Option<usize>,
}

#[derive(Debug)]
struct DynamicForestData<Item> {
    nodes: Vec<DynamicTreeNode<Item>>,
}

/// An arena containing QPA dynamic-tree nodes that may be linked together.
///
/// QPA allocates each node as an independently garbage-collected object. A
/// Rust arena gives the same nodes stable identity while safely owning their
/// cyclic parent/child graph. All nodes passed to one operation must therefore
/// have been created by the same forest.
pub struct DynamicForest<Item> {
    data: Rc<RefCell<DynamicForestData<Item>>>,
}

impl<Item> DynamicForest<Item> {
    /// Create an empty dynamic-tree arena.
    #[must_use]
    pub fn new() -> Self {
        Self {
            data: Rc::new(RefCell::new(DynamicForestData { nodes: Vec::new() })),
        }
    }

    /// Translate QPA's `CreateDynamicTree()` or `CreateDynamicTree(item)`.
    #[must_use]
    pub fn create(&self, item: Option<Item>) -> DynamicTree<Item> {
        let mut data = self.data.borrow_mut();
        let index = data.nodes.len();
        data.nodes.push(DynamicTreeNode {
            item,
            parent: None,
            left: None,
            right: None,
        });
        DynamicTree {
            forest: Rc::clone(&self.data),
            index,
        }
    }

    /// Number of nodes allocated in this arena.
    #[must_use]
    pub fn len(&self) -> usize {
        self.data.borrow().nodes.len()
    }

    /// Whether this arena contains no nodes.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl<Item> Clone for DynamicForest<Item> {
    fn clone(&self) -> Self {
        Self {
            data: Rc::clone(&self.data),
        }
    }
}

impl<Item> Default for DynamicForest<Item> {
    fn default() -> Self {
        Self::new()
    }
}

impl<Item> fmt::Debug for DynamicForest<Item> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("DynamicForest")
            .field("nodes", &self.len())
            .finish()
    }
}

/// An identity-bearing node in a [`DynamicForest`].
pub struct DynamicTree<Item> {
    forest: Rc<RefCell<DynamicForestData<Item>>>,
    index: usize,
}

impl<Item> DynamicTree<Item> {
    /// Whether the optional `item` component was supplied at construction.
    #[must_use]
    pub fn has_item(&self) -> bool {
        self.forest.borrow().nodes[self.index].item.is_some()
    }

    /// Clone the optional item stored in this node.
    #[must_use]
    pub fn item(&self) -> Option<Item>
    where
        Item: Clone,
    {
        self.forest.borrow().nodes[self.index].item.clone()
    }

    /// Replace the mutable `item` component and return its old value.
    pub fn replace_item(&self, item: Option<Item>) -> Option<Item> {
        std::mem::replace(&mut self.forest.borrow_mut().nodes[self.index].item, item)
    }

    /// Current represented/auxiliary parent, or `None` for QPA's `0`.
    #[must_use]
    pub fn parent(&self) -> Option<Self> {
        self.related(self.parent_index())
    }

    /// Current left auxiliary child, or `None` for QPA's `0`.
    #[must_use]
    pub fn left(&self) -> Option<Self> {
        self.related(self.left_index())
    }

    /// Current right auxiliary child, or `None` for QPA's `0`.
    #[must_use]
    pub fn right(&self) -> Option<Self> {
        self.related(self.right_index())
    }

    /// Translate QPA's `RightRotateDynamicTree`.
    ///
    /// # Errors
    ///
    /// Returns [`DynamicTreeError::MissingLeftChild`] when the QPA method
    /// would attempt to access a component of `0`.
    pub fn right_rotate(&self) -> Result<(), DynamicTreeError> {
        Self::right_rotate_index(&self.forest, self.index)
    }

    /// Translate QPA's `LeftRotateDynamicTree`.
    ///
    /// # Errors
    ///
    /// Returns [`DynamicTreeError::MissingRightChild`] when the QPA method
    /// would attempt to access a component of `0`.
    pub fn left_rotate(&self) -> Result<(), DynamicTreeError> {
        Self::left_rotate_index(&self.forest, self.index)
    }

    /// Translate QPA's `SpliceDynamicTree`.
    ///
    /// # Errors
    ///
    /// Returns [`DynamicTreeError::NotMiddleChild`] for QPA's level-one
    /// assertion that the node is neither auxiliary child of its parent.
    pub fn splice(&self) -> Result<(), DynamicTreeError> {
        Self::splice_index(&self.forest, self.index)
    }

    /// Translate QPA's two-stage `SplayDynamicTree` method.
    ///
    /// # Errors
    ///
    /// Propagates an invalid rotation or splice precondition if the stored
    /// dynamic-tree links are inconsistent.
    pub fn splay(&self) -> Result<(), DynamicTreeError> {
        let mut current = self.index;
        while Self::parent_index_in(&self.forest, current).is_some() {
            Self::splay_auxiliary(&self.forest, current)?;
            let Some(parent) = Self::parent_index_in(&self.forest, current) else {
                // QPA assigns its `0` sentinel to `w` here and terminates the
                // next loop test.
                break;
            };
            current = parent;
        }

        current = self.index;
        while Self::parent_index_in(&self.forest, current).is_some() {
            Self::splice_index(&self.forest, current)?;
            let Some(parent) = Self::parent_index_in(&self.forest, current) else {
                break;
            };
            current = parent;
        }
        Self::splay_auxiliary(&self.forest, self.index)
    }

    /// Translate QPA's `RootOfDynamicTree`.
    ///
    /// # Errors
    ///
    /// Propagates an invalid splay precondition if the stored links are
    /// inconsistent.
    pub fn root(&self) -> Result<Self, DynamicTreeError> {
        self.splay()?;
        let mut root = self.index;
        while let Some(right) = Self::right_index_in(&self.forest, root) {
            root = right;
        }
        let root = Self {
            forest: Rc::clone(&self.forest),
            index: root,
        };
        root.splay()?;
        Ok(root)
    }

    /// Translate QPA's `LinkDynamicTrees(v, w)`, with `self` as `v`.
    ///
    /// QPA requires callers not to introduce represented-tree cycles. Like the
    /// source method, this operation does not perform a cycle test.
    ///
    /// # Errors
    ///
    /// Returns [`DynamicTreeError::DifferentForests`] when the nodes do not
    /// share the arena required by Rust's cycle-free ownership model.
    pub fn link(&self, parent: &Self) -> Result<(), DynamicTreeError> {
        self.ensure_same_forest(parent)?;
        self.splay()?;
        parent.splay()?;
        self.forest.borrow_mut().nodes[self.index].parent = Some(parent.index);
        Ok(())
    }

    /// Translate QPA's `CutDynamicTree`.
    ///
    /// # Errors
    ///
    /// Propagates an invalid splay precondition if the stored links are
    /// inconsistent.
    pub fn cut(&self) -> Result<(), DynamicTreeError> {
        self.splay()?;
        let right = self.right_index();
        if let Some(right) = right {
            let mut data = self.forest.borrow_mut();
            data.nodes[right].parent = None;
            data.nodes[self.index].right = None;
        }
        Ok(())
    }

    fn ensure_same_forest(&self, other: &Self) -> Result<(), DynamicTreeError> {
        if Rc::ptr_eq(&self.forest, &other.forest) {
            Ok(())
        } else {
            Err(DynamicTreeError::DifferentForests)
        }
    }

    fn related(&self, index: Option<usize>) -> Option<Self> {
        index.map(|index| Self {
            forest: Rc::clone(&self.forest),
            index,
        })
    }

    fn parent_index(&self) -> Option<usize> {
        Self::parent_index_in(&self.forest, self.index)
    }

    fn left_index(&self) -> Option<usize> {
        Self::left_index_in(&self.forest, self.index)
    }

    fn right_index(&self) -> Option<usize> {
        Self::right_index_in(&self.forest, self.index)
    }

    fn parent_index_in(
        forest: &Rc<RefCell<DynamicForestData<Item>>>,
        index: usize,
    ) -> Option<usize> {
        forest.borrow().nodes[index].parent
    }

    fn left_index_in(forest: &Rc<RefCell<DynamicForestData<Item>>>, index: usize) -> Option<usize> {
        forest.borrow().nodes[index].left
    }

    fn right_index_in(
        forest: &Rc<RefCell<DynamicForestData<Item>>>,
        index: usize,
    ) -> Option<usize> {
        forest.borrow().nodes[index].right
    }

    fn right_rotate_index(
        forest: &Rc<RefCell<DynamicForestData<Item>>>,
        y: usize,
    ) -> Result<(), DynamicTreeError> {
        let mut data = forest.borrow_mut();
        let x = data.nodes[y]
            .left
            .ok_or(DynamicTreeError::MissingLeftChild)?;
        let z = data.nodes[y].parent;
        if let Some(z) = z {
            if data.nodes[z].left == Some(y) {
                data.nodes[z].left = Some(x);
            } else if data.nodes[z].right == Some(y) {
                data.nodes[z].right = Some(x);
            }
        }
        let x_right = data.nodes[x].right;
        data.nodes[y].left = x_right;
        data.nodes[x].right = Some(y);
        data.nodes[x].parent = z;
        data.nodes[y].parent = Some(x);
        if let Some(child) = x_right {
            data.nodes[child].parent = Some(y);
        }
        Ok(())
    }

    fn left_rotate_index(
        forest: &Rc<RefCell<DynamicForestData<Item>>>,
        y: usize,
    ) -> Result<(), DynamicTreeError> {
        let mut data = forest.borrow_mut();
        let x = data.nodes[y]
            .right
            .ok_or(DynamicTreeError::MissingRightChild)?;
        let z = data.nodes[y].parent;
        if let Some(z) = z {
            if data.nodes[z].left == Some(y) {
                data.nodes[z].left = Some(x);
            } else if data.nodes[z].right == Some(y) {
                data.nodes[z].right = Some(x);
            }
        }
        let x_left = data.nodes[x].left;
        data.nodes[y].right = x_left;
        data.nodes[x].left = Some(y);
        data.nodes[x].parent = z;
        data.nodes[y].parent = Some(x);
        if let Some(child) = x_left {
            data.nodes[child].parent = Some(y);
        }
        Ok(())
    }

    fn splice_index(
        forest: &Rc<RefCell<DynamicForestData<Item>>>,
        node: usize,
    ) -> Result<(), DynamicTreeError> {
        let mut data = forest.borrow_mut();
        if let Some(parent) = data.nodes[node].parent {
            if data.nodes[parent].left == Some(node) || data.nodes[parent].right == Some(node) {
                return Err(DynamicTreeError::NotMiddleChild);
            }
            data.nodes[parent].left = Some(node);
        }
        Ok(())
    }

    fn splay_auxiliary(
        forest: &Rc<RefCell<DynamicForestData<Item>>>,
        node: usize,
    ) -> Result<(), DynamicTreeError> {
        loop {
            let Some(parent) = Self::parent_index_in(forest, node) else {
                return Ok(());
            };
            let is_left = Self::left_index_in(forest, parent) == Some(node);
            let is_right = Self::right_index_in(forest, parent) == Some(node);
            if !is_left && !is_right {
                return Ok(());
            }

            let grandparent = Self::parent_index_in(forest, parent);
            if is_left {
                if grandparent.is_some_and(|grandparent| {
                    Self::left_index_in(forest, grandparent) == Some(parent)
                }) {
                    Self::right_rotate_index(
                        forest,
                        grandparent.expect("the source comparison found a grandparent"),
                    )?;
                }
                let current_parent =
                    Self::parent_index_in(forest, node).expect("an auxiliary child has a parent");
                Self::right_rotate_index(forest, current_parent)?;
            } else if is_right {
                if grandparent.is_some_and(|grandparent| {
                    Self::right_index_in(forest, grandparent) == Some(parent)
                }) {
                    Self::left_rotate_index(
                        forest,
                        grandparent.expect("the source comparison found a grandparent"),
                    )?;
                }
                let current_parent =
                    Self::parent_index_in(forest, node).expect("an auxiliary child has a parent");
                Self::left_rotate_index(forest, current_parent)?;
            }
        }
    }
}

impl<Item> Clone for DynamicTree<Item> {
    fn clone(&self) -> Self {
        Self {
            forest: Rc::clone(&self.forest),
            index: self.index,
        }
    }
}

impl<Item> PartialEq for DynamicTree<Item> {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index && Rc::ptr_eq(&self.forest, &other.forest)
    }
}

impl<Item> Eq for DynamicTree<Item> {}

impl<Item> fmt::Debug for DynamicTree<Item> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("DynamicTree")
            .field("index", &self.index)
            .finish_non_exhaustive()
    }
}

/// Typed failures replacing QPA record-component errors and assertions.
#[derive(Clone, Copy, Debug, Eq, Error, PartialEq)]
pub enum DynamicTreeError {
    /// A right rotation requires `y!.left` to be a node.
    #[error("a right rotation requires a left child")]
    MissingLeftChild,
    /// A left rotation requires `y!.right` to be a node.
    #[error("a left rotation requires a right child")]
    MissingRightChild,
    /// `SpliceDynamicTree` requires a represented-tree middle child.
    #[error("the node is not a middle child")]
    NotMiddleChild,
    /// Rust arena handles from different ownership domains cannot be linked.
    #[error("dynamic-tree nodes belong to different forests")]
    DifferentForests,
}

#[cfg(test)]
mod tests {
    use serde_json::{Value, json};

    use super::{DynamicForest, DynamicTree, DynamicTreeError};

    fn set_auxiliary_children<Item>(
        parent: &DynamicTree<Item>,
        left: &DynamicTree<Item>,
        right: &DynamicTree<Item>,
    ) {
        parent.ensure_same_forest(left).unwrap();
        parent.ensure_same_forest(right).unwrap();
        let mut data = parent.forest.borrow_mut();
        data.nodes[parent.index].left = Some(left.index);
        data.nodes[parent.index].right = Some(right.index);
        data.nodes[left.index].parent = Some(parent.index);
        data.nodes[right.index].parent = Some(parent.index);
    }

    fn set_parent<Item>(node: &DynamicTree<Item>, parent: &DynamicTree<Item>) {
        node.ensure_same_forest(parent).unwrap();
        node.forest.borrow_mut().nodes[node.index].parent = Some(parent.index);
    }

    fn label(node: Option<DynamicTree<&'static str>>) -> Value {
        node.and_then(|node| node.item())
            .map_or(Value::Null, Value::from)
    }

    fn snapshot(nodes: &[DynamicTree<&'static str>]) -> Value {
        json!(
            nodes
                .iter()
                .map(|node| [
                    label(Some(node.clone())),
                    label(node.parent()),
                    label(node.left()),
                    label(node.right()),
                ])
                .collect::<Vec<_>>()
        )
    }

    fn root_labels(nodes: &[DynamicTree<&'static str>]) -> Value {
        json!(
            nodes
                .iter()
                .map(|node| node.root().unwrap().item().unwrap())
                .collect::<Vec<_>>()
        )
    }

    #[test]
    fn pinned_qpa_dynamic_tree_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/dynamic_tree.schema-v1.json"
        ))
        .unwrap();

        let forest = DynamicForest::new();
        let left = forest.create(Some("left"));
        let middle = forest.create(Some("middle"));
        let right = forest.create(Some("right"));
        set_auxiliary_children(&middle, &left, &right);
        middle.right_rotate().unwrap();
        assert_eq!(snapshot(&[left, middle, right]), expected["right_rotation"]);

        let left = forest.create(Some("left"));
        let middle = forest.create(Some("middle"));
        let right = forest.create(Some("right"));
        set_auxiliary_children(&middle, &left, &right);
        middle.left_rotate().unwrap();
        assert_eq!(snapshot(&[left, middle, right]), expected["left_rotation"]);

        let splice_parent = forest.create(Some("parent"));
        let splice_child = forest.create(Some("child"));
        splice_child.link(&splice_parent).unwrap();
        splice_child.splice().unwrap();
        assert_eq!(snapshot(&[splice_parent, splice_child]), expected["splice"]);

        let a = forest.create(Some("a"));
        let b = forest.create(Some("b"));
        let c = forest.create(Some("c"));
        let d = forest.create(Some("d"));
        a.link(&b).unwrap();
        b.link(&c).unwrap();
        d.link(&c).unwrap();
        let nodes = [a, b, c, d];
        assert_eq!(root_labels(&nodes), expected["roots_before"]);
        nodes[1].cut().unwrap();
        assert_eq!(root_labels(&nodes), expected["roots_after_cut"]);
        nodes[1].link(&nodes[3]).unwrap();
        assert_eq!(root_labels(&nodes), expected["roots_after_relink"]);
    }

    #[test]
    fn qpa_identity_optional_item_and_failures_are_explicit() {
        let forest = DynamicForest::<String>::new();
        assert!(forest.is_empty());
        let empty = forest.create(None);
        let clone = empty.clone();
        let other = forest.create(None);
        assert_eq!(forest.len(), 2);
        assert!(!empty.has_item());
        assert_eq!(empty.replace_item(Some("item".into())), None);
        assert_eq!(empty.item().as_deref(), Some("item"));
        assert_eq!(empty, clone);
        assert_ne!(empty, other);

        assert_eq!(
            empty.right_rotate(),
            Err(DynamicTreeError::MissingLeftChild)
        );
        assert_eq!(
            empty.left_rotate(),
            Err(DynamicTreeError::MissingRightChild)
        );

        let other_forest = DynamicForest::new();
        let foreign = other_forest.create(Some("foreign".to_owned()));
        assert_eq!(
            empty.link(&foreign),
            Err(DynamicTreeError::DifferentForests)
        );
    }

    #[test]
    fn splice_enforces_qpa_middle_child_assertion() {
        let forest = DynamicForest::new();
        let left = forest.create(Some("left"));
        let parent = forest.create(Some("parent"));
        let right = forest.create(Some("right"));
        set_auxiliary_children(&parent, &left, &right);
        assert_eq!(left.splice(), Err(DynamicTreeError::NotMiddleChild));
        assert_eq!(right.splice(), Err(DynamicTreeError::NotMiddleChild));
    }

    #[test]
    fn rotations_preserve_every_qpa_parent_and_transfer_branch() {
        let forest = DynamicForest::new();
        let z = forest.create(Some("z"));
        let y = forest.create(Some("y"));
        let sibling = forest.create(Some("sibling"));
        let x = forest.create(Some("x"));
        let y_right = forest.create(Some("y-right"));
        let x_left = forest.create(Some("x-left"));
        let transfer = forest.create(Some("transfer"));
        set_auxiliary_children(&z, &y, &sibling);
        set_auxiliary_children(&y, &x, &y_right);
        set_auxiliary_children(&x, &x_left, &transfer);

        y.right_rotate().unwrap();
        assert_eq!(z.left(), Some(x.clone()));
        assert_eq!(x.parent(), Some(z));
        assert_eq!(x.right(), Some(y.clone()));
        assert_eq!(y.parent(), Some(x));
        assert_eq!(y.left(), Some(transfer.clone()));
        assert_eq!(transfer.parent(), Some(y));

        let z = forest.create(Some("z2"));
        let sibling = forest.create(Some("sibling2"));
        let y = forest.create(Some("y2"));
        let y_left = forest.create(Some("y-left2"));
        let x = forest.create(Some("x2"));
        let transfer = forest.create(Some("transfer2"));
        let x_right = forest.create(Some("x-right2"));
        set_auxiliary_children(&z, &sibling, &y);
        set_auxiliary_children(&y, &y_left, &x);
        set_auxiliary_children(&x, &transfer, &x_right);

        y.left_rotate().unwrap();
        assert_eq!(z.right(), Some(x.clone()));
        assert_eq!(x.parent(), Some(z));
        assert_eq!(x.left(), Some(y.clone()));
        assert_eq!(y.parent(), Some(x));
        assert_eq!(y.right(), Some(transfer.clone()));
        assert_eq!(transfer.parent(), Some(y));

        let middle_parent = forest.create(Some("middle-parent"));
        let middle = forest.create(Some("middle"));
        let middle_left = forest.create(Some("middle-left"));
        let middle_right = forest.create(Some("middle-right"));
        set_auxiliary_children(&middle, &middle_left, &middle_right);
        set_parent(&middle, &middle_parent);
        middle.right_rotate().unwrap();
        assert_eq!(middle_parent.left(), None);
        assert_eq!(middle_parent.right(), None);
        assert_eq!(middle_left.parent(), Some(middle_parent));
    }
}
