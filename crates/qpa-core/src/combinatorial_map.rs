//! QPA combinatorial maps and their foundational algorithms.
//!
//! This module translates the field-independent first half of
//! `lib/combinatorialmap.gi`.  In particular, permutation products use GAP's
//! right action: under `left * right`, a point is acted on by `left` first and
//! by `right` second.  This convention is observable in QPA's face and dual
//! permutations.

use std::{
    collections::VecDeque,
    fmt,
    sync::{Arc, OnceLock},
};

use numerica::domains::{integer::Integer, rational::Rational};
use thiserror::Error;

/// A finite permutation in GAP's one-based point convention.
///
/// Trailing fixed points are omitted, just as GAP permutations can act on
/// larger symmetric groups without changing their value.  Consequently an
/// empty image list represents the identity permutation.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Permutation {
    images: Box<[usize]>,
}

impl Permutation {
    /// Return the identity permutation.
    #[must_use]
    pub fn identity() -> Self {
        Self {
            images: Box::new([]),
        }
    }

    /// Construct a permutation from the images of `1, ..., n`.
    ///
    /// Trailing fixed points are removed from the stored representation.
    ///
    /// # Errors
    ///
    /// Returns an error unless the entries are a bijection of `1, ..., n`.
    pub fn from_images(mut images: Vec<usize>) -> Result<Self, PermutationError> {
        let size = images.len();
        let mut seen = vec![false; size];
        for &image in &images {
            if image == 0 || image > size {
                return Err(PermutationError::ImageOutsideDomain { image, size });
            }
            if seen[image - 1] {
                return Err(PermutationError::RepeatedImage(image));
            }
            seen[image - 1] = true;
        }
        while images.last().is_some_and(|&image| image == images.len()) {
            images.pop();
        }
        Ok(Self {
            images: images.into_boxed_slice(),
        })
    }

    /// Construct disjoint cycle notation on `1, ..., size`.
    ///
    /// Empty and one-point cycles have no effect.  GAP cycle notation used by
    /// QPA's fixtures is disjoint, so repeated points are rejected explicitly.
    ///
    /// # Errors
    ///
    /// Returns an error if a cycle contains zero, a point greater than
    /// `size`, or a point already used by another nontrivial cycle.
    pub fn from_cycles(size: usize, cycles: &[&[usize]]) -> Result<Self, PermutationError> {
        let mut images = (1..=size).collect::<Vec<_>>();
        let mut used = vec![false; size];
        for cycle in cycles {
            if cycle.len() < 2 {
                continue;
            }
            for &point in *cycle {
                if point == 0 || point > size {
                    return Err(PermutationError::CyclePointOutsideDomain { point, size });
                }
                if used[point - 1] {
                    return Err(PermutationError::RepeatedCyclePoint(point));
                }
                used[point - 1] = true;
            }
            for (&point, &image) in cycle
                .iter()
                .zip(cycle.iter().cycle().skip(1))
                .take(cycle.len())
            {
                images[point - 1] = image;
            }
        }
        Self::from_images(images)
    }

    /// Return the largest moved point, or zero for the identity.
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.images.len()
    }

    /// Return the image of a point.
    ///
    /// GAP permutations fix zero and every point beyond their degree.
    #[must_use]
    pub fn image(&self, point: usize) -> usize {
        point
            .checked_sub(1)
            .and_then(|index| self.images.get(index).copied())
            .unwrap_or(point)
    }

    /// Return images on `1, ..., size`, including implicit fixed points.
    ///
    /// # Errors
    ///
    /// Returns an error if `size` is smaller than the largest moved point.
    pub fn images_on(&self, size: usize) -> Result<Vec<usize>, PermutationError> {
        if self.degree() > size {
            return Err(PermutationError::DegreeExceedsSize {
                degree: self.degree(),
                size,
            });
        }
        Ok((1..=size).map(|point| self.image(point)).collect())
    }

    /// Return `self * right` using GAP's right-action convention.
    #[must_use]
    pub fn gap_product(&self, right: &Self) -> Self {
        let degree = self.degree().max(right.degree());
        let images = (1..=degree)
            .map(|point| right.image(self.image(point)))
            .collect();
        Self::from_valid_images(images)
    }

    /// Return the inverse permutation.
    #[must_use]
    pub fn inverse(&self) -> Self {
        let mut images = vec![0; self.degree()];
        for point in 1..=self.degree() {
            images[self.image(point) - 1] = point;
        }
        Self::from_valid_images(images)
    }

    /// Return the orbits meeting the supplied points in GAP domain order.
    ///
    /// This is QPA's single-generator use of GAP `OrbitsPerms`.  Orbit entries
    /// follow repeated forward action by this permutation.
    #[must_use]
    pub fn orbits_on(&self, points: &[usize]) -> Vec<Vec<usize>> {
        let mut points = points.to_vec();
        points.sort_unstable();
        points.dedup();
        let mut visited = Vec::new();
        let mut result = Vec::new();
        for start in points {
            if visited.contains(&start) {
                continue;
            }
            let mut orbit = Vec::new();
            let mut point = start;
            while !orbit.contains(&point) {
                orbit.push(point);
                visited.push(point);
                point = self.image(point);
            }
            result.push(orbit);
        }
        result
    }

    fn transposition(left: usize, right: usize) -> Self {
        if left == right {
            return Self::identity();
        }
        let degree = left.max(right);
        let mut images = (1..=degree).collect::<Vec<_>>();
        images[left - 1] = right;
        images[right - 1] = left;
        Self::from_valid_images(images)
    }

    fn from_valid_images(mut images: Vec<usize>) -> Self {
        while images.last().is_some_and(|&image| image == images.len()) {
            images.pop();
        }
        Self {
            images: images.into_boxed_slice(),
        }
    }

    fn has_order_two(&self) -> bool {
        !self.images.is_empty()
            && (1..=self.degree()).all(|point| self.image(self.image(point)) == point)
    }
}

/// Failure to construct or inspect a finite permutation.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PermutationError {
    /// An image is not in the represented one-based domain.
    #[error("permutation image {image} is outside 1..={size}")]
    ImageOutsideDomain {
        /// The invalid image.
        image: usize,
        /// The represented domain size.
        size: usize,
    },
    /// Two domain points have the same image.
    #[error("permutation image {0} occurs more than once")]
    RepeatedImage(usize),
    /// A cycle point is not in the requested one-based domain.
    #[error("cycle point {point} is outside 1..={size}")]
    CyclePointOutsideDomain {
        /// The invalid point.
        point: usize,
        /// The requested domain size.
        size: usize,
    },
    /// Disjoint cycle notation reused a point.
    #[error("cycle point {0} occurs in more than one nontrivial cycle")]
    RepeatedCyclePoint(usize),
    /// A requested image list is too short for the permutation.
    #[error("permutation degree {degree} exceeds requested size {size}")]
    DegreeExceedsSize {
        /// The largest moved point.
        degree: usize,
        /// The requested domain size.
        size: usize,
    },
}

struct CombinatorialMapData {
    size: usize,
    ordering: Permutation,
    pairing: Permutation,
    marked_edges: Box<[usize]>,
    faces: OnceLock<Vec<Vec<usize>>>,
    genus: OnceLock<Rational>,
    dual: OnceLock<CombinatorialMap>,
    marked_boundaries: OnceLock<Vec<(Vec<usize>, usize)>>,
    boundary_curves: OnceLock<Vec<(Vec<usize>, usize)>>,
}

/// A QPA combinatorial map on the half-edges `1, ..., size`.
///
/// Cloning preserves GAP object identity.  A fresh constructor call creates a
/// distinct object even when all four constructor arguments are equal.
#[derive(Clone)]
pub struct CombinatorialMap(Arc<CombinatorialMapData>);

impl CombinatorialMap {
    /// Translate QPA's four-argument `CombinatorialMap` constructor.
    ///
    /// The pairing must have order exactly two.  Thus fixed half-edges are
    /// allowed, but the identity permutation is rejected.  QPA does not
    /// validate the marked-edge list, and neither does this constructor.
    ///
    /// # Errors
    ///
    /// Returns an error when `size` is zero, either permutation moves a point
    /// outside `1, ..., size`, or the pairing does not have order two.
    pub fn new(
        size: usize,
        ordering: Permutation,
        pairing: Permutation,
        marked_edges: Vec<usize>,
    ) -> Result<Self, CombinatorialMapError> {
        if size == 0 {
            return Err(CombinatorialMapError::SizeMustBePositive);
        }
        if ordering.degree() > size {
            return Err(CombinatorialMapError::OrderingOutsideUnderlyingSet);
        }
        if pairing.degree() > size {
            return Err(CombinatorialMapError::PairingOutsideUnderlyingSet);
        }
        if !pairing.has_order_two() {
            return Err(CombinatorialMapError::PairingNotInvolution);
        }
        Ok(Self::from_valid_parts(
            size,
            ordering,
            pairing,
            marked_edges,
        ))
    }

    fn from_valid_parts(
        size: usize,
        ordering: Permutation,
        pairing: Permutation,
        marked_edges: Vec<usize>,
    ) -> Self {
        Self(Arc::new(CombinatorialMapData {
            size,
            ordering,
            pairing,
            marked_edges: marked_edges.into_boxed_slice(),
            faces: OnceLock::new(),
            genus: OnceLock::new(),
            dual: OnceLock::new(),
            marked_boundaries: OnceLock::new(),
            boundary_curves: OnceLock::new(),
        }))
    }

    /// Return the number of half-edges.
    #[must_use]
    pub fn size(&self) -> usize {
        self.0.size
    }

    /// Return QPA's ordering permutation.
    #[must_use]
    pub fn ordering(&self) -> &Permutation {
        &self.0.ordering
    }

    /// Return QPA's pairing involution.
    #[must_use]
    pub fn pairing(&self) -> &Permutation {
        &self.0.pairing
    }

    /// Return QPA's unvalidated marked-half-edge list.
    #[must_use]
    pub fn marked_edges(&self) -> &[usize] {
        &self.0.marked_edges
    }

    /// Return the face orbits of `pairing * ordering`.
    #[must_use]
    pub fn faces(&self) -> &[Vec<usize>] {
        self.0.faces.get_or_init(|| {
            let face_permutation = self.pairing().gap_product(self.ordering());
            face_permutation.orbits_on(&(1..=self.size()).collect::<Vec<_>>())
        })
    }

    /// Return the QPA genus as an exact Numerica rational.
    ///
    /// QPA uses `size / 2` even when fixed pairing points make the input
    /// topologically degenerate, so fractional and negative values are
    /// intentionally preserved.
    #[must_use]
    pub fn genus(&self) -> &Rational {
        self.0.genus.get_or_init(|| {
            let vertices = self
                .ordering()
                .orbits_on(&(1..=self.size()).collect::<Vec<_>>())
                .len();
            let mut numerator = Integer::from(4);
            numerator += Integer::from(self.size());
            numerator -= Integer::from(vertices) * Integer::from(2);
            numerator -= Integer::from(self.faces().len()) * Integer::from(2);
            Rational::from((numerator, Integer::from(4)))
        })
    }

    /// Return QPA's cached dual combinatorial map.
    #[must_use]
    pub fn dual(&self) -> &Self {
        self.0.dual.get_or_init(|| {
            Self::from_valid_parts(
                self.size(),
                self.pairing().gap_product(self.ordering()),
                self.pairing().clone(),
                self.marked_edges().to_vec(),
            )
        })
    }

    /// Return whether two handles are the same QPA combinatorial-map object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Translate QPA's `RemoveEdgeOfCombinatorialMap` permutation formulas.
    ///
    /// As in QPA, this checks only that `edge` has two in-range entries; it
    /// does not check that those entries are paired.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid edge list or if QPA's resulting pairing
    /// fails the combinatorial-map constructor.
    pub fn remove_edge(&self, edge: &[usize]) -> Result<Self, CombinatorialMapError> {
        if edge.len() != 2 || edge.iter().any(|&point| point == 0 || point > self.size()) {
            return Err(CombinatorialMapError::InvalidEdge);
        }
        let first = edge[0];
        let second = edge[1];
        let ordering = self.ordering();
        let new_ordering = if second == ordering.image(first) {
            ordering
                .gap_product(&Permutation::transposition(first, second))
                .gap_product(&Permutation::transposition(second, ordering.image(second)))
        } else if first == ordering.image(second) {
            ordering
                .gap_product(&Permutation::transposition(second, first))
                .gap_product(&Permutation::transposition(first, ordering.image(first)))
        } else {
            ordering
                .gap_product(&Permutation::transposition(second, ordering.image(second)))
                .gap_product(&Permutation::transposition(first, ordering.image(first)))
        };
        let new_pairing = Permutation::transposition(first, second).gap_product(self.pairing());
        Self::new(self.size(), new_ordering, new_pairing, Vec::new())
    }

    /// Return each face together with QPA's number of marked boundary points.
    #[must_use]
    pub fn marked_boundaries(&self) -> &[(Vec<usize>, usize)] {
        self.0.marked_boundaries.get_or_init(|| {
            let inverse_ordering = self.ordering().inverse();
            let mut result = self
                .faces()
                .iter()
                .map(|face| {
                    let count = face
                        .iter()
                        .filter(|&&half_edge| {
                            self.marked_edges()
                                .contains(&inverse_ordering.image(half_edge))
                        })
                        .count();
                    (face.clone(), count)
                })
                .collect::<Vec<_>>();
            result.sort();
            result
        })
    }

    /// Translate QPA's combinatorial winding number.
    ///
    /// # Errors
    ///
    /// Returns an error for an empty or odd-length curve.  QPA's loop does not
    /// terminate when a requested pair lies in different ordering orbits;
    /// Rust reports that violated documented precondition instead.
    pub fn winding_number(&self, curve: &[usize]) -> Result<i64, CombinatorialMapError> {
        if curve.is_empty() {
            return Err(CombinatorialMapError::TrivialCurve);
        }
        if !curve.len().is_multiple_of(2) {
            return Err(CombinatorialMapError::OddCurveLength);
        }
        let mut curve = curve.to_vec();
        if curve[0] == self.pairing().image(curve[1]) {
            curve.rotate_left(1);
        }

        let mut winding = 0_i64;
        for pair in curve.chunks_exact(2) {
            let mut point = pair[0];
            let target = pair[1];
            if point == target {
                continue;
            }
            let mut local = 1_i64;
            let start = point;
            loop {
                if self.marked_edges().contains(&point) {
                    local = -1;
                }
                point = self.ordering().image(point);
                if point == target {
                    break;
                }
                if point == start {
                    return Err(CombinatorialMapError::CurvePairInDifferentVertices);
                }
            }
            winding += local;
        }
        Ok(winding)
    }

    /// Return QPA's reduced closed curves around the boundary components.
    ///
    /// # Errors
    ///
    /// Returns an explicit error if QPA's four-entry cancellation would index
    /// outside the list.  Valid surface maps satisfy the cancellation
    /// precondition; degenerate pairings accepted by the constructor need not.
    pub fn boundary_curves(&self) -> Result<&[(Vec<usize>, usize)], CombinatorialMapError> {
        if let Some(curves) = self.0.boundary_curves.get() {
            return Ok(curves);
        }
        let curves = self.compute_boundary_curves()?;
        Ok(self.0.boundary_curves.get_or_init(|| curves))
    }

    /// Return QPA's depth-first cover tree rooted at the vertex containing
    /// `half_edge`.
    ///
    /// Roots larger than the map size are accepted and produce the same empty
    /// result as GAP.  Zero has no applicable QPA `IsPosInt` method.
    ///
    /// # Errors
    ///
    /// Returns an error when `half_edge` is zero.
    pub fn depth_search(&self, half_edge: usize) -> Result<Vec<Vec<usize>>, CombinatorialMapError> {
        self.cover_tree_search(half_edge, SearchKind::Depth)
    }

    /// Return QPA's breadth-first cover tree rooted at the vertex containing
    /// `half_edge`.
    ///
    /// # Errors
    ///
    /// Returns an error when `half_edge` is zero.
    pub fn width_search(&self, half_edge: usize) -> Result<Vec<Vec<usize>>, CombinatorialMapError> {
        self.cover_tree_search(half_edge, SearchKind::Width)
    }

    fn compute_boundary_curves(&self) -> Result<Vec<(Vec<usize>, usize)>, CombinatorialMapError> {
        let mut curves = self
            .marked_boundaries()
            .iter()
            .map(|(face, count)| {
                let curve: Vec<usize> = face
                    .iter()
                    .flat_map(|&half_edge| [half_edge, self.pairing().image(half_edge)])
                    .collect();
                (curve, *count)
            })
            .collect::<Vec<_>>();

        for (stored_curve, _) in &mut curves {
            let mut curve = std::mem::take(stored_curve);
            loop {
                if curve.is_empty() {
                    break;
                }
                let mut changed = false;
                if curve.first() == curve.last() {
                    if curve.len() < 4 {
                        return Err(CombinatorialMapError::BoundaryCancellationUndefined);
                    }
                    curve.drain(0..2);
                    curve.truncate(curve.len() - 2);
                    changed = true;
                } else {
                    for right in 1..curve.len() {
                        if curve[right - 1] == curve[right] {
                            let Some(start) = right.checked_sub(2) else {
                                return Err(CombinatorialMapError::BoundaryCancellationUndefined);
                            };
                            if start + 4 > curve.len() {
                                return Err(CombinatorialMapError::BoundaryCancellationUndefined);
                            }
                            curve.drain(start..start + 4);
                            changed = true;
                            break;
                        }
                    }
                }
                if !changed {
                    break;
                }
            }
            *stored_curve = curve;
        }
        Ok(curves)
    }

    fn cover_tree_search(
        &self,
        root: usize,
        kind: SearchKind,
    ) -> Result<Vec<Vec<usize>>, CombinatorialMapError> {
        if root == 0 {
            return Err(CombinatorialMapError::RootMustBePositive);
        }
        let mut root_vertex = self.ordering().orbits_on(&[root]).remove(0);
        root_vertex.sort_unstable();
        let mut visited = Vec::new();
        let mut result = Vec::new();
        let mut pending = VecDeque::from([(root_vertex, root)]);

        while let Some((current_vertex, current_edge)) = match kind {
            SearchKind::Depth => pending.pop_back(),
            SearchKind::Width => pending.pop_front(),
        } {
            if visited.contains(&current_vertex) {
                continue;
            }
            visited.push(current_vertex.clone());
            result.push(current_edge);
            for half_edge in current_vertex {
                let paired = self.pairing().image(half_edge);
                let mut vertex = self.ordering().orbits_on(&[paired]).remove(0);
                vertex.sort_unstable();
                if !visited.contains(&vertex) {
                    pending.push_back((vertex, paired));
                }
            }
        }
        let position = result
            .iter()
            .position(|&half_edge| half_edge == root)
            .expect("the search result begins with its root");
        result.remove(position);
        Ok(self.pairing().orbits_on(&result))
    }
}

impl fmt::Debug for CombinatorialMap {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("CombinatorialMap")
            .field("size", &self.size())
            .field("ordering", self.ordering())
            .field("pairing", self.pairing())
            .field("marked_edges", &self.marked_edges())
            .finish_non_exhaustive()
    }
}

#[derive(Clone, Copy)]
enum SearchKind {
    Depth,
    Width,
}

/// Failure to construct or operate on a QPA combinatorial map.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum CombinatorialMapError {
    /// QPA dispatch requires a positive underlying-set size.
    #[error("the combinatorial-map size must be positive")]
    SizeMustBePositive,
    /// The ordering is not in `SymmetricGroup(size)`.
    #[error("the ordering is not a permutation of proper form")]
    OrderingOutsideUnderlyingSet,
    /// The pairing is not in `SymmetricGroup(size)`.
    #[error("the pairing is not a permutation of proper form")]
    PairingOutsideUnderlyingSet,
    /// The pairing permutation does not have order exactly two.
    #[error("the pairing is not an involution")]
    PairingNotInvolution,
    /// The edge argument is not a two-entry in-range list.
    #[error("the entered list does not correspond to an edge")]
    InvalidEdge,
    /// QPA rejects an empty winding-number curve.
    #[error("the chosen curve is trivial")]
    TrivialCurve,
    /// QPA's pair indexing requires an even curve length.
    #[error("the chosen curve must have even length")]
    OddCurveLength,
    /// A documented curve pair lies in two different ordering orbits.
    #[error("adjacent curve half-edges lie at different vertices")]
    CurvePairInDifferentVertices,
    /// QPA's boundary cancellation would remove outside its list.
    #[error("QPA boundary cancellation is undefined for this degenerate map")]
    BoundaryCancellationUndefined,
    /// Search dispatch requires a positive half-edge integer.
    #[error("the search root must be positive")]
    RootMustBePositive,
}

#[cfg(test)]
mod tests {
    use numerica::domains::rational::Rational;
    use serde_json::{Value, json};

    use super::{CombinatorialMap, CombinatorialMapError, Permutation, PermutationError};

    fn permutation(size: usize, cycles: &[&[usize]]) -> Permutation {
        Permutation::from_cycles(size, cycles).unwrap()
    }

    fn rational_json(value: &Rational) -> Value {
        json!([
            value.numerator().to_string().parse::<i64>().unwrap(),
            value.denominator().to_string().parse::<i64>().unwrap()
        ])
    }

    fn map_json(map: &CombinatorialMap) -> Value {
        json!({
            "size": map.size(),
            "ordering": map.ordering().images_on(map.size()).unwrap(),
            "pairing": map.pairing().images_on(map.size()).unwrap(),
            "marked": map.marked_edges(),
            "faces": map.faces(),
            "genus": rational_json(map.genus()),
            "marked_boundaries": map.marked_boundaries(),
        })
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_combinatorial_map_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/combinatorialmap_foundation.schema-v1.json"
        ))
        .unwrap();

        let torus = CombinatorialMap::new(
            4,
            permutation(4, &[&[1, 2, 3, 4]]),
            permutation(4, &[&[1, 3], &[2, 4]]),
            vec![1, 3],
        )
        .unwrap();
        let removed = torus.remove_edge(&[1, 3]).unwrap();
        let search_map = CombinatorialMap::new(
            8,
            permutation(8, &[&[1, 2, 3], &[4, 5], &[6, 7]]),
            permutation(8, &[&[1, 4], &[2, 6], &[3, 8], &[5, 7]]),
            vec![1, 5, 8],
        )
        .unwrap();
        let fixed_pairing = CombinatorialMap::new(
            3,
            Permutation::identity(),
            permutation(3, &[&[1, 2]]),
            vec![0, 4],
        )
        .unwrap();

        assert_eq!(map_json(&torus), expected["torus"]);
        assert_eq!(map_json(torus.dual()), expected["dual"]);
        assert_eq!(map_json(&removed), expected["removed"]);
        assert_eq!(map_json(&search_map), expected["search_map"]);
        assert_eq!(map_json(&fixed_pairing), expected["fixed_pairing"]);
        assert_eq!(
            json!(search_map.depth_search(1).unwrap()),
            expected["searches"]["depth_1"]
        );
        assert_eq!(
            json!(search_map.width_search(1).unwrap()),
            expected["searches"]["width_1"]
        );
        assert_eq!(
            json!(search_map.depth_search(9).unwrap()),
            expected["searches"]["depth_9"]
        );
        assert_eq!(
            json!(torus.boundary_curves().unwrap()),
            expected["boundary_curves"]["torus"]
        );
        assert_eq!(
            json!(search_map.boundary_curves().unwrap()),
            expected["boundary_curves"]["search_map"]
        );
        assert_eq!(
            json!({
                "equal": torus.winding_number(&[1, 1]).unwrap(),
                "marked": torus.winding_number(&[1, 2]).unwrap(),
                "rotated": torus.winding_number(&[2, 4]).unwrap(),
            }),
            expected["winding"]
        );
    }

    #[test]
    fn qpa_constructor_boundaries_and_attribute_identity_are_explicit() {
        assert!(matches!(
            Permutation::from_images(vec![2, 2]),
            Err(PermutationError::RepeatedImage(2))
        ));
        assert!(matches!(
            CombinatorialMap::new(
                3,
                permutation(4, &[&[1, 4]]),
                permutation(3, &[&[1, 2]]),
                vec![]
            ),
            Err(CombinatorialMapError::OrderingOutsideUnderlyingSet)
        ));
        assert!(matches!(
            CombinatorialMap::new(3, Permutation::identity(), Permutation::identity(), vec![]),
            Err(CombinatorialMapError::PairingNotInvolution)
        ));

        let map = CombinatorialMap::new(
            4,
            permutation(4, &[&[1, 2, 3, 4]]),
            permutation(4, &[&[1, 3], &[2, 4]]),
            vec![1],
        )
        .unwrap();
        assert!(map.is_identical(&map.clone()));
        assert!(map.dual().is_identical(map.dual()));
        let fresh = CombinatorialMap::new(
            4,
            map.ordering().clone(),
            map.pairing().clone(),
            map.marked_edges().to_vec(),
        )
        .unwrap();
        assert!(!map.is_identical(&fresh));
        assert!(matches!(
            map.remove_edge(&[1]),
            Err(CombinatorialMapError::InvalidEdge)
        ));
        assert!(matches!(
            map.winding_number(&[]),
            Err(CombinatorialMapError::TrivialCurve)
        ));
        assert!(matches!(
            map.depth_search(0),
            Err(CombinatorialMapError::RootMustBePositive)
        ));

        let degenerate = CombinatorialMap::new(
            3,
            Permutation::identity(),
            permutation(3, &[&[1, 2]]),
            vec![0, 4],
        )
        .unwrap();
        assert!(matches!(
            degenerate.winding_number(&[1, 2]),
            Err(CombinatorialMapError::CurvePairInDifferentVertices)
        ));
        assert!(matches!(
            degenerate.boundary_curves(),
            Err(CombinatorialMapError::BoundaryCancellationUndefined)
        ));
    }
}
