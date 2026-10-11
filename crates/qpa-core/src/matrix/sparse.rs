//! Sparse storage for GAP's insertion-ordered semi-echelon algorithm.
//!
//! The input-row order, increasing pivot-column sweep, normalization and
//! relation-vector order match `DenseMatrix` exactly. Only absent zeros and
//! unused transformation data are omitted.

use std::collections::BTreeMap;

use numerica::domains::Field;

use super::{DenseMatrix, adaptive_row::AdaptiveRow, binary::PackedSemiEchelon};

#[derive(Default, Debug)]
struct FillProfile {
    promotions: usize,
    peak_sparse_entries: usize,
}

pub(super) type Row<F> = BTreeMap<usize, <F as numerica::domains::Set>::Element>;

/// Incremental exact span membership; the accepted original rows span the
/// same space as all inputs. Final canonicalization is left to the caller.
pub(crate) struct IncrementalRowSpan<F: Field> {
    heads: BTreeMap<usize, Row<F>>,
    field: F,
}

impl<F: Field> IncrementalRowSpan<F> {
    pub(crate) fn new(field: F) -> Self {
        Self {
            heads: BTreeMap::new(),
            field,
        }
    }

    pub(crate) fn insert(&mut self, values: &[F::Element]) -> bool {
        let mut row: Row<F> = values
            .iter()
            .enumerate()
            .filter(|(_, value)| !self.field.is_zero(value))
            .map(|(i, value)| (i, value.clone()))
            .collect();
        while let Some((&pivot, value)) = row.first_key_value() {
            if let Some(head) = self.heads.get(&pivot) {
                let factor = self.field.neg(value);
                add_scaled(&mut row, head, &factor, &self.field);
            } else {
                let inverse = self.field.inv(value);
                scale(&mut row, &inverse, &self.field);
                self.heads.insert(pivot, row);
                return true;
            }
        }
        false
    }
}

pub(crate) struct SparseMatrix<F: Field> {
    rows: Vec<Row<F>>,
    nrows: usize,
    ncols: usize,
    field: F,
}

impl<F: Field> SparseMatrix<F> {
    pub(crate) fn zero(nrows: usize, ncols: usize, field: F) -> Self {
        Self {
            // Unconstrained Hom dimensions must remain constant-space even
            // when there are very many possible vertex-map entries.
            rows: if ncols == 0 {
                Vec::new()
            } else {
                vec![BTreeMap::new(); nrows]
            },
            nrows,
            ncols,
            field,
        }
    }

    pub(crate) fn nrows(&self) -> usize {
        self.nrows
    }
    pub(crate) fn ncols(&self) -> usize {
        self.ncols
    }
    pub(crate) fn field(&self) -> &F {
        &self.field
    }

    pub(crate) fn set(&mut self, row: usize, column: usize, value: F::Element) {
        debug_assert!(column < self.ncols);
        if self.field.is_zero(&value) {
            self.rows[row].remove(&column);
        } else {
            self.rows[row].insert(column, value);
        }
    }

    pub(crate) fn add(&mut self, row: usize, column: usize, value: &F::Element) {
        if self.field.is_zero(value) {
            return;
        }
        let value = self.rows[row]
            .get(&column)
            .map_or_else(|| value.clone(), |old| self.field.add(old, value));
        self.set(row, column, value);
    }

    // Keep the existing packed/dense kernels for small systems. Large Hom
    // systems never allocate the Cartesian product merely to test its rank.
    fn small_dense(&self) -> Option<DenseMatrix<F>> {
        if self.nrows().checked_mul(self.ncols)? > 65_536 {
            return None;
        }
        let mut matrix = DenseMatrix::zero(self.nrows(), self.ncols, self.field.clone());
        for (i, row) in self.rows.iter().enumerate() {
            for (&j, value) in row {
                matrix[(i, j)] = value.clone();
            }
        }
        Some(matrix)
    }

    // Pack sparse GF(2) input directly; do not allocate a dense scalar matrix
    // first. The bound includes transformation/relation bits when requested.
    fn packed<const TRANSFORM: bool>(&self) -> Option<PackedSemiEchelon> {
        if self.nrows < 32 || self.ncols < 64 || self.field.size() != Some(2.into()) {
            return None;
        }
        let words = self.ncols.div_ceil(64).checked_add(if TRANSFORM {
            self.nrows.div_ceil(64)
        } else {
            0
        })?;
        if self.nrows.checked_mul(words)?.checked_mul(8)? > 64 * 1024 * 1024 {
            return None;
        }
        // A packed pivot sweep loses to sparse lookup on very sparse inputs.
        // Keep those sparse (individual rows may still adapt after fill-in).
        let entries = self
            .rows
            .iter()
            .try_fold(0_usize, |n, row| n.checked_add(row.len()))?;
        if entries < self.nrows.checked_mul(self.ncols)?.div_ceil(32) {
            return None;
        }
        Some(PackedSemiEchelon::from_rows::<TRANSFORM>(
            self.nrows,
            self.ncols,
            self.rows.iter().map(|row| {
                let mut bits = vec![0_u64; self.ncols.div_ceil(64)];
                for &column in row.keys() {
                    bits[column / 64] |= 1 << (column % 64);
                }
                bits
            }),
        ))
    }

    pub(crate) fn rank(&self) -> usize {
        if self.ncols == 0 {
            return 0;
        }
        if let Some(packed) = self.packed::<false>() {
            return packed.vectors.len();
        }
        self.small_dense()
            .map_or_else(|| self.eliminate::<false>().0, |matrix| matrix.rank())
    }

    pub(crate) fn nullspace_rows(&self) -> DenseMatrix<F> {
        if let Some(packed) = self.packed::<true>() {
            return packed.relation_matrix(self.nrows, &self.field);
        }
        self.small_dense().map_or_else(
            || self.sparse_nullspace_rows(),
            |matrix| matrix.nullspace_rows(),
        )
    }

    fn sparse_nullspace_rows(&self) -> DenseMatrix<F> {
        let (_, relations, _) = self.eliminate::<true>();
        let mut result = DenseMatrix::zero(relations.len(), self.nrows(), self.field.clone());
        for (row, relation) in relations.into_iter().enumerate() {
            for (column, value) in relation {
                result[(row, column)] = value;
            }
        }
        result
    }

    fn eliminate<const TRANSFORM: bool>(&self) -> (usize, Vec<Row<F>>, FillProfile) {
        let mut profile = FillProfile::default();
        // Bound retained inline dense storage, independent of matrix shape.
        let mut budget = (64 * 1024 * 1024) / std::mem::size_of::<F::Element>().max(1);
        if self.ncols == 0 {
            let relations = if TRANSFORM {
                (0..self.nrows)
                    .map(|i| BTreeMap::from([(i, self.field.one())]))
                    .collect()
            } else {
                Vec::new()
            };
            return (0, relations, profile);
        }
        let mut heads = vec![None; self.ncols];
        let mut vectors: Vec<AdaptiveRow<F>> = Vec::new();
        let mut coefficients: Vec<Row<F>> = Vec::new();
        let mut relations = Vec::new();
        for (index, input) in self.rows.iter().enumerate() {
            let mut row = AdaptiveRow::Sparse(input.clone());
            profile.peak_sparse_entries = profile.peak_sparse_entries.max(row.sparse_len());
            profile.promotions += usize::from(row.promote(self.ncols, &mut budget, &self.field));
            let mut combination = BTreeMap::new();
            if TRANSFORM {
                combination.insert(index, self.field.one());
            }
            let mut start = 0;
            while let Some((column, coefficient)) = row.next(start, &self.field) {
                if let Some(head) = heads[column] {
                    let scale = self.field.neg(&coefficient);
                    row.add_scaled(&vectors[head], &scale, column, &self.field);
                    profile.peak_sparse_entries = profile.peak_sparse_entries.max(row.sparse_len());
                    profile.promotions +=
                        usize::from(row.promote(self.ncols, &mut budget, &self.field));
                    if TRANSFORM {
                        add_scaled(&mut combination, &coefficients[head], &scale, &self.field);
                    }
                }
                start = column + 1;
            }
            if let Some((pivot, value)) = row.next(0, &self.field) {
                let inverse = self.field.inv(&value);
                row.scale(&inverse, &self.field);
                if TRANSFORM {
                    scale(&mut combination, &inverse, &self.field);
                    coefficients.push(combination);
                }
                heads[pivot] = Some(vectors.len());
                vectors.push(row);
            } else {
                budget += row.dense_len();
                if TRANSFORM {
                    relations.push(combination);
                }
            }
        }
        (vectors.len(), relations, profile)
    }
}

pub(super) fn add_scaled<F: Field>(
    row: &mut Row<F>,
    other: &Row<F>,
    scalar: &F::Element,
    field: &F,
) {
    for (&column, value) in other {
        let value = if field.is_one(scalar) {
            value.clone()
        } else {
            field.mul(value, scalar)
        };
        let value = row
            .get(&column)
            .map_or_else(|| value.clone(), |old| field.add(old, &value));
        if field.is_zero(&value) {
            row.remove(&column);
        } else {
            row.insert(column, value);
        }
    }
}

pub(super) fn scale<F: Field>(row: &mut Row<F>, scalar: &F::Element, field: &F) {
    if field.is_one(scalar) {
        return;
    }
    for value in row.values_mut() {
        *value = field.mul(&*value, scalar);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use numerica::domains::{
        finite_field::{Z2, Zp},
        rational::Q,
    };

    fn compare<F: Field>(field: &F) {
        // Includes new early pivots with already-existing later pivots, zero
        // matrices, rectangular shapes, dependent rows and cancellation.
        for nrows in 0..9 {
            for ncols in 0..10 {
                for seed in 0..9 {
                    let mut dense = DenseMatrix::zero(nrows, ncols, field.clone());
                    let mut sparse = SparseMatrix::zero(nrows, ncols, field.clone());
                    for i in 0..nrows {
                        for j in 0..ncols {
                            let value = match (i * 17 + j * 13 + seed * 7 + i * j) % 7 {
                                0..=2 => field.zero(),
                                3 => field.one(),
                                4 => field.neg(&field.one()),
                                _ => field.add(&field.one(), &field.one()),
                            };
                            dense[(i, j)] = value.clone();
                            sparse.set(i, j, value);
                        }
                    }
                    let mut span = IncrementalRowSpan::new(field.clone());
                    let mut accepted = Vec::new();
                    let mut rank = 0;
                    for i in 0..nrows {
                        if span.insert(dense.row(i)) {
                            rank += 1;
                            accepted.extend_from_slice(dense.row(i));
                        }
                    }
                    let accepted =
                        DenseMatrix::from_flat(rank, ncols, accepted, field.clone()).unwrap();
                    assert_eq!(rank, dense.rank());
                    assert_eq!(
                        accepted.canonical_row_space_basis(),
                        dense.canonical_row_space_basis()
                    );
                    assert_eq!(sparse.eliminate::<false>().0, dense.rank());
                    assert_eq!(sparse.sparse_nullspace_rows(), dense.nullspace_rows());
                    assert_eq!(sparse.rank(), dense.rank());
                    assert_eq!(sparse.nullspace_rows(), dense.nullspace_rows());
                }
            }
        }
        let mut sparse = SparseMatrix::zero(1, 1, field.clone());
        sparse.add(0, 0, &field.one());
        sparse.add(0, 0, &field.neg(&field.one()));
        assert_eq!(sparse.rank(), 0);
    }

    #[test]
    fn unconstrained_shape_does_not_allocate_empty_rows() {
        let matrix = SparseMatrix::zero(1_000_000_000, 0, Q);
        assert_eq!(matrix.nrows(), 1_000_000_000);
        assert!(matrix.rows.is_empty());
        assert_eq!(matrix.rank(), 0);
    }

    fn compare_fill<F: Field>(field: &F) {
        let (nrows, ncols) = (96, 768);
        let mut dense = DenseMatrix::zero(nrows, ncols, field.clone());
        let mut state = 1937_u64;
        for row in 0..48 {
            let pivot = 47 - row;
            dense[(row, pivot)] = field.one();
            for column in pivot + 1..ncols {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                if state >> 59 == 0 || (column < 48 && state >> 62 == 0) {
                    dense[(row, column)] = field.one();
                }
            }
        }
        for row in 48..nrows {
            for column in 0..ncols {
                dense[(row, column)] = field.add(
                    &dense[(row - 48, column)],
                    &dense[((row - 47) % 48, column)],
                );
            }
        }
        let mut sparse = SparseMatrix::zero(nrows, ncols, field.clone());
        for row in 0..nrows {
            for column in 0..ncols {
                sparse.set(row, column, dense[(row, column)].clone());
            }
        }
        let (rank, _, profile) = sparse.eliminate::<true>();
        assert_eq!(rank, dense.rank());
        assert_eq!(sparse.sparse_nullspace_rows(), dense.nullspace_rows());
        assert_eq!(sparse.nullspace_rows(), dense.nullspace_rows());
        assert_eq!(sparse.rank(), dense.rank());
        assert!(profile.promotions > 0, "{profile:?}");
        eprintln!(
            "fill profile {}: input_max={}, peak_sparse={}, promotions={}",
            std::any::type_name::<F>(),
            sparse.rows.iter().map(Row::<F>::len).max().unwrap(),
            profile.peak_sparse_entries,
            profile.promotions
        );
        let mut no_budget = 0;
        let mut row = AdaptiveRow::<F>::Sparse((0..ncols).map(|i| (i, field.one())).collect());
        assert!(!row.promote(ncols, &mut no_budget, field));
        if field.size() == Some(2.into()) {
            assert!(sparse.packed::<true>().is_some());
        } else {
            assert!(sparse.packed::<true>().is_none());
        }
    }

    #[test]
    fn binary_packing_respects_density_and_memory_guards() {
        let mut diagonal = SparseMatrix::zero(324, 648, Z2);
        for i in 0..324 {
            diagonal.set(i, i, 1);
        }
        assert!(diagonal.packed::<false>().is_none());
        assert_eq!(diagonal.rank(), 324);
        let huge = SparseMatrix::zero(32, usize::MAX, Z2);
        assert!(huge.packed::<true>().is_none());
    }

    #[test]
    fn fill_in_adaptation_and_direct_binary_packing_preserve_ordered_nullspaces() {
        compare_fill(&Q);
        compare_fill(&Z2);
        compare_fill(&Zp::new(101));
        compare_fill(&crate::GaloisField::conway(2, 2).unwrap());
    }

    #[test]
    fn sparse_elimination_matches_dense_rank_and_ordered_nullspace() {
        compare(&Q);
        compare(&Z2);
        compare(&Zp::new(101));
    }
}
