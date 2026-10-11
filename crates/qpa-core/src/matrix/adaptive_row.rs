use super::sparse::{Row, add_scaled, scale};
use super::{DenseMatrix, Field, add_scaled_row, scale_row};

// Dense storage is per row, never the Cartesian product of a large Hom system.
// The caller accounts for all simultaneously retained dense coefficient slots.
pub(super) enum AdaptiveRow<F: Field> {
    Sparse(Row<F>),
    Dense(Vec<F::Element>),
}

impl<F: Field> AdaptiveRow<F> {
    pub(super) fn promote(&mut self, width: usize, budget: &mut usize, field: &F) -> bool {
        let Self::Sparse(row) = self else {
            return false;
        };
        if width == 0 || row.len() < width.div_ceil(4) || width > *budget {
            return false;
        }
        let mut dense = vec![field.zero(); width];
        for (column, value) in std::mem::take(row) {
            dense[column] = value;
        }
        *self = Self::Dense(dense);
        *budget -= width;
        true
    }

    pub(super) fn sparse_len(&self) -> usize {
        match self {
            Self::Sparse(row) => row.len(),
            Self::Dense(_) => 0,
        }
    }

    pub(super) fn dense_len(&self) -> usize {
        match self {
            Self::Sparse(_) => 0,
            Self::Dense(row) => row.len(),
        }
    }

    pub(super) fn next(&self, start: usize, field: &F) -> Option<(usize, F::Element)> {
        match self {
            Self::Sparse(row) => row.range(start..).next().map(|(&i, v)| (i, v.clone())),
            Self::Dense(row) => row[start..]
                .iter()
                .position(|v| !field.is_zero(v))
                .map(|offset| (start + offset, row[start + offset].clone())),
        }
    }

    pub(super) fn scale(&mut self, factor: &F::Element, field: &F) {
        match self {
            Self::Sparse(row) => scale(row, factor, field),
            Self::Dense(row) => scale_row(field, row, factor),
        }
    }

    pub(super) fn add_scaled(
        &mut self,
        other: &Self,
        factor: &F::Element,
        start: usize,
        field: &F,
    ) {
        match (self, other) {
            (Self::Sparse(row), Self::Sparse(other)) => add_scaled(row, other, factor, field),
            (Self::Dense(row), Self::Dense(other)) => {
                add_scaled_row(field, &mut row[start..], &other[start..], factor);
            }
            (Self::Dense(row), Self::Sparse(other)) => {
                for (&column, value) in other {
                    if field.is_one(factor) {
                        field.add_assign(&mut row[column], value);
                    } else {
                        field.add_mul_assign(&mut row[column], value, factor);
                    }
                }
            }
            (Self::Sparse(row), Self::Dense(other)) => {
                for (column, value) in other.iter().enumerate().skip(start) {
                    if field.is_zero(value) {
                        continue;
                    }
                    let value = if field.is_one(factor) {
                        value.clone()
                    } else {
                        field.mul(value, factor)
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
        }
    }
}

// Keep the dense rank path from exporting a second, unused copy of its basis.
impl<F: Field> DenseMatrix<F> {
    pub(super) fn rank_without_export(&self) -> usize {
        self.semi_echelon::<false, false>()
            .heads
            .iter()
            .filter(|head| head.is_some())
            .count()
    }
}
