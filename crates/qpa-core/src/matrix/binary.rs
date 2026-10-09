//! The same column-ordered elimination with GF(2) rows stored in machine words.

use super::{DenseMatrix, Field, SemiEchelonTransformation};

pub(super) struct PackedSemiEchelon {
    heads: Vec<Option<usize>>,
    pub(super) vectors: Vec<Vec<u64>>,
    coefficients: Vec<Vec<u64>>,
    relations: Vec<Vec<u64>>,
}

impl<F: Field> DenseMatrix<F> {
    pub(super) fn use_packed_binary(&self) -> bool {
        // Packing is amortized by larger reductions. Test cardinality, not
        // characteristic: extension fields of characteristic two cannot use XOR
        // for multiplication by an arbitrary nonzero pivot coefficient.
        self.nrows >= 32 && self.ncols >= 64 && self.field.size() == Some(2.into())
    }

    pub(super) fn packed_binary<const TRANSFORM: bool>(&self) -> PackedSemiEchelon {
        let mut result = PackedSemiEchelon {
            heads: vec![None; self.ncols],
            vectors: Vec::new(),
            coefficients: Vec::new(),
            relations: Vec::new(),
        };
        for row_index in 0..self.nrows {
            let mut row = vec![0_u64; self.ncols.div_ceil(64)];
            for (column, entry) in self.row(row_index).iter().enumerate() {
                if !self.field.is_zero(entry) {
                    row[column / 64] |= 1 << (column % 64);
                }
            }
            let mut combination = if TRANSFORM {
                let mut row = vec![0_u64; self.nrows.div_ceil(64)];
                row[row_index / 64] = 1 << (row_index % 64);
                row
            } else {
                Vec::new()
            };
            // Preserve GAP's column order even when pivots were discovered in
            // another order. Earlier bits of each basis vector are zero.
            for (column, &head) in result.heads.iter().enumerate() {
                if let Some(head) = head
                    && row[column / 64] & (1 << (column % 64)) != 0
                {
                    xor_row(
                        &mut row[column / 64..],
                        &result.vectors[head][column / 64..],
                    );
                    if TRANSFORM {
                        let end = row_index / 64 + 1;
                        xor_row(&mut combination[..end], &result.coefficients[head][..end]);
                    }
                }
            }
            let pivot = row.iter().enumerate().find_map(|(word, &bits)| {
                (bits != 0).then(|| word * 64 + bits.trailing_zeros() as usize)
            });
            if let Some(column) = pivot {
                result.heads[column] = Some(result.vectors.len());
                result.vectors.push(row);
                if TRANSFORM {
                    result.coefficients.push(combination);
                }
            } else if TRANSFORM {
                result.relations.push(combination);
            }
        }
        result
    }
}

impl PackedSemiEchelon {
    pub(super) fn unpack_relations<F: Field>(&self, matrix: &DenseMatrix<F>) -> DenseMatrix<F> {
        unpack_rows(&self.relations, matrix.nrows, &matrix.field)
    }

    pub(super) fn unpack<F: Field>(self, matrix: &DenseMatrix<F>) -> SemiEchelonTransformation<F> {
        SemiEchelonTransformation {
            heads: self.heads,
            vectors: unpack_rows(&self.vectors, matrix.ncols, &matrix.field),
            coefficients: unpack_rows(&self.coefficients, matrix.nrows, &matrix.field),
            relations: unpack_rows(&self.relations, matrix.nrows, &matrix.field),
        }
    }
}

fn xor_row(target: &mut [u64], source: &[u64]) {
    for (target, source) in target.iter_mut().zip(source) {
        *target ^= source;
    }
}

fn unpack_rows<F: Field>(rows: &[Vec<u64>], width: usize, field: &F) -> DenseMatrix<F> {
    let mut data = Vec::with_capacity(rows.len() * width);
    for row in rows {
        data.extend((0..width).map(|column| {
            if row[column / 64] & (1 << (column % 64)) == 0 {
                field.zero()
            } else {
                field.one()
            }
        }));
    }
    DenseMatrix {
        nrows: rows.len(),
        ncols: width,
        data,
        field: field.clone(),
    }
}
