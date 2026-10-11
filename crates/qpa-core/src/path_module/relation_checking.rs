// Local to one validation call: no cross-module or algebra reference cycles.
// Bound both prefix metadata and retained coefficient slots. Arbitrary-precision
// scalar payloads are additional; this is not a hard RSS bound.
#[derive(Clone)]
enum RelationMatrix<'a, F: Field> {
    Arrow(&'a DenseMatrix<F>),
    Product(Arc<DenseMatrix<F>>),
}

impl<F: Field> RelationMatrix<'_, F> {
    fn matrix(&self) -> &DenseMatrix<F> {
        match self {
            Self::Arrow(matrix) => matrix,
            Self::Product(matrix) => matrix,
        }
    }
}

#[derive(Clone)]
struct RelationPathAction<'a, F: Field> {
    shape: (usize, usize),
    field: F,
    matrix: Option<RelationMatrix<'a, F>>, // None denotes an exactly zero matrix.
}

impl<'a, F: Field> RelationPathAction<'a, F> {
    fn from_arrow(matrix: &'a DenseMatrix<F>, zero: bool) -> Self {
        Self {
            shape: (matrix.nrows(), matrix.ncols()),
            field: matrix.field().clone(),
            matrix: (!zero).then_some(RelationMatrix::Arrow(matrix)),
        }
    }

    fn from_matrix(matrix: DenseMatrix<F>) -> Self {
        Self {
            shape: (matrix.nrows(), matrix.ncols()),
            field: matrix.field().clone(),
            matrix: (!matrix.is_zero()).then(|| RelationMatrix::Product(Arc::new(matrix))),
        }
    }

    fn multiply(&self, right: &DenseMatrix<F>, right_zero: bool) -> Result<Self, MatrixError> {
        // Keep even NC constructor field/shape failures, including zero products.
        if &self.field != right.field() {
            return Err(MatrixError::DifferentFields);
        }
        if self.shape.1 != right.nrows() {
            return Err(MatrixError::IncompatibleProduct {
                left: self.shape,
                right: (right.nrows(), right.ncols()),
            });
        }
        if let Some(left) = &self.matrix
            && !right_zero
        {
            return Ok(Self::from_matrix(left.matrix().multiply(right)?));
        }
        Ok(Self {
            shape: (self.shape.0, right.ncols()),
            field: self.field.clone(),
            matrix: None,
        })
    }
}

#[allow(clippy::too_many_lines)]
fn relations_are_satisfied<F: Field>(
    algebra: &QuiverAlgebra<F>,
    dimension_vector: &[usize],
    arrow_matrices: &[DenseMatrix<F>],
    field: &F,
) -> Result<bool, MatrixError> {
    use std::collections::HashMap;
    const MAX_DEPTH: usize = 16;
    const MAX_PREFIXES: usize = 4096;
    const MAX_VALUES: usize = 1024;
    if algebra.relations().is_empty() {
        return Ok(true);
    }
    let mut frequency = HashMap::<Vec<ArrowIndex>, usize>::new();
    for relation in algebra.relations() {
        for term in relation.terms() {
            if let Some(walk) = term.path().walk() {
                for end in 2..=walk.len().min(MAX_DEPTH) {
                    if let Some(count) = frequency.get_mut(&walk[..end]) {
                        *count += 1;
                    } else if frequency.len() < MAX_PREFIXES {
                        frequency.insert(walk[..end].to_vec(), 1);
                    }
                }
            }
        }
    }
    let zero_arrows = arrow_matrices
        .iter()
        .map(DenseMatrix::is_zero)
        .collect::<Vec<_>>();
    let mut cache = HashMap::<Vec<ArrowIndex>, RelationPathAction<'_, F>>::new();
    let mut remaining = (8 * 1024 * 1024) / std::mem::size_of::<F::Element>().max(1);
    for relation in algebra.relations() {
        let mut shape_and_field: Option<((usize, usize), F)> = None;
        let mut value: Option<DenseMatrix<F>> = None;
        for term in relation.terms() {
            let walk = term
                .path()
                .walk()
                .expect("relation terms are nonzero paths");
            let action = if let Some(first) = walk.first() {
                let cached = (2..=walk.len().min(MAX_DEPTH))
                    .rev()
                    .find_map(|end| cache.get(&walk[..end]).map(|action| (end, action.clone())));
                let (start, mut action) = cached.unwrap_or_else(|| {
                    (
                        1,
                        RelationPathAction::from_arrow(
                            &arrow_matrices[first.zero_based()],
                            zero_arrows[first.zero_based()],
                        ),
                    )
                });
                for end in start + 1..=walk.len() {
                    let index = walk[end - 1].zero_based();
                    action = action.multiply(&arrow_matrices[index], zero_arrows[index])?;
                    if end <= MAX_DEPTH
                        && cache.len() < MAX_VALUES
                        && frequency.get(&walk[..end]).is_some_and(|&n| n > 1)
                    {
                        let slots = if action.matrix.is_some() {
                            action.shape.0.saturating_mul(action.shape.1)
                        } else {
                            0
                        };
                        if slots <= remaining {
                            remaining -= slots;
                            cache.insert(walk[..end].to_vec(), action.clone());
                        }
                    }
                }
                action
            } else {
                let PathEndpoint::Vertex(vertex) = term.path().source() else {
                    unreachable!()
                };
                RelationPathAction::from_matrix(DenseMatrix::identity(
                    physical_dimension(dimension_vector[vertex.index().zero_based()]),
                    field.clone(),
                ))
            };
            if let Some((shape, previous_field)) = &shape_and_field {
                if previous_field != &action.field {
                    return Err(MatrixError::DifferentFields);
                }
                if *shape != action.shape {
                    return Err(MatrixError::DifferentShapes {
                        left: *shape,
                        right: action.shape,
                    });
                }
            } else {
                shape_and_field = Some((action.shape, action.field.clone()));
            }
            if let Some(matrix) = action.matrix {
                let sum = value.get_or_insert_with(|| {
                    DenseMatrix::zero(action.shape.0, action.shape.1, action.field)
                });
                sum.add_scaled_assign(matrix.matrix(), term.coefficient())?;
            }
        }
        if value.is_some_and(|matrix| !matrix.is_zero()) {
            return Ok(false);
        }
    }
    Ok(true)
}
