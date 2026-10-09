//! Regression checks against the original full GAP transformation translation.

use qpa_core::api::coefficients::{Field, Integer, Q, Ring, Z2, Zp};
use qpa_core::{DenseMatrix, GaloisField, MatrixError, SemiEchelonTransformation};

fn rows<F: Field>(values: Vec<Vec<F::Element>>, width: usize, field: F) -> DenseMatrix<F> {
    DenseMatrix::from_flat(
        values.len(),
        width,
        values.into_iter().flatten().collect(),
        field,
    )
    .unwrap()
}

// Intentionally retain full rows and transformations here. This is the original
// algorithm, including its column-order reductions and choice of dependencies.
fn reference<F: Field>(input: &DenseMatrix<F>) -> SemiEchelonTransformation<F> {
    let field = input.field();
    let mut matrix = input.clone();
    let mut transform = DenseMatrix::identity(input.nrows(), field.clone());
    let mut heads: Vec<Option<usize>> = vec![None; input.ncols()];
    let mut vectors: Vec<Vec<F::Element>> = Vec::new();
    let mut coefficients: Vec<Vec<F::Element>> = Vec::new();
    let mut relations = Vec::new();
    for row in 0..input.nrows() {
        for column in 0..input.ncols() {
            if let Some(head) = heads[column] {
                let scale = field.neg(&matrix[(row, column)]);
                if !field.is_zero(&scale) {
                    for (index, value) in coefficients[head].iter().enumerate() {
                        field.add_assign(&mut transform[(row, index)], &field.mul(value, &scale));
                    }
                    for (index, value) in vectors[head].iter().enumerate() {
                        field.add_assign(&mut matrix[(row, index)], &field.mul(value, &scale));
                    }
                }
            }
        }
        if let Some(column) = matrix.row(row).iter().position(|x| !field.is_zero(x)) {
            let inverse = field.inv(&matrix[(row, column)]);
            for index in 0..input.nrows() {
                field.mul_assign(&mut transform[(row, index)], &inverse);
            }
            for index in 0..input.ncols() {
                field.mul_assign(&mut matrix[(row, index)], &inverse);
            }
            coefficients.push(transform.row(row).to_vec());
            vectors.push(matrix.row(row).to_vec());
            heads[column] = Some(vectors.len() - 1);
        } else {
            relations.push(transform.row(row).to_vec());
        }
    }
    SemiEchelonTransformation {
        heads,
        vectors: rows(vectors, input.ncols(), field.clone()),
        coefficients: rows(coefficients, input.nrows(), field.clone()),
        relations: rows(relations, input.nrows(), field.clone()),
    }
}

fn check_field<F: Field>(field: &F, generator: &F::Element) {
    for height in 0..7 {
        for width in 0..7 {
            for seed in 0..4 {
                let data = (0..height)
                    .flat_map(|row| {
                        (0..width).map(move |column| {
                            field.add(
                                &field.nth(Integer::from(
                                    (row * 13 + column * 7 + seed * 11 + row * column * 3) % 9,
                                )),
                                &field.mul(
                                    generator,
                                    &field.nth(Integer::from((row + column + seed) % 3)),
                                ),
                            )
                        })
                    })
                    .collect();
                let matrix = DenseMatrix::from_flat(height, width, data, field.clone()).unwrap();
                let expected = reference(&matrix);
                let actual = matrix.semi_echelon_transformation();
                assert_eq!(actual, expected);
                assert_eq!(matrix.semi_echelon_basis(), expected.vectors);
                assert_eq!(matrix.rank(), expected.vectors.nrows());
                assert_eq!(
                    actual.coefficients.multiply(&matrix).unwrap(),
                    actual.vectors
                );
                assert_eq!(
                    actual.relations.multiply(&matrix).unwrap(),
                    DenseMatrix::zero(actual.relations.nrows(), width, field.clone())
                );
                for target in
                    (0..height)
                        .map(|row| matrix.row(row).to_vec())
                        .chain([vec![field.zero(); width]])
                {
                    let solution = actual.solution_row(&target).unwrap().unwrap();
                    let row = DenseMatrix::from_flat(1, height, solution, field.clone()).unwrap();
                    assert_eq!(row.multiply(&matrix).unwrap().row(0), target);
                }
                assert!(matches!(
                    actual.solution_row(&vec![field.zero(); width + 1]),
                    Err(MatrixError::TargetLength { .. })
                ));
            }
        }
    }
}

#[test]
fn elimination_preserves_exact_bases_across_fields_and_empty_shapes() {
    check_field(&Q, &Q.one());
    check_field(&Z2, &Z2.one());
    let prime = Zp::new(5);
    check_field(&prime, &prime.one());
    let extension = GaloisField::conway(2, 2).unwrap();
    check_field(&extension, &extension.primitive_element());
}

#[test]
fn prepared_solutions_preserve_gap_choice_and_failure() {
    let matrix = DenseMatrix::from_rows(
        vec![
            vec![Q.zero(), Q.nth(2.into()), Q.nth(2.into())],
            vec![Q.nth(3.into()), Q.zero(), Q.zero()],
            vec![Q.zero(), Q.nth(4.into()), Q.nth(4.into())],
        ],
        Q,
    )
    .unwrap();
    let semi = matrix.semi_echelon_transformation();
    assert_eq!(semi.heads, vec![Some(1), Some(0), None]);
    assert_eq!(semi, reference(&matrix));
    assert_eq!(matrix.semi_echelon_basis(), semi.vectors);
    assert_eq!(
        semi.solution_row(&[Q.nth(3.into()), Q.nth(2.into()), Q.nth(2.into())])
            .unwrap(),
        Some(vec![Q.one(), Q.one(), Q.zero()])
    );
    assert_eq!(
        semi.solution_row(&[Q.zero(), Q.zero(), Q.one()]).unwrap(),
        None
    );
}
