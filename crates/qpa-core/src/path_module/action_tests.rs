use super::*;
use crate::{PathAlgebra, Quiver};
use numerica::domains::{
    finite_field::{Z2, Zp},
    rational::Q,
};

#[allow(clippy::many_single_char_names)]
fn compare_actions<F: Field>(field: &F) {
    let v = crate::VertexIndex::new(1).unwrap();
    let quiver = Quiver::from_indexed(
        1,
        vec![
            crate::ArrowSpec::named(v, v, "x"),
            crate::ArrowSpec::named(v, v, "y"),
        ],
    )
    .unwrap();
    let parent = PathAlgebra::new(field.clone(), quiver).unwrap();
    let z = field.zero();
    let o = field.one();
    let x = DenseMatrix::from_rows(
        vec![
            vec![z.clone(), o.clone(), z.clone()],
            vec![z.clone(), z.clone(), o.clone()],
            vec![z.clone(), z.clone(), z.clone()],
        ],
        field.clone(),
    )
    .unwrap();
    let y = DenseMatrix::from_rows(
        vec![
            vec![o.clone(), z.clone(), o.clone()],
            vec![z.clone(), o.clone(), z.clone()],
            vec![o.clone(), o.clone(), z.clone()],
        ],
        field.clone(),
    )
    .unwrap();
    let module = right_module_over_path_algebra(
        QuiverAlgebra::Path(parent.clone()),
        vec![3],
        &[NamedArrowMatrix::new("x", x), NamedArrowMatrix::new("y", y)],
    )
    .unwrap()
    .unwrap();
    let row = DenseMatrix::from_rows(
        vec![vec![o.clone(), field.add(&o, &o), field.neg(&o)]],
        field.clone(),
    )
    .unwrap();
    let vector = module
        .element_from_vertex_coordinates(vec![row.row(0).to_vec()])
        .unwrap();
    let arrows = parent
        .quiver()
        .arrows()
        .map(|a| parent.element(a.path()).unwrap())
        .collect::<Vec<_>>();
    for degree in 0..8 {
        let word = arrows[0].multiply(&arrows[1]).unwrap().pow(degree).unwrap();
        let polynomial = word
            .sub(&arrows[1].pow(degree).unwrap())
            .unwrap()
            .add(&parent.one())
            .unwrap();
        let mut expected = DenseMatrix::zero(1, 3, field.clone());
        // Frozen pre-optimization matrix-chain evaluation.
        for term in polynomial.terms() {
            let action = path_action_matrix(
                term.path(),
                module.dimension_vector(),
                module.arrow_matrices(),
                field.clone(),
            )
            .unwrap();
            expected = expected
                .add(
                    &row.multiply(&action)
                        .unwrap()
                        .scale_left(term.coefficient()),
                )
                .unwrap();
        }
        let actual = vector.act(&QuiverAlgebraElement::Path(polynomial)).unwrap();
        assert_eq!(actual.external_representation()[0], expected.row(0));
    }
}

#[test]
fn vector_first_path_action_matches_full_matrix_chain() {
    compare_actions(&Q);
    compare_actions(&Z2);
    compare_actions(&Zp::new(101));
}
