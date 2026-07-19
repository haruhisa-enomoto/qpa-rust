//! Construct a quiver, its rational path algebra, and a right matrix module.

use std::{error::Error, io};

use qpa_core::api::prelude::*;

fn vertex(index: usize) -> VertexIndex {
    VertexIndex::new(index).expect("example uses positive one-based indices")
}

fn main() -> Result<(), Box<dyn Error>> {
    // QPA convention: vertices and arrows are numbered from one.
    let quiver = Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")])?;
    let path_algebra = PathAlgebra::new(Q, quiver.clone())?;

    // Path products are written in traversal order: a * b means first a,
    // then b. Here the identity acts on the arrow element as expected.
    let arrow = quiver
        .arrows()
        .next()
        .expect("the example quiver has one arrow");
    let a = path_algebra.element(arrow.path())?;
    assert_eq!(a.multiply(&path_algebra.one())?, a);

    // Right modules use row vectors, so an arrow i -> j has a
    // dim(M_i)-by-dim(M_j) matrix.
    let action = DenseMatrix::from_rows(vec![vec![Q.one()]], Q)?;
    let module = right_module_over_path_algebra(
        QuiverAlgebra::Path(path_algebra),
        vec![1, 1],
        &[NamedArrowMatrix::new("a", action)],
    )?
    .ok_or_else(|| io::Error::other("the entered actions violated a relation"))?;

    assert_eq!(module.dimension_vector(), &[1, 1]);
    println!("module dimension = {}", module.dimension());
    Ok(())
}
