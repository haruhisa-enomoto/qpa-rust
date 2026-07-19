//! Compute the dual numbers as a quotient through QPA's GBNP call chain.

use std::error::Error;

use qpa_gbnp::api::prelude::*;

fn vertex(index: usize) -> VertexIndex {
    VertexIndex::new(index).expect("example uses positive one-based indices")
}

fn main() -> Result<(), Box<dyn Error>> {
    let quiver = Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")])?;
    let path_algebra = PathAlgebra::new(Z2, quiver.clone())?;

    let x_path = quiver
        .arrows()
        .next()
        .expect("the example quiver has one loop")
        .path();
    let x = path_algebra.element(x_path)?;
    let x_squared = x.multiply(&x)?;

    // This is QPA's KQ / [x^2] call through its configured GBNP backend.
    let quotient = quotient_by_relations(path_algebra, &[x_squared])?;
    let dimension = quotient_dimension(&quotient)?;

    let AlgebraDimension::Finite(dimension) = dimension else {
        panic!("the dual numbers are finite-dimensional");
    };
    assert_eq!(dimension.to_string(), "2");
    println!("dimension of F_2[x]/(x^2) = {dimension}");
    Ok(())
}
