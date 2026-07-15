//! Centers of QPA path algebras and finite-dimensional quotients.
//!
//! The two methods deliberately follow the separate implementations in
//! `pathalg.gi`: parent path algebras use connected components and oriented
//! cycles, while quotients solve the arrow-commutator coefficient matrix.

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    AlgebraDimension, DenseMatrix, MatrixError, PathAlgebra, PathAlgebraElement, PathAlgebraError,
    PathEndpoint, QuiverError, QuotientError, QuotientPathAlgebra, QuotientPathAlgebraElement,
};

/// QPA's center subalgebra for an unfactored path algebra.
#[derive(Clone, Debug)]
pub struct PathAlgebraCentre<F: Field> {
    ambient: PathAlgebra<F>,
    generators: Vec<PathAlgebraElement<F>>,
    dimension: AlgebraDimension,
}

impl<F: Field> PathAlgebraCentre<F> {
    /// Return the ambient path algebra.
    #[must_use]
    pub const fn ambient(&self) -> &PathAlgebra<F> {
        &self.ambient
    }

    /// Return QPA's `GeneratorsOfAlgebraWithOne` list for the center.
    #[must_use]
    pub fn generators(&self) -> &[PathAlgebraElement<F>] {
        &self.generators
    }

    /// Return the center's vector-space dimension, including infinity for an
    /// oriented-cycle component.
    #[must_use]
    pub fn dimension(&self) -> &AlgebraDimension {
        &self.dimension
    }

    /// Test membership in the center by commuting with QPA's ambient
    /// vertex-and-arrow generator list.
    ///
    /// # Errors
    ///
    /// Returns an error for an element from another path-algebra family or a
    /// failed path product.
    pub fn contains(&self, element: &PathAlgebraElement<F>) -> Result<bool, CentreError> {
        if !self.ambient().is_identical(element.algebra()) {
            return Err(CentreError::DifferentAmbientAlgebra);
        }
        for generator in self.ambient().generators() {
            if element.multiply(&generator)? != generator.multiply(element)? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

/// QPA's known-basis center subalgebra for a quotient path algebra.
#[derive(Clone, Debug)]
pub struct QuotientPathAlgebraCentre<F: Field> {
    ambient: QuotientPathAlgebra<F>,
    basis: Vec<QuotientPathAlgebraElement<F>>,
}

impl<F: Field> QuotientPathAlgebraCentre<F> {
    /// Return the ambient quotient path algebra.
    #[must_use]
    pub const fn ambient(&self) -> &QuotientPathAlgebra<F> {
        &self.ambient
    }

    /// Return QPA's center basis supplied to `SubalgebraWithOne(...,
    /// "basis")`.
    #[must_use]
    pub fn basis(&self) -> &[QuotientPathAlgebraElement<F>] {
        &self.basis
    }

    /// Return the exact finite center dimension.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.basis.len()
    }

    /// Test membership in the center by commuting with QPA's ambient
    /// vertex-and-arrow generator list.
    ///
    /// # Errors
    ///
    /// Returns an error for another quotient family or failed normal-form
    /// arithmetic.
    pub fn contains(&self, element: &QuotientPathAlgebraElement<F>) -> Result<bool, CentreError> {
        if !self.ambient().is_identical(element.algebra()) {
            return Err(CentreError::DifferentAmbientAlgebra);
        }
        for generator in self.ambient().generators()? {
            if element.multiply(&generator)? != generator.multiply(element)? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

impl<F: Field> PathAlgebra<F> {
    /// Translate QPA's component/oriented-cycle `Centre` method for a path
    /// algebra.
    ///
    /// # Errors
    ///
    /// Propagates QPA connected-component, name-resolution, and path-product
    /// failures.
    pub fn centre(&self) -> Result<PathAlgebraCentre<F>, CentreError> {
        let quiver = self.quiver();
        let original_vertices = quiver.vertices().collect::<Vec<_>>();
        let original_arrows = quiver.arrows().collect::<Vec<_>>();
        let mut generators = Vec::new();
        let mut has_oriented_cycle = false;

        for component in quiver.connected_components()? {
            let component_vertices = component.vertices().collect::<Vec<_>>();
            let mut block_identity = self.zero();
            for component_vertex in &component_vertices {
                let original = original_vertices
                    .iter()
                    .find(|vertex| vertex.name() == component_vertex.name())
                    .ok_or(CentreError::ComponentNameMapping)?;
                block_identity = block_identity.add(&self.element(original.path())?)?;
            }
            generators.push(block_identity);

            let oriented_cycle = component_vertices
                .iter()
                .all(|vertex| vertex.in_degree() == 1 && vertex.out_degree() == 1);
            if !oriented_cycle {
                continue;
            }
            let component_arrows = component.arrows().collect::<Vec<_>>();
            let comparison_arrows = component_arrows
                .iter()
                .map(|component_arrow| {
                    original_arrows
                        .iter()
                        .find(|arrow| arrow.name() == component_arrow.name())
                        .cloned()
                        .ok_or(CentreError::ComponentNameMapping)
                })
                .collect::<Result<Vec<_>, CentreError>>()?;
            let Some(first) = comparison_arrows.first().cloned() else {
                continue;
            };

            let mut cycle = vec![first.clone()];
            let mut last = first;
            for _ in 1..comparison_arrows.len() {
                let target = last.target();
                last = comparison_arrows
                    .iter()
                    .find(|arrow| arrow.source() == target)
                    .cloned()
                    .ok_or(CentreError::BrokenOrientedCycle)?;
                cycle.push(last.clone());
            }

            let mut cycle_sum = self.zero();
            let length = cycle.len();
            for rotation in 0..length {
                let mut rotated_cycle = self.one();
                for offset in 0..length {
                    let arrow = &cycle[(rotation + offset) % length];
                    rotated_cycle = rotated_cycle.multiply(&self.element(arrow.path())?)?;
                }
                cycle_sum = cycle_sum.add(&rotated_cycle)?;
            }
            generators.push(cycle_sum);
            has_oriented_cycle = true;
        }

        let dimension = if has_oriented_cycle {
            AlgebraDimension::Infinite
        } else {
            AlgebraDimension::Finite(generators.len().into())
        };
        Ok(PathAlgebraCentre {
            ambient: self.clone(),
            generators,
            dimension,
        })
    }
}

impl<F: Field> QuotientPathAlgebra<F> {
    /// Translate QPA's finite-dimensional quotient `Centre` method using its
    /// arrow-commutator coefficient matrix and GAP `NullspaceMat` convention.
    ///
    /// # Errors
    ///
    /// The quotient must be finite dimensional and its cached canonical basis
    /// must contain normal-form vectors. Arithmetic, coordinate, matrix, and
    /// nontip failures are propagated.
    pub fn centre(&self) -> Result<QuotientPathAlgebraCentre<F>, CentreError> {
        if !self.is_finite_dimensional()? {
            return Err(CentreError::InfiniteDimensionalQuotient);
        }
        let canonical_basis = self.canonical_basis_object()?;
        let basis = canonical_basis
            .basis_vectors()
            .ok_or(QuotientError::CanonicalBasisVectorsUnavailable)?;
        let cycles = basis
            .iter()
            .filter(|element| {
                matches!(
                    (
                        element.leading_monomial().source(),
                        element.leading_monomial().target()
                    ),
                    (PathEndpoint::Vertex(source), PathEndpoint::Vertex(target)) if source == target
                )
            })
            .cloned()
            .collect::<Vec<_>>();
        let generators = self.generators()?;
        let vertex_count = self.original_path_algebra().quiver().number_of_vertices();
        let arrow_count = self.original_path_algebra().quiver().number_of_arrows();
        let arrows = &generators[vertex_count..vertex_count + arrow_count];
        let column_count = basis.len().saturating_mul(arrow_count);
        let mut entries = Vec::with_capacity(cycles.len().saturating_mul(column_count));

        for cycle in &cycles {
            for arrow in arrows {
                let commutator = cycle.multiply(arrow)?.sub(&arrow.multiply(cycle)?)?;
                entries.extend(canonical_basis.coefficients(&commutator)?);
            }
        }
        let matrix = DenseMatrix::from_flat(
            cycles.len(),
            column_count,
            entries,
            self.original_path_algebra().field().clone(),
        )?;
        let solutions = matrix.nullspace_rows();
        let field = self.original_path_algebra().field();
        let mut center_basis = Vec::with_capacity(solutions.nrows());
        for row in 0..solutions.nrows() {
            let mut element = self.zero();
            for (coefficient, cycle) in solutions.row(row).iter().zip(&cycles) {
                if !field.is_zero(coefficient) {
                    element = element.add(&cycle.scale_left(coefficient)?)?;
                }
            }
            center_basis.push(element);
        }
        Ok(QuotientPathAlgebraCentre {
            ambient: self.clone(),
            basis: center_basis,
        })
    }
}

/// Failure in QPA's two path-algebra center constructions.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum CentreError {
    /// A membership candidate belongs to another ambient algebra family.
    #[error("the element belongs to a different ambient algebra")]
    DifferentAmbientAlgebra,
    /// The quotient method rejects an infinite-dimensional algebra.
    #[error("the entered quotient algebra is not finite dimensional")]
    InfiniteDimensionalQuotient,
    /// QPA's component-name `Position` mapping had no original generator.
    #[error("a connected-component generator name has no ambient match")]
    ComponentNameMapping,
    /// A component passed QPA's degree-one test but cycle following failed.
    #[error("the degree-one component does not contain the expected oriented cycle")]
    BrokenOrientedCycle,
    /// Parent path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Quotient arithmetic or basis conversion failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Connected-component construction failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// GAP-compatible nullspace matrix construction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
}
