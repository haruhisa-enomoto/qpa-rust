//! QPA's legacy `present2` projective-presentation object.
//!
//! The constructor in QPA 1.37 is observably broken: after finding the source
//! vertex of a nonzero row it calls the list form of `IsLeftUniform`, whose
//! implementation applies `CoefficientsAndMagmaElements` to that quiver
//! vertex. GAP has no such method, so no ordinary nonempty presentation matrix
//! reaches object creation. This module preserves that boundary and separately
//! exposes the declared attribute-storing representation.

use numerica::domains::Field;
use thiserror::Error;

use crate::{PathAlgebra, PathAlgebraElement, PathEndpoint, Vertex};

/// QPA's `IsProjectivePresentationDefaultRep` positional/attribute object.
#[derive(Clone, Debug)]
pub struct ProjectivePresentation<F: Field> {
    parent_algebra: PathAlgebra<F>,
    matrix_representation: Vec<Vec<PathAlgebraElement<F>>>,
    matrix_num_rows: usize,
    matrix_num_cols: usize,
    p0_vertex_list: Vec<Vertex>,
    p1_vertex_list: Option<Vec<Vertex>>,
}

impl<F: Field> ProjectivePresentation<F> {
    /// Translate QPA's public `ProjectivePresentation(A, mat)` method.
    ///
    /// In pinned QPA this returns no object: empty/empty-row lists fail method
    /// dispatch, ragged matrices reach the explicit row-length error, an
    /// all-zero first row passes `fail` to `CoefficientsAndMagmaElements`, and
    /// every nonzero first row hits the pinned vertex-uniformity dispatch bug.
    ///
    /// # Errors
    ///
    /// Returns the exact source boundary selected by the matrix shape and first
    /// row.
    pub fn from_matrix(
        _parent_algebra: &PathAlgebra<F>,
        matrix: &[Vec<PathAlgebraElement<F>>],
    ) -> Result<Self, ProjectivePresentationError> {
        let Some(first_row) = matrix.first() else {
            return Err(ProjectivePresentationError::NotMatrixDispatch);
        };
        if first_row.is_empty() || matrix.iter().any(Vec::is_empty) {
            return Err(ProjectivePresentationError::NotMatrixDispatch);
        }
        let number_of_columns = first_row.len();
        if matrix
            .iter()
            .skip(1)
            .any(|row| row.len() != number_of_columns)
        {
            return Err(ProjectivePresentationError::RowsOfUnequalLength);
        }

        let first_nonzero = first_row
            .iter()
            .find(|entry| !entry.is_zero())
            .ok_or(ProjectivePresentationError::AllZeroFirstRow)?;
        let first_term = first_nonzero
            .terms()
            .first()
            .ok_or(ProjectivePresentationError::AllZeroFirstRow)?;
        if !matches!(first_term.path().source(), PathEndpoint::Vertex(_)) {
            return Err(ProjectivePresentationError::AllZeroFirstRow);
        }

        Err(ProjectivePresentationError::VertexUniformityDispatch)
    }

    /// Construct the declared attribute-storing representation directly.
    ///
    /// This mirrors QPA's public `ObjectifyWithAttributes` route and therefore
    /// stores the supplied row/column counts independently of the matrix. The
    /// broken public constructor supplies no reachable object of its own.
    #[must_use]
    pub const fn from_attributes(
        parent_algebra: PathAlgebra<F>,
        matrix_representation: Vec<Vec<PathAlgebraElement<F>>>,
        matrix_num_rows: usize,
        matrix_num_cols: usize,
        p0_vertex_list: Vec<Vertex>,
        p1_vertex_list: Option<Vec<Vertex>>,
    ) -> Self {
        Self {
            parent_algebra,
            matrix_representation,
            matrix_num_rows,
            matrix_num_cols,
            p0_vertex_list,
            p1_vertex_list,
        }
    }

    /// Return QPA's `ParentAlgebra` attribute.
    #[must_use]
    pub const fn parent_algebra(&self) -> &PathAlgebra<F> {
        &self.parent_algebra
    }

    /// Return QPA's `MatrixRepresentation` attribute.
    #[must_use]
    pub fn matrix_representation(&self) -> &[Vec<PathAlgebraElement<F>>] {
        &self.matrix_representation
    }

    /// Return QPA's independently stored `MatrixNumRows` attribute.
    #[must_use]
    pub const fn matrix_num_rows(&self) -> usize {
        self.matrix_num_rows
    }

    /// Return QPA's independently stored `MatrixNumCols` attribute.
    #[must_use]
    pub const fn matrix_num_cols(&self) -> usize {
        self.matrix_num_cols
    }

    /// Return QPA's `P0VertexList` attribute.
    #[must_use]
    pub fn p0_vertex_list(&self) -> &[Vertex] {
        &self.p0_vertex_list
    }

    /// Return QPA's optional `P1VertexList` attribute.
    ///
    /// `None` corresponds to `HasP1VertexList(obj) = false`; the pinned
    /// constructor never sets this declared attribute.
    #[must_use]
    pub fn p1_vertex_list(&self) -> Option<&[Vertex]> {
        self.p1_vertex_list.as_deref()
    }
}

/// Failure in QPA's legacy `ProjectivePresentation` constructor.
#[derive(Clone, Copy, Debug, Eq, Error, PartialEq)]
pub enum ProjectivePresentationError {
    /// Empty matrices and matrices with an empty row do not satisfy `IsMatrix`.
    #[error("no ProjectivePresentation method applies to this matrix list")]
    NotMatrixDispatch,
    /// QPA's explicit rectangularity check.
    #[error("Matrix must have rows of equal length.")]
    RowsOfUnequalLength,
    /// `First(row, x -> not IsZero(x))` returned `fail`.
    #[error("QPA passes fail to CoefficientsAndMagmaElements for an all-zero row")]
    AllZeroFirstRow,
    /// The pinned list-uniformity method applies coefficients to a quiver vertex.
    #[error("QPA has no CoefficientsAndMagmaElements method for a quiver vertex")]
    VertexUniformityDispatch,
}

#[cfg(test)]
mod tests {
    use numerica::domains::{integer::Integer, rational::Q};
    use serde_json::{Value, json};

    use crate::{ArrowIndex, ArrowSpec, PathAlgebra, PathAlgebraElement, Quiver, VertexIndex};

    use super::{ProjectivePresentation, ProjectivePresentationError};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn arrow(index: usize) -> ArrowIndex {
        ArrowIndex::new(index).unwrap()
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the present2 fixture uses only small integers")
    }

    fn element_json(
        element: &PathAlgebraElement<numerica::domains::rational::RationalField>,
    ) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().to_string(),
                        "coefficient": [
                            small_integer(term.coefficient().numerator_ref()),
                            small_integer(term.coefficient().denominator_ref())
                        ]
                    })
                })
                .collect(),
        )
    }

    fn fixture_algebra() -> (
        Quiver,
        PathAlgebra<numerica::domains::rational::RationalField>,
    ) {
        let quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(1), vertex(3), "b"),
                ArrowSpec::named(vertex(2), vertex(3), "c"),
            ],
        )
        .unwrap();
        let algebra = PathAlgebra::new(Q, quiver.clone()).unwrap();
        (quiver, algebra)
    }

    #[test]
    fn pinned_qpa_present2_representation_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/present2.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1"
            })
        );

        let (quiver, algebra) = fixture_algebra();
        let a = algebra
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let b = algebra
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();
        let c = algebra
            .element(quiver.arrow(arrow(3)).unwrap().path())
            .unwrap();
        let v2 = algebra
            .element(quiver.vertex(vertex(2)).unwrap().path())
            .unwrap();
        let matrix = vec![vec![a.add(&b).unwrap(), algebra.zero()], vec![c, v2]];
        let presentation = ProjectivePresentation::from_attributes(
            algebra.clone(),
            matrix.clone(),
            2,
            2,
            vec![
                quiver.vertex(vertex(1)).unwrap(),
                quiver.vertex(vertex(2)).unwrap(),
            ],
            None,
        );
        assert_eq!(expected["category"], true);
        assert_eq!(expected["default_rep"], true);
        assert_eq!(expected["global_family"], true);
        assert_eq!(
            json!({
                "parent_identical": presentation.parent_algebra().is_identical(&algebra),
                "matrix_equal": presentation.matrix_representation() == matrix,
                "matrix": presentation.matrix_representation().iter().map(|row| row.iter().map(element_json).collect::<Vec<_>>()).collect::<Vec<_>>(),
                "shape": [presentation.matrix_num_rows(), presentation.matrix_num_cols()],
                "p0": presentation.p0_vertex_list().iter().map(|vertex| vertex.name().to_owned()).collect::<Vec<_>>(),
                "has_p1": presentation.p1_vertex_list().is_some()
            }),
            json!({
                "parent_identical": expected["parent_identical"],
                "matrix_equal": expected["matrix_equal"],
                "matrix": expected["matrix"],
                "shape": expected["shape"],
                "p0": expected["p0"],
                "has_p1": expected["has_p1"]
            })
        );
    }

    #[test]
    fn pinned_constructor_failures_remain_unrepaired() {
        let (quiver, algebra) = fixture_algebra();
        let a = algebra
            .element(quiver.arrow(arrow(1)).unwrap().path())
            .unwrap();
        let b = algebra
            .element(quiver.arrow(arrow(2)).unwrap().path())
            .unwrap();

        assert_eq!(
            ProjectivePresentation::from_matrix(&algebra, &[vec![a.clone(), b.clone()]])
                .unwrap_err(),
            ProjectivePresentationError::VertexUniformityDispatch
        );
        assert_eq!(
            ProjectivePresentation::from_matrix(&algebra, &[vec![algebra.zero()]]).unwrap_err(),
            ProjectivePresentationError::AllZeroFirstRow
        );
        assert_eq!(
            ProjectivePresentation::from_matrix(&algebra, &[vec![a.clone()], vec![a, b]])
                .unwrap_err(),
            ProjectivePresentationError::RowsOfUnequalLength
        );
        assert_eq!(
            ProjectivePresentation::from_matrix(&algebra, &[]).unwrap_err(),
            ProjectivePresentationError::NotMatrixDispatch
        );
        assert_eq!(
            ProjectivePresentation::from_matrix(&algebra, &[vec![]]).unwrap_err(),
            ProjectivePresentationError::NotMatrixDispatch
        );
        assert_eq!(
            ProjectivePresentationError::RowsOfUnequalLength.to_string(),
            "Matrix must have rows of equal length."
        );
    }
}
