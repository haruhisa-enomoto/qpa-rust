//! Projective path-algebra presentations of finite-dimensional modules.
//!
//! This is the source-order translation of QPA's
//! `ProjectivePathAlgebraPresentation` method in
//! `lib/projpathalgmodule.gi`. It deliberately computes the relation basis
//! through QPA's configured GBNP adapter and retains both the entered kernel
//! generators and the completely reduced module Groebner basis.

use numerica::domains::Field;
use qpa_core::{
    DenseMatrix, MatrixError, PathAlgebraElement, PathAlgebraIdeal, PathAlgebraModule,
    PathAlgebraModuleError, PathAlgebraVector, ProjectivePathModuleError, QuiverAlgebra,
    QuiverAlgebraElement, QuiverAlgebraElementError, QuotientError,
    RightProjectivePathAlgebraModule, RightProjectivePathAlgebraSubmodule,
    completely_reduce_path_algebra_module_groebner_basis, minimal_generating_set_of_module,
    right_projective_path_algebra_module,
};
use thiserror::Error;

use crate::{
    module::{EnvelopingModuleError, basis_of_projectives},
    qpa_conversion::{QpaGbnpError, gbnp_groebner_basis},
};

/// QPA's list-valued `ProjectivePathAlgebraPresentation(M)` result.
///
/// The source returns a list of length five when the kernel is empty and a
/// list of length six otherwise. Rust gives those positions names while
/// preserving the exceptional empty-list values through [`Self::result_length`].
#[derive(Clone, Debug)]
pub struct ProjectivePathAlgebraPresentation<F: Field> {
    projective: RightProjectivePathAlgebraModule<F>,
    submodule: Option<RightProjectivePathAlgebraSubmodule<F>>,
    projective_generators: Vec<PathAlgebraVector<F>>,
    syzygy_generators: Vec<PathAlgebraVector<F>>,
    matrix: Vec<Vec<PathAlgebraElement<F>>>,
    reduced_groebner_basis: Vec<PathAlgebraVector<F>>,
    relation_ideal: PathAlgebraIdeal<F>,
}

impl<F: Field> ProjectivePathAlgebraPresentation<F> {
    /// Return the inducing projective module over the original path algebra.
    #[must_use]
    pub const fn projective(&self) -> &RightProjectivePathAlgebraModule<F> {
        &self.projective
    }

    /// Return QPA's entered kernel submodule, absent in its empty-kernel branch.
    #[must_use]
    pub const fn submodule(&self) -> Option<&RightProjectivePathAlgebraSubmodule<F>> {
        self.submodule.as_ref()
    }

    /// Return QPA's exposed `f0` list.
    ///
    /// This is deliberately empty in QPA's length-five result, even though
    /// [`Self::projective`] still has its canonical coordinate generators.
    #[must_use]
    pub fn projective_generators(&self) -> &[PathAlgebraVector<F>] {
        &self.projective_generators
    }

    /// Return QPA's `f1`, reconstructed from its reduced raw vectors.
    #[must_use]
    pub fn syzygy_generators(&self) -> &[PathAlgebraVector<F>] {
        &self.syzygy_generators
    }

    /// Return QPA's transposed path-algebra matrix from `f1` to `f0`.
    #[must_use]
    pub fn matrix(&self) -> &[Vec<PathAlgebraElement<F>>] {
        &self.matrix
    }

    /// Return QPA's raw `redRG` positional vectors.
    #[must_use]
    pub fn reduced_groebner_basis(&self) -> &[PathAlgebraVector<F>] {
        &self.reduced_groebner_basis
    }

    /// Return the fresh defining ideal and attached GBNP basis used by QPA.
    #[must_use]
    pub const fn relation_ideal(&self) -> &PathAlgebraIdeal<F> {
        &self.relation_ideal
    }

    /// Return the observable GAP result-list length, either five or six.
    #[must_use]
    pub const fn result_length(&self) -> usize {
        if self.submodule.is_some() { 6 } else { 5 }
    }
}

/// Translate QPA's `ProjectivePathAlgebraPresentation(M)` method.
///
/// The implementation follows the installed GAP source: it recomputes a
/// complete defining-relation basis through GBNP, constructs the vector-space
/// kernel of the projective cover from `BasisOfProjectives`, appends the
/// relation-ideal generators, and completely reduces the resulting right
/// module Groebner basis.
///
/// # Errors
///
/// QPA 1.37 reads the quotient ideal before its nominal path-algebra branch,
/// so actual path-algebra input fails; Rust exposes that pinned boundary as
/// [`ProjectivePathAlgebraPresentationError::PathAlgebraHasNoDefiningIdeal`].
/// Other failures propagate the configured GBNP, finite projective-basis,
/// exact matrix/module arithmetic, quotient, and projective-vector machinery.
#[allow(clippy::too_many_lines)]
pub fn projective_path_algebra_presentation<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<ProjectivePathAlgebraPresentation<F>, ProjectivePathAlgebraPresentationError> {
    let algebra = module.right_acting_algebra();
    if matches!(algebra, QuiverAlgebra::Path(_)) {
        return Err(ProjectivePathAlgebraPresentationError::PathAlgebraHasNoDefiningIdeal);
    }

    let parent = algebra.original_path_algebra().clone();
    let relations = algebra.relations().to_vec();
    let complete_relations = gbnp_groebner_basis(&relations, &parent)?;
    let relation_ideal = PathAlgebraIdeal::new(parent.clone(), relations)?;
    relation_ideal.with_complete_groebner_basis(&complete_relations)?;
    let right_relations = relation_ideal.right_groebner_basis()?.relations().to_vec();

    let projective_bases = basis_of_projectives(algebra)?
        .ok_or(ProjectivePathAlgebraPresentationError::ProjectiveBasisUnavailable)?;
    let projective_bases = projective_bases
        .into_iter()
        .map(|projective| projective.into_iter().flatten().collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let module_generators = minimal_generating_set_of_module(module.clone())?;
    let projective_cover = module_generators
        .iter()
        .enumerate()
        .map(|(position, generator)| {
            let support = generator.support_vertices();
            let [vertex] = support.as_slice() else {
                return Err(
                    ProjectivePathAlgebraPresentationError::GeneratorNotVertexSupported {
                        position: position + 1,
                        support,
                    },
                );
            };
            Ok(*vertex)
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut cover_images = Vec::new();
    for (generator, &cover_vertex) in module_generators.iter().zip(&projective_cover) {
        let basis = &projective_bases[cover_vertex - 1];
        for basis_element in basis {
            cover_images.push(generator.act(basis_element)?);
        }
    }
    let cover_matrix = DenseMatrix::from_flat(
        cover_images.len(),
        module.dimension(),
        cover_images
            .iter()
            .flat_map(qpa_core::PathAlgebraModuleElement::basis_coefficients)
            .collect(),
        algebra.field().clone(),
    )?;
    let solutions = cover_matrix.nullspace_rows();

    let mut lifted_solutions = Vec::with_capacity(solutions.nrows());
    for solution_index in 0..solutions.nrows() {
        let solution = solutions.row(solution_index);
        let mut offset = 0;
        let mut components = Vec::with_capacity(projective_cover.len());
        for &cover_vertex in &projective_cover {
            let basis = &projective_bases[cover_vertex - 1];
            let end = offset + basis.len();
            let coefficients = solution
                .get(offset..end)
                .ok_or(ProjectivePathAlgebraPresentationError::KernelCoefficientCountMismatch)?;
            let quotient_element = linear_combination(algebra, basis, coefficients)?;
            components.push(parent_representative(quotient_element));
            offset = end;
        }
        if offset != solution.len() {
            return Err(ProjectivePathAlgebraPresentationError::KernelCoefficientCountMismatch);
        }
        lifted_solutions.push(components);
    }

    let parent_generators = parent.generators();
    let vertex_list = projective_cover
        .iter()
        .map(|&vertex| parent_generators[vertex - 1].clone())
        .collect::<Vec<_>>();
    let projective = right_projective_path_algebra_module(&parent, &vertex_list)?;
    let f0 = projective.generators()?;
    let mut syzygies = lifted_solutions
        .into_iter()
        .map(|components| projective.vector(components))
        .collect::<Result<Vec<_>, _>>()?;

    for generator in &f0 {
        for relation in &right_relations {
            let product = generator.act(relation)?;
            if !product.is_zero() {
                syzygies.push(product);
            }
        }
    }

    if syzygies.is_empty() {
        return Ok(ProjectivePathAlgebraPresentation {
            projective,
            submodule: None,
            projective_generators: Vec::new(),
            syzygy_generators: Vec::new(),
            matrix: Vec::new(),
            reduced_groebner_basis: Vec::new(),
            relation_ideal,
        });
    }

    let submodule = projective.submodule(syzygies)?;
    let right_basis = submodule.right_groebner_basis()?;
    let reduced_groebner_basis =
        completely_reduce_path_algebra_module_groebner_basis(&right_basis)?;
    let syzygy_generators = reduced_groebner_basis.clone();
    let matrix = (0..projective.vertex_indices().len())
        .map(|position| {
            syzygy_generators
                .iter()
                .map(|generator| generator.components()[position].clone())
                .collect()
        })
        .collect();

    Ok(ProjectivePathAlgebraPresentation {
        projective,
        submodule: Some(submodule),
        projective_generators: f0,
        syzygy_generators,
        matrix,
        reduced_groebner_basis,
        relation_ideal,
    })
}

fn linear_combination<F: Field>(
    algebra: &QuiverAlgebra<F>,
    basis: &[QuiverAlgebraElement<F>],
    coefficients: &[F::Element],
) -> Result<QuiverAlgebraElement<F>, ProjectivePathAlgebraPresentationError> {
    if basis.len() != coefficients.len() {
        return Err(ProjectivePathAlgebraPresentationError::KernelCoefficientCountMismatch);
    }
    let mut result = algebra.zero();
    for (element, coefficient) in basis.iter().zip(coefficients) {
        result = result.add(&element.scale_left(coefficient)?)?;
    }
    Ok(result)
}

fn parent_representative<F: Field>(element: QuiverAlgebraElement<F>) -> PathAlgebraElement<F> {
    match element {
        QuiverAlgebraElement::Path(element) => element,
        QuiverAlgebraElement::Quotient(element) => element.representative().clone(),
    }
}

/// Failure in QPA's path-algebra presentation construction.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ProjectivePathAlgebraPresentationError {
    /// The installed QPA 1.37 method reads `fam!.ideal` even for path algebras.
    #[error(
        "QPA's ProjectivePathAlgebraPresentation method cannot read a defining ideal from a path algebra"
    )]
    PathAlgebraHasNoDefiningIdeal,
    /// `BasisOfProjectives` returned QPA's printed `fail` value.
    #[error("BasisOfProjectives returned fail while constructing the projective presentation")]
    ProjectiveBasisUnavailable,
    /// QPA indexes the first support vertex of each minimal generator.
    #[error(
        "minimal module generator {position} has vertex support {support:?}, expected one vertex"
    )]
    GeneratorNotVertexSupported {
        /// One-based generator position.
        position: usize,
        /// Rejected one-based support positions.
        support: Vec<usize>,
    },
    /// The nullspace row did not split into the selected projective bases.
    #[error("the projective-cover nullspace row has the wrong coefficient count")]
    KernelCoefficientCountMismatch,
    /// QPA's configured GBNP adapter failed.
    #[error(transparent)]
    Gbnp(#[from] QpaGbnpError),
    /// QPA's finite projective-basis machinery failed.
    #[error(transparent)]
    Module(#[from] EnvelopingModuleError),
    /// Exact module construction or action failed.
    #[error(transparent)]
    PathModule(#[from] PathAlgebraModuleError),
    /// Exact vector-space kernel arithmetic failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// Defining-ideal or right-basis construction failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Arithmetic in the path-algebra/quotient union failed.
    #[error(transparent)]
    QuiverElement(#[from] QuiverAlgebraElementError),
    /// Projective path-vector construction or reduction failed.
    #[error(transparent)]
    ProjectiveModule(#[from] ProjectivePathModuleError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, Rational},
    };
    use qpa_core::{
        ArrowSpec, DenseMatrix, NamedArrowMatrix, PathAlgebra, PathAlgebraElement,
        PathAlgebraVector, Quiver, QuiverAlgebra, VertexIndex, right_module_over_path_algebra,
    };
    use serde_json::{Value, json};

    use super::{
        ProjectivePathAlgebraPresentation, ProjectivePathAlgebraPresentationError,
        projective_path_algebra_presentation,
    };
    use crate::{indec_projective_modules, quotient_by_relations};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn small_integer(value: &Integer) -> i64 {
        value.to_string().parse().unwrap()
    }

    fn element_json<F: Field>(
        element: &PathAlgebraElement<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!(
            element
                .external_representation()
                .terms()
                .iter()
                .map(|(path, scalar)| json!({
                    "path": path,
                    "coefficient": coefficient(scalar),
                }))
                .collect::<Vec<_>>()
        )
    }

    fn vector_json<F: Field>(
        vector: &PathAlgebraVector<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!(
            vector
                .components()
                .iter()
                .map(|element| element_json(element, coefficient))
                .collect::<Vec<_>>()
        )
    }

    fn vector_list_json<F: Field>(
        vectors: &[PathAlgebraVector<F>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!(
            vectors
                .iter()
                .map(|vector| vector_json(vector, coefficient))
                .collect::<Vec<_>>()
        )
    }

    fn element_matrix_json<F: Field>(
        matrix: &[Vec<PathAlgebraElement<F>>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!(
            matrix
                .iter()
                .map(|row| row
                    .iter()
                    .map(|element| element_json(element, coefficient))
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>()
        )
    }

    fn presentation_json<F: Field>(
        presentation: &ProjectivePathAlgebraPresentation<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let entered = presentation
            .submodule()
            .map_or_else(Vec::new, |submodule| submodule.generators().to_vec());
        json!({
            "result_length": presentation.result_length(),
            "projective_vertices": presentation
                .projective()
                .vertex_indices()
                .iter()
                .map(|vertex| vec![vertex.get()])
                .collect::<Vec<_>>(),
            "entered_submodule_generators": vector_list_json(&entered, coefficient),
            "f0": vector_list_json(presentation.projective_generators(), coefficient),
            "f1": vector_list_json(presentation.syzygy_generators(), coefficient),
            "matrix": element_matrix_json(presentation.matrix(), coefficient),
            "reduced_basis": vector_list_json(
                presentation.reduced_groebner_basis(),
                coefficient,
            ),
        })
    }

    fn presentation_case<F: Field>(
        field: &F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(
                3,
                vec![
                    ArrowSpec::named(vertex(1), vertex(2), "a"),
                    ArrowSpec::named(vertex(2), vertex(3), "b"),
                ],
            )
            .unwrap(),
        )
        .unwrap();
        let generators = parent.generators();
        let relation = generators[3].multiply(&generators[4]).unwrap();
        let algebra = QuiverAlgebra::Quotient(quotient_by_relations(parent, &[relation]).unwrap());
        let signed_module = right_module_over_path_algebra(
            algebra.clone(),
            vec![2, 1, 0],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![field.one()], vec![field.one()]], field.clone())
                    .unwrap(),
            )],
        )
        .unwrap()
        .unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();

        let relation_projective = projective_path_algebra_presentation(&projectives[0]).unwrap();
        let signed_kernel = projective_path_algebra_presentation(&signed_module).unwrap();
        let sink_projective = projective_path_algebra_presentation(&projectives[2]).unwrap();
        json!({
            "relation_projective": presentation_json(&relation_projective, coefficient),
            "signed_kernel": presentation_json(&signed_kernel, coefficient),
            "sink_projective": presentation_json(&sink_projective, coefficient),
        })
    }

    #[test]
    fn pinned_qpa_projective_path_presentation_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/projpathalgmodule_presentation.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let actual = json!({
            "schema_version": 1,
            "kind": "projpathalgmodule_presentation",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654",
            },
            "binary": presentation_case(&Z2, |entry| {
                small_integer(&Z2.to_integer(entry)).into()
            }),
            "quinary": presentation_case(&field5, |entry| {
                small_integer(&field5.to_integer(entry)).into()
            }),
            "rational": presentation_case(&Q, |entry: &Rational| {
                json!([
                    small_integer(entry.numerator_ref()),
                    small_integer(entry.denominator_ref()),
                ])
            }),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn path_algebra_input_preserves_qpa_source_failure() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let module = right_module_over_path_algebra(algebra, vec![1, 0], &[])
            .unwrap()
            .unwrap();
        assert_eq!(
            projective_path_algebra_presentation(&module).unwrap_err(),
            ProjectivePathAlgebraPresentationError::PathAlgebraHasNoDefiningIdeal,
        );
    }
}
