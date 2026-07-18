//! QPA's Green--Solberg--Zacharia projective resolutions.
//!
//! This module translates the explicit-presentation constructor from
//! `lib/moduleprojres.gi`. The resolution is computed over the original path
//! algebra and retains QPA's parallel quotient-projective, path-projective,
//! differential, `f`, and `fprime` lists.

use std::{cmp::Ordering, fmt, sync::Arc};

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    PathAlgebraElement, PathAlgebraIdeal, PathAlgebraVector, PathEndpoint,
    ProjectivePathModuleError, QuotientError, QuotientPathAlgebra, QuotientPathAlgebraElement,
    RightProjectivePathAlgebraModule, RightProjectivePathAlgebraSubmodule, VertexIndex, first_part,
    right_projective_path_algebra_module, x_set_of_path_algebra_vector,
};

struct RightProjectiveQuotientPathAlgebraModuleData<F: Field> {
    algebra: QuotientPathAlgebra<F>,
    vertices: Box<[VertexIndex]>,
}

/// QPA's quotient-algebra `RightProjectiveModule(A, vertices)`.
///
/// The GSZ constructor only observes the parent quotient and ordered vertex
/// summands of these modules. This identity-preserving descriptor retains that
/// exact source data while the parallel `RProjectives` list stores the full
/// path-algebra vector families used by the computation.
#[derive(Clone)]
pub struct RightProjectiveQuotientPathAlgebraModule<F: Field>(
    Arc<RightProjectiveQuotientPathAlgebraModuleData<F>>,
);

impl<F: Field> RightProjectiveQuotientPathAlgebraModule<F> {
    fn from_vertex_indices(algebra: QuotientPathAlgebra<F>, vertices: Vec<VertexIndex>) -> Self {
        Self(Arc::new(RightProjectiveQuotientPathAlgebraModuleData {
            algebra,
            vertices: vertices.into_boxed_slice(),
        }))
    }

    /// Return the parent quotient path algebra.
    #[must_use]
    pub fn algebra(&self) -> &QuotientPathAlgebra<F> {
        &self.0.algebra
    }

    /// Return the summand vertices in QPA input order.
    #[must_use]
    pub fn vertex_indices(&self) -> &[VertexIndex] {
        &self.0.vertices
    }

    /// Return whether two handles denote the identical QPA projective object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Translate QPA's `LiftPathAlgebraModule` for a quotient vertex projective.
    ///
    /// QPA lifts the summand vertices to the original path algebra and adds,
    /// in vertex-major order, every nonzero product of such a vertex with the
    /// defining ideal's right Gröbner-basis enumerator. Since this constructor
    /// marks its input vertex-projective, no entered quotient-module generators
    /// precede those ideal relations.
    ///
    /// # Errors
    ///
    /// Requires QPA's attached complete two-sided Gröbner basis and propagates
    /// right-basis, path-projective, multiplication, and vector failures.
    pub fn lift_path_algebra_module(
        &self,
    ) -> Result<LiftedProjectivePathAlgebraModule<F>, ProjectiveResolutionError> {
        let ideal = self.algebra().ideal();
        let basis = ideal
            .groebner_basis()
            .ok_or(ProjectiveResolutionError::MissingGroebnerBasis)?;
        if !basis.is_complete() {
            return Err(ProjectiveResolutionError::IncompleteGroebnerBasis);
        }
        let parent = self.algebra().original_path_algebra();
        let vertices = self
            .vertex_indices()
            .iter()
            .copied()
            .map(|vertex| embedded_vertex(parent, vertex))
            .collect::<Result<Vec<_>, _>>()?;
        let lifted_parent = right_projective_path_algebra_module(parent, &vertices)?;
        let right_basis = ideal.right_groebner_basis()?;
        let mut generators = Vec::new();
        for (position, vertex) in vertices.iter().enumerate() {
            for relation in right_basis.enumerator() {
                let component = vertex.multiply(relation)?;
                if component.is_zero() {
                    continue;
                }
                let mut components = vec![parent.zero(); vertices.len()];
                components[position] = component;
                generators.push(lifted_parent.vector(components)?);
            }
        }
        if generators.is_empty() {
            Ok(LiftedProjectivePathAlgebraModule::Projective(lifted_parent))
        } else {
            let submodule = lifted_parent.submodule(generators)?;
            Ok(LiftedProjectivePathAlgebraModule::Submodule(submodule))
        }
    }
}

impl<F: Field> fmt::Debug for RightProjectiveQuotientPathAlgebraModule<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("RightProjectiveQuotientPathAlgebraModule")
            .field("algebra", self.algebra())
            .field("vertices", &self.vertex_indices())
            .finish()
    }
}

impl<F: Field> PartialEq for RightProjectiveQuotientPathAlgebraModule<F> {
    fn eq(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}

impl<F: Field> Eq for RightProjectiveQuotientPathAlgebraModule<F> {}

/// QPA's path-algebra lift of a quotient vertex projective.
#[derive(Clone, Debug)]
pub enum LiftedProjectivePathAlgebraModule<F: Field> {
    /// No right ideal relation starts at a selected vertex, so QPA returns the
    /// whole lifted projective.
    Projective(RightProjectivePathAlgebraModule<F>),
    /// QPA returns the submodule generated by the lifted ideal relations.
    Submodule(RightProjectivePathAlgebraSubmodule<F>),
}

impl<F: Field> LiftedProjectivePathAlgebraModule<F> {
    /// Return the lifted whole projective/vector family.
    #[must_use]
    pub fn parent_module(&self) -> &RightProjectivePathAlgebraModule<F> {
        match self {
            Self::Projective(module) => module,
            Self::Submodule(module) => module.parent(),
        }
    }

    /// Return QPA's submodule result, when ideal generators were nonempty.
    #[must_use]
    pub const fn submodule(&self) -> Option<&RightProjectivePathAlgebraSubmodule<F>> {
        match self {
            Self::Projective(_) => None,
            Self::Submodule(module) => Some(module),
        }
    }

    /// Return whether QPA selected its submodule branch.
    #[must_use]
    pub const fn is_submodule(&self) -> bool {
        matches!(self, Self::Submodule(_))
    }
}

/// Construct QPA's quotient-path-algebra vertex projective.
///
/// # Errors
///
/// Every entered element must be a coefficient-one embedded vertex in the
/// identical quotient algebra.
pub fn right_projective_quotient_path_algebra_module<F: Field>(
    algebra: &QuotientPathAlgebra<F>,
    vertices: &[QuotientPathAlgebraElement<F>],
) -> Result<RightProjectiveQuotientPathAlgebraModule<F>, ProjectiveResolutionError> {
    if vertices.is_empty() {
        return Err(ProjectiveResolutionError::EmptyQuotientProjectiveVertices);
    }
    let mut indices = Vec::with_capacity(vertices.len());
    for (position, vertex) in vertices.iter().enumerate() {
        if !vertex.algebra().is_identical(algebra) {
            return Err(ProjectiveResolutionError::ForeignQuotientVertex {
                position: position + 1,
            });
        }
        let [term] = vertex.representative().terms() else {
            return Err(ProjectiveResolutionError::NotEmbeddedQuotientVertex {
                position: position + 1,
            });
        };
        let PathEndpoint::Vertex(source) = term.path().source() else {
            return Err(ProjectiveResolutionError::NotEmbeddedQuotientVertex {
                position: position + 1,
            });
        };
        if !term.path().is_vertex()
            || !algebra
                .original_path_algebra()
                .field()
                .is_one(term.coefficient())
        {
            return Err(ProjectiveResolutionError::NotEmbeddedQuotientVertex {
                position: position + 1,
            });
        }
        indices.push(source.index());
    }
    Ok(RightProjectiveQuotientPathAlgebraModule::from_vertex_indices(algebra.clone(), indices))
}

/// QPA's `ProjectiveResolutionFpPathAlgebraModule` attribute object.
#[derive(Clone, Debug)]
pub struct ProjectiveResolutionFpPathAlgebraModule<F: Field> {
    parent_algebra: QuotientPathAlgebra<F>,
    ring_ideal: PathAlgebraIdeal<F>,
    projectives_attribute_is_bound: bool,
    projectives: Vec<RightProjectiveQuotientPathAlgebraModule<F>>,
    r_projectives: Vec<RightProjectivePathAlgebraModule<F>>,
    maps: Vec<Vec<Vec<PathAlgebraElement<F>>>>,
    f_lists: Vec<Vec<PathAlgebraVector<F>>>,
    f_prime_lists: Vec<Vec<PathAlgebraVector<F>>>,
    finite_at_projective: Option<usize>,
}

impl<F: Field> ProjectiveResolutionFpPathAlgebraModule<F> {
    /// Return QPA's `ParentAlgebra` attribute.
    #[must_use]
    pub const fn parent_algebra(&self) -> &QuotientPathAlgebra<F> {
        &self.parent_algebra
    }

    /// Return QPA's `RingIdeal` attribute.
    #[must_use]
    pub const fn ring_ideal(&self) -> &PathAlgebraIdeal<F> {
        &self.ring_ideal
    }

    /// Return QPA's quotient-algebra `Projectives` list.
    #[must_use]
    pub fn projectives(&self) -> &[RightProjectiveQuotientPathAlgebraModule<F>] {
        &self.projectives
    }

    /// Return QPA's `HasProjectives` value.
    ///
    /// `ProjectiveResolutionOfPathAlgebraModule` deliberately leaves this
    /// attribute unbound, whereas the explicit-presentation constructor binds
    /// its parallel quotient-projective list.
    #[must_use]
    pub const fn projectives_attribute_is_bound(&self) -> bool {
        self.projectives_attribute_is_bound
    }

    /// Return QPA's original-path-algebra `RProjectives` list.
    #[must_use]
    pub fn r_projectives(&self) -> &[RightProjectivePathAlgebraModule<F>] {
        &self.r_projectives
    }

    /// Return QPA's ordered differential matrices over the original path algebra.
    #[must_use]
    pub fn maps(&self) -> &[Vec<Vec<PathAlgebraElement<F>>>] {
        &self.maps
    }

    /// Return the first half of QPA's `ProjectivesFList` pair.
    #[must_use]
    pub fn f_lists(&self) -> &[Vec<PathAlgebraVector<F>>] {
        &self.f_lists
    }

    /// Return the second half of QPA's `ProjectivesFList` pair.
    #[must_use]
    pub fn f_prime_lists(&self) -> &[Vec<PathAlgebraVector<F>>] {
        &self.f_prime_lists
    }

    /// Return the degree at which QPA prints `finite at projective`, if reached.
    #[must_use]
    pub const fn finite_at_projective(&self) -> Option<usize> {
        self.finite_at_projective
    }
}

/// Translate QPA's four-argument explicit-presentation GSZ constructor.
///
/// `presentation_map` is entered over `algebra`; QPA lifts its external
/// representations into the original path algebra before computing module
/// Gröbner bases. As in the source, the degree-two data is constructed even
/// when `termint` is one.
///
/// # Errors
///
/// Enforces QPA's positive-term, attached-basis, nonempty-table, uniformity,
/// vertex-consistency, and nonzero-column preconditions. Exact quotient,
/// path-algebra, projective-vector, and GSZ failures are propagated.
#[allow(clippy::too_many_lines)]
pub fn projective_resolution_fp_path_algebra_module<F: Field>(
    algebra: &QuotientPathAlgebra<F>,
    ideal: &PathAlgebraIdeal<F>,
    presentation_map: &[Vec<QuotientPathAlgebraElement<F>>],
    termint: usize,
) -> Result<ProjectiveResolutionFpPathAlgebraModule<F>, ProjectiveResolutionError> {
    if termint == 0 {
        return Err(ProjectiveResolutionError::NonpositiveTerm);
    }
    if ideal.groebner_basis().is_none() {
        return Err(ProjectiveResolutionError::MissingGroebnerBasis);
    }
    if !ideal
        .algebra()
        .is_identical(algebra.original_path_algebra())
    {
        return Err(ProjectiveResolutionError::ForeignResolutionIdeal);
    }
    let Some(first_row) = presentation_map.first() else {
        return Err(ProjectiveResolutionError::EmptyPresentationMap);
    };
    let column_count = first_row.len();
    for (row_position, row) in presentation_map.iter().enumerate() {
        if row.len() != column_count {
            return Err(ProjectiveResolutionError::RaggedPresentationMap {
                row: row_position + 1,
                actual: row.len(),
                expected: column_count,
            });
        }
        for (column_position, entry) in row.iter().enumerate() {
            if !entry.algebra().is_identical(algebra) {
                return Err(ProjectiveResolutionError::ForeignPresentationEntry {
                    row: row_position + 1,
                    column: column_position + 1,
                });
            }
            if !entry.is_left_uniform() {
                return Err(ProjectiveResolutionError::PresentationEntryNotLeftUniform {
                    row: row_position + 1,
                    column: column_position + 1,
                });
            }
        }
    }

    let mut source_vertices = first_row
        .iter()
        .map(quotient_element_source_vertex)
        .collect::<Result<Vec<_>, _>>()?;
    for (row_position, row) in presentation_map.iter().enumerate() {
        for (column_position, entry) in row.iter().enumerate() {
            let source = quotient_element_source_vertex(entry)?;
            match (source_vertices[column_position], source) {
                (None, Some(vertex)) => source_vertices[column_position] = Some(vertex),
                (Some(expected), Some(actual)) if expected != actual => {
                    return Err(ProjectiveResolutionError::MismatchedStartingVertices {
                        row: row_position + 1,
                        column: column_position + 1,
                        expected,
                        actual,
                    });
                }
                _ => {}
            }
        }
    }
    let source_vertices = source_vertices
        .into_iter()
        .enumerate()
        .map(|(column, vertex)| {
            vertex.ok_or(ProjectiveResolutionError::ZeroPresentationColumn { column: column + 1 })
        })
        .collect::<Result<Vec<_>, _>>()?;
    if source_vertices.is_empty() {
        return Err(ProjectiveResolutionError::EmptyResolutionProjective { degree: 0 });
    }

    let parent = ideal.algebra();
    let first_r_projective = path_projective_from_indices(parent, &source_vertices)?;
    let first_projective = RightProjectiveQuotientPathAlgebraModule::from_vertex_indices(
        algebra.clone(),
        source_vertices,
    );
    let f_zero = first_r_projective.generators()?;
    let lifted_presentation = presentation_map
        .iter()
        .map(|row| {
            row.iter()
                .map(|entry| entry.representative().clone())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let presentation_vectors = lifted_presentation
        .iter()
        .cloned()
        .map(|row| first_r_projective.vector(row))
        .collect::<Result<Vec<_>, _>>()?;

    let right_relations = ideal.right_groebner_basis()?.relations().to_vec();
    let mut ideal_vectors = Vec::new();
    for (position, vertex) in first_r_projective
        .vertex_indices()
        .iter()
        .copied()
        .enumerate()
    {
        let vertex_element = embedded_vertex(parent, vertex)?;
        for relation in &right_relations {
            let product = vertex_element.multiply(relation)?;
            if product.is_zero() {
                continue;
            }
            let mut components = vec![parent.zero(); first_r_projective.vertex_indices().len()];
            components[position] = product;
            ideal_vectors.push(first_r_projective.vector(components)?);
        }
    }
    if ideal_vectors.is_empty() {
        return Err(ProjectiveResolutionError::EmptyProjectiveIdealGenerators);
    }

    let mut combined = presentation_vectors;
    combined.extend(ideal_vectors.iter().cloned());
    let presentation_submodule = first_r_projective.submodule(combined)?;
    let ideal_submodule = first_r_projective.submodule(ideal_vectors)?;
    let presentation_basis = presentation_submodule.right_groebner_basis()?;
    let ideal_basis = ideal_submodule.right_groebner_basis()?;

    let mut f_one = Vec::new();
    let mut f_prime_one = Vec::new();
    for vector in presentation_basis.basis_vectors() {
        if ideal_basis.completely_reduce(vector)?.is_zero() {
            f_prime_one.push(vector.clone());
        } else {
            f_one.push(vector.clone());
        }
    }
    sort_vectors_descending(&mut f_one)?;
    sort_vectors_descending(&mut f_prime_one)?;

    let f_one_vertices = right_uniform_target_vertices(&f_one, 1)?;
    if f_one_vertices.is_empty() {
        return Err(ProjectiveResolutionError::EmptyResolutionProjective { degree: 1 });
    }
    let second_r_projective = path_projective_from_indices(parent, &f_one_vertices)?;
    let second_projective = RightProjectiveQuotientPathAlgebraModule::from_vertex_indices(
        algebra.clone(),
        f_one_vertices,
    );

    let mut projectives = vec![first_projective, second_projective];
    let mut r_projectives = vec![first_r_projective.clone(), second_r_projective];
    let mut maps = vec![lifted_presentation];
    let mut f_lists = vec![f_zero, f_one];
    let mut f_prime_lists = vec![Vec::new(), f_prime_one];

    let degree_two = next_resolution_step(
        &f_lists[0],
        &f_prime_lists[0],
        &f_lists[1],
        &f_prime_lists[1],
        ideal,
        &first_r_projective,
        false,
    )?;
    maps.push(degree_two.map);
    f_lists.push(degree_two.f);
    f_prime_lists.push(degree_two.f_prime);
    let degree_two_vertices = right_uniform_target_vertices(&f_lists[2], 2)?;
    if degree_two_vertices.is_empty() {
        return Err(ProjectiveResolutionError::EmptyResolutionProjective { degree: 2 });
    }
    r_projectives.push(path_projective_from_indices(parent, &degree_two_vertices)?);
    projectives.push(
        RightProjectiveQuotientPathAlgebraModule::from_vertex_indices(
            algebra.clone(),
            degree_two_vertices,
        ),
    );

    let mut finite_at_projective = None;
    for degree in 3..=termint {
        let step = next_resolution_step(
            &f_lists[degree - 2],
            &f_prime_lists[degree - 2],
            &f_lists[degree - 1],
            &f_prime_lists[degree - 1],
            ideal,
            &first_r_projective,
            true,
        )?;
        let map_is_empty = step.map.is_empty();
        maps.push(step.map);
        f_lists.push(step.f);
        f_prime_lists.push(step.f_prime);
        if map_is_empty {
            finite_at_projective = Some(degree);
            break;
        }
        let vertices = right_uniform_target_vertices(&f_lists[degree], degree)?;
        r_projectives.push(path_projective_from_indices(parent, &vertices)?);
        projectives.push(
            RightProjectiveQuotientPathAlgebraModule::from_vertex_indices(
                algebra.clone(),
                vertices,
            ),
        );
    }

    Ok(ProjectiveResolutionFpPathAlgebraModule {
        parent_algebra: algebra.clone(),
        ring_ideal: ideal.clone(),
        projectives_attribute_is_bound: true,
        projectives,
        r_projectives,
        maps,
        f_lists,
        f_prime_lists,
        finite_at_projective,
    })
}

/// Complete QPA's `ProjectiveResolutionOfPathAlgebraModule` after its
/// GBNP-backed projective presentation has been constructed.
///
/// This cross-crate hook starts with QPA's `ProjPres[1]`, `ProjPres[3]`,
/// `ProjPres[4]`, and `ProjPres[5]` values. It splits the first module
/// Groebner basis into `f^1` and `f^1prime` against `P^0 I`, then applies the
/// same GSZ X-set recurrence as the explicit-presentation constructor.
/// Unlike that constructor, QPA leaves the quotient `Projectives` attribute
/// unbound and begins its repeated loop at degree two.
///
/// # Errors
///
/// Enforces the positive term, identical parent algebra, attached defining
/// basis, projective-vector family, nonempty relation-submodule, and GSZ
/// uniformity boundaries of the installed source method.
#[doc(hidden)]
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn projective_resolution_of_path_algebra_module_from_presentation<F: Field>(
    algebra: &QuotientPathAlgebra<F>,
    ideal: &PathAlgebraIdeal<F>,
    first_projective: &RightProjectivePathAlgebraModule<F>,
    f_zero: &[PathAlgebraVector<F>],
    presentation_syzygies: &[PathAlgebraVector<F>],
    presentation_map: &[Vec<PathAlgebraElement<F>>],
    termint: usize,
) -> Result<ProjectiveResolutionFpPathAlgebraModule<F>, ProjectiveResolutionError> {
    if termint == 0 {
        return Err(ProjectiveResolutionError::NonpositiveTerm);
    }
    if ideal.groebner_basis().is_none() {
        return Err(ProjectiveResolutionError::MissingGroebnerBasis);
    }
    if !ideal
        .algebra()
        .is_identical(algebra.original_path_algebra())
    {
        return Err(ProjectiveResolutionError::ForeignResolutionIdeal);
    }
    if !first_projective.algebra().is_identical(ideal.algebra()) {
        return Err(ProjectiveResolutionError::ForeignInitialPathProjective);
    }
    if f_zero
        .iter()
        .chain(presentation_syzygies)
        .any(|vector| !vector.module().is_identical(first_projective))
    {
        return Err(ProjectiveResolutionError::ForeignInitialResolutionVector);
    }
    if presentation_map
        .iter()
        .flatten()
        .any(|entry| !entry.algebra().is_identical(ideal.algebra()))
    {
        return Err(ProjectiveResolutionError::ForeignInitialResolutionMap);
    }

    let right_relations = ideal.right_groebner_basis()?.relations().to_vec();
    let mut projective_ideal_generators = Vec::new();
    for generator in f_zero {
        for relation in &right_relations {
            let product = generator.act(relation)?;
            if !product.is_zero() {
                projective_ideal_generators.push(product);
            }
        }
    }
    if presentation_syzygies.is_empty() {
        return Err(ProjectiveResolutionError::EmptyPresentationSyzygies);
    }
    if projective_ideal_generators.is_empty() {
        return Err(ProjectiveResolutionError::EmptyProjectiveIdealGenerators);
    }

    let presentation_submodule = first_projective.submodule(presentation_syzygies.to_vec())?;
    let ideal_submodule = first_projective.submodule(projective_ideal_generators)?;
    let presentation_basis = presentation_submodule.right_groebner_basis()?;
    let ideal_basis = ideal_submodule.right_groebner_basis()?;

    let mut f_one = Vec::new();
    let mut f_prime_one = Vec::new();
    for vector in presentation_basis.basis_vectors() {
        if ideal_basis.completely_reduce(vector)?.is_zero() {
            f_prime_one.push(vector.clone());
        } else {
            f_one.push(vector.clone());
        }
    }
    sort_vectors_descending(&mut f_one)?;
    sort_vectors_descending(&mut f_prime_one)?;

    let f_one_vertices = right_uniform_target_vertices(&f_one, 1)?;
    if f_one_vertices.is_empty() {
        return Err(ProjectiveResolutionError::EmptyResolutionProjective { degree: 1 });
    }
    let second_r_projective = path_projective_from_indices(ideal.algebra(), &f_one_vertices)?;
    let mut r_projectives = vec![first_projective.clone(), second_r_projective];
    let mut maps = vec![presentation_map.to_vec()];
    let mut f_lists = vec![f_zero.to_vec(), f_one];
    let mut f_prime_lists = vec![Vec::new(), f_prime_one];
    let mut finite_at_projective = None;

    for degree in 2..=termint {
        let step = next_resolution_step(
            &f_lists[degree - 2],
            &f_prime_lists[degree - 2],
            &f_lists[degree - 1],
            &f_prime_lists[degree - 1],
            ideal,
            first_projective,
            true,
        )?;
        let map_is_empty = step.map.is_empty();
        maps.push(step.map);
        f_lists.push(step.f);
        f_prime_lists.push(step.f_prime);
        if map_is_empty {
            finite_at_projective = Some(degree);
            break;
        }
        let vertices = right_uniform_target_vertices(&f_lists[degree], degree)?;
        if vertices.is_empty() {
            return Err(ProjectiveResolutionError::EmptyResolutionProjective { degree });
        }
        r_projectives.push(path_projective_from_indices(ideal.algebra(), &vertices)?);
    }

    Ok(ProjectiveResolutionFpPathAlgebraModule {
        parent_algebra: algebra.clone(),
        ring_ideal: ideal.clone(),
        projectives_attribute_is_bound: false,
        projectives: Vec::new(),
        r_projectives,
        maps,
        f_lists,
        f_prime_lists,
        finite_at_projective,
    })
}

fn quotient_element_source_vertex<F: Field>(
    element: &QuotientPathAlgebraElement<F>,
) -> Result<Option<VertexIndex>, ProjectiveResolutionError> {
    if element.is_zero() {
        return Ok(None);
    }
    match element.leading_monomial().source() {
        PathEndpoint::Vertex(vertex) => Ok(Some(vertex.index())),
        PathEndpoint::Zero(_) => Err(ProjectiveResolutionError::ZeroLeadingPath),
    }
}

fn embedded_vertex<F: Field>(
    algebra: &crate::PathAlgebra<F>,
    vertex: VertexIndex,
) -> Result<PathAlgebraElement<F>, ProjectiveResolutionError> {
    let path = algebra
        .quiver()
        .vertex(vertex)
        .ok_or(ProjectiveResolutionError::VertexOutOfRange { vertex })?
        .path();
    Ok(algebra.element(path)?)
}

fn path_projective_from_indices<F: Field>(
    algebra: &crate::PathAlgebra<F>,
    vertices: &[VertexIndex],
) -> Result<RightProjectivePathAlgebraModule<F>, ProjectiveResolutionError> {
    let elements = vertices
        .iter()
        .copied()
        .map(|vertex| embedded_vertex(algebra, vertex))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(right_projective_path_algebra_module(algebra, &elements)?)
}

fn sort_vectors_descending<F: Field>(
    vectors: &mut Vec<PathAlgebraVector<F>>,
) -> Result<(), ProjectiveResolutionError> {
    let entered = std::mem::take(vectors);
    let mut ascending = Vec::with_capacity(entered.len());
    for vector in entered {
        let mut insertion = ascending.len();
        for (position, existing) in ascending.iter().enumerate() {
            if vector.compare(existing)? == Ordering::Less {
                insertion = position;
                break;
            }
        }
        ascending.insert(insertion, vector);
    }
    ascending.reverse();
    *vectors = ascending;
    Ok(())
}

fn right_uniform_target_vertices<F: Field>(
    vectors: &[PathAlgebraVector<F>],
    degree: usize,
) -> Result<Vec<VertexIndex>, ProjectiveResolutionError> {
    vectors
        .iter()
        .enumerate()
        .map(|(position, vector)| {
            if !vector.is_right_uniform() {
                return Err(ProjectiveResolutionError::ResolutionVectorNotRightUniform {
                    degree,
                    position: position + 1,
                });
            }
            let target = vector.target_vertex().ok_or(
                ProjectiveResolutionError::ResolutionVectorNotRightUniform {
                    degree,
                    position: position + 1,
                },
            )?;
            match target.source() {
                PathEndpoint::Vertex(vertex) => Ok(vertex.index()),
                PathEndpoint::Zero(_) => Err(ProjectiveResolutionError::ZeroResolutionTarget {
                    degree,
                    position: position + 1,
                }),
            }
        })
        .collect()
}

struct ResolutionStep<F: Field> {
    map: Vec<Vec<PathAlgebraElement<F>>>,
    f: Vec<PathAlgebraVector<F>>,
    f_prime: Vec<PathAlgebraVector<F>>,
}

#[allow(clippy::too_many_arguments)]
fn next_resolution_step<F: Field>(
    previous_f: &[PathAlgebraVector<F>],
    previous_f_prime: &[PathAlgebraVector<F>],
    current_f: &[PathAlgebraVector<F>],
    current_f_prime: &[PathAlgebraVector<F>],
    ideal: &PathAlgebraIdeal<F>,
    first_projective: &RightProjectivePathAlgebraModule<F>,
    skip_zero_overlap: bool,
) -> Result<ResolutionStep<F>, ProjectiveResolutionError> {
    let mut map = Vec::new();
    let mut next_f = Vec::new();
    let mut next_f_prime = Vec::new();
    for current in current_f {
        let x_set = x_set_of_path_algebra_vector(previous_f, previous_f_prime, ideal, current)?;
        for overlap in x_set.o_set {
            let generator = previous_f.get(overlap.generator_position - 1).ok_or(
                ProjectiveResolutionError::InvalidXSetGeneratorPosition {
                    position: overlap.generator_position,
                    length: previous_f.len(),
                },
            )?;
            let factor = overlap.q_prime.multiply(&overlap.relation)?;
            let product = generator.act(&factor)?;
            if skip_zero_overlap && product.is_zero() {
                continue;
            }
            let coefficients = first_part(current_f, current_f_prime, &product)?;
            let mut next = first_projective.zero();
            for (basis_vector, coefficient) in current_f.iter().zip(&coefficients) {
                next = next.add(&basis_vector.act(coefficient)?)?;
            }
            map.push(coefficients);
            next_f.push(next);
        }
        for nonoverlap in x_set.n_set {
            let extended = current.act(&nonoverlap.z)?;
            next_f_prime.push(extended.act(&nonoverlap.relation)?);
        }
    }
    Ok(ResolutionStep {
        map,
        f: next_f,
        f_prime: next_f_prime,
    })
}

/// Errors from QPA's explicit GSZ projective-resolution constructor.
#[derive(Clone, Debug, Error, Eq, PartialEq)]
pub enum ProjectiveResolutionError {
    /// QPA dispatches the constructor only on a positive integer.
    #[error("the projective-resolution term must be positive")]
    NonpositiveTerm,
    /// The ideal fails QPA's `HasGroebnerBasisOfIdeal` method filter.
    #[error("the projective-resolution ideal has no attached Groebner basis")]
    MissingGroebnerBasis,
    /// QPA's lift method rejects an attached basis not marked complete.
    #[error("the quotient defining ideal's Groebner basis is not complete")]
    IncompleteGroebnerBasis,
    /// The quotient and ideal do not share the identical original path algebra.
    #[error("the projective-resolution ideal belongs to another path algebra")]
    ForeignResolutionIdeal,
    /// QPA explicitly rejects an empty outer presentation list.
    #[error("the projective presentation map must be nonempty")]
    EmptyPresentationMap,
    /// Rust makes the `IsRingElementTable` rectangularity filter explicit.
    #[error("presentation row {row} has length {actual}, expected {expected}")]
    RaggedPresentationMap {
        /// One-based row position.
        row: usize,
        /// Entered row length.
        actual: usize,
        /// Required row length.
        expected: usize,
    },
    /// An entry belongs to another quotient family.
    #[error("presentation entry ({row}, {column}) belongs to another quotient algebra")]
    ForeignPresentationEntry {
        /// One-based row position.
        row: usize,
        /// One-based column position.
        column: usize,
    },
    /// QPA explicitly requires every entry to be left uniform.
    #[error("presentation entry ({row}, {column}) is not left uniform")]
    PresentationEntryNotLeftUniform {
        /// One-based row position.
        row: usize,
        /// One-based column position.
        column: usize,
    },
    /// Nonzero entries in one column have inconsistent source vertices.
    #[error(
        "presentation entry ({row}, {column}) starts at vertex {actual:?}, expected {expected:?}"
    )]
    MismatchedStartingVertices {
        /// One-based row position.
        row: usize,
        /// One-based column position.
        column: usize,
        /// Vertex selected by an earlier row.
        expected: VertexIndex,
        /// Vertex in this row.
        actual: VertexIndex,
    },
    /// QPA explicitly rejects an all-zero column.
    #[error("presentation column {column} contains only zeroes")]
    ZeroPresentationColumn {
        /// One-based column position.
        column: usize,
    },
    /// A resolution vertex index was not in the parent quiver.
    #[error("resolution vertex {vertex:?} is outside the parent quiver")]
    VertexOutOfRange {
        /// Offending vertex index.
        vertex: VertexIndex,
    },
    /// The ideal contributes no generators to the initial projective.
    #[error("the initial projective times the right ideal basis is empty")]
    EmptyProjectiveIdealGenerators,
    /// QPA's module-resolution wrapper cannot form `N` from an empty `ProjPres[4]`.
    #[error("the projective path presentation exposes no first syzygy generators")]
    EmptyPresentationSyzygies,
    /// The presentation projective must use the quotient's original path algebra.
    #[error("the initial path projective belongs to another path algebra")]
    ForeignInitialPathProjective,
    /// Every initial `f` vector must belong to the presentation projective.
    #[error("an initial resolution vector belongs to another projective family")]
    ForeignInitialResolutionVector,
    /// Every first-map entry must lie in the quotient's original path algebra.
    #[error("the initial resolution map contains an element from another path algebra")]
    ForeignInitialResolutionMap,
    /// QPA's quotient-projective method does not dispatch on an empty list.
    #[error("a quotient vertex projective requires a nonempty vertex list")]
    EmptyQuotientProjectiveVertices,
    /// QPA's constructor reached an empty vertex list before its finite-stop branch.
    #[error("the degree-{degree} resolution projective has no vertex summands")]
    EmptyResolutionProjective {
        /// Projective degree whose vertex list is empty.
        degree: usize,
    },
    /// QPA's source or target helper unexpectedly returned the zero path.
    #[error("a nonzero presentation entry has a zero leading path")]
    ZeroLeadingPath,
    /// A computed `f` vector violates the GSZ right-uniform invariant.
    #[error("f^{degree} vector {position} is not right uniform")]
    ResolutionVectorNotRightUniform {
        /// GSZ superscript.
        degree: usize,
        /// One-based vector position.
        position: usize,
    },
    /// A computed zero `f` vector has no projective target vertex.
    #[error("f^{degree} vector {position} has zero target vertex")]
    ZeroResolutionTarget {
        /// GSZ superscript.
        degree: usize,
        /// One-based vector position.
        position: usize,
    },
    /// An X-set returned an invalid one-based generator position.
    #[error("X-set generator position {position} is outside a list of length {length}")]
    InvalidXSetGeneratorPosition {
        /// One-based position from the X-set.
        position: usize,
        /// Available generator count.
        length: usize,
    },
    /// A quotient-projective vertex belongs to another quotient.
    #[error("quotient-projective vertex {position} belongs to another algebra")]
    ForeignQuotientVertex {
        /// One-based list position.
        position: usize,
    },
    /// A quotient-projective entry is not a coefficient-one vertex.
    #[error("quotient-projective entry {position} is not an embedded vertex")]
    NotEmbeddedQuotientVertex {
        /// One-based list position.
        position: usize,
    },
    /// Original path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] crate::PathAlgebraError),
    /// Quotient or ideal arithmetic failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Projective path-vector arithmetic failed.
    #[error(transparent)]
    Projective(#[from] ProjectivePathModuleError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field,
        finite_field::{FiniteFieldCore, Z2, Zp},
        rational::{Q, Rational},
    };
    use serde_json::{Value, json};

    use crate::{
        ArrowSpec, PathAlgebra, PathAlgebraElement, PathAlgebraIdeal, PathAlgebraVector, Quiver,
        QuotientPathAlgebra, VertexIndex,
    };

    use super::{
        ProjectiveResolutionError, ProjectiveResolutionFpPathAlgebraModule,
        projective_resolution_fp_path_algebra_module,
        right_projective_quotient_path_algebra_module,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn element_json<F: Field>(
        element: &PathAlgebraElement<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().external_representation(),
                        "coefficient": coefficient(term.coefficient()),
                    })
                })
                .collect(),
        )
    }

    fn element_list_json<F: Field>(
        elements: &[PathAlgebraElement<F>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            elements
                .iter()
                .map(|element| element_json(element, coefficient))
                .collect(),
        )
    }

    fn vector_json<F: Field>(
        vector: &PathAlgebraVector<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        element_list_json(vector.components(), coefficient)
    }

    fn vector_lists_json<F: Field>(
        lists: &[Vec<PathAlgebraVector<F>>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            lists
                .iter()
                .map(|vectors| {
                    Value::Array(
                        vectors
                            .iter()
                            .map(|vector| vector_json(vector, coefficient))
                            .collect(),
                    )
                })
                .collect(),
        )
    }

    fn vector_list_json<F: Field>(
        vectors: &[PathAlgebraVector<F>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            vectors
                .iter()
                .map(|vector| vector_json(vector, coefficient))
                .collect(),
        )
    }

    fn maps_json<F: Field>(
        maps: &[Vec<Vec<PathAlgebraElement<F>>>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            maps.iter()
                .map(|matrix| {
                    Value::Array(
                        matrix
                            .iter()
                            .map(|row| element_list_json(row, coefficient))
                            .collect(),
                    )
                })
                .collect(),
        )
    }

    fn vertex_path_json<F: Field>(algebra: &PathAlgebra<F>, vertices: &[VertexIndex]) -> Value {
        Value::Array(
            vertices
                .iter()
                .map(|vertex| {
                    json!(
                        algebra
                            .quiver()
                            .vertex(*vertex)
                            .unwrap()
                            .path()
                            .external_representation()
                    )
                })
                .collect(),
        )
    }

    fn resolution_json<F: Field>(
        resolution: &ProjectiveResolutionFpPathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let parent = resolution.ring_ideal().algebra();
        json!({
            "maps": maps_json(resolution.maps(), coefficient),
            "f": vector_lists_json(resolution.f_lists(), coefficient),
            "f_prime": vector_lists_json(resolution.f_prime_lists(), coefficient),
            "r_projective_vertices": resolution.r_projectives().iter().map(|projective| {
                vertex_path_json(parent, projective.vertex_indices())
            }).collect::<Vec<_>>(),
            "projective_vertices": resolution.projectives().iter().map(|projective| {
                vertex_path_json(parent, projective.vertex_indices())
            }).collect::<Vec<_>>(),
        })
    }

    fn resolution_case<F: Field>(
        field: &F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let dual_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let dual_generators = dual_parent.generators();
        let dual_relation = dual_generators[1].multiply(&dual_generators[1]).unwrap();
        let dual_ideal = PathAlgebraIdeal::new(dual_parent, vec![dual_relation.clone()])
            .unwrap()
            .with_complete_groebner_basis(std::slice::from_ref(&dual_relation))
            .unwrap();
        let dual_algebra = QuotientPathAlgebra::by_ideal(&dual_ideal).unwrap();
        let dual_x = dual_algebra.generators().unwrap()[1].clone();
        let term_four = projective_resolution_fp_path_algebra_module(
            &dual_algebra,
            &dual_ideal,
            &[vec![dual_x.clone()]],
            4,
        )
        .unwrap();
        let term_one = projective_resolution_fp_path_algebra_module(
            &dual_algebra,
            &dual_ideal,
            &[vec![dual_x]],
            1,
        )
        .unwrap();

        let overlap_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(
                1,
                vec![
                    ArrowSpec::named(vertex(1), vertex(1), "x"),
                    ArrowSpec::named(vertex(1), vertex(1), "y"),
                ],
            )
            .unwrap(),
        )
        .unwrap();
        let overlap_generators = overlap_parent.generators();
        let overlap_relations = vec![
            overlap_generators[1]
                .multiply(&overlap_generators[1])
                .unwrap(),
            overlap_generators[1]
                .multiply(&overlap_generators[2])
                .unwrap(),
            overlap_generators[2]
                .multiply(&overlap_generators[2])
                .unwrap(),
        ];
        let overlap_ideal = PathAlgebraIdeal::new(overlap_parent, overlap_relations.clone())
            .unwrap()
            .with_complete_groebner_basis(&overlap_relations)
            .unwrap();
        let overlap_algebra = QuotientPathAlgebra::by_ideal(&overlap_ideal).unwrap();
        let overlap_generators = overlap_algebra.generators().unwrap();
        let overlap = projective_resolution_fp_path_algebra_module(
            &overlap_algebra,
            &overlap_ideal,
            &[
                vec![overlap_generators[1].clone()],
                vec![overlap_generators[2].clone()],
            ],
            3,
        )
        .unwrap();

        json!({
            "term_four": resolution_json(&term_four, coefficient),
            "term_one": resolution_json(&term_one, coefficient),
            "overlap": resolution_json(&overlap, coefficient),
        })
    }

    fn lift_case<F: Field>(field: &F, coefficient: impl Fn(&F::Element) -> Value + Copy) -> Value {
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
        let ideal = PathAlgebraIdeal::new(parent, vec![relation.clone()])
            .unwrap()
            .with_complete_groebner_basis(std::slice::from_ref(&relation))
            .unwrap();
        let algebra = QuotientPathAlgebra::by_ideal(&ideal).unwrap();
        let quotient_generators = algebra.generators().unwrap();
        let source = right_projective_quotient_path_algebra_module(
            &algebra,
            &[
                quotient_generators[0].clone(),
                quotient_generators[1].clone(),
            ],
        )
        .unwrap();
        let source_lift = source.lift_path_algebra_module().unwrap();
        let sink = right_projective_quotient_path_algebra_module(
            &algebra,
            std::slice::from_ref(&quotient_generators[2]),
        )
        .unwrap();
        let sink_lift = sink.lift_path_algebra_module().unwrap();
        let lift_json = |lift: &super::LiftedProjectivePathAlgebraModule<F>| {
            json!({
                "is_submodule": lift.is_submodule(),
                "parent_vertices": lift
                    .parent_module()
                    .vertex_indices()
                    .iter()
                    .map(|vertex| vertex.get())
                    .collect::<Vec<_>>(),
                "entered_generators": lift.submodule().map_or_else(
                    || Value::Array(Vec::new()),
                    |submodule| vector_list_json(submodule.generators(), coefficient),
                ),
            })
        };
        json!({
            "source_vertices": source.vertex_indices().iter().map(|vertex| vertex.get()).collect::<Vec<_>>(),
            "source_lift": lift_json(&source_lift),
            "sink_vertices": sink.vertex_indices().iter().map(|vertex| vertex.get()).collect::<Vec<_>>(),
            "sink_lift": lift_json(&sink_lift),
        })
    }

    #[test]
    fn pinned_qpa_moduleprojres_resolution_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/moduleprojres_resolution.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "moduleprojres_resolution",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": resolution_case(&Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": resolution_case(&field5, |entry| json!(
                field5_coefficients
                    .to_integer(entry)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
            )),
            "rational": resolution_case(&Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn pinned_qpa_projective_path_module_lift_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/projpathalgmodule_lift.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "projpathalgmodule_lift",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": lift_case(&Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": lift_case(&field5, |entry| json!(
                field5_coefficients
                    .to_integer(entry)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
            )),
            "rational": lift_case(&Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn resolution_boundary_checks_are_typed() {
        let parent = PathAlgebra::new(
            Z2,
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let generators = parent.generators();
        let relation = generators[1].multiply(&generators[1]).unwrap();
        let unattached = PathAlgebraIdeal::new(parent.clone(), vec![relation.clone()]).unwrap();
        let ideal = unattached
            .with_complete_groebner_basis(std::slice::from_ref(&relation))
            .unwrap();
        let algebra = QuotientPathAlgebra::by_ideal(&ideal).unwrap();
        let quotient_generators = algebra.generators().unwrap();
        let x = quotient_generators[1].clone();

        assert!(matches!(
            projective_resolution_fp_path_algebra_module(&algebra, &ideal, &[vec![x.clone()]], 0),
            Err(ProjectiveResolutionError::NonpositiveTerm)
        ));
        assert!(matches!(
            projective_resolution_fp_path_algebra_module(&algebra, &ideal, &[], 1),
            Err(ProjectiveResolutionError::EmptyPresentationMap)
        ));
        assert!(matches!(
            projective_resolution_fp_path_algebra_module(
                &algebra,
                &ideal,
                &[vec![x.clone()], vec![]],
                1,
            ),
            Err(ProjectiveResolutionError::RaggedPresentationMap { row: 2, .. })
        ));
        assert!(matches!(
            projective_resolution_fp_path_algebra_module(
                &algebra,
                &ideal,
                &[vec![algebra.zero()]],
                1,
            ),
            Err(ProjectiveResolutionError::ZeroPresentationColumn { column: 1 })
        ));

        let quotient_projective = right_projective_quotient_path_algebra_module(
            &algebra,
            std::slice::from_ref(&quotient_generators[0]),
        )
        .unwrap();
        assert_eq!(quotient_projective.vertex_indices(), &[vertex(1)]);
        assert!(matches!(
            right_projective_quotient_path_algebra_module(&algebra, std::slice::from_ref(&x)),
            Err(ProjectiveResolutionError::NotEmbeddedQuotientVertex { position: 1 })
        ));
        assert!(matches!(
            right_projective_quotient_path_algebra_module(&algebra, &[]),
            Err(ProjectiveResolutionError::EmptyQuotientProjectiveVertices)
        ));
        assert!(matches!(
            projective_resolution_fp_path_algebra_module(
                &algebra,
                &ideal,
                &[vec![quotient_generators[0].clone()]],
                3,
            ),
            Err(ProjectiveResolutionError::EmptyResolutionProjective { degree: 2 })
        ));
    }
}
