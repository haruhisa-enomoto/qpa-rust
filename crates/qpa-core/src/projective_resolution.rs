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
    CachedAttribute, DenseMatrix, Path, PathAlgebra, PathAlgebraElement, PathAlgebraIdeal,
    PathAlgebraVector, PathEndpoint, ProjectivePathModuleError, QuotientError, QuotientPathAlgebra,
    QuotientPathAlgebraElement, RightProjectivePathAlgebraModule,
    RightProjectivePathAlgebraSubmodule, VertexIndex, first_part,
    right_projective_path_algebra_module, x_set_of_path_algebra_vector,
};

/// Rust logging target corresponding to QPA's projective-resolution info class.
pub const INFO_PROJECTIVE_RESOLUTION_FP_PATH_ALGEBRA_MODULE: &str =
    "InfoProjectiveResolutionFpPathAlgebraModule";

/// QPA leaves its projective-resolution information class at level zero.
pub const INFO_PROJECTIVE_RESOLUTION_FP_PATH_ALGEBRA_MODULE_DEFAULT_LEVEL: usize = 0;

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

    /// Return QPA's canonical module zero.
    #[must_use]
    pub fn zero(&self) -> QuotientPathAlgebraVector<F> {
        QuotientPathAlgebraVector::new_unchecked(
            self.clone(),
            vec![self.algebra().zero(); self.vertex_indices().len()],
            self.vertex_indices().len() + 1,
        )
    }

    /// Return QPA's coordinate embedded-vertex generators.
    ///
    /// # Errors
    ///
    /// Propagates quotient vertex embedding failures.
    pub fn generators(
        &self,
    ) -> Result<Vec<QuotientPathAlgebraVector<F>>, ProjectiveResolutionError> {
        self.vertex_indices()
            .iter()
            .copied()
            .enumerate()
            .map(|(position, vertex)| {
                let path = self
                    .algebra()
                    .original_path_algebra()
                    .quiver()
                    .vertex(vertex)
                    .ok_or(ProjectiveResolutionError::VertexOutOfRange { vertex })?
                    .path();
                let mut components = vec![self.algebra().zero(); self.vertex_indices().len()];
                components[position] = self.algebra().path(path)?;
                Ok(QuotientPathAlgebraVector::new_unchecked(
                    self.clone(),
                    components,
                    position + 1,
                ))
            })
            .collect()
    }

    /// Translate QPA's checked `PathAlgebraVector`/`Vectorize` constructor for
    /// quotient-projective modules.
    ///
    /// # Errors
    ///
    /// Requires the fixed component count, identical quotient family, and
    /// QPA's left-uniform source-vertex conditions.
    pub fn vector(
        &self,
        components: Vec<QuotientPathAlgebraElement<F>>,
    ) -> Result<QuotientPathAlgebraVector<F>, ProjectiveResolutionError> {
        if components.len() != self.vertex_indices().len() {
            return Err(ProjectiveResolutionError::QuotientVectorComponentCount {
                actual: components.len(),
                expected: self.vertex_indices().len(),
            });
        }
        if components
            .iter()
            .any(|component| !component.algebra().is_identical(self.algebra()))
        {
            return Err(ProjectiveResolutionError::ForeignQuotientVectorComponent);
        }
        for (position, component) in components.iter().enumerate().rev() {
            if !component.is_left_uniform() {
                return Err(
                    ProjectiveResolutionError::QuotientVectorComponentNotLeftUniform {
                        position: position + 1,
                    },
                );
            }
            if component.is_zero() {
                continue;
            }
            let PathEndpoint::Vertex(source) = component.leading_monomial().source() else {
                return Err(
                    ProjectiveResolutionError::QuotientVectorComponentSourceMismatch {
                        position: position + 1,
                    },
                );
            };
            if source.index() != self.vertex_indices()[position] {
                return Err(
                    ProjectiveResolutionError::QuotientVectorComponentSourceMismatch {
                        position: position + 1,
                    },
                );
            }
        }
        let leading_position = quotient_vector_leading_position(self, &components)?;
        Ok(QuotientPathAlgebraVector::new_unchecked(
            self.clone(),
            components,
            leading_position,
        ))
    }

    /// Translate QPA's unchecked positional constructor for this vector family.
    #[must_use]
    pub fn vector_nc(
        &self,
        components: Vec<QuotientPathAlgebraElement<F>>,
        leading_position: usize,
        _normal: bool,
    ) -> QuotientPathAlgebraVector<F> {
        QuotientPathAlgebraVector::new_unchecked(self.clone(), components, leading_position)
    }

    /// Translate the canonical basis of a quotient vertex-projective module.
    ///
    /// `None` is QPA's canonical-basis object without basis vectors.
    ///
    /// # Errors
    ///
    /// Propagates quotient basis, multiplication, vector, and ordering errors.
    pub fn canonical_basis(
        &self,
    ) -> Result<Option<Vec<QuotientPathAlgebraVector<F>>>, ProjectiveResolutionError> {
        let Some(algebra_basis) = self.algebra().canonical_basis()? else {
            return Ok(None);
        };
        let generators = self.generators()?;
        let mut basis = Vec::new();
        for (position, generator) in generators.iter().enumerate() {
            let vertex = &generator.components()[position];
            for element in &algebra_basis {
                if vertex.multiply(element)?.is_zero() {
                    continue;
                }
                let mut components = vec![self.algebra().zero(); self.vertex_indices().len()];
                components[position] = element.clone();
                basis.push(self.vector(components)?);
            }
        }
        sort_quotient_path_algebra_vectors(&mut basis)?;
        Ok(Some(basis))
    }

    /// Return the finite QPA module dimension when its canonical basis exists.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::canonical_basis`] failures.
    pub fn dimension(&self) -> Result<Option<usize>, ProjectiveResolutionError> {
        Ok(self.canonical_basis()?.map(|basis| basis.len()))
    }

    /// Return QPA's installed vertex-projective property for this constructor.
    #[must_use]
    pub const fn is_vertex_projective_module(&self) -> bool {
        true
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

/// QPA's positional vector for a quotient vertex-projective module.
#[derive(Clone, Debug)]
pub struct QuotientPathAlgebraVector<F: Field> {
    module: RightProjectiveQuotientPathAlgebraModule<F>,
    components: Vec<QuotientPathAlgebraElement<F>>,
    leading_position: usize,
}

impl<F: Field> QuotientPathAlgebraVector<F> {
    fn new_unchecked(
        module: RightProjectiveQuotientPathAlgebraModule<F>,
        components: Vec<QuotientPathAlgebraElement<F>>,
        leading_position: usize,
    ) -> Self {
        Self {
            module,
            components,
            leading_position,
        }
    }

    /// Return the identical quotient-projective vector family.
    #[must_use]
    pub const fn module(&self) -> &RightProjectiveQuotientPathAlgebraModule<F> {
        &self.module
    }

    /// Return the quotient-algebra coordinates.
    #[must_use]
    pub fn components(&self) -> &[QuotientPathAlgebraElement<F>] {
        &self.components
    }

    /// Translate QPA's generic path-vector `PrintObj` method.
    ///
    /// GAP prints the stored coordinate list; quotient coordinates therefore
    /// retain the brackets supplied by their own `PrintObj` overload.
    #[must_use]
    pub fn print_obj_with(&self, coefficient: impl Fn(&F::Element) -> String + Copy) -> String {
        let entries = self
            .components()
            .iter()
            .map(|component| component.print_obj_with(coefficient))
            .collect::<Vec<_>>()
            .join(", ");
        format!("[ {entries} ]")
    }

    /// Return QPA's cached one-based leading position.
    #[must_use]
    pub const fn leading_position(&self) -> usize {
        self.leading_position
    }

    /// Return the vector family's stored `IsNormalForm` filter.
    #[must_use]
    pub fn is_normal_form(&self) -> bool {
        self.module.algebra().has_normal_form_function()
    }

    /// Return whether all quotient coordinates vanish.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.components
            .iter()
            .all(QuotientPathAlgebraElement::is_zero)
    }

    /// Add vectors from an identical QPA family and recompute the tip.
    ///
    /// # Errors
    ///
    /// Rejects distinct module families and propagates quotient arithmetic.
    pub fn add(&self, right: &Self) -> Result<Self, ProjectiveResolutionError> {
        self.ensure_same_family(right)?;
        self.module.vector(
            self.components
                .iter()
                .zip(&right.components)
                .map(|(left, right)| left.add(right))
                .collect::<Result<Vec<_>, _>>()?,
        )
    }

    /// Return the additive inverse while preserving QPA's cached position.
    ///
    /// # Errors
    ///
    /// Propagates quotient negation.
    pub fn neg(&self) -> Result<Self, ProjectiveResolutionError> {
        Ok(Self::new_unchecked(
            self.module.clone(),
            self.components
                .iter()
                .map(QuotientPathAlgebraElement::neg)
                .collect::<Result<Vec<_>, _>>()?,
            self.leading_position,
        ))
    }

    /// Subtract vectors from an identical family and recompute the tip.
    ///
    /// # Errors
    ///
    /// Rejects distinct families and propagates quotient arithmetic.
    pub fn sub(&self, right: &Self) -> Result<Self, ProjectiveResolutionError> {
        self.ensure_same_family(right)?;
        self.module.vector(
            self.components
                .iter()
                .zip(&right.components)
                .map(|(left, right)| left.sub(right))
                .collect::<Result<Vec<_>, _>>()?,
        )
    }

    /// Scale every coordinate on the left.
    ///
    /// # Errors
    ///
    /// Propagates quotient scaling and normalization.
    pub fn scale_left(&self, scalar: &F::Element) -> Result<Self, ProjectiveResolutionError> {
        if self
            .module
            .algebra()
            .original_path_algebra()
            .field()
            .is_zero(scalar)
        {
            return Ok(self.module.zero());
        }
        Ok(Self::new_unchecked(
            self.module.clone(),
            self.components
                .iter()
                .map(|component| component.scale_left(scalar))
                .collect::<Result<Vec<_>, _>>()?,
            self.leading_position,
        ))
    }

    /// Scale every coordinate on the right.
    ///
    /// QPA implements both scalar actions by the same componentwise field
    /// multiplication, so this delegates to [`Self::scale_left`].
    ///
    /// # Errors
    ///
    /// Propagates quotient scaling and normalization.
    pub fn scale_right(&self, scalar: &F::Element) -> Result<Self, ProjectiveResolutionError> {
        self.scale_left(scalar)
    }

    /// Apply QPA's coordinatewise right quotient-algebra action.
    ///
    /// # Errors
    ///
    /// Rejects a foreign quotient and propagates multiplication and checked
    /// vector construction.
    pub fn act(
        &self,
        element: &QuotientPathAlgebraElement<F>,
    ) -> Result<Self, ProjectiveResolutionError> {
        if !element.algebra().is_identical(self.module.algebra()) {
            return Err(ProjectiveResolutionError::ForeignQuotientVectorAction);
        }
        self.module.vector(
            self.components
                .iter()
                .map(|component| component.multiply(element))
                .collect::<Result<Vec<_>, _>>()?,
        )
    }

    /// Return QPA's sparse quotient-vector leading term.
    #[must_use]
    pub fn leading_term(&self) -> Self {
        let Some(position) = self.leading_position.checked_sub(1) else {
            return self.clone();
        };
        let Some(component) = self.components.get(position) else {
            return self.clone();
        };
        let mut components = vec![self.module.algebra().zero(); self.components.len()];
        components[position] = component.leading_term();
        Self::new_unchecked(self.module.clone(), components, self.leading_position)
    }

    /// Return QPA's vector-tip coefficient.
    ///
    /// # Errors
    ///
    /// QPA indexes coordinate one for zero and fails on an empty vector.
    pub fn leading_coefficient(&self) -> Result<F::Element, ProjectiveResolutionError> {
        let position = if self.leading_position <= self.components.len() {
            self.leading_position - 1
        } else {
            0
        };
        self.components
            .get(position)
            .map(QuotientPathAlgebraElement::leading_coefficient)
            .ok_or(ProjectiveResolutionError::EmptyQuotientVectorLeadingComponent)
    }

    /// Return the coordinate containing QPA's vector tip.
    ///
    /// # Errors
    ///
    /// QPA indexes coordinate one for zero and fails on an empty vector.
    pub fn leading_component(
        &self,
    ) -> Result<&QuotientPathAlgebraElement<F>, ProjectiveResolutionError> {
        let position = if self.is_zero() {
            0
        } else {
            self.leading_position - 1
        };
        self.components
            .get(position)
            .ok_or(ProjectiveResolutionError::EmptyQuotientVectorLeadingComponent)
    }

    /// Translate QPA's vector `TargetVertex`; `None` is GAP's `fail`.
    #[must_use]
    pub fn target_vertex(&self) -> Option<Path> {
        if !self.is_right_uniform() {
            return None;
        }
        if self.is_zero() {
            return Some(
                self.module
                    .algebra()
                    .original_path_algebra()
                    .quiver()
                    .zero(),
            );
        }
        match self.leading_path().ok()?.target() {
            PathEndpoint::Vertex(vertex) => Some(vertex.path()),
            PathEndpoint::Zero(path) => Some(path),
        }
    }

    /// Translate QPA's vector `SourceVertex`; `None` is GAP's `fail`.
    ///
    /// The source method surprisingly checks right uniformity, and that
    /// condition is retained.
    #[must_use]
    pub fn source_vertex(&self) -> Option<Path> {
        if !self.is_right_uniform() {
            return None;
        }
        if self.is_zero() {
            return Some(
                self.module
                    .algebra()
                    .original_path_algebra()
                    .quiver()
                    .zero(),
            );
        }
        match self.leading_path().ok()?.source() {
            PathEndpoint::Vertex(vertex) => Some(vertex.path()),
            PathEndpoint::Zero(path) => Some(path),
        }
    }

    /// Return whether all components vanish from one common source vertex.
    #[must_use]
    pub fn is_left_uniform(&self) -> bool {
        self.is_uniform_at(Path::source, QuotientPathAlgebraElement::is_left_uniform)
    }

    /// Return whether all components vanish at one common target vertex.
    #[must_use]
    pub fn is_right_uniform(&self) -> bool {
        self.is_uniform_at(Path::target, QuotientPathAlgebraElement::is_right_uniform)
    }

    /// Return whether the vector is both left and right uniform.
    #[must_use]
    pub fn is_uniform(&self) -> bool {
        self.is_left_uniform() && self.is_right_uniform()
    }

    /// Compare vectors by tip monomial and descending component position.
    ///
    /// # Errors
    ///
    /// Rejects distinct vector families and propagates path ordering.
    pub fn compare(&self, right: &Self) -> Result<Ordering, ProjectiveResolutionError> {
        self.ensure_same_family(right)?;
        if !self.is_normal_form() || !right.is_normal_form() {
            return Err(ProjectiveResolutionError::QuotientVectorsNotInNormalForm);
        }
        let path_order = self
            .module
            .algebra()
            .original_path_algebra()
            .ordering()
            .compare(&self.leading_path()?, &right.leading_path()?)?;
        if path_order.is_eq() {
            Ok(right.leading_position.cmp(&self.leading_position))
        } else {
            Ok(path_order)
        }
    }

    /// Translate QPA's `IsLeftDivisible(x, y)` predicate.
    ///
    /// The source additionally requires the dividend tip walk to be nonempty,
    /// so a vertex tip is never left-divisible here.
    ///
    /// # Errors
    ///
    /// Rejects distinct vector families and explicit zero inputs.
    pub fn is_left_divisible(&self, divisor: &Self) -> Result<bool, ProjectiveResolutionError> {
        self.ensure_same_family(divisor)?;
        if self.is_zero() || divisor.is_zero() {
            return Err(ProjectiveResolutionError::ZeroQuotientVectorDivision);
        }
        let dividend_path = self.leading_path()?;
        let divisor_path = divisor.leading_path()?;
        let dividend_walk = dividend_path
            .walk()
            .ok_or(ProjectivePathModuleError::PathHasNoArrowWalk)?;
        let divisor_walk = divisor_path
            .walk()
            .ok_or(ProjectivePathModuleError::PathHasNoArrowWalk)?;
        Ok(self.leading_position == divisor.leading_position
            && !dividend_walk.is_empty()
            && dividend_walk.starts_with(divisor_walk))
    }

    fn ensure_same_family(&self, right: &Self) -> Result<(), ProjectiveResolutionError> {
        if self.module.is_identical(&right.module) {
            Ok(())
        } else {
            Err(ProjectiveResolutionError::DifferentQuotientVectorFamilies)
        }
    }

    fn leading_path(&self) -> Result<crate::Path, ProjectiveResolutionError> {
        if self.is_zero() {
            return Err(
                ProjectiveResolutionError::InvalidQuotientVectorLeadingPosition {
                    position: self.leading_position,
                    length: self.components.len(),
                },
            );
        }
        self.components
            .get(self.leading_position - 1)
            .map(QuotientPathAlgebraElement::leading_monomial)
            .ok_or(
                ProjectiveResolutionError::InvalidQuotientVectorLeadingPosition {
                    position: self.leading_position,
                    length: self.components.len(),
                },
            )
    }

    fn is_uniform_at(
        &self,
        endpoint: impl Fn(&crate::Path) -> PathEndpoint,
        component_is_uniform: impl Fn(&QuotientPathAlgebraElement<F>) -> bool,
    ) -> bool {
        if self.is_zero() {
            return true;
        }
        if self
            .components
            .iter()
            .any(|component| !component_is_uniform(component))
        {
            return false;
        }
        let mut endpoints = self
            .components
            .iter()
            .filter(|component| !component.is_zero())
            .filter_map(|component| {
                component.representative().terms().first().and_then(|term| {
                    match endpoint(term.path()) {
                        PathEndpoint::Vertex(vertex) => Some(vertex),
                        PathEndpoint::Zero(_) => None,
                    }
                })
            });
        let Some(first) = endpoints.next() else {
            return true;
        };
        endpoints.all(|candidate| candidate == first)
    }
}

impl<F: Field> PartialEq for QuotientPathAlgebraVector<F> {
    fn eq(&self, other: &Self) -> bool {
        self.module.is_identical(&other.module) && self.components == other.components
    }
}

impl<F: Field> Eq for QuotientPathAlgebraVector<F> {}

fn quotient_vector_leading_position<F: Field>(
    module: &RightProjectiveQuotientPathAlgebraModule<F>,
    components: &[QuotientPathAlgebraElement<F>],
) -> Result<usize, ProjectiveResolutionError> {
    let Some(first) = components.first() else {
        return Ok(1);
    };
    let ordering = module.algebra().original_path_algebra().ordering();
    let mut position = 1;
    let mut greatest = first.leading_monomial();
    for (index, component) in components.iter().enumerate().skip(1) {
        let candidate = component.leading_monomial();
        if ordering.compare(&candidate, &greatest)?.is_gt() {
            greatest = candidate;
            position = index + 1;
        }
    }
    Ok(position)
}

fn sort_quotient_path_algebra_vectors<F: Field>(
    vectors: &mut Vec<QuotientPathAlgebraVector<F>>,
) -> Result<(), ProjectiveResolutionError> {
    let entered = std::mem::take(vectors);
    for vector in entered {
        let mut insertion = vectors.len();
        for (position, existing) in vectors.iter().enumerate() {
            if vector.compare(existing)?.is_lt() {
                insertion = position;
                break;
            }
        }
        vectors.insert(insertion, vector);
    }
    Ok(())
}

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
    r_projectives_vertex_list: CachedAttribute<Vec<Vec<VertexIndex>>>,
    projectives_f_prime_list: CachedAttribute<Vec<Vec<PathAlgebraVector<F>>>>,
    finite_at_projective: Option<usize>,
}

/// Homology after tensoring a projective resolution with the vertex
/// semisimple quotient of an admissible path-algebra quotient.
///
/// Rows are ordered by homological degree and columns by the one-based quiver
/// vertex order.  `differential_vertex_ranks[0]` is the zero differential;
/// row `n` for `n > 0` records the vertexwise rank of `d_n`.  Homology is
/// returned only through degrees for which both adjacent differentials were
/// computed.  Its dimensions are the vertex multiplicities in the minimal
/// projective resolution.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SemisimpleProjectiveResolutionHomology {
    projective_vertex_multiplicities: Vec<Vec<usize>>,
    differential_vertex_ranks: Vec<Vec<usize>>,
    homology_vertex_dimensions: Vec<Vec<usize>>,
}

impl SemisimpleProjectiveResolutionHomology {
    /// Return the vertex multiplicities of every available projective term.
    #[must_use]
    pub fn projective_vertex_multiplicities(&self) -> &[Vec<usize>] {
        &self.projective_vertex_multiplicities
    }

    /// Return the vertexwise ranks of the reduced differentials.
    #[must_use]
    pub fn differential_vertex_ranks(&self) -> &[Vec<usize>] {
        &self.differential_vertex_ranks
    }

    /// Return the vertexwise semisimple homology dimensions.
    ///
    /// For an admissible quotient these are the projective multiplicities in
    /// the minimal resolution.
    #[must_use]
    pub fn homology_vertex_dimensions(&self) -> &[Vec<usize>] {
        &self.homology_vertex_dimensions
    }

    /// Return the minimal projective Betti multiplicities.
    ///
    /// This is the same data as [`Self::homology_vertex_dimensions`], named
    /// for its use in bounded projective-dimension decisions.
    #[must_use]
    pub fn minimal_projective_betti_multiplicities(&self) -> &[Vec<usize>] {
        &self.homology_vertex_dimensions
    }
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

    /// Return the bound state of QPA's mutable `RProjectivesVertexList` attribute.
    ///
    /// QPA's constructors never set this separately declared legacy attribute.
    #[must_use]
    pub const fn r_projectives_vertex_list(&self) -> &CachedAttribute<Vec<Vec<VertexIndex>>> {
        &self.r_projectives_vertex_list
    }

    /// Bind QPA's mutable `RProjectivesVertexList` attribute.
    pub fn set_r_projectives_vertex_list(&mut self, value: Vec<Vec<VertexIndex>>) {
        self.r_projectives_vertex_list = CachedAttribute::Value(value);
    }

    /// Return the bound state of QPA's mutable `ProjectivesFPrimeList` attribute.
    ///
    /// This is distinct from the populated second component of
    /// `ProjectivesFList`; QPA's constructors leave this legacy attribute
    /// unbound.
    #[must_use]
    pub const fn projectives_f_prime_list(
        &self,
    ) -> &CachedAttribute<Vec<Vec<PathAlgebraVector<F>>>> {
        &self.projectives_f_prime_list
    }

    /// Bind QPA's mutable `ProjectivesFPrimeList` attribute.
    pub fn set_projectives_f_prime_list(&mut self, value: Vec<Vec<PathAlgebraVector<F>>>) {
        self.projectives_f_prime_list = CachedAttribute::Value(value);
    }

    /// Return the degree at which QPA prints `finite at projective`, if reached.
    #[must_use]
    pub const fn finite_at_projective(&self) -> Option<usize> {
        self.finite_at_projective
    }

    /// Reduce the differential matrices modulo the arrow ideal and compute
    /// their vertexwise homology.
    ///
    /// If `P` is any projective resolution over an admissible quotient
    /// `A = kQ/I`, then `H_n(P tensor_A A/rad(A))` records the multiplicities
    /// of the indecomposable projectives in the minimal resolution.  This
    /// therefore extracts exact minimal Betti data even when the GSZ
    /// resolution stored here is not itself minimal.
    ///
    /// # Errors
    ///
    /// Requires an admissible quotient and internally consistent projective
    /// and differential shapes.  Exact admissibility failures are propagated.
    #[allow(clippy::too_many_lines)]
    pub fn semisimple_homology(
        &self,
    ) -> Result<SemisimpleProjectiveResolutionHomology, ProjectiveResolutionError> {
        if !self.parent_algebra.is_admissible()? {
            return Err(ProjectiveResolutionError::NonadmissibleSemisimpleHomology);
        }

        let quiver = self.parent_algebra.original_path_algebra().quiver();
        let vertex_count = quiver.number_of_vertices();
        let field = self.parent_algebra.original_path_algebra().field();
        let projective_vertices = self
            .r_projectives
            .iter()
            .map(RightProjectivePathAlgebraModule::vertex_indices)
            .collect::<Vec<_>>();
        let projective_vertex_multiplicities = projective_vertices
            .iter()
            .map(|vertices| {
                let mut multiplicities = vec![0_usize; vertex_count];
                for vertex in *vertices {
                    multiplicities[vertex.get() - 1] += 1;
                }
                multiplicities
            })
            .collect::<Vec<_>>();

        let mut differential_vertex_ranks = vec![vec![0_usize; vertex_count]];
        for (map_position, stored_map) in self.maps.iter().enumerate() {
            let degree = map_position + 1;
            let target_vertices = projective_vertices
                .get(degree - 1)
                .ok_or(ProjectiveResolutionError::SemisimpleDifferentialMissingTarget { degree })?;
            // QPA stores the entered projective-presentation matrix at
            // `maps[0]`, while it rebuilds P1 from the Gröbner basis `f^1`.
            // Thus d1 is the component matrix of `f^1`; only later `maps`
            // entries are already the actual differentials between the
            // stored projectives.
            let differential_rows = if degree == 1 {
                self.f_lists
                    .get(1)
                    .ok_or(
                        ProjectiveResolutionError::SemisimpleDifferentialMissingSource { degree },
                    )?
                    .iter()
                    .map(PathAlgebraVector::components)
                    .collect::<Vec<_>>()
            } else {
                stored_map.iter().map(Vec::as_slice).collect::<Vec<_>>()
            };

            if differential_rows.is_empty() {
                if degree < projective_vertices.len() {
                    return Err(ProjectiveResolutionError::SemisimpleDifferentialRowCount {
                        degree,
                        actual: 0,
                        expected: projective_vertices[degree].len(),
                    });
                }
                differential_vertex_ranks.push(vec![0_usize; vertex_count]);
                continue;
            }

            let source_vertices = projective_vertices
                .get(degree)
                .ok_or(ProjectiveResolutionError::SemisimpleDifferentialMissingSource { degree })?;
            if differential_rows.len() != source_vertices.len() {
                return Err(ProjectiveResolutionError::SemisimpleDifferentialRowCount {
                    degree,
                    actual: differential_rows.len(),
                    expected: source_vertices.len(),
                });
            }
            for (row, entries) in differential_rows.iter().enumerate() {
                if entries.len() != target_vertices.len() {
                    return Err(
                        ProjectiveResolutionError::SemisimpleDifferentialColumnCount {
                            degree,
                            row: row + 1,
                            actual: entries.len(),
                            expected: target_vertices.len(),
                        },
                    );
                }
            }

            let mut row_counts = vec![0; vertex_count];
            let row_positions = source_vertices
                .iter()
                .map(|vertex| {
                    let index = vertex.get() - 1;
                    let position = row_counts[index];
                    row_counts[index] += 1;
                    position
                })
                .collect::<Vec<_>>();
            let mut column_counts = vec![0; vertex_count];
            let column_positions = target_vertices
                .iter()
                .map(|vertex| {
                    let index = vertex.get() - 1;
                    let position = column_counts[index];
                    column_counts[index] += 1;
                    position
                })
                .collect::<Vec<_>>();
            let mut vertex_blocks = row_counts
                .iter()
                .zip(&column_counts)
                .map(|(&rows, &columns)| DenseMatrix::zero(rows, columns, field.clone()))
                .collect::<Vec<_>>();

            for (row, entries) in differential_rows.iter().enumerate() {
                for (column, entry) in entries.iter().enumerate() {
                    for term in entry.terms().iter().filter(|term| term.path().is_vertex()) {
                        let PathEndpoint::Vertex(term_vertex) = term.path().source() else {
                            unreachable!("a vertex path has a vertex source");
                        };
                        let term_vertex = term_vertex.index();
                        let source_vertex = source_vertices[row];
                        let target_vertex = target_vertices[column];
                        if term_vertex != source_vertex || term_vertex != target_vertex {
                            return Err(
                                ProjectiveResolutionError::SemisimpleDifferentialVertexMismatch {
                                    degree,
                                    row: row + 1,
                                    column: column + 1,
                                    term_vertex,
                                    source_vertex,
                                    target_vertex,
                                },
                            );
                        }
                        let vertex = term_vertex.get() - 1;
                        field.add_assign(
                            &mut vertex_blocks[vertex]
                                [(row_positions[row], column_positions[column])],
                            term.coefficient(),
                        );
                    }
                }
            }
            differential_vertex_ranks.push(vertex_blocks.iter().map(DenseMatrix::rank).collect());
        }

        let homology_degree_count = projective_vertex_multiplicities.len().min(self.maps.len());
        let mut homology_vertex_dimensions = Vec::with_capacity(homology_degree_count);
        for degree in 0..homology_degree_count {
            let mut dimensions = Vec::with_capacity(vertex_count);
            for vertex in 0..vertex_count {
                let projective = projective_vertex_multiplicities[degree][vertex];
                let differential_rank = differential_vertex_ranks[degree][vertex];
                let next_differential_rank = differential_vertex_ranks[degree + 1][vertex];
                let Some(dimension) = projective
                    .checked_sub(differential_rank)
                    .and_then(|value| value.checked_sub(next_differential_rank))
                else {
                    return Err(
                        ProjectiveResolutionError::InvalidSemisimpleHomologyDimension {
                            degree,
                            vertex: vertex + 1,
                            projective,
                            differential_rank,
                            next_differential_rank,
                        },
                    );
                };
                dimensions.push(dimension);
            }
            homology_vertex_dimensions.push(dimensions);
        }

        Ok(SemisimpleProjectiveResolutionHomology {
            projective_vertex_multiplicities,
            differential_vertex_ranks,
            homology_vertex_dimensions,
        })
    }
}

/// Represent QPA's declared but uninstalled three-argument
/// `ProjectiveResolutionFpPathAlgebraModule(A, I, presentation_map)` overload.
///
/// # Errors
///
/// Pinned QPA installs only the four-argument overload and therefore always
/// reaches its missing-method boundary here.
pub const fn projective_resolution_fp_path_algebra_module_without_term<F: Field>(
    _algebra: &QuotientPathAlgebra<F>,
    _ideal: &PathAlgebraIdeal<F>,
    _presentation_map: &[Vec<QuotientPathAlgebraElement<F>>],
) -> Result<ProjectiveResolutionFpPathAlgebraModule<F>, ProjectiveResolutionError> {
    Err(ProjectiveResolutionError::UninstalledThreeArgumentConstructor)
}

/// Represent QPA's declared but uninstalled `FindNextRProjective` operation.
///
/// # Errors
///
/// Always returns QPA's explicit missing-method boundary.
pub const fn find_next_r_projective<F: Field>(
    _resolution: &ProjectiveResolutionFpPathAlgebraModule<F>,
    _ring: &PathAlgebra<F>,
    _term: usize,
) -> Result<(), ProjectiveResolutionError> {
    Err(ProjectiveResolutionError::UninstalledFindNextRProjective)
}

/// Represent QPA's declared but uninstalled `FindNextSyzygy` operation.
///
/// # Errors
///
/// Always returns QPA's explicit missing-method boundary.
pub const fn find_next_syzygy<F: Field>(
    _resolution: &ProjectiveResolutionFpPathAlgebraModule<F>,
    _term: usize,
) -> Result<(), ProjectiveResolutionError> {
    Err(ProjectiveResolutionError::UninstalledFindNextSyzygy)
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
        r_projectives_vertex_list: CachedAttribute::Unbound,
        projectives_f_prime_list: CachedAttribute::Unbound,
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
        r_projectives_vertex_list: CachedAttribute::Unbound,
        projectives_f_prime_list: CachedAttribute::Unbound,
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
    /// QPA declares the three-argument overload but installs no method for it.
    #[error("QPA installs no three-argument ProjectiveResolutionFpPathAlgebraModule method")]
    UninstalledThreeArgumentConstructor,
    /// QPA declares `FindNextRProjective` but installs no method for it.
    #[error("QPA installs no method for FindNextRProjective")]
    UninstalledFindNextRProjective,
    /// QPA declares `FindNextSyzygy` but installs no method for it.
    #[error("QPA installs no method for FindNextSyzygy")]
    UninstalledFindNextSyzygy,
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
    /// Vertex semisimple homology represents the radical quotient only for an
    /// admissible path-algebra quotient.
    #[error("semisimple resolution homology requires an admissible quotient")]
    NonadmissibleSemisimpleHomology,
    /// A differential has no stored target projective.
    #[error("degree-{degree} differential has no stored target projective")]
    SemisimpleDifferentialMissingTarget {
        /// Differential degree.
        degree: usize,
    },
    /// A nonzero differential has no stored source projective.
    #[error("degree-{degree} differential has no stored source projective")]
    SemisimpleDifferentialMissingSource {
        /// Differential degree.
        degree: usize,
    },
    /// A differential row count disagrees with its source projective.
    #[error(
        "degree-{degree} differential has {actual} rows, expected {expected} from its source projective"
    )]
    SemisimpleDifferentialRowCount {
        /// Differential degree.
        degree: usize,
        /// Stored matrix row count.
        actual: usize,
        /// Source-projective summand count.
        expected: usize,
    },
    /// A differential column count disagrees with its target projective.
    #[error(
        "degree-{degree} differential row {row} has {actual} columns, expected {expected} from its target projective"
    )]
    SemisimpleDifferentialColumnCount {
        /// Differential degree.
        degree: usize,
        /// One-based matrix row.
        row: usize,
        /// Stored matrix column count.
        actual: usize,
        /// Target-projective summand count.
        expected: usize,
    },
    /// A degree-zero path term must lie in the same vertex block on source and
    /// target.
    #[error(
        "degree-{degree} differential entry ({row}, {column}) has vertex term {term_vertex:?}, source {source_vertex:?}, target {target_vertex:?}"
    )]
    SemisimpleDifferentialVertexMismatch {
        /// Differential degree.
        degree: usize,
        /// One-based matrix row.
        row: usize,
        /// One-based matrix column.
        column: usize,
        /// Vertex of the degree-zero path term.
        term_vertex: VertexIndex,
        /// Vertex of the source projective summand.
        source_vertex: VertexIndex,
        /// Vertex of the target projective summand.
        target_vertex: VertexIndex,
    },
    /// Adjacent reduced differential ranks exceeded the corresponding
    /// projective multiplicity.
    #[error(
        "degree-{degree} vertex {vertex} has projective multiplicity {projective}, adjacent ranks {differential_rank} and {next_differential_rank}"
    )]
    InvalidSemisimpleHomologyDimension {
        /// Homological degree.
        degree: usize,
        /// Quiver vertex.
        vertex: usize,
        /// Projective summand multiplicity.
        projective: usize,
        /// Rank of `d_degree`.
        differential_rank: usize,
        /// Rank of `d_(degree+1)`.
        next_differential_rank: usize,
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
    /// Checked quotient vectors use the project's fixed coordinate count.
    #[error("quotient path vector has {actual} components, expected {expected}")]
    QuotientVectorComponentCount {
        /// Entered coordinate count.
        actual: usize,
        /// Required coordinate count.
        expected: usize,
    },
    /// Every coordinate must belong to the module's quotient algebra.
    #[error("a quotient path-vector component belongs to another quotient algebra")]
    ForeignQuotientVectorComponent,
    /// QPA checks left uniformity from the last coordinate backwards.
    #[error("quotient path-vector component {position} is not left uniform")]
    QuotientVectorComponentNotLeftUniform {
        /// One-based coordinate.
        position: usize,
    },
    /// A nonzero coordinate must start at its projective summand vertex.
    #[error("quotient path-vector component {position} starts at the wrong vertex")]
    QuotientVectorComponentSourceMismatch {
        /// One-based coordinate.
        position: usize,
    },
    /// Vector arithmetic requires an identical quotient-projective family.
    #[error("the quotient path vectors belong to different module families")]
    DifferentQuotientVectorFamilies,
    /// QPA only installs quotient-vector ordering for normal-form vectors.
    #[error("quotient path-vector ordering requires normal forms")]
    QuotientVectorsNotInNormalForm,
    /// Right action requires the identical quotient algebra.
    #[error("the quotient vector action belongs to another quotient algebra")]
    ForeignQuotientVectorAction,
    /// QPA's divisibility method has no tip for a zero vector.
    #[error("zero quotient path vectors have no divisibility tip")]
    ZeroQuotientVectorDivision,
    /// QPA indexes coordinate one for zero and fails on a zero-rank vector.
    #[error("an empty quotient path vector has no leading component")]
    EmptyQuotientVectorLeadingComponent,
    /// An unchecked positional vector stored an invalid one-based tip index.
    #[error("leading position {position} is outside a quotient vector of length {length}")]
    InvalidQuotientVectorLeadingPosition {
        /// Stored one-based position.
        position: usize,
        /// Coordinate count.
        length: usize,
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
    /// QPA path ordering failed.
    #[error(transparent)]
    Ordering(#[from] crate::OrderingError),
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
        INFO_PROJECTIVE_RESOLUTION_FP_PATH_ALGEBRA_MODULE,
        INFO_PROJECTIVE_RESOLUTION_FP_PATH_ALGEBRA_MODULE_DEFAULT_LEVEL, ProjectiveResolutionError,
        ProjectiveResolutionFpPathAlgebraModule, find_next_r_projective, find_next_syzygy,
        projective_resolution_fp_path_algebra_module,
        projective_resolution_fp_path_algebra_module_without_term,
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

    fn quotient_vector_json<F: Field>(
        vector: &super::QuotientPathAlgebraVector<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            vector
                .components()
                .iter()
                .map(|component| element_json(component.representative(), coefficient))
                .collect(),
        )
    }

    fn quotient_vector_list_json<F: Field>(
        vectors: &[super::QuotientPathAlgebraVector<F>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            vectors
                .iter()
                .map(|vector| quotient_vector_json(vector, coefficient))
                .collect(),
        )
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

    fn quotient_vector_case<F: Field>(
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
        let ideal = PathAlgebraIdeal::new(parent, vec![relation.clone()])
            .unwrap()
            .with_complete_groebner_basis(std::slice::from_ref(&relation))
            .unwrap();
        let algebra = QuotientPathAlgebra::by_ideal(&ideal).unwrap();
        let quotient_generators = algebra.generators().unwrap();
        let projective = right_projective_quotient_path_algebra_module(
            &algebra,
            &[
                quotient_generators[0].clone(),
                quotient_generators[1].clone(),
            ],
        )
        .unwrap();
        let vector = projective
            .vector(vec![
                quotient_generators[3].clone(),
                quotient_generators[1].clone(),
            ])
            .unwrap();
        let other = projective
            .vector(vec![quotient_generators[0].clone(), algebra.zero()])
            .unwrap();
        let checked_zero = projective
            .vector(vec![algebra.zero(), algebra.zero()])
            .unwrap();
        let nc = projective.vector_nc(vector.components().to_vec(), 2, false);
        let canonical_basis = projective.canonical_basis().unwrap().unwrap();
        let scalar = field.add(&field.one(), &field.one());
        json!({
            "dimension": projective.dimension().unwrap().unwrap(),
            "generators": quotient_vector_list_json(&projective.generators().unwrap(), coefficient),
            "canonical_basis": quotient_vector_list_json(&canonical_basis, coefficient),
            "vector": quotient_vector_json(&vector, coefficient),
            "other": quotient_vector_json(&other, coefficient),
            "sum": quotient_vector_json(&vector.add(&other).unwrap(), coefficient),
            "difference": quotient_vector_json(&vector.sub(&other).unwrap(), coefficient),
            "negated": quotient_vector_json(&vector.neg().unwrap(), coefficient),
            "scaled_left": quotient_vector_json(
                &vector.scale_left(&scalar).unwrap(),
                coefficient,
            ),
            "scaled_right": quotient_vector_json(
                &vector.scale_right(&scalar).unwrap(),
                coefficient,
            ),
            "action_b": quotient_vector_json(
                &vector.act(&quotient_generators[4]).unwrap(),
                coefficient,
            ),
            "leading_term": quotient_vector_json(&vector.leading_term(), coefficient),
            "leading_component": element_json(
                vector.leading_component().unwrap().representative(),
                coefficient,
            ),
            "leading_position": vector.leading_position(),
            "leading_coefficient": coefficient(&vector.leading_coefficient().unwrap()),
            "target_vertex": vector.target_vertex().unwrap().external_representation(),
            "source_vertex": vector.source_vertex().unwrap().external_representation(),
            "left_divisibility": [
                vector.is_left_divisible(&other).unwrap(),
                other.is_left_divisible(&vector).unwrap(),
            ],
            "uniformity": [
                vector.is_left_uniform(),
                vector.is_right_uniform(),
                vector.is_uniform(),
            ],
            "normal_form": vector.is_normal_form(),
            "zero_positions": [
                projective.zero().leading_position(),
                checked_zero.leading_position(),
            ],
            "nc": {
                "vector": quotient_vector_json(&nc, coefficient),
                "leading_position": nc.leading_position(),
            },
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
    fn pinned_qpa_projective_quotient_vector_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/projpathalgmodule_quotient_vectors.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "projpathalgmodule_quotient_vectors",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": quotient_vector_case(&Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": quotient_vector_case(&field5, |entry| json!(
                field5_coefficients
                    .to_integer(entry)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
            )),
            "rational": quotient_vector_case(&Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn quotient_projective_vector_boundaries_are_typed() {
        let parent = PathAlgebra::new(
            Z2,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap();
        let generators = parent.generators();
        let relation = generators[2].clone();
        let unattached = PathAlgebraIdeal::new(parent, vec![relation.clone()]).unwrap();
        let ideal = unattached
            .with_complete_groebner_basis(std::slice::from_ref(&relation))
            .unwrap();
        let algebra = QuotientPathAlgebra::by_ideal(&ideal).unwrap();
        let quotient_generators = algebra.generators().unwrap();
        let projective = right_projective_quotient_path_algebra_module(
            &algebra,
            &[
                quotient_generators[0].clone(),
                quotient_generators[1].clone(),
            ],
        )
        .unwrap();

        assert!(matches!(
            projective.vector(vec![quotient_generators[0].clone()]),
            Err(ProjectiveResolutionError::QuotientVectorComponentCount {
                actual: 1,
                expected: 2,
            })
        ));
        assert!(matches!(
            projective.vector(vec![quotient_generators[1].clone(), algebra.zero()]),
            Err(ProjectiveResolutionError::QuotientVectorComponentSourceMismatch { position: 1 })
        ));
        assert!(matches!(
            projective.vector(vec![
                quotient_generators[0].add(&quotient_generators[1]).unwrap(),
                algebra.zero(),
            ]),
            Err(ProjectiveResolutionError::QuotientVectorComponentNotLeftUniform { position: 1 })
        ));

        let other_projective = right_projective_quotient_path_algebra_module(
            &algebra,
            &[
                quotient_generators[0].clone(),
                quotient_generators[1].clone(),
            ],
        )
        .unwrap();
        let vector = projective.generators().unwrap()[0].clone();
        let other_vector = other_projective.generators().unwrap()[0].clone();
        assert!(matches!(
            vector.add(&other_vector),
            Err(ProjectiveResolutionError::DifferentQuotientVectorFamilies)
        ));
        assert!(matches!(
            vector.is_left_divisible(&projective.zero()),
            Err(ProjectiveResolutionError::ZeroQuotientVectorDivision)
        ));

        let foreign_algebra = QuotientPathAlgebra::by_ideal(&ideal).unwrap();
        let foreign_vertex = foreign_algebra.generators().unwrap()[0].clone();
        assert!(matches!(
            projective.vector(vec![foreign_vertex.clone(), foreign_algebra.zero()]),
            Err(ProjectiveResolutionError::ForeignQuotientVectorComponent)
        ));
        assert!(matches!(
            vector.act(&foreign_vertex),
            Err(ProjectiveResolutionError::ForeignQuotientVectorAction)
        ));

        let incomplete_ideal =
            PathAlgebraIdeal::new(ideal.algebra().clone(), vec![relation]).unwrap();
        let incomplete_algebra = QuotientPathAlgebra::by_ideal(&incomplete_ideal).unwrap();
        let incomplete_generators = incomplete_algebra.generators().unwrap();
        let incomplete_projective = right_projective_quotient_path_algebra_module(
            &incomplete_algebra,
            std::slice::from_ref(&incomplete_generators[0]),
        )
        .unwrap();
        let incomplete_vector = incomplete_projective.generators().unwrap()[0].clone();
        assert!(matches!(
            incomplete_vector.compare(&incomplete_vector),
            Err(ProjectiveResolutionError::QuotientVectorsNotInNormalForm)
        ));
    }

    #[test]
    fn semisimple_homology_cancels_a_nonminimal_gsz_pair() {
        let parent = PathAlgebra::new(
            Z2,
            Quiver::from_indexed(
                5,
                vec![
                    ArrowSpec::named(vertex(1), vertex(2), "a"),
                    ArrowSpec::named(vertex(2), vertex(4), "b"),
                    ArrowSpec::named(vertex(1), vertex(3), "c"),
                    ArrowSpec::named(vertex(3), vertex(4), "d"),
                    ArrowSpec::named(vertex(4), vertex(5), "e"),
                    ArrowSpec::named(vertex(4), vertex(5), "f"),
                ],
            )
            .unwrap(),
        )
        .unwrap();
        let generators = parent.generators();
        let relations = vec![
            generators[5]
                .multiply(&generators[6])
                .unwrap()
                .sub(&generators[7].multiply(&generators[8]).unwrap())
                .unwrap(),
            generators[6].multiply(&generators[10]).unwrap(),
            generators[8].multiply(&generators[9]).unwrap(),
        ];
        let ideal = PathAlgebraIdeal::new(parent, relations)
            .unwrap()
            .with_high_level_groebner_basis()
            .unwrap();
        let algebra = QuotientPathAlgebra::by_ideal(&ideal).unwrap();
        let quotient_generators = algebra.generators().unwrap();
        let resolution = projective_resolution_fp_path_algebra_module(
            &algebra,
            &ideal,
            &[
                vec![quotient_generators[5].clone()],
                vec![quotient_generators[7].clone()],
            ],
            5,
        )
        .unwrap();

        assert_eq!(resolution.finite_at_projective(), Some(4));
        let homology = resolution.semisimple_homology().unwrap();
        assert_eq!(
            homology.projective_vertex_multiplicities(),
            &[
                vec![1, 0, 0, 0, 0],
                vec![0, 1, 1, 0, 0],
                vec![0, 0, 0, 1, 1],
                vec![0, 0, 0, 0, 1],
            ]
        );
        assert_eq!(
            homology.differential_vertex_ranks(),
            &[
                vec![0, 0, 0, 0, 0],
                vec![0, 0, 0, 0, 0],
                vec![0, 0, 0, 0, 0],
                vec![0, 0, 0, 0, 1],
                vec![0, 0, 0, 0, 0],
            ]
        );
        assert_eq!(
            homology.minimal_projective_betti_multiplicities(),
            &[
                vec![1, 0, 0, 0, 0],
                vec![0, 1, 1, 0, 0],
                vec![0, 0, 0, 1, 0],
                vec![0, 0, 0, 0, 0],
            ]
        );
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
            projective_resolution_fp_path_algebra_module_without_term(
                &algebra,
                &ideal,
                &[vec![x.clone()]],
            ),
            Err(ProjectiveResolutionError::UninstalledThreeArgumentConstructor)
        ));

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

        let mut resolution =
            projective_resolution_fp_path_algebra_module(&algebra, &ideal, &[vec![x]], 1).unwrap();
        assert!(!resolution.r_projectives_vertex_list().is_bound());
        assert!(!resolution.projectives_f_prime_list().is_bound());
        assert!(matches!(
            find_next_r_projective(&resolution, &parent, 1),
            Err(ProjectiveResolutionError::UninstalledFindNextRProjective)
        ));
        assert!(matches!(
            find_next_syzygy(&resolution, 1),
            Err(ProjectiveResolutionError::UninstalledFindNextSyzygy)
        ));
        resolution.set_r_projectives_vertex_list(vec![vec![vertex(1)]]);
        resolution.set_projectives_f_prime_list(vec![Vec::new()]);
        assert!(resolution.r_projectives_vertex_list().is_bound());
        assert!(resolution.projectives_f_prime_list().is_bound());
        assert_eq!(
            INFO_PROJECTIVE_RESOLUTION_FP_PATH_ALGEBRA_MODULE,
            "InfoProjectiveResolutionFpPathAlgebraModule"
        );
        assert_eq!(
            INFO_PROJECTIVE_RESOLUTION_FP_PATH_ALGEBRA_MODULE_DEFAULT_LEVEL,
            0
        );
    }
}
