//! QPA's finitely presented path-algebra modules.
//!
//! This is the representation constructed by
//! `NaturalHomomorphismBySubAlgebraModule` and used by both overloads of
//! `VertexProjectivePresentation`. Elements retain QPA's path-vector external
//! representation and are normalized by the right module Groebner basis of
//! the entered submodule.

use std::{cmp::Ordering, fmt, sync::Arc};

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    Path, PathAlgebra, PathAlgebraElement, PathEndpoint, ProjectivePathModuleError,
    RightPathAlgebraModuleGroebnerBasis, RightProjectivePathAlgebraModule,
    RightProjectivePathAlgebraSubmodule, VertexIndex, right_projective_path_algebra_module,
};

#[derive(Clone, Debug)]
struct RawFpPathAlgebraVector<F: Field> {
    components: Vec<PathAlgebraElement<F>>,
    leading_position: usize,
}

struct RightFpPathAlgebraModuleData<F: Field> {
    parent: RightProjectivePathAlgebraModule<F>,
    bottom: RightProjectivePathAlgebraSubmodule<F>,
    groebner_basis: RightPathAlgebraModuleGroebnerBasis<F>,
    basis: Box<[RawFpPathAlgebraVector<F>]>,
}

/// QPA's quotient `P/N` represented as an `IsFpPathAlgebraModule`.
///
/// Every construction creates a fresh family. Clones retain that exact family
/// identity, as GAP's `IsIdenticalObj` filters require.
#[derive(Clone)]
pub struct RightFpPathAlgebraModule<F: Field>(Arc<RightFpPathAlgebraModuleData<F>>);

impl<F: Field> RightFpPathAlgebraModule<F> {
    /// Return the acting path algebra.
    #[must_use]
    pub fn algebra(&self) -> &PathAlgebra<F> {
        self.0.parent.algebra()
    }

    /// Return QPA's projective preimage module.
    #[must_use]
    pub fn parent_projective(&self) -> &RightProjectivePathAlgebraModule<F> {
        &self.0.parent
    }

    /// Return QPA's `bottomModule` before its pure-path-algebra lift.
    #[must_use]
    pub fn bottom_submodule(&self) -> &RightProjectivePathAlgebraSubmodule<F> {
        &self.0.bottom
    }

    /// Return the right module Groebner basis used by the normal-form function.
    #[must_use]
    pub fn groebner_basis(&self) -> &RightPathAlgebraModuleGroebnerBasis<F> {
        &self.0.groebner_basis
    }

    /// Return whether two handles denote the identical QPA module/vector family.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Return QPA's unchecked family zero with position `rank + 1`.
    #[must_use]
    pub fn zero(&self) -> FpPathAlgebraVector<F> {
        self.raw_vector(RawFpPathAlgebraVector {
            components: vec![
                self.algebra().zero();
                self.parent_projective().vertex_indices().len()
            ],
            leading_position: self.parent_projective().vertex_indices().len() + 1,
        })
    }

    /// Return QPA's unchecked coordinate vertex generators.
    ///
    /// # Errors
    ///
    /// Propagates the parent module's vertex embedding boundary.
    pub fn generators(&self) -> Result<Vec<FpPathAlgebraVector<F>>, FpPathModuleError> {
        self.parent_projective()
            .generators()?
            .into_iter()
            .map(|vector| {
                Ok(self.raw_vector(RawFpPathAlgebraVector {
                    components: vector.components().to_vec(),
                    leading_position: vector.leading_position(),
                }))
            })
            .collect()
    }

    /// Return the normal monomial basis installed by QPA.
    #[must_use]
    pub fn basis(&self) -> Vec<FpPathAlgebraVector<F>> {
        self.0
            .basis
            .iter()
            .cloned()
            .map(|raw| self.raw_vector(raw))
            .collect()
    }

    /// Return QPA's installed finite dimension.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.0.basis.len()
    }

    /// Return QPA's installed `IsFpPathAlgebraModule` filter.
    #[must_use]
    pub const fn is_fp_path_algebra_module(&self) -> bool {
        true
    }

    /// Return QPA's installed `IsPathAlgebraModule` filter.
    #[must_use]
    pub const fn is_path_algebra_module(&self) -> bool {
        true
    }

    /// Return QPA's installed whole-family filter.
    #[must_use]
    pub const fn is_whole_family(&self) -> bool {
        true
    }

    /// Translate the checked `PathAlgebraVector`/`Vectorize` route.
    ///
    /// QPA first checks the projective coordinate conditions and chooses a tip,
    /// then calls the family normal-form function.
    ///
    /// # Errors
    ///
    /// Propagates parent-vector validation and Groebner reduction failures.
    pub fn vector(
        &self,
        components: Vec<PathAlgebraElement<F>>,
    ) -> Result<FpPathAlgebraVector<F>, FpPathModuleError> {
        let preimage = self.parent_projective().vector(components)?;
        self.normal_form(&preimage)
    }

    /// Translate QPA's nonchecking constructor.
    ///
    /// With `normal = true`, QPA trusts all supplied data. Otherwise it invokes
    /// the same normal-form function as the checked constructor.
    ///
    /// # Errors
    ///
    /// A requested reduction propagates malformed-position and Groebner errors.
    pub fn vector_nc(
        &self,
        components: Vec<PathAlgebraElement<F>>,
        leading_position: usize,
        normal: bool,
    ) -> Result<FpPathAlgebraVector<F>, FpPathModuleError> {
        if normal {
            return Ok(self.raw_vector(RawFpPathAlgebraVector {
                components,
                leading_position,
            }));
        }
        self.normal_form(
            &self
                .parent_projective()
                .vector_nc(components, leading_position, false),
        )
    }

    fn normal_form(
        &self,
        preimage: &crate::PathAlgebraVector<F>,
    ) -> Result<FpPathAlgebraVector<F>, FpPathModuleError> {
        let reduced = self.groebner_basis().completely_reduce(preimage)?;
        Ok(self.raw_vector(RawFpPathAlgebraVector {
            components: reduced.components().to_vec(),
            leading_position: reduced.leading_position(),
        }))
    }

    fn raw_vector(&self, raw: RawFpPathAlgebraVector<F>) -> FpPathAlgebraVector<F> {
        FpPathAlgebraVector {
            module: self.clone(),
            components: raw.components,
            leading_position: raw.leading_position,
        }
    }
}

impl<F: Field> fmt::Debug for RightFpPathAlgebraModule<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("RightFpPathAlgebraModule")
            .field("parent_projective", self.parent_projective())
            .field("bottom_submodule", self.bottom_submodule())
            .field("dimension", &self.dimension())
            .finish()
    }
}

impl<F: Field> PartialEq for RightFpPathAlgebraModule<F> {
    fn eq(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}

impl<F: Field> Eq for RightFpPathAlgebraModule<F> {}

/// A normal-form vector in QPA's finitely presented module family.
#[derive(Clone, Debug)]
pub struct FpPathAlgebraVector<F: Field> {
    module: RightFpPathAlgebraModule<F>,
    components: Vec<PathAlgebraElement<F>>,
    leading_position: usize,
}

impl<F: Field> FpPathAlgebraVector<F> {
    /// Return the identical owning finitely presented module.
    #[must_use]
    pub const fn module(&self) -> &RightFpPathAlgebraModule<F> {
        &self.module
    }

    /// Return QPA's stored path-algebra coordinates.
    #[must_use]
    pub fn components(&self) -> &[PathAlgebraElement<F>] {
        &self.components
    }

    /// Return QPA's cached one-based leading position.
    #[must_use]
    pub const fn leading_position(&self) -> usize {
        self.leading_position
    }

    /// Return the default type's installed normal-form filter.
    #[must_use]
    pub const fn is_normal_form(&self) -> bool {
        true
    }

    /// Return whether every coordinate vanishes.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.components.iter().all(PathAlgebraElement::is_zero)
    }

    /// Add identical-family normal forms through QPA's `PathAlgebraVectorNC`.
    ///
    /// # Errors
    ///
    /// Rejects distinct quotient families and propagates coordinate arithmetic.
    pub fn add(&self, right: &Self) -> Result<Self, FpPathModuleError> {
        self.ensure_same_family(right)?;
        self.normal_combination(right, PathAlgebraElement::add)
    }

    /// Subtract identical-family normal forms through QPA's unchecked route.
    ///
    /// # Errors
    ///
    /// Rejects distinct quotient families and propagates coordinate arithmetic.
    pub fn sub(&self, right: &Self) -> Result<Self, FpPathModuleError> {
        self.ensure_same_family(right)?;
        self.normal_combination(right, PathAlgebraElement::sub)
    }

    /// Return QPA's unchecked additive inverse.
    ///
    /// # Errors
    ///
    /// Propagates coordinate negation.
    pub fn neg(&self) -> Result<Self, FpPathModuleError> {
        self.module.vector_nc(
            self.components
                .iter()
                .map(PathAlgebraElement::neg)
                .collect::<Result<Vec<_>, _>>()?,
            self.leading_position,
            true,
        )
    }

    /// Return QPA's unchecked left scalar multiple.
    ///
    /// # Errors
    ///
    /// Propagates coordinate scaling.
    pub fn scale_left(&self, scalar: &F::Element) -> Result<Self, FpPathModuleError> {
        if self.module.algebra().field().is_zero(scalar) {
            return Ok(self.module.zero());
        }
        self.module.vector_nc(
            self.components
                .iter()
                .map(|component| component.scale_left(scalar))
                .collect::<Result<Vec<_>, _>>()?,
            self.leading_position,
            true,
        )
    }

    /// Return QPA's identical right scalar action.
    ///
    /// # Errors
    ///
    /// Propagates coordinate scaling.
    pub fn scale_right(&self, scalar: &F::Element) -> Result<Self, FpPathModuleError> {
        self.scale_left(scalar)
    }

    /// Apply the checked coordinatewise right path-algebra action and normalize.
    ///
    /// # Errors
    ///
    /// Rejects a foreign action algebra and propagates exact multiplication.
    pub fn act(&self, element: &PathAlgebraElement<F>) -> Result<Self, FpPathModuleError> {
        if !element.algebra().is_identical(self.module.algebra()) {
            return Err(FpPathModuleError::ForeignActionAlgebra);
        }
        self.module.vector(
            self.components
                .iter()
                .map(|component| component.multiply(element))
                .collect::<Result<Vec<_>, _>>()?,
        )
    }

    /// Return QPA's sparse unchecked leading term.
    ///
    /// # Errors
    ///
    /// The trusted NC constructor itself cannot fail in this branch.
    pub fn leading_term(&self) -> Result<Self, FpPathModuleError> {
        let Some(position) = self.leading_position.checked_sub(1) else {
            return Ok(self.clone());
        };
        let Some(component) = self.components.get(position) else {
            return Ok(self.clone());
        };
        let mut components = vec![self.module.algebra().zero(); self.components.len()];
        components[position] = component.leading_term();
        self.module
            .vector_nc(components, self.leading_position, true)
    }

    /// Return QPA's vector-tip coefficient.
    ///
    /// # Errors
    ///
    /// QPA indexes coordinate one for zero and fails for rank zero.
    pub fn leading_coefficient(&self) -> Result<F::Element, FpPathModuleError> {
        let position = if self.leading_position <= self.components.len() {
            self.leading_position - 1
        } else {
            0
        };
        self.components
            .get(position)
            .map(PathAlgebraElement::leading_coefficient)
            .ok_or(FpPathModuleError::EmptyLeadingComponent)
    }

    /// Return the coordinate containing QPA's vector tip.
    ///
    /// # Errors
    ///
    /// QPA indexes coordinate one for zero and fails for rank zero.
    pub fn leading_component(&self) -> Result<&PathAlgebraElement<F>, FpPathModuleError> {
        let position = if self.is_zero() {
            0
        } else {
            self.leading_position - 1
        };
        self.components
            .get(position)
            .ok_or(FpPathModuleError::EmptyLeadingComponent)
    }

    /// Translate QPA's vector target; `None` is GAP's `fail`.
    #[must_use]
    pub fn target_vertex(&self) -> Option<Path> {
        if !self.is_right_uniform() {
            return None;
        }
        if self.is_zero() {
            return Some(self.module.algebra().quiver().zero());
        }
        match self.leading_path().ok()?.target() {
            PathEndpoint::Vertex(vertex) => Some(vertex.path()),
            PathEndpoint::Zero(path) => Some(path),
        }
    }

    /// Translate QPA's surprising right-uniform-gated vector source.
    #[must_use]
    pub fn source_vertex(&self) -> Option<Path> {
        if !self.is_right_uniform() {
            return None;
        }
        if self.is_zero() {
            return Some(self.module.algebra().quiver().zero());
        }
        match self.leading_path().ok()?.source() {
            PathEndpoint::Vertex(vertex) => Some(vertex.path()),
            PathEndpoint::Zero(path) => Some(path),
        }
    }

    /// Return QPA's left-uniformity predicate.
    #[must_use]
    pub fn is_left_uniform(&self) -> bool {
        self.is_uniform_at(Path::source, PathAlgebraElement::is_left_uniform)
    }

    /// Return QPA's right-uniformity predicate.
    #[must_use]
    pub fn is_right_uniform(&self) -> bool {
        self.is_uniform_at(Path::target, PathAlgebraElement::is_right_uniform)
    }

    /// Return whether the vector is both left and right uniform.
    #[must_use]
    pub fn is_uniform(&self) -> bool {
        self.is_left_uniform() && self.is_right_uniform()
    }

    /// Compare normal forms by QPA's tip/descending-position order.
    ///
    /// # Errors
    ///
    /// Rejects distinct families, zero tips, and ordering failures.
    pub fn compare(&self, right: &Self) -> Result<Ordering, FpPathModuleError> {
        self.ensure_same_family(right)?;
        let path_order = self
            .module
            .algebra()
            .ordering()
            .compare(&self.leading_path()?, &right.leading_path()?)?;
        if path_order.is_eq() {
            Ok(right.leading_position.cmp(&self.leading_position))
        } else {
            Ok(path_order)
        }
    }

    /// Translate QPA's vector left-divisibility predicate.
    ///
    /// # Errors
    ///
    /// Rejects distinct families and zero or non-walk tips.
    pub fn is_left_divisible(&self, divisor: &Self) -> Result<bool, FpPathModuleError> {
        self.ensure_same_family(divisor)?;
        if self.is_zero() || divisor.is_zero() {
            return Err(FpPathModuleError::ZeroVectorDivision);
        }
        let dividend_path = self.leading_path()?;
        let divisor_path = divisor.leading_path()?;
        let dividend_walk = dividend_path
            .walk()
            .ok_or(FpPathModuleError::PathHasNoArrowWalk)?;
        let divisor_walk = divisor_path
            .walk()
            .ok_or(FpPathModuleError::PathHasNoArrowWalk)?;
        Ok(self.leading_position == divisor.leading_position
            && !dividend_walk.is_empty()
            && dividend_walk.starts_with(divisor_walk))
    }

    fn normal_combination(
        &self,
        right: &Self,
        operation: impl Fn(
            &PathAlgebraElement<F>,
            &PathAlgebraElement<F>,
        ) -> Result<PathAlgebraElement<F>, crate::PathAlgebraError>,
    ) -> Result<Self, FpPathModuleError> {
        let components = self
            .components
            .iter()
            .zip(&right.components)
            .map(|(left, right)| operation(left, right))
            .collect::<Result<Vec<_>, _>>()?;
        let positioned = self.module.parent_projective().vector(components)?;
        self.module.vector_nc(
            positioned.components().to_vec(),
            positioned.leading_position(),
            true,
        )
    }

    fn ensure_same_family(&self, right: &Self) -> Result<(), FpPathModuleError> {
        if self.module.is_identical(&right.module) {
            Ok(())
        } else {
            Err(FpPathModuleError::DifferentVectorFamilies)
        }
    }

    fn leading_path(&self) -> Result<Path, FpPathModuleError> {
        if self.is_zero() {
            return Err(FpPathModuleError::ZeroVectorDivision);
        }
        self.components
            .get(self.leading_position - 1)
            .map(PathAlgebraElement::leading_monomial)
            .ok_or(FpPathModuleError::InvalidLeadingPosition {
                position: self.leading_position,
                length: self.components.len(),
            })
    }

    fn is_uniform_at(
        &self,
        endpoint: impl Fn(&Path) -> PathEndpoint,
        component_is_uniform: impl Fn(&PathAlgebraElement<F>) -> bool,
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
                component
                    .terms()
                    .first()
                    .and_then(|term| match endpoint(term.path()) {
                        PathEndpoint::Vertex(vertex) => Some(vertex),
                        PathEndpoint::Zero(_) => None,
                    })
            });
        let Some(first) = endpoints.next() else {
            return true;
        };
        endpoints.all(|candidate| candidate == first)
    }
}

impl<F: Field> PartialEq for FpPathAlgebraVector<F> {
    fn eq(&self, other: &Self) -> bool {
        self.module.is_identical(&other.module) && self.components == other.components
    }
}

impl<F: Field> Eq for FpPathAlgebraVector<F> {}

/// QPA's surjective natural map from `P` to the finitely presented `P/N`.
#[derive(Clone, Debug)]
pub struct NaturalSubmoduleHomomorphism<F: Field> {
    source: RightProjectivePathAlgebraModule<F>,
    range: RightFpPathAlgebraModule<F>,
}

impl<F: Field> NaturalSubmoduleHomomorphism<F> {
    /// Return the source vertex-projective module.
    #[must_use]
    pub const fn source(&self) -> &RightProjectivePathAlgebraModule<F> {
        &self.source
    }

    /// Return the quotient-range module.
    #[must_use]
    pub const fn range(&self) -> &RightFpPathAlgebraModule<F> {
        &self.range
    }

    /// Return QPA's installed surjectivity property.
    #[must_use]
    pub const fn is_surjective(&self) -> bool {
        true
    }

    /// Return the entered bottom module, which is the natural map's kernel.
    #[must_use]
    pub fn kernel(&self) -> &RightProjectivePathAlgebraSubmodule<F> {
        self.range.bottom_submodule()
    }

    /// Map a source vector by retaining its external coordinates and reducing.
    ///
    /// # Errors
    ///
    /// Rejects a vector from another projective family and propagates reduction.
    pub fn image(
        &self,
        vector: &crate::PathAlgebraVector<F>,
    ) -> Result<FpPathAlgebraVector<F>, FpPathModuleError> {
        if !vector.module().is_identical(&self.source) {
            return Err(FpPathModuleError::ForeignHomomorphismSourceVector);
        }
        self.range.vector(vector.components().to_vec())
    }
}

/// Translate QPA's `NaturalHomomorphismBySubAlgebraModule(P, N)`.
///
/// # Errors
///
/// Requires `N` to use `P`'s identical vector family and requires a finite
/// canonical basis, just as QPA's basis-installation route does.
pub fn natural_homomorphism_by_projective_submodule<F: Field>(
    parent: &RightProjectivePathAlgebraModule<F>,
    submodule: &RightProjectivePathAlgebraSubmodule<F>,
) -> Result<NaturalSubmoduleHomomorphism<F>, FpPathModuleError> {
    if !submodule.parent().is_identical(parent) {
        return Err(FpPathModuleError::ForeignBottomSubmodule);
    }
    let groebner_basis = submodule.right_groebner_basis()?;
    let canonical_basis = parent
        .canonical_basis()?
        .ok_or(FpPathModuleError::CanonicalBasisUnavailable)?;
    let mut basis = Vec::new();
    for vector in canonical_basis {
        let normal = groebner_basis.completely_reduce(&vector)?;
        if vector.components() == normal.components() {
            basis.push(RawFpPathAlgebraVector {
                components: vector.components().to_vec(),
                leading_position: vector.leading_position(),
            });
        }
    }
    let range = RightFpPathAlgebraModule(Arc::new(RightFpPathAlgebraModuleData {
        parent: parent.clone(),
        bottom: submodule.clone(),
        groebner_basis,
        basis: basis.into_boxed_slice(),
    }));
    Ok(NaturalSubmoduleHomomorphism {
        source: parent.clone(),
        range,
    })
}

/// Translate `VertexProjectivePresentation(A, map)` with inferred vertices.
///
/// # Errors
///
/// Preserves QPA's table-shape, left-uniformity, source compatibility, and
/// nonzero-column checks before constructing `P/N` through the natural map.
pub fn vertex_projective_presentation<F: Field>(
    algebra: &PathAlgebra<F>,
    map: &[Vec<PathAlgebraElement<F>>],
) -> Result<RightFpPathAlgebraModule<F>, FpPathModuleError> {
    let first = map.first().ok_or(FpPathModuleError::EmptyPresentationMap)?;
    ensure_rectangular_map(algebra, map, first.len())?;
    ensure_left_uniform_map(map)?;

    let mut vertices = first.iter().map(source_vertex_index).collect::<Vec<_>>();
    for row in map {
        for (position, component) in row.iter().enumerate() {
            let candidate = source_vertex_index(component);
            match (vertices[position], candidate) {
                (None, Some(vertex)) => vertices[position] = Some(vertex),
                (Some(expected), Some(actual)) if expected != actual => {
                    return Err(FpPathModuleError::MismatchedStartingVertex {
                        position: position + 1,
                    });
                }
                _ => {}
            }
        }
    }
    let vertices = vertices
        .into_iter()
        .enumerate()
        .map(|(position, vertex)| {
            vertex.ok_or(FpPathModuleError::AllZeroPresentationColumn {
                position: position + 1,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    construct_vertex_presentation(algebra, map, &vertices)
}

/// Translate `VertexProjectivePresentation(A, map, verts)`.
///
/// # Errors
///
/// Requires coefficient-one embedded vertices and checks every nonzero map
/// component against the corresponding explicit source vertex.
pub fn vertex_projective_presentation_with_vertices<F: Field>(
    algebra: &PathAlgebra<F>,
    map: &[Vec<PathAlgebraElement<F>>],
    vertices: &[PathAlgebraElement<F>],
) -> Result<RightFpPathAlgebraModule<F>, FpPathModuleError> {
    map.first().ok_or(FpPathModuleError::EmptyPresentationMap)?;
    ensure_rectangular_map(algebra, map, vertices.len())?;
    ensure_left_uniform_map(map)?;
    let parent = right_projective_path_algebra_module(algebra, vertices)?;
    for row in map {
        for (position, component) in row.iter().enumerate() {
            if let Some(actual) = source_vertex_index(component)
                && actual != parent.vertex_indices()[position]
            {
                return Err(FpPathModuleError::MismatchedStartingVertex {
                    position: position + 1,
                });
            }
        }
    }
    quotient_by_map(&parent, map)
}

fn construct_vertex_presentation<F: Field>(
    algebra: &PathAlgebra<F>,
    map: &[Vec<PathAlgebraElement<F>>],
    vertices: &[VertexIndex],
) -> Result<RightFpPathAlgebraModule<F>, FpPathModuleError> {
    let embedded = vertices
        .iter()
        .map(|vertex| {
            algebra
                .quiver()
                .vertex(*vertex)
                .ok_or(ProjectivePathModuleError::VertexOutOfRange { vertex: *vertex })
                .and_then(|vertex| algebra.element(vertex.path()).map_err(Into::into))
        })
        .collect::<Result<Vec<_>, ProjectivePathModuleError>>()?;
    let parent = right_projective_path_algebra_module(algebra, &embedded)?;
    quotient_by_map(&parent, map)
}

fn quotient_by_map<F: Field>(
    parent: &RightProjectivePathAlgebraModule<F>,
    map: &[Vec<PathAlgebraElement<F>>],
) -> Result<RightFpPathAlgebraModule<F>, FpPathModuleError> {
    let relations = map
        .iter()
        .map(|row| parent.vector(row.clone()))
        .collect::<Result<Vec<_>, _>>()?;
    let submodule = parent.submodule(relations)?;
    Ok(
        natural_homomorphism_by_projective_submodule(parent, &submodule)?
            .range()
            .clone(),
    )
}

fn ensure_rectangular_map<F: Field>(
    algebra: &PathAlgebra<F>,
    map: &[Vec<PathAlgebraElement<F>>],
    columns: usize,
) -> Result<(), FpPathModuleError> {
    for (row, entries) in map.iter().enumerate() {
        if entries.len() != columns {
            return Err(FpPathModuleError::PresentationRowLength {
                row: row + 1,
                actual: entries.len(),
                expected: columns,
            });
        }
        if entries
            .iter()
            .any(|entry| !entry.algebra().is_identical(algebra))
        {
            return Err(FpPathModuleError::ForeignPresentationEntry { row: row + 1 });
        }
    }
    Ok(())
}

fn ensure_left_uniform_map<F: Field>(
    map: &[Vec<PathAlgebraElement<F>>],
) -> Result<(), FpPathModuleError> {
    for (row, entries) in map.iter().enumerate() {
        if entries.iter().any(|entry| !entry.is_left_uniform()) {
            return Err(FpPathModuleError::NonuniformPresentationRow { row: row + 1 });
        }
    }
    Ok(())
}

fn source_vertex_index<F: Field>(element: &PathAlgebraElement<F>) -> Option<VertexIndex> {
    if element.is_zero() {
        return None;
    }
    match element.leading_monomial().source() {
        PathEndpoint::Vertex(vertex) => Some(vertex.index()),
        PathEndpoint::Zero(_) => None,
    }
}

/// Failure in QPA's finitely presented path-module route.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum FpPathModuleError {
    /// The natural map requires an identical ambient projective family.
    #[error("the bottom submodule does not belong to the entered projective module")]
    ForeignBottomSubmodule,
    /// A cyclic path algebra has no canonical module basis to quotient.
    #[error("the projective module has no available canonical basis")]
    CanonicalBasisUnavailable,
    /// Arithmetic requires one identical finitely presented vector family.
    #[error("the finitely presented vectors belong to different module families")]
    DifferentVectorFamilies,
    /// Right action requires the identical parent path algebra.
    #[error("the finitely presented vector action belongs to another path algebra")]
    ForeignActionAlgebra,
    /// Natural-map images require a vector in the map's exact source family.
    #[error("the vector does not belong to the natural homomorphism source")]
    ForeignHomomorphismSourceVector,
    /// QPA indexes coordinate one for zero and fails for rank zero.
    #[error("an empty finitely presented vector has no leading component")]
    EmptyLeadingComponent,
    /// An unchecked vector retained an invalid one-based position.
    #[error("leading position {position} is outside a vector of length {length}")]
    InvalidLeadingPosition {
        /// Stored one-based position.
        position: usize,
        /// Coordinate count.
        length: usize,
    },
    /// QPA's divisibility operation has no zero-vector tip.
    #[error("finitely presented path-vector divisibility does not accept zero inputs")]
    ZeroVectorDivision,
    /// `WalkOfPath` is unavailable for the entered tip.
    #[error("a finitely presented path-vector tip has no arrow walk")]
    PathHasNoArrowWalk,
    /// Empty GAP lists do not satisfy the installed table/list dispatch.
    #[error("no vertex-projective presentation method applies to an empty map")]
    EmptyPresentationMap,
    /// The GAP table/explicit-length boundary requires a rectangular map.
    #[error("presentation row {row} has length {actual}, expected {expected}")]
    PresentationRowLength {
        /// One-based row.
        row: usize,
        /// Actual row length.
        actual: usize,
        /// Required row length.
        expected: usize,
    },
    /// Every map entry must belong to the identical acting algebra.
    #[error("presentation row {row} contains an entry from another algebra")]
    ForeignPresentationEntry {
        /// One-based row.
        row: usize,
    },
    /// QPA requires every entry in each row to be left uniform.
    #[error("presentation row {row} contains a non-left-uniform entry")]
    NonuniformPresentationRow {
        /// One-based row.
        row: usize,
    },
    /// Nonzero components in a column must have one common source vertex.
    #[error("presentation column {position} contains mismatched starting vertices")]
    MismatchedStartingVertex {
        /// One-based column.
        position: usize,
    },
    /// Inferred presentations require a nonzero entry in every column.
    #[error("presentation column {position} contains only zeroes")]
    AllZeroPresentationColumn {
        /// One-based column.
        position: usize,
    },
    /// The underlying projective-vector implementation failed.
    #[error(transparent)]
    Projective(#[from] ProjectivePathModuleError),
    /// Parent path-algebra arithmetic failed directly.
    #[error(transparent)]
    PathAlgebra(#[from] crate::PathAlgebraError),
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

    use crate::{ArrowSpec, PathAlgebra, PathAlgebraElement, Quiver, VertexIndex};

    use super::{
        FpPathAlgebraVector, FpPathModuleError, natural_homomorphism_by_projective_submodule,
        vertex_projective_presentation, vertex_projective_presentation_with_vertices,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    type PresentationData<F> = (
        PathAlgebra<F>,
        Vec<PathAlgebraElement<F>>,
        Vec<Vec<PathAlgebraElement<F>>>,
    );

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

    fn vector_json<F: Field>(
        vector: &FpPathAlgebraVector<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            vector
                .components()
                .iter()
                .map(|component| element_json(component, coefficient))
                .collect(),
        )
    }

    fn vector_list_json<F: Field>(
        vectors: &[FpPathAlgebraVector<F>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            vectors
                .iter()
                .map(|vector| vector_json(vector, coefficient))
                .collect(),
        )
    }

    fn presentation_data<F: Field>(field: F) -> PresentationData<F> {
        let algebra = PathAlgebra::new(
            field,
            Quiver::from_indexed(
                3,
                vec![
                    ArrowSpec::named(vertex(1), vertex(2), "a"),
                    ArrowSpec::named(vertex(2), vertex(3), "b"),
                    ArrowSpec::named(vertex(1), vertex(3), "c"),
                ],
            )
            .unwrap(),
        )
        .unwrap();
        let generators = algebra.generators();
        let map = vec![
            vec![generators[3].clone(), generators[1].neg().unwrap()],
            vec![generators[5].clone(), generators[4].neg().unwrap()],
        ];
        (algebra, generators, map)
    }

    fn presentation_case<F: Field>(
        field: F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let (algebra, generators, map) = presentation_data(field);
        let module = vertex_projective_presentation(&algebra, &map).unwrap();
        let explicit = vertex_projective_presentation_with_vertices(
            &algebra,
            &map,
            &[generators[0].clone(), generators[1].clone()],
        )
        .unwrap();
        let module_generators = module.generators().unwrap();
        let relation = module.vector(map[0].clone()).unwrap();
        let probe_a = module
            .vector(vec![generators[3].clone(), algebra.zero()])
            .unwrap();
        let probe_c = module
            .vector(vec![generators[5].clone(), algebra.zero()])
            .unwrap();
        let probe_b = module
            .vector(vec![algebra.zero(), generators[4].clone()])
            .unwrap();
        let action_a = module_generators[0].act(&generators[3]).unwrap();
        let action_b = module_generators[1].act(&generators[4]).unwrap();
        let generator_sum = module_generators[0].add(&module_generators[1]).unwrap();
        let negative_b = probe_b.neg().unwrap();
        let two = algebra
            .field()
            .add(&algebra.field().one(), &algebra.field().one());
        let scaled_b = probe_b.scale_left(&two).unwrap();
        let cancelled_b = probe_b.add(&negative_b).unwrap();
        json!({
            "dimension": module.dimension(),
            "properties": [
                module.is_fp_path_algebra_module(),
                module.is_path_algebra_module(),
                module.is_whole_family(),
            ],
            "generators": vector_list_json(&module_generators, coefficient),
            "basis": vector_list_json(&module.basis(), coefficient),
            "zero": vector_json(&module.zero(), coefficient),
            "positions": [
                module.zero().leading_position(),
                relation.leading_position(),
                module_generators[0].leading_position(),
                module_generators[1].leading_position(),
            ],
            "normal_forms": vector_list_json(
                &[relation, probe_a, probe_c, probe_b],
                coefficient,
            ),
            "actions": vector_list_json(&[action_a, action_b], coefficient),
            "arithmetic": vector_list_json(
                &[generator_sum, negative_b, scaled_b, cancelled_b],
                coefficient,
            ),
            "explicit": {
                "dimension": explicit.dimension(),
                "basis": vector_list_json(&explicit.basis(), coefficient),
            },
        })
    }

    fn presentation_error_success<F: Field>(
        algebra: &PathAlgebra<F>,
        generators: &[PathAlgebraElement<F>],
        map: &[Vec<PathAlgebraElement<F>>],
    ) -> Value {
        let nonuniform = generators[3].add(&generators[4]).unwrap();
        json!({
            "success": [
                vertex_projective_presentation(algebra, map).is_ok(),
                vertex_projective_presentation_with_vertices(
                    algebra,
                    map,
                    &[generators[0].clone(), generators[1].clone()],
                ).is_ok(),
            ],
            "failure_success": [
                vertex_projective_presentation(algebra, &[]).is_ok(),
                vertex_projective_presentation(
                    algebra,
                    &[map[0].clone(), vec![generators[5].clone()]],
                ).is_ok(),
                vertex_projective_presentation(
                    algebra,
                    &[vec![nonuniform, generators[1].clone()]],
                ).is_ok(),
                vertex_projective_presentation(
                    algebra,
                    &[
                        map[0].clone(),
                        vec![generators[4].clone(), algebra.zero()],
                    ],
                ).is_ok(),
                vertex_projective_presentation(
                    algebra,
                    &[vec![algebra.zero(), generators[1].clone()]],
                ).is_ok(),
                vertex_projective_presentation_with_vertices(
                    algebra,
                    &[],
                    &[generators[0].clone(), generators[1].clone()],
                ).is_ok(),
                vertex_projective_presentation_with_vertices(
                    algebra,
                    &[vec![generators[3].clone()]],
                    &[generators[0].clone(), generators[1].clone()],
                ).is_ok(),
                vertex_projective_presentation_with_vertices(
                    algebra,
                    map,
                    &[generators[3].clone(), generators[1].clone()],
                ).is_ok(),
                vertex_projective_presentation_with_vertices(
                    algebra,
                    map,
                    &[generators[1].clone(), generators[1].clone()],
                ).is_ok(),
            ],
        })
    }

    #[test]
    fn pinned_qpa_vertex_presentation_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/projpathalgmodule_vertex_presentation.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "projpathalgmodule_vertex_presentation",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654",
            },
            "binary": presentation_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": presentation_case(field5, |entry| json!(
                field5_coefficients
                    .to_integer(entry)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
            )),
            "rational": presentation_case(Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn fp_path_module_boundaries_are_typed() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/projpathalgmodule_vertex_presentation_errors.schema-v1.json"
        ))
        .unwrap();
        let (algebra, generators, map) = presentation_data(Z2);
        assert_eq!(
            vertex_projective_presentation(&algebra, &[]),
            Err(FpPathModuleError::EmptyPresentationMap)
        );
        assert!(matches!(
            vertex_projective_presentation(
                &algebra,
                &[map[0].clone(), vec![generators[5].clone()]],
            ),
            Err(FpPathModuleError::PresentationRowLength { row: 2, .. })
        ));
        let nonuniform = generators[3].add(&generators[4]).unwrap();
        assert!(matches!(
            vertex_projective_presentation(&algebra, &[vec![nonuniform, generators[1].clone()]],),
            Err(FpPathModuleError::NonuniformPresentationRow { row: 1 })
        ));
        assert!(matches!(
            vertex_projective_presentation(
                &algebra,
                &[map[0].clone(), vec![generators[4].clone(), algebra.zero()],],
            ),
            Err(FpPathModuleError::MismatchedStartingVertex { position: 1 })
        ));
        assert!(matches!(
            vertex_projective_presentation(
                &algebra,
                &[vec![algebra.zero(), generators[1].clone()]],
            ),
            Err(FpPathModuleError::AllZeroPresentationColumn { position: 1 })
        ));

        let module = vertex_projective_presentation(&algebra, &map).unwrap();
        let other = vertex_projective_presentation(&algebra, &map).unwrap();
        assert!(matches!(
            module.generators().unwrap()[0].add(&other.generators().unwrap()[0]),
            Err(FpPathModuleError::DifferentVectorFamilies)
        ));
        let foreign = PathAlgebra::new(
            Z2,
            Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
        )
        .unwrap();
        assert!(matches!(
            module.generators().unwrap()[0].act(&foreign.generators()[0]),
            Err(FpPathModuleError::ForeignActionAlgebra)
        ));

        let parent = module.parent_projective();
        let natural =
            natural_homomorphism_by_projective_submodule(parent, module.bottom_submodule())
                .unwrap();
        assert!(natural.is_surjective());
        assert!(natural.kernel().is_identical(module.bottom_submodule()));
        let relation = parent.vector(map[0].clone()).unwrap();
        assert!(natural.image(&relation).unwrap().is_zero());

        let actual = presentation_error_success(&algebra, &generators, &map);
        assert_eq!(actual["success"], expected["success"]);
        assert_eq!(actual["failure_success"], expected["failure_success"]);
    }
}
