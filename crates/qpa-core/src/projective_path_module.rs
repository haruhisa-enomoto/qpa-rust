//! QPA's vertex-projective modules over a parent path algebra.
//!
//! This translates the foundational `PathAlgebraVector` representation from
//! `lib/projpathalgmodule.gi` together with the vector reductions used by
//! `lib/moduleprojres.gi`. These are polynomial path-algebra vectors, not the
//! finite-dimensional matrix representations exposed by [`crate::PathAlgebraModule`].

use std::{cmp::Ordering, fmt, sync::Arc};

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    Path, PathAlgebra, PathAlgebraElement, PathAlgebraError, PathEndpoint, PathError, VertexIndex,
};

struct RightProjectivePathAlgebraModuleData<F: Field> {
    algebra: PathAlgebra<F>,
    vertices: Box<[VertexIndex]>,
}

/// QPA's vertex-projective `RightProjectiveModule(KQ, vertices)`.
///
/// Each construction creates a fresh vector family, even when its algebra and
/// vertex list equal an earlier construction. Clones retain family identity.
#[derive(Clone)]
pub struct RightProjectivePathAlgebraModule<F: Field>(Arc<RightProjectivePathAlgebraModuleData<F>>);

impl<F: Field> RightProjectivePathAlgebraModule<F> {
    /// Return the parent path algebra.
    #[must_use]
    pub fn algebra(&self) -> &PathAlgebra<F> {
        &self.0.algebra
    }

    /// Return the summand vertices in QPA input order.
    #[must_use]
    pub fn vertex_indices(&self) -> &[VertexIndex] {
        &self.0.vertices
    }

    /// Return whether two handles denote the identical QPA module/vector family.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Return the canonical module zero.
    ///
    /// QPA constructs this through `PathAlgebraVectorNC` with leading position
    /// `Length(vertices) + 1`, which differs observably from a checked zero
    /// built by [`Self::vector`].
    #[must_use]
    pub fn zero(&self) -> PathAlgebraVector<F> {
        PathAlgebraVector::new_unchecked(
            self.clone(),
            vec![self.algebra().zero(); self.vertex_indices().len()],
            self.vertex_indices().len() + 1,
        )
    }

    /// Return QPA's coordinate vertex generators.
    ///
    /// # Errors
    ///
    /// Propagates the internal vertex embedding check.
    pub fn generators(&self) -> Result<Vec<PathAlgebraVector<F>>, ProjectivePathModuleError> {
        self.vertex_indices()
            .iter()
            .copied()
            .enumerate()
            .map(|(position, vertex)| {
                let mut components = vec![self.algebra().zero(); self.vertex_indices().len()];
                components[position] = self.algebra().element(
                    self.algebra()
                        .quiver()
                        .vertex(vertex)
                        .ok_or(ProjectivePathModuleError::VertexOutOfRange { vertex })?
                        .path(),
                )?;
                Ok(PathAlgebraVector::new_unchecked(
                    self.clone(),
                    components,
                    position + 1,
                ))
            })
            .collect()
    }

    /// Translate QPA's checked `PathAlgebraVector`/`Vectorize` constructor.
    ///
    /// The leading position uses the greatest component tip monomial and the
    /// smallest one-based coordinate in a tie. Consequently an all-zero
    /// checked vector has position one.
    ///
    /// # Errors
    ///
    /// Requires exactly one component per summand, all in the parent algebra,
    /// with QPA's left-uniform and source-vertex conditions.
    pub fn vector(
        &self,
        components: Vec<PathAlgebraElement<F>>,
    ) -> Result<PathAlgebraVector<F>, ProjectivePathModuleError> {
        if components.len() != self.vertex_indices().len() {
            return Err(ProjectivePathModuleError::ComponentCount {
                actual: components.len(),
                expected: self.vertex_indices().len(),
            });
        }
        if components
            .iter()
            .any(|component| !component.algebra().is_identical(self.algebra()))
        {
            return Err(ProjectivePathModuleError::ForeignComponentAlgebra);
        }
        for (position, component) in components.iter().enumerate().rev() {
            if !component.is_left_uniform() {
                return Err(ProjectivePathModuleError::ComponentNotLeftUniform {
                    position: position + 1,
                });
            }
            if component.is_zero() {
                continue;
            }
            let PathEndpoint::Vertex(source) = component.leading_monomial().source() else {
                return Err(ProjectivePathModuleError::ComponentSourceMismatch {
                    position: position + 1,
                });
            };
            if source.index() != self.vertex_indices()[position] {
                return Err(ProjectivePathModuleError::ComponentSourceMismatch {
                    position: position + 1,
                });
            }
        }
        let leading_position = leading_position(self.algebra(), &components)?;
        Ok(PathAlgebraVector::new_unchecked(
            self.clone(),
            components,
            leading_position,
        ))
    }

    /// Return the finite dimension when QPA's canonical path basis is finite.
    #[must_use]
    pub fn dimension(&self) -> Option<usize> {
        let basis = self.algebra().canonical_basis()?;
        Some(
            self.vertex_indices()
                .iter()
                .map(|vertex| {
                    basis
                        .iter()
                        .filter(|element| {
                            matches!(
                                element.leading_monomial().source(),
                                PathEndpoint::Vertex(source) if source.index() == *vertex
                            )
                        })
                        .count()
                })
                .sum(),
        )
    }
}

impl<F: Field> fmt::Debug for RightProjectivePathAlgebraModule<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("RightProjectivePathAlgebraModule")
            .field("algebra", self.algebra())
            .field("vertices", &self.vertex_indices())
            .finish()
    }
}

impl<F: Field> PartialEq for RightProjectivePathAlgebraModule<F> {
    fn eq(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}

impl<F: Field> Eq for RightProjectivePathAlgebraModule<F> {}

/// Construct QPA's vertex-list overload of `RightProjectiveModule`.
///
/// # Errors
///
/// Every entry must be a coefficient-one embedded vertex in the identical
/// parent path algebra.
pub fn right_projective_path_algebra_module<F: Field>(
    algebra: &PathAlgebra<F>,
    vertices: &[PathAlgebraElement<F>],
) -> Result<RightProjectivePathAlgebraModule<F>, ProjectivePathModuleError> {
    let mut indices = Vec::with_capacity(vertices.len());
    for (position, vertex) in vertices.iter().enumerate() {
        if !vertex.algebra().is_identical(algebra) {
            return Err(ProjectivePathModuleError::ForeignVertexAlgebra {
                position: position + 1,
            });
        }
        let [term] = vertex.terms() else {
            return Err(ProjectivePathModuleError::NotEmbeddedVertex {
                position: position + 1,
            });
        };
        let PathEndpoint::Vertex(source) = term.path().source() else {
            return Err(ProjectivePathModuleError::NotEmbeddedVertex {
                position: position + 1,
            });
        };
        if !term.path().is_vertex() || !algebra.field().is_one(term.coefficient()) {
            return Err(ProjectivePathModuleError::NotEmbeddedVertex {
                position: position + 1,
            });
        }
        indices.push(source.index());
    }
    Ok(RightProjectivePathAlgebraModule(Arc::new(
        RightProjectivePathAlgebraModuleData {
            algebra: algebra.clone(),
            vertices: indices.into_boxed_slice(),
        },
    )))
}

/// QPA's positional `PathAlgebraVector` external representation.
#[derive(Clone, Debug)]
pub struct PathAlgebraVector<F: Field> {
    module: RightProjectivePathAlgebraModule<F>,
    components: Vec<PathAlgebraElement<F>>,
    leading_position: usize,
}

impl<F: Field> PathAlgebraVector<F> {
    fn new_unchecked(
        module: RightProjectivePathAlgebraModule<F>,
        components: Vec<PathAlgebraElement<F>>,
        leading_position: usize,
    ) -> Self {
        Self {
            module,
            components,
            leading_position,
        }
    }

    /// Return the owning vertex-projective module/vector family.
    #[must_use]
    pub const fn module(&self) -> &RightProjectivePathAlgebraModule<F> {
        &self.module
    }

    /// Return the path-algebra coordinates.
    #[must_use]
    pub fn components(&self) -> &[PathAlgebraElement<F>] {
        &self.components
    }

    /// Return the cached one-based QPA leading position.
    #[must_use]
    pub const fn leading_position(&self) -> usize {
        self.leading_position
    }

    /// Return whether all components vanish.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.components.iter().all(PathAlgebraElement::is_zero)
    }

    /// Add two vectors in an identical QPA family and recompute the tip.
    ///
    /// # Errors
    ///
    /// Rejects distinct module families and propagates component arithmetic.
    pub fn add(&self, right: &Self) -> Result<Self, ProjectivePathModuleError> {
        self.ensure_same_family(right)?;
        self.module.vector(
            self.components
                .iter()
                .zip(&right.components)
                .map(|(left, right)| left.add(right))
                .collect::<Result<Vec<_>, _>>()?,
        )
    }

    /// Subtract two vectors in an identical QPA family and recompute the tip.
    ///
    /// # Errors
    ///
    /// Rejects distinct module families and propagates component arithmetic.
    pub fn sub(&self, right: &Self) -> Result<Self, ProjectivePathModuleError> {
        self.ensure_same_family(right)?;
        self.module.vector(
            self.components
                .iter()
                .zip(&right.components)
                .map(|(left, right)| left.sub(right))
                .collect::<Result<Vec<_>, _>>()?,
        )
    }

    /// Return the additive inverse while preserving QPA's cached position.
    ///
    /// # Errors
    ///
    /// Propagates component negation.
    pub fn neg(&self) -> Result<Self, ProjectivePathModuleError> {
        Ok(Self::new_unchecked(
            self.module.clone(),
            self.components
                .iter()
                .map(PathAlgebraElement::neg)
                .collect::<Result<Vec<_>, _>>()?,
            self.leading_position,
        ))
    }

    /// Multiply by a scalar on the left as QPA does.
    ///
    /// # Errors
    ///
    /// Propagates exact component scaling.
    pub fn scale_left(&self, scalar: &F::Element) -> Result<Self, ProjectivePathModuleError> {
        if self.module.algebra().field().is_zero(scalar) {
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

    /// Multiply every coordinate on the right by a path-algebra element.
    ///
    /// QPA routes the result through the checked constructor, including its
    /// leading-position recomputation.
    ///
    /// # Errors
    ///
    /// Rejects a foreign algebra element and propagates multiplication and
    /// checked vector construction.
    pub fn act(&self, element: &PathAlgebraElement<F>) -> Result<Self, ProjectivePathModuleError> {
        if !element.algebra().is_identical(self.module.algebra()) {
            return Err(ProjectivePathModuleError::ForeignActionAlgebra);
        }
        self.module.vector(
            self.components
                .iter()
                .map(|component| component.multiply(element))
                .collect::<Result<Vec<_>, _>>()?,
        )
    }

    /// Return QPA's sparse vector leading term.
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

    /// Return the coefficient of QPA's vector tip.
    ///
    /// # Errors
    ///
    /// QPA indexes the first coordinate for its canonical zero and therefore
    /// fails on the zero-rank projective; Rust reports that case explicitly.
    pub fn leading_coefficient(&self) -> Result<F::Element, ProjectivePathModuleError> {
        let position = if self.leading_position <= self.components.len() {
            self.leading_position - 1
        } else {
            0
        };
        self.components
            .get(position)
            .map(PathAlgebraElement::leading_coefficient)
            .ok_or(ProjectivePathModuleError::EmptyVectorLeadingComponent)
    }

    /// Return the coordinate containing QPA's vector tip.
    ///
    /// # Errors
    ///
    /// QPA indexes coordinate one for zero and fails on an empty vector.
    pub fn leading_component(&self) -> Result<&PathAlgebraElement<F>, ProjectivePathModuleError> {
        let position = if self.is_zero() {
            0
        } else {
            self.leading_position - 1
        };
        self.components
            .get(position)
            .ok_or(ProjectivePathModuleError::EmptyVectorLeadingComponent)
    }

    /// Translate QPA's vector `TargetVertex`; `None` is GAP's `fail`.
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
            return Some(self.module.algebra().quiver().zero());
        }
        match self.leading_path().ok()?.source() {
            PathEndpoint::Vertex(vertex) => Some(vertex.path()),
            PathEndpoint::Zero(path) => Some(path),
        }
    }

    /// Translate QPA's vector `IsLeftUniform`.
    #[must_use]
    pub fn is_left_uniform(&self) -> bool {
        if self.is_zero() {
            return true;
        }
        if self
            .components
            .iter()
            .any(|component| !component.is_left_uniform())
        {
            return false;
        }
        common_nonzero_endpoint(&self.components, Path::source)
    }

    /// Translate QPA's vector `IsRightUniform`.
    #[must_use]
    pub fn is_right_uniform(&self) -> bool {
        if self.is_zero() {
            return true;
        }
        if self
            .components
            .iter()
            .any(|component| !component.is_right_uniform())
        {
            return false;
        }
        common_nonzero_endpoint(&self.components, Path::target)
    }

    /// Return whether the vector is both left and right uniform.
    #[must_use]
    pub fn is_uniform(&self) -> bool {
        self.is_left_uniform() && self.is_right_uniform()
    }

    /// Compare two vectors with QPA's tip-monomial/descending-position order.
    ///
    /// # Errors
    ///
    /// Rejects vectors from distinct module families and propagates path
    /// ordering failures.
    pub fn compare(&self, right: &Self) -> Result<Ordering, ProjectivePathModuleError> {
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

    /// Translate QPA's `IsLeftDivisible(x, y)` predicate.
    ///
    /// The source additionally requires the dividend tip walk to be nonempty,
    /// so a vertex tip is never left-divisible here.
    ///
    /// # Errors
    ///
    /// Rejects distinct vector families and explicit zero inputs.
    pub fn is_left_divisible(&self, divisor: &Self) -> Result<bool, ProjectivePathModuleError> {
        self.ensure_same_family(divisor)?;
        if self.is_zero() || divisor.is_zero() {
            return Err(ProjectivePathModuleError::ZeroVectorDivision);
        }
        let dividend_path = self.leading_path()?;
        let divisor_path = divisor.leading_path()?;
        let dividend_walk = path_walk(&dividend_path)?;
        let divisor_walk = path_walk(&divisor_path)?;
        Ok(self.leading_position == divisor.leading_position
            && !dividend_walk.is_empty()
            && dividend_walk.starts_with(divisor_walk))
    }

    fn ensure_same_family(&self, other: &Self) -> Result<(), ProjectivePathModuleError> {
        if self.module.is_identical(&other.module) {
            Ok(())
        } else {
            Err(ProjectivePathModuleError::DifferentVectorFamilies)
        }
    }

    fn leading_path(&self) -> Result<Path, ProjectivePathModuleError> {
        if self.is_zero() {
            return Err(ProjectivePathModuleError::ZeroVectorDivision);
        }
        self.components
            .get(self.leading_position - 1)
            .map(PathAlgebraElement::leading_monomial)
            .ok_or(ProjectivePathModuleError::InvalidLeadingPosition {
                position: self.leading_position,
                length: self.components.len(),
            })
    }
}

impl<F: Field> PartialEq for PathAlgebraVector<F> {
    fn eq(&self, other: &Self) -> bool {
        self.module.is_identical(&other.module) && self.components == other.components
    }
}

impl<F: Field> Eq for PathAlgebraVector<F> {}

/// Translate `moduleprojres.gi`'s `LeftDivision(x, y)`.
///
/// `None` is QPA's `false`; success returns the coefficient-bearing right path
/// factor whose multiplication makes the two vector tips equal. Unlike
/// [`PathAlgebraVector::is_left_divisible`], equal vertex tips are accepted.
///
/// # Errors
///
/// QPA explicitly rejects zero inputs. Rust also rejects different vector
/// families and propagates exact path/field arithmetic.
pub fn left_division<F: Field>(
    dividend: &PathAlgebraVector<F>,
    divisor: &PathAlgebraVector<F>,
) -> Result<Option<PathAlgebraElement<F>>, ProjectivePathModuleError> {
    dividend.ensure_same_family(divisor)?;
    if dividend.is_zero() || divisor.is_zero() {
        return Err(ProjectivePathModuleError::ZeroVectorDivision);
    }
    let dividend_path = dividend.leading_path()?;
    let divisor_path = divisor.leading_path()?;
    let dividend_walk = path_walk(&dividend_path)?;
    let divisor_walk = path_walk(&divisor_path)?;
    if dividend.leading_position != divisor.leading_position
        || !dividend_walk.starts_with(divisor_walk)
    {
        return Ok(None);
    }
    let field = dividend.module.algebra().field();
    let coefficient = field.div(
        &dividend.leading_coefficient()?,
        &divisor.leading_coefficient()?,
    );
    let factor = right_path_factor(
        dividend.module.algebra(),
        &dividend_path,
        divisor_walk.len(),
    )?
    .scale_left(&coefficient)?;
    Ok(Some(factor))
}

/// Translate QPA's one-list `TipReduce(M)` overload.
///
/// The ordered restart loop removes zeros, scans each vector against every
/// other current vector, performs only the first available reduction, and
/// restarts from the beginning.
///
/// # Errors
///
/// Propagates exact left-divisibility and reduction failures.
pub fn tip_reduce_path_algebra_vectors<F: Field>(
    vectors: &[PathAlgebraVector<F>],
) -> Result<Vec<PathAlgebraVector<F>>, ProjectivePathModuleError> {
    let mut reduced = vectors.to_vec();
    loop {
        reduced.retain(|vector| !vector.is_zero());
        let mut changed = false;
        'scan: for left in 0..reduced.len() {
            for right in 0..reduced.len() {
                if left == right {
                    continue;
                }
                if reduced[left].is_left_divisible(&reduced[right])? {
                    reduced[left] = reduce_right_module_element(&reduced[left], &reduced[right])?;
                    changed = true;
                    break 'scan;
                }
            }
        }
        if !changed {
            return Ok(reduced);
        }
    }
}

/// Translate QPA's two-argument `TipReduce(H, el)` overload.
///
/// Reducers are tried in entered order and the scan restarts after every
/// successful leading-term cancellation.
///
/// # Errors
///
/// Propagates exact left-divisibility and reduction failures.
pub fn tip_reduce_path_algebra_vector<F: Field>(
    reducers: &[PathAlgebraVector<F>],
    element: &PathAlgebraVector<F>,
) -> Result<PathAlgebraVector<F>, ProjectivePathModuleError> {
    let mut reduced = element.clone();
    while !reduced.is_zero() {
        let mut changed = false;
        for reducer in reducers {
            if reduced.is_left_divisible(reducer)? {
                reduced = reduce_right_module_element(&reduced, reducer)?;
                changed = true;
                break;
            }
        }
        if !changed {
            break;
        }
    }
    Ok(reduced)
}

fn reduce_right_module_element<F: Field>(
    dividend: &PathAlgebraVector<F>,
    divisor: &PathAlgebraVector<F>,
) -> Result<PathAlgebraVector<F>, ProjectivePathModuleError> {
    dividend.ensure_same_family(divisor)?;
    let dividend_path = dividend.leading_path()?;
    let divisor_path = divisor.leading_path()?;
    let dividend_walk = path_walk(&dividend_path)?;
    let divisor_walk = path_walk(&divisor_path)?;
    if dividend.leading_position != divisor.leading_position
        || !dividend_walk.starts_with(divisor_walk)
    {
        return Err(ProjectivePathModuleError::ReductionByNondivisor);
    }
    let field = dividend.module.algebra().field();
    let coefficient = field.div(
        &dividend.leading_coefficient()?,
        &divisor.leading_coefficient()?,
    );
    let factor = if divisor_walk.len() == dividend_walk.len() {
        dividend.module.algebra().one()
    } else {
        right_path_factor(
            dividend.module.algebra(),
            &dividend_path,
            divisor_walk.len(),
        )?
    }
    .scale_left(&coefficient)?;
    dividend.sub(&divisor.act(&factor)?)
}

fn leading_position<F: Field>(
    algebra: &PathAlgebra<F>,
    components: &[PathAlgebraElement<F>],
) -> Result<usize, ProjectivePathModuleError> {
    let Some(first) = components.first() else {
        return Ok(1);
    };
    let mut position = 1;
    let mut greatest = first.leading_monomial();
    for (index, component) in components.iter().enumerate().skip(1) {
        let candidate = component.leading_monomial();
        if algebra.ordering().compare(&candidate, &greatest)?.is_gt() {
            greatest = candidate;
            position = index + 1;
        }
    }
    Ok(position)
}

fn common_nonzero_endpoint<F: Field>(
    components: &[PathAlgebraElement<F>],
    endpoint: impl Fn(&Path) -> PathEndpoint,
) -> bool {
    let mut endpoints = components
        .iter()
        .filter(|component| !component.is_zero())
        .filter_map(|component| match endpoint(component.terms()[0].path()) {
            PathEndpoint::Vertex(vertex) => Some(vertex),
            PathEndpoint::Zero(_) => None,
        });
    let Some(first) = endpoints.next() else {
        return true;
    };
    endpoints.all(|candidate| candidate == first)
}

fn path_walk(path: &Path) -> Result<&[crate::ArrowIndex], ProjectivePathModuleError> {
    path.walk()
        .ok_or(ProjectivePathModuleError::PathHasNoArrowWalk)
}

fn right_path_factor<F: Field>(
    algebra: &PathAlgebra<F>,
    dividend_path: &Path,
    divisor_length: usize,
) -> Result<PathAlgebraElement<F>, ProjectivePathModuleError> {
    let walk = path_walk(dividend_path)?;
    if divisor_length == walk.len() {
        let PathEndpoint::Vertex(target) = dividend_path.target() else {
            return Err(ProjectivePathModuleError::PathHasNoVertexTarget);
        };
        return Ok(algebra.element(target.path())?);
    }
    let arrows = walk[divisor_length..]
        .iter()
        .map(|index| {
            algebra
                .quiver()
                .arrow(*index)
                .ok_or(ProjectivePathModuleError::MissingPathArrow)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let path = Path::from_arrow_product(algebra.quiver(), &arrows)?
        .ok_or(ProjectivePathModuleError::EmptyRightFactor)?;
    Ok(algebra.element(path)?)
}

/// Failure in QPA's projective path-algebra vector machinery.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ProjectivePathModuleError {
    /// The vertex element belongs to another path algebra.
    #[error("projective vertex {position} belongs to another path algebra")]
    ForeignVertexAlgebra {
        /// One-based entered position.
        position: usize,
    },
    /// QPA's vertex overload accepts only coefficient-one trivial paths.
    #[error("projective entry {position} is not an embedded vertex")]
    NotEmbeddedVertex {
        /// One-based entered position.
        position: usize,
    },
    /// An internal vertex index is outside the parent quiver.
    #[error("projective vertex {vertex:?} is outside the parent quiver")]
    VertexOutOfRange {
        /// Rejected vertex index.
        vertex: VertexIndex,
    },
    /// Checked vectors must have the family's fixed number of coordinates.
    #[error("path-algebra vector has {actual} components, expected {expected}")]
    ComponentCount {
        /// Entered count.
        actual: usize,
        /// Family count.
        expected: usize,
    },
    /// Every component must use the family's path algebra.
    #[error("a vector component belongs to another path algebra")]
    ForeignComponentAlgebra,
    /// QPA checks left uniformity from the final coordinate backwards.
    #[error("path-algebra vector component {position} is not left uniform")]
    ComponentNotLeftUniform {
        /// One-based coordinate.
        position: usize,
    },
    /// A nonzero coordinate must start at its family vertex.
    #[error("path-algebra vector component {position} starts at the wrong vertex")]
    ComponentSourceMismatch {
        /// One-based coordinate.
        position: usize,
    },
    /// QPA installs vector arithmetic only for identical families.
    #[error("the path-algebra vectors belong to different module families")]
    DifferentVectorFamilies,
    /// A right action element must belong to the module's parent algebra.
    #[error("the right action element belongs to another path algebra")]
    ForeignActionAlgebra,
    /// QPA indexes a nonexistent coordinate for this empty projective.
    #[error("an empty path-algebra vector has no leading component")]
    EmptyVectorLeadingComponent,
    /// The stored positional representation violated its family length.
    #[error("leading position {position} is outside a vector of length {length}")]
    InvalidLeadingPosition {
        /// Stored one-based position.
        position: usize,
        /// Component count.
        length: usize,
    },
    /// QPA's `LeftDivision` explicitly rejects zeroes.
    #[error("path-algebra vector division does not accept zero inputs")]
    ZeroVectorDivision,
    /// `WalkOfPath` is unavailable for this tip path.
    #[error("a vector tip has no arrow walk")]
    PathHasNoArrowWalk,
    /// A nonzero path factor must have a vertex target.
    #[error("a vector tip has no vertex target")]
    PathHasNoVertexTarget,
    /// A stored path walk referenced no quiver arrow.
    #[error("a vector tip references a missing quiver arrow")]
    MissingPathArrow,
    /// The successful division branch cannot have an empty generic path product.
    #[error("a nonempty vector right factor produced an empty path product")]
    EmptyRightFactor,
    /// Internal tip reduction called its reducer on a nondivisible vector.
    #[error("tip reduction selected a vector that does not divide the current tip")]
    ReductionByNondivisor,
    /// Parent path-algebra arithmetic failed.
    #[error(transparent)]
    Algebra(#[from] PathAlgebraError),
    /// Path construction failed.
    #[error(transparent)]
    Path(#[from] PathError),
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
        PathAlgebraVector, ProjectivePathModuleError, left_division,
        right_projective_path_algebra_module, tip_reduce_path_algebra_vector,
        tip_reduce_path_algebra_vectors,
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

    fn vector_json<F: Field>(
        vector: &PathAlgebraVector<F>,
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

    #[allow(clippy::too_many_lines)]
    fn projective_path_vector_case<F: Field>(
        field: &F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let algebra = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(
                3,
                vec![
                    ArrowSpec::named(vertex(1), vertex(2), "a"),
                    ArrowSpec::named(vertex(1), vertex(2), "b"),
                    ArrowSpec::named(vertex(2), vertex(3), "c"),
                ],
            )
            .unwrap(),
        )
        .unwrap();
        let generators = algebra.generators();
        let (vertex_v, vertex_w, arrow_a, arrow_b, arrow_c) = (
            &generators[0],
            &generators[1],
            &generators[3],
            &generators[4],
            &generators[5],
        );
        let zero = algebra.zero();
        let projective = right_projective_path_algebra_module(
            &algebra,
            &[vertex_v.clone(), vertex_v.clone(), vertex_w.clone()],
        )
        .unwrap();
        let bc = arrow_b.multiply(arrow_c).unwrap();
        let ac = arrow_a.multiply(arrow_c).unwrap();
        let p1 = projective
            .vector(vec![bc.clone(), ac, arrow_c.clone()])
            .unwrap();
        let p2 = projective
            .vector(vec![arrow_a.clone(), arrow_b.clone(), vertex_w.clone()])
            .unwrap();
        let reducible = projective
            .vector(vec![bc.add(arrow_a).unwrap(), zero.clone(), zero.clone()])
            .unwrap();
        let divisor = projective
            .vector(vec![arrow_b.clone(), zero.clone(), zero.clone()])
            .unwrap();
        let exact = projective
            .vector(vec![bc, zero.clone(), zero.clone()])
            .unwrap();
        let nonuniform = projective
            .vector(vec![arrow_b.clone(), vertex_v.clone(), zero.clone()])
            .unwrap();
        let vertex_vector = projective
            .vector(vec![vertex_v.clone(), zero.clone(), zero.clone()])
            .unwrap();
        let checked_zero = projective
            .vector(vec![zero.clone(), zero.clone(), zero])
            .unwrap();
        json!({
            "dimension": projective.dimension().unwrap(),
            "generators": vector_list_json(&projective.generators().unwrap(), coefficient),
            "p1": vector_json(&p1, coefficient),
            "p2": vector_json(&p2, coefficient),
            "sum": vector_json(&p1.add(&p2).unwrap(), coefficient),
            "difference": vector_json(&p1.sub(&p2).unwrap(), coefficient),
            "negated": vector_json(&p2.neg().unwrap(), coefficient),
            "action_c": vector_json(&p2.act(arrow_c).unwrap(), coefficient),
            "leading": {
                "position": p1.leading_position(),
                "coefficient": coefficient(&p1.leading_coefficient().unwrap()),
                "component": element_json(p1.leading_component().unwrap(), coefficient),
                "term": vector_json(&p1.leading_term(), coefficient),
                "target": p1.target_vertex().unwrap().external_representation(),
                "source": p1.source_vertex().unwrap().external_representation(),
            },
            "uniformity": {
                "p1": [p1.is_left_uniform(), p1.is_right_uniform(), p1.is_uniform()],
                "divisor": [
                    divisor.is_left_uniform(),
                    divisor.is_right_uniform(),
                    divisor.is_uniform(),
                ],
            },
            "ordering": {
                "p1_lt_p2": p1.compare(&p2).unwrap().is_lt(),
                "p2_lt_p1": p2.compare(&p1).unwrap().is_lt(),
            },
            "zero_positions": {
                "module_zero": projective.zero().leading_position(),
                "checked_zero": checked_zero.leading_position(),
            },
            "division": {
                "divisible": exact.is_left_divisible(&divisor).unwrap(),
                "vertex_divisible": vertex_vector.is_left_divisible(&vertex_vector).unwrap(),
                "factor": element_json(
                    &left_division(&exact, &divisor).unwrap().unwrap(),
                    coefficient,
                ),
                "vertex_factor": element_json(
                    &left_division(&vertex_vector, &vertex_vector).unwrap().unwrap(),
                    coefficient,
                ),
                "wrong_position": p1.is_left_divisible(&p2).unwrap(),
            },
            "tip_reduce_list": vector_list_json(
                &tip_reduce_path_algebra_vectors(&[reducible.clone(), divisor.clone()]).unwrap(),
                coefficient,
            ),
            "tip_reduce_exact_list": vector_list_json(
                &tip_reduce_path_algebra_vectors(&[exact, divisor.clone()]).unwrap(),
                coefficient,
            ),
            "tip_reduce_element": vector_json(
                &tip_reduce_path_algebra_vector(&[divisor], &reducible).unwrap(),
                coefficient,
            ),
            "tip_reduce_nonuniform_exact": vector_json(
                &tip_reduce_path_algebra_vector(
                    std::slice::from_ref(&nonuniform),
                    &nonuniform,
                )
                .unwrap(),
                coefficient,
            ),
        })
    }

    #[test]
    fn pinned_qpa_projective_path_vectors_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/projpathalgmodule_vectors.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "projpathalgmodule_vectors",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": projective_path_vector_case(&Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": projective_path_vector_case(&field5, |entry| json!(
                field5_coefficients
                    .to_integer(entry)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
            )),
            "rational": projective_path_vector_case(&Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn projective_path_vectors_preserve_family_and_validation_boundaries() {
        let algebra = PathAlgebra::new(
            Z2,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap();
        let generators = algebra.generators();
        let first =
            right_projective_path_algebra_module(&algebra, &[generators[0].clone()]).unwrap();
        let second =
            right_projective_path_algebra_module(&algebra, &[generators[0].clone()]).unwrap();
        let first_vector = first.vector(vec![generators[0].clone()]).unwrap();
        let second_vector = second.vector(vec![generators[0].clone()]).unwrap();
        assert!(matches!(
            first_vector.add(&second_vector),
            Err(ProjectivePathModuleError::DifferentVectorFamilies)
        ));
        assert!(matches!(
            first.vector(vec![]),
            Err(ProjectivePathModuleError::ComponentCount { .. })
        ));
        assert!(matches!(
            first.vector(vec![generators[1].clone()]),
            Err(ProjectivePathModuleError::ComponentSourceMismatch { position: 1 })
        ));
        assert!(matches!(
            left_division(&first.zero(), &first_vector),
            Err(ProjectivePathModuleError::ZeroVectorDivision)
        ));
    }
}
