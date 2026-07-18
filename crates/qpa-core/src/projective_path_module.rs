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
    Path, PathAlgebra, PathAlgebraElement, PathAlgebraError, PathAlgebraIdeal, PathEndpoint,
    PathError, QuiverStaticDictionary, StaticDictionaryError, VertexIndex,
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

    /// Translate QPA's canonical basis of a vertex-projective path module.
    ///
    /// QPA first iterates summands and the canonical path-algebra basis, then
    /// sorts the resulting vectors by path tip and descending component
    /// position. `None` is QPA's unavailable basis for a cyclic quiver.
    ///
    /// # Errors
    ///
    /// Propagates vector construction and path-ordering failures.
    pub fn canonical_basis(
        &self,
    ) -> Result<Option<Vec<PathAlgebraVector<F>>>, ProjectivePathModuleError> {
        let Some(algebra_basis) = self.algebra().canonical_basis() else {
            return Ok(None);
        };
        let mut basis = Vec::new();
        for (position, vertex) in self.vertex_indices().iter().copied().enumerate() {
            for element in &algebra_basis {
                if !matches!(
                    element.leading_monomial().source(),
                    PathEndpoint::Vertex(source) if source.index() == vertex
                ) {
                    continue;
                }
                let mut components = vec![self.algebra().zero(); self.vertex_indices().len()];
                components[position] = element.clone();
                basis.push(self.vector(components)?);
            }
        }
        sort_path_algebra_vectors(&mut basis)?;
        Ok(Some(basis))
    }

    /// Translate QPA's structural `NewBasis(self, vectors)` overload.
    ///
    /// # Errors
    ///
    /// Rust makes QPA's homogeneous-family dispatch boundary explicit.
    pub fn new_basis(
        &self,
        vectors: Vec<PathAlgebraVector<F>>,
    ) -> Result<PathAlgebraVectorBasis<F>, ProjectivePathModuleError> {
        if vectors
            .iter()
            .any(|vector| !self.is_identical(vector.module()))
        {
            return Err(ProjectivePathModuleError::ForeignVectorBasisVector);
        }
        Ok(PathAlgebraVectorBasis::structural_for_projective(
            self.clone(),
            vectors,
        ))
    }

    /// Construct QPA's `SubAlgebraModule` with entered generators retained.
    ///
    /// # Errors
    ///
    /// QPA's installed method does not dispatch on an empty GAP list. Every
    /// entered generator must belong to this identical projective-vector family.
    pub fn submodule(
        &self,
        generators: Vec<PathAlgebraVector<F>>,
    ) -> Result<RightProjectivePathAlgebraSubmodule<F>, ProjectivePathModuleError> {
        RightProjectivePathAlgebraSubmodule::new(self, generators, false)
    }

    /// Translate QPA's `SubAlgebraModule(self, generators, "basis")` overload.
    ///
    /// QPA trusts the entered nonempty list as a field basis, stores it as the
    /// left-module generators, and sets the dimension to its length without an
    /// independence test. Rust exposes that literal branch as a named method
    /// instead of accepting an arbitrary sentinel string.
    ///
    /// # Errors
    ///
    /// Uses the same family and nonempty-list checks as [`Self::submodule`].
    pub fn submodule_with_basis(
        &self,
        generators: Vec<PathAlgebraVector<F>>,
    ) -> Result<RightProjectivePathAlgebraSubmodule<F>, ProjectivePathModuleError> {
        RightProjectivePathAlgebraSubmodule::new(self, generators, true)
    }

    /// Translate `UniformGeneratorsOfModule` for the whole projective module.
    ///
    /// # Errors
    ///
    /// Propagates generator construction, right action, and ordering failures.
    pub fn uniform_generators(
        &self,
    ) -> Result<Vec<PathAlgebraVector<F>>, ProjectivePathModuleError> {
        uniform_generators_of_module(self, &self.generators()?)
    }

    /// Translate `RightGroebnerBasisOfModule` for the whole projective module.
    ///
    /// # Errors
    ///
    /// Propagates uniformization, tip reduction, and dictionary failures.
    pub fn right_groebner_basis(
        &self,
    ) -> Result<RightPathAlgebraModuleGroebnerBasis<F>, ProjectivePathModuleError> {
        let uniform_generators = self.uniform_generators()?;
        RightPathAlgebraModuleGroebnerBasis::new(
            RightPathAlgebraModuleSource::Whole(self.clone()),
            &uniform_generators,
        )
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

    /// Translate QPA's pure-path-algebra `PathAlgebraVectorNC` constructor.
    ///
    /// The source performs no component, family, length, or position checks,
    /// and ignores its `normal` Boolean when no quotient normal-form function
    /// is installed. This deliberately retains that unsafe semantic boundary.
    #[must_use]
    pub fn vector_nc(
        &self,
        components: Vec<PathAlgebraElement<F>>,
        leading_position: usize,
        _normal: bool,
    ) -> PathAlgebraVector<F> {
        PathAlgebraVector::new_unchecked(self.clone(), components, leading_position)
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

struct RightProjectivePathAlgebraSubmoduleData<F: Field> {
    parent: RightProjectivePathAlgebraModule<F>,
    generators: Box<[PathAlgebraVector<F>]>,
    entered_generators_are_basis: bool,
}

/// QPA's `SubAlgebraModule` inside a vertex-projective path module.
///
/// The entered generating list is retained verbatim, including zero and
/// redundant generators. Clones preserve QPA object identity.
#[derive(Clone)]
pub struct RightProjectivePathAlgebraSubmodule<F: Field>(
    Arc<RightProjectivePathAlgebraSubmoduleData<F>>,
);

impl<F: Field> RightProjectivePathAlgebraSubmodule<F> {
    fn new(
        parent: &RightProjectivePathAlgebraModule<F>,
        generators: Vec<PathAlgebraVector<F>>,
        entered_generators_are_basis: bool,
    ) -> Result<Self, ProjectivePathModuleError> {
        if generators.is_empty() {
            return Err(ProjectivePathModuleError::EmptySubmoduleGenerators);
        }
        if generators
            .iter()
            .any(|generator| !generator.module().is_identical(parent))
        {
            return Err(ProjectivePathModuleError::ForeignSubmoduleGenerator);
        }
        Ok(Self(Arc::new(RightProjectivePathAlgebraSubmoduleData {
            parent: parent.clone(),
            generators: generators.into_boxed_slice(),
            entered_generators_are_basis,
        })))
    }

    /// Return the ambient vertex-projective module.
    #[must_use]
    pub fn parent(&self) -> &RightProjectivePathAlgebraModule<F> {
        &self.0.parent
    }

    /// Return QPA's entered `GeneratorsOfAlgebraModule` list.
    #[must_use]
    pub fn generators(&self) -> &[PathAlgebraVector<F>] {
        &self.0.generators
    }

    /// Return whether two handles denote the identical QPA submodule object.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Translate QPA's structural `NewBasis(self, vectors)` overload.
    ///
    /// This does not test that the vectors span or belong to this submodule;
    /// QPA's method only relies on their common ambient element family.
    ///
    /// # Errors
    ///
    /// Rejects vectors from another projective-vector family.
    pub fn new_basis(
        &self,
        vectors: Vec<PathAlgebraVector<F>>,
    ) -> Result<PathAlgebraVectorBasis<F>, ProjectivePathModuleError> {
        if vectors
            .iter()
            .any(|vector| !self.parent().is_identical(vector.module()))
        {
            return Err(ProjectivePathModuleError::ForeignVectorBasisVector);
        }
        Ok(PathAlgebraVectorBasis::structural_for_submodule(
            self.clone(),
            vectors,
        ))
    }

    /// Translate QPA's vertex-idempotent splitting into uniform generators.
    ///
    /// # Errors
    ///
    /// Propagates exact right action and vector-ordering failures.
    pub fn uniform_generators(
        &self,
    ) -> Result<Vec<PathAlgebraVector<F>>, ProjectivePathModuleError> {
        uniform_generators_of_module(self.parent(), self.generators())
    }

    /// Translate QPA's right module Gröbner-basis construction.
    ///
    /// # Errors
    ///
    /// Propagates uniformization, tip reduction, and dictionary failures.
    pub fn right_groebner_basis(
        &self,
    ) -> Result<RightPathAlgebraModuleGroebnerBasis<F>, ProjectivePathModuleError> {
        let uniform_generators = self.uniform_generators()?;
        RightPathAlgebraModuleGroebnerBasis::new(
            RightPathAlgebraModuleSource::Submodule(self.clone()),
            &uniform_generators,
        )
    }

    /// Translate QPA's `Basis` method for a projective path submodule.
    ///
    /// For the three-argument `"basis"` constructor this is the exact entered
    /// list. Otherwise QPA subtracts each ambient canonical vector's complete
    /// normal form modulo a right Gröbner basis, removes zeroes and duplicates
    /// with `AsSet`, and stores the resulting structural basis.
    /// `None` is QPA's unavailable canonical basis for a cyclic path algebra.
    ///
    /// # Errors
    ///
    /// Propagates canonical-basis, Gröbner, reduction, and ordering failures.
    pub fn basis(&self) -> Result<Option<PathAlgebraVectorBasis<F>>, ProjectivePathModuleError> {
        if self.0.entered_generators_are_basis {
            return Ok(Some(self.new_basis(self.generators().to_vec())?));
        }
        let Some(ambient_basis) = self.parent().canonical_basis()? else {
            return Ok(None);
        };
        let groebner_basis = self.right_groebner_basis()?;
        let mut vectors = Vec::new();
        for ambient_vector in ambient_basis {
            let normal_form = groebner_basis.completely_reduce(&ambient_vector)?;
            let vector = ambient_vector.sub(&normal_form)?;
            if !vector.is_zero() {
                add_set_path_algebra_vector(&mut vectors, vector)?;
            }
        }
        Ok(Some(self.new_basis(vectors)?))
    }

    /// Return the dimension installed by QPA's projective-submodule basis.
    ///
    /// `None` means the ambient canonical basis is unavailable.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::basis`] failures.
    pub fn dimension(&self) -> Result<Option<usize>, ProjectivePathModuleError> {
        Ok(self.basis()?.map(|basis| basis.basis_vectors().len()))
    }
}

impl<F: Field> fmt::Debug for RightProjectivePathAlgebraSubmodule<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("RightProjectivePathAlgebraSubmodule")
            .field("parent", self.parent())
            .field("generators", &self.generators())
            .finish()
    }
}

impl<F: Field> PartialEq for RightProjectivePathAlgebraSubmodule<F> {
    fn eq(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}

impl<F: Field> Eq for RightProjectivePathAlgebraSubmodule<F> {}

#[derive(Clone, Debug)]
enum RightPathAlgebraModuleSource<F: Field> {
    Whole(RightProjectivePathAlgebraModule<F>),
    Submodule(RightProjectivePathAlgebraSubmodule<F>),
}

impl<F: Field> RightPathAlgebraModuleSource<F> {
    fn parent(&self) -> &RightProjectivePathAlgebraModule<F> {
        match self {
            Self::Whole(module) => module,
            Self::Submodule(module) => module.parent(),
        }
    }

    fn is_whole(&self) -> bool {
        matches!(self, Self::Whole(_))
    }
}

/// QPA's right Gröbner-basis object for a projective path module or submodule.
#[derive(Clone, Debug)]
pub struct RightPathAlgebraModuleGroebnerBasis<F: Field> {
    underlying: RightPathAlgebraModuleSource<F>,
    grouped_basis: Vec<Vec<PathAlgebraVector<F>>>,
    static_dictionaries: Vec<Option<QuiverStaticDictionary>>,
    basis_vectors: Vec<PathAlgebraVector<F>>,
}

impl<F: Field> RightPathAlgebraModuleGroebnerBasis<F> {
    fn new(
        underlying: RightPathAlgebraModuleSource<F>,
        uniform_generators: &[PathAlgebraVector<F>],
    ) -> Result<Self, ProjectivePathModuleError> {
        let parent = underlying.parent();
        let reduced = tip_reduce_path_algebra_vectors(uniform_generators)?;
        let mut grouped_basis = vec![Vec::new(); parent.vertex_indices().len()];
        for vector in reduced {
            let position = vector.leading_position().checked_sub(1).ok_or(
                ProjectivePathModuleError::InvalidLeadingPosition {
                    position: vector.leading_position(),
                    length: parent.vertex_indices().len(),
                },
            )?;
            let Some(group) = grouped_basis.get_mut(position) else {
                return Err(ProjectivePathModuleError::InvalidLeadingPosition {
                    position: vector.leading_position(),
                    length: parent.vertex_indices().len(),
                });
            };
            group.push(vector);
        }

        let static_dictionaries = grouped_basis
            .iter()
            .map(|group| {
                if group.is_empty() {
                    return Ok(None);
                }
                let patterns = group
                    .iter()
                    .map(PathAlgebraVector::leading_path)
                    .collect::<Result<Vec<_>, _>>()?;
                Ok(Some(QuiverStaticDictionary::new(
                    parent.algebra().quiver(),
                    &patterns,
                )?))
            })
            .collect::<Result<Vec<_>, ProjectivePathModuleError>>()?;
        let basis_vectors = grouped_basis.iter().flatten().cloned().collect();
        Ok(Self {
            underlying,
            grouped_basis,
            static_dictionaries,
            basis_vectors,
        })
    }

    /// Return the ambient vertex-projective module/vector family.
    #[must_use]
    pub fn parent_module(&self) -> &RightProjectivePathAlgebraModule<F> {
        self.underlying.parent()
    }

    /// Return whether QPA constructed this basis for the whole projective.
    #[must_use]
    pub fn is_for_whole_module(&self) -> bool {
        self.underlying.is_whole()
    }

    /// Return QPA's exact underlying whole module, when applicable.
    #[must_use]
    pub fn underlying_projective_module(&self) -> Option<&RightProjectivePathAlgebraModule<F>> {
        match &self.underlying {
            RightPathAlgebraModuleSource::Whole(module) => Some(module),
            RightPathAlgebraModuleSource::Submodule(_) => None,
        }
    }

    /// Return QPA's exact underlying submodule, when applicable.
    #[must_use]
    pub fn underlying_submodule(&self) -> Option<&RightProjectivePathAlgebraSubmodule<F>> {
        match &self.underlying {
            RightPathAlgebraModuleSource::Whole(_) => None,
            RightPathAlgebraModuleSource::Submodule(module) => Some(module),
        }
    }

    /// Return QPA's `BasisVectors`, flattened by leading coordinate.
    #[must_use]
    pub fn basis_vectors(&self) -> &[PathAlgebraVector<F>] {
        &self.basis_vectors
    }

    /// Translate `CompletelyReduce(gb, vector)`.
    ///
    /// QPA selects the first shortest prefix match in the static dictionary;
    /// unreducible leading terms are moved to the result one at a time.
    ///
    /// # Errors
    ///
    /// The vector must belong to the identical ambient family. Exact path,
    /// field, and static-dictionary failures are propagated.
    pub fn completely_reduce(
        &self,
        vector: &PathAlgebraVector<F>,
    ) -> Result<PathAlgebraVector<F>, ProjectivePathModuleError> {
        if !vector.module().is_identical(self.parent_module()) {
            return Err(ProjectivePathModuleError::DifferentVectorFamilies);
        }
        if vector.is_zero() {
            return Ok(vector.clone());
        }

        let mut current = vector.clone();
        let mut normal_form = vector.module().zero();
        loop {
            let position = current.leading_position() - 1;
            let path = current.leading_path()?;
            let matches = match &self.static_dictionaries[position] {
                Some(dictionary) => dictionary.prefix_search(&path)?,
                None => Vec::new(),
            };
            if let Some(pattern) = matches
                .first()
                .and_then(|matched| matched.pattern_indices.first())
            {
                let reducer = &self.grouped_basis[position][pattern - 1];
                current = reduce_right_module_element(&current, reducer)?;
            } else {
                let leading_term = current.leading_term();
                normal_form = normal_form.add(&leading_term)?;
                current = current.sub(&leading_term)?;
            }
            if current.is_zero() {
                return Ok(normal_form);
            }
        }
    }
}

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

    /// Return QPA's stored `IsNormalForm` filter.
    ///
    /// Pure path-projective families create their default vector type with
    /// this filter set, so the generic `ReturnFalse` method is not consulted.
    #[must_use]
    pub const fn is_normal_form(&self) -> bool {
        true
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

/// A free left field-module whose elements are QPA path-algebra vectors.
///
/// This is the typed counterpart of the `LeftModuleByGenerators` objects on
/// which `projpathalgmodule.gi` installs `NewBasis`, `BasisNC`, `Basis`,
/// `BasisOfDomain`, `Coefficients`, and membership.
#[derive(Clone, Debug)]
pub struct PathAlgebraVectorSpace<F: Field> {
    module: RightProjectivePathAlgebraModule<F>,
    generators: Vec<PathAlgebraVector<F>>,
}

impl<F: Field> PathAlgebraVectorSpace<F> {
    /// Construct a field-space from entered path-vector generators.
    ///
    /// # Errors
    ///
    /// Every generator must belong to the identical projective-vector family.
    pub fn by_generators(
        module: RightProjectivePathAlgebraModule<F>,
        generators: Vec<PathAlgebraVector<F>>,
    ) -> Result<Self, ProjectivePathModuleError> {
        if generators
            .iter()
            .any(|generator| !module.is_identical(generator.module()))
        {
            return Err(ProjectivePathModuleError::ForeignVectorSpaceGenerator);
        }
        Ok(Self { module, generators })
    }

    /// Return the identical ambient projective-vector family.
    #[must_use]
    pub const fn module(&self) -> &RightProjectivePathAlgebraModule<F> {
        &self.module
    }

    /// Return `GeneratorsOfLeftModule` in entered order.
    #[must_use]
    pub fn generators(&self) -> &[PathAlgebraVector<F>] {
        &self.generators
    }

    /// Translate QPA's structural `NewBasis(V, vectors)` method.
    ///
    /// No echelon data is attached, so `Coefficients` on the result is not
    /// available—exactly as in QPA.
    ///
    /// # Errors
    ///
    /// Rust makes QPA's implicit homogeneous-family boundary explicit.
    pub fn new_basis(
        &self,
        vectors: Vec<PathAlgebraVector<F>>,
    ) -> Result<PathAlgebraVectorBasis<F>, ProjectivePathModuleError> {
        self.ensure_ambient_family(&vectors)?;
        Ok(PathAlgebraVectorBasis::structural_for_space(
            self.clone(),
            vectors,
        ))
    }

    /// Translate QPA's unchecked `BasisNC(V, vectors)` method.
    ///
    /// `Ok(None)` is GAP's `fail` if QPA's own triangularization removes a
    /// dependent nonzero vector. The ordered basis vectors remain the entered
    /// nonzero vectors; echelon and sparse base-change data are stored apart.
    ///
    /// # Errors
    ///
    /// Rejects a foreign vector family and propagates exact arithmetic and
    /// path-ordering failures.
    pub fn basis_nc(
        &self,
        vectors: Vec<PathAlgebraVector<F>>,
    ) -> Result<Option<PathAlgebraVectorBasis<F>>, ProjectivePathModuleError> {
        self.ensure_ambient_family(&vectors)?;
        let vectors = vectors
            .into_iter()
            .filter(|vector| !vector.is_zero())
            .collect::<Vec<_>>();
        let triangular = triangularize_path_algebra_vector_list(&vectors)?;
        if triangular.echelon_basis.len() != vectors.len() {
            return Ok(None);
        }
        Ok(Some(PathAlgebraVectorBasis::echelon_for_space(
            self.clone(),
            vectors,
            triangular,
        )))
    }

    /// Translate QPA's checked `Basis(V, vectors)` overload.
    ///
    /// `Ok(None)` is GAP's `fail` when membership or independence fails.
    /// Membership intentionally uses QPA's tip-scanning coefficient method,
    /// including its observable early-stop behavior.
    ///
    /// # Errors
    ///
    /// Rejects a foreign vector family and propagates exact basis arithmetic.
    pub fn basis_from_vectors(
        &self,
        vectors: Vec<PathAlgebraVector<F>>,
    ) -> Result<Option<PathAlgebraVectorBasis<F>>, ProjectivePathModuleError> {
        self.ensure_ambient_family(&vectors)?;
        for vector in &vectors {
            if !self.contains(vector)? {
                return Ok(None);
            }
        }
        self.basis_nc(vectors)
    }

    /// Translate QPA's `BasisOfDomain(V)` method.
    ///
    /// QPA triangularizes the left-module generators, uses those echelon
    /// vectors themselves as `BasisVectors`, and resets the sparse base-change
    /// matrix to the identity.
    ///
    /// # Errors
    ///
    /// Propagates exact triangularization failures.
    pub fn basis_of_domain(&self) -> Result<PathAlgebraVectorBasis<F>, ProjectivePathModuleError> {
        let triangular = triangularize_path_algebra_vector_list(&self.generators)?;
        let vectors = triangular.echelon_basis.clone();
        let field = self.module.algebra().field();
        let base_change = (0..vectors.len())
            .map(|position| vec![(position, field.one())])
            .collect();
        Ok(PathAlgebraVectorBasis::echelon_for_space(
            self.clone(),
            vectors,
            PathAlgebraVectorTriangularization {
                echelon_basis: triangular.echelon_basis,
                heads: triangular.heads,
                base_change,
            },
        ))
    }

    /// Translate QPA's default `Basis(V)` method.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::basis_of_domain`] failures.
    pub fn basis(&self) -> Result<PathAlgebraVectorBasis<F>, ProjectivePathModuleError> {
        self.basis_of_domain()
    }

    /// Translate QPA's `vector in V` method.
    ///
    /// A different projective-vector family is outside this typed space.
    ///
    /// # Errors
    ///
    /// Propagates QPA's basis construction and coefficient scan.
    pub fn contains(
        &self,
        vector: &PathAlgebraVector<F>,
    ) -> Result<bool, ProjectivePathModuleError> {
        if !self.module.is_identical(vector.module()) {
            return Ok(false);
        }
        Ok(self.basis_of_domain()?.coefficients(vector)?.is_some())
    }

    fn ensure_ambient_family(
        &self,
        vectors: &[PathAlgebraVector<F>],
    ) -> Result<(), ProjectivePathModuleError> {
        if vectors
            .iter()
            .any(|vector| !self.module.is_identical(vector.module()))
        {
            Err(ProjectivePathModuleError::ForeignVectorBasisVector)
        } else {
            Ok(())
        }
    }
}

#[derive(Clone, Debug)]
enum PathAlgebraVectorBasisUnderlying<F: Field> {
    Space(PathAlgebraVectorSpace<F>),
    Projective(RightProjectivePathAlgebraModule<F>),
    Submodule(RightProjectivePathAlgebraSubmodule<F>),
}

type PathAlgebraVectorBaseChange<E> = Vec<Vec<(usize, E)>>;

/// QPA's `IsBasisOfPathAlgebraVectorSpace` object.
#[derive(Clone, Debug)]
pub struct PathAlgebraVectorBasis<F: Field> {
    underlying: PathAlgebraVectorBasisUnderlying<F>,
    basis_vectors: Vec<PathAlgebraVector<F>>,
    echelon_basis: Option<Vec<PathAlgebraVector<F>>>,
    heads: Option<Vec<PathAlgebraVector<F>>>,
    base_change: Option<PathAlgebraVectorBaseChange<F::Element>>,
}

impl<F: Field> PathAlgebraVectorBasis<F> {
    fn structural_for_space(
        space: PathAlgebraVectorSpace<F>,
        basis_vectors: Vec<PathAlgebraVector<F>>,
    ) -> Self {
        Self {
            underlying: PathAlgebraVectorBasisUnderlying::Space(space),
            basis_vectors,
            echelon_basis: None,
            heads: None,
            base_change: None,
        }
    }

    fn structural_for_submodule(
        submodule: RightProjectivePathAlgebraSubmodule<F>,
        basis_vectors: Vec<PathAlgebraVector<F>>,
    ) -> Self {
        Self {
            underlying: PathAlgebraVectorBasisUnderlying::Submodule(submodule),
            basis_vectors,
            echelon_basis: None,
            heads: None,
            base_change: None,
        }
    }

    fn structural_for_projective(
        module: RightProjectivePathAlgebraModule<F>,
        basis_vectors: Vec<PathAlgebraVector<F>>,
    ) -> Self {
        Self {
            underlying: PathAlgebraVectorBasisUnderlying::Projective(module),
            basis_vectors,
            echelon_basis: None,
            heads: None,
            base_change: None,
        }
    }

    fn echelon_for_space(
        space: PathAlgebraVectorSpace<F>,
        basis_vectors: Vec<PathAlgebraVector<F>>,
        triangular: PathAlgebraVectorTriangularization<F>,
    ) -> Self {
        Self {
            underlying: PathAlgebraVectorBasisUnderlying::Space(space),
            basis_vectors,
            echelon_basis: Some(triangular.echelon_basis),
            heads: Some(triangular.heads),
            base_change: Some(triangular.base_change),
        }
    }

    /// Return QPA's ordered `BasisVectors` list.
    #[must_use]
    pub fn basis_vectors(&self) -> &[PathAlgebraVector<F>] {
        &self.basis_vectors
    }

    /// Return the underlying free field-space, when this basis came from one.
    #[must_use]
    pub fn underlying_space(&self) -> Option<&PathAlgebraVectorSpace<F>> {
        match &self.underlying {
            PathAlgebraVectorBasisUnderlying::Space(space) => Some(space),
            PathAlgebraVectorBasisUnderlying::Projective(_)
            | PathAlgebraVectorBasisUnderlying::Submodule(_) => None,
        }
    }

    /// Return the underlying whole projective module, when applicable.
    #[must_use]
    pub fn underlying_projective_module(&self) -> Option<&RightProjectivePathAlgebraModule<F>> {
        match &self.underlying {
            PathAlgebraVectorBasisUnderlying::Projective(module) => Some(module),
            PathAlgebraVectorBasisUnderlying::Space(_)
            | PathAlgebraVectorBasisUnderlying::Submodule(_) => None,
        }
    }

    /// Return the underlying projective submodule, when applicable.
    #[must_use]
    pub fn underlying_submodule(&self) -> Option<&RightProjectivePathAlgebraSubmodule<F>> {
        match &self.underlying {
            PathAlgebraVectorBasisUnderlying::Space(_)
            | PathAlgebraVectorBasisUnderlying::Projective(_) => None,
            PathAlgebraVectorBasisUnderlying::Submodule(submodule) => Some(submodule),
        }
    }

    /// Translate QPA's `Coefficients(B, vector)` tip scan.
    ///
    /// `Ok(None)` is GAP's `fail`. In particular, this preserves QPA's scan
    /// stopping at a head in a larger component position, even when a later
    /// head could span the query.
    ///
    /// # Errors
    ///
    /// Structural `NewBasis` values have no echelon data. Foreign vector
    /// families and exact arithmetic/order failures are reported explicitly.
    pub fn coefficients(
        &self,
        vector: &PathAlgebraVector<F>,
    ) -> Result<Option<Vec<F::Element>>, ProjectivePathModuleError> {
        let module = match &self.underlying {
            PathAlgebraVectorBasisUnderlying::Space(space) => space.module(),
            PathAlgebraVectorBasisUnderlying::Projective(module) => module,
            PathAlgebraVectorBasisUnderlying::Submodule(submodule) => submodule.parent(),
        };
        if !module.is_identical(vector.module()) {
            return Err(ProjectivePathModuleError::ForeignVectorBasisVector);
        }
        let echelon_basis = self
            .echelon_basis
            .as_ref()
            .ok_or(ProjectivePathModuleError::MissingVectorBasisEchelonData)?;
        let heads = self
            .heads
            .as_ref()
            .ok_or(ProjectivePathModuleError::MissingVectorBasisEchelonData)?;
        let base_change = self
            .base_change
            .as_ref()
            .ok_or(ProjectivePathModuleError::MissingVectorBasisEchelonData)?;
        let field = module.algebra().field();
        let mut current = vector.clone();
        let mut coefficients = vec![field.zero(); self.basis_vectors.len()];
        let mut index = heads.len();

        while index > 0 && !current.is_zero() {
            let leading_monomial = current.leading_path()?;
            let position = current.leading_position();
            let head_monomial = loop {
                let head = &heads[index - 1];
                let head_monomial = head.leading_path()?;
                let head_position = head.leading_position();
                index -= 1;
                let monomial_not_greater = !module
                    .algebra()
                    .ordering()
                    .compare(&head_monomial, &leading_monomial)?
                    .is_gt();
                if (head_position == position && monomial_not_greater)
                    || position < head_position
                    || index == 0
                {
                    break head_monomial;
                }
            };
            if leading_monomial != head_monomial {
                return Ok(None);
            }
            let coefficient = current.leading_coefficient()?;
            current = current.sub(&echelon_basis[index].scale_left(&coefficient)?)?;
            for (basis_position, base_coefficient) in &base_change[index] {
                field.add_mul_assign(
                    &mut coefficients[*basis_position],
                    base_coefficient,
                    &coefficient,
                );
            }
        }

        Ok(current.is_zero().then_some(coefficients))
    }
}

#[derive(Clone, Debug)]
struct PathAlgebraVectorTriangularization<F: Field> {
    echelon_basis: Vec<PathAlgebraVector<F>>,
    heads: Vec<PathAlgebraVector<F>>,
    base_change: PathAlgebraVectorBaseChange<F::Element>,
}

fn triangularize_path_algebra_vector_list<F: Field>(
    entered: &[PathAlgebraVector<F>],
) -> Result<PathAlgebraVectorTriangularization<F>, ProjectivePathModuleError> {
    let mut vectors = entered
        .iter()
        .filter(|vector| !vector.is_zero())
        .cloned()
        .collect::<Vec<_>>();
    if vectors.is_empty() {
        return Ok(PathAlgebraVectorTriangularization {
            echelon_basis: Vec::new(),
            heads: Vec::new(),
            base_change: Vec::new(),
        });
    }
    let module = vectors[0].module().clone();
    let field = module.algebra().field();
    let mut base_change = (0..vectors.len())
        .map(|position| vec![(position, field.one())])
        .collect::<Vec<_>>();
    sort_parallel_path_algebra_vectors(&mut vectors, &mut base_change)?;
    let mut heads = Vec::new();
    let mut k = vectors.len();

    while k > 0 {
        let reducer = k - 1;
        let leading_coefficient = vectors[reducer].leading_coefficient()?;
        let inverse = field.div(&field.one(), &leading_coefficient);
        vectors[reducer] = vectors[reducer].scale_left(&inverse)?;
        for (_, coefficient) in &mut base_change[reducer] {
            *coefficient = field.div(coefficient, &leading_coefficient);
        }
        let head = vectors[reducer].leading_term();
        let head_position = head.leading_position();
        heads.push(head.clone());
        let mut sort = false;

        for position in (0..reducer).rev() {
            let test_head = vectors[position].leading_term();
            let same_component_head = head_position == test_head.leading_position()
                && head.components()[head_position - 1]
                    == test_head.components()[head_position - 1];
            if !same_component_head {
                continue;
            }
            sort = true;
            let coefficient = vectors[position].leading_coefficient()?;
            vectors[position] =
                vectors[position].sub(&vectors[reducer].scale_left(&coefficient)?)?;
            let reducer_change = base_change[reducer].clone();
            for (basis_position, reducer_coefficient) in reducer_change {
                let correction = field.neg(&field.mul(&coefficient, &reducer_coefficient));
                match base_change[position].binary_search_by_key(&basis_position, |entry| entry.0) {
                    Ok(existing) => {
                        field.add_assign(&mut base_change[position][existing].1, &correction);
                    }
                    Err(insertion) => {
                        base_change[position].insert(insertion, (basis_position, correction));
                    }
                }
            }
        }

        if sort {
            let old_k = k;
            let mut position = 0;
            while position < old_k - 1 {
                if vectors[position].is_zero() {
                    vectors.remove(position);
                    base_change.remove(position);
                    k -= 1;
                } else {
                    position += 1;
                }
            }
            sort_parallel_path_algebra_vectors(&mut vectors, &mut base_change)?;
        }
        k -= 1;
    }
    heads.reverse();
    Ok(PathAlgebraVectorTriangularization {
        echelon_basis: vectors,
        heads,
        base_change,
    })
}

fn sort_parallel_path_algebra_vectors<F: Field>(
    vectors: &mut Vec<PathAlgebraVector<F>>,
    parallel: &mut Vec<Vec<(usize, F::Element)>>,
) -> Result<(), ProjectivePathModuleError> {
    let entered_vectors = std::mem::take(vectors);
    let entered_parallel = std::mem::take(parallel);
    for (vector, data) in entered_vectors.into_iter().zip(entered_parallel) {
        let mut insertion = vectors.len();
        for (position, existing) in vectors.iter().enumerate() {
            if vector.compare(existing)?.is_lt() {
                insertion = position;
                break;
            }
        }
        vectors.insert(insertion, vector);
        parallel.insert(insertion, data);
    }
    Ok(())
}

fn uniform_generators_of_module<F: Field>(
    module: &RightProjectivePathAlgebraModule<F>,
    generators: &[PathAlgebraVector<F>],
) -> Result<Vec<PathAlgebraVector<F>>, ProjectivePathModuleError> {
    let mut uniform = Vec::new();
    for vertex in module.algebra().quiver().vertices() {
        let idempotent = module.algebra().element(vertex.path())?;
        for generator in generators {
            let component = generator.act(&idempotent)?;
            if !component.is_zero() {
                add_set_path_algebra_vector(&mut uniform, component)?;
            }
        }
    }
    Ok(uniform)
}

fn sort_path_algebra_vectors<F: Field>(
    vectors: &mut Vec<PathAlgebraVector<F>>,
) -> Result<(), ProjectivePathModuleError> {
    let entered = std::mem::take(vectors);
    for vector in entered {
        insert_path_algebra_vector(vectors, vector, false)?;
    }
    Ok(())
}

fn add_set_path_algebra_vector<F: Field>(
    vectors: &mut Vec<PathAlgebraVector<F>>,
    vector: PathAlgebraVector<F>,
) -> Result<(), ProjectivePathModuleError> {
    insert_path_algebra_vector(vectors, vector, true)
}

fn insert_path_algebra_vector<F: Field>(
    vectors: &mut Vec<PathAlgebraVector<F>>,
    vector: PathAlgebraVector<F>,
    remove_equal: bool,
) -> Result<(), ProjectivePathModuleError> {
    for (position, existing) in vectors.iter().enumerate() {
        if remove_equal && existing == &vector {
            return Ok(());
        }
        if vector.compare(existing)?.is_lt() {
            vectors.insert(position, vector);
            return Ok(());
        }
    }
    vectors.push(vector);
    Ok(())
}

/// One QPA O-set entry returned by `XSetOfPathAlgebraVector`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PathAlgebraVectorOSetEntry<F: Field> {
    /// Tip monomial of the selected right Gröbner-basis relation.
    pub right_tip: Path,
    /// One-based position of the selected `f` generator.
    pub generator_position: usize,
    /// QPA's coefficient-bearing prefix `q'`.
    pub q_prime: PathAlgebraElement<F>,
    /// Selected relation from the attached two-sided Gröbner basis.
    pub relation: PathAlgebraElement<F>,
}

/// One QPA N-set entry returned by `XSetOfPathAlgebraVector`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PathAlgebraVectorNSetEntry<F: Field> {
    /// Tip monomial of the selected right Gröbner-basis relation.
    pub right_tip: Path,
    /// QPA's intervening path `z`.
    pub z: PathAlgebraElement<F>,
    /// Selected relation from the attached two-sided Gröbner basis.
    pub relation: PathAlgebraElement<F>,
}

/// QPA's ordered O-set and N-set pair.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PathAlgebraVectorXSet<F: Field> {
    /// Green--Solberg--Zacharia overlap entries.
    pub o_set: Vec<PathAlgebraVectorOSetEntry<F>>,
    /// Green--Solberg--Zacharia nonoverlap entries.
    pub n_set: Vec<PathAlgebraVectorNSetEntry<F>>,
}

/// Translate QPA's `FirstPart(f, fprime, v)` scan.
///
/// Reducers are concatenated in entered order. Each complete pass tries every
/// nonzero reducer against the current tip, and successful divisions are
/// accumulated. Only the coefficients belonging to `f` are returned.
///
/// # Errors
///
/// Propagates exact vector division and arithmetic failures. QPA loops forever
/// when the entered reducers do not span the vector; Rust reports that stalled
/// source branch explicitly.
pub fn first_part<F: Field>(
    f: &[PathAlgebraVector<F>],
    f_prime: &[PathAlgebraVector<F>],
    vector: &PathAlgebraVector<F>,
) -> Result<Vec<PathAlgebraElement<F>>, ProjectivePathModuleError> {
    let algebra = vector.module().algebra();
    let reducers = f.iter().chain(f_prime).collect::<Vec<_>>();
    let mut coefficients = vec![algebra.zero(); reducers.len()];
    let mut current = vector.clone();
    while !current.is_zero() {
        let before = current.clone();
        let mut reduced = false;
        for (position, reducer) in reducers.iter().enumerate() {
            if current.is_zero() || reducer.is_zero() {
                continue;
            }
            if let Some(factor) = left_division(&current, reducer)? {
                coefficients[position] = coefficients[position].add(&factor)?;
                current = current.sub(&reducer.act(&factor)?)?;
                reduced = true;
            }
        }
        if !reduced || current == before {
            return Err(ProjectivePathModuleError::FirstPartNotSpanned);
        }
    }
    coefficients.truncate(f.len());
    Ok(coefficients)
}

/// Translate QPA's `XSetOfPathAlgebraVector(f, fprime, I, fn)`.
///
/// This retains the source's right-basis pruning, maximal coefficient-tip
/// choice, first suffix relation, and entered relation order exactly.
///
/// # Errors
///
/// Requires an attached finite-nontip Gröbner basis for `ideal`, a nonzero
/// `fn` spanned by the entered vectors, and the identical parent path algebra.
/// All exact Gröbner, path, field, and vector failures are propagated.
#[allow(clippy::too_many_lines)]
pub fn x_set_of_path_algebra_vector<F: Field>(
    f: &[PathAlgebraVector<F>],
    f_prime: &[PathAlgebraVector<F>],
    ideal: &PathAlgebraIdeal<F>,
    f_n: &PathAlgebraVector<F>,
) -> Result<PathAlgebraVectorXSet<F>, ProjectivePathModuleError> {
    let algebra = f_n.module().algebra();
    if !ideal.algebra().is_identical(algebra) {
        return Err(ProjectivePathModuleError::ForeignXSetIdeal);
    }
    if f_n.is_zero() {
        return Err(ProjectivePathModuleError::ZeroXSetVector);
    }

    let right_basis = ideal.right_groebner_basis()?;
    let right_relations = right_basis.relations();
    if right_relations.is_empty() {
        return Err(ProjectivePathModuleError::EmptyXSetRightBasis);
    }
    let mut right_tips = right_relations
        .iter()
        .enumerate()
        .map(|(position, relation)| (relation.leading_monomial(), position))
        .collect::<Vec<_>>();
    let mut divided = vec![false; right_tips.len()];
    for left in 0..right_tips.len() {
        let left_walk = path_walk(&right_tips[left].0)?;
        for right in 0..right_tips.len() {
            if left != right
                && !divided[right]
                && path_walk(&right_tips[right].0)?.starts_with(left_walk)
            {
                divided[right] = true;
            }
        }
    }
    right_tips = right_tips
        .into_iter()
        .enumerate()
        .filter_map(|(position, tip)| (!divided[position]).then_some(tip))
        .collect();

    let f_n_tip = f_n.leading_path()?;
    let f_n_walk = path_walk(&f_n_tip)?;
    let f_n_coefficient = f_n.leading_coefficient()?;
    let coefficients = first_part(f, f_prime, f_n)?;
    let mut selected: Option<(usize, Path)> = None;
    for (position, generator) in f.iter().enumerate() {
        if generator.is_zero() || !f_n_walk.starts_with(path_walk(&generator.leading_path()?)?) {
            continue;
        }
        let coefficient = &coefficients[position];
        if coefficient.is_zero() {
            continue;
        }
        let tip = coefficient.leading_monomial();
        let replace = match &selected {
            None => true,
            Some((_, greatest)) => algebra.ordering().compare(&tip, greatest)?.is_gt(),
        };
        if replace {
            selected = Some((position, tip));
        }
    }
    let (generator_position, reduced_tip) =
        selected.ok_or(ProjectivePathModuleError::MissingFirstPartTip)?;
    let reduced_walk = path_walk(&reduced_tip)?.to_vec();
    let reduced_length = reduced_walk.len();
    let mut test_set = Vec::new();
    for tip in right_tips {
        if path_walk(&tip.0)?.starts_with(&reduced_walk) {
            test_set.push(tip.0);
        }
    }

    let two_sided_basis = ideal
        .groebner_basis()
        .ok_or(ProjectivePathModuleError::MissingXSetGroebnerBasis)?;
    let two_sided_relations = two_sided_basis.relations();
    let two_sided_walks = two_sided_relations
        .iter()
        .map(|relation| path_walk(&relation.leading_monomial()).map(<[_]>::to_vec))
        .collect::<Result<Vec<_>, _>>()?;

    let mut o_set = Vec::new();
    let mut n_set = Vec::new();
    for right_tip in test_set {
        let right_walk = path_walk(&right_tip)?;
        let right_length = right_walk.len();
        let relation_position = two_sided_walks
            .iter()
            .position(|walk| right_walk.ends_with(walk))
            .ok_or(ProjectivePathModuleError::MissingXSetSuffixRelation)?;
        let relation = &two_sided_relations[relation_position];
        let relation_length = two_sided_walks[relation_position].len();
        let prefix_length = right_length - relation_length;

        if reduced_length <= prefix_length {
            let z = if reduced_length == prefix_length {
                algebra.one()
            } else {
                path_segment_element(algebra, &right_tip, reduced_length, prefix_length)?
            };
            n_set.push(PathAlgebraVectorNSetEntry {
                right_tip,
                z,
                relation: relation.clone(),
            });
        } else {
            let coefficient = algebra
                .field()
                .div(&f_n_coefficient, &relation.leading_coefficient());
            let q_prime = if prefix_length == 0 {
                algebra.one()
            } else {
                path_segment_element(algebra, &right_tip, 0, prefix_length)?
            }
            .scale_left(&coefficient)?;
            o_set.push(PathAlgebraVectorOSetEntry {
                right_tip,
                generator_position: generator_position + 1,
                q_prime,
                relation: relation.clone(),
            });
        }
    }

    Ok(PathAlgebraVectorXSet { o_set, n_set })
}

fn path_segment_element<F: Field>(
    algebra: &PathAlgebra<F>,
    path: &Path,
    start: usize,
    end: usize,
) -> Result<PathAlgebraElement<F>, ProjectivePathModuleError> {
    let walk = path_walk(path)?;
    let arrows = walk[start..end]
        .iter()
        .map(|index| {
            algebra
                .quiver()
                .arrow(*index)
                .ok_or(ProjectivePathModuleError::MissingPathArrow)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let segment = Path::from_arrow_product(algebra.quiver(), &arrows)?
        .ok_or(ProjectivePathModuleError::EmptyRightFactor)?;
    Ok(algebra.element(segment)?)
}

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

/// Translate QPA's `CompletelyReduceGroebnerBasisForModule(GB)`.
///
/// QPA first tip-reduces the flattened basis. It then preserves each leading
/// term and completely reduces only that relation's tail against the original
/// entered Gröbner-basis object.
///
/// # Errors
///
/// Propagates exact tip reduction, vector arithmetic, and complete-reduction
/// failures.
pub fn completely_reduce_path_algebra_module_groebner_basis<F: Field>(
    basis: &RightPathAlgebraModuleGroebnerBasis<F>,
) -> Result<Vec<PathAlgebraVector<F>>, ProjectivePathModuleError> {
    let mut reduced = tip_reduce_path_algebra_vectors(basis.basis_vectors())?;
    for relation in &mut reduced {
        let leading_term = relation.leading_term();
        let tail = relation.sub(&leading_term)?;
        *relation = leading_term.add(&basis.completely_reduce(&tail)?)?;
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
    /// Submodule generators must belong to their ambient vector family.
    #[error("a submodule generator belongs to another projective-vector family")]
    ForeignSubmoduleGenerator,
    /// A generated field-space has one fixed projective-vector family.
    #[error("a path-vector-space generator belongs to another projective-vector family")]
    ForeignVectorSpaceGenerator,
    /// Basis construction and coordinates use one fixed vector family.
    #[error("a path-vector basis received a vector from another projective-vector family")]
    ForeignVectorBasisVector,
    /// QPA's structural `NewBasis` does not install coefficient data.
    #[error("the path-vector basis has no echelon data for coefficients")]
    MissingVectorBasisEchelonData,
    /// GAP's empty list misses QPA's element-collection method dispatch.
    #[error("QPA does not dispatch its projective SubAlgebraModule method on an empty list")]
    EmptySubmoduleGenerators,
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
    /// QPA's `FirstPart` source loop stalls when the vector is not spanned.
    #[error("the entered f and f-prime vectors do not reduce the FirstPart vector to zero")]
    FirstPartNotSpanned,
    /// The X-set ideal and projective vector must use the identical algebra.
    #[error("the X-set ideal belongs to another path algebra")]
    ForeignXSetIdeal,
    /// QPA indexes the leading path of its X-set vector unconditionally.
    #[error("XSetOfPathAlgebraVector does not accept a zero vector")]
    ZeroXSetVector,
    /// QPA indexes the first relation of its right basis unconditionally.
    #[error("the X-set right Groebner basis is empty")]
    EmptyXSetRightBasis,
    /// No nonzero first-part coefficient supplied QPA's maximal reduced tip.
    #[error("XSetOfPathAlgebraVector found no usable FirstPart coefficient tip")]
    MissingFirstPartTip,
    /// The ideal lost the attached two-sided basis required after right-basis construction.
    #[error("the X-set ideal has no attached two-sided Groebner basis")]
    MissingXSetGroebnerBasis,
    /// QPA's unbounded suffix scan found no originating two-sided relation.
    #[error("a right Groebner-basis tip has no two-sided relation suffix")]
    MissingXSetSuffixRelation,
    /// Parent path-algebra arithmetic failed.
    #[error(transparent)]
    Algebra(#[from] PathAlgebraError),
    /// Path construction failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// QPA path ordering failed.
    #[error(transparent)]
    Ordering(#[from] crate::OrderingError),
    /// Ideal-level Gröbner construction or reduction failed.
    #[error(transparent)]
    Quotient(#[from] crate::QuotientError),
    /// QPA's module Gröbner basis uses quiver static dictionaries.
    #[error(transparent)]
    StaticDictionary(#[from] StaticDictionaryError),
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
        ArrowSpec, PathAlgebra, PathAlgebraElement, PathAlgebraIdeal, Quiver, VertexIndex,
    };

    use super::{
        PathAlgebraVector, PathAlgebraVectorNSetEntry, PathAlgebraVectorOSetEntry,
        PathAlgebraVectorSpace, ProjectivePathModuleError, first_part, left_division,
        right_projective_path_algebra_module, tip_reduce_path_algebra_vector,
        tip_reduce_path_algebra_vectors, x_set_of_path_algebra_vector,
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

    fn o_set_json<F: Field>(
        entries: &[PathAlgebraVectorOSetEntry<F>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            entries
                .iter()
                .map(|entry| {
                    json!({
                        "right_tip": entry.right_tip.external_representation(),
                        "generator_position": entry.generator_position,
                        "q_prime": element_json(&entry.q_prime, coefficient),
                        "relation": element_json(&entry.relation, coefficient),
                    })
                })
                .collect(),
        )
    }

    fn n_set_json<F: Field>(
        entries: &[PathAlgebraVectorNSetEntry<F>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            entries
                .iter()
                .map(|entry| {
                    json!({
                        "right_tip": entry.right_tip.external_representation(),
                        "z": element_json(&entry.z, coefficient),
                        "relation": element_json(&entry.relation, coefficient),
                    })
                })
                .collect(),
        )
    }

    fn moduleprojres_xset_firstpart_case<F: Field>(
        field: &F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let algebra = PathAlgebra::new(
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
        let generators = algebra.generators();
        let (vertex_element, arrow_x, arrow_y) = (&generators[0], &generators[1], &generators[2]);
        let relations = vec![
            arrow_x.multiply(arrow_x).unwrap(),
            arrow_x.multiply(arrow_y).unwrap(),
            arrow_y.multiply(arrow_y).unwrap(),
        ];
        let ideal = PathAlgebraIdeal::new(algebra.clone(), relations.clone())
            .unwrap()
            .with_complete_groebner_basis(&relations)
            .unwrap();
        let projective =
            right_projective_path_algebra_module(&algebra, std::slice::from_ref(vertex_element))
                .unwrap();
        let generator = projective.generators().unwrap().into_iter().next().unwrap();
        let f_x = projective.vector(vec![arrow_x.clone()]).unwrap();
        let f_y = projective.vector(vec![arrow_y.clone()]).unwrap();
        let mixed = projective
            .vector(vec![
                arrow_y
                    .multiply(arrow_x)
                    .unwrap()
                    .add(&arrow_x.multiply(arrow_y).unwrap())
                    .unwrap(),
            ])
            .unwrap();
        let prime_only = projective
            .vector(vec![arrow_y.multiply(arrow_x).unwrap()])
            .unwrap();
        let mixed_first = first_part(
            std::slice::from_ref(&f_x),
            std::slice::from_ref(&f_y),
            &mixed,
        )
        .unwrap();
        let prime_first = first_part(
            std::slice::from_ref(&f_x),
            std::slice::from_ref(&f_y),
            &prime_only,
        )
        .unwrap();
        let f_n = projective.vector(vec![arrow_y.neg().unwrap()]).unwrap();
        let x_set =
            x_set_of_path_algebra_vector(std::slice::from_ref(&generator), &[], &ideal, &f_n)
                .unwrap();
        let vertex_x_set =
            x_set_of_path_algebra_vector(std::slice::from_ref(&generator), &[], &ideal, &generator)
                .unwrap();
        let right_relations = ideal.right_groebner_basis().unwrap();
        let two_sided_relations = ideal.groebner_basis().unwrap();
        json!({
            "first_part_mixed": element_list_json(&mixed_first, coefficient),
            "first_part_prime_only": element_list_json(&prime_first, coefficient),
            "two_sided_relations": element_list_json(two_sided_relations.relations(), coefficient),
            "right_relations": element_list_json(right_relations.relations(), coefficient),
            "o_set": o_set_json(&x_set.o_set, coefficient),
            "n_set": n_set_json(&x_set.n_set, coefficient),
            "vertex_o_set": o_set_json(&vertex_x_set.o_set, coefficient),
            "vertex_n_set": n_set_json(&vertex_x_set.n_set, coefficient),
        })
    }

    #[test]
    fn pinned_qpa_moduleprojres_xset_firstpart_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/moduleprojres_xset_firstpart.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "moduleprojres_xset_firstpart",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": moduleprojres_xset_firstpart_case(&Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": moduleprojres_xset_firstpart_case(&field5, |entry| json!(
                field5_coefficients
                    .to_integer(entry)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
            )),
            "rational": moduleprojres_xset_firstpart_case(&Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
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

    #[allow(clippy::too_many_lines)]
    fn projective_path_basis_case<F: Field>(
        field: &F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let algebra = PathAlgebra::new(
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
        let generators = algebra.generators();
        let (vertex_v, arrow_a, arrow_b) = (&generators[0], &generators[3], &generators[4]);
        let zero = algebra.zero();
        let projective =
            right_projective_path_algebra_module(&algebra, &[vertex_v.clone(), vertex_v.clone()])
                .unwrap();
        let ab = arrow_a.multiply(arrow_b).unwrap();
        let u = projective
            .vector(vec![arrow_a.clone(), zero.clone()])
            .unwrap();
        let two = field.add(&field.one(), &field.one());
        let v = projective
            .vector(vec![
                arrow_a.scale_left(&two).unwrap().add(vertex_v).unwrap(),
                zero.clone(),
            ])
            .unwrap();
        let w = projective.vector(vec![zero.clone(), ab.clone()]).unwrap();
        let raw = projective
            .vector(vec![zero.clone(), vertex_v.clone()])
            .unwrap();
        let entered = vec![projective.zero(), v.clone(), u.clone(), w.clone()];
        let space =
            PathAlgebraVectorSpace::by_generators(projective.clone(), entered.clone()).unwrap();
        let new_basis = space.new_basis(entered.clone()).unwrap();
        let unchecked = space.basis_nc(entered.clone()).unwrap().unwrap();
        let checked = space.basis_from_vectors(entered).unwrap();
        let domain_basis = space.basis_of_domain().unwrap();
        let dependent = space
            .basis_from_vectors(vec![u.clone(), u.clone()])
            .unwrap();
        let outside =
            PathAlgebraVectorSpace::by_generators(projective.clone(), vec![u.clone()]).unwrap();
        let probe = v.add(&w).unwrap();
        let nc = projective.vector_nc(raw.components().to_vec(), 2, false);
        assert!(new_basis.underlying_space().is_some());
        assert!(matches!(
            new_basis.coefficients(&probe),
            Err(ProjectivePathModuleError::MissingVectorBasisEchelonData)
        ));
        let projective_new_basis = projective.new_basis(vec![u.clone()]).unwrap();
        assert!(
            projective_new_basis
                .underlying_projective_module()
                .unwrap()
                .is_identical(&projective)
        );

        let g1 = projective
            .vector(vec![arrow_a.add(&ab).unwrap(), arrow_a.clone()])
            .unwrap();
        let g2 = projective.vector(vec![ab, zero]).unwrap();
        let submodule = projective.submodule(vec![g1.clone(), g2.clone()]).unwrap();
        let submodule_basis = submodule.basis().unwrap().unwrap();
        let based_submodule = projective.submodule_with_basis(vec![g1, g2]).unwrap();
        let based_basis = based_submodule.basis().unwrap().unwrap();
        assert!(
            based_basis
                .underlying_submodule()
                .unwrap()
                .is_identical(&based_submodule)
        );
        let coefficient_list = |entries: Option<Vec<F::Element>>| match entries {
            Some(entries) => Value::Array(entries.iter().map(coefficient).collect()),
            None => Value::Null,
        };
        json!({
            "new_basis": vector_list_json(new_basis.basis_vectors(), coefficient),
            "unchecked_basis": vector_list_json(unchecked.basis_vectors(), coefficient),
            "checked_basis": checked
                .as_ref()
                .map_or(Value::Null, |basis| vector_list_json(basis.basis_vectors(), coefficient)),
            "domain_basis": vector_list_json(domain_basis.basis_vectors(), coefficient),
            "probe_coefficients": coefficient_list(domain_basis.coefficients(&probe).unwrap()),
            "outside_coefficients": coefficient_list(domain_basis.coefficients(&raw).unwrap()),
            "dependent_failed": dependent.is_none(),
            "membership": [
                space.contains(&u).unwrap(),
                space.contains(&raw).unwrap(),
                space.contains(&projective.zero()).unwrap(),
                outside.contains(&raw).unwrap(),
            ],
            "nc": {
                "vector": vector_json(&nc, coefficient),
                "leading_position": nc.leading_position(),
                "is_normal_form": nc.is_normal_form(),
            },
            "submodule_basis": vector_list_json(submodule_basis.basis_vectors(), coefficient),
            "submodule_dimension": submodule.dimension().unwrap().unwrap(),
            "based_basis": vector_list_json(based_basis.basis_vectors(), coefficient),
            "based_dimension": based_submodule.dimension().unwrap().unwrap(),
        })
    }

    #[allow(clippy::too_many_lines)]
    fn projective_path_groebner_case<F: Field>(
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
        let ac = arrow_a.multiply(arrow_c).unwrap();
        let bc = arrow_b.multiply(arrow_c).unwrap();
        let a_plus_ac = arrow_a.add(&ac).unwrap();
        let g1 = projective
            .vector(vec![a_plus_ac.clone(), arrow_b.clone(), zero.clone()])
            .unwrap();
        let g2 = projective
            .vector(vec![bc.add(&ac).unwrap(), zero.clone(), arrow_c.clone()])
            .unwrap();
        let g3 = projective
            .vector(vec![a_plus_ac, zero.clone(), zero.clone()])
            .unwrap();
        let submodule = projective
            .submodule(vec![g1.clone(), g2.clone(), g3.clone()])
            .unwrap();
        let zero_submodule = projective.submodule(vec![projective.zero()]).unwrap();
        let groebner_basis = submodule.right_groebner_basis().unwrap();
        assert!(
            groebner_basis
                .underlying_submodule()
                .unwrap()
                .is_identical(&submodule)
        );

        let probe1 = projective
            .vector(vec![ac.add(arrow_a).unwrap(), zero.clone(), zero.clone()])
            .unwrap();
        let probe2 = projective
            .vector(vec![
                bc.add(&ac).unwrap().add(arrow_a).unwrap(),
                zero.clone(),
                arrow_c.clone(),
            ])
            .unwrap();
        let probe3 = projective.vector(vec![bc, ac, zero]).unwrap();
        json!({
            "projective_dimension": projective.dimension().unwrap(),
            "projective_basis": vector_list_json(
                &projective.canonical_basis().unwrap().unwrap(),
                coefficient,
            ),
            "whole_uniform_generators": vector_list_json(
                &projective.uniform_generators().unwrap(),
                coefficient,
            ),
            "entered_generators": vector_list_json(submodule.generators(), coefficient),
            "uniform_generators": vector_list_json(
                &submodule.uniform_generators().unwrap(),
                coefficient,
            ),
            "zero_uniform_generators": vector_list_json(
                &zero_submodule.uniform_generators().unwrap(),
                coefficient,
            ),
            "groebner_basis": vector_list_json(groebner_basis.basis_vectors(), coefficient),
            "reduced_probes": vector_list_json(
                &[
                    groebner_basis.completely_reduce(&probe1).unwrap(),
                    groebner_basis.completely_reduce(&probe2).unwrap(),
                    groebner_basis.completely_reduce(&probe3).unwrap(),
                ],
                coefficient,
            ),
        })
    }

    #[test]
    fn pinned_qpa_projective_path_module_groebner_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/projpathalgmodule_groebner.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "projpathalgmodule_groebner",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": projective_path_groebner_case(&Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": projective_path_groebner_case(&field5, |entry| json!(
                field5_coefficients
                    .to_integer(entry)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
            )),
            "rational": projective_path_groebner_case(&Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn pinned_qpa_projective_path_module_basis_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/projpathalgmodule_basis.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "projpathalgmodule_basis",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": projective_path_basis_case(&Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": projective_path_basis_case(&field5, |entry| json!(
                field5_coefficients
                    .to_integer(entry)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
            )),
            "rational": projective_path_basis_case(&Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
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
        assert!(matches!(
            first_part(&[], &[], &first_vector),
            Err(ProjectivePathModuleError::FirstPartNotSpanned)
        ));
        assert!(first_part(&[], &[], &first.zero()).unwrap().is_empty());
        let unattached_ideal =
            PathAlgebraIdeal::new(algebra.clone(), vec![generators[1].clone()]).unwrap();
        assert!(matches!(
            x_set_of_path_algebra_vector(
                std::slice::from_ref(&first_vector),
                &[],
                &unattached_ideal,
                &first_vector,
            ),
            Err(ProjectivePathModuleError::Quotient(
                crate::QuotientError::MissingGroebnerBasis
            ))
        ));
        assert!(matches!(
            x_set_of_path_algebra_vector(
                std::slice::from_ref(&first_vector),
                &[],
                &unattached_ideal,
                &first.zero(),
            ),
            Err(ProjectivePathModuleError::ZeroXSetVector)
        ));
        assert!(matches!(
            first.submodule(vec![second_vector.clone()]),
            Err(ProjectivePathModuleError::ForeignSubmoduleGenerator)
        ));
        assert!(matches!(
            first.submodule(Vec::new()),
            Err(ProjectivePathModuleError::EmptySubmoduleGenerators)
        ));
        let whole_basis = first.right_groebner_basis().unwrap();
        assert!(whole_basis.is_for_whole_module());
        assert!(
            whole_basis
                .completely_reduce(&first_vector)
                .unwrap()
                .is_zero()
        );
        assert!(matches!(
            whole_basis.completely_reduce(&second_vector),
            Err(ProjectivePathModuleError::DifferentVectorFamilies)
        ));

        let cyclic = PathAlgebra::new(
            Z2,
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let cyclic_projective =
            right_projective_path_algebra_module(&cyclic, &[cyclic.generators()[0].clone()])
                .unwrap();
        assert!(cyclic_projective.canonical_basis().unwrap().is_none());
    }
}
