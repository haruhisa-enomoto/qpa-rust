//! QPA's uniform-element overload of `RightProjectiveModule`.
//!
//! This older overload stores the entered uniform elements in its vector
//! family but does not install the `vertices` component expected by QPA's
//! checked vector constructor. The resulting module has usable unchecked
//! generators, zero, negation, and nonzero scalar multiples, while every
//! checked nonzero result reaches the pinned missing-`vertices` error.

use std::{cmp::Ordering, fmt, sync::Arc};

use numerica::domains::Field;
use thiserror::Error;

use crate::projective_path_module::qpa_path_algebra_components_string;
use crate::{Path, PathAlgebra, PathAlgebraElement, PathEndpoint};

struct RightUniformPathAlgebraModuleData<F: Field> {
    algebra: PathAlgebra<F>,
    elements: Box<[PathAlgebraElement<F>]>,
}

/// QPA's `RightProjectiveModule(KQ, uniform_elements)` overload.
///
/// The name distinguishes this observably broken source overload from QPA's
/// working coefficient-one vertex overload.
#[derive(Clone)]
pub struct RightUniformPathAlgebraModule<F: Field>(Arc<RightUniformPathAlgebraModuleData<F>>);

impl<F: Field> RightUniformPathAlgebraModule<F> {
    /// Return the identical parent path algebra.
    #[must_use]
    pub fn algebra(&self) -> &PathAlgebra<F> {
        &self.0.algebra
    }

    /// Return QPA's immutable `elements` family component.
    #[must_use]
    pub fn elements(&self) -> &[PathAlgebraElement<F>] {
        &self.0.elements
    }

    /// Return whether two handles denote the same QPA module/vector family.
    #[must_use]
    pub fn is_identical(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Return QPA's unchecked module zero with cached position `len + 1`.
    #[must_use]
    pub fn zero(&self) -> RightUniformPathAlgebraVector<F> {
        RightUniformPathAlgebraVector::new_unchecked(
            self.clone(),
            vec![self.algebra().zero(); self.elements().len()],
            self.elements().len() + 1,
        )
    }

    /// Return QPA's unchecked coordinate generators containing the entered
    /// uniform elements themselves.
    #[must_use]
    pub fn generators(&self) -> Vec<RightUniformPathAlgebraVector<F>> {
        self.elements()
            .iter()
            .cloned()
            .enumerate()
            .map(|(position, element)| {
                let mut components = vec![self.algebra().zero(); self.elements().len()];
                components[position] = element;
                RightUniformPathAlgebraVector::new_unchecked(self.clone(), components, position + 1)
            })
            .collect()
    }

    /// Translate QPA's checked `PathAlgebraVector` constructor for this family.
    ///
    /// # Errors
    ///
    /// Checks component family, length, and left uniformity before preserving
    /// QPA's unbound `vertices` access for the highest nonzero coordinate.
    pub fn vector(
        &self,
        components: Vec<PathAlgebraElement<F>>,
    ) -> Result<RightUniformPathAlgebraVector<F>, UniformProjectiveModuleError> {
        if components.len() != self.elements().len() {
            return Err(UniformProjectiveModuleError::ComponentCount {
                actual: components.len(),
                expected: self.elements().len(),
            });
        }
        if components
            .iter()
            .any(|component| !component.algebra().is_identical(self.algebra()))
        {
            return Err(UniformProjectiveModuleError::ForeignComponentAlgebra);
        }
        for (position, component) in components.iter().enumerate().rev() {
            if !component.is_left_uniform() {
                return Err(UniformProjectiveModuleError::ComponentNotLeftUniform {
                    position: position + 1,
                });
            }
            if !component.is_zero() {
                return Err(UniformProjectiveModuleError::MissingVertices {
                    position: position + 1,
                });
            }
        }
        Ok(RightUniformPathAlgebraVector::new_unchecked(
            self.clone(),
            components,
            1,
        ))
    }

    /// Translate QPA's nonchecking positional constructor.
    #[must_use]
    pub fn vector_nc(
        &self,
        components: Vec<PathAlgebraElement<F>>,
        leading_position: usize,
        _normal: bool,
    ) -> RightUniformPathAlgebraVector<F> {
        RightUniformPathAlgebraVector::new_unchecked(self.clone(), components, leading_position)
    }

    /// Return QPA's constructor-installed vertex-projective property.
    #[must_use]
    pub const fn is_vertex_projective_module(&self) -> bool {
        true
    }

    /// Translate the canonical-basis method's missing-`vertices` boundary.
    ///
    /// `None` remains QPA's unavailable canonical algebra basis for a cyclic
    /// quiver; an empty-rank module has the empty basis.
    ///
    /// # Errors
    ///
    /// Every nonempty finite-rank family lacks QPA's required `vertices`
    /// component.
    pub fn canonical_basis(
        &self,
    ) -> Result<Option<Vec<RightUniformPathAlgebraVector<F>>>, UniformProjectiveModuleError> {
        if self.algebra().canonical_basis().is_none() {
            return Ok(None);
        }
        if self.elements().is_empty() {
            Ok(Some(Vec::new()))
        } else {
            Err(UniformProjectiveModuleError::MissingVertices { position: 1 })
        }
    }

    /// Translate QPA's vertex-major uniform-generator loop.
    ///
    /// # Errors
    ///
    /// For every nonzero entered uniform element, its target vertex produces
    /// a nonzero checked action and therefore QPA's missing-`vertices` error.
    pub fn uniform_generators(
        &self,
    ) -> Result<Vec<RightUniformPathAlgebraVector<F>>, UniformProjectiveModuleError> {
        let generators = self.generators();
        let mut result = Vec::new();
        for vertex in self.algebra().quiver().vertices() {
            let embedded = self.algebra().element(vertex.path())?;
            for generator in &generators {
                let acted = generator.act(&embedded)?;
                if !acted.is_zero() && !result.contains(&acted) {
                    result.push(acted);
                }
            }
        }
        Ok(result)
    }
}

impl<F: Field> fmt::Debug for RightUniformPathAlgebraModule<F> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("RightUniformPathAlgebraModule")
            .field("algebra", self.algebra())
            .field("elements", &self.elements())
            .finish()
    }
}

impl<F: Field> PartialEq for RightUniformPathAlgebraModule<F> {
    fn eq(&self, other: &Self) -> bool {
        self.is_identical(other)
    }
}

impl<F: Field> Eq for RightUniformPathAlgebraModule<F> {}

/// QPA's positional vector in a uniform-element projective family.
#[derive(Clone, Debug)]
pub struct RightUniformPathAlgebraVector<F: Field> {
    module: RightUniformPathAlgebraModule<F>,
    components: Vec<PathAlgebraElement<F>>,
    leading_position: usize,
}

impl<F: Field> RightUniformPathAlgebraVector<F> {
    fn new_unchecked(
        module: RightUniformPathAlgebraModule<F>,
        components: Vec<PathAlgebraElement<F>>,
        leading_position: usize,
    ) -> Self {
        Self {
            module,
            components,
            leading_position,
        }
    }

    /// Return the identical owning module/vector family.
    #[must_use]
    pub const fn module(&self) -> &RightUniformPathAlgebraModule<F> {
        &self.module
    }

    /// Return QPA's stored component list.
    #[must_use]
    pub fn components(&self) -> &[PathAlgebraElement<F>] {
        &self.components
    }

    /// Translate QPA's generic path-vector `PrintObj` method.
    ///
    /// The caller supplies GAP-facing coefficient spelling because Numerica
    /// owns exact arithmetic values rather than GAP presentation syntax.
    #[must_use]
    pub fn print_obj_with(&self, coefficient: impl Fn(&F::Element) -> String + Copy) -> String {
        qpa_path_algebra_components_string(self.components(), coefficient)
    }

    /// Return QPA's cached one-based leading position.
    #[must_use]
    pub const fn leading_position(&self) -> usize {
        self.leading_position
    }

    /// Return the default type's installed normal-form property.
    #[must_use]
    pub const fn is_normal_form(&self) -> bool {
        true
    }

    /// Return whether every component vanishes.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.components.iter().all(PathAlgebraElement::is_zero)
    }

    /// Add identical-family vectors through QPA's checked constructor.
    ///
    /// # Errors
    ///
    /// A nonzero result reaches the missing-`vertices` source bug.
    pub fn add(&self, right: &Self) -> Result<Self, UniformProjectiveModuleError> {
        self.ensure_same_family(right)?;
        self.module.vector(
            self.components
                .iter()
                .zip(&right.components)
                .map(|(left, right)| left.add(right))
                .collect::<Result<Vec<_>, _>>()?,
        )
    }

    /// Subtract identical-family vectors through QPA's checked constructor.
    ///
    /// # Errors
    ///
    /// A nonzero result reaches the missing-`vertices` source bug.
    pub fn sub(&self, right: &Self) -> Result<Self, UniformProjectiveModuleError> {
        self.ensure_same_family(right)?;
        self.module.vector(
            self.components
                .iter()
                .zip(&right.components)
                .map(|(left, right)| left.sub(right))
                .collect::<Result<Vec<_>, _>>()?,
        )
    }

    /// Return QPA's unchecked additive inverse.
    ///
    /// # Errors
    ///
    /// Propagates component negation.
    pub fn neg(&self) -> Result<Self, UniformProjectiveModuleError> {
        Ok(Self::new_unchecked(
            self.module.clone(),
            self.components
                .iter()
                .map(PathAlgebraElement::neg)
                .collect::<Result<Vec<_>, _>>()?,
            self.leading_position,
        ))
    }

    /// Return QPA's unchecked left scalar multiple.
    ///
    /// # Errors
    ///
    /// Propagates component scaling.
    pub fn scale_left(&self, scalar: &F::Element) -> Result<Self, UniformProjectiveModuleError> {
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

    /// Return QPA's identical right scalar action.
    ///
    /// # Errors
    ///
    /// Propagates component scaling.
    pub fn scale_right(&self, scalar: &F::Element) -> Result<Self, UniformProjectiveModuleError> {
        self.scale_left(scalar)
    }

    /// Apply QPA's checked coordinatewise right path-algebra action.
    ///
    /// # Errors
    ///
    /// Rejects a foreign algebra; a nonzero result reaches the missing-family
    /// `vertices` component.
    pub fn act(
        &self,
        element: &PathAlgebraElement<F>,
    ) -> Result<Self, UniformProjectiveModuleError> {
        if !element.algebra().is_identical(self.module.algebra()) {
            return Err(UniformProjectiveModuleError::ForeignActionAlgebra);
        }
        self.module.vector(
            self.components
                .iter()
                .map(|component| component.multiply(element))
                .collect::<Result<Vec<_>, _>>()?,
        )
    }

    /// Return QPA's sparse unchecked leading term.
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
    /// QPA indexes coordinate one for zero and fails for rank zero.
    pub fn leading_coefficient(&self) -> Result<F::Element, UniformProjectiveModuleError> {
        let position = if self.leading_position <= self.components.len() {
            self.leading_position - 1
        } else {
            0
        };
        self.components
            .get(position)
            .map(PathAlgebraElement::leading_coefficient)
            .ok_or(UniformProjectiveModuleError::EmptyLeadingComponent)
    }

    /// Return the coordinate containing QPA's vector tip.
    ///
    /// # Errors
    ///
    /// QPA indexes coordinate one for zero and fails for rank zero.
    pub fn leading_component(
        &self,
    ) -> Result<&PathAlgebraElement<F>, UniformProjectiveModuleError> {
        let position = if self.is_zero() {
            0
        } else {
            self.leading_position - 1
        };
        self.components
            .get(position)
            .ok_or(UniformProjectiveModuleError::EmptyLeadingComponent)
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

    /// Translate QPA's surprising right-uniform-gated `SourceVertex`.
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

    /// Translate QPA's vector left-uniformity predicate.
    #[must_use]
    pub fn is_left_uniform(&self) -> bool {
        self.is_uniform_at(Path::source, PathAlgebraElement::is_left_uniform)
    }

    /// Translate QPA's vector right-uniformity predicate.
    #[must_use]
    pub fn is_right_uniform(&self) -> bool {
        self.is_uniform_at(Path::target, PathAlgebraElement::is_right_uniform)
    }

    /// Return whether the vector is both left and right uniform.
    #[must_use]
    pub fn is_uniform(&self) -> bool {
        self.is_left_uniform() && self.is_right_uniform()
    }

    /// Compare vectors by tip and descending component position.
    ///
    /// # Errors
    ///
    /// Rejects different vector families, zero tips, and ordering failures.
    pub fn compare(&self, right: &Self) -> Result<Ordering, UniformProjectiveModuleError> {
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
    /// Rejects different families and zero or non-walk tips.
    pub fn is_left_divisible(&self, divisor: &Self) -> Result<bool, UniformProjectiveModuleError> {
        self.ensure_same_family(divisor)?;
        if self.is_zero() || divisor.is_zero() {
            return Err(UniformProjectiveModuleError::ZeroVectorDivision);
        }
        let dividend_path = self.leading_path()?;
        let divisor_path = divisor.leading_path()?;
        let dividend_walk = dividend_path
            .walk()
            .ok_or(UniformProjectiveModuleError::PathHasNoArrowWalk)?;
        let divisor_walk = divisor_path
            .walk()
            .ok_or(UniformProjectiveModuleError::PathHasNoArrowWalk)?;
        Ok(self.leading_position == divisor.leading_position
            && !dividend_walk.is_empty()
            && dividend_walk.starts_with(divisor_walk))
    }

    fn ensure_same_family(&self, right: &Self) -> Result<(), UniformProjectiveModuleError> {
        if self.module.is_identical(&right.module) {
            Ok(())
        } else {
            Err(UniformProjectiveModuleError::DifferentVectorFamilies)
        }
    }

    fn leading_path(&self) -> Result<Path, UniformProjectiveModuleError> {
        if self.is_zero() {
            return Err(UniformProjectiveModuleError::ZeroVectorDivision);
        }
        self.components
            .get(self.leading_position - 1)
            .map(PathAlgebraElement::leading_monomial)
            .ok_or(UniformProjectiveModuleError::InvalidLeadingPosition {
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

impl<F: Field> PartialEq for RightUniformPathAlgebraVector<F> {
    fn eq(&self, other: &Self) -> bool {
        self.module.is_identical(&other.module) && self.components == other.components
    }
}

impl<F: Field> Eq for RightUniformPathAlgebraVector<F> {}

/// Construct QPA's uniform-element overload of `RightProjectiveModule`.
///
/// # Errors
///
/// Every element must lie in the identical parent algebra and be uniform.
pub fn right_uniform_path_algebra_module<F: Field>(
    algebra: &PathAlgebra<F>,
    elements: &[PathAlgebraElement<F>],
) -> Result<RightUniformPathAlgebraModule<F>, UniformProjectiveModuleError> {
    for (position, element) in elements.iter().enumerate() {
        if !element.algebra().is_identical(algebra) {
            return Err(UniformProjectiveModuleError::ForeignSummand {
                position: position + 1,
            });
        }
        if !element.is_uniform() {
            return Err(UniformProjectiveModuleError::NonuniformSummand {
                position: position + 1,
            });
        }
    }
    Ok(RightUniformPathAlgebraModule(Arc::new(
        RightUniformPathAlgebraModuleData {
            algebra: algebra.clone(),
            elements: elements.to_vec().into_boxed_slice(),
        },
    )))
}

/// Failure in QPA's uniform-element projective overload.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum UniformProjectiveModuleError {
    /// An entered summand belongs to another path algebra.
    #[error("uniform projective summand {position} belongs to another path algebra")]
    ForeignSummand {
        /// One-based list position.
        position: usize,
    },
    /// QPA's first overload requires every entered element to be uniform.
    #[error("uniform projective summand {position} is not uniform")]
    NonuniformSummand {
        /// One-based list position.
        position: usize,
    },
    /// Checked vectors use the family's fixed rank.
    #[error("uniform path vector has {actual} components, expected {expected}")]
    ComponentCount {
        /// Entered coordinate count.
        actual: usize,
        /// Required coordinate count.
        expected: usize,
    },
    /// Every component must belong to the module's parent algebra.
    #[error("a uniform path-vector component belongs to another path algebra")]
    ForeignComponentAlgebra,
    /// QPA checks component left uniformity before reading `vertices`.
    #[error("uniform path-vector component {position} is not left uniform")]
    ComponentNotLeftUniform {
        /// One-based coordinate.
        position: usize,
    },
    /// The source overload never installs the family component it reads.
    #[error(
        "QPA's uniform-element vector family has no vertices component at coordinate {position}"
    )]
    MissingVertices {
        /// First offending one-based coordinate in QPA's reverse scan.
        position: usize,
    },
    /// Arithmetic requires the identical source family.
    #[error("the uniform path vectors belong to different module families")]
    DifferentVectorFamilies,
    /// Right action requires the identical parent algebra.
    #[error("the uniform vector action belongs to another path algebra")]
    ForeignActionAlgebra,
    /// QPA indexes coordinate one for zero and fails at rank zero.
    #[error("an empty uniform path vector has no leading component")]
    EmptyLeadingComponent,
    /// An unchecked vector stored an invalid one-based leading position.
    #[error("leading position {position} is outside a uniform vector of length {length}")]
    InvalidLeadingPosition {
        /// Stored one-based position.
        position: usize,
        /// Coordinate count.
        length: usize,
    },
    /// QPA's divisibility operation has no zero-vector tip.
    #[error("uniform path-vector divisibility does not accept zero inputs")]
    ZeroVectorDivision,
    /// `WalkOfPath` is unavailable for the entered tip.
    #[error("a uniform path-vector tip has no arrow walk")]
    PathHasNoArrowWalk,
    /// Parent path-algebra arithmetic failed.
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
        RightUniformPathAlgebraVector, UniformProjectiveModuleError,
        right_uniform_path_algebra_module,
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
        vector: &RightUniformPathAlgebraVector<F>,
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
        vectors: &[RightUniformPathAlgebraVector<F>],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        Value::Array(
            vectors
                .iter()
                .map(|vector| vector_json(vector, coefficient))
                .collect(),
        )
    }

    fn uniform_projective_case<F: Field>(
        field: &F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let algebra = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap();
        let generators = algebra.generators();
        let module = right_uniform_path_algebra_module(&algebra, &[generators[2].clone()]).unwrap();
        let generator = module.generators()[0].clone();
        let source_action = generator.act(&generators[0]).unwrap();
        let negative = generator.neg().unwrap();
        let two = field.add(&field.one(), &field.one());
        let scaled_two = generator.scale_left(&two).unwrap();
        let cancelled = generator.add(&negative).unwrap();
        json!({
            "generators": vector_list_json(&module.generators(), coefficient),
            "module_zero": vector_json(&module.zero(), coefficient),
            "source_action": vector_json(&source_action, coefficient),
            "negative": vector_json(&negative, coefficient),
            "scaled_two": vector_json(&scaled_two, coefficient),
            "cancelled": vector_json(&cancelled, coefficient),
            "positions": [
                module.zero().leading_position(),
                cancelled.leading_position(),
                generator.leading_position(),
            ],
            "properties": [
                module.is_vertex_projective_module(),
                generator.is_normal_form(),
            ],
            "failure_success": [
                generator.act(&generators[1]).is_ok(),
                module.vector(vec![generators[2].clone()]).is_ok(),
                module.canonical_basis().is_ok(),
                module.uniform_generators().is_ok(),
            ],
        })
    }

    #[test]
    fn pinned_qpa_uniform_projective_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/projpathalgmodule_uniform_projective.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "projpathalgmodule_uniform_projective",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": uniform_projective_case(&Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": uniform_projective_case(&field5, |entry| json!(
                field5_coefficients
                    .to_integer(entry)
                    .to_string()
                    .parse::<i64>()
                    .unwrap()
            )),
            "rational": uniform_projective_case(&Q, |entry: &Rational| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });
        assert_eq!(actual, expected);
    }

    #[test]
    fn uniform_projective_boundaries_are_typed() {
        let algebra = PathAlgebra::new(
            Z2,
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
        )
        .unwrap();
        let generators = algebra.generators();
        let module = right_uniform_path_algebra_module(&algebra, &[generators[2].clone()]).unwrap();
        let generator = module.generators()[0].clone();
        assert!(matches!(
            generator.act(&generators[1]),
            Err(UniformProjectiveModuleError::MissingVertices { position: 1 })
        ));
        assert!(matches!(
            module.vector(Vec::new()),
            Err(UniformProjectiveModuleError::ComponentCount { .. })
        ));
        let nonuniform = generators[0].add(&generators[1]).unwrap();
        assert!(matches!(
            right_uniform_path_algebra_module(&algebra, &[nonuniform]),
            Err(UniformProjectiveModuleError::NonuniformSummand { position: 1 })
        ));
        let foreign = PathAlgebra::new(
            Z2,
            Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
        )
        .unwrap();
        assert!(matches!(
            right_uniform_path_algebra_module(&algebra, &[foreign.generators()[0].clone()]),
            Err(UniformProjectiveModuleError::ForeignSummand { position: 1 })
        ));
        assert!(matches!(
            generator.is_left_divisible(&module.zero()),
            Err(UniformProjectiveModuleError::ZeroVectorDivision)
        ));
    }
}
