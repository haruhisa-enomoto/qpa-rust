//! QPA's combinatorial substrate for tensor products of path algebras.
//!
//! This module translates the dependency-free first part of
//! `lib/pathalgtensor.gi`: product quivers, their retained decomposition,
//! inclusion and projection of paths, and `ReversePath`.  Algebra tensor
//! products are layered on these operations separately.

use std::{collections::HashMap, fmt};

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    Generator, NamedArrowSpec, Path, PathEndpoint, PathError, Quiver, QuiverAlgebra,
    QuiverAlgebraElement, QuiverAlgebraElementError, QuiverError, Vertex,
};

/// QPA's two-factor decomposition object attached to a product quiver.
#[derive(Clone)]
pub struct QuiverProductDecomposition {
    factors: [Quiver; 2],
    lookup_table: HashMap<Box<str>, [Generator; 2]>,
}

impl QuiverProductDecomposition {
    /// Return whether QPA's list-like decomposition binds this positive index.
    #[must_use]
    pub const fn is_bound(&self, index: usize) -> bool {
        index > 0 && index <= self.factors.len()
    }

    /// Extract a factor using QPA's one-based `dec[index]` convention.
    ///
    /// # Errors
    ///
    /// Returns an explicit error where GAP's list access would fail.
    pub fn factor(&self, index: usize) -> Result<&Quiver, QuiverProductError> {
        self.factors
            .get(index.wrapping_sub(1))
            .ok_or(QuiverProductError::FactorIndexOutOfRange { index })
    }
}

impl fmt::Debug for QuiverProductDecomposition {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("<quiver product decomposition>")
    }
}

impl fmt::Display for QuiverProductDecomposition {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("<quiver product decomposition>")
    }
}

/// A quiver constructed by QPA's `QuiverProduct`, together with its attached
/// decomposition attribute.
#[derive(Clone, Debug)]
pub struct QuiverProduct {
    quiver: Quiver,
    decomposition: QuiverProductDecomposition,
}

impl QuiverProduct {
    /// Return the product quiver itself.
    #[must_use]
    pub const fn quiver(&self) -> &Quiver {
        &self.quiver
    }

    /// Return QPA's attached `QuiverProductDecomposition` value.
    #[must_use]
    pub const fn decomposition(&self) -> &QuiverProductDecomposition {
        &self.decomposition
    }
}

/// Construct QPA's Cartesian product of two quivers.
///
/// Vertices are ordered as `Cartesian(Vertices(q1), Vertices(q2))`. Arrows
/// first hold the first coordinate fixed and traverse `q2`, then hold the
/// second coordinate fixed and traverse `q1`. Names are joined by `_` exactly
/// as in the source.
///
/// # Errors
///
/// Propagates QPA-compatible named-quiver construction failures.
pub fn quiver_product(q1: &Quiver, q2: &Quiver) -> Result<QuiverProduct, QuiverProductError> {
    let mut vertex_pairs = Vec::new();
    let mut vertex_names = Vec::new();
    for left in q1.vertices() {
        for right in q2.vertices() {
            let pair = [
                Generator::Vertex(left.clone()),
                Generator::Vertex(right.clone()),
            ];
            vertex_names.push(compose_name(&pair));
            vertex_pairs.push(pair);
        }
    }

    let mut arrow_pairs = Vec::new();
    for left in q1.vertices() {
        for right in q2.arrows() {
            arrow_pairs.push([Generator::Vertex(left.clone()), Generator::Arrow(right)]);
        }
    }
    for left in q1.arrows() {
        for right in q2.vertices() {
            arrow_pairs.push([Generator::Arrow(left.clone()), Generator::Vertex(right)]);
        }
    }

    let arrow_definitions = arrow_pairs
        .iter()
        .map(|pair| {
            NamedArrowSpec::named(
                compose_name(&[
                    Generator::Vertex(generator_source(&pair[0])),
                    Generator::Vertex(generator_source(&pair[1])),
                ]),
                compose_name(&[
                    Generator::Vertex(generator_target(&pair[0])),
                    Generator::Vertex(generator_target(&pair[1])),
                ]),
                compose_name(pair),
            )
        })
        .collect();

    let mut lookup_table = HashMap::new();
    for (name, pair) in vertex_names.iter().cloned().zip(vertex_pairs) {
        lookup_table.insert(name, pair);
    }
    for pair in &arrow_pairs {
        lookup_table.insert(compose_name(pair), pair.clone());
    }

    let quiver = Quiver::from_named(vertex_names, arrow_definitions)?;
    Ok(QuiverProduct {
        quiver,
        decomposition: QuiverProductDecomposition {
            factors: [q1.clone(), q2.clone()],
            lookup_table,
        },
    })
}

/// Translate QPA's two overloads of `WalkOfPathOrVertex`.
///
/// A vertex becomes a singleton generator list; a nontrivial path becomes its
/// arrow walk. QPA has no `WalkOfPath` value for its zero path.
///
/// # Errors
///
/// Returns [`PathError::ZeroHasNoWalk`] for the distinguished zero path.
pub fn walk_of_path_or_vertex(path: &Path) -> Result<Vec<Generator>, PathError> {
    if path.is_zero() {
        return Err(PathError::ZeroHasNoWalk);
    }
    if path.is_vertex() {
        let PathEndpoint::Vertex(vertex) = path.source() else {
            unreachable!("a nonzero vertex path has a vertex endpoint");
        };
        return Ok(vec![Generator::Vertex(vertex)]);
    }
    path.walk()
        .ok_or(PathError::ZeroHasNoWalk)?
        .iter()
        .map(|&index| {
            path.quiver().arrow(index).map(Generator::Arrow).ok_or(
                PathError::GeneratorIndexOutOfRange {
                    index: path.quiver().number_of_vertices() + index.get(),
                    number_of_generators: path.quiver().number_of_vertices()
                        + path.quiver().number_of_arrows(),
                },
            )
        })
        .collect()
}

/// Translate QPA's `ReversePath(path) = Product(Reversed(WalkOfPath(path)))`.
///
/// This is deliberately not the opposite path. Reversing a composable walk
/// commonly makes adjacent arrows noncomposable, in which case QPA's product
/// is the zero path and Rust returns that same value.
///
/// # Errors
///
/// Propagates QPA's unbound zero walk, empty generic product, and path-family
/// failures.
pub fn reverse_path(path: &Path) -> Result<Path, PathError> {
    let walk = path.walk().ok_or(PathError::ZeroHasNoWalk)?;
    let reversed = walk
        .iter()
        .rev()
        .map(|&index| {
            path.quiver().arrow(index).map(|arrow| arrow.path()).ok_or(
                PathError::GeneratorIndexOutOfRange {
                    index: path.quiver().number_of_vertices() + index.get(),
                    number_of_generators: path.quiver().number_of_vertices()
                        + path.quiver().number_of_arrows(),
                },
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let Some((first, rest)) = reversed.split_first() else {
        return Err(PathError::EmptyProduct);
    };
    let mut product = first.clone();
    for arrow in rest {
        product = product.multiply(arrow)?;
    }
    Ok(product)
}

/// Include two paths or vertices in QPA's selected product-quiver route.
///
/// The first path is traversed at the source of the second; the second is then
/// traversed at the target of the first. Generator lookup is name-based, as in
/// the GAP record-component implementation.
///
/// # Errors
///
/// Reports zero-path walks, missing composed record components, and path
/// multiplication failures.
pub fn include_in_product_quiver(
    factors: [&Path; 2],
    product: &QuiverProduct,
) -> Result<Path, QuiverProductError> {
    let walks = [
        walk_of_path_or_vertex(factors[0])?,
        walk_of_path_or_vertex(factors[1])?,
    ];
    let join = [
        endpoint_vertex(factors[0].target())?,
        endpoint_vertex(factors[1].source())?,
    ];

    let mut product_walk = Vec::with_capacity(walks[0].len().saturating_add(walks[1].len()));
    for generator in &walks[0] {
        product_walk.push(product_generator_path(
            product,
            &compose_name(&[generator.clone(), Generator::Vertex(join[1].clone())]),
        )?);
    }
    for generator in &walks[1] {
        product_walk.push(product_generator_path(
            product,
            &compose_name(&[Generator::Vertex(join[0].clone()), generator.clone()]),
        )?);
    }
    product_of_paths(&product_walk).map_err(Into::into)
}

/// Project a product-quiver path to one of its two retained factors.
///
/// Projection maps each product generator through QPA's name-keyed lookup
/// table and then multiplies the resulting factor generators in walk order.
///
/// # Errors
///
/// The path must belong to this product, and `index` must be one or two.
pub fn project_from_product_quiver(
    index: usize,
    path: &Path,
    product: &QuiverProduct,
) -> Result<Path, QuiverProductError> {
    let factor_position = index
        .checked_sub(1)
        .filter(|&position| position < 2)
        .ok_or(QuiverProductError::FactorIndexOutOfRange { index })?;
    if !path.quiver().is_identical(product.quiver()) {
        return Err(QuiverProductError::PathNotInProductQuiver);
    }

    let mut projected = Vec::new();
    for generator in walk_of_path_or_vertex(path)? {
        let name = generator_name(&generator);
        let pair = product
            .decomposition
            .lookup_table
            .get(name)
            .ok_or_else(|| QuiverProductError::MissingDecompositionEntry { name: name.into() })?;
        projected.push(pair[factor_position].path());
    }
    product_of_paths(&projected).map_err(Into::into)
}

/// Include a path or vertex in a path algebra or quotient using QPA's
/// name-based record-component lookup.
///
/// The source maps every entry of `WalkOfPathOrVertex(path)` to the algebra
/// component with the same printed name, then multiplies those components.
/// This intentionally permits an input from another quiver when all required
/// names exist, just as GAP record lookup does.
///
/// # Errors
///
/// Reports the zero-path walk boundary, a missing algebra record component,
/// or path/quotient arithmetic failure.
pub fn include_in_path_algebra<F: Field>(
    path: &Path,
    algebra: &QuiverAlgebra<F>,
) -> Result<QuiverAlgebraElement<F>, PathAlgebraTensorError> {
    let mut included = Vec::new();
    for generator in walk_of_path_or_vertex(path)? {
        let name = generator_name(&generator);
        let element = algebra
            .generator_element_by_name(name)?
            .ok_or_else(|| PathAlgebraTensorError::MissingAlgebraGenerator { name: name.into() })?;
        included.push(element);
    }
    product_of_algebra_elements(&included).map_err(Into::into)
}

/// Return QPA's ordered list of vertex idempotents for a path algebra or
/// quotient.
///
/// Each source-quiver vertex is passed through [`include_in_path_algebra`]
/// rather than embedded positionally, preserving record-name collision
/// behavior.
///
/// # Errors
///
/// Propagates source-translated inclusion and quotient-reduction failures.
pub fn vertices_of_path_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Vec<QuiverAlgebraElement<F>>, PathAlgebraTensorError> {
    algebra
        .quiver()
        .vertices()
        .map(|vertex| include_in_path_algebra(&vertex.path(), algebra))
        .collect()
}

/// One record returned by QPA's `PathAlgebraElementTerms` helper.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PathAlgebraElementTermRecord<F: Field> {
    monomial: Path,
    coefficient: F::Element,
}

impl<F: Field> PathAlgebraElementTermRecord<F> {
    /// Return the source record's `monom` component.
    #[must_use]
    pub const fn monomial(&self) -> &Path {
        &self.monomial
    }

    /// Return the source record's `coeff` component.
    #[must_use]
    pub const fn coefficient(&self) -> &F::Element {
        &self.coefficient
    }
}

/// Translate QPA's tip-peeling `PathAlgebraElementTerms` loop.
///
/// Terms are returned from greatest to least monomial. The implementation
/// deliberately repeats `Tip`, `TipMonomial`, `TipCoefficient`, and
/// subtraction until the tip coefficient is zero, matching the source rather
/// than merely reversing Rust's sparse storage.
///
/// # Errors
///
/// Propagates path-algebra family and quotient normal-form failures from each
/// source subtraction.
pub fn path_algebra_element_terms<F: Field>(
    element: &QuiverAlgebraElement<F>,
) -> Result<Vec<PathAlgebraElementTermRecord<F>>, PathAlgebraTensorError> {
    let mut terms = Vec::new();
    let mut rest = element.clone();
    loop {
        let tip = rest.leading_term();
        let monomial = rest.leading_monomial();
        let coefficient = rest.leading_coefficient();
        if rest.field().is_zero(&coefficient) {
            break;
        }
        rest = rest.sub(&tip)?;
        terms.push(PathAlgebraElementTermRecord {
            monomial,
            coefficient,
        });
    }
    Ok(terms)
}

/// Translate QPA's bilinear `SimpleTensor([x, y], target)` expansion.
///
/// QPA first returns target zero if either factor is zero. Otherwise it peels
/// both elements into descending tip records, traverses their Cartesian
/// product in GAP order, includes each monomial pair through the product
/// quiver, multiplies the two coefficients, and sums the resulting target
/// elements in that order.
///
/// Rust receives the retained [`QuiverProduct`] explicitly because product
/// decomposition is structured metadata rather than a dynamically attached
/// GAP attribute.
///
/// # Errors
///
/// Nonzero factors and the target must use the same exact Numerica field, and
/// the supplied product metadata must belong to the target quiver. All
/// source-translated term, inclusion, and arithmetic failures are propagated.
pub fn simple_tensor<F: Field>(
    factors: [&QuiverAlgebraElement<F>; 2],
    target: &QuiverAlgebra<F>,
    product: &QuiverProduct,
) -> Result<QuiverAlgebraElement<F>, PathAlgebraTensorError> {
    if factors[0].is_zero() || factors[1].is_zero() {
        return Ok(target.zero());
    }
    if factors
        .iter()
        .any(|factor| factor.field() != target.field())
    {
        return Err(PathAlgebraTensorError::DifferentCoefficientFields);
    }
    if !target.quiver().is_identical(product.quiver()) {
        return Err(PathAlgebraTensorError::TargetProductQuiverMismatch);
    }

    let parts = [
        path_algebra_element_terms(factors[0])?,
        path_algebra_element_terms(factors[1])?,
    ];
    let mut summands = Vec::with_capacity(parts[0].len().saturating_mul(parts[1].len()));
    for left in &parts[0] {
        for right in &parts[1] {
            let included_path =
                include_in_product_quiver([left.monomial(), right.monomial()], product)?;
            let included = include_in_path_algebra(&included_path, target)?;
            let coefficient = target.field().mul(left.coefficient(), right.coefficient());
            summands.push(included.scale_left(&coefficient)?);
        }
    }
    sum_algebra_elements(&summands).map_err(Into::into)
}

fn product_of_algebra_elements<F: Field>(
    elements: &[QuiverAlgebraElement<F>],
) -> Result<QuiverAlgebraElement<F>, QuiverAlgebraElementError> {
    let Some((first, rest)) = elements.split_first() else {
        return Err(QuiverAlgebraElementError::EmptyProduct);
    };
    let mut product = first.clone();
    for element in rest {
        product = product.multiply(element)?;
    }
    Ok(product)
}

fn sum_algebra_elements<F: Field>(
    elements: &[QuiverAlgebraElement<F>],
) -> Result<QuiverAlgebraElement<F>, QuiverAlgebraElementError> {
    let Some((first, rest)) = elements.split_first() else {
        return Err(QuiverAlgebraElementError::EmptySum);
    };
    let mut sum = first.clone();
    for element in rest {
        sum = sum.add(element)?;
    }
    Ok(sum)
}

fn product_of_paths(paths: &[Path]) -> Result<Path, PathError> {
    let Some((first, rest)) = paths.split_first() else {
        return Err(PathError::EmptyProduct);
    };
    let mut product = first.clone();
    for path in rest {
        product = product.multiply(path)?;
    }
    Ok(product)
}

fn product_generator_path(product: &QuiverProduct, name: &str) -> Result<Path, QuiverProductError> {
    product
        .quiver
        .generator_by_name(name)
        .map(|generator| generator.path())
        .ok_or_else(|| QuiverProductError::MissingProductGenerator { name: name.into() })
}

fn endpoint_vertex(endpoint: PathEndpoint) -> Result<Vertex, PathError> {
    match endpoint {
        PathEndpoint::Vertex(vertex) => Ok(vertex),
        PathEndpoint::Zero(_) => Err(PathError::ZeroHasNoWalk),
    }
}

fn generator_source(generator: &Generator) -> Vertex {
    match generator {
        Generator::Vertex(vertex) => vertex.clone(),
        Generator::Arrow(arrow) => arrow.source(),
    }
}

fn generator_target(generator: &Generator) -> Vertex {
    match generator {
        Generator::Vertex(vertex) => vertex.clone(),
        Generator::Arrow(arrow) => arrow.target(),
    }
}

fn generator_name(generator: &Generator) -> &str {
    match generator {
        Generator::Vertex(vertex) => vertex.name(),
        Generator::Arrow(arrow) => arrow.name(),
    }
}

fn compose_name(pair: &[Generator; 2]) -> Box<str> {
    format!("{}_{}", generator_name(&pair[0]), generator_name(&pair[1])).into_boxed_str()
}

/// Failure in QPA's product-quiver utility layer.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum QuiverProductError {
    /// QPA's decomposition has exactly two one-based factor positions.
    #[error("product-quiver factor index {index} is outside 1..=2")]
    FactorIndexOutOfRange {
        /// Rejected one-based index.
        index: usize,
    },
    /// Inclusion composed a name absent from the product quiver record.
    #[error("cannot find product-quiver generator named {name}")]
    MissingProductGenerator {
        /// Missing composed name.
        name: Box<str>,
    },
    /// Projection found no decomposition slot for a product generator name.
    #[error("cannot decompose product-quiver generator named {name}")]
    MissingDecompositionEntry {
        /// Missing product-generator name.
        name: Box<str>,
    },
    /// Projection was paired with a different product-quiver object.
    #[error("the path does not belong to this product quiver")]
    PathNotInProductQuiver,
    /// Product quiver construction failed.
    #[error(transparent)]
    Quiver(#[from] QuiverError),
    /// A source path operation failed.
    #[error(transparent)]
    Path(#[from] PathError),
}

/// Failure in QPA's path-algebra portion of the tensor helper layer.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PathAlgebraTensorError {
    /// QPA's algebra record had no component with the required generator name.
    #[error("cannot find path-algebra generator named {name}")]
    MissingAlgebraGenerator {
        /// Missing source generator name.
        name: Box<str>,
    },
    /// Factor coefficients cannot act on the requested target algebra.
    #[error("simple-tensor factors and target use different coefficient fields")]
    DifferentCoefficientFields,
    /// The explicit Rust decomposition metadata belongs to another quiver.
    #[error("the target algebra does not use the supplied product quiver")]
    TargetProductQuiverMismatch,
    /// A source path operation failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// Quiver-algebra element arithmetic or reduction failed.
    #[error(transparent)]
    Element(#[from] QuiverAlgebraElementError),
    /// Product-quiver inclusion failed.
    #[error(transparent)]
    Product(#[from] QuiverProductError),
}

#[cfg(test)]
mod tests {
    use crate::{
        ArrowSpec, NamedArrowSpec, PathAlgebra, PathError, Quiver, QuiverAlgebra,
        QuiverAlgebraElement, QuotientPathAlgebra, VertexIndex,
    };
    use numerica::domains::{
        Field, Ring,
        finite_field::Z2,
        integer::Integer,
        rational::{Q, Rational, RationalField},
    };
    use serde_json::{Value, json};

    use super::{
        PathAlgebraElementTermRecord, QuiverProductError, include_in_path_algebra,
        include_in_product_quiver, path_algebra_element_terms, project_from_product_quiver,
        quiver_product, reverse_path, simple_tensor, vertices_of_path_algebra,
        walk_of_path_or_vertex,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn factors() -> [Quiver; 2] {
        [
            Quiver::from_named(
                vec!["u1".into(), "u2".into()],
                vec![NamedArrowSpec::named("u1", "u2", "a")],
            )
            .unwrap(),
            Quiver::from_named(
                vec!["v1".into(), "v2".into(), "v3".into()],
                vec![
                    NamedArrowSpec::named("v1", "v2", "b"),
                    NamedArrowSpec::named("v2", "v3", "c"),
                ],
            )
            .unwrap(),
        ]
    }

    fn q(numerator: i64) -> Rational {
        Q.to_element(Integer::from(numerator), Integer::from(1), true)
    }

    fn terms_json<F: Field>(
        terms: &[PathAlgebraElementTermRecord<F>],
        coefficient_text: impl Fn(&F::Element) -> String,
    ) -> Value {
        Value::Array(
            terms
                .iter()
                .map(|term| {
                    json!({
                        "monomial": term.monomial().to_string(),
                        "coefficient": coefficient_text(term.coefficient()),
                    })
                })
                .collect(),
        )
    }

    fn rational_element_json(element: &QuiverAlgebraElement<RationalField>) -> Value {
        terms_json(
            &path_algebra_element_terms(element).unwrap(),
            ToString::to_string,
        )
    }

    fn binary_element_json(element: &QuiverAlgebraElement<Z2>) -> Value {
        terms_json(&path_algebra_element_terms(element).unwrap(), |_| {
            "Z(2)^0".to_owned()
        })
    }

    #[test]
    fn product_construction_and_decomposition_follow_cartesian_order() {
        let [q1, q2] = factors();
        let product = quiver_product(&q1, &q2).unwrap();
        assert_eq!(product.quiver().number_of_vertices(), 6);
        assert_eq!(product.quiver().number_of_arrows(), 7);
        assert_eq!(
            product
                .quiver()
                .vertices()
                .map(|vertex| vertex.name().to_owned())
                .collect::<Vec<_>>(),
            ["u1_v1", "u1_v2", "u1_v3", "u2_v1", "u2_v2", "u2_v3"]
        );
        assert_eq!(
            product
                .quiver()
                .arrows()
                .map(|arrow| arrow.name().to_owned())
                .collect::<Vec<_>>(),
            ["u1_b", "u1_c", "u2_b", "u2_c", "a_v1", "a_v2", "a_v3"]
        );
        assert!(product.decomposition().is_bound(1));
        assert!(product.decomposition().is_bound(2));
        assert!(!product.decomposition().is_bound(3));
        assert!(product.decomposition().factor(1).unwrap().is_identical(&q1));
        assert!(product.decomposition().factor(2).unwrap().is_identical(&q2));
        assert_eq!(
            product.decomposition().factor(3),
            Err(QuiverProductError::FactorIndexOutOfRange { index: 3 })
        );
        assert_eq!(
            product.decomposition().to_string(),
            "<quiver product decomposition>"
        );
    }

    #[test]
    fn inclusion_and_projection_follow_qpa_name_lookup_and_products() {
        let [q1, q2] = factors();
        let product = quiver_product(&q1, &q2).unwrap();
        let a = q1.generator_by_name("a").unwrap().path();
        let b = q2.generator_by_name("b").unwrap().path();
        let c = q2.generator_by_name("c").unwrap().path();
        let bc = b.multiply(&c).unwrap();

        let included = include_in_product_quiver([&a, &bc], &product).unwrap();
        assert_eq!(included.to_string(), "a_v1*u2_b*u2_c");
        assert_eq!(
            project_from_product_quiver(1, &included, &product)
                .unwrap()
                .to_string(),
            "a"
        );
        assert_eq!(
            project_from_product_quiver(2, &included, &product)
                .unwrap()
                .to_string(),
            "b*c"
        );

        let u1 = q1.generator_by_name("u1").unwrap().path();
        let v2 = q2.generator_by_name("v2").unwrap().path();
        assert_eq!(
            include_in_product_quiver([&u1, &v2], &product)
                .unwrap()
                .to_string(),
            "u1_v2"
        );
    }

    #[test]
    fn walk_and_reverse_preserve_source_boundaries() {
        let quiver = Quiver::from_indexed(
            3,
            vec![
                ArrowSpec::named(vertex(1), vertex(2), "a"),
                ArrowSpec::named(vertex(2), vertex(3), "b"),
            ],
        )
        .unwrap();
        let a = quiver.generator_by_name("a").unwrap().path();
        let b = quiver.generator_by_name("b").unwrap().path();
        let ab = a.multiply(&b).unwrap();
        assert_eq!(walk_of_path_or_vertex(&ab).unwrap().len(), 2);
        assert_eq!(reverse_path(&a).unwrap(), a);
        assert!(reverse_path(&ab).unwrap().is_zero());
        assert_eq!(
            reverse_path(&quiver.vertices().next().unwrap().path()),
            Err(PathError::EmptyProduct)
        );
        assert_eq!(
            walk_of_path_or_vertex(&quiver.zero()),
            Err(PathError::ZeroHasNoWalk)
        );
    }

    #[test]
    fn pinned_qpa_product_quiver_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalgtensor_quiver.schema-v1.json"
        ))
        .unwrap();
        let [q1, q2] = factors();
        let product = quiver_product(&q1, &q2).unwrap();
        let a = q1.generator_by_name("a").unwrap().path();
        let b = q2.generator_by_name("b").unwrap().path();
        let c = q2.generator_by_name("c").unwrap().path();
        let bc = b.multiply(&c).unwrap();
        let included = include_in_product_quiver([&a, &bc], &product).unwrap();
        let u1 = q1.generator_by_name("u1").unwrap().path();
        let v2 = q2.generator_by_name("v2").unwrap().path();

        let endpoints = product
            .quiver()
            .arrows()
            .map(|arrow| json!([arrow.source().name(), arrow.target().name()]))
            .collect::<Vec<_>>();
        let actual = json!({
            "product": {
                "vertices": product.quiver().vertices().map(|vertex| vertex.name().to_owned()).collect::<Vec<_>>(),
                "arrows": product.quiver().arrows().map(|arrow| arrow.name().to_owned()).collect::<Vec<_>>(),
                "arrow_endpoints": endpoints,
            },
            "decomposition": {
                "length": 2,
                "bound": [
                    product.decomposition().is_bound(1),
                    product.decomposition().is_bound(2),
                    product.decomposition().is_bound(3),
                ],
                "factor_identity": [
                    product.decomposition().factor(1).unwrap().is_identical(&q1),
                    product.decomposition().factor(2).unwrap().is_identical(&q2),
                ],
            },
            "walks": {
                "vertex": walk_of_path_or_vertex(&u1).unwrap().iter().map(super::generator_name).collect::<Vec<_>>(),
                "arrow": walk_of_path_or_vertex(&a).unwrap().iter().map(super::generator_name).collect::<Vec<_>>(),
                "path": walk_of_path_or_vertex(&bc).unwrap().iter().map(super::generator_name).collect::<Vec<_>>(),
            },
            "include": {
                "path": included.to_string(),
                "vertices": include_in_product_quiver([&u1, &v2], &product).unwrap().to_string(),
            },
            "project": {
                "first": project_from_product_quiver(1, &included, &product).unwrap().to_string(),
                "second": project_from_product_quiver(2, &included, &product).unwrap().to_string(),
            }
        });
        for key in ["product", "decomposition", "walks", "include", "project"] {
            assert_eq!(actual[key], expected[key], "fixture section {key}");
        }

        let loops = Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "s"),
                ArrowSpec::named(vertex(1), vertex(1), "t"),
            ],
        )
        .unwrap();
        let st = loops
            .generator_by_name("s")
            .unwrap()
            .path()
            .multiply(&loops.generator_by_name("t").unwrap().path())
            .unwrap();
        assert_eq!(
            reverse_path(&st).unwrap().to_string(),
            expected["reverse"]["loop_path"]
        );
    }

    #[test]
    fn pinned_qpa_path_algebra_element_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalgtensor_elements.schema-v1.json"
        ))
        .unwrap();

        let quiver = Quiver::from_named(
            vec!["u".into(), "v".into()],
            vec![
                NamedArrowSpec::named("u", "u", "a"),
                NamedArrowSpec::named("u", "v", "b"),
            ],
        )
        .unwrap();
        let parent = PathAlgebra::new(Q, quiver.clone()).unwrap();
        let u = quiver.generator_by_name("u").unwrap().path();
        let a = quiver.generator_by_name("a").unwrap().path();
        let b = quiver.generator_by_name("b").unwrap().path();
        let aa = a.multiply(&a).unwrap();
        let ab = a.multiply(&b).unwrap();
        let aa_relation = parent.element(aa.clone()).unwrap();
        let quotient = QuotientPathAlgebra::by_relations_high_level(
            parent.clone(),
            std::slice::from_ref(&aa_relation),
        )
        .unwrap();
        let parent_algebra = QuiverAlgebra::Path(parent.clone());
        let quotient_algebra = QuiverAlgebra::Quotient(quotient.clone());
        let parent_element = QuiverAlgebraElement::Path(
            parent
                .from_terms([(u, q(2)), (a.clone(), q(3)), (ab, q(5))])
                .unwrap(),
        );
        let quotient_element = QuiverAlgebraElement::Quotient(
            quotient
                .element(match &parent_element {
                    QuiverAlgebraElement::Path(element) => element.clone(),
                    QuiverAlgebraElement::Quotient(_) => unreachable!(),
                })
                .unwrap(),
        );
        let parent_vertices = vertices_of_path_algebra(&parent_algebra).unwrap();
        let quotient_vertices = vertices_of_path_algebra(&quotient_algebra).unwrap();
        let rational = json!({
            "include_parent": rational_element_json(
                &include_in_path_algebra(&a.multiply(&b).unwrap(), &parent_algebra).unwrap()
            ),
            "include_quotient_zero": rational_element_json(
                &include_in_path_algebra(&aa, &quotient_algebra).unwrap()
            ),
            "parent_vertices": parent_vertices.iter().map(rational_element_json).collect::<Vec<_>>(),
            "quotient_vertices": quotient_vertices.iter().map(rational_element_json).collect::<Vec<_>>(),
            "parent_terms": rational_element_json(&parent_element),
            "quotient_terms": rational_element_json(&quotient_element),
        });
        assert_eq!(rational, expected["rational"]);

        let binary_quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "c")]).unwrap();
        let binary_parent = PathAlgebra::new(Z2, binary_quiver.clone()).unwrap();
        let v1 = binary_quiver.generator_by_name("v1").unwrap().path();
        let c = binary_quiver.generator_by_name("c").unwrap().path();
        let binary_element = QuiverAlgebraElement::Path(
            binary_parent
                .from_terms([(v1, Z2.one()), (c.clone(), Z2.one())])
                .unwrap(),
        );
        let binary_algebra = QuiverAlgebra::Path(binary_parent);
        let binary_vertices = vertices_of_path_algebra(&binary_algebra).unwrap();
        let binary = json!({
            "include": binary_element_json(&include_in_path_algebra(&c, &binary_algebra).unwrap()),
            "vertices": binary_vertices.iter().map(binary_element_json).collect::<Vec<_>>(),
            "terms": binary_element_json(&binary_element),
        });
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn pinned_qpa_simple_tensor_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalgtensor_simple.schema-v1.json"
        ))
        .unwrap();
        let q1 = Quiver::from_named(
            vec!["u1".into(), "u2".into()],
            vec![NamedArrowSpec::named("u1", "u2", "a")],
        )
        .unwrap();
        let q2 = Quiver::from_named(vec!["v".into()], vec![NamedArrowSpec::named("v", "v", "b")])
            .unwrap();
        let product = quiver_product(&q1, &q2).unwrap();

        let left_parent = PathAlgebra::new(Q, q1.clone()).unwrap();
        let right_parent = PathAlgebra::new(Q, q2.clone()).unwrap();
        let b_path = q2.generator_by_name("b").unwrap().path();
        let bb_relation = right_parent
            .element(b_path.multiply(&b_path).unwrap())
            .unwrap();
        let right_quotient = QuotientPathAlgebra::by_relations_high_level(
            right_parent.clone(),
            std::slice::from_ref(&bb_relation),
        )
        .unwrap();
        let target_parent = PathAlgebra::new(Q, product.quiver().clone()).unwrap();
        let a_v = product.quiver().generator_by_name("a_v").unwrap().path();
        let u2_b = product.quiver().generator_by_name("u2_b").unwrap().path();
        let killed = target_parent.element(a_v.multiply(&u2_b).unwrap()).unwrap();
        let target_quotient = QuotientPathAlgebra::by_relations_high_level(
            target_parent.clone(),
            std::slice::from_ref(&killed),
        )
        .unwrap();

        let u1 = q1.generator_by_name("u1").unwrap().path();
        let a_path = q1.generator_by_name("a").unwrap().path();
        let vertex_v = q2.generator_by_name("v").unwrap().path();
        let left_element = QuiverAlgebraElement::Path(
            left_parent
                .from_terms([(u1, q(2)), (a_path, q(3))])
                .unwrap(),
        );
        let right_representative = right_parent
            .from_terms([(vertex_v, q(5)), (b_path, q(7))])
            .unwrap();
        let right_element =
            QuiverAlgebraElement::Quotient(right_quotient.element(right_representative).unwrap());
        let left_algebra = QuiverAlgebra::Path(left_parent.clone());
        let target_algebra = QuiverAlgebra::Path(target_parent);
        let quotient_target_algebra = QuiverAlgebra::Quotient(target_quotient);
        let zero = left_algebra.zero();
        let a_element = left_algebra
            .generator_element_by_name("a")
            .unwrap()
            .unwrap();
        let b_element = QuiverAlgebra::Quotient(right_quotient)
            .generator_element_by_name("b")
            .unwrap()
            .unwrap();
        let rational = json!({
            "quotient_factor": rational_element_json(
                &simple_tensor([&left_element, &right_element], &target_algebra, &product).unwrap()
            ),
            "zero_short_circuit": rational_element_json(
                &simple_tensor([&zero, &right_element], &target_algebra, &product).unwrap()
            ),
            "quotient_target_reduction": rational_element_json(
                &simple_tensor([&a_element, &b_element], &quotient_target_algebra, &product).unwrap()
            ),
        });
        assert_eq!(rational, expected["rational"]);

        let binary_left_parent = PathAlgebra::new(Z2, q1).unwrap();
        let binary_right_parent = PathAlgebra::new(Z2, q2).unwrap();
        let binary_target_parent = PathAlgebra::new(Z2, product.quiver().clone()).unwrap();
        let binary_left = QuiverAlgebra::Path(binary_left_parent);
        let binary_right = QuiverAlgebra::Path(binary_right_parent);
        let binary_target = QuiverAlgebra::Path(binary_target_parent);
        let binary_x = binary_left
            .generator_element_by_name("u1")
            .unwrap()
            .unwrap()
            .add(&binary_left.generator_element_by_name("a").unwrap().unwrap())
            .unwrap();
        let binary_y = binary_right
            .generator_element_by_name("v")
            .unwrap()
            .unwrap()
            .add(
                &binary_right
                    .generator_element_by_name("b")
                    .unwrap()
                    .unwrap(),
            )
            .unwrap();
        let binary = json!({
            "expansion": binary_element_json(
                &simple_tensor([&binary_x, &binary_y], &binary_target, &product).unwrap()
            ),
        });
        assert_eq!(binary, expected["binary"]);
    }
}
