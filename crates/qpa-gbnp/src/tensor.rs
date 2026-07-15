//! QPA tensor products of path algebras and their quotients.
//!
//! This translates `TensorProductOfAlgebras` and
//! `TensorProductOfPathAlgebras` from `lib/pathalgtensor.gi`.  The complete
//! tensor Gröbner list is the source theorem's induced factor bases followed
//! by commutativity relations; it is attached directly rather than recomputed
//! by an independent algorithm.

use numerica::domains::Field;
use qpa_core::{
    PathAlgebra, PathAlgebraElement, PathAlgebraError, PathAlgebraIdeal, PathAlgebraTensorError,
    PathError, QuiverAlgebra, QuiverAlgebraElement, QuiverAlgebraElementError,
    QuiverAlgebraHomomorphism, QuiverAlgebraHomomorphismError, QuiverProduct, QuiverProductError,
    QuotientError, QuotientPathAlgebra, include_in_path_algebra, include_in_product_quiver,
    quiver_product, simple_tensor, vertices_of_path_algebra,
};
use thiserror::Error;

use crate::{PathalgGbnpError, ensure_quotient_groebner_basis};

/// QPA's error text reached when the tensor presentation has no generators
/// for its defining ideal.
pub const EMPTY_TENSOR_RELATIONS_ERROR: &str =
    "first argument must be a ring,\nsecond argument must be a list of generators";

/// A tensor-product quotient together with QPA's retained product quiver and
/// original algebra factors.
#[derive(Clone, Debug)]
pub struct TensorProductAlgebra<F: Field> {
    algebra: QuiverAlgebra<F>,
    product: QuiverProduct,
    factors: [QuiverAlgebra<F>; 2],
}

impl<F: Field> TensorProductAlgebra<F> {
    /// Return the tensor product as QPA's quiver-algebra union.
    #[must_use]
    pub const fn algebra(&self) -> &QuiverAlgebra<F> {
        &self.algebra
    }

    /// Return the quotient constructed by the source method.
    #[must_use]
    pub const fn quotient(&self) -> &QuotientPathAlgebra<F> {
        match &self.algebra {
            QuiverAlgebra::Quotient(quotient) => quotient,
            QuiverAlgebra::Path(_) => unreachable!(),
        }
    }

    /// Return QPA's product quiver with its decomposition metadata.
    #[must_use]
    pub const fn product_quiver(&self) -> &QuiverProduct {
        &self.product
    }

    /// Return QPA's retained `TensorProductDecomposition` list.
    #[must_use]
    pub const fn factors(&self) -> &[QuiverAlgebra<F>; 2] {
        &self.factors
    }

    /// Extract one retained factor using QPA's one-based indexing.
    ///
    /// # Errors
    ///
    /// Returns an explicit range error where GAP list access would fail.
    pub fn factor(&self, index: usize) -> Result<&QuiverAlgebra<F>, TensorProductError> {
        self.factors
            .get(index.wrapping_sub(1))
            .ok_or(TensorProductError::FactorIndexOutOfRange { index })
    }

    /// Apply QPA's `SimpleTensor` using this tensor product's retained product
    /// metadata.
    ///
    /// # Errors
    ///
    /// Propagates source-translated simple-tensor failures.
    pub fn simple_tensor(
        &self,
        factors: [&QuiverAlgebraElement<F>; 2],
    ) -> Result<QuiverAlgebraElement<F>, TensorProductError> {
        Ok(simple_tensor(
            factors,
            self.algebra(),
            self.product_quiver(),
        )?)
    }

    /// Translate `TensorAlgebraInclusion(T, n)`.
    ///
    /// QPA selects the first factor only for the integer `1`; every other
    /// integer selects the second factor, including zero and negative values.
    ///
    /// # Errors
    ///
    /// Propagates generator construction, simple-tensor, and homomorphism
    /// construction failures.
    pub fn inclusion(
        &self,
        factor_number: i64,
    ) -> Result<QuiverAlgebraHomomorphism<F>, TensorProductError> {
        tensor_algebra_inclusion(self, factor_number)
    }
}

/// Translate QPA's binary `TensorProductOfAlgebras` method.
///
/// This is the method-level delegator to [`tensor_product_of_path_algebras`].
///
/// # Errors
///
/// Propagates the global constructor's source-translated failures.
pub fn tensor_product_of_algebras<F: Field>(
    left: &QuiverAlgebra<F>,
    right: &QuiverAlgebra<F>,
) -> Result<TensorProductAlgebra<F>, TensorProductError> {
    tensor_product_of_path_algebras([left, right])
}

/// Translate QPA's `TensorAlgebraInclusion(T, n)` method.
///
/// The mapping sends `x` to `x tensor 1` when `n = 1`, and to `1 tensor x`
/// for every other integer. Its explicitly retained generators and images
/// are QPA's `GeneratorsOfAlgebra` list, including the leading identity that
/// GAP adds for quotient path algebras.
///
/// # Errors
///
/// Propagates generator construction, simple-tensor, and homomorphism
/// construction failures.
pub fn tensor_algebra_inclusion<F: Field>(
    tensor: &TensorProductAlgebra<F>,
    factor_number: i64,
) -> Result<QuiverAlgebraHomomorphism<F>, TensorProductError> {
    let selected = usize::from(factor_number != 1);
    let source = tensor.factors[selected].clone();
    let generators = source.generators()?;
    let other_one = tensor.factors[1 - selected].one()?;
    let images = generators
        .iter()
        .map(|generator| {
            if factor_number == 1 {
                tensor.simple_tensor([generator, &other_one])
            } else {
                tensor.simple_tensor([&other_one, generator])
            }
        })
        .collect::<Result<Vec<_>, TensorProductError>>()?;
    Ok(QuiverAlgebraHomomorphism::by_generator_images(
        source,
        tensor.algebra.clone(),
        images,
    )?)
}

/// Translate QPA's `TensorProductOfPathAlgebras([A, B])` global function.
///
/// Relation order is exactly `comm_rels` followed by the relations induced
/// from the first and second factors. The attached complete list is exactly
/// the induced factor Gröbner bases followed by `comm_rels`, matching the
/// Bjørshol theorem path in QPA 1.37.
///
/// # Errors
///
/// Rejects unequal fields before constructing a product. It preserves QPA's
/// empty-ideal generator failure and propagates product, GBNP, and quotient
/// construction errors.
pub fn tensor_product_of_path_algebras<F: Field>(
    factors: [&QuiverAlgebra<F>; 2],
) -> Result<TensorProductAlgebra<F>, TensorProductError> {
    if factors[0].field() != factors[1].field() {
        return Err(TensorProductError::DifferentCoefficientFields);
    }

    let product = quiver_product(factors[0].quiver(), factors[1].quiver())?;
    let product_parent = PathAlgebra::new(factors[0].field().clone(), product.quiver().clone())?;
    let product_path_algebra = QuiverAlgebra::Path(product_parent.clone());

    let commutativity =
        commutativity_relations(factors[0], factors[1], &product_path_algebra, &product)?;
    let factor_relations = [
        relation_elements(factors[0].relations()),
        relation_elements(factors[1].relations()),
    ];
    let factor_vertices = [
        vertices_of_path_algebra(factors[0])?,
        vertices_of_path_algebra(factors[1])?,
    ];
    let induced_relations = induced_elements(
        &factor_relations,
        &factor_vertices,
        &product_path_algebra,
        &product,
    )?;

    let mut tensor_relations = commutativity.clone();
    tensor_relations.extend(induced_relations);
    if tensor_relations.is_empty() {
        return Err(TensorProductError::EmptyTensorRelations);
    }

    let factor_bases = [
        groebner_elements(factors[0])?,
        groebner_elements(factors[1])?,
    ];
    let induced_basis = induced_elements(
        &factor_bases,
        &factor_vertices,
        &product_path_algebra,
        &product,
    )?;
    let mut tensor_basis = induced_basis;
    tensor_basis.extend(commutativity);
    if tensor_basis.is_empty() {
        return Err(TensorProductError::EmptyTensorGroebnerBasis);
    }

    let ideal = PathAlgebraIdeal::new(product_parent, tensor_relations)?;
    ideal.with_complete_groebner_basis(&tensor_basis)?;
    let quotient = QuotientPathAlgebra::by_ideal(&ideal)?;
    Ok(TensorProductAlgebra {
        algebra: QuiverAlgebra::Quotient(quotient),
        product,
        factors: [factors[0].clone(), factors[1].clone()],
    })
}

fn commutativity_relations<F: Field>(
    left: &QuiverAlgebra<F>,
    right: &QuiverAlgebra<F>,
    target: &QuiverAlgebra<F>,
    product: &QuiverProduct,
) -> Result<Vec<PathAlgebraElement<F>>, TensorProductError> {
    let mut relations = Vec::new();
    for left_arrow in left.quiver().arrows() {
        for right_arrow in right.quiver().arrows() {
            let left_arrow_path = left_arrow.path();
            let right_arrow_path = right_arrow.path();
            let right_source = right_arrow.source().path();
            let right_target = right_arrow.target().path();
            let left_source = left_arrow.source().path();
            let left_target = left_arrow.target().path();

            let first =
                include_in_product_quiver([&left_arrow_path, &right_source], product)?.multiply(
                    &include_in_product_quiver([&left_target, &right_arrow_path], product)?,
                )?;
            let second =
                include_in_product_quiver([&left_source, &right_arrow_path], product)?.multiply(
                    &include_in_product_quiver([&left_arrow_path, &right_target], product)?,
                )?;
            let first = include_in_path_algebra(&first, target)?;
            let second = include_in_path_algebra(&second, target)?;
            relations.push(path_element(first.sub(&second)?)?);
        }
    }
    Ok(relations)
}

fn relation_elements<F: Field>(
    relations: &[PathAlgebraElement<F>],
) -> Vec<QuiverAlgebraElement<F>> {
    relations
        .iter()
        .cloned()
        .map(QuiverAlgebraElement::Path)
        .collect()
}

fn groebner_elements<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Vec<QuiverAlgebraElement<F>>, TensorProductError> {
    match algebra {
        QuiverAlgebra::Path(_) => Ok(Vec::new()),
        QuiverAlgebra::Quotient(quotient) => Ok(ensure_quotient_groebner_basis(quotient)?
            .relations()
            .iter()
            .cloned()
            .map(QuiverAlgebraElement::Path)
            .collect()),
    }
}

fn induced_elements<F: Field>(
    factor_elements: &[Vec<QuiverAlgebraElement<F>>; 2],
    factor_vertices: &[Vec<QuiverAlgebraElement<F>>; 2],
    target: &QuiverAlgebra<F>,
    product: &QuiverProduct,
) -> Result<Vec<PathAlgebraElement<F>>, TensorProductError> {
    let mut induced = Vec::new();
    for relation in &factor_elements[0] {
        for vertex in &factor_vertices[1] {
            induced.push(path_element(simple_tensor(
                [relation, vertex],
                target,
                product,
            )?)?);
        }
    }
    for vertex in &factor_vertices[0] {
        for relation in &factor_elements[1] {
            induced.push(path_element(simple_tensor(
                [vertex, relation],
                target,
                product,
            )?)?);
        }
    }
    Ok(induced)
}

fn path_element<F: Field>(
    element: QuiverAlgebraElement<F>,
) -> Result<PathAlgebraElement<F>, TensorProductError> {
    match element {
        QuiverAlgebraElement::Path(element) => Ok(element),
        QuiverAlgebraElement::Quotient(_) => Err(TensorProductError::UnexpectedTargetVariant),
    }
}

/// Failure in QPA's full tensor-product constructor.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum TensorProductError {
    /// QPA checks equality of the two left-acting domains first.
    #[error("path algebras must be over the same field")]
    DifferentCoefficientFields,
    /// QPA's `Ideal(product_pa, [])` dispatch rejects an empty generator list.
    #[error("{EMPTY_TENSOR_RELATIONS_ERROR}")]
    EmptyTensorRelations,
    /// A nonempty defining list unexpectedly had no supplied complete basis.
    #[error("QPA cannot attach an empty tensor-product Groebner basis")]
    EmptyTensorGroebnerBasis,
    /// The internal product target must remain the path algebra constructed by
    /// the source method.
    #[error("a tensor presentation helper returned a quotient target element")]
    UnexpectedTargetVariant,
    /// QPA's retained decomposition contains exactly two factors.
    #[error("tensor-product factor index {index} is outside 1..=2")]
    FactorIndexOutOfRange {
        /// Rejected one-based index.
        index: usize,
    },
    /// Product-quiver construction or inclusion failed.
    #[error(transparent)]
    Product(#[from] QuiverProductError),
    /// Product-path multiplication failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// Path-algebra construction failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Core path-algebra tensor helpers failed.
    #[error(transparent)]
    Tensor(#[from] PathAlgebraTensorError),
    /// Arithmetic over QPA's quiver-algebra element union failed.
    #[error(transparent)]
    Element(#[from] QuiverAlgebraElementError),
    /// The configured QPA/GBNP basis path failed.
    #[error(transparent)]
    Backend(#[from] PathalgGbnpError),
    /// Ideal attachment or quotient construction failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Construction or application of the stored algebra mapping failed.
    #[error(transparent)]
    Homomorphism(#[from] QuiverAlgebraHomomorphismError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field,
        finite_field::Z2,
        integer::Integer,
        rational::{Q, Rational, RationalField},
    };
    use qpa_core::{
        NamedArrowSpec, Path, PathAlgebra, PathAlgebraElement, PathEndpoint, Quiver, QuiverAlgebra,
        QuiverAlgebraElement, path_algebra_element_terms,
    };
    use serde_json::{Value, json};

    use super::{TensorProductError, tensor_algebra_inclusion, tensor_product_of_algebras};
    use crate::quotient_by_relations;

    fn q(numerator: i64) -> Rational {
        Q.to_element(Integer::from(numerator), Integer::from(1), true)
    }

    fn path_text(path: &Path) -> String {
        if path.is_zero() {
            return "0".to_owned();
        }
        if let PathEndpoint::Vertex(vertex) = path.source()
            && path.is_vertex()
        {
            return vertex.name().to_owned();
        }
        path.walk()
            .unwrap()
            .iter()
            .map(|&index| path.quiver().arrow(index).unwrap().name().to_owned())
            .collect::<Vec<_>>()
            .join("*")
    }

    fn terms_json<F: Field>(
        element: &PathAlgebraElement<F>,
        coefficient_text: impl Fn(&F::Element) -> String,
    ) -> Value {
        let wrapped = QuiverAlgebraElement::Path(element.clone());
        Value::Array(
            path_algebra_element_terms(&wrapped)
                .unwrap()
                .iter()
                .map(|term| {
                    json!({
                        "monomial": path_text(term.monomial()),
                        "coefficient": coefficient_text(term.coefficient()),
                    })
                })
                .collect(),
        )
    }

    fn rational_terms(element: &PathAlgebraElement<RationalField>) -> Value {
        terms_json(element, ToString::to_string)
    }

    fn binary_terms(element: &PathAlgebraElement<Z2>) -> Value {
        terms_json(element, |_| "Z(2)^0".to_owned())
    }

    fn rational_union_terms(element: &QuiverAlgebraElement<RationalField>) -> Value {
        union_terms(element, ToString::to_string)
    }

    fn binary_union_terms(element: &QuiverAlgebraElement<Z2>) -> Value {
        union_terms(element, |_| "Z(2)^0".to_owned())
    }

    fn union_terms<F: Field>(
        element: &QuiverAlgebraElement<F>,
        coefficient_text: impl Fn(&F::Element) -> String,
    ) -> Value {
        Value::Array(
            path_algebra_element_terms(element)
                .unwrap()
                .iter()
                .map(|term| {
                    json!({
                        "monomial": path_text(term.monomial()),
                        "coefficient": coefficient_text(term.coefficient()),
                    })
                })
                .collect(),
        )
    }

    fn quivers() -> [Quiver; 2] {
        [
            Quiver::from_named(
                vec!["u1".into(), "u2".into()],
                vec![NamedArrowSpec::named("u1", "u2", "a")],
            )
            .unwrap(),
            Quiver::from_named(vec!["v".into()], vec![NamedArrowSpec::named("v", "v", "b")])
                .unwrap(),
        ]
    }

    #[test]
    fn pinned_qpa_full_tensor_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalgtensor_full.schema-v1.json"
        ))
        .unwrap();
        let [q1, q2] = quivers();
        let left_parent = PathAlgebra::new(Q, q1.clone()).unwrap();
        let right_parent = PathAlgebra::new(Q, q2.clone()).unwrap();
        let b_path = q2.generator_by_name("b").unwrap().path();
        let right_relation = right_parent
            .element(b_path.multiply(&b_path).unwrap())
            .unwrap();
        let right_quotient =
            quotient_by_relations(right_parent, std::slice::from_ref(&right_relation)).unwrap();
        let left = QuiverAlgebra::Path(left_parent.clone());
        let right = QuiverAlgebra::Quotient(right_quotient);
        let tensor = tensor_product_of_algebras(&left, &right).unwrap();
        let quotient = tensor.quotient();
        let basis = quotient.ideal().groebner_basis().unwrap();
        let left_arrow = left.generator_element_by_name("a").unwrap().unwrap();
        let right_arrow = right.generator_element_by_name("b").unwrap().unwrap();
        let rational = json!({
            "vertices": tensor.product_quiver().quiver().vertices().map(|vertex| vertex.name().to_owned()).collect::<Vec<_>>(),
            "arrows": tensor.product_quiver().quiver().arrows().map(|arrow| arrow.name().to_owned()).collect::<Vec<_>>(),
            "relations": quotient.relations().iter().map(rational_terms).collect::<Vec<_>>(),
            "groebner_basis": basis.relations().iter().map(rational_terms).collect::<Vec<_>>(),
            "factor_identity": [
                matches!(tensor.factor(1).unwrap(), QuiverAlgebra::Path(algebra) if algebra.is_identical(&left_parent)),
                tensor.factor(2).unwrap() == &right,
            ],
            "simple_tensor": rational_union_terms(&tensor.simple_tensor([&left_arrow, &right_arrow]).unwrap()),
        });
        assert_eq!(rational, expected["rational"]);

        let binary_left = QuiverAlgebra::Path(PathAlgebra::new(Z2, q1).unwrap());
        let binary_right = QuiverAlgebra::Path(PathAlgebra::new(Z2, q2).unwrap());
        let binary_tensor = tensor_product_of_algebras(&binary_left, &binary_right).unwrap();
        let binary_basis = binary_tensor.quotient().ideal().groebner_basis().unwrap();
        let binary = json!({
            "relations": binary_tensor.quotient().relations().iter().map(binary_terms).collect::<Vec<_>>(),
            "groebner_basis": binary_basis.relations().iter().map(binary_terms).collect::<Vec<_>>(),
        });
        assert_eq!(binary, expected["binary"]);
    }

    #[test]
    fn empty_source_presentation_failure_is_explicit() {
        let quiver = Quiver::from_named(vec!["v".into()], Vec::new()).unwrap();
        let left = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver.clone()).unwrap());
        let right = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        assert!(matches!(
            tensor_product_of_algebras(&left, &right),
            Err(TensorProductError::EmptyTensorRelations)
        ));
        assert_eq!(
            TensorProductError::EmptyTensorRelations.to_string(),
            super::EMPTY_TENSOR_RELATIONS_ERROR
        );
    }

    #[test]
    fn scalar_coefficients_are_multiplied_by_the_target_field() {
        let [q1, q2] = quivers();
        let left = QuiverAlgebra::Path(PathAlgebra::new(Q, q1).unwrap());
        let right = QuiverAlgebra::Path(PathAlgebra::new(Q, q2).unwrap());
        let tensor = tensor_product_of_algebras(&left, &right).unwrap();
        let left_arrow = left
            .generator_element_by_name("a")
            .unwrap()
            .unwrap()
            .scale_left(&q(2))
            .unwrap();
        let right_arrow = right
            .generator_element_by_name("b")
            .unwrap()
            .unwrap()
            .scale_left(&q(3))
            .unwrap();
        let result = tensor.simple_tensor([&left_arrow, &right_arrow]).unwrap();
        assert_eq!(result.leading_coefficient(), q(6));
    }

    #[test]
    fn tensor_inclusions_store_qpa_generators_images_and_else_branch() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalgtensor_inclusion.schema-v1.json"
        ))
        .unwrap();
        let [q1, q2] = quivers();
        let left_parent = PathAlgebra::new(Q, q1.clone()).unwrap();
        let right_parent = PathAlgebra::new(Q, q2.clone()).unwrap();
        let b = right_parent.generators()[1].clone();
        let right = QuiverAlgebra::Quotient(
            quotient_by_relations(right_parent, &[b.multiply(&b).unwrap()]).unwrap(),
        );
        let left = QuiverAlgebra::Path(left_parent);
        let tensor = tensor_product_of_algebras(&left, &right).unwrap();

        let left_inclusion = tensor_algebra_inclusion(&tensor, 1).unwrap();
        let right_inclusion = tensor.inclusion(2).unwrap();
        let zero_inclusion = tensor.inclusion(0).unwrap();
        assert_eq!(left_inclusion.source(), &left);
        assert_eq!(right_inclusion.source(), &right);
        assert_eq!(left_inclusion.range(), tensor.algebra());
        assert_eq!(left_inclusion.generators(), left.generators().unwrap());
        assert_eq!(right_inclusion.generators(), right.generators().unwrap());
        assert_eq!(zero_inclusion.generators(), right_inclusion.generators());
        assert_eq!(
            zero_inclusion.generator_images(),
            right_inclusion.generator_images()
        );

        let left_vertex = left.generator_element_by_name("u1").unwrap().unwrap();
        let left_arrow = left.generator_element_by_name("a").unwrap().unwrap();
        let left_element = left_vertex
            .add(&left_arrow.scale_left(&q(2)).unwrap())
            .unwrap();
        let right_vertex = right.generator_element_by_name("v").unwrap().unwrap();
        let right_arrow = right.generator_element_by_name("b").unwrap().unwrap();
        let right_element = right_vertex
            .add(&right_arrow.scale_left(&q(3)).unwrap())
            .unwrap();
        let rational = json!({
            "left_generators": left_inclusion.generators().iter().map(rational_union_terms).collect::<Vec<_>>(),
            "left_images": left_inclusion.generator_images().iter().map(rational_union_terms).collect::<Vec<_>>(),
            "right_generators": right_inclusion.generators().iter().map(rational_union_terms).collect::<Vec<_>>(),
            "right_images": right_inclusion.generator_images().iter().map(rational_union_terms).collect::<Vec<_>>(),
            "zero_selects_second": [
                zero_inclusion.generators() == right_inclusion.generators(),
                zero_inclusion.generator_images() == right_inclusion.generator_images(),
            ],
            "source_range_identity": [
                left_inclusion.source() == &left,
                right_inclusion.source() == &right,
                left_inclusion.range() == tensor.algebra(),
                right_inclusion.range() == tensor.algebra(),
            ],
            "left_element_image": rational_union_terms(&left_inclusion.image(&left_element).unwrap()),
            "right_element_image": rational_union_terms(&right_inclusion.image(&right_element).unwrap()),
            "zero_index_image": rational_union_terms(&zero_inclusion.image(&right_element).unwrap()),
        });
        assert_eq!(rational, expected["rational"]);

        let binary_left = QuiverAlgebra::Path(PathAlgebra::new(Z2, q1).unwrap());
        let binary_right_parent = PathAlgebra::new(Z2, q2).unwrap();
        let binary_b = binary_right_parent.generators()[1].clone();
        let binary_right = QuiverAlgebra::Quotient(
            quotient_by_relations(
                binary_right_parent,
                &[binary_b.multiply(&binary_b).unwrap()],
            )
            .unwrap(),
        );
        let binary_tensor = tensor_product_of_algebras(&binary_left, &binary_right).unwrap();
        let binary_left_inclusion = binary_tensor.inclusion(1).unwrap();
        let binary_right_inclusion = binary_tensor.inclusion(2).unwrap();
        let binary = json!({
            "left_images": binary_left_inclusion.generator_images().iter().map(binary_union_terms).collect::<Vec<_>>(),
            "right_images": binary_right_inclusion.generator_images().iter().map(binary_union_terms).collect::<Vec<_>>(),
        });
        assert_eq!(binary, expected["binary"]);
    }
}
