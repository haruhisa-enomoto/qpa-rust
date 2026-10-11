//! QPA's path-algebra/GBNP conversion layer.
//!
//! This translates `lib/gbnp.gi`: `QPA_Path2Cohen`, its one-vertex free
//! variant, `QPA_Cohen2Path`, `QPA_RelationsForPathAlgebra`, and the checked
//! and unchecked `GBNPGroebnerBasis` methods.

use numerica::domains::Field;
use qpa_core::{
    GroebnerError, PathAlgebra, PathAlgebraElement, PathAlgebraError, PathAlgebraIdeal,
    make_uniform_list, qpa_in_arrow_ideal,
};
use thiserror::Error;

use crate::{NpError, NpMonomial, NpPolynomial, SgrobnerError, clean_np, sgrobner};

mod components;

/// Raw two-list Cohen/NP representation returned by QPA conversion.
///
/// Unlike [`NpPolynomial`], this preserves QPA's source term order until GBNP
/// applies `CleanNP` during preprocessing.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CohenPolynomial<F: Field> {
    field: F,
    monomials: Vec<NpMonomial>,
    coefficients: Vec<F::Element>,
}

impl<F: Field> CohenPolynomial<F> {
    /// Construct the raw two-list representation accepted by
    /// `QPA_Cohen2Path`.
    ///
    /// QPA indexes one coefficient for every monomial and ignores trailing
    /// coefficients. Rust reports the corresponding missing-list-entry
    /// failure before conversion instead of allowing an out-of-bounds access.
    ///
    /// # Errors
    ///
    /// Returns [`QpaGbnpError::MissingCohenCoefficient`] when the coefficient
    /// list is shorter than the monomial list.
    pub fn from_raw(
        field: F,
        monomials: Vec<Vec<i32>>,
        coefficients: Vec<F::Element>,
    ) -> Result<Self, QpaGbnpError> {
        if coefficients.len() < monomials.len() {
            return Err(QpaGbnpError::MissingCohenCoefficient {
                monomials: monomials.len(),
                coefficients: coefficients.len(),
            });
        }
        Ok(Self {
            field,
            monomials: monomials.into_iter().map(NpMonomial::from).collect(),
            coefficients,
        })
    }

    /// Coefficient field carried by this Cohen representation.
    #[must_use]
    pub const fn field(&self) -> &F {
        &self.field
    }

    /// Raw monomial list in QPA exterior-representation order.
    #[must_use]
    pub fn monomials(&self) -> &[NpMonomial] {
        &self.monomials
    }

    /// Raw coefficient list parallel to [`Self::monomials`].
    #[must_use]
    pub fn coefficients(&self) -> &[F::Element] {
        &self.coefficients
    }

    /// Apply GBNP `CleanNP` to this raw representation.
    ///
    /// # Errors
    ///
    /// Propagates the NP raw-list validation error.
    pub fn clean(&self) -> Result<NpPolynomial<F>, NpError> {
        clean_np(
            self.field.clone(),
            self.monomials.clone(),
            self.coefficients.clone(),
        )
    }
}

/// Common read-only view of raw QPA Cohen values and cleaned GBNP NP values.
pub trait CohenRepresentation<F: Field> {
    /// Coefficient field carried by the representation.
    fn coefficient_field(&self) -> &F;

    /// Monomial/coefficient pairs in the representation's current order.
    fn term_pairs(&self) -> Vec<(&NpMonomial, &F::Element)>;
}

impl<F: Field> CohenRepresentation<F> for CohenPolynomial<F> {
    fn coefficient_field(&self) -> &F {
        &self.field
    }

    fn term_pairs(&self) -> Vec<(&NpMonomial, &F::Element)> {
        self.monomials.iter().zip(&self.coefficients).collect()
    }
}

impl<F: Field> CohenRepresentation<F> for NpPolynomial<F> {
    fn coefficient_field(&self) -> &F {
        self.field()
    }

    fn term_pairs(&self) -> Vec<(&NpMonomial, &F::Element)> {
        self.terms()
            .iter()
            .map(|term| (term.monomial(), term.coefficient()))
            .collect()
    }
}

/// Translate `QPA_Path2Cohen` while preserving its raw term order.
///
/// # Errors
///
/// Returns an error if an exterior generator index cannot fit GBNP's signed
/// generator representation.
pub fn path_to_cohen<F: Field>(
    elements: &[PathAlgebraElement<F>],
) -> Result<Vec<CohenPolynomial<F>>, QpaGbnpError> {
    elements.iter().map(path_element_to_cohen).collect()
}

/// Translate `QPA_Path2CohenFree` for a one-vertex path algebra.
///
/// The sole vertex word `[1]` becomes GBNP's empty identity word.
///
/// # Errors
///
/// Returns an error unless every element belongs to a one-vertex path algebra,
/// or if a generator index does not fit GBNP's representation.
pub fn path_to_cohen_free<F: Field>(
    elements: &[PathAlgebraElement<F>],
) -> Result<Vec<CohenPolynomial<F>>, QpaGbnpError> {
    let mut result = Vec::with_capacity(elements.len());
    for element in elements {
        if element.algebra().quiver().number_of_vertices() != 1 {
            return Err(QpaGbnpError::NotFreeOneVertexAlgebra);
        }
        let mut representation = path_element_to_cohen(element)?;
        for monomial in &mut representation.monomials {
            if monomial.generators() == [1] {
                *monomial = NpMonomial::from(&[][..]);
            }
        }
        result.push(representation);
    }
    Ok(result)
}

fn path_element_to_cohen<F: Field>(
    element: &PathAlgebraElement<F>,
) -> Result<CohenPolynomial<F>, QpaGbnpError> {
    let mut monomials = Vec::with_capacity(element.terms().len());
    let mut coefficients = Vec::with_capacity(element.terms().len());
    for term in element.terms() {
        monomials.push(NpMonomial::from(
            term.path()
                .external_representation()
                .into_iter()
                .map(|generator| {
                    i32::try_from(generator)
                        .map_err(|_| QpaGbnpError::GeneratorIndexTooLarge(generator))
                })
                .collect::<Result<Vec<_>, _>>()?,
        ));
        coefficients.push(term.coefficient().clone());
    }
    Ok(CohenPolynomial {
        field: element.algebra().field().clone(),
        monomials,
        coefficients,
    })
}

/// Translate `QPA_Cohen2Path`.
///
/// # Errors
///
/// Returns an error for an invalid NP generator number or path-algebra
/// arithmetic failure.
pub fn cohen_to_path<F: Field, R: CohenRepresentation<F>>(
    representations: &[R],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, QpaGbnpError> {
    let generators = algebra.generators();
    let mut elements = Vec::with_capacity(representations.len());
    for representation in representations {
        if representation.coefficient_field() != algebra.field() {
            return Err(QpaGbnpError::DifferentFields);
        }
        let mut polynomial = algebra.zero();
        for (monomial, coefficient) in representation.term_pairs() {
            let mut word = algebra.one().scale_left(coefficient)?;
            for &generator in monomial.generators() {
                let one_based = usize::try_from(generator)
                    .map_err(|_| QpaGbnpError::InvalidGenerator(generator))?;
                if one_based == 0 {
                    return Err(QpaGbnpError::InvalidGenerator(generator));
                }
                let factor = generators
                    .get(one_based - 1)
                    .ok_or(QpaGbnpError::InvalidGenerator(generator))?;
                word = word.multiply(factor)?;
            }
            polynomial = polynomial.add(&word)?;
        }
        elements.push(polynomial);
    }
    Ok(elements)
}

/// Translate `QPA_RelationsForPathAlgebra` into raw Cohen polynomials.
///
/// # Errors
///
/// Returns an error if a quiver generator index cannot fit NP format.
pub fn relations_for_path_algebra<F: Field>(
    algebra: &PathAlgebra<F>,
) -> Result<Vec<CohenPolynomial<F>>, QpaGbnpError> {
    let quiver = algebra.quiver();
    let vertices = quiver.number_of_vertices();
    let mut relations = Vec::new();
    let one = algebra.field().one();
    let negative_one = algebra.field().neg(&one);

    let mut identity_words = (1..=vertices)
        .map(|index| word(&[index]))
        .collect::<Result<Vec<_>, _>>()?;
    identity_words.push(NpMonomial::from(&[][..]));
    let mut identity_coefficients = vec![one.clone(); vertices];
    identity_coefficients.push(negative_one.clone());
    relations.push(relation(algebra, identity_words, identity_coefficients));

    for vertex in 1..=vertices {
        relations.push(relation(
            algebra,
            vec![word(&[vertex, vertex])?, word(&[vertex])?],
            vec![one.clone(), negative_one.clone()],
        ));
    }
    for left in 1..=vertices {
        for right in left + 1..=vertices {
            relations.push(relation(
                algebra,
                vec![word(&[left, right])?],
                vec![one.clone()],
            ));
            relations.push(relation(
                algebra,
                vec![word(&[right, left])?],
                vec![one.clone()],
            ));
        }
    }

    for arrow in quiver.arrows() {
        let generator = arrow.generator_index();
        let source = arrow.source().index().get();
        let target = arrow.target().index().get();
        relations.push(relation(
            algebra,
            vec![word(&[source, generator])?, word(&[generator])?],
            vec![one.clone(), negative_one.clone()],
        ));
        relations.push(relation(
            algebra,
            vec![word(&[generator, target])?, word(&[generator])?],
            vec![one.clone(), negative_one.clone()],
        ));
    }

    let arrows = quiver.arrows().collect::<Vec<_>>();
    for left in &arrows {
        for right in &arrows {
            if left.target() != right.source() {
                relations.push(relation(
                    algebra,
                    vec![word(&[left.generator_index(), right.generator_index()])?],
                    vec![one.clone()],
                ));
            }
        }
    }
    Ok(relations)
}

fn word(generators: &[usize]) -> Result<NpMonomial, QpaGbnpError> {
    Ok(NpMonomial::from(
        generators
            .iter()
            .map(|&generator| {
                i32::try_from(generator)
                    .map_err(|_| QpaGbnpError::GeneratorIndexTooLarge(generator))
            })
            .collect::<Result<Vec<_>, _>>()?,
    ))
}

fn relation<F: Field>(
    algebra: &PathAlgebra<F>,
    monomials: Vec<NpMonomial>,
    coefficients: Vec<F::Element>,
) -> CohenPolynomial<F> {
    CohenPolynomial {
        field: algebra.field().clone(),
        monomials,
        coefficients,
    }
}

/// Translate checked QPA `GBNPGroebnerBasis` for arrow-ideal relations.
///
/// When the quiver has both internal and crossing component relations,
/// complete the internal ideal first and reduce the crossing relations by it.
/// This structural preconditioning preserves the ideal and the final reduced
/// basis, and requires no caller-supplied partition or presentation metadata.
///
/// # Errors
///
/// Returns an error for elements outside the algebra/arrow ideal or any GBNP
/// conversion and computation failure.
pub fn gbnp_groebner_basis<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, QpaGbnpError> {
    if !qpa_in_arrow_ideal(elements, algebra)? {
        return Err(QpaGbnpError::OutsideArrowIdeal);
    }
    let uniform = make_uniform_list(elements, algebra)?;
    let input = if let Some(input) = components::precondition(&uniform, algebra)? {
        input
    } else {
        path_to_cohen(&uniform)?
            .iter()
            .map(CohenPolynomial::clean)
            .collect::<Result<Vec<_>, _>>()?
    };
    cohen_to_path(&sgrobner(&input)?, algebra)
}

/// Translate QPA's `GBNPGroebnerBasis(I)` ideal overload.
///
/// The source reads `LeftActingRingOfIdeal(I)` and `GeneratorsOfIdeal(I)`, then
/// dispatches unchanged to the checked list/path-algebra overload.
///
/// # Errors
///
/// Propagates the checked QPA/GBNP conversion and computation failures.
pub fn gbnp_groebner_basis_of_ideal<F: Field>(
    ideal: &PathAlgebraIdeal<F>,
) -> Result<Vec<PathAlgebraElement<F>>, QpaGbnpError> {
    gbnp_groebner_basis(ideal.generators(), ideal.algebra())
}

/// Translate QPA `GBNPGroebnerBasisNC`, including the source's structural path
/// algebra relations for quivers with more than one vertex.
///
/// # Errors
///
/// Returns an error for elements outside the algebra or any GBNP conversion
/// and computation failure.
pub fn gbnp_groebner_basis_nc<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<Vec<PathAlgebraElement<F>>, QpaGbnpError> {
    ensure_elements(elements, algebra)?;
    let uniform = make_uniform_list(elements, algebra)?;
    let input = if algebra.quiver().number_of_vertices() > 1 {
        let mut converted = path_to_cohen(&uniform)?
            .iter()
            .map(CohenPolynomial::clean)
            .collect::<Result<Vec<_>, _>>()?;
        converted.extend(
            relations_for_path_algebra(algebra)?
                .iter()
                .map(CohenPolynomial::clean)
                .collect::<Result<Vec<_>, _>>()?,
        );
        converted
    } else {
        path_to_cohen_free(&uniform)?
            .iter()
            .map(CohenPolynomial::clean)
            .collect::<Result<Vec<_>, _>>()?
    };
    let result = cohen_to_path(&sgrobner(&input)?, algebra)?;
    Ok(result
        .into_iter()
        .filter(|element| !element.is_zero())
        .collect())
}

fn ensure_elements<F: Field>(
    elements: &[PathAlgebraElement<F>],
    algebra: &PathAlgebra<F>,
) -> Result<(), QpaGbnpError> {
    if elements
        .iter()
        .all(|element| element.algebra().is_identical(algebra))
    {
        Ok(())
    } else {
        Err(QpaGbnpError::DifferentPathAlgebras)
    }
}

/// Failure in QPA's source-translated GBNP adapter.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum QpaGbnpError {
    /// A QPA generator index does not fit GBNP's signed word entries.
    #[error("QPA generator index {0} does not fit GBNP's NP representation")]
    GeneratorIndexTooLarge(usize),
    /// An NP word contains zero, a negative module generator, or no such QPA
    /// algebra generator.
    #[error("NP generator {0} is not a generator of the target path algebra")]
    InvalidGenerator(i32),
    /// `QPA_Path2CohenFree` requires exactly one quiver vertex.
    #[error("QPA_Path2CohenFree requires a one-vertex path algebra")]
    NotFreeOneVertexAlgebra,
    /// NP and path-algebra coefficients use different Numerica fields.
    #[error("NP and path-algebra coefficients use different fields")]
    DifferentFields,
    /// Checked `GBNPGroebnerBasis` accepts only arrow-ideal elements.
    #[error("the relations are not all in QPA's arrow ideal")]
    OutsideArrowIdeal,
    /// Relations do not all belong to the supplied path algebra family.
    #[error("the relations do not all belong to the supplied path algebra")]
    DifferentPathAlgebras,
    /// QPA would index an unbound coefficient-list position.
    #[error("Cohen representation has {monomials} monomials but only {coefficients} coefficients")]
    MissingCohenCoefficient {
        /// Number of monomial entries.
        monomials: usize,
        /// Number of coefficient entries.
        coefficients: usize,
    },
    /// QPA high-level validation/uniformization failed.
    #[error(transparent)]
    Groebner(#[from] GroebnerError),
    /// Path-algebra arithmetic failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// NP cleaning failed.
    #[error(transparent)]
    Np(#[from] NpError),
    /// GBNP `SGrobner` failed.
    #[error(transparent)]
    Sgrobner(#[from] SgrobnerError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Set,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::{Q, RationalField},
    };
    use qpa_core::{
        ArrowSpec, PathAlgebra, PathAlgebraElement, PathAlgebraIdeal, Quiver, VertexIndex,
        make_uniform, make_uniform_list, qpa_in_arrow_ideal,
    };
    use serde_json::{Value, json};

    use super::{
        CohenPolynomial, cohen_to_path, gbnp_groebner_basis, gbnp_groebner_basis_nc,
        gbnp_groebner_basis_of_ideal, path_to_cohen, path_to_cohen_free,
        relations_for_path_algebra,
    };

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    struct OneVertexCase<F: Field> {
        algebra: PathAlgebra<F>,
        vertex: PathAlgebraElement<F>,
        a: PathAlgebraElement<F>,
        b: PathAlgebraElement<F>,
        first: PathAlgebraElement<F>,
        second: PathAlgebraElement<F>,
    }

    fn one_vertex_case<F: Field>(field: F) -> OneVertexCase<F> {
        let quiver = Quiver::from_indexed(
            1,
            vec![
                ArrowSpec::named(vertex(1), vertex(1), "a"),
                ArrowSpec::named(vertex(1), vertex(1), "b"),
            ],
        )
        .unwrap();
        let algebra = PathAlgebra::new(field, quiver).unwrap();
        let generators = algebra.generators();
        let vertex = generators[0].clone();
        let a = generators[1].clone();
        let b = generators[2].clone();
        let first = a
            .multiply(&a)
            .unwrap()
            .multiply(&b)
            .unwrap()
            .sub(&vertex)
            .unwrap();
        let second = a
            .multiply(&b)
            .unwrap()
            .multiply(&b)
            .unwrap()
            .sub(&vertex)
            .unwrap();
        OneVertexCase {
            algebra,
            vertex,
            a,
            b,
            first,
            second,
        }
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the QPA/GBNP fixture uses only small integers")
    }

    fn cohen_json<F: Field>(
        polynomial: &CohenPolynomial<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            polynomial
                .monomials()
                .iter()
                .zip(polynomial.coefficients())
                .map(|(monomial, coefficient)| {
                    json!({
                        "monomial": monomial.generators(),
                        "coefficient": coefficient_json(coefficient)
                    })
                })
                .collect(),
        )
    }

    fn cohen_list_json<F: Field>(
        polynomials: &[CohenPolynomial<F>],
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            polynomials
                .iter()
                .map(|polynomial| cohen_json(polynomial, coefficient_json))
                .collect(),
        )
    }

    fn element_json<F: Field>(
        element: &PathAlgebraElement<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            element
                .terms()
                .iter()
                .map(|term| {
                    json!({
                        "path": term.path().external_representation(),
                        "coefficient": coefficient_json(term.coefficient())
                    })
                })
                .collect(),
        )
    }

    fn element_list_json<F: Field>(
        elements: &[PathAlgebraElement<F>],
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        Value::Array(
            elements
                .iter()
                .map(|element| element_json(element, coefficient_json))
                .collect(),
        )
    }

    fn rational_coefficient_json(coefficient: &<RationalField as Set>::Element) -> Value {
        json!([
            small_integer(coefficient.numerator_ref()),
            small_integer(coefficient.denominator_ref())
        ])
    }

    #[test]
    fn one_vertex_nc_adapter_reproduces_the_gbnp_example_in_qpa_generators() {
        let case = one_vertex_case(Q);

        let raw = path_to_cohen(std::slice::from_ref(&case.first)).unwrap();
        assert_eq!(raw[0].monomials()[0].generators(), [1]);
        let free = path_to_cohen_free(std::slice::from_ref(&case.first)).unwrap();
        assert!(free[0].monomials()[0].generators().is_empty());

        let basis = gbnp_groebner_basis_nc(&[case.first, case.second], &case.algebra).unwrap();
        assert_eq!(basis.len(), 2);
        assert_eq!(basis[0].leading_monomial().to_string(), "b");
        assert_eq!(basis[1].leading_monomial().to_string(), "a^3");
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn pinned_qpa_gbnp_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/qpa_gbnp.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(
            expected["oracle"],
            json!({
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            })
        );

        let rational = one_vertex_case(Q);
        let rational_relations = [rational.first.clone(), rational.second.clone()];
        let rational_raw = path_to_cohen(&rational_relations).unwrap();
        let rational_free = path_to_cohen_free(&rational_relations).unwrap();
        let rational_nc = gbnp_groebner_basis_nc(&rational_relations, &rational.algebra).unwrap();
        let checked_input = [
            rational.a.multiply(&rational.a).unwrap(),
            rational.a.multiply(&rational.b).unwrap(),
        ];
        let rational_checked = gbnp_groebner_basis(&checked_input, &rational.algebra).unwrap();
        let checked_ideal =
            PathAlgebraIdeal::new(rational.algebra.clone(), checked_input.to_vec()).unwrap();
        let rational_ideal = gbnp_groebner_basis_of_ideal(&checked_ideal).unwrap();
        let vertex_is_in_arrow_ideal =
            qpa_in_arrow_ideal(std::slice::from_ref(&rational.vertex), &rational.algebra).unwrap();

        let field5 = Zp::new(5);
        let f5 = one_vertex_case(field5.clone());
        let f5_relations = [f5.first.clone(), f5.second.clone()];
        let f5_free = path_to_cohen_free(&f5_relations).unwrap();
        let f5_nc = gbnp_groebner_basis_nc(&f5_relations, &f5.algebra).unwrap();
        let f5_coefficient_json = |coefficient: &<Zp as Set>::Element| {
            small_integer(&field5.to_integer(coefficient)).into()
        };

        let f2 = one_vertex_case(Z2);
        let f2_relations = [f2.first.clone(), f2.second.clone()];
        let f2_free = path_to_cohen_free(&f2_relations).unwrap();
        let f2_nc = gbnp_groebner_basis_nc(&f2_relations, &f2.algebra).unwrap();
        let f2_coefficient_json =
            |coefficient: &<Z2 as Set>::Element| small_integer(&Z2.to_integer(coefficient)).into();

        let line =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "c")]).unwrap();
        let line_algebra = PathAlgebra::new(Q, line).unwrap();
        let line_generators = line_algebra.generators();
        let c = line_generators[2].clone();
        let line_nonuniform = line_generators[0]
            .add(&c)
            .unwrap()
            .add(&line_generators[1])
            .unwrap();
        let uniform_element = make_uniform(&line_nonuniform, &line_algebra).unwrap();
        let uniform_list = make_uniform_list(&[line_nonuniform, c.clone()], &line_algebra).unwrap();
        let structure = relations_for_path_algebra(&line_algebra).unwrap();
        let structure_as_paths = cohen_to_path(&structure, &line_algebra).unwrap();
        let line_nc = gbnp_groebner_basis_nc(std::slice::from_ref(&c), &line_algebra).unwrap();

        assert_eq!(
            json!({
                "schema_version": 1,
                "kind": "qpa_gbnp",
                "oracle": expected["oracle"].clone(),
                "rational": {
                    "raw": cohen_list_json(&rational_raw, &rational_coefficient_json),
                    "free": cohen_list_json(&rational_free, &rational_coefficient_json),
                    "nc_basis": element_list_json(&rational_nc, &rational_coefficient_json),
                    "checked_basis": element_list_json(
                        &rational_checked,
                        &rational_coefficient_json
                    ),
                    "ideal_basis": element_list_json(
                        &rational_ideal,
                        &rational_coefficient_json
                    ),
                    "vertex_is_in_arrow_ideal": vertex_is_in_arrow_ideal
                },
                "f5": {
                    "free": cohen_list_json(&f5_free, &f5_coefficient_json),
                    "nc_basis": element_list_json(&f5_nc, &f5_coefficient_json)
                },
                "f2": {
                    "free": cohen_list_json(&f2_free, &f2_coefficient_json),
                    "nc_basis": element_list_json(&f2_nc, &f2_coefficient_json)
                },
                "multi_vertex": {
                    "structure_relations": cohen_list_json(
                        &structure,
                        &rational_coefficient_json
                    ),
                    "uniform_element": element_list_json(
                        &uniform_element,
                        &rational_coefficient_json
                    ),
                    "uniform_list": element_list_json(
                        &uniform_list,
                        &rational_coefficient_json
                    ),
                    "structure_as_paths": element_list_json(
                        &structure_as_paths,
                        &rational_coefficient_json
                    ),
                    "nc_basis": element_list_json(&line_nc, &rational_coefficient_json)
                }
            }),
            expected
        );
    }
}
