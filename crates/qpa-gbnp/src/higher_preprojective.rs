//! QPA's module-based higher preprojective-algebra construction.
//!
//! The public QPA method in `pathalgpredef.gi` first builds a graded algebra
//! by structure constants and then passes it to `AlgebraAsQuiverAlgebra`.
//! This module translates that first stage literally so the transitive
//! conversion dependency can be implemented and tested separately.

use std::ops::Range;

use numerica::domains::{Field, integer::Integer};
use qpa_core::{
    MatrixError, PathAlgebraModule, PathAlgebraModuleError, QpaCoefficientField,
    StructureConstantsAlgebra, StructureConstantsError, hom_over_algebra_with_basis_function,
};
use thiserror::Error;

use crate::{
    AlgebraAsQuiverAlgebraError, AlgebraAsQuiverAlgebraResult, EnvelopingModuleError,
    algebra_as_quiver_algebra, tr_d_homomorphism_iterated, tr_d_iterated,
};

/// The graded structure-constant stage of QPA's
/// `PreprojectiveAlgebra(M, n)` method.
#[derive(Clone, Debug)]
pub struct ModulePreprojectiveStructure<F: Field> {
    algebra: StructureConstantsAlgebra<F>,
    degrees: Vec<Range<usize>>,
}

impl<F: Field> ModulePreprojectiveStructure<F> {
    /// Return the algebra produced immediately before QPA calls
    /// `AlgebraAsQuiverAlgebra`.
    #[must_use]
    pub const fn algebra(&self) -> &StructureConstantsAlgebra<F> {
        &self.algebra
    }

    /// Return the zero-based half-open coordinate blocks corresponding to
    /// `Hom(M, TrD^i(M))` after QPA removes trailing empty degrees.
    #[must_use]
    pub fn degrees(&self) -> &[Range<usize>] {
        &self.degrees
    }
}

/// Translate the structure-constant stage of QPA's
/// `PreprojectiveAlgebra(M, n)` method.
///
/// The caller must establish QPA's `IsHereditaryAlgebra(RightActingAlgebra(M))`
/// dispatch condition. The method builds the ordered Hom bases for degrees
/// `0..n`, checks that degree `n+1` vanishes, removes trailing empty Hom
/// spaces, and uses QPA's repeated `DualOfModuleHomomorphism` followed by
/// `TransposeOfModuleHomomorphism` to define every graded product.
///
/// # Errors
///
/// Reproduces the source's insufficient-bound and empty-Hom failures, rejects
/// an integer outside the iterated functors' current machine range, and
/// propagates exact Hom, functor, composition, coordinate, and tensor errors.
#[allow(clippy::too_many_lines)]
pub fn module_preprojective_structure<F: Field>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<ModulePreprojectiveStructure<F>, ModulePreprojectiveError> {
    let bound_i64 = bound.to_string().parse::<i64>().map_err(|_| {
        ModulePreprojectiveError::BoundOutsideFunctorRange {
            value: bound.clone(),
        }
    })?;

    let mut hom_spaces = Vec::new();
    if bound_i64 >= 0 {
        let degree_count = usize::try_from(bound_i64)
            .ok()
            .and_then(|value| value.checked_add(1))
            .ok_or_else(|| ModulePreprojectiveError::BoundOutsideFunctorRange {
                value: bound.clone(),
            })?;
        hom_spaces.reserve(degree_count);
        for degree in 0..degree_count {
            let exponent = i64::try_from(degree).map_err(|_| {
                ModulePreprojectiveError::BoundOutsideFunctorRange {
                    value: bound.clone(),
                }
            })?;
            let target = tr_d_iterated(module, exponent)?;
            let hom = hom_over_algebra_with_basis_function(module, &target)?
                .ok_or(ModulePreprojectiveError::InternalHomRejected)?;
            hom_spaces.push(hom);
        }
    }

    let terminal_exponent = bound_i64.checked_add(1).ok_or_else(|| {
        ModulePreprojectiveError::BoundOutsideFunctorRange {
            value: bound.clone(),
        }
    })?;
    if tr_d_iterated(module, terminal_exponent)?.dimension() > 0 {
        return Err(ModulePreprojectiveError::MustIncreaseSecondArgument);
    }

    let Some(top) = hom_spaces
        .iter()
        .rposition(|hom| !hom.basis().is_empty())
        .map(|index| index + 1)
    else {
        return Err(ModulePreprojectiveError::NoNonzeroHomDegree);
    };
    hom_spaces.truncate(top);

    let dimension = hom_spaces.iter().try_fold(0usize, |total, hom| {
        total
            .checked_add(hom.basis().len())
            .ok_or(ModulePreprojectiveError::AlgebraDimensionOverflow)
    })?;
    let mut degrees = Vec::with_capacity(top);
    let mut origin = 0usize;
    for hom in &hom_spaces {
        let target = origin
            .checked_add(hom.basis().len())
            .ok_or(ModulePreprojectiveError::AlgebraDimensionOverflow)?;
        degrees.push(origin..target);
        origin = target;
    }

    let tensor_length = dimension
        .checked_mul(dimension)
        .and_then(|square| square.checked_mul(dimension))
        .ok_or(ModulePreprojectiveError::MultiplicationTensorOverflow { dimension })?;
    let field = module.right_acting_algebra().field().clone();
    let mut products = vec![field.zero(); tensor_length];

    for left_degree in 0..top {
        for (left_local, left_map) in hom_spaces[left_degree].basis().iter().enumerate() {
            let left_index = degrees[left_degree].start + left_local;
            for right_degree in 0..top {
                if left_degree + right_degree >= top {
                    continue;
                }
                for (right_local, right_map) in hom_spaces[right_degree].basis().iter().enumerate()
                {
                    let right_index = degrees[right_degree].start + right_local;
                    let exponent = i64::try_from(left_degree).map_err(|_| {
                        ModulePreprojectiveError::BoundOutsideFunctorRange {
                            value: bound.clone(),
                        }
                    })?;
                    let translated_right = tr_d_homomorphism_iterated(right_map, exponent)?;
                    let product = left_map.compose(&translated_right)?;
                    let coefficients = hom_spaces[left_degree + right_degree]
                        .coefficients(&product)?
                        .ok_or(ModulePreprojectiveError::ProductOutsideHomBasis)?;
                    let output_degree = &degrees[left_degree + right_degree];
                    if coefficients.len() != output_degree.len() {
                        return Err(ModulePreprojectiveError::ProductCoordinateLength {
                            actual: coefficients.len(),
                            expected: output_degree.len(),
                        });
                    }
                    let product_start = (left_index * dimension + right_index) * dimension;
                    for (local, coefficient) in coefficients.into_iter().enumerate() {
                        products[product_start + output_degree.start + local] = coefficient;
                    }
                }
            }
        }
    }

    Ok(ModulePreprojectiveStructure {
        algebra: StructureConstantsAlgebra::new(field, dimension, products)?,
        degrees,
    })
}

/// Translate QPA's complete module overload `PreprojectiveAlgebra(M, n)`.
///
/// This applies [`module_preprojective_structure`] and then QPA's
/// `AlgebraAsQuiverAlgebra`, returning the latter method's quiver algebra and
/// generator images. As in the source method's dispatch, the caller must have
/// established that the module's right-acting algebra is hereditary.
///
/// # Errors
///
/// Propagates the graded Hom construction and the finite split-prime-field
/// algebra conversion.
pub fn module_preprojective_algebra<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
    bound: &Integer,
) -> Result<AlgebraAsQuiverAlgebraResult<F>, ModulePreprojectiveError> {
    let structure = module_preprojective_structure(module, bound)?;
    Ok(algebra_as_quiver_algebra(structure.algebra())?)
}

/// Failure in the structure-constant stage of QPA's module preprojective
/// algebra method.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ModulePreprojectiveError {
    /// The current iterated functor API uses a signed machine exponent.
    #[error("the preprojective degree bound is outside the supported functor range: {value}")]
    BoundOutsideFunctorRange {
        /// Rejected QPA integer.
        value: Integer,
    },
    /// QPA explicitly asks for a larger truncation bound.
    #[error("Must increase the second argument.")]
    MustIncreaseSecondArgument,
    /// QPA's `PositionNonZero` returns `fail` when every retained Hom space is empty.
    #[error("PositionNonZero returned fail for the preprojective Hom dimensions")]
    NoNonzeroHomDegree,
    /// Every internal Hom call uses modules over the same acting algebra.
    #[error("HomOverAlgebra rejected an internal preprojective Hom space")]
    InternalHomRejected,
    /// The total graded dimension overflowed Rust's addressable size.
    #[error("the preprojective structure-algebra dimension overflows usize")]
    AlgebraDimensionOverflow,
    /// The flattened cubic multiplication tensor overflowed Rust's addressable size.
    #[error("the multiplication tensor for dimension {dimension} overflows usize")]
    MultiplicationTensorOverflow {
        /// Structure-algebra dimension.
        dimension: usize,
    },
    /// A composed map must remain in QPA's corresponding Hom basis.
    #[error("a graded preprojective product lies outside its Hom basis")]
    ProductOutsideHomBasis,
    /// The coordinate closure returned a row inconsistent with its Hom basis.
    #[error("a graded product has {actual} coordinates, expected {expected}")]
    ProductCoordinateLength {
        /// Returned coordinate count.
        actual: usize,
        /// Hom-basis dimension.
        expected: usize,
    },
    /// Module dual/transpose construction failed.
    #[error(transparent)]
    Enveloping(#[from] EnvelopingModuleError),
    /// Hom construction or map composition failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// Hom-basis coordinate elimination failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// Structure-constant algebra construction failed.
    #[error(transparent)]
    Structure(#[from] StructureConstantsError),
    /// QPA's final finite-dimensional-algebra conversion failed.
    #[error(transparent)]
    Conversion(#[from] AlgebraAsQuiverAlgebraError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field,
        finite_field::{FiniteFieldCore, Z2, Zp},
        integer::Integer,
        rational::Q,
    };
    use qpa_core::{
        ArrowSpec, PathAlgebra, QpaCoefficientField, Quiver, QuiverAlgebra, VertexIndex,
        direct_sum_of_qpa_modules,
    };
    use serde_json::{Value, json};

    use super::{
        ModulePreprojectiveError, ModulePreprojectiveStructure, module_preprojective_algebra,
        module_preprojective_structure,
    };
    use crate::{algebra_as_quiver_algebra, indec_projective_modules, simple_modules};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn additive_generator_case<F: Field>(field: F) -> ModulePreprojectiveStructure<F> {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(field, quiver).unwrap());
        let simples = simple_modules(&algebra).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let module = direct_sum_of_qpa_modules(&[
            simples[0].clone(),
            simples[1].clone(),
            projectives[0].clone(),
        ])
        .unwrap()
        .unwrap();
        module_preprojective_structure(&module, &Integer::from(1)).unwrap()
    }

    fn structure_json<F: QpaCoefficientField>(
        structure: &ModulePreprojectiveStructure<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        let algebra = structure.algebra();
        let mut products = Vec::new();
        for left in 0..algebra.dimension() {
            for right in 0..algebra.dimension() {
                let coordinates = algebra.basis_product(left, right);
                if coordinates
                    .iter()
                    .all(|coefficient| algebra.field().is_zero(coefficient))
                {
                    continue;
                }
                let sparse = coordinates
                    .iter()
                    .enumerate()
                    .filter(|(_, coefficient)| !algebra.field().is_zero(coefficient))
                    .map(|(coordinate, coefficient)| {
                        json!([coordinate + 1, coefficient_json(coefficient)])
                    })
                    .collect::<Vec<_>>();
                products.push(json!({
                    "left": left + 1,
                    "right": right + 1,
                    "coordinates": sparse
                }));
            }
        }
        let sparse_basis = |basis: &qpa_core::DenseMatrix<F>| {
            (0..basis.nrows())
                .map(|row| {
                    basis
                        .row(row)
                        .iter()
                        .enumerate()
                        .filter(|(_, coefficient)| !algebra.field().is_zero(coefficient))
                        .map(|(coordinate, coefficient)| {
                            json!([coordinate + 1, coefficient_json(coefficient)])
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };
        let sparse_rows = |rows: &[Vec<F::Element>]| {
            rows.iter()
                .map(|row| {
                    row.iter()
                        .enumerate()
                        .filter(|(_, coefficient)| !algebra.field().is_zero(coefficient))
                        .map(|(coordinate, coefficient)| {
                            json!([coordinate + 1, coefficient_json(coefficient)])
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };
        let centre = algebra.centre().unwrap();
        let radical = algebra.radical().unwrap();
        let mut result = json!({
            "dimension": algebra.dimension(),
            "hom_dimensions": structure.degrees().iter().map(std::ops::Range::len).collect::<Vec<_>>(),
            "degrees": structure.degrees().iter().map(|range| (range.start + 1..=range.end).collect::<Vec<_>>()).collect::<Vec<_>>(),
            "nonzero_products": products,
            "centre_basis": sparse_basis(centre.basis()),
            "radical_basis": sparse_basis(radical.basis())
        });
        if algebra.field().qpa_elements().unwrap().is_some() {
            let quotient = algebra
                .quotient_by_ideal_basis(radical.basis().clone())
                .unwrap();
            let idempotents = quotient
                .algebra()
                .central_idempotents_split_prime()
                .unwrap();
            let lifted = quotient
                .lift_complete_set_of_orthogonal_idempotents(&idempotents)
                .unwrap();
            let object = result.as_object_mut().unwrap();
            object.insert(
                "semisimple_lifts".into(),
                json!(sparse_basis(quotient.lifts())),
            );
            object.insert(
                "semisimple_idempotents".into(),
                json!(sparse_rows(&idempotents)),
            );
            object.insert("lifted_vertices".into(), json!(sparse_rows(&lifted)));
        }
        result
    }

    fn small_integer(value: &Integer) -> i64 {
        value
            .to_string()
            .parse()
            .expect("the preprojective fixture uses only small integers")
    }

    fn final_conversion_json<F: QpaCoefficientField>(
        structure: &ModulePreprojectiveStructure<F>,
        coefficient_json: &impl Fn(&F::Element) -> Value,
    ) -> Value {
        let converted = algebra_as_quiver_algebra(structure.algebra()).unwrap();
        let quiver = converted.algebra().quiver();
        let arrows = quiver
            .arrows()
            .map(|arrow| json!([arrow.source().index().get(), arrow.target().index().get()]))
            .collect::<Vec<_>>();
        let generator_images = converted
            .generator_images()
            .iter()
            .map(|coordinates| {
                coordinates
                    .iter()
                    .enumerate()
                    .filter(|(_, coefficient)| !structure.algebra().field().is_zero(coefficient))
                    .map(|(coordinate, coefficient)| {
                        json!([coordinate + 1, coefficient_json(coefficient)])
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let relations = converted
            .algebra()
            .relations()
            .iter()
            .map(|relation| {
                relation
                    .terms()
                    .iter()
                    .map(|term| {
                        json!({
                            "walk": term.path().walk().unwrap().iter().map(|arrow| arrow.get()).collect::<Vec<_>>(),
                            "coefficient": coefficient_json(term.coefficient())
                        })
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let dimension = converted
            .algebra()
            .as_quotient()
            .unwrap()
            .canonical_basis()
            .unwrap()
            .unwrap()
            .len();
        json!({
            "dimension": dimension,
            "vertices": quiver.number_of_vertices(),
            "arrows": arrows,
            "generator_images": generator_images,
            "loewy_length": converted.loewy_length(),
            "relations": relations
        })
    }

    #[test]
    fn pinned_qpa_graded_hom_products_match() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalg_module_preprojective_structure.schema-v1.json"
        ))
        .unwrap();

        let rational = additive_generator_case(Q);
        assert_eq!(
            structure_json(&rational, &|coefficient| {
                json!([
                    small_integer(coefficient.numerator_ref()),
                    small_integer(coefficient.denominator_ref())
                ])
            }),
            expected["rational"]
        );

        let binary = additive_generator_case(Z2);
        assert_eq!(
            structure_json(&binary, &|coefficient| {
                small_integer(&Z2.to_integer(coefficient)).into()
            }),
            expected["f2"]
        );
        assert_eq!(
            final_conversion_json(&binary, &|coefficient| {
                small_integer(&Z2.to_integer(coefficient)).into()
            }),
            expected["final"]["f2"]
        );

        let field5 = Zp::new(5);
        let prime_five = additive_generator_case(field5.clone());
        assert_eq!(
            structure_json(&prime_five, &|coefficient| {
                small_integer(&field5.to_integer(coefficient)).into()
            }),
            expected["f5"]
        );
        assert_eq!(
            final_conversion_json(&prime_five, &|coefficient| {
                small_integer(&field5.to_integer(coefficient)).into()
            }),
            expected["final"]["f5"]
        );
    }

    #[test]
    fn public_module_preprojective_method_reaches_qpa_conversion() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let simples = simple_modules(&algebra).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let module = direct_sum_of_qpa_modules(&[
            simples[0].clone(),
            simples[1].clone(),
            projectives[0].clone(),
        ])
        .unwrap()
        .unwrap();
        let converted = module_preprojective_algebra(&module, &Integer::from(1)).unwrap();
        assert_eq!(converted.algebra().quiver().number_of_arrows(), 3);
        assert_eq!(converted.algebra().relations().len(), 2);
    }

    #[test]
    fn public_module_preprojective_method_preserves_qpa_finite_field_error() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());
        let simples = simple_modules(&algebra).unwrap();
        let projectives = indec_projective_modules(&algebra).unwrap();
        let module = direct_sum_of_qpa_modules(&[
            simples[0].clone(),
            simples[1].clone(),
            projectives[0].clone(),
        ])
        .unwrap()
        .unwrap();
        assert!(matches!(
            module_preprojective_algebra(&module, &Integer::from(1)),
            Err(ModulePreprojectiveError::Conversion(
                crate::AlgebraAsQuiverAlgebraError::InfiniteCoefficientField
            ))
        ));
    }

    #[test]
    fn insufficient_bound_reproduces_qpa_error() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let simple = simple_modules(&algebra).unwrap().remove(1);
        assert!(matches!(
            module_preprojective_structure(&simple, &Integer::from(0)),
            Err(ModulePreprojectiveError::MustIncreaseSecondArgument)
        ));
    }
}
