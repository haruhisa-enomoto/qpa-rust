//! QPA modules constructed from enveloping algebras.
//!
//! This begins the source-order translation of `lib/pathalgtensor.gi`'s
//! module section. The regular bimodule method deliberately uses QPA's
//! canonical-basis partition and its checked path-module constructor.

use numerica::domains::Field;
use qpa_core::{
    DenseMatrix, MatrixError, NamedArrowMatrix, PathAlgebraError, PathAlgebraModule,
    PathAlgebraModuleError, PathAlgebraTensorError, PathEndpoint, PathError, QuiverAlgebra,
    QuiverAlgebraElement, QuiverAlgebraElementError, QuiverProduct, QuiverProductError,
    QuotientError, include_in_path_algebra, project_from_product_quiver,
    right_module_over_path_algebra, right_module_over_path_algebra_with_tensor,
};
use thiserror::Error;

use crate::{
    PathalgGbnpError, TensorProductError, ensure_quotient_groebner_basis, enveloping_algebra,
    opposite_path_algebra, quotient_is_finite_dimensional,
};

/// Translate QPA's cached `AlgebraAsModuleOverEnvelopingAlgebra(A)`
/// attribute (and its `AlgebraAsModuleOfEnvelopingAlgebra` synonym).
///
/// The canonical basis is partitioned into the spaces `e_i A e_j`. For each
/// product-quiver arrow, the source method acts on those basis vectors by the
/// opposite of its first projection on the left and its second projection on
/// the right, then selects canonical-basis coordinates in the target space.
/// The resulting row-action matrices are passed through QPA's checked
/// `RightModuleOverPathAlgebra` translation.
///
/// `Ok(None)` is QPA's `fail` for an infinite-dimensional algebra (or the
/// checked constructor's theoretically possible relation failure).
///
/// # Errors
///
/// Propagates the configured GBNP finite-dimensionality/basis path,
/// enveloping-algebra construction, product projections, exact algebra and
/// matrix arithmetic, and checked module construction.
pub fn algebra_as_module_over_enveloping_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_regular_enveloping_module() {
        return Ok(Some(cached));
    }
    if algebra.regular_enveloping_module_attribute_is_fail() {
        return Ok(None);
    }

    let Some(basis_vectors) = finite_canonical_basis(algebra)? else {
        algebra.mark_regular_enveloping_module_attribute_fail();
        return Ok(None);
    };

    let enveloping = enveloping_algebra(algebra)?;
    let product = enveloping.product_quiver();
    let vertices = algebra.quiver().vertices().collect::<Vec<_>>();
    let mut vector_spaces = vec![vec![Vec::<usize>::new(); vertices.len()]; vertices.len()];

    let vertex_elements = vertices
        .iter()
        .map(|vertex| include_in_path_algebra(&vertex.path(), algebra))
        .collect::<Result<Vec<_>, PathAlgebraTensorError>>()?;
    for (left_index, left) in vertex_elements.iter().enumerate() {
        for (right_index, right) in vertex_elements.iter().enumerate() {
            for (basis_index, basis_vector) in basis_vectors.iter().enumerate() {
                let component = left.multiply(basis_vector)?.multiply(right)?;
                if !component.is_zero() {
                    vector_spaces[left_index][right_index].push(basis_index);
                }
            }
        }
    }
    let dimension_vector = vector_spaces
        .iter()
        .flat_map(|row| row.iter().map(Vec::len))
        .collect::<Vec<_>>();

    let mut module_specification = Vec::new();
    for arrow in enveloping.algebra().quiver().arrows() {
        let components = [
            project_from_product_quiver(1, &arrow.path(), product)?,
            project_from_product_quiver(2, &arrow.path(), product)?,
        ];
        let source_indices = components.each_ref().map(|component| {
            let PathEndpoint::Vertex(vertex) = component.source() else {
                unreachable!("a nonzero projected product arrow has a vertex source")
            };
            vertex.index().get() - 1
        });
        let target_indices = components.each_ref().map(|component| {
            let PathEndpoint::Vertex(vertex) = component.target() else {
                unreachable!("a nonzero projected product arrow has a vertex target")
            };
            vertex.index().get() - 1
        });
        let source_space = &vector_spaces[source_indices[0]][source_indices[1]];
        let target_space = &vector_spaces[target_indices[0]][target_indices[1]];
        if source_space.is_empty() || target_space.is_empty() {
            continue;
        }

        let left_path = components[0].quiver().opposite_path(&components[0])?;
        let left = include_in_path_algebra(&left_path, algebra)?;
        let right = include_in_path_algebra(&components[1], algebra)?;
        let rows = source_space
            .iter()
            .map(|&basis_index| {
                let image = left
                    .multiply(&basis_vectors[basis_index])?
                    .multiply(&right)?;
                let coefficients = canonical_basis_coefficients(algebra, &basis_vectors, &image)?;
                Ok(target_space
                    .iter()
                    .map(|&target_index| coefficients[target_index].clone())
                    .collect::<Vec<_>>())
            })
            .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
        let matrix = DenseMatrix::from_rows(rows, algebra.field().clone())?;
        module_specification.push(NamedArrowMatrix::new(arrow.name(), matrix));
    }

    let Some(module) = right_module_over_path_algebra_with_tensor(
        enveloping.algebra().clone(),
        dimension_vector,
        &module_specification,
        enveloping,
    )?
    else {
        algebra.mark_regular_enveloping_module_attribute_fail();
        return Ok(None);
    };
    Ok(Some(algebra.cache_regular_enveloping_module(&module)))
}

/// QPA's spelling synonym for
/// [`algebra_as_module_over_enveloping_algebra`].
///
/// # Errors
///
/// Propagates the identical source method failures.
pub fn algebra_as_module_of_enveloping_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    algebra_as_module_over_enveloping_algebra(algebra)
}

/// Translate QPA's cached `DualOfModule(M)` for path-algebra matrix modules.
///
/// Arrow matrices are transposed without reordering, the dimension vector is
/// retained, and the result acts over QPA's opposite algebra. For a nonzero
/// module, QPA explicitly installs the reciprocal dual attribute on the
/// result; its zero-module early return does not.
///
/// # Errors
///
/// Propagates opposite-algebra construction and the checked right-module
/// constructor. A generated action failing the opposite relations is reported
/// as an internal source-translation failure.
pub fn dual_of_module<F: Field>(
    module: &PathAlgebraModule<F>,
) -> Result<PathAlgebraModule<F>, ModuleFunctorError> {
    if let Some(cached) = module.cached_dual_module() {
        return Ok(cached);
    }
    let opposite = match module.right_acting_algebra() {
        QuiverAlgebra::Path(algebra) => QuiverAlgebra::Path(algebra.opposite()?),
        QuiverAlgebra::Quotient(algebra) => {
            QuiverAlgebra::Quotient(opposite_path_algebra(algebra)?)
        }
    };
    let dimensions = module.dimension_vector().to_vec();
    let mut specification = Vec::new();
    if module.dimension() != 0 {
        for (position, arrow) in opposite.quiver().arrows().enumerate() {
            let source = dimensions[arrow.source().index().get() - 1];
            let target = dimensions[arrow.target().index().get() - 1];
            if source != 0 && target != 0 {
                specification.push(NamedArrowMatrix::new(
                    arrow.name(),
                    module.arrow_matrices()[position].transpose(),
                ));
            }
        }
    }
    let dual = right_module_over_path_algebra(opposite, dimensions, &specification)?
        .ok_or(ModuleFunctorError::GeneratedDualViolatesRelations)?;
    if module.dimension() == 0 {
        Ok(module.cache_dual_module(&dual))
    } else {
        Ok(module.install_dual_module(&dual)?)
    }
}

/// Translate QPA's cached
/// `DualOfAlgebraAsModuleOverEnvelopingAlgebra(A)` attribute.
///
/// QPA first dualizes the regular bimodule over the opposite enveloping
/// algebra. It then identifies that opposite with `A^op tensor A` by swapping
/// the two product-quiver coordinates, permutes vertex spaces and arrow
/// matrices through the generated names, and reconstructs a checked module
/// over the original enveloping algebra.
///
/// `Ok(None)` is the cached GAP `fail` value for an infinite-dimensional
/// algebra or an unexpected checked-relation failure.
///
/// # Errors
///
/// Propagates finite-dimensionality, enveloping/dual construction, product
/// projection and name-map failures, and checked module construction.
pub fn dual_of_algebra_as_module_over_enveloping_algebra<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<PathAlgebraModule<F>>, EnvelopingModuleError> {
    if let Some(cached) = algebra.cached_dual_enveloping_module() {
        return Ok(Some(cached));
    }
    if algebra.dual_enveloping_module_attribute_is_fail() {
        return Ok(None);
    }
    if !is_finite_dimensional(algebra)? {
        algebra.mark_dual_enveloping_module_attribute_fail();
        return Ok(None);
    }

    let enveloping = enveloping_algebra(algebra)?;
    let Some(regular) = algebra_as_module_over_enveloping_algebra(algebra)? else {
        algebra.mark_dual_enveloping_module_attribute_fail();
        return Ok(None);
    };
    let dual = dual_of_module(&regular)?;
    let product = enveloping.product_quiver();
    let vertices = enveloping.algebra().quiver().vertices().collect::<Vec<_>>();
    let arrows = enveloping.algebra().quiver().arrows().collect::<Vec<_>>();
    let vertex_names = vertices
        .iter()
        .map(qpa_core::Vertex::name)
        .collect::<Vec<_>>();
    let arrow_names = arrows.iter().map(qpa_core::Arrow::name).collect::<Vec<_>>();

    let vertex_positions = vertices
        .iter()
        .map(|vertex| {
            let name = swapped_product_generator_name(&vertex.path(), product, algebra)?;
            vertex_names
                .iter()
                .position(|candidate| *candidate == name.as_ref())
                .ok_or(EnvelopingModuleError::MissingSwappedProductGenerator { name })
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let arrow_positions = arrows
        .iter()
        .map(|arrow| {
            let name = swapped_product_generator_name(&arrow.path(), product, algebra)?;
            arrow_names
                .iter()
                .position(|candidate| *candidate == name.as_ref())
                .ok_or(EnvelopingModuleError::MissingSwappedProductGenerator { name })
        })
        .collect::<Result<Vec<_>, EnvelopingModuleError>>()?;
    let dimensions = vertex_positions
        .iter()
        .map(|&position| dual.dimension_vector()[position])
        .collect::<Vec<_>>();
    let matrices = arrow_positions
        .iter()
        .map(|&position| dual.arrow_matrices()[position].clone())
        .collect::<Vec<_>>();
    let specification = arrows
        .iter()
        .enumerate()
        .filter(|(_, arrow)| {
            dimensions[arrow.source().index().get() - 1] != 0
                && dimensions[arrow.target().index().get() - 1] != 0
        })
        .map(|(position, arrow)| NamedArrowMatrix::new(arrow.name(), matrices[position].clone()))
        .collect::<Vec<_>>();

    let Some(module) = right_module_over_path_algebra_with_tensor(
        enveloping.algebra().clone(),
        dimensions,
        &specification,
        enveloping,
    )?
    else {
        algebra.mark_dual_enveloping_module_attribute_fail();
        return Ok(None);
    };
    Ok(Some(algebra.cache_dual_enveloping_module(&module)))
}

fn is_finite_dimensional<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<bool, EnvelopingModuleError> {
    match algebra {
        QuiverAlgebra::Path(algebra) => Ok(algebra.has_finite_canonical_basis()),
        QuiverAlgebra::Quotient(algebra) => Ok(quotient_is_finite_dimensional(algebra)?),
    }
}

fn swapped_product_generator_name<F: Field>(
    generator: &qpa_core::Path,
    product: &QuiverProduct,
    algebra: &QuiverAlgebra<F>,
) -> Result<Box<str>, EnvelopingModuleError> {
    let first = project_from_product_quiver(1, generator, product)?;
    let second = project_from_product_quiver(2, generator, product)?;
    let first_name = first.quiver().opposite_quiver_name(&first.to_string())?;
    let second_name = algebra.quiver().opposite_quiver_name(&second.to_string())?;
    Ok(format!("{second_name}_{first_name}").into_boxed_str())
}

fn finite_canonical_basis<F: Field>(
    algebra: &QuiverAlgebra<F>,
) -> Result<Option<Vec<QuiverAlgebraElement<F>>>, EnvelopingModuleError> {
    match algebra {
        QuiverAlgebra::Path(path_algebra) => Ok(path_algebra
            .canonical_basis()
            .map(|basis| basis.into_iter().map(QuiverAlgebraElement::Path).collect())),
        QuiverAlgebra::Quotient(quotient) => {
            if !quotient_is_finite_dimensional(quotient)? {
                return Ok(None);
            }
            ensure_quotient_groebner_basis(quotient)?;
            Ok(Some(
                quotient
                    .canonical_basis()?
                    .ok_or(EnvelopingModuleError::MissingFiniteCanonicalBasis)?
                    .into_iter()
                    .map(QuiverAlgebraElement::Quotient)
                    .collect(),
            ))
        }
    }
}

fn canonical_basis_coefficients<F: Field>(
    algebra: &QuiverAlgebra<F>,
    basis_vectors: &[QuiverAlgebraElement<F>],
    element: &QuiverAlgebraElement<F>,
) -> Result<Vec<F::Element>, EnvelopingModuleError> {
    match (algebra, element) {
        (QuiverAlgebra::Path(path_algebra), QuiverAlgebraElement::Path(element)) => {
            let mut coordinates = vec![path_algebra.field().zero(); basis_vectors.len()];
            for term in element.terms() {
                let position = basis_vectors
                    .iter()
                    .position(|basis_vector| {
                        matches!(basis_vector, QuiverAlgebraElement::Path(basis)
                            if basis.terms().len() == 1
                                && basis.terms()[0].path() == term.path())
                    })
                    .ok_or(EnvelopingModuleError::ElementOutsideCanonicalBasis)?;
                coordinates[position] = term.coefficient().clone();
            }
            Ok(coordinates)
        }
        (QuiverAlgebra::Quotient(quotient), QuiverAlgebraElement::Quotient(element)) => {
            Ok(quotient.canonical_basis_coefficients(element)?)
        }
        _ => Err(EnvelopingModuleError::UnexpectedElementVariant),
    }
}

/// Failure while constructing QPA's regular bimodule over its enveloping
/// algebra.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum EnvelopingModuleError {
    /// A finite quotient must expose the canonical basis prepared by GBNP.
    #[error("a finite-dimensional quotient has no canonical basis vectors")]
    MissingFiniteCanonicalBasis,
    /// A path-algebra image contained a path absent from its canonical basis.
    #[error("an element path is absent from the finite canonical basis")]
    ElementOutsideCanonicalBasis,
    /// The algebra and element variants must remain paired.
    #[error("canonical-basis arithmetic returned another quiver-algebra variant")]
    UnexpectedElementVariant,
    /// The coordinate-swap name must resolve in the original product quiver.
    #[error("cannot find swapped product-quiver generator {name:?}")]
    MissingSwappedProductGenerator {
        /// Generated QPA product name.
        name: Box<str>,
    },
    /// The configured QPA/GBNP finite-basis path failed.
    #[error(transparent)]
    Backend(#[from] PathalgGbnpError),
    /// Enveloping-algebra construction failed.
    #[error(transparent)]
    Tensor(#[from] TensorProductError),
    /// Product-quiver projection failed.
    #[error(transparent)]
    Product(#[from] QuiverProductError),
    /// Opposite-path reconstruction failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// Path inclusion in the source algebra failed.
    #[error(transparent)]
    PathTensor(#[from] PathAlgebraTensorError),
    /// Source-algebra multiplication failed.
    #[error(transparent)]
    Element(#[from] QuiverAlgebraElementError),
    /// Quotient canonical-basis coordinates failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Exact matrix construction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// The checked QPA right-module constructor rejected its generated data.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// Generic QPA module dualization failed.
    #[error(transparent)]
    Functor(#[from] ModuleFunctorError),
}

/// Failure in QPA's generic module-dual functor.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ModuleFunctorError {
    /// The generated transposed actions should satisfy the opposite relations.
    #[error("the generated dual module does not satisfy the opposite-algebra relations")]
    GeneratedDualViolatesRelations,
    /// Parent path-algebra opposite construction failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Quotient opposite construction through QPA's GBNP backend failed.
    #[error(transparent)]
    Backend(#[from] PathalgGbnpError),
    /// The checked module constructor or reciprocal cache failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
}

#[cfg(test)]
mod tests {
    use numerica::domains::{
        Field, Ring,
        finite_field::{FiniteFieldCore, Z2},
        integer::Integer,
        rational::{Q, Rational},
    };
    use qpa_core::{
        ArrowSpec, DenseMatrix, NamedArrowSpec, PathAlgebra, PathAlgebraModule, Quiver,
        QuiverAlgebra, VertexIndex,
    };
    use serde_json::{Value, json};

    use super::{
        algebra_as_module_of_enveloping_algebra, algebra_as_module_over_enveloping_algebra,
        dual_of_algebra_as_module_over_enveloping_algebra, dual_of_module,
    };
    use crate::{enveloping_algebra, quotient_by_relations};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn q(value: i64) -> Rational {
        Q.nth(Integer::from(value))
    }

    #[test]
    fn regular_bimodule_for_dual_numbers_has_qpa_row_actions() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let parent = PathAlgebra::new(Q, quiver).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );

        let module = algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        let again = algebra_as_module_of_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        let enveloping = enveloping_algebra(&algebra).unwrap();

        assert_eq!(module.dimension_vector(), &[2]);
        assert_eq!(module.arrow_matrices().len(), 2);
        assert_eq!(
            module.arrow_matrices()[0].as_flat(),
            &[q(0), q(1), q(0), q(0)]
        );
        assert_eq!(
            module.arrow_matrices()[1].as_flat(),
            &[q(0), q(1), q(0), q(0)]
        );
        assert!(module.is_identical(&again));
        assert!(
            module
                .retained_tensor_product()
                .unwrap()
                .is_identical(&enveloping)
        );

        let generic_dual = dual_of_module(&module).unwrap();
        let generic_back = dual_of_module(&generic_dual).unwrap();
        let bimodule_dual = dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        assert_eq!(
            generic_dual.arrow_matrices()[0].as_flat(),
            &[q(0), q(0), q(1), q(0)]
        );
        assert!(generic_back.is_identical(&module));
        assert_eq!(
            bimodule_dual.arrow_matrices()[0].as_flat(),
            &[q(0), q(0), q(1), q(0)]
        );
    }

    #[test]
    fn regular_bimodule_supports_characteristic_two_and_zero_spaces() {
        let quiver =
            Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, quiver).unwrap());
        let module = algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();

        assert_eq!(module.dimension_vector(), &[1, 1, 0, 1]);
        assert_eq!(module.dimension(), 3);
        assert_eq!(module.arrow_matrices().len(), 4);
        assert!(module.arrow_matrices().iter().all(|matrix| {
            matrix
                .as_flat()
                .iter()
                .all(|entry| Z2.to_integer(entry) == 0 || Z2.to_integer(entry) == 1)
        }));

        let dual = dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
            .unwrap()
            .unwrap();
        assert_eq!(dual.dimension_vector(), &[1, 0, 1, 1]);
    }

    #[test]
    fn cyclic_path_algebra_returns_qpa_fail_before_enveloping_construction() {
        let quiver =
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap();
        let algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, quiver).unwrap());

        assert!(
            algebra_as_module_over_enveloping_algebra(&algebra)
                .unwrap()
                .is_none()
        );
        assert!(algebra.cached_enveloping_algebra().is_none());
        assert!(
            dual_of_algebra_as_module_over_enveloping_algebra(&algebra)
                .unwrap()
                .is_none()
        );
        assert!(algebra.has_dual_enveloping_module_attribute());
    }

    fn matrix_json<F: Field>(
        matrix: &DenseMatrix<F>,
        coefficient: impl Fn(&F::Element) -> String,
    ) -> Value {
        json!({
            "shape": [matrix.nrows(), matrix.ncols()],
            "rows": (0..matrix.nrows())
                .map(|row| matrix.row(row).iter().map(&coefficient).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
        })
    }

    fn module_json<F: Field>(
        module: &PathAlgebraModule<F>,
        coefficient: impl Fn(&F::Element) -> String + Copy,
    ) -> Value {
        json!({
            "dimension_vector": module.dimension_vector(),
            "dimension": module.dimension(),
            "arrows": module
                .right_acting_algebra()
                .quiver()
                .arrows()
                .map(|arrow| arrow.name().to_owned())
                .collect::<Vec<_>>(),
            "matrices": module
                .arrow_matrices()
                .iter()
                .map(|matrix| matrix_json(matrix, coefficient))
                .collect::<Vec<_>>(),
        })
    }

    #[test]
    fn pinned_qpa_regular_enveloping_module_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalgtensor_regular_module.schema-v1.json"
        ))
        .unwrap();

        let loop_quiver =
            Quiver::from_named(vec!["v".into()], vec![NamedArrowSpec::named("v", "v", "x")])
                .unwrap();
        let parent = PathAlgebra::new(Q, loop_quiver.clone()).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let rational_algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let rational_cached_before = rational_algebra.has_regular_enveloping_module_attribute();
        let rational_module = algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let rational_again = algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let rational_synonym = algebra_as_module_of_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let rational_enveloping = enveloping_algebra(&rational_algebra).unwrap();
        let rational_json = json!({
            "attribute_identity": [
                rational_cached_before,
                rational_algebra.has_regular_enveloping_module_attribute(),
                rational_module.is_identical(&rational_again),
                rational_module.is_identical(&rational_synonym),
                rational_module.right_acting_algebra() == rational_enveloping.algebra(),
            ],
            "module": module_json(&rational_module, ToString::to_string),
        });

        let path_quiver = Quiver::from_named(
            vec!["u1".into(), "u2".into()],
            vec![NamedArrowSpec::named("u1", "u2", "a")],
        )
        .unwrap();
        let binary_algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, path_quiver).unwrap());
        let binary_module = algebra_as_module_over_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_again = algebra_as_module_of_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_enveloping = enveloping_algebra(&binary_algebra).unwrap();
        let binary_json = json!({
            "attribute_identity": [
                binary_module.is_identical(&binary_again),
                binary_module.right_acting_algebra() == binary_enveloping.algebra(),
            ],
            "module": module_json(&binary_module, |entry| {
                if Z2.is_zero(entry) {
                    "0*Z(2)".to_owned()
                } else {
                    "Z(2)^0".to_owned()
                }
            }),
        });

        let infinite_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, loop_quiver).unwrap());
        let infinite_result = algebra_as_module_over_enveloping_algebra(&infinite_algebra).unwrap();
        let infinite_json = json!({
            "is_fail": infinite_result.is_none(),
            "attribute_bound": infinite_algebra.has_regular_enveloping_module_attribute(),
            "enveloping_bound": infinite_algebra.cached_enveloping_algebra().is_some(),
        });

        assert_eq!(rational_json, expected["rational_quotient"]);
        assert_eq!(binary_json, expected["binary_path"]);
        assert_eq!(infinite_json, expected["infinite_path"]);
    }

    #[test]
    fn pinned_qpa_dual_enveloping_module_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/pathalgtensor_dual_module.schema-v1.json"
        ))
        .unwrap();

        let loop_quiver =
            Quiver::from_named(vec!["v".into()], vec![NamedArrowSpec::named("v", "v", "x")])
                .unwrap();
        let parent = PathAlgebra::new(Q, loop_quiver.clone()).unwrap();
        let x = parent
            .element(parent.quiver().arrows().next().unwrap().path())
            .unwrap();
        let rational_algebra = QuiverAlgebra::Quotient(
            quotient_by_relations(parent, &[x.multiply(&x).unwrap()]).unwrap(),
        );
        let cached_before = rational_algebra.has_dual_enveloping_module_attribute();
        let bimodule_dual = dual_of_algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let bimodule_again = dual_of_algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let regular = algebra_as_module_over_enveloping_algebra(&rational_algebra)
            .unwrap()
            .unwrap();
        let generic_dual = dual_of_module(&regular).unwrap();
        let generic_again = dual_of_module(&regular).unwrap();
        let generic_back = dual_of_module(&generic_dual).unwrap();
        let rational_enveloping = enveloping_algebra(&rational_algebra).unwrap();
        let rational_json = json!({
            "attribute_identity": [
                cached_before,
                rational_algebra.has_dual_enveloping_module_attribute(),
                bimodule_dual.is_identical(&bimodule_again),
                bimodule_dual.right_acting_algebra() == rational_enveloping.algebra(),
            ],
            "generic_identity": [
                generic_dual.is_identical(&generic_again),
                generic_back.is_identical(&regular),
            ],
            "generic_dual": module_json(&generic_dual, ToString::to_string),
            "bimodule_dual": module_json(&bimodule_dual, ToString::to_string),
        });

        let path_quiver = Quiver::from_named(
            vec!["u1".into(), "u2".into()],
            vec![NamedArrowSpec::named("u1", "u2", "a")],
        )
        .unwrap();
        let binary_algebra = QuiverAlgebra::Path(PathAlgebra::new(Z2, path_quiver).unwrap());
        let binary_dual = dual_of_algebra_as_module_over_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_again = dual_of_algebra_as_module_over_enveloping_algebra(&binary_algebra)
            .unwrap()
            .unwrap();
        let binary_enveloping = enveloping_algebra(&binary_algebra).unwrap();
        let binary_json = json!({
            "attribute_identity": [
                binary_dual.is_identical(&binary_again),
                binary_dual.right_acting_algebra() == binary_enveloping.algebra(),
            ],
            "bimodule_dual": module_json(&binary_dual, |entry| {
                if Z2.is_zero(entry) {
                    "0*Z(2)".to_owned()
                } else {
                    "Z(2)^0".to_owned()
                }
            }),
        });

        let infinite_algebra = QuiverAlgebra::Path(PathAlgebra::new(Q, loop_quiver).unwrap());
        let infinite_result =
            dual_of_algebra_as_module_over_enveloping_algebra(&infinite_algebra).unwrap();
        let infinite_json = json!({
            "is_fail": infinite_result.is_none(),
            "attribute_bound": infinite_algebra.has_dual_enveloping_module_attribute(),
            "regular_bound": infinite_algebra.has_regular_enveloping_module_attribute(),
            "enveloping_bound": infinite_algebra.cached_enveloping_algebra().is_some(),
        });

        assert_eq!(rational_json, expected["rational_quotient"]);
        assert_eq!(binary_json, expected["binary_path"]);
        assert_eq!(infinite_json, expected["infinite_path"]);
    }
}
