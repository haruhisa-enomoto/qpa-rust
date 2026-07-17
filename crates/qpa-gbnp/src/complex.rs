//! QPA complexes and homological resolutions.
//!
//! This starts with the lazy projective-resolution specialization used by
//! `modulehomalg.gi`. Its generated differentials follow `complex.gi`'s
//! `"next/repeat"` closure exactly.

use std::sync::{Arc, Mutex};

use numerica::domains::{Field, integer::Integer};
use qpa_core::{
    PathAlgebraModule, PathAlgebraModuleHomomorphism, QpaCoefficientField, kernel_inclusion,
    zero_module, zero_module_mapping,
};

use crate::module::{EnvelopingModuleError, projective_cover};

#[derive(Debug)]
struct ProjectiveResolutionData<F: Field> {
    module: PathAlgebraModule<F>,
    zero: PathAlgebraModule<F>,
    zero_complex: bool,
    nonnegative_differentials: Mutex<Vec<PathAlgebraModuleHomomorphism<F>>>,
}

/// QPA's lazy `ProjectiveResolution(M)` complex, including `M` in degree -1.
///
/// The differential in degree zero is `ProjectiveCover(M)`. Positive
/// differentials are generated on demand from the preceding differential by
/// QPA's exact closure
/// `ProjectiveCover(Kernel(d)) * KernelInclusion(d)`. Negative degrees use the
/// zero part installed by `Complex`, with `M -> 0` in degree -1.
#[derive(Clone, Debug)]
pub struct ProjectiveResolution<F: Field>(Arc<ProjectiveResolutionData<F>>);

impl<F: QpaCoefficientField> ProjectiveResolution<F> {
    /// Return QPA's differential in an arbitrary GAP integer degree.
    ///
    /// # Errors
    ///
    /// Reports a nonnegative degree outside Rust's address space and
    /// propagates projective-cover, kernel-inclusion, composition, and lazy
    /// cache failures.
    pub fn differential(
        &self,
        degree: &Integer,
    ) -> Result<PathAlgebraModuleHomomorphism<F>, EnvelopingModuleError> {
        if self.0.zero_complex {
            return Ok(zero_module_mapping(
                self.0.zero.clone(),
                self.0.zero.clone(),
            )?);
        }
        if degree < &Integer::from(0) {
            if degree == &Integer::from(-1) {
                return Ok(zero_module_mapping(
                    self.0.module.clone(),
                    self.0.zero.clone(),
                )?);
            }
            return Ok(zero_module_mapping(
                self.0.zero.clone(),
                self.0.zero.clone(),
            )?);
        }
        let index = usize::try_from(degree.clone()).map_err(|_| {
            EnvelopingModuleError::ResolutionDegreeNotRepresentable {
                degree: degree.clone(),
            }
        })?;
        let mut differentials = self
            .0
            .nonnegative_differentials
            .lock()
            .map_err(|_| EnvelopingModuleError::ProjectiveResolutionCachePoisoned)?;
        while differentials.len() <= index {
            let previous = differentials
                .last()
                .ok_or(EnvelopingModuleError::EmptyProjectiveResolutionCache)?;
            let inclusion = kernel_inclusion(previous)?;
            let cover = projective_cover(inclusion.source())?;
            differentials.push(cover.compose(&inclusion)?);
        }
        Ok(differentials[index].clone())
    }

    /// Return QPA's `ObjectOfComplex` value, the source of the differential in
    /// the entered degree.
    ///
    /// # Errors
    ///
    /// Propagates [`Self::differential`].
    pub fn object(&self, degree: &Integer) -> Result<PathAlgebraModule<F>, EnvelopingModuleError> {
        Ok(self.differential(degree)?.source().clone())
    }
}

/// Translate QPA's `ProjectiveResolution(M)` constructor.
///
/// QPA computes the degree-zero cover eagerly, then installs the positive
/// `"next/repeat"` closure. For a zero module, `Complex` evaluates one further
/// differential and collapses the result to its canonical zero complex; this
/// branch is reproduced here.
///
/// # Errors
///
/// Propagates projective-cover, kernel-inclusion, composition, and zero-module
/// construction failures.
pub fn projective_resolution<F: QpaCoefficientField>(
    module: &PathAlgebraModule<F>,
) -> Result<ProjectiveResolution<F>, EnvelopingModuleError> {
    let cover = projective_cover(module)?;
    let zero = zero_module(module.right_acting_algebra())?;
    let zero_complex = if cover.is_zero() {
        let inclusion = kernel_inclusion(&cover)?;
        let next_cover = projective_cover(inclusion.source())?;
        next_cover.compose(&inclusion)?.is_zero()
    } else {
        false
    };
    let differentials = if zero_complex {
        Vec::new()
    } else {
        vec![cover]
    };
    Ok(ProjectiveResolution(Arc::new(ProjectiveResolutionData {
        module: module.clone(),
        zero,
        zero_complex,
        nonnegative_differentials: Mutex::new(differentials),
    })))
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
        ArrowSpec, DenseMatrix, NamedArrowMatrix, PathAlgebra, PathAlgebraModule,
        PathAlgebraModuleHomomorphism, QpaCoefficientField, Quiver, QuiverAlgebra, VertexIndex,
        right_module_over_path_algebra, zero_module,
    };
    use serde_json::{Value, json};

    use super::{ProjectiveResolution, projective_resolution};
    use crate::{EnvelopingModuleError, nth_syzygy, quotient_by_relations};

    fn vertex(index: usize) -> VertexIndex {
        VertexIndex::new(index).unwrap()
    }

    fn map_json<F: Field>(
        map: &PathAlgebraModuleHomomorphism<F>,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!({
            "source_dimension": map.source().dimension_vector(),
            "range_dimension": map.range().dimension_vector(),
            "matrices": map.vertex_matrices().iter().map(|matrix| json!(
                (0..matrix.nrows()).map(|row| matrix.row(row).iter()
                    .map(coefficient).collect::<Vec<_>>()).collect::<Vec<_>>()
            )).collect::<Vec<_>>(),
        })
    }

    fn resolution_json<F: QpaCoefficientField>(
        resolution: &ProjectiveResolution<F>,
        degrees: &[i64],
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        json!({
            "objects": degrees.iter().map(|&degree| resolution
                .object(&Integer::from(degree)).unwrap().dimension_vector().to_vec())
                .collect::<Vec<_>>(),
            "differentials": degrees.iter().map(|&degree| map_json(
                &resolution.differential(&Integer::from(degree)).unwrap(), coefficient))
                .collect::<Vec<_>>(),
        })
    }

    #[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
    fn projective_resolution_case<F: QpaCoefficientField>(
        field: F,
        coefficient: impl Fn(&F::Element) -> Value + Copy,
    ) -> Value {
        let zero = field.zero();
        let one = field.one();
        let a2 = QuiverAlgebra::Path(
            PathAlgebra::new(
                field.clone(),
                Quiver::from_indexed(2, vec![ArrowSpec::named(vertex(1), vertex(2), "a")]).unwrap(),
            )
            .unwrap(),
        );
        let a2_zero = zero_module(&a2).unwrap();
        let simple_one = right_module_over_path_algebra(a2.clone(), vec![1, 0], &[])
            .unwrap()
            .unwrap();
        let projective = right_module_over_path_algebra(
            a2,
            vec![1, 1],
            &[NamedArrowMatrix::new(
                "a",
                DenseMatrix::from_rows(vec![vec![one]], field.clone()).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        let loop_parent = PathAlgebra::new(
            field.clone(),
            Quiver::from_indexed(1, vec![ArrowSpec::named(vertex(1), vertex(1), "x")]).unwrap(),
        )
        .unwrap();
        let x = loop_parent.generators()[1].clone();
        let dual = QuiverAlgebra::Quotient(
            quotient_by_relations(loop_parent, &[x.pow(2).unwrap()]).unwrap(),
        );
        let dual_simple = right_module_over_path_algebra(
            dual,
            vec![1],
            &[NamedArrowMatrix::new(
                "x",
                DenseMatrix::from_rows(vec![vec![zero]], field).unwrap(),
            )],
        )
        .unwrap()
        .unwrap();

        let nth_zero = nth_syzygy(&simple_one, &Integer::from(0)).unwrap();
        let module_dimension = |module: PathAlgebraModule<F>| {
            json!({
                "dimension": module.dimension_vector(),
            })
        };
        json!({
            "zero_resolution": resolution_json(
                &projective_resolution(&a2_zero).unwrap(), &[-2, -1, 0, 1, 2], coefficient),
            "a2_simple_resolution": resolution_json(
                &projective_resolution(&simple_one).unwrap(), &[-2, -1, 0, 1, 2], coefficient),
            "a2_projective_resolution": resolution_json(
                &projective_resolution(&projective).unwrap(), &[-2, -1, 0, 1, 2], coefficient),
            "dual_simple_resolution": resolution_json(
                &projective_resolution(&dual_simple).unwrap(), &[-1, 0, 1, 2, 3], coefficient),
            "nth_zero_identical": nth_zero.is_identical(&simple_one),
            "a2_nth_0": module_dimension(nth_zero),
            "a2_nth_1": module_dimension(nth_syzygy(&simple_one, &Integer::from(1)).unwrap()),
            "a2_nth_2": module_dimension(nth_syzygy(&simple_one, &Integer::from(2)).unwrap()),
            "a2_nth_3": module_dimension(nth_syzygy(&simple_one, &Integer::from(3)).unwrap()),
            "dual_nth_1": module_dimension(nth_syzygy(&dual_simple, &Integer::from(1)).unwrap()),
            "dual_nth_2": module_dimension(nth_syzygy(&dual_simple, &Integer::from(2)).unwrap()),
            "dual_nth_3": module_dimension(nth_syzygy(&dual_simple, &Integer::from(3)).unwrap()),
        })
    }

    #[test]
    fn pinned_qpa_projective_resolution_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/complex_projective_resolution.schema-v1.json"
        ))
        .unwrap();
        let field5 = Zp::new(5);
        let field5_coefficients = field5.clone();
        let actual = json!({
            "schema_version": 1,
            "kind": "complex_projective_resolution",
            "oracle": {
                "gap_version": "4.17dev",
                "gap_revision": "4efe877859ab0bc0ac041bf635a996fb897ee0c9",
                "qpa_version": "1.37",
                "qpa_revision": "8d5f847258d4853700ce3a5d947db2a1441109f1",
                "gbnp_version": "1.1.0",
                "gbnp_revision": "e0b46c1aac16f3fbabe743e3169dc6bd86bdc654"
            },
            "binary": projective_resolution_case(Z2, |entry| json!(
                Z2.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "quinary": projective_resolution_case(field5, |entry| json!(
                field5_coefficients.to_integer(entry).to_string().parse::<i64>().unwrap()
            )),
            "rational": projective_resolution_case(Q, |entry| json!([
                entry.numerator_ref().to_string().parse::<i64>().unwrap(),
                entry.denominator_ref().to_string().parse::<i64>().unwrap(),
            ])),
        });

        assert_eq!(actual, expected);
    }

    #[test]
    fn nth_syzygy_rejects_negative_degree_before_constructing_resolution() {
        let algebra = QuiverAlgebra::Path(
            PathAlgebra::new(
                Z2,
                Quiver::from_indexed(1, Vec::<ArrowSpec>::new()).unwrap(),
            )
            .unwrap(),
        );
        let simple = right_module_over_path_algebra(algebra, vec![1], &[])
            .unwrap()
            .unwrap();
        assert!(matches!(
            nth_syzygy(&simple, &Integer::from(-1)),
            Err(EnvelopingModuleError::NegativeSyzygyDegree { .. })
        ));
        assert!(simple.cached_projective_cover().is_none());
    }
}
