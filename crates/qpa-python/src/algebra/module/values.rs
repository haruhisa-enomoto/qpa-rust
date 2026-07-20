//! Python methods for module elements and basis-coordinate mappings.

use pyo3::pymethods;

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

#[pymethods]
impl PyModuleElement {
    #[getter]
    fn module(&self) -> PyMatrixModule {
        let kind = map_ref!(
            &self.kind,
            ModuleElementKind => MatrixModuleKind,
            |element| element.module().clone()
        );
        PyMatrixModule { kind }
    }

    #[getter]
    fn external_representation(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, ModuleElementKind, |element| scalar_rows(
            element.module().right_acting_algebra().field(),
            element.external_representation(),
        ))
    }

    #[getter]
    fn basis_coefficients(&self) -> Vec<PyScalar> {
        dispatch_ref!(&self.kind, ModuleElementKind, |element| {
            let field = element.module().right_acting_algebra().field();
            element
                .basis_coefficients()
                .iter()
                .map(|value| field.wrap_scalar(value))
                .collect()
        })
    }

    fn hom_from_projective(&self) -> PyResult<Option<PyModuleHomomorphism>> {
        dispatch_ref!(&self.kind, ModuleElementKind, |element| Ok(
            qpa_gbnp::hom_from_projective(element, element.module())
                .map_err(invalid)?
                .map(wrap_module_homomorphism)
        ))
    }

    #[getter]
    fn support_vertices(&self) -> Vec<usize> {
        dispatch_ref!(&self.kind, ModuleElementKind, |element| element
            .support_vertices())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, ModuleElementKind, |element| element.is_zero())
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleElementKind => ModuleElementKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("module elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn subtract(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleElementKind => ModuleElementKind,
            |left, right| left.subtract(right).map_err(invalid)?,
            invalid("module elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn negate(&self) -> Self {
        let kind = map_ref!(
            &self.kind,
            ModuleElementKind => ModuleElementKind,
            |element| element.additive_inverse()
        );
        Self { kind }
    }

    fn scale_left(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleElementKind => ModuleElementKind,
            |element| element.scale_left(
                &element
                    .module()
                    .right_acting_algebra()
                    .field()
                    .scalar_value(scalar)?,
            )
        );
        Ok(Self { kind })
    }

    fn scale_right(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            ModuleElementKind => ModuleElementKind,
            |element| element.scale_right(
                &element
                    .module()
                    .right_acting_algebra()
                    .field()
                    .scalar_value(scalar)?,
            )
        );
        Ok(Self { kind })
    }

    fn act_by_algebra_element(&self, element: &PyAlgebraElement) -> PyResult<Self> {
        let kind = match (&self.kind, &element.kind) {
            (ModuleElementKind::Rational(left), AlgebraElementKind::Rational(right)) => {
                ModuleElementKind::Rational(
                    left.act(&QuiverAlgebraElement::Path(right.clone()))
                        .map_err(invalid)?,
                )
            }
            (ModuleElementKind::Binary(left), AlgebraElementKind::Binary(right)) => {
                ModuleElementKind::Binary(
                    left.act(&QuiverAlgebraElement::Path(right.clone()))
                        .map_err(invalid)?,
                )
            }
            (ModuleElementKind::Prime(left), AlgebraElementKind::Prime(right)) => {
                ModuleElementKind::Prime(
                    left.act(&QuiverAlgebraElement::Path(right.clone()))
                        .map_err(invalid)?,
                )
            }
            _ => {
                return Err(invalid(
                    "acting element belongs to another coefficient field",
                ));
            }
        };
        Ok(Self { kind })
    }

    fn act_by_quotient_element(&self, element: &PyQuotientElement) -> PyResult<Self> {
        let kind = match (&self.kind, &element.kind) {
            (ModuleElementKind::Rational(left), QuotientElementKind::Rational(right)) => {
                ModuleElementKind::Rational(
                    left.act(&QuiverAlgebraElement::Quotient(right.clone()))
                        .map_err(invalid)?,
                )
            }
            (ModuleElementKind::Binary(left), QuotientElementKind::Binary(right)) => {
                ModuleElementKind::Binary(
                    left.act(&QuiverAlgebraElement::Quotient(right.clone()))
                        .map_err(invalid)?,
                )
            }
            (ModuleElementKind::Prime(left), QuotientElementKind::Prime(right)) => {
                ModuleElementKind::Prime(
                    left.act(&QuiverAlgebraElement::Quotient(right.clone()))
                        .map_err(invalid)?,
                )
            }
            _ => {
                return Err(invalid(
                    "acting element belongs to another coefficient field",
                ));
            }
        };
        Ok(Self { kind })
    }

    fn __add__(&self, right: &Self) -> PyResult<Self> {
        self.add(right)
    }

    fn __sub__(&self, right: &Self) -> PyResult<Self> {
        self.subtract(right)
    }

    fn __neg__(&self) -> Self {
        self.negate()
    }

    fn __eq__(&self, right: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            ModuleElementKind,
            |left, right| left == right,
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, ModuleElementKind, |element| {
            let field = element.module().right_acting_algebra().field();
            element.qpa_print_string_with(|value| field.scalar_text(value))
        })
    }
}

#[pymethods]
impl PyModuleBasisMapping {
    #[getter]
    fn source(&self) -> PyMatrixModule {
        let kind = map_ref!(
            &self.kind,
            ModuleBasisMappingKind => MatrixModuleKind,
            |mapping| mapping.source().clone()
        );
        PyMatrixModule { kind }
    }

    #[getter]
    fn range(&self) -> PyMatrixModule {
        let kind = map_ref!(
            &self.kind,
            ModuleBasisMappingKind => MatrixModuleKind,
            |mapping| mapping.range().clone()
        );
        PyMatrixModule { kind }
    }

    #[getter]
    fn basis_images(&self) -> Vec<PyModuleElement> {
        dispatch_ref!(&self.kind, ModuleBasisMappingKind, |mapping| {
            wrap_module_elements(mapping.basis_images().to_vec())
        })
    }

    #[getter]
    fn is_surjective(&self) -> bool {
        dispatch_ref!(&self.kind, ModuleBasisMappingKind, |mapping| mapping
            .is_surjective())
    }

    fn image(&self, element: &PyModuleElement) -> PyResult<PyModuleElement> {
        match (&self.kind, &element.kind) {
            (ModuleBasisMappingKind::Rational(mapping), ModuleElementKind::Rational(element)) => {
                Ok(wrap_module_element(
                    mapping.image(element).map_err(invalid)?,
                ))
            }
            (ModuleBasisMappingKind::Binary(mapping), ModuleElementKind::Binary(element)) => Ok(
                wrap_module_element(mapping.image(element).map_err(invalid)?),
            ),
            (ModuleBasisMappingKind::Prime(mapping), ModuleElementKind::Prime(element)) => Ok(
                wrap_module_element(mapping.image(element).map_err(invalid)?),
            ),
            _ => Err(invalid(
                "basis mapping and element belong to different coefficient fields",
            )),
        }
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            ModuleBasisMappingKind,
            |left, right| left.is_identical(right),
            false
        )
    }
}
