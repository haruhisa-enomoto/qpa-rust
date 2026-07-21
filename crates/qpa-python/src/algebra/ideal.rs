//! Methods on QPA path-algebra ideals.

use pyo3::pymethods;

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

#[pymethods]
impl PyPathAlgebraIdeal {
    #[getter]
    fn algebra(&self) -> PyPathAlgebra {
        let kind = map_ref!(
            &self.kind,
            IdealKind => AlgebraKind,
            |ideal| ideal.algebra().clone()
        );
        PyPathAlgebra { kind }
    }

    #[getter]
    fn generators(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| wrap_algebra_elements(
            ideal.generators().to_vec()
        ))
    }

    #[getter]
    fn is_trivial(&self) -> bool {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal.is_trivial())
    }

    #[getter]
    fn generators_are_monomial(&self) -> bool {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .generators_are_monomial())
    }

    #[getter]
    fn attached_groebner_basis(&self) -> Option<PyGroebnerBasis> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .groebner_basis()
            .map(wrap_groebner_basis))
    }

    fn with_high_level_groebner_basis(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            IdealKind => IdealKind,
            |ideal| ideal.with_high_level_groebner_basis().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn with_complete_groebner_basis(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<Self> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(&self.kind, IdealKind => IdealKind, |ideal| {
            let relations = runtime::extract_all(relations)?;
            ideal
                .with_complete_groebner_basis(&relations)
                .map_err(invalid)?
        });
        Ok(Self { kind })
    }

    /// Invoke QPA's `GBNPGroebnerBasis` ideal overload.
    fn gbnp_groebner_basis(&self) -> PyResult<Vec<PyAlgebraElement>> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| {
            Ok(wrap_algebra_elements(
                qpa_gbnp::gbnp_groebner_basis_of_ideal(ideal).map_err(invalid)?,
            ))
        })
    }

    fn right_groebner_basis(&self) -> PyResult<PyGroebnerBasis> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| Ok(wrap_groebner_basis(
            ideal.right_groebner_basis().map_err(invalid)?,
        )))
    }

    fn contains(&self, element: &PyAlgebraElement) -> PyResult<bool> {
        let element = element.kind.clone();
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .contains(runtime::Extract::extract(element)?)
            .map_err(invalid))
    }

    fn is_monomial(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .is_monomial()
            .map_err(invalid))
    }

    fn is_admissible(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .is_admissible()
            .map_err(invalid))
    }

    fn product(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            IdealKind => IdealKind,
            |left, right| left.product(right).map_err(invalid)?,
            invalid("ideals belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn minimal_generating_set(&self) -> PyResult<Vec<PyAlgebraElement>> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| Ok(wrap_algebra_elements(
            qpa_gbnp::ideal_minimal_generating_set(ideal).map_err(invalid)?,
        )))
    }

    fn quotient(&self) -> PyResult<PyQuotientAlgebra> {
        let kind = map_ref!(
            &self.kind,
            IdealKind => QuotientKind,
            |ideal| build_ideal_quotient(ideal)?
        );
        Ok(PyQuotientAlgebra { kind })
    }

    fn quotient_without_groebner_basis(&self) -> PyResult<PyQuotientAlgebra> {
        let kind = map_ref!(
            &self.kind,
            IdealKind => QuotientKind,
            |ideal| QuotientPathAlgebra::by_ideal(ideal).map_err(invalid)?
        );
        Ok(PyQuotientAlgebra { kind })
    }

    fn natural_homomorphism(&self) -> PyResult<PyNaturalQuotientMap> {
        dispatch_ref!(&self.kind, IdealKind, |ideal| {
            mapping::natural_quotient_map_for(ideal)
        })
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            IdealKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __contains__(&self, element: &PyAlgebraElement) -> PyResult<bool> {
        self.contains(element)
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, IdealKind, |ideal| ideal
            .view_obj_with(&ideal.algebra().field().python_kind().display()))
    }

    fn view(&self) -> String {
        self.__repr__()
    }
}
