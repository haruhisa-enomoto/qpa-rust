//! Methods on QPA Gröbner-basis objects.

use pyo3::{Py, PyRef, PyRefMut, Python, pymethods};

#[allow(clippy::wildcard_imports)]
// Private runtime enums and helpers form one binding domain.
use super::*;

#[pymethods]
impl PyGroebnerBasis {
    #[staticmethod]
    fn unreduced_right(
        algebra: &PyPathAlgebra,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<Self> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &algebra.kind,
            AlgebraKind => GroebnerKind,
            |algebra| GroebnerBasis::from_unreduced_right_relations(
                algebra,
                &runtime::extract_all(relations)?,
            )
            .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    #[getter]
    fn algebra(&self) -> PyPathAlgebra {
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => AlgebraKind,
            |basis| basis.algebra().clone()
        );
        PyPathAlgebra { kind }
    }

    #[getter]
    fn relations(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| wrap_algebra_elements(
            basis.relations().to_vec()
        ))
    }

    #[getter]
    fn enumerator(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| wrap_algebra_elements(
            basis.enumerator().to_vec()
        ))
    }

    #[getter]
    fn as_list(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| wrap_algebra_elements(
            basis.as_list().to_vec()
        ))
    }

    #[getter]
    fn has_cached_enumerator(&self) -> bool {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis
            .has_cached_enumerator())
    }

    #[getter]
    fn has_cached_as_list(&self) -> bool {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis.has_cached_as_list())
    }

    #[getter]
    fn kind(&self) -> &'static str {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| match basis.kind() {
            GroebnerBasisKind::TwoSided => "two-sided",
            GroebnerBasisKind::Right => "right",
        })
    }

    #[getter]
    fn is_right(&self) -> bool {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis.kind()
            == GroebnerBasisKind::Right)
    }

    #[getter]
    fn is_complete(&self) -> bool {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis.is_complete())
    }

    #[getter]
    fn tip_reduced_property(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis
            .tip_reduced_property())
    }

    #[getter]
    fn completely_reduced_property(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis
            .completely_reduced_property())
    }

    #[getter]
    fn homogeneous_property(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis
            .homogeneous_property())
    }

    fn set_homogeneous_property(&mut self, value: bool) {
        dispatch_ref!(&mut self.kind, GroebnerKind, |basis| basis
            .set_homogeneous_property(value));
    }

    fn set_enumerator(
        &mut self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<bool> {
        let relations = clone_relations(py, relations);
        dispatch_ref!(&mut self.kind, GroebnerKind, |basis| Ok(
            basis.set_enumerator(runtime::extract_all(relations)?)
        ))
    }

    fn set_as_list(
        &mut self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<bool> {
        let relations = clone_relations(py, relations);
        dispatch_ref!(&mut self.kind, GroebnerKind, |basis| Ok(
            basis.set_as_list(runtime::extract_all(relations)?)
        ))
    }

    fn tip_reduced_basis(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => GroebnerKind,
            |basis| basis.tip_reduced_basis().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn completely_reduced_basis(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => GroebnerKind,
            |basis| basis.completely_reduced_basis().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn right_basis(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => GroebnerKind,
            |basis| basis.right_basis().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn tip_reduce(&self, element: &PyAlgebraElement) -> PyResult<PyAlgebraElement> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => AlgebraElementKind,
            |basis| basis
                .tip_reduce(runtime::Extract::extract(element)?)
                .map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn completely_reduce(&self, element: &PyAlgebraElement) -> PyResult<PyAlgebraElement> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => AlgebraElementKind,
            |basis| basis
                .completely_reduce(runtime::Extract::extract(element)?)
                .map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn remainder(&self, element: &PyAlgebraElement) -> PyResult<PyAlgebraElement> {
        let element = element.kind.clone();
        let kind = map_ref!(
            &self.kind,
            GroebnerKind => AlgebraElementKind,
            |basis| basis
                .remainder(runtime::Extract::extract(element)?)
                .map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn is_prefix_of_tip_in_tip_ideal(&self, element: &PyAlgebraElement) -> PyResult<bool> {
        let element = element.kind.clone();
        dispatch_ref!(&self.kind, GroebnerKind, |basis| {
            let element = runtime::Extract::extract(element)?;
            basis
                .is_prefix_of_tip_in_tip_ideal(&element)
                .map_err(invalid)
        })
    }

    #[getter]
    fn admits_finitely_many_nontips(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis
            .admits_finitely_many_nontips()
            .map_err(invalid))
    }

    #[getter]
    fn nontip_size(&self) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| Ok(basis
            .nontip_size()
            .map_err(invalid)?
            .map(|value| value.to_string())))
    }

    fn nontips(&self) -> PyResult<Option<Vec<PyPath>>> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| Ok(basis
            .nontips()
            .map_err(invalid)?
            .map(|paths| paths
                .into_iter()
                .map(|inner| PyPath { inner })
                .collect())))
    }

    fn iterator(&self) -> PyGroebnerBasisIterator {
        PyGroebnerBasisIterator {
            relations: self.relations(),
            position: 0,
        }
    }

    fn print(&self) -> PyResult<String> {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis
            .print_string_with(
                |algebra| {
                    algebra.print_obj_with(&algebra.field().python_kind().display(), None)
                },
                |element| {
                    element.print_obj_with(|value| element.algebra().field().scalar_text(value))
                },
            )
            .map_err(invalid))
    }

    fn __len__(&self) -> usize {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis.relations().len())
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, GroebnerKind, |basis| basis.view_string())
    }

    fn __iter__(&self) -> PyGroebnerBasisIterator {
        self.iterator()
    }
}

#[pymethods]
impl PyGroebnerBasisIterator {
    #[getter]
    fn is_done(&self) -> bool {
        self.position >= self.relations.len()
    }

    fn next_element(&mut self) -> Option<PyAlgebraElement> {
        let element = self.relations.get(self.position).cloned();
        if element.is_some() {
            self.position += 1;
        }
        element
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyAlgebraElement> {
        slf.next_element()
    }
}
