//! Python wrappers for structure-constant quotients and algebra mappings.

use pyo3::prelude::*;
use qpa_core::{
    DenseMatrix, StructureConstantsAlgebra, StructureConstantsAlgebraElement,
    StructureConstantsAlgebraHomomorphism, StructureConstantsQuotient, lift_idempotent,
    lift_idempotents_for_decomposition, lift_two_orthogonal_idempotents,
    lifting_complete_set_of_orthogonal_idempotents,
};

use super::*;

#[derive(Clone)]
pub(super) enum StructureConstantsQuotientKind {
    Rational(StructureConstantsQuotient<RationalField>),
    Binary(StructureConstantsQuotient<Z2>),
    Prime(StructureConstantsQuotient<Zp>),
}

/// A natural quotient of a structure-constants algebra by an ideal basis.
#[pyclass(name = "StructureConstantsQuotient", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyStructureConstantsQuotient {
    kind: StructureConstantsQuotientKind,
}

#[derive(Clone)]
enum StructureConstantsHomomorphismKind {
    Rational(StructureConstantsAlgebraHomomorphism<RationalField>),
    Binary(StructureConstantsAlgebraHomomorphism<Z2>),
    Prime(StructureConstantsAlgebraHomomorphism<Zp>),
}

/// A structure-constants algebra homomorphism stored by ordered basis images.
#[pyclass(name = "StructureConstantsHomomorphism", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyStructureConstantsHomomorphism {
    kind: StructureConstantsHomomorphismKind,
}

macro_rules! impl_quotient_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<StructureConstantsQuotient<$field>> for StructureConstantsQuotientKind {
            fn wrap(value: StructureConstantsQuotient<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<StructureConstantsAlgebraHomomorphism<$field>>
            for StructureConstantsHomomorphismKind
        {
            fn wrap(value: StructureConstantsAlgebraHomomorphism<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_quotient_runtime!(RationalField, Rational);
impl_quotient_runtime!(Z2, Binary);
impl_quotient_runtime!(Zp, Prime);

fn scalar_matrix<F: runtime::PythonField>(
    py: Python<'_>,
    field: &F,
    expected_rows: usize,
    expected_columns: usize,
    rows: Vec<Vec<Py<PyScalar>>>,
) -> PyResult<DenseMatrix<F>> {
    if rows.len() != expected_rows || rows.iter().any(|row| row.len() != expected_columns) {
        return Err(invalid(format!(
            "expected a {expected_rows} x {expected_columns} coefficient matrix",
        )));
    }
    let entries = rows
        .into_iter()
        .flatten()
        .map(|coefficient| field.scalar_value(&coefficient.borrow(py)))
        .collect::<PyResult<Vec<_>>>()?;
    DenseMatrix::from_flat(expected_rows, expected_columns, entries, field.clone()).map_err(invalid)
}

fn dense_rows<F: runtime::PythonField>(matrix: &DenseMatrix<F>) -> Vec<Vec<PyScalar>> {
    structure_matrix_rows(matrix)
}

fn coordinates_in<F: runtime::PythonField>(
    kind: StructureConstantsElementKind,
    algebra: &StructureConstantsAlgebra<F>,
) -> PyResult<Vec<F::Element>>
where
    StructureConstantsElementKind: runtime::Extract<StructureConstantsAlgebraElement<F>>,
{
    let element: StructureConstantsAlgebraElement<F> = runtime::Extract::extract(kind)?;
    if element.algebra() != algebra {
        return Err(invalid(
            "structure-constants element belongs to a different algebra",
        ));
    }
    Ok(element.coordinates().to_vec())
}

fn coordinate_rows_in<F: runtime::PythonField>(
    kinds: Vec<StructureConstantsElementKind>,
    algebra: &StructureConstantsAlgebra<F>,
) -> PyResult<Vec<Vec<F::Element>>>
where
    StructureConstantsElementKind: runtime::Extract<StructureConstantsAlgebraElement<F>>,
{
    kinds
        .into_iter()
        .map(|kind| coordinates_in(kind, algebra))
        .collect()
}

fn element_kinds(
    py: Python<'_>,
    elements: Vec<Py<PyStructureConstantsElement>>,
) -> Vec<StructureConstantsElementKind> {
    elements
        .into_iter()
        .map(|element| element.borrow(py).kind.clone())
        .collect()
}

pub(super) fn structure_constants_quotient_for<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: &StructureConstantsAlgebra<F>,
    ideal_basis: Vec<Vec<Py<PyScalar>>>,
) -> PyResult<PyStructureConstantsQuotient>
where
    StructureConstantsQuotientKind: runtime::Wrap<StructureConstantsQuotient<F>>,
{
    let row_count = ideal_basis.len();
    let ideal_basis = scalar_matrix(
        py,
        algebra.field(),
        row_count,
        algebra.dimension(),
        ideal_basis,
    )?;
    Ok(PyStructureConstantsQuotient {
        kind: runtime::Wrap::wrap(
            algebra
                .quotient_by_ideal_basis(ideal_basis)
                .map_err(invalid)?,
        ),
    })
}

fn homomorphism_for<F: runtime::PythonField>(
    py: Python<'_>,
    source: &StructureConstantsAlgebra<F>,
    range: &StructureConstantsAlgebra<F>,
    basis_images: Vec<Vec<Py<PyScalar>>>,
) -> PyResult<PyStructureConstantsHomomorphism>
where
    StructureConstantsHomomorphismKind: runtime::Wrap<StructureConstantsAlgebraHomomorphism<F>>,
{
    let basis_images = scalar_matrix(
        py,
        source.field(),
        source.dimension(),
        range.dimension(),
        basis_images,
    )?;
    Ok(PyStructureConstantsHomomorphism {
        kind: runtime::Wrap::wrap(
            StructureConstantsAlgebraHomomorphism::by_basis_images(
                source.clone(),
                range.clone(),
                basis_images,
            )
            .map_err(invalid)?,
        ),
    })
}

#[pymethods]
impl PyStructureConstantsQuotient {
    #[getter]
    fn original(&self) -> PyStructureConstantsAlgebra {
        let kind = map_ref!(
            &self.kind,
            StructureConstantsQuotientKind => StructureConstantsAlgebraKind,
            |quotient| quotient.original().clone()
        );
        PyStructureConstantsAlgebra { kind }
    }

    #[getter]
    fn ideal_basis(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, StructureConstantsQuotientKind, |quotient| {
            dense_rows(quotient.ideal_basis())
        })
    }

    #[getter]
    fn algebra(&self) -> PyStructureConstantsAlgebra {
        let kind = map_ref!(
            &self.kind,
            StructureConstantsQuotientKind => StructureConstantsAlgebraKind,
            |quotient| quotient.algebra().clone()
        );
        PyStructureConstantsAlgebra { kind }
    }

    #[getter]
    fn lifts(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(&self.kind, StructureConstantsQuotientKind, |quotient| {
            dense_rows(quotient.lifts())
        })
    }

    fn project(
        &self,
        element: &PyStructureConstantsElement,
    ) -> PyResult<PyStructureConstantsElement> {
        dispatch_ref!(&self.kind, StructureConstantsQuotientKind, |quotient| {
            let coordinates = coordinates_in(element.kind.clone(), quotient.original())?;
            Ok(wrap_structure_element(
                quotient
                    .algebra()
                    .element(
                        quotient
                            .project_coordinates(&coordinates)
                            .map_err(invalid)?,
                    )
                    .map_err(invalid)?,
            ))
        })
    }

    fn lift(&self, element: &PyStructureConstantsElement) -> PyResult<PyStructureConstantsElement> {
        dispatch_ref!(&self.kind, StructureConstantsQuotientKind, |quotient| {
            let coordinates = coordinates_in(element.kind.clone(), quotient.algebra())?;
            Ok(wrap_structure_element(
                quotient
                    .original()
                    .element(quotient.lift_coordinates(&coordinates).map_err(invalid)?)
                    .map_err(invalid)?,
            ))
        })
    }

    fn natural_homomorphism(&self) -> PyResult<PyStructureConstantsHomomorphism> {
        dispatch_ref!(&self.kind, StructureConstantsQuotientKind, |quotient| {
            Ok(PyStructureConstantsHomomorphism {
                kind: runtime::Wrap::wrap(quotient.natural_homomorphism().map_err(invalid)?),
            })
        })
    }

    fn lifting_idempotent(
        &self,
        idempotent: &PyStructureConstantsElement,
    ) -> PyResult<Option<PyStructureConstantsElement>> {
        dispatch_ref!(&self.kind, StructureConstantsQuotientKind, |quotient| {
            let coordinates = coordinates_in(idempotent.kind.clone(), quotient.algebra())?;
            quotient
                .lifting_idempotent(&coordinates)
                .map_err(invalid)?
                .map(|coordinates| {
                    quotient
                        .original()
                        .element(coordinates)
                        .map(wrap_structure_element)
                        .map_err(invalid)
                })
                .transpose()
        })
    }

    fn lift_idempotent(
        &self,
        idempotent: &PyStructureConstantsElement,
    ) -> PyResult<PyStructureConstantsElement> {
        dispatch_ref!(&self.kind, StructureConstantsQuotientKind, |quotient| {
            let coordinates = coordinates_in(idempotent.kind.clone(), quotient.algebra())?;
            Ok(wrap_structure_element(
                quotient
                    .original()
                    .element(quotient.lift_idempotent(&coordinates).map_err(invalid)?)
                    .map_err(invalid)?,
            ))
        })
    }

    fn lift_complete_set_of_orthogonal_idempotents(
        &self,
        py: Python<'_>,
        idempotents: Vec<Py<PyStructureConstantsElement>>,
    ) -> PyResult<Vec<PyStructureConstantsElement>> {
        let kinds = element_kinds(py, idempotents);
        dispatch_ref!(&self.kind, StructureConstantsQuotientKind, |quotient| {
            let coordinates = coordinate_rows_in(kinds, quotient.algebra())?;
            quotient
                .lift_complete_set_of_orthogonal_idempotents(&coordinates)
                .map_err(invalid)?
                .into_iter()
                .map(|coordinates| {
                    quotient
                        .original()
                        .element(coordinates)
                        .map(wrap_structure_element)
                        .map_err(invalid)
                })
                .collect()
        })
    }

    fn lift_two_orthogonal_idempotents(
        &self,
        first: &PyStructureConstantsElement,
        second: &PyStructureConstantsElement,
    ) -> PyResult<Vec<PyStructureConstantsElement>> {
        dispatch_ref!(&self.kind, StructureConstantsQuotientKind, |quotient| {
            let first = coordinates_in(first.kind.clone(), quotient.original())?;
            let second = coordinates_in(second.kind.clone(), quotient.algebra())?;
            quotient
                .lift_two_orthogonal_idempotents(&first, &second)
                .map_err(invalid)?
                .into_iter()
                .map(|coordinates| {
                    quotient
                        .original()
                        .element(coordinates)
                        .map(wrap_structure_element)
                        .map_err(invalid)
                })
                .collect()
        })
    }

    fn lift_idempotents_for_decomposition(
        &self,
        py: Python<'_>,
        idempotents: Vec<Py<PyStructureConstantsElement>>,
    ) -> PyResult<Vec<PyStructureConstantsElement>> {
        let kinds = element_kinds(py, idempotents);
        dispatch_ref!(&self.kind, StructureConstantsQuotientKind, |quotient| {
            let coordinates = coordinate_rows_in(kinds, quotient.algebra())?;
            quotient
                .lift_idempotents_for_decomposition(&coordinates)
                .map_err(invalid)?
                .into_iter()
                .map(|coordinates| {
                    quotient
                        .original()
                        .element(coordinates)
                        .map(wrap_structure_element)
                        .map_err(invalid)
                })
                .collect()
        })
    }
}

#[pymethods]
impl PyStructureConstantsHomomorphism {
    #[staticmethod]
    fn by_basis_images(
        py: Python<'_>,
        source: &PyStructureConstantsAlgebra,
        range: &PyStructureConstantsAlgebra,
        basis_images: Vec<Vec<Py<PyScalar>>>,
    ) -> PyResult<Self> {
        dispatch_pair_ref!(
            &source.kind,
            &range.kind,
            StructureConstantsAlgebraKind,
            |source, range| homomorphism_for(py, source, range, basis_images),
            Err(invalid(
                "source and range algebras use different coefficient fields"
            ))
        )
    }

    #[getter]
    fn source(&self) -> PyStructureConstantsAlgebra {
        let kind = map_ref!(
            &self.kind,
            StructureConstantsHomomorphismKind => StructureConstantsAlgebraKind,
            |homomorphism| homomorphism.source().clone()
        );
        PyStructureConstantsAlgebra { kind }
    }

    #[getter]
    fn range(&self) -> PyStructureConstantsAlgebra {
        let kind = map_ref!(
            &self.kind,
            StructureConstantsHomomorphismKind => StructureConstantsAlgebraKind,
            |homomorphism| homomorphism.range().clone()
        );
        PyStructureConstantsAlgebra { kind }
    }

    #[getter]
    fn basis_images(&self) -> Vec<Vec<PyScalar>> {
        dispatch_ref!(
            &self.kind,
            StructureConstantsHomomorphismKind,
            |homomorphism| dense_rows(homomorphism.basis_images())
        )
    }

    fn image(
        &self,
        element: &PyStructureConstantsElement,
    ) -> PyResult<PyStructureConstantsElement> {
        dispatch_ref!(
            &self.kind,
            StructureConstantsHomomorphismKind,
            |homomorphism| {
                let coordinates = coordinates_in(element.kind.clone(), homomorphism.source())?;
                Ok(wrap_structure_element(
                    homomorphism
                        .range()
                        .element(
                            homomorphism
                                .image_coordinates(&coordinates)
                                .map_err(invalid)?,
                        )
                        .map_err(invalid)?,
                ))
            }
        )
    }

    fn preimage_representative(
        &self,
        element: &PyStructureConstantsElement,
    ) -> PyResult<Option<PyStructureConstantsElement>> {
        dispatch_ref!(
            &self.kind,
            StructureConstantsHomomorphismKind,
            |homomorphism| {
                let coordinates = coordinates_in(element.kind.clone(), homomorphism.range())?;
                homomorphism
                    .preimage_representative(&coordinates)
                    .map_err(invalid)?
                    .map(|coordinates| {
                        homomorphism
                            .source()
                            .element(coordinates)
                            .map(wrap_structure_element)
                            .map_err(invalid)
                    })
                    .transpose()
            }
        )
    }

    fn lift_idempotent(
        &self,
        idempotent: &PyStructureConstantsElement,
    ) -> PyResult<PyStructureConstantsElement> {
        dispatch_ref!(
            &self.kind,
            StructureConstantsHomomorphismKind,
            |homomorphism| {
                let idempotent = coordinates_in(idempotent.kind.clone(), homomorphism.range())?;
                Ok(wrap_structure_element(
                    homomorphism
                        .source()
                        .element(lift_idempotent(homomorphism, &idempotent).map_err(invalid)?)
                        .map_err(invalid)?,
                ))
            }
        )
    }

    fn lift_idempotents_for_decomposition(
        &self,
        py: Python<'_>,
        idempotents: Vec<Py<PyStructureConstantsElement>>,
    ) -> PyResult<Vec<PyStructureConstantsElement>> {
        let kinds = element_kinds(py, idempotents);
        dispatch_ref!(
            &self.kind,
            StructureConstantsHomomorphismKind,
            |homomorphism| {
                let idempotents = coordinate_rows_in(kinds, homomorphism.range())?;
                lift_idempotents_for_decomposition(homomorphism, &idempotents)
                    .map_err(invalid)?
                    .into_iter()
                    .map(|coordinates| {
                        homomorphism
                            .source()
                            .element(coordinates)
                            .map(wrap_structure_element)
                            .map_err(invalid)
                    })
                    .collect()
            }
        )
    }

    fn lift_two_orthogonal_idempotents(
        &self,
        first: &PyStructureConstantsElement,
        second: &PyStructureConstantsElement,
    ) -> PyResult<Vec<PyStructureConstantsElement>> {
        dispatch_ref!(
            &self.kind,
            StructureConstantsHomomorphismKind,
            |homomorphism| {
                let first = coordinates_in(first.kind.clone(), homomorphism.source())?;
                let second = coordinates_in(second.kind.clone(), homomorphism.range())?;
                lift_two_orthogonal_idempotents(homomorphism, &first, &second)
                    .map_err(invalid)?
                    .into_iter()
                    .map(|coordinates| {
                        homomorphism
                            .source()
                            .element(coordinates)
                            .map(wrap_structure_element)
                            .map_err(invalid)
                    })
                    .collect()
            }
        )
    }

    fn lifting_complete_set_of_orthogonal_idempotents(
        &self,
        py: Python<'_>,
        idempotents: Vec<Py<PyStructureConstantsElement>>,
    ) -> PyResult<Vec<PyStructureConstantsElement>> {
        let kinds = element_kinds(py, idempotents);
        dispatch_ref!(
            &self.kind,
            StructureConstantsHomomorphismKind,
            |homomorphism| {
                let idempotents = coordinate_rows_in(kinds, homomorphism.range())?;
                lifting_complete_set_of_orthogonal_idempotents(homomorphism, &idempotents)
                    .map_err(invalid)?
                    .into_iter()
                    .map(|coordinates| {
                        homomorphism
                            .source()
                            .element(coordinates)
                            .map(wrap_structure_element)
                            .map_err(invalid)
                    })
                    .collect()
            }
        )
    }
}
