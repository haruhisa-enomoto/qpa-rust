//! `PyO3` bindings over the source-faithful QPA Rust surface.

#![allow(clippy::needless_pass_by_value)]
#![allow(missing_docs)] // PyO3 exception structs cannot receive rustdoc attributes.

mod algebra;
mod field;
mod quiver;
mod runtime;

use algebra::{
    PyAlgebraElement, PyGroebnerBasis, PyMatrixModule, PyModuleElement, PyModuleHomomorphism,
    PyModuleSubmodule, PyPathAlgebra, PyPathAlgebraIdeal, PyQuotientAlgebra, PyQuotientElement,
};
pub(crate) use field::{FieldKind, PyField, PyScalar, ScalarKind};
use quiver::{PyPath, PyQuiver};

use pyo3::{
    PyErr, PyResult, create_exception, exceptions::PyException, prelude::*, types::PyModule,
};

// Base exception for failures reported by the Rust QPA binding.
create_exception!(_qpa, QpaError, PyException);
// Invalid mathematical input, field ownership, or object-family exception.
create_exception!(_qpa, InvalidInputError, QpaError);
// QPA's explicit `fail` result for actions that do not satisfy relations.
create_exception!(_qpa, RelationNotSatisfiedError, QpaError);

fn invalid(error: impl std::fmt::Display) -> PyErr {
    InvalidInputError::new_err(error.to_string())
}

#[pymodule]
fn _qpa(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("QpaError", module.py().get_type::<QpaError>())?;
    module.add(
        "InvalidInputError",
        module.py().get_type::<InvalidInputError>(),
    )?;
    module.add(
        "RelationNotSatisfiedError",
        module.py().get_type::<RelationNotSatisfiedError>(),
    )?;
    module.add_class::<PyField>()?;
    module.add_class::<PyScalar>()?;
    module.add_class::<PyQuiver>()?;
    module.add_class::<PyPath>()?;
    module.add_class::<PyPathAlgebra>()?;
    module.add_class::<PyAlgebraElement>()?;
    module.add_class::<PyPathAlgebraIdeal>()?;
    module.add_class::<PyGroebnerBasis>()?;
    module.add_class::<PyQuotientAlgebra>()?;
    module.add_class::<PyQuotientElement>()?;
    module.add_class::<PyMatrixModule>()?;
    module.add_class::<PyModuleElement>()?;
    module.add_class::<PyModuleHomomorphism>()?;
    module.add_class::<PyModuleSubmodule>()?;
    Ok(())
}
