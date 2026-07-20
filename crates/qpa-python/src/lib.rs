//! `PyO3` bindings over the source-faithful QPA Rust surface.

#![allow(clippy::needless_pass_by_value)]
#![allow(missing_docs)] // PyO3 exception structs cannot receive rustdoc attributes.

mod algebra;
mod field;
mod quiver;
mod runtime;

use algebra::{
    PyAddMCoresolution, PyAddMResolution, PyAdrAlgebraOfAlgebraResult,
    PyAlgebraAsQuiverAlgebraResult, PyAlgebraElement, PyAlgebraHomomorphism,
    PyAllComplementsResult, PyAllSubmodulesOfModuleResult, PyCommonDirectSummandResult,
    PyCotiltingModuleResult, PyEndOfBasicModuleAsQuiverAlgebraResult,
    PyEndOfModuleAsQuiverAlgebraResult, PyExtAlgebraGeneratorsResult, PyExtOverAlgebraResult,
    PyFiniteAddMResolutionResult, PyFpPathAlgebraByImagesHomomorphism,
    PyFpPathAlgebraMatrixHomomorphism, PyFpPathAlgebraVector, PyGroebnerBasis,
    PyGroebnerBasisIterator, PyGszProjectiveResolution, PyInjectiveResolution,
    PyLiftedProjectivePathAlgebraModule, PyMatrixModule, PyMaximalCommonDirectSummand,
    PyModuleBasisMapping, PyModuleChainMap, PyModuleComplex,
    PyModuleDecompositionWithMultiplicities, PyModuleElement, PyModuleEndomorphismAlgebra,
    PyModuleEndomorphismQuotientAlgebra, PyModuleEndomorphismQuotientElement,
    PyModuleEndomorphismQuotientMap, PyModuleHomomorphism, PyModuleHomomorphismBasis,
    PyModuleSubmodule, PyMoreMinimalVersionResult, PyNaturalSubmoduleHomomorphism,
    PyNonisomorphicDirectSummandData, PyOppositeStructureConstantsAlgebra,
    PyOppositeStructureConstantsElement, PyPathAlgebra, PyPathAlgebraIdeal,
    PyPathAlgebraModuleGroebnerBasis, PyPathAlgebraVector, PyPathAlgebraVectorBasis,
    PyPathAlgebraVectorSpace, PyPathAlgebraVectorXSet, PyPredecessorValuation,
    PyPredecessorsOfModuleResult, PyProjectiveIsomorphismComponent,
    PyProjectivePathAlgebraPresentation, PyProjectiveResolution, PyQuiverAlgebra,
    PyQuotientAlgebra, PyQuotientElement, PyQuotientPathAlgebraVector, PyRightFpPathAlgebraModule,
    PyRightFpPathAlgebraSubmodule, PyRightProjectiveModule, PyRightProjectiveQuotientModule,
    PyRightProjectiveSubmodule, PyStructureConstantsAlgebra, PyStructureConstantsElement,
    PyStructureConstantsHomomorphism, PyStructureConstantsQuotient, PyStructureConstantsSubalgebra,
    PyTiltingModuleResult,
};
pub(crate) use field::{FieldKind, PyField, PyScalar, ScalarKind};
use quiver::{
    PyPath, PyQuiver, PyQuiverArrow, PyQuiverOrdering, PyQuiverPathEnumerator,
    PyQuiverPathIterator, PyQuiverVertex, py_unique_quiver_name,
};

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
    module.add_class::<PyQuiverVertex>()?;
    module.add_class::<PyQuiverArrow>()?;
    module.add_class::<PyQuiverPathIterator>()?;
    module.add_class::<PyQuiverPathEnumerator>()?;
    module.add_class::<PyQuiverOrdering>()?;
    module.add_class::<PyPath>()?;
    module.add_function(wrap_pyfunction!(py_unique_quiver_name, module)?)?;
    module.add_class::<PyPathAlgebra>()?;
    module.add_class::<PyAlgebraElement>()?;
    module.add_class::<PyAlgebraHomomorphism>()?;
    module.add_class::<PyStructureConstantsAlgebra>()?;
    module.add_class::<PyStructureConstantsElement>()?;
    module.add_class::<PyStructureConstantsSubalgebra>()?;
    module.add_class::<PyStructureConstantsQuotient>()?;
    module.add_class::<PyStructureConstantsHomomorphism>()?;
    module.add_class::<PyOppositeStructureConstantsAlgebra>()?;
    module.add_class::<PyOppositeStructureConstantsElement>()?;
    module.add_class::<PyQuiverAlgebra>()?;
    module.add_class::<PyAlgebraAsQuiverAlgebraResult>()?;
    module.add_class::<PyAdrAlgebraOfAlgebraResult>()?;
    module.add_class::<PyPathAlgebraIdeal>()?;
    module.add_class::<PyGroebnerBasis>()?;
    module.add_class::<PyGroebnerBasisIterator>()?;
    module.add_class::<PyRightProjectiveModule>()?;
    module.add_class::<PyRightProjectiveSubmodule>()?;
    module.add_class::<PyPathAlgebraVector>()?;
    module.add_class::<PyPathAlgebraVectorXSet>()?;
    module.add_class::<PyPathAlgebraVectorSpace>()?;
    module.add_class::<PyPathAlgebraVectorBasis>()?;
    module.add_class::<PyPathAlgebraModuleGroebnerBasis>()?;
    module.add_class::<PyRightProjectiveQuotientModule>()?;
    module.add_class::<PyQuotientPathAlgebraVector>()?;
    module.add_class::<PyLiftedProjectivePathAlgebraModule>()?;
    module.add_class::<PyProjectivePathAlgebraPresentation>()?;
    module.add_class::<PyGszProjectiveResolution>()?;
    module.add_class::<PyRightFpPathAlgebraModule>()?;
    module.add_class::<PyFpPathAlgebraVector>()?;
    module.add_class::<PyRightFpPathAlgebraSubmodule>()?;
    module.add_class::<PyNaturalSubmoduleHomomorphism>()?;
    module.add_class::<PyFpPathAlgebraByImagesHomomorphism>()?;
    module.add_class::<PyFpPathAlgebraMatrixHomomorphism>()?;
    module.add_class::<PyQuotientAlgebra>()?;
    module.add_class::<PyQuotientElement>()?;
    module.add_class::<PyMatrixModule>()?;
    module.add_class::<PyModuleElement>()?;
    module.add_class::<PyModuleHomomorphism>()?;
    module.add_class::<PyModuleBasisMapping>()?;
    module.add_class::<PyModuleHomomorphismBasis>()?;
    module.add_class::<PyModuleEndomorphismAlgebra>()?;
    module.add_class::<PyEndOfModuleAsQuiverAlgebraResult>()?;
    module.add_class::<PyEndOfBasicModuleAsQuiverAlgebraResult>()?;
    module.add_class::<PyModuleEndomorphismQuotientAlgebra>()?;
    module.add_class::<PyModuleEndomorphismQuotientElement>()?;
    module.add_class::<PyModuleEndomorphismQuotientMap>()?;
    module.add_class::<PyModuleSubmodule>()?;
    module.add_class::<PyProjectiveResolution>()?;
    module.add_class::<PyInjectiveResolution>()?;
    module.add_class::<PyModuleComplex>()?;
    module.add_class::<PyModuleChainMap>()?;
    module.add_class::<PyProjectiveIsomorphismComponent>()?;
    module.add_class::<PyExtOverAlgebraResult>()?;
    module.add_class::<PyExtAlgebraGeneratorsResult>()?;
    module.add_class::<PyAddMCoresolution>()?;
    module.add_class::<PyFiniteAddMResolutionResult>()?;
    module.add_class::<PyTiltingModuleResult>()?;
    module.add_class::<PyCotiltingModuleResult>()?;
    module.add_class::<PyAddMResolution>()?;
    module.add_class::<PyAllComplementsResult>()?;
    module.add_class::<PyPredecessorValuation>()?;
    module.add_class::<PyPredecessorsOfModuleResult>()?;
    module.add_class::<PyModuleDecompositionWithMultiplicities>()?;
    module.add_class::<PyAllSubmodulesOfModuleResult>()?;
    module.add_class::<PyNonisomorphicDirectSummandData>()?;
    module.add_class::<PyMoreMinimalVersionResult>()?;
    module.add_class::<PyCommonDirectSummandResult>()?;
    module.add_class::<PyMaximalCommonDirectSummand>()?;
    Ok(())
}
