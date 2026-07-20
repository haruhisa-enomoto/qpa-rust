//! `PyO3` bindings over the source-faithful QPA Rust surface.

#![allow(clippy::needless_pass_by_value)]
#![allow(missing_docs)] // PyO3 exception structs cannot receive rustdoc attributes.

mod algebra;
mod dynamic_tree;
mod field;
mod poset;
mod quiver;
mod runtime;
mod static_dictionary;
mod suffix_tree;
mod util;

use algebra::{
    PyAddMCoresolution, PyAddMResolution, PyAdrAlgebraOfAlgebraResult,
    PyAlgebraAsQuiverAlgebraResult, PyAlgebraElement, PyAlgebraHomomorphism,
    PyAllComplementsResult, PyAllSubmodulesOfModuleResult, PyArQuiverNumerical,
    PyBrauerConfiguration, PyBrauerConfigurationAlgebra, PyBrauerConfigurationAlgebraResult,
    PyCodimensionConjectureReport, PyCombinatorialMap, PyCommonDirectSummandResult,
    PyCotiltingModuleResult, PyEndOfBasicModuleAsQuiverAlgebraResult,
    PyEndOfModuleAsQuiverAlgebraResult, PyEulerBilinearForm, PyExtAlgebraGeneratorsResult,
    PyExtOverAlgebraResult, PyFiniteAddMResolutionResult, PyFpPathAlgebraByImagesHomomorphism,
    PyFpPathAlgebraMatrixHomomorphism, PyFpPathAlgebraVector, PyFrobeniusForm,
    PyFrobeniusLinearFunctional, PyGroebnerBasis, PyGroebnerBasisIterator,
    PyGszProjectiveResolution, PyInjectiveResolution, PyLiftedProjectivePathAlgebraModule,
    PyMatrixModule, PyMaximalCommonDirectSummand, PyModuleBasisMapping, PyModuleCategory,
    PyModuleChainMap, PyModuleComplex, PyModuleDecompositionWithMultiplicities, PyModuleElement,
    PyModuleEndomorphismAlgebra, PyModuleEndomorphismQuotientAlgebra,
    PyModuleEndomorphismQuotientElement, PyModuleEndomorphismQuotientMap, PyModuleHomomorphism,
    PyModuleHomomorphismBasis, PyModuleMorphismList, PyModuleSubmodule, PyMoreMinimalVersionResult,
    PyNakayamaAlgebraResult, PyNakayamaAutomorphism, PyNakayamaPermutation, PyNaturalQuotientMap,
    PyNaturalSubmoduleHomomorphism, PyNonisomorphicDirectSummandData,
    PyOppositeStructureConstantsAlgebra, PyOppositeStructureConstantsElement, PyPathAlgebra,
    PyPathAlgebraCentre, PyPathAlgebraIdeal, PyPathAlgebraModuleGroebnerBasis, PyPathAlgebraVector,
    PyPathAlgebraVectorBasis, PyPathAlgebraVectorSpace, PyPathAlgebraVectorXSet,
    PyPathVectorCoefficientsResult, PyPosetAlgebra, PyPredecessorValuation,
    PyPredecessorsOfModuleResult, PyPredefinedAlgebra, PyPredefinedArQuiverData,
    PyProjectiveIsomorphismComponent, PyProjectivePathAlgebraPresentation,
    PyProjectivePresentation, PyProjectiveResolution, PyQuiverAlgebra, PyQuiverAlgebraElement,
    PyQuiverAlgebraUnitGroup, PyQuotientAlgebra, PyQuotientElement, PyQuotientNiceBasis,
    PyQuotientPathAlgebraCentre, PyQuotientPathAlgebraVector, PyQuotientVectorCoefficientsResult,
    PyRightFpPathAlgebraModule, PyRightFpPathAlgebraSubmodule, PyRightProjectiveModule,
    PyRightProjectiveQuotientModule, PyRightProjectiveSubmodule,
    PySimpleEndoProjectiveResolutionResult, PyStringSigmaEpsilon, PyStructureConstantsAlgebra,
    PyStructureConstantsElement, PyStructureConstantsHomomorphism, PyStructureConstantsQuotient,
    PyStructureConstantsSubalgebra, PyTensorProductAlgebra, PyTiltingModuleResult, PyUnitForm,
    py_admissible_sequence_generator, py_brauer_configuration_algebra, py_canonical_algebra,
    py_complex_and_chain_maps, py_compute_lps_array, py_cyclic_permutations_of_band,
    py_finite_global_dim_green, py_finite_global_dim_kirk_kuz, py_format_multiplicity_vector,
    py_format_multiplicity_vectors, py_kronecker_algebra, py_leftmost_occurrence,
    py_leftmost_occurrence_inv_len_inv_lex, py_nakayama_algebra, py_opposite_relations,
    py_poset_algebra, py_predefined_ar_quiver_catalog, py_predefined_ar_quiver_data,
    py_read_algebra, py_string_quiver_rho,
};
use dynamic_tree::{PyDynamicForest, PyDynamicTree};
pub(crate) use field::{FieldKind, PyField, PyScalar, ScalarKind};
use quiver::{
    PyPath, PyQuiver, PyQuiverArrow, PyQuiverOrdering, PyQuiverPathEnumerator,
    PyQuiverPathIterator, PyQuiverProduct, PyQuiverProductDecomposition, PyQuiverVertex,
    py_unique_quiver_name,
};
use static_dictionary::{
    PyQuiverStaticDictionary, PyStaticDictionary, py_static_dictionary_info_messages,
};
use suffix_tree::{
    PySuffixTree, PySuffixTreeEdge, PySuffixTreeNode, py_info_dynamic_dictionary_name,
};

use poset::{PyPoset, PyPosetOrder};
use pyo3::{
    PyErr, PyResult, create_exception, exceptions::PyException, prelude::*, types::PyModule,
};
use util::py_positions_nonzero;

// Base exception for failures reported by the Rust QPA binding.
create_exception!(_qpa, QpaError, PyException);
// Invalid mathematical input, field ownership, or object-family exception.
create_exception!(_qpa, InvalidInputError, QpaError);
// QPA's explicit `fail` result for actions that do not satisfy relations.
create_exception!(_qpa, RelationNotSatisfiedError, QpaError);

fn invalid(error: impl std::fmt::Display) -> PyErr {
    InvalidInputError::new_err(error.to_string())
}

fn register_foundation(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDynamicForest>()?;
    module.add_class::<PyDynamicTree>()?;
    module.add_class::<PySuffixTree>()?;
    module.add_class::<PySuffixTreeNode>()?;
    module.add_class::<PySuffixTreeEdge>()?;
    module.add_class::<PyField>()?;
    module.add_class::<PyScalar>()?;
    module.add_class::<PyQuiver>()?;
    module.add_class::<PyQuiverVertex>()?;
    module.add_class::<PyQuiverArrow>()?;
    module.add_class::<PyQuiverPathIterator>()?;
    module.add_class::<PyQuiverPathEnumerator>()?;
    module.add_class::<PyQuiverOrdering>()?;
    module.add_class::<PyQuiverProduct>()?;
    module.add_class::<PyQuiverProductDecomposition>()?;
    module.add_class::<PyPath>()?;
    module.add_class::<PyPoset>()?;
    module.add_class::<PyPosetOrder>()?;
    module.add_function(wrap_pyfunction!(py_unique_quiver_name, module)?)?;
    module.add_function(wrap_pyfunction!(py_complex_and_chain_maps, module)?)?;
    module.add_function(wrap_pyfunction!(py_compute_lps_array, module)?)?;
    module.add_function(wrap_pyfunction!(py_cyclic_permutations_of_band, module)?)?;
    module.add_function(wrap_pyfunction!(py_predefined_ar_quiver_data, module)?)?;
    module.add_function(wrap_pyfunction!(py_predefined_ar_quiver_catalog, module)?)?;
    module.add_function(wrap_pyfunction!(py_format_multiplicity_vector, module)?)?;
    module.add_function(wrap_pyfunction!(py_format_multiplicity_vectors, module)?)?;
    module.add_function(wrap_pyfunction!(py_string_quiver_rho, module)?)?;
    module.add_function(wrap_pyfunction!(py_admissible_sequence_generator, module)?)?;
    module.add_function(wrap_pyfunction!(py_nakayama_algebra, module)?)?;
    module.add_function(wrap_pyfunction!(py_canonical_algebra, module)?)?;
    module.add_function(wrap_pyfunction!(py_kronecker_algebra, module)?)?;
    module.add_function(wrap_pyfunction!(py_finite_global_dim_green, module)?)?;
    module.add_function(wrap_pyfunction!(py_finite_global_dim_kirk_kuz, module)?)?;
    module.add_function(wrap_pyfunction!(py_poset_algebra, module)?)?;
    module.add_function(wrap_pyfunction!(py_brauer_configuration_algebra, module)?)?;
    module.add_function(wrap_pyfunction!(py_opposite_relations, module)?)?;
    module.add_function(wrap_pyfunction!(py_read_algebra, module)?)?;
    module.add_function(wrap_pyfunction!(py_leftmost_occurrence, module)?)?;
    module.add_function(wrap_pyfunction!(
        py_leftmost_occurrence_inv_len_inv_lex,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(py_positions_nonzero, module)?)?;
    module.add_function(wrap_pyfunction!(py_info_dynamic_dictionary_name, module)?)?;
    module.add_function(wrap_pyfunction!(
        py_static_dictionary_info_messages,
        module
    )?)?;
    Ok(())
}

fn register_algebra(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyStaticDictionary>()?;
    module.add_class::<PyQuiverStaticDictionary>()?;
    module.add_class::<PyPathAlgebra>()?;
    module.add_class::<PyProjectivePresentation>()?;
    module.add_class::<PyAlgebraElement>()?;
    module.add_class::<PyAlgebraHomomorphism>()?;
    module.add_class::<PyNaturalQuotientMap>()?;
    module.add_class::<PyPathAlgebraCentre>()?;
    module.add_class::<PyQuotientPathAlgebraCentre>()?;
    module.add_class::<PyQuiverAlgebraUnitGroup>()?;
    module.add_class::<PyQuotientNiceBasis>()?;
    module.add_class::<PyPathVectorCoefficientsResult>()?;
    module.add_class::<PyQuotientVectorCoefficientsResult>()?;
    module.add_class::<PyStructureConstantsAlgebra>()?;
    module.add_class::<PyStructureConstantsElement>()?;
    module.add_class::<PyStructureConstantsSubalgebra>()?;
    module.add_class::<PyStructureConstantsQuotient>()?;
    module.add_class::<PyStructureConstantsHomomorphism>()?;
    module.add_class::<PyOppositeStructureConstantsAlgebra>()?;
    module.add_class::<PyOppositeStructureConstantsElement>()?;
    module.add_class::<PyQuiverAlgebra>()?;
    module.add_class::<PyQuiverAlgebraElement>()?;
    module.add_class::<PyNakayamaPermutation>()?;
    module.add_class::<PyNakayamaAutomorphism>()?;
    module.add_class::<PyFrobeniusForm>()?;
    module.add_class::<PyFrobeniusLinearFunctional>()?;
    module.add_class::<PyTensorProductAlgebra>()?;
    module.add_class::<PyCombinatorialMap>()?;
    module.add_class::<PyUnitForm>()?;
    module.add_class::<PyEulerBilinearForm>()?;
    module.add_class::<PyArQuiverNumerical>()?;
    module.add_class::<PyPredefinedArQuiverData>()?;
    module.add_class::<PyCodimensionConjectureReport>()?;
    module.add_class::<PyStringSigmaEpsilon>()?;
    module.add_class::<PyPredefinedAlgebra>()?;
    module.add_class::<PyNakayamaAlgebraResult>()?;
    module.add_class::<PyPosetAlgebra>()?;
    module.add_class::<PyBrauerConfiguration>()?;
    module.add_class::<PyBrauerConfigurationAlgebra>()?;
    module.add_class::<PyBrauerConfigurationAlgebraResult>()?;
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
    Ok(())
}

fn register_modules(module: &Bound<'_, PyModule>) -> PyResult<()> {
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
    module.add_class::<PySimpleEndoProjectiveResolutionResult>()?;
    module.add_class::<PyModuleComplex>()?;
    module.add_class::<PyModuleCategory>()?;
    module.add_class::<PyModuleMorphismList>()?;
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
    register_foundation(module)?;
    register_algebra(module)?;
    register_modules(module)
}
