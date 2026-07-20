//! Python wrappers for path algebras, ideals, Gröbner bases, and quotients.

mod associated;
mod brauer;
mod combinatorial;
mod combinatorial_rep;
mod complex_foundation;
mod conversion;
mod corner;
mod degeneration;
mod element_utilities;
mod finite_structures;
mod fp_projective;
mod groebner;
mod homological;
mod ideal;
mod invariants;
mod mapping;
mod module;
mod nakayama;
mod predefined;
mod projective;
mod projective_gsz;
mod quotient;
mod quotient_projective;
mod string;
mod tensor;
mod validity;
mod vector_reduction;

pub(crate) use brauer::{
    PyBrauerConfiguration, PyBrauerConfigurationAlgebra, PyBrauerConfigurationAlgebraResult,
    py_brauer_configuration_algebra,
};
pub(crate) use combinatorial::PyCombinatorialMap;
pub(crate) use combinatorial_rep::{PyEulerBilinearForm, PyUnitForm};
pub(crate) use complex_foundation::{PyModuleCategory, PyModuleMorphismList};
use conversion::{
    AdrAlgebraOfAlgebraResultKind, StructureConstantsAlgebraKind, adr_algebra_of_basic_module_for,
    wrap_quiver_algebra,
};
pub(crate) use conversion::{
    PyAdrAlgebraOfAlgebraResult, PyAlgebraAsQuiverAlgebraResult,
    PyOppositeStructureConstantsAlgebra, PyOppositeStructureConstantsElement, PyQuiverAlgebra,
    PyStructureConstantsAlgebra, PyStructureConstantsElement, PyStructureConstantsHomomorphism,
    PyStructureConstantsQuotient, PyStructureConstantsSubalgebra,
};
pub(crate) use degeneration::{
    PyArQuiverNumerical, PyCodimensionConjectureReport, PyPredefinedArQuiverData,
    py_format_multiplicity_vector, py_format_multiplicity_vectors, py_predefined_ar_quiver_catalog,
    py_predefined_ar_quiver_data,
};
pub(crate) use element_utilities::py_opposite_relations;
pub(crate) use finite_structures::{
    PyPathAlgebraCentre, PyQuiverAlgebraUnitGroup, PyQuotientPathAlgebraCentre,
};
pub(crate) use fp_projective::{
    PyFpPathAlgebraByImagesHomomorphism, PyFpPathAlgebraMatrixHomomorphism, PyFpPathAlgebraVector,
    PyNaturalSubmoduleHomomorphism, PyRightFpPathAlgebraModule, PyRightFpPathAlgebraSubmodule,
};
pub(crate) use homological::py_complex_and_chain_maps;
pub(crate) use mapping::PyNaturalQuotientMap;
pub(crate) use module::{PyGszProjectiveResolution, PyProjectivePathAlgebraPresentation};
pub(crate) use nakayama::{
    PyFrobeniusForm, PyFrobeniusLinearFunctional, PyNakayamaAutomorphism, PyNakayamaPermutation,
};
pub(crate) use predefined::{
    PyNakayamaAlgebraResult, PyPosetAlgebra, PyPredefinedAlgebra, py_admissible_sequence_generator,
    py_canonical_algebra, py_finite_global_dim_green, py_finite_global_dim_kirk_kuz,
    py_kronecker_algebra, py_nakayama_algebra, py_poset_algebra,
};
pub(crate) use projective::{
    PyPathAlgebraModuleGroebnerBasis, PyPathAlgebraVector, PyPathAlgebraVectorBasis,
    PyPathAlgebraVectorSpace, PyRightProjectiveModule, PyRightProjectiveSubmodule,
};
pub(crate) use projective_gsz::PyPathAlgebraVectorXSet;
pub(crate) use quotient_projective::{
    PyLiftedProjectivePathAlgebraModule, PyQuotientPathAlgebraVector,
    PyRightProjectiveQuotientModule,
};
pub(crate) use string::{py_compute_lps_array, py_cyclic_permutations_of_band};
pub(crate) use tensor::{PyQuiverAlgebraElement, PyTensorProductAlgebra};
pub(crate) use validity::{PyStringSigmaEpsilon, py_string_quiver_rho};
pub(crate) use vector_reduction::{
    PyPathVectorCoefficientsResult, PyQuotientVectorCoefficientsResult,
};

use numerica::domains::{
    finite_field::{Z2, Zp},
    integer::Integer,
    rational::{Q, RationalField},
};
use pyo3::{Py, PyRef, PyResult, Python, prelude::*, types::PyList};
use qpa_core::{
    AlgebraDimension, DenseMatrix, GroebnerBasis, GroebnerBasisKind, NamedArrowMatrix, PathAlgebra,
    PathAlgebraElement, PathAlgebraExternalRepresentation, PathAlgebraIdeal, PathAlgebraModule,
    PathAlgebraModuleElement, PathAlgebraModuleHomomorphism, PathAlgebraModuleSubmodule,
    QuiverAlgebra, QuiverAlgebraElement, QuotientPathAlgebra, QuotientPathAlgebraElement,
    right_module_over_path_algebra,
};

use crate::field::parse_python_integer;
use crate::{
    FieldKind, PyField, PyScalar, RelationNotSatisfiedError, invalid,
    quiver::{PyPath, PyQuiver, PyQuiverOrdering, ordering_value},
    runtime::{self, PythonField as _, dispatch_pair_ref, dispatch_ref, map_pair_ref, map_ref},
};
#[derive(Clone)]
enum AlgebraKind {
    Rational(PathAlgebra<RationalField>),
    Binary(PathAlgebra<Z2>),
    Prime(PathAlgebra<Zp>),
}

/// Runtime-dispatched path algebra backed by the existing generic Rust type.
#[pyclass(name = "PathAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPathAlgebra {
    kind: AlgebraKind,
}

#[derive(Clone)]
enum AlgebraElementKind {
    Rational(PathAlgebraElement<RationalField>),
    Binary(PathAlgebraElement<Z2>),
    Prime(PathAlgebraElement<Zp>),
}

macro_rules! impl_algebra_element_extract {
    ($field:ty, $variant:ident) => {
        impl runtime::Extract<PathAlgebraElement<$field>> for AlgebraElementKind {
            fn extract(self) -> PyResult<PathAlgebraElement<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid("relation belongs to another coefficient field")),
                }
            }
        }

        impl runtime::Wrap<PathAlgebraElement<$field>> for AlgebraElementKind {
            fn wrap(value: PathAlgebraElement<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_algebra_element_extract!(RationalField, Rational);
impl_algebra_element_extract!(Z2, Binary);
impl_algebra_element_extract!(Zp, Prime);

/// Sparse element of a runtime-dispatched Python path algebra.
#[pyclass(name = "AlgebraElement", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAlgebraElement {
    kind: AlgebraElementKind,
}

#[derive(Clone)]
enum IdealKind {
    Rational(PathAlgebraIdeal<RationalField>),
    Binary(PathAlgebraIdeal<Z2>),
    Prime(PathAlgebraIdeal<Zp>),
}

/// Two-sided path-algebra ideal using QPA's lazy Gröbner-basis behavior.
#[pyclass(name = "PathAlgebraIdeal", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPathAlgebraIdeal {
    kind: IdealKind,
}

#[derive(Clone)]
enum GroebnerKind {
    Rational(GroebnerBasis<RationalField>),
    Binary(GroebnerBasis<Z2>),
    Prime(GroebnerBasis<Zp>),
}

/// QPA's two-sided or right Gröbner-basis object.
#[pyclass(name = "GroebnerBasis", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyGroebnerBasis {
    kind: GroebnerKind,
}

/// Mutable iterator over a snapshot of QPA's stored Gröbner relations.
#[pyclass(name = "GroebnerBasisIterator", skip_from_py_object)]
pub(crate) struct PyGroebnerBasisIterator {
    relations: Vec<PyAlgebraElement>,
    position: usize,
}

macro_rules! impl_runtime_wrap {
    ($wrapper:ty, $value:ty, $variant:ident) => {
        impl runtime::Wrap<$value> for $wrapper {
            fn wrap(value: $value) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_runtime_wrap!(IdealKind, PathAlgebraIdeal<RationalField>, Rational);
impl_runtime_wrap!(IdealKind, PathAlgebraIdeal<Z2>, Binary);
impl_runtime_wrap!(IdealKind, PathAlgebraIdeal<Zp>, Prime);
impl_runtime_wrap!(GroebnerKind, GroebnerBasis<RationalField>, Rational);
impl_runtime_wrap!(GroebnerKind, GroebnerBasis<Z2>, Binary);
impl_runtime_wrap!(GroebnerKind, GroebnerBasis<Zp>, Prime);

#[derive(Clone)]
enum QuotientKind {
    Rational(QuotientPathAlgebra<RationalField>),
    Binary(QuotientPathAlgebra<Z2>),
    Prime(QuotientPathAlgebra<Zp>),
}

/// Quotient path algebra constructed by QPA's high-level Groebner route.
#[pyclass(name = "QuotientAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuotientAlgebra {
    kind: QuotientKind,
}

#[derive(Clone)]
enum QuotientElementKind {
    Rational(QuotientPathAlgebraElement<RationalField>),
    Binary(QuotientPathAlgebraElement<Z2>),
    Prime(QuotientPathAlgebraElement<Zp>),
}

macro_rules! impl_quotient_element_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Extract<QuotientPathAlgebraElement<$field>> for QuotientElementKind {
            fn extract(self) -> PyResult<QuotientPathAlgebraElement<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "quotient element belongs to another coefficient field",
                    )),
                }
            }
        }

        impl runtime::Wrap<QuotientPathAlgebraElement<$field>> for QuotientElementKind {
            fn wrap(value: QuotientPathAlgebraElement<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_quotient_element_runtime!(RationalField, Rational);
impl_quotient_element_runtime!(Z2, Binary);
impl_quotient_element_runtime!(Zp, Prime);

/// Element of a runtime-dispatched quotient path algebra.
#[pyclass(name = "QuotientElement", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyQuotientElement {
    kind: QuotientElementKind,
}

#[derive(Clone)]
enum AlgebraHomomorphismKind {
    Rational(qpa_core::QuiverAlgebraHomomorphism<RationalField>),
    Binary(qpa_core::QuiverAlgebraHomomorphism<Z2>),
    Prime(qpa_core::QuiverAlgebraHomomorphism<Zp>),
}

/// A QPA quiver-algebra homomorphism stored by ordered generator images.
#[pyclass(name = "AlgebraHomomorphism", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAlgebraHomomorphism {
    kind: AlgebraHomomorphismKind,
}

macro_rules! impl_algebra_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<PathAlgebra<$field>> for AlgebraKind {
            fn wrap(value: PathAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<QuotientPathAlgebra<$field>> for QuotientKind {
            fn wrap(value: QuotientPathAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::QuiverAlgebraHomomorphism<$field>>
            for AlgebraHomomorphismKind
        {
            fn wrap(value: qpa_core::QuiverAlgebraHomomorphism<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Extract<qpa_core::QuiverAlgebraHomomorphism<$field>>
            for AlgebraHomomorphismKind
        {
            fn extract(self) -> PyResult<qpa_core::QuiverAlgebraHomomorphism<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "algebra homomorphism belongs to another coefficient field",
                    )),
                }
            }
        }
    };
}

impl_algebra_runtime!(RationalField, Rational);
impl_algebra_runtime!(Z2, Binary);
impl_algebra_runtime!(Zp, Prime);

#[derive(Clone)]
enum MatrixModuleKind {
    Rational(PathAlgebraModule<RationalField>),
    Binary(PathAlgebraModule<Z2>),
    Prime(PathAlgebraModule<Zp>),
}

macro_rules! impl_module_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Extract<PathAlgebraModule<$field>> for MatrixModuleKind {
            fn extract(self) -> PyResult<PathAlgebraModule<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid("module belongs to another coefficient field")),
                }
            }
        }

        impl runtime::Wrap<PathAlgebraModule<$field>> for MatrixModuleKind {
            fn wrap(value: PathAlgebraModule<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_module_runtime!(RationalField, Rational);
impl_module_runtime!(Z2, Binary);
impl_module_runtime!(Zp, Prime);

/// Checked QPA right matrix module over a path algebra or quotient.
#[pyclass(name = "MatrixModule", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyMatrixModule {
    kind: MatrixModuleKind,
}

#[derive(Clone)]
enum ModuleElementKind {
    Rational(PathAlgebraModuleElement<RationalField>),
    Binary(PathAlgebraModuleElement<Z2>),
    Prime(PathAlgebraModuleElement<Zp>),
}

/// Element of a QPA matrix module over a path algebra or quotient.
#[pyclass(name = "ModuleElement", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleElement {
    kind: ModuleElementKind,
}

#[derive(Clone)]
enum ModuleHomomorphismKind {
    Rational(PathAlgebraModuleHomomorphism<RationalField>),
    Binary(PathAlgebraModuleHomomorphism<Z2>),
    Prime(PathAlgebraModuleHomomorphism<Zp>),
}

/// QPA homomorphism between right matrix modules.
#[pyclass(name = "ModuleHomomorphism", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleHomomorphism {
    kind: ModuleHomomorphismKind,
}

#[derive(Clone)]
enum ModuleBasisMappingKind {
    Rational(qpa_core::PathAlgebraModuleBasisMapping<RationalField>),
    Binary(qpa_core::PathAlgebraModuleBasisMapping<Z2>),
    Prime(qpa_core::PathAlgebraModuleBasisMapping<Zp>),
}

/// QPA's basis-image mapping returned for a quotient by a subalgebra module.
#[pyclass(name = "ModuleBasisMapping", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleBasisMapping {
    kind: ModuleBasisMappingKind,
}

#[derive(Clone)]
enum ModuleHomomorphismBasisKind {
    Rational(qpa_core::PathAlgebraModuleHomomorphismBasis<RationalField>),
    Binary(qpa_core::PathAlgebraModuleHomomorphismBasis<Z2>),
    Prime(qpa_core::PathAlgebraModuleHomomorphismBasis<Zp>),
}

/// QPA's ordered Hom basis together with its coefficient function.
#[pyclass(name = "ModuleHomomorphismBasis", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleHomomorphismBasis {
    kind: ModuleHomomorphismBasisKind,
}

#[derive(Clone)]
enum ModuleEndomorphismAlgebraKind {
    Rational(qpa_core::ModuleEndomorphismAlgebra<RationalField>),
    Binary(qpa_core::ModuleEndomorphismAlgebra<Z2>),
    Prime(qpa_core::ModuleEndomorphismAlgebra<Zp>),
}

/// QPA's matrix-algebra realization of a module endomorphism algebra.
#[pyclass(name = "ModuleEndomorphismAlgebra", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleEndomorphismAlgebra {
    kind: ModuleEndomorphismAlgebraKind,
}

#[derive(Clone)]
enum ModuleEndomorphismQuotientAlgebraKind {
    Rational(qpa_core::ModuleEndomorphismQuotientAlgebra<RationalField>),
    Binary(qpa_core::ModuleEndomorphismQuotientAlgebra<Z2>),
    Prime(qpa_core::ModuleEndomorphismQuotientAlgebra<Zp>),
}

/// QPA's stable endomorphism algebra in natural-quotient coordinates.
#[pyclass(
    name = "ModuleEndomorphismQuotientAlgebra",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PyModuleEndomorphismQuotientAlgebra {
    kind: ModuleEndomorphismQuotientAlgebraKind,
}

#[derive(Clone)]
enum ModuleEndomorphismQuotientElementKind {
    Rational(qpa_core::ModuleEndomorphismQuotientElement<RationalField>),
    Binary(qpa_core::ModuleEndomorphismQuotientElement<Z2>),
    Prime(qpa_core::ModuleEndomorphismQuotientElement<Zp>),
}

/// One exact element of a stable endomorphism quotient algebra.
#[pyclass(
    name = "ModuleEndomorphismQuotientElement",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PyModuleEndomorphismQuotientElement {
    kind: ModuleEndomorphismQuotientElementKind,
}

#[derive(Clone)]
enum ModuleEndomorphismQuotientMapKind {
    Rational(qpa_core::ModuleEndomorphismQuotientHomomorphism<RationalField>),
    Binary(qpa_core::ModuleEndomorphismQuotientHomomorphism<Z2>),
    Prime(qpa_core::ModuleEndomorphismQuotientHomomorphism<Zp>),
}

/// QPA's natural map from `End(M)` to its quotient modulo projectives.
#[pyclass(name = "ModuleEndomorphismQuotientMap", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleEndomorphismQuotientMap {
    kind: ModuleEndomorphismQuotientMapKind,
}

#[derive(Clone)]
enum ModuleSubmoduleKind {
    Rational(PathAlgebraModuleSubmodule<RationalField>),
    Binary(PathAlgebraModuleSubmodule<Z2>),
    Prime(PathAlgebraModuleSubmodule<Zp>),
}

/// QPA's generator-defined concrete submodule of a path-algebra matrix module.
#[pyclass(name = "ModuleSubmodule", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleSubmodule {
    kind: ModuleSubmoduleKind,
}

#[derive(Clone)]
enum ProjectiveResolutionKind {
    Rational(qpa_gbnp::ProjectiveResolution<RationalField>),
    Binary(qpa_gbnp::ProjectiveResolution<Z2>),
    Prime(qpa_gbnp::ProjectiveResolution<Zp>),
}

/// QPA's lazy projective resolution, including the resolved module in degree -1.
#[pyclass(name = "ProjectiveResolution", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyProjectiveResolution {
    kind: ProjectiveResolutionKind,
}

#[derive(Clone)]
enum InjectiveResolutionKind {
    Rational(qpa_gbnp::InjectiveResolution<RationalField>),
    Binary(qpa_gbnp::InjectiveResolution<Z2>),
    Prime(qpa_gbnp::InjectiveResolution<Zp>),
}

/// QPA's lazy injective resolution, including the resolved module in degree zero.
#[pyclass(name = "InjectiveResolution", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyInjectiveResolution {
    kind: InjectiveResolutionKind,
}

#[derive(Clone)]
enum SimpleEndoProjectiveResolutionKind {
    Rational(qpa_gbnp::ProjectiveResolutionOfSimpleModuleOverEndoResult<RationalField>),
    Binary(qpa_gbnp::ProjectiveResolutionOfSimpleModuleOverEndoResult<Z2>),
    Prime(qpa_gbnp::ProjectiveResolutionOfSimpleModuleOverEndoResult<Zp>),
}

/// QPA's projective-dimension status and stripped syzygies over an endomorphism ring.
#[pyclass(
    name = "ProjectiveResolutionOfSimpleModuleOverEndoResult",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PySimpleEndoProjectiveResolutionResult {
    kind: SimpleEndoProjectiveResolutionKind,
}

#[derive(Clone)]
enum ModuleComplexKind {
    Rational(qpa_gbnp::ModuleComplex<RationalField>),
    Binary(qpa_gbnp::ModuleComplex<Z2>),
    Prime(qpa_gbnp::ModuleComplex<Zp>),
}

/// A QPA complex in the category of right modules over one algebra.
#[pyclass(name = "ModuleComplex", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleComplex {
    kind: ModuleComplexKind,
}

#[derive(Clone)]
enum ModuleChainMapKind {
    Rational(qpa_gbnp::ModuleChainMap<RationalField>),
    Binary(qpa_gbnp::ModuleChainMap<Z2>),
    Prime(qpa_gbnp::ModuleChainMap<Zp>),
}

/// A QPA chain map between complexes of right modules.
#[pyclass(name = "ModuleChainMap", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyModuleChainMap {
    kind: ModuleChainMapKind,
}

#[derive(Clone)]
enum ProjectiveIsomorphismComponentKind {
    Rational(qpa_gbnp::ProjectiveIsomorphismComponent<RationalField>),
    Binary(qpa_gbnp::ProjectiveIsomorphismComponent<Z2>),
    Prime(qpa_gbnp::ProjectiveIsomorphismComponent<Zp>),
}

/// QPA's selected isomorphism component in a map of decomposed projectives.
#[pyclass(name = "ProjectiveIsomorphismComponent", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyProjectiveIsomorphismComponent {
    kind: ProjectiveIsomorphismComponentKind,
}

#[derive(Clone)]
enum ExtOverAlgebraKind {
    Rational(qpa_gbnp::ExtOverAlgebraResult<RationalField>),
    Binary(qpa_gbnp::ExtOverAlgebraResult<Z2>),
    Prime(qpa_gbnp::ExtOverAlgebraResult<Zp>),
}

/// QPA's syzygy inclusion, ordered basis, and coefficient closure for Ext^1.
#[pyclass(name = "ExtOverAlgebraResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyExtOverAlgebraResult {
    kind: ExtOverAlgebraKind,
}

#[derive(Clone)]
enum ExtAlgebraGeneratorsKind {
    Rational(qpa_gbnp::ExtAlgebraGeneratorsResult<RationalField>),
    Binary(qpa_gbnp::ExtAlgebraGeneratorsResult<Z2>),
    Prime(qpa_gbnp::ExtAlgebraGeneratorsResult<Zp>),
}

/// QPA's degreewise dimensions and retained Ext-algebra generators.
#[pyclass(name = "ExtAlgebraGeneratorsResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyExtAlgebraGeneratorsResult {
    kind: ExtAlgebraGeneratorsKind,
}

#[derive(Clone)]
enum AddMCoresolutionKind {
    Rational(qpa_gbnp::AddMCoresolution<RationalField>),
    Binary(qpa_gbnp::AddMCoresolution<Z2>),
    Prime(qpa_gbnp::AddMCoresolution<Zp>),
}

/// A finite add(M)-coresolution returned by QPA.
#[pyclass(name = "AddMCoresolution", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAddMCoresolution {
    kind: AddMCoresolutionKind,
}

#[derive(Clone)]
enum FiniteAddMResolutionKind {
    Rational(qpa_gbnp::FiniteAddMResolutionResult<RationalField>),
    Binary(qpa_gbnp::FiniteAddMResolutionResult<Z2>),
    Prime(qpa_gbnp::FiniteAddMResolutionResult<Zp>),
}

/// QPA's literal-true or finite-complex add(M)-resolution result.
#[pyclass(name = "FiniteAddMResolutionResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyFiniteAddMResolutionResult {
    kind: FiniteAddMResolutionKind,
}

#[derive(Clone)]
enum TiltingModuleResultKind {
    Rational(qpa_gbnp::TiltingModuleResult<RationalField>),
    Binary(qpa_gbnp::TiltingModuleResult<Z2>),
    Prime(qpa_gbnp::TiltingModuleResult<Zp>),
}

/// QPA's bounded tilting-module success record.
#[pyclass(name = "TiltingModuleResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyTiltingModuleResult {
    kind: TiltingModuleResultKind,
}

#[derive(Clone)]
enum CotiltingModuleResultKind {
    Rational(qpa_gbnp::CotiltingModuleResult<RationalField>),
    Binary(qpa_gbnp::CotiltingModuleResult<Z2>),
    Prime(qpa_gbnp::CotiltingModuleResult<Zp>),
}

/// QPA's bounded cotilting-module success record.
#[pyclass(name = "CotiltingModuleResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyCotiltingModuleResult {
    kind: CotiltingModuleResultKind,
}

#[derive(Clone)]
enum AddMResolutionKind {
    Rational(qpa_gbnp::AddMResolution<RationalField>),
    Binary(qpa_gbnp::AddMResolution<Z2>),
    Prime(qpa_gbnp::AddMResolution<Zp>),
}

/// A finite add(M)-resolution returned by QPA's complement computation.
#[pyclass(name = "AddMResolution", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAddMResolution {
    kind: AddMResolutionKind,
}

#[derive(Clone)]
enum AllComplementsResultKind {
    Rational(qpa_gbnp::AllComplementsResult<RationalField>),
    Binary(qpa_gbnp::AllComplementsResult<Z2>),
    Prime(qpa_gbnp::AllComplementsResult<Zp>),
}

/// QPA's possibly empty resolution and coresolution complement pair.
#[pyclass(name = "AllComplementsResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAllComplementsResult {
    kind: AllComplementsResultKind,
}

#[derive(Clone)]
enum PredecessorsOfModuleResultKind {
    Rational(qpa_gbnp::PredecessorsOfModuleResult<RationalField>),
    Binary(qpa_gbnp::PredecessorsOfModuleResult<Z2>),
    Prime(qpa_gbnp::PredecessorsOfModuleResult<Zp>),
}

/// QPA's module-or-layered predecessor result.
#[pyclass(name = "PredecessorsOfModuleResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPredecessorsOfModuleResult {
    kind: PredecessorsOfModuleResultKind,
}

/// One one-based QPA predecessor valuation entry.
#[pyclass(name = "PredecessorValuation", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyPredecessorValuation {
    value: qpa_gbnp::PredecessorValuation,
}

#[derive(Clone)]
enum ModuleDecompositionWithMultiplicitiesKind {
    Rational(qpa_core::ModuleDecompositionWithMultiplicities<RationalField>),
    Binary(qpa_core::ModuleDecompositionWithMultiplicities<Z2>),
    Prime(qpa_core::ModuleDecompositionWithMultiplicities<Zp>),
}

/// QPA's basic summands paired with their multiplicities.
#[pyclass(
    name = "ModuleDecompositionWithMultiplicities",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PyModuleDecompositionWithMultiplicities {
    kind: ModuleDecompositionWithMultiplicitiesKind,
}

#[derive(Clone)]
enum AllSubmodulesOfModuleResultKind {
    Rational(qpa_core::AllSubmodulesOfModuleResult<RationalField>),
    Binary(qpa_core::AllSubmodulesOfModuleResult<Z2>),
    Prime(qpa_core::AllSubmodulesOfModuleResult<Zp>),
}

/// QPA's flat or dimension-indexed list of submodule inclusions.
#[pyclass(name = "AllSubmodulesOfModuleResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAllSubmodulesOfModuleResult {
    kind: AllSubmodulesOfModuleResultKind,
}

/// QPA's number and semisimple-block dimensions of nonisomorphic summands.
#[pyclass(name = "NonisomorphicDirectSummandData", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyNonisomorphicDirectSummandData {
    value: qpa_core::NonisomorphicDirectSummandData,
}

#[derive(Clone)]
enum MoreMinimalVersionResultKind {
    Rational(qpa_core::MoreMinimalVersionResult<RationalField>),
    Binary(qpa_core::MoreMinimalVersionResult<Z2>),
    Prime(qpa_core::MoreMinimalVersionResult<Zp>),
}

/// QPA's unchanged-map or one-step-reduced minimal-version result.
#[pyclass(name = "MoreMinimalVersionResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyMoreMinimalVersionResult {
    kind: MoreMinimalVersionResultKind,
}

#[derive(Clone)]
enum CommonDirectSummandResultKind {
    Rational(qpa_core::CommonDirectSummandResult<RationalField>),
    Binary(qpa_core::CommonDirectSummandResult<Z2>),
    Prime(qpa_core::CommonDirectSummandResult<Zp>),
}

/// QPA's different-algebra, false, or four-module common-summand result.
#[pyclass(name = "CommonDirectSummandResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyCommonDirectSummandResult {
    kind: CommonDirectSummandResultKind,
}

#[derive(Clone)]
enum MaximalCommonDirectSummandKind {
    Rational(qpa_core::MaximalCommonDirectSummand<RationalField>),
    Binary(qpa_core::MaximalCommonDirectSummand<Z2>),
    Prime(qpa_core::MaximalCommonDirectSummand<Zp>),
}

/// QPA's ordered common summands and two remaining complements.
#[pyclass(name = "MaximalCommonDirectSummand", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyMaximalCommonDirectSummand {
    kind: MaximalCommonDirectSummandKind,
}

#[derive(Clone)]
enum EndOfModuleAsQuiverAlgebraResultKind {
    Rational(qpa_gbnp::EndOfModuleAsQuiverAlgebraResult<RationalField>),
    Binary(qpa_gbnp::EndOfModuleAsQuiverAlgebraResult<Z2>),
    Prime(qpa_gbnp::EndOfModuleAsQuiverAlgebraResult<Zp>),
}

/// QPA's endomorphism algebra, adjacency matrix, and quiver-algebra conversion.
#[pyclass(name = "EndOfModuleAsQuiverAlgebraResult", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyEndOfModuleAsQuiverAlgebraResult {
    kind: EndOfModuleAsQuiverAlgebraResultKind,
}

#[derive(Clone)]
enum EndOfBasicModuleAsQuiverAlgebraResultKind {
    Rational(qpa_gbnp::EndOfBasicModuleAsQuiverAlgebraResult<RationalField>),
    Binary(qpa_gbnp::EndOfBasicModuleAsQuiverAlgebraResult<Z2>),
    Prime(qpa_gbnp::EndOfBasicModuleAsQuiverAlgebraResult<Zp>),
}

/// QPA's quiver presentation of the endomorphism ring of a basic module list.
#[pyclass(
    name = "EndOfBasicModuleAsQuiverAlgebraResult",
    frozen,
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PyEndOfBasicModuleAsQuiverAlgebraResult {
    kind: EndOfBasicModuleAsQuiverAlgebraResultKind,
}

macro_rules! impl_module_value_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Extract<PathAlgebraModuleElement<$field>> for ModuleElementKind {
            fn extract(self) -> PyResult<PathAlgebraModuleElement<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "module element belongs to another coefficient field",
                    )),
                }
            }
        }

        impl runtime::Wrap<PathAlgebraModuleElement<$field>> for ModuleElementKind {
            fn wrap(value: PathAlgebraModuleElement<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<PathAlgebraModuleHomomorphism<$field>> for ModuleHomomorphismKind {
            fn wrap(value: PathAlgebraModuleHomomorphism<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Extract<PathAlgebraModuleHomomorphism<$field>> for ModuleHomomorphismKind {
            fn extract(self) -> PyResult<PathAlgebraModuleHomomorphism<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "module homomorphism belongs to another coefficient field",
                    )),
                }
            }
        }

        impl runtime::Wrap<PathAlgebraModuleSubmodule<$field>> for ModuleSubmoduleKind {
            fn wrap(value: PathAlgebraModuleSubmodule<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::PathAlgebraModuleBasisMapping<$field>>
            for ModuleBasisMappingKind
        {
            fn wrap(value: qpa_core::PathAlgebraModuleBasisMapping<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::PathAlgebraModuleHomomorphismBasis<$field>>
            for ModuleHomomorphismBasisKind
        {
            fn wrap(value: qpa_core::PathAlgebraModuleHomomorphismBasis<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::ModuleEndomorphismAlgebra<$field>>
            for ModuleEndomorphismAlgebraKind
        {
            fn wrap(value: qpa_core::ModuleEndomorphismAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::ModuleEndomorphismQuotientAlgebra<$field>>
            for ModuleEndomorphismQuotientAlgebraKind
        {
            fn wrap(value: qpa_core::ModuleEndomorphismQuotientAlgebra<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::ModuleEndomorphismQuotientElement<$field>>
            for ModuleEndomorphismQuotientElementKind
        {
            fn wrap(value: qpa_core::ModuleEndomorphismQuotientElement<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::ModuleEndomorphismQuotientHomomorphism<$field>>
            for ModuleEndomorphismQuotientMapKind
        {
            fn wrap(value: qpa_core::ModuleEndomorphismQuotientHomomorphism<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_module_value_runtime!(RationalField, Rational);
impl_module_value_runtime!(Z2, Binary);
impl_module_value_runtime!(Zp, Prime);

macro_rules! impl_resolution_runtime {
    ($field:ty, $variant:ident) => {
        impl runtime::Wrap<qpa_gbnp::ProjectiveResolution<$field>> for ProjectiveResolutionKind {
            fn wrap(value: qpa_gbnp::ProjectiveResolution<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::InjectiveResolution<$field>> for InjectiveResolutionKind {
            fn wrap(value: qpa_gbnp::InjectiveResolution<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::ProjectiveResolutionOfSimpleModuleOverEndoResult<$field>>
            for SimpleEndoProjectiveResolutionKind
        {
            fn wrap(
                value: qpa_gbnp::ProjectiveResolutionOfSimpleModuleOverEndoResult<$field>,
            ) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::ModuleComplex<$field>> for ModuleComplexKind {
            fn wrap(value: qpa_gbnp::ModuleComplex<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Extract<qpa_gbnp::ModuleComplex<$field>> for ModuleComplexKind {
            fn extract(self) -> PyResult<qpa_gbnp::ModuleComplex<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid("complex belongs to another coefficient field")),
                }
            }
        }

        impl runtime::Wrap<qpa_gbnp::ModuleChainMap<$field>> for ModuleChainMapKind {
            fn wrap(value: qpa_gbnp::ModuleChainMap<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::ProjectiveIsomorphismComponent<$field>>
            for ProjectiveIsomorphismComponentKind
        {
            fn wrap(value: qpa_gbnp::ProjectiveIsomorphismComponent<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Extract<qpa_gbnp::ProjectiveIsomorphismComponent<$field>>
            for ProjectiveIsomorphismComponentKind
        {
            fn extract(self) -> PyResult<qpa_gbnp::ProjectiveIsomorphismComponent<$field>> {
                match self {
                    Self::$variant(value) => Ok(value),
                    _ => Err(invalid(
                        "isomorphism component belongs to another coefficient field",
                    )),
                }
            }
        }

        impl runtime::Wrap<qpa_gbnp::ExtOverAlgebraResult<$field>> for ExtOverAlgebraKind {
            fn wrap(value: qpa_gbnp::ExtOverAlgebraResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::ExtAlgebraGeneratorsResult<$field>>
            for ExtAlgebraGeneratorsKind
        {
            fn wrap(value: qpa_gbnp::ExtAlgebraGeneratorsResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::AddMCoresolution<$field>> for AddMCoresolutionKind {
            fn wrap(value: qpa_gbnp::AddMCoresolution<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::FiniteAddMResolutionResult<$field>>
            for FiniteAddMResolutionKind
        {
            fn wrap(value: qpa_gbnp::FiniteAddMResolutionResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::TiltingModuleResult<$field>> for TiltingModuleResultKind {
            fn wrap(value: qpa_gbnp::TiltingModuleResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::CotiltingModuleResult<$field>> for CotiltingModuleResultKind {
            fn wrap(value: qpa_gbnp::CotiltingModuleResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::AddMResolution<$field>> for AddMResolutionKind {
            fn wrap(value: qpa_gbnp::AddMResolution<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::AllComplementsResult<$field>> for AllComplementsResultKind {
            fn wrap(value: qpa_gbnp::AllComplementsResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::PredecessorsOfModuleResult<$field>>
            for PredecessorsOfModuleResultKind
        {
            fn wrap(value: qpa_gbnp::PredecessorsOfModuleResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::ModuleDecompositionWithMultiplicities<$field>>
            for ModuleDecompositionWithMultiplicitiesKind
        {
            fn wrap(value: qpa_core::ModuleDecompositionWithMultiplicities<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::AllSubmodulesOfModuleResult<$field>>
            for AllSubmodulesOfModuleResultKind
        {
            fn wrap(value: qpa_core::AllSubmodulesOfModuleResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::MoreMinimalVersionResult<$field>>
            for MoreMinimalVersionResultKind
        {
            fn wrap(value: qpa_core::MoreMinimalVersionResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::CommonDirectSummandResult<$field>>
            for CommonDirectSummandResultKind
        {
            fn wrap(value: qpa_core::CommonDirectSummandResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_core::MaximalCommonDirectSummand<$field>>
            for MaximalCommonDirectSummandKind
        {
            fn wrap(value: qpa_core::MaximalCommonDirectSummand<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::EndOfModuleAsQuiverAlgebraResult<$field>>
            for EndOfModuleAsQuiverAlgebraResultKind
        {
            fn wrap(value: qpa_gbnp::EndOfModuleAsQuiverAlgebraResult<$field>) -> Self {
                Self::$variant(value)
            }
        }

        impl runtime::Wrap<qpa_gbnp::EndOfBasicModuleAsQuiverAlgebraResult<$field>>
            for EndOfBasicModuleAsQuiverAlgebraResultKind
        {
            fn wrap(value: qpa_gbnp::EndOfBasicModuleAsQuiverAlgebraResult<$field>) -> Self {
                Self::$variant(value)
            }
        }
    };
}

impl_resolution_runtime!(RationalField, Rational);
impl_resolution_runtime!(Z2, Binary);
impl_resolution_runtime!(Zp, Prime);

fn clone_module_homomorphisms(
    py: Python<'_>,
    homomorphisms: Vec<Py<PyModuleHomomorphism>>,
) -> Vec<ModuleHomomorphismKind> {
    homomorphisms
        .into_iter()
        .map(|homomorphism| homomorphism.borrow(py).kind.clone())
        .collect()
}

fn global_dimension_text(
    result: Option<qpa_gbnp::GlobalDimensionOfAlgebraResult>,
) -> Option<String> {
    result.map(|result| match result {
        qpa_gbnp::GlobalDimensionOfAlgebraResult::Finite(dimension) => dimension.to_string(),
        qpa_gbnp::GlobalDimensionOfAlgebraResult::Infinite => "infinity".into(),
        qpa_gbnp::GlobalDimensionOfAlgebraResult::ExceedsBound => "exceeds_bound".into(),
    })
}

fn dominant_dimension_text(result: Option<qpa_gbnp::DominantDimensionResult>) -> Option<String> {
    result.map(|result| match result {
        qpa_gbnp::DominantDimensionResult::Finite(dimension) => dimension.to_string(),
        qpa_gbnp::DominantDimensionResult::Infinite => "infinity".into(),
        qpa_gbnp::DominantDimensionResult::ExceedsBound => "exceeds_bound".into(),
    })
}

fn clone_relations(
    py: Python<'_>,
    relations: Vec<Py<PyAlgebraElement>>,
) -> Vec<AlgebraElementKind> {
    relations
        .into_iter()
        .map(|relation| relation.borrow(py).kind.clone())
        .collect()
}

fn wrap_algebra_elements<F: runtime::PythonField>(
    elements: Vec<PathAlgebraElement<F>>,
) -> Vec<PyAlgebraElement>
where
    AlgebraElementKind: runtime::Wrap<PathAlgebraElement<F>>,
{
    runtime::wrap_all::<AlgebraElementKind, _>(elements)
        .into_iter()
        .map(|kind| PyAlgebraElement { kind })
        .collect()
}

fn wrap_groebner_basis<F: runtime::PythonField>(basis: GroebnerBasis<F>) -> PyGroebnerBasis
where
    GroebnerKind: runtime::Wrap<GroebnerBasis<F>>,
{
    PyGroebnerBasis {
        kind: runtime::Wrap::wrap(basis),
    }
}

fn wrap_quotient_elements<F: runtime::PythonField>(
    elements: Vec<QuotientPathAlgebraElement<F>>,
) -> Vec<PyQuotientElement>
where
    QuotientElementKind: runtime::Wrap<QuotientPathAlgebraElement<F>>,
{
    runtime::wrap_all::<QuotientElementKind, _>(elements)
        .into_iter()
        .map(|kind| PyQuotientElement { kind })
        .collect()
}

fn wrap_matrix_modules<F: runtime::PythonField>(
    modules: Vec<PathAlgebraModule<F>>,
) -> Vec<PyMatrixModule>
where
    MatrixModuleKind: runtime::Wrap<PathAlgebraModule<F>>,
{
    runtime::wrap_all::<MatrixModuleKind, _>(modules)
        .into_iter()
        .map(|kind| PyMatrixModule { kind })
        .collect()
}

fn homomorphism_of_projectives_from_path_matrix<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    matrix: Vec<Vec<AlgebraElementKind>>,
    source_vertices: &[usize],
    range_vertices: &[usize],
) -> PyResult<PyModuleHomomorphism>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
    ModuleHomomorphismKind: runtime::Wrap<PathAlgebraModuleHomomorphism<F>>,
{
    let matrix = matrix
        .into_iter()
        .map(|row| {
            Ok(runtime::extract_all(row)?
                .into_iter()
                .map(QuiverAlgebraElement::Path)
                .collect::<Vec<_>>())
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PyModuleHomomorphism {
        kind: runtime::Wrap::wrap(
            qpa_gbnp::from_matrix_to_homomorphism_of_projectives(
                &QuiverAlgebra::Path(algebra.clone()),
                &matrix,
                source_vertices,
                range_vertices,
            )
            .map_err(invalid)?,
        ),
    })
}

fn build_algebra_homomorphism<F: runtime::PythonField>(
    source: QuiverAlgebra<F>,
    range: QuiverAlgebra<F>,
    generator_images: Vec<QuiverAlgebraElement<F>>,
    generators_with_one: bool,
) -> PyResult<PyAlgebraHomomorphism>
where
    AlgebraHomomorphismKind: runtime::Wrap<qpa_core::QuiverAlgebraHomomorphism<F>>,
{
    if source.field() != range.field() {
        return Err(invalid(
            "source and range algebras have different coefficient fields",
        ));
    }
    let homomorphism = if generators_with_one {
        qpa_core::QuiverAlgebraHomomorphism::by_generators_with_one_images(
            source,
            range,
            generator_images,
        )
    } else {
        qpa_core::QuiverAlgebraHomomorphism::by_generator_images(source, range, generator_images)
    }
    .map_err(invalid)?;
    Ok(PyAlgebraHomomorphism {
        kind: runtime::Wrap::wrap(homomorphism),
    })
}

fn path_generator_images<F: runtime::PythonField>(
    images: Vec<AlgebraElementKind>,
) -> PyResult<Vec<QuiverAlgebraElement<F>>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    Ok(runtime::extract_all(images)?
        .into_iter()
        .map(QuiverAlgebraElement::Path)
        .collect())
}

fn quotient_generator_images<F: runtime::PythonField>(
    images: Vec<QuotientElementKind>,
) -> PyResult<Vec<QuiverAlgebraElement<F>>>
where
    QuotientElementKind: runtime::Extract<QuotientPathAlgebraElement<F>>,
{
    Ok(runtime::extract_all(images)?
        .into_iter()
        .map(QuiverAlgebraElement::Quotient)
        .collect())
}

fn build_ideal<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    generators: Vec<AlgebraElementKind>,
) -> PyResult<PathAlgebraIdeal<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    PathAlgebraIdeal::new(algebra.clone(), runtime::extract_all(generators)?).map_err(invalid)
}

fn build_groebner_basis<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<GroebnerBasis<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    GroebnerBasis::high_level(algebra, &runtime::extract_all(relations)?).map_err(invalid)
}

#[pymethods]
impl PyAlgebraHomomorphism {
    fn enveloping_algebra_homomorphism(&self) -> PyResult<Self> {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            tensor::enveloping_homomorphism_for(homomorphism)
        })
    }

    #[getter]
    fn source(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            conversion::wrap_quiver_algebra(homomorphism.source().clone())
        })
    }

    #[getter]
    fn range(&self) -> PyQuiverAlgebra {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            conversion::wrap_quiver_algebra(homomorphism.range().clone())
        })
    }

    #[getter]
    fn generators(&self) -> Vec<PyQuiverAlgebraElement> {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            tensor::wrap_quiver_elements(homomorphism.generators())
        })
    }

    #[getter]
    fn generator_images(&self) -> Vec<PyQuiverAlgebraElement> {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            tensor::wrap_quiver_elements(homomorphism.generator_images())
        })
    }

    fn image(&self, element: &PyQuiverAlgebraElement) -> PyResult<PyQuiverAlgebraElement> {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            tensor::homomorphism_image_for(homomorphism, element.kind.clone())
        })
    }

    /// Python spelling of QPA's retained-generator `MappedExpression` route.
    fn mapped_expression(
        &self,
        element: &PyQuiverAlgebraElement,
    ) -> PyResult<PyQuiverAlgebraElement> {
        self.image(element)
    }

    /// Python spelling of QPA's `ImagesRepresentative` method.
    fn images_representative(
        &self,
        element: &PyQuiverAlgebraElement,
    ) -> PyResult<PyQuiverAlgebraElement> {
        self.image(element)
    }

    #[staticmethod]
    #[pyo3(signature = (source, range, generator_images, generators_with_one=false))]
    fn path_to_path(
        py: Python<'_>,
        source: &PyPathAlgebra,
        range: &PyPathAlgebra,
        generator_images: Vec<Py<PyAlgebraElement>>,
        generators_with_one: bool,
    ) -> PyResult<Self> {
        let generator_images = clone_relations(py, generator_images);
        dispatch_pair_ref!(
            &source.kind,
            &range.kind,
            AlgebraKind,
            |source, range| build_algebra_homomorphism(
                QuiverAlgebra::Path(source.clone()),
                QuiverAlgebra::Path(range.clone()),
                path_generator_images(generator_images)?,
                generators_with_one,
            ),
            Err(invalid(
                "source and range algebras have different coefficient fields"
            ))
        )
    }

    #[staticmethod]
    #[pyo3(signature = (source, range, generator_images, generators_with_one=false))]
    fn path_to_quotient(
        py: Python<'_>,
        source: &PyPathAlgebra,
        range: &PyQuotientAlgebra,
        generator_images: Vec<Py<PyQuotientElement>>,
        generators_with_one: bool,
    ) -> PyResult<Self> {
        let generator_images = generator_images
            .into_iter()
            .map(|element| element.borrow(py).kind.clone())
            .collect();
        match (&source.kind, &range.kind) {
            (AlgebraKind::Rational(source), QuotientKind::Rational(range)) => {
                build_algebra_homomorphism(
                    QuiverAlgebra::Path(source.clone()),
                    QuiverAlgebra::Quotient(range.clone()),
                    quotient_generator_images(generator_images)?,
                    generators_with_one,
                )
            }
            (AlgebraKind::Binary(source), QuotientKind::Binary(range)) => {
                build_algebra_homomorphism(
                    QuiverAlgebra::Path(source.clone()),
                    QuiverAlgebra::Quotient(range.clone()),
                    quotient_generator_images(generator_images)?,
                    generators_with_one,
                )
            }
            (AlgebraKind::Prime(source), QuotientKind::Prime(range)) => build_algebra_homomorphism(
                QuiverAlgebra::Path(source.clone()),
                QuiverAlgebra::Quotient(range.clone()),
                quotient_generator_images(generator_images)?,
                generators_with_one,
            ),
            _ => Err(invalid(
                "source and range algebras have different coefficient fields",
            )),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (source, range, generator_images, generators_with_one=false))]
    fn quotient_to_path(
        py: Python<'_>,
        source: &PyQuotientAlgebra,
        range: &PyPathAlgebra,
        generator_images: Vec<Py<PyAlgebraElement>>,
        generators_with_one: bool,
    ) -> PyResult<Self> {
        let generator_images = clone_relations(py, generator_images);
        match (&source.kind, &range.kind) {
            (QuotientKind::Rational(source), AlgebraKind::Rational(range)) => {
                build_algebra_homomorphism(
                    QuiverAlgebra::Quotient(source.clone()),
                    QuiverAlgebra::Path(range.clone()),
                    path_generator_images(generator_images)?,
                    generators_with_one,
                )
            }
            (QuotientKind::Binary(source), AlgebraKind::Binary(range)) => {
                build_algebra_homomorphism(
                    QuiverAlgebra::Quotient(source.clone()),
                    QuiverAlgebra::Path(range.clone()),
                    path_generator_images(generator_images)?,
                    generators_with_one,
                )
            }
            (QuotientKind::Prime(source), AlgebraKind::Prime(range)) => build_algebra_homomorphism(
                QuiverAlgebra::Quotient(source.clone()),
                QuiverAlgebra::Path(range.clone()),
                path_generator_images(generator_images)?,
                generators_with_one,
            ),
            _ => Err(invalid(
                "source and range algebras have different coefficient fields",
            )),
        }
    }

    #[staticmethod]
    #[pyo3(signature = (source, range, generator_images, generators_with_one=false))]
    fn quotient_to_quotient(
        py: Python<'_>,
        source: &PyQuotientAlgebra,
        range: &PyQuotientAlgebra,
        generator_images: Vec<Py<PyQuotientElement>>,
        generators_with_one: bool,
    ) -> PyResult<Self> {
        let generator_images = generator_images
            .into_iter()
            .map(|element| element.borrow(py).kind.clone())
            .collect();
        dispatch_pair_ref!(
            &source.kind,
            &range.kind,
            QuotientKind,
            |source, range| build_algebra_homomorphism(
                QuiverAlgebra::Quotient(source.clone()),
                QuiverAlgebra::Quotient(range.clone()),
                quotient_generator_images(generator_images)?,
                generators_with_one,
            ),
            Err(invalid(
                "source and range algebras have different coefficient fields"
            ))
        )
    }

    #[getter]
    fn source_kind(&self) -> &'static str {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            match homomorphism.source() {
                QuiverAlgebra::Path(_) => "path",
                QuiverAlgebra::Quotient(_) => "quotient",
            }
        })
    }

    #[getter]
    fn range_kind(&self) -> &'static str {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            match homomorphism.range() {
                QuiverAlgebra::Path(_) => "path",
                QuiverAlgebra::Quotient(_) => "quotient",
            }
        })
    }

    #[getter]
    fn source_path_algebra(&self) -> Option<PyPathAlgebra> {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            match homomorphism.source() {
                QuiverAlgebra::Path(algebra) => Some(PyPathAlgebra {
                    kind: runtime::Wrap::wrap(algebra.clone()),
                }),
                QuiverAlgebra::Quotient(_) => None,
            }
        })
    }

    #[getter]
    fn source_quotient_algebra(&self) -> Option<PyQuotientAlgebra> {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            match homomorphism.source() {
                QuiverAlgebra::Path(_) => None,
                QuiverAlgebra::Quotient(algebra) => Some(PyQuotientAlgebra {
                    kind: runtime::Wrap::wrap(algebra.clone()),
                }),
            }
        })
    }

    #[getter]
    fn range_path_algebra(&self) -> Option<PyPathAlgebra> {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            match homomorphism.range() {
                QuiverAlgebra::Path(algebra) => Some(PyPathAlgebra {
                    kind: runtime::Wrap::wrap(algebra.clone()),
                }),
                QuiverAlgebra::Quotient(_) => None,
            }
        })
    }

    #[getter]
    fn range_quotient_algebra(&self) -> Option<PyQuotientAlgebra> {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            match homomorphism.range() {
                QuiverAlgebra::Path(_) => None,
                QuiverAlgebra::Quotient(algebra) => Some(PyQuotientAlgebra {
                    kind: runtime::Wrap::wrap(algebra.clone()),
                }),
            }
        })
    }

    #[getter]
    fn generator_count(&self) -> usize {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            homomorphism.generators().len()
        })
    }

    #[getter]
    fn uses_generators_with_one(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraHomomorphismKind, |homomorphism| {
            homomorphism.uses_generators_with_one()
        })
    }

    fn opposite(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraHomomorphismKind => AlgebraHomomorphismKind,
            |homomorphism| qpa_gbnp::opposite_algebra_homomorphism(homomorphism)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn __eq__(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            AlgebraHomomorphismKind,
            |left, right| left == right,
            false
        )
    }
}

#[pymethods]
impl PyPathAlgebra {
    #[new]
    fn new(field: PyRef<'_, PyField>, quiver: PyRef<'_, PyQuiver>) -> PyResult<Self> {
        let kind = match field.kind {
            FieldKind::Rational => AlgebraKind::Rational(match quiver.rust_ordering() {
                Some(ordering) => {
                    PathAlgebra::with_ordering(Q, quiver.rust_quiver().clone(), ordering.clone())
                        .map_err(invalid)?
                }
                None => PathAlgebra::new(Q, quiver.rust_quiver().clone()).map_err(invalid)?,
            }),
            FieldKind::Binary => AlgebraKind::Binary(match quiver.rust_ordering() {
                Some(ordering) => {
                    PathAlgebra::with_ordering(Z2, quiver.rust_quiver().clone(), ordering.clone())
                        .map_err(invalid)?
                }
                None => PathAlgebra::new(Z2, quiver.rust_quiver().clone()).map_err(invalid)?,
            }),
            FieldKind::Prime(characteristic) => AlgebraKind::Prime(match quiver.rust_ordering() {
                Some(ordering) => PathAlgebra::with_ordering(
                    Zp::new(characteristic),
                    quiver.rust_quiver().clone(),
                    ordering.clone(),
                )
                .map_err(invalid)?,
                None => PathAlgebra::new(Zp::new(characteristic), quiver.rust_quiver().clone())
                    .map_err(invalid)?,
            }),
        };
        Ok(Self { kind })
    }

    fn right_projective_module(&self, vertices: Vec<usize>) -> PyResult<PyRightProjectiveModule> {
        projective::right_projective_module(self, vertices)
    }

    fn preprojective_algebra(&self) -> PyResult<PyQuiverAlgebra> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            predefined::preprojective_for(algebra)
        })
    }

    fn complexity(&self, degree: i64) -> PyResult<usize> {
        let degree = Integer::from(degree);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::complexity_of_algebra(&QuiverAlgebra::Path(algebra.clone()), &degree)
                .map_err(invalid)
        })
    }

    fn centre(&self) -> PyResult<PyPathAlgebraCentre> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            finite_structures::path_centre_for(algebra)
        })
    }

    fn units(&self) -> PyResult<PyQuiverAlgebraUnitGroup> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            finite_structures::units_for(&QuiverAlgebra::Path(algebra.clone()))
        })
    }

    #[getter]
    fn has_elements_of_units(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            QuiverAlgebra::Path(algebra.clone()).has_elements_of_units()
        })
    }

    fn elements_of_units(&self) -> PyResult<Vec<PyQuiverAlgebraElement>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            finite_structures::elements_of_units_for(&QuiverAlgebra::Path(algebra.clone()))
        })
    }

    fn right_groebner_basis_from_generators(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyPathAlgebraModuleGroebnerBasis> {
        projective::right_groebner_basis_from_generators(self, py, generators)
    }

    fn vertex_projective_presentation(
        &self,
        py: Python<'_>,
        rows: Vec<Vec<Py<PyAlgebraElement>>>,
    ) -> PyResult<PyRightFpPathAlgebraModule> {
        fp_projective::vertex_presentation(self, py, rows)
    }

    fn vertex_projective_presentation_with_vertices(
        &self,
        py: Python<'_>,
        rows: Vec<Vec<Py<PyAlgebraElement>>>,
        vertices: Vec<usize>,
    ) -> PyResult<PyRightFpPathAlgebraModule> {
        fp_projective::vertex_presentation_with_vertices(self, py, rows, vertices)
    }

    #[getter]
    fn field(&self) -> PyField {
        let kind = dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .field()
            .python_kind());
        PyField { kind }
    }

    #[getter]
    fn quiver(&self) -> PyQuiver {
        let inner = dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra.quiver().clone());
        let ordering = dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .ordering()
            .clone());
        PyQuiver::from_inner_with_ordering(inner, ordering)
    }

    /// QPA's `OriginalPathAlgebra(A)` returns a path algebra unchanged.
    #[getter]
    fn original_path_algebra(&self) -> Self {
        self.clone()
    }

    #[getter]
    fn ordering(&self) -> PyQuiverOrdering {
        let inner = dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .ordering()
            .clone());
        PyQuiverOrdering::from_inner(inner)
    }

    fn adr_algebra(&self) -> PyResult<PyAdrAlgebraOfAlgebraResult> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AdrAlgebraOfAlgebraResultKind,
            |algebra| qpa_gbnp::adr_algebra_of_algebra(&QuiverAlgebra::Path(algebra.clone()))
                .map_err(invalid)?
        );
        Ok(PyAdrAlgebraOfAlgebraResult { kind })
    }

    #[getter]
    fn global_dimension(&self) -> usize {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .qpa_global_dimension())
    }

    fn global_dimension_bounded(&self, bound: usize) -> PyResult<Option<String>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            Ok(global_dimension_text(
                qpa_gbnp::global_dimension_of_algebra(&QuiverAlgebra::Path(algebra.clone()), bound)
                    .map_err(invalid)?,
            ))
        })
    }

    #[getter]
    fn cached_gorenstein_dimension(&self) -> Option<usize> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::gorenstein_dimension(&QuiverAlgebra::Path(algebra.clone()))
        })
    }

    fn gorenstein_dimension_bounded(&self, bound: i64) -> PyResult<Option<usize>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::gorenstein_dimension_of_algebra(&QuiverAlgebra::Path(algebra.clone()), &bound)
                .map_err(invalid)
        })
    }

    fn dominant_dimension_bounded(&self, bound: i64) -> PyResult<Option<String>> {
        let bound = Integer::from(bound);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| Ok(
            dominant_dimension_text(
                qpa_gbnp::dominant_dimension_of_algebra(
                    &QuiverAlgebra::Path(algebra.clone()),
                    &bound,
                )
                .map_err(invalid)?,
            )
        ))
    }

    #[getter]
    fn has_finite_canonical_basis(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .has_finite_canonical_basis())
    }

    /// Return QPA's `IsPathRing` property.
    #[getter]
    fn is_path_ring(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra.is_path_ring())
    }

    #[getter]
    fn is_monomial(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_monomial(algebra)
        })
    }

    #[getter]
    fn is_special_biserial(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_special_biserial(algebra)
        })
    }

    #[getter]
    fn is_string(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_string(algebra)
        })
    }

    #[getter]
    fn is_gentle(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::path_algebra_is_gentle(algebra)
        })
    }

    #[getter]
    fn is_gorenstein(&self) -> PyResult<bool> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::gentle_algebra_is_gorenstein(&QuiverAlgebra::Path(algebra.clone()))
                .map_err(invalid)
        })
    }

    #[getter]
    fn is_hereditary(&self) -> PyResult<Option<bool>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::is_hereditary_algebra(&QuiverAlgebra::Path(algebra.clone())).map_err(invalid)
        })
    }

    #[getter]
    fn is_finite_global_dimension(&self) -> Option<bool> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            qpa_gbnp::is_finite_global_dimension_from_known_properties(&QuiverAlgebra::Path(
                algebra.clone(),
            ))
        })
    }

    #[getter]
    fn generators(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| wrap_algebra_elements(
            algebra.generators()
        ))
    }

    #[getter]
    fn relations(&self) -> Vec<PyAlgebraElement> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| wrap_algebra_elements(
            algebra.relations()
        ))
    }

    fn canonical_basis(&self) -> Option<Vec<PyAlgebraElement>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .canonical_basis()
            .map(wrap_algebra_elements))
    }

    fn opposite(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraKind,
            |algebra| algebra.opposite().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn generator(&self, name: &str) -> PyResult<PyAlgebraElement> {
        let path = dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .quiver()
            .generator_by_name(name)
            .ok_or_else(|| invalid(format!("unknown quiver generator {name:?}")))?
            .path());
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra.element(path).map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    /// Embed a path through QPA's `ElementOfPathAlgebra` operation.
    fn element(&self, path: &PyPath) -> PyResult<PyAlgebraElement> {
        let path = path.inner.clone();
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra.element(path).map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    /// Rebuild a path-algebra element from QPA's canonical external form.
    fn element_from_external_representation(
        &self,
        py: Python<'_>,
        representation: Vec<(Vec<usize>, Py<PyScalar>)>,
    ) -> PyResult<PyAlgebraElement> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra
                .element_from_external_representation(&path_external_representation(
                    py,
                    algebra,
                    representation,
                )?)
                .map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn term(&self, path: &PyPath, coefficient: &PyScalar) -> PyResult<PyAlgebraElement> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra
                .from_terms([(
                    path.inner.clone(),
                    coefficient.value_for(algebra.field())?,
                )])
                .map_err(invalid)?
        );
        Ok(PyAlgebraElement { kind })
    }

    fn zero(&self) -> PyAlgebraElement {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra.zero()
        );
        PyAlgebraElement { kind }
    }

    fn one(&self) -> PyAlgebraElement {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => AlgebraElementKind,
            |algebra| algebra.one()
        );
        PyAlgebraElement { kind }
    }

    fn quotient(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyQuotientAlgebra> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => QuotientKind,
            |algebra| build_quotient(algebra, relations)?
        );
        Ok(PyQuotientAlgebra { kind })
    }

    fn quotient_high_level(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyQuotientAlgebra> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => QuotientKind,
            |algebra| build_high_level_quotient(algebra, relations)?
        );
        Ok(PyQuotientAlgebra { kind })
    }

    fn ideal(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyPathAlgebraIdeal> {
        let generators = clone_relations(py, generators);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => IdealKind,
            |algebra| build_ideal(algebra, generators)?
        );
        Ok(PyPathAlgebraIdeal { kind })
    }

    fn groebner_basis(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyGroebnerBasis> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => GroebnerKind,
            |algebra| build_groebner_basis(algebra, relations)?
        );
        Ok(PyGroebnerBasis { kind })
    }

    fn gbnp_groebner_basis(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyGroebnerBasis> {
        let relations = clone_relations(py, relations);
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => GroebnerKind,
            |algebra| build_gbnp_groebner_basis(algebra, relations)?
        );
        Ok(PyGroebnerBasis { kind })
    }

    /// Invoke the `GroebnerBasisFunction` selected by QPA's default
    /// `PathAlgebra` constructor, namely `GBNPGroebnerBasis`.
    fn groebner_basis_function(
        &self,
        py: Python<'_>,
        relations: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyGroebnerBasis> {
        self.gbnp_groebner_basis(py, relations)
    }

    fn make_uniform_on_right(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<Vec<PyAlgebraElement>> {
        let generators = clone_relations(py, generators);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            let generators = runtime::extract_all(generators)?;
            Ok(wrap_algebra_elements(
                algebra
                    .make_uniform_on_right(&generators)
                    .map_err(invalid)?,
            ))
        })
    }

    fn generators_times_arrows_on_right(
        &self,
        py: Python<'_>,
        generators: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<Vec<PyAlgebraElement>> {
        let generators = clone_relations(py, generators);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            let generators = runtime::extract_all(generators)?;
            Ok(wrap_algebra_elements(
                algebra
                    .generators_times_arrows_on_right(&generators)
                    .map_err(invalid)?,
            ))
        })
    }

    fn tip_monomial_and_coefficient_of_vector(
        &self,
        py: Python<'_>,
        vector: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<(usize, PyPath, PyScalar)> {
        let vector = clone_relations(py, vector);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            vector_reduction::path_tip_for(algebra, vector)
        })
    }

    /// Mutate and return the identical outer Python list, matching QPA.
    fn tip_reduce_vectors(&self, vectors: &Bound<'_, PyList>) -> PyResult<Option<Py<PyList>>> {
        let vector_kinds = vector_reduction::path_rows(vectors)?;
        let reduced = dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            vector_reduction::path_reduce_for(algebra, vector_kinds)
        })?;
        let Some(reduced) = reduced else {
            return Ok(None);
        };
        vector_reduction::replace_path_rows(vectors, reduced)?;
        Ok(Some(vectors.clone().unbind()))
    }

    fn coefficients_of_vectors(
        &self,
        py: Python<'_>,
        basis: Vec<Vec<Py<PyAlgebraElement>>>,
        vector: Vec<Py<PyAlgebraElement>>,
    ) -> PyResult<PyPathVectorCoefficientsResult> {
        let basis = basis
            .into_iter()
            .map(|row| clone_relations(py, row))
            .collect();
        let vector = clone_relations(py, vector);
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            vector_reduction::path_coefficients_for(algebra, basis, vector)
        })
    }

    /// Return QPA's generator-global assignment plan without mutating Python globals.
    #[pyo3(signature = (read_only_names=Vec::new()))]
    fn generator_variable_assignments(
        &self,
        read_only_names: Vec<String>,
    ) -> PyResult<Vec<(String, PyQuiverAlgebraElement)>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            mapping::generator_variable_assignments_for(
                &QuiverAlgebra::Path(algebra.clone()),
                &read_only_names,
            )
        })
    }

    fn nth_power_of_arrow_ideal(&self, exponent: i64) -> PyResult<Vec<PyAlgebraElement>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| Ok(
            wrap_algebra_elements(
                algebra
                    .nth_power_of_arrow_ideal(exponent)
                    .map_err(invalid)?,
            )
        ))
    }

    fn left_support_of_quiver_algebra_element(
        &self,
        element: &PyAlgebraElement,
    ) -> PyResult<Vec<usize>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            Ok(element_utilities::path_support_for(
                algebra,
                element.kind.clone(),
                element_utilities::SupportSide::Left,
            )?
            .pop()
            .unwrap_or_default())
        })
    }

    fn right_support_of_quiver_algebra_element(
        &self,
        element: &PyAlgebraElement,
    ) -> PyResult<Vec<usize>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            Ok(element_utilities::path_support_for(
                algebra,
                element.kind.clone(),
                element_utilities::SupportSide::Right,
            )?
            .pop()
            .unwrap_or_default())
        })
    }

    fn support_of_quiver_algebra_element(
        &self,
        element: &PyAlgebraElement,
    ) -> PyResult<Vec<Vec<usize>>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            element_utilities::path_support_for(
                algebra,
                element.kind.clone(),
                element_utilities::SupportSide::Both,
            )
        })
    }

    /// Mutate and return the identical Python list, as QPA mutates and returns `rels`.
    fn add_nth_power_to_relations(
        &self,
        relations: &Bound<'_, PyList>,
        exponent: i64,
    ) -> PyResult<Py<PyList>> {
        let relation_kinds = relations
            .iter()
            .map(|relation| {
                relation
                    .extract::<PyRef<'_, PyAlgebraElement>>()
                    .map(|relation| relation.kind.clone())
                    .map_err(Into::into)
            })
            .collect::<PyResult<Vec<_>>>()?;
        let appended = dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            element_utilities::add_power_for(algebra, relation_kinds, exponent)
        })?;
        for relation in appended {
            relations.append(Py::new(relations.py(), relation)?)?;
        }
        Ok(relations.clone().unbind())
    }

    fn module_from_integer_matrices(
        &self,
        dimension_vector: Vec<usize>,
        actions: Vec<(String, Vec<Vec<i64>>)>,
    ) -> PyResult<PyMatrixModule> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => MatrixModuleKind,
            |algebra| build_path_module(algebra, dimension_vector, actions)?
        );
        Ok(PyMatrixModule { kind })
    }

    fn hom_of_fp_path_algebra_modules(
        &self,
        source: &PyRightFpPathAlgebraModule,
        range: &PyRightFpPathAlgebraModule,
    ) -> PyResult<()> {
        fp_projective::hom(self, source, range)
    }

    fn end_of_fp_path_algebra_module(&self, module: &PyRightFpPathAlgebraModule) -> PyResult<()> {
        fp_projective::end(self, module)
    }

    fn homomorphism_of_projectives_from_matrix(
        &self,
        py: Python<'_>,
        matrix: Vec<Vec<Py<PyAlgebraElement>>>,
        source_vertices: Vec<usize>,
        range_vertices: Vec<usize>,
    ) -> PyResult<PyModuleHomomorphism> {
        let matrix = matrix
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|element| element.borrow(py).kind.clone())
                    .collect()
            })
            .collect();
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| {
            homomorphism_of_projectives_from_path_matrix(
                algebra,
                matrix,
                &source_vertices,
                &range_vertices,
            )
        })
    }

    fn zero_module(&self) -> PyResult<PyMatrixModule> {
        let kind = map_ref!(
            &self.kind,
            AlgebraKind => MatrixModuleKind,
            |algebra| qpa_core::zero_module(&QuiverAlgebra::Path(algebra.clone()))
                .map_err(invalid)?
        );
        Ok(PyMatrixModule { kind })
    }

    fn simple_modules(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| Ok(wrap_matrix_modules(
            qpa_gbnp::simple_modules(&QuiverAlgebra::Path(algebra.clone())).map_err(invalid)?,
        )))
    }

    fn indecomposable_projectives(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| Ok(wrap_matrix_modules(
            qpa_gbnp::indec_projective_modules(&QuiverAlgebra::Path(algebra.clone()))
                .map_err(invalid)?,
        )))
    }

    fn indecomposable_injectives(&self) -> PyResult<Vec<PyMatrixModule>> {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| Ok(wrap_matrix_modules(
            qpa_gbnp::indec_injective_modules(&QuiverAlgebra::Path(algebra.clone()))
                .map_err(invalid)?,
        )))
    }

    fn is_identical(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            AlgebraKind,
            |left, right| left.is_identical(right),
            false
        )
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, AlgebraKind, |algebra| algebra
            .view_obj_with(&algebra.field().python_kind().display()))
    }
}

fn build_quotient<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<QuotientPathAlgebra<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let relations = runtime::extract_all(relations)?;
    qpa_gbnp::quotient_by_relations(algebra.clone(), &relations).map_err(invalid)
}

fn path_external_representation<F: runtime::PythonField>(
    py: Python<'_>,
    algebra: &PathAlgebra<F>,
    representation: Vec<(Vec<usize>, Py<PyScalar>)>,
) -> PyResult<PathAlgebraExternalRepresentation<F::Element>> {
    let terms = representation
        .into_iter()
        .map(|(path, coefficient)| {
            Ok((path, algebra.field().scalar_value(&coefficient.borrow(py))?))
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(PathAlgebraExternalRepresentation::new(
        algebra.field().integer_value(&Integer::from(0)),
        terms,
    ))
}

fn build_high_level_quotient<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<QuotientPathAlgebra<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let relations = runtime::extract_all(relations)?;
    QuotientPathAlgebra::by_relations_high_level(algebra.clone(), &relations).map_err(invalid)
}

fn build_gbnp_groebner_basis<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    relations: Vec<AlgebraElementKind>,
) -> PyResult<GroebnerBasis<F>>
where
    AlgebraElementKind: runtime::Extract<PathAlgebraElement<F>>,
{
    let relations = runtime::extract_all(relations)?;
    let complete = qpa_gbnp::gbnp_groebner_basis(&relations, algebra).map_err(invalid)?;
    GroebnerBasis::from_complete_relations(algebra, &complete).map_err(invalid)
}

fn build_ideal_quotient<F: runtime::PythonField>(
    ideal: &PathAlgebraIdeal<F>,
) -> PyResult<QuotientPathAlgebra<F>> {
    if ideal.groebner_basis().is_none() {
        let complete = qpa_gbnp::gbnp_groebner_basis_of_ideal(ideal).map_err(invalid)?;
        ideal
            .with_complete_groebner_basis(&complete)
            .map_err(invalid)?;
    }
    QuotientPathAlgebra::by_ideal(ideal).map_err(invalid)
}

fn build_path_module<F: runtime::PythonField>(
    algebra: &PathAlgebra<F>,
    dimension_vector: Vec<usize>,
    actions: Vec<(String, Vec<Vec<i64>>)>,
) -> PyResult<PathAlgebraModule<F>> {
    build_matrix_module(
        QuiverAlgebra::Path(algebra.clone()),
        dimension_vector,
        actions,
    )
}

fn build_matrix_module<F: runtime::PythonField>(
    algebra: QuiverAlgebra<F>,
    dimension_vector: Vec<usize>,
    actions: Vec<(String, Vec<Vec<i64>>)>,
) -> PyResult<PathAlgebraModule<F>> {
    let actions = actions
        .into_iter()
        .map(|(label, rows)| {
            let rows = rows
                .into_iter()
                .map(|row| {
                    row.into_iter()
                        .map(|value| algebra.field().integer_value(&Integer::from(value)))
                        .collect()
                })
                .collect();
            let matrix = DenseMatrix::from_rows(rows, algebra.field().clone()).map_err(invalid)?;
            Ok(NamedArrowMatrix::new(label, matrix))
        })
        .collect::<PyResult<Vec<_>>>()?;
    right_module_over_path_algebra(algebra, dimension_vector, &actions)
        .map_err(invalid)?
        .ok_or_else(|| {
            RelationNotSatisfiedError::new_err("module actions do not satisfy relations")
        })
}

#[pymethods]
impl PyAlgebraElement {
    #[getter]
    fn algebra(&self) -> PyPathAlgebra {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraKind,
            |element| element.algebra().clone()
        );
        PyPathAlgebra { kind }
    }

    #[getter]
    fn terms(&self) -> Vec<(PyPath, PyScalar)> {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .terms()
            .iter()
            .map(|term| (
                PyPath {
                    inner: term.path().clone(),
                },
                element.algebra().field().wrap_scalar(term.coefficient()),
            ))
            .collect())
    }

    #[getter]
    fn external_representation(&self) -> Vec<(Vec<usize>, PyScalar)> {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .external_representation()
            .terms()
            .iter()
            .map(|(path, coefficient)| (
                path.clone(),
                element.algebra().field().wrap_scalar(coefficient),
            ))
            .collect())
    }

    #[getter]
    fn is_zero(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element.is_zero())
    }

    #[getter]
    fn is_left_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .is_left_uniform())
    }

    #[getter]
    fn is_right_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .is_right_uniform())
    }

    #[getter]
    fn is_uniform(&self) -> bool {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .is_uniform())
    }

    #[getter]
    fn leading_term(&self) -> Self {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.leading_term()
        );
        Self { kind }
    }

    #[getter]
    fn leading_coefficient(&self) -> PyScalar {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .algebra()
            .field()
            .wrap_scalar(&element.leading_coefficient()))
    }

    #[getter]
    fn leading_monomial(&self) -> PyPath {
        let inner = dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .leading_monomial());
        PyPath { inner }
    }

    #[getter]
    fn tip(&self) -> Self {
        self.leading_term()
    }

    #[getter]
    fn tip_coefficient(&self) -> PyScalar {
        self.leading_coefficient()
    }

    #[getter]
    fn tip_monomial(&self) -> PyPath {
        self.leading_monomial()
    }

    fn path_removal(&self, py: Python<'_>, vertex_positions: Vec<Py<PyAny>>) -> PyResult<Self> {
        let positions = vertex_positions
            .into_iter()
            .map(|position| parse_python_integer(position.bind(py)))
            .collect::<PyResult<Vec<_>>>()?
            .into_iter()
            .filter_map(|position| position.to_string().parse::<usize>().ok())
            .collect::<Vec<_>>();
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| {
            element_utilities::path_removal_for(element, &positions)
        })
    }

    #[getter]
    fn vertex_position(&self) -> PyResult<usize> {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .vertex_position()
            .map_err(invalid))
    }

    fn add(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind => AlgebraElementKind,
            |left, right| left.add(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn negate(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.neg().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn subtract(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind => AlgebraElementKind,
            |left, right| left.sub(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn multiply(&self, right: &Self) -> PyResult<Self> {
        let kind = map_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind => AlgebraElementKind,
            |left, right| left.multiply(right).map_err(invalid)?,
            invalid("elements belong to different coefficient fields")
        )?;
        Ok(Self { kind })
    }

    fn pow(&self, exponent: u64) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.pow(exponent).map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_left(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element
                .scale_left(&scalar.value_for(element.algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn scale_right(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element
                .scale_right(&scalar.value_for(element.algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn divide_scalar(&self, scalar: &PyScalar) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element
                .divide_scalar(&scalar.value_for(element.algebra().field())?)
                .map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn compare(&self, right: &Self) -> PyResult<i8> {
        dispatch_pair_ref!(
            &self.kind,
            &right.kind,
            AlgebraElementKind,
            |left, right| Ok(ordering_value(left.compare(right).map_err(invalid)?)),
            Err(invalid("elements belong to different coefficient fields"))
        )
    }

    fn opposite(&self) -> PyResult<Self> {
        let kind = map_ref!(
            &self.kind,
            AlgebraElementKind => AlgebraElementKind,
            |element| element.opposite().map_err(invalid)?
        );
        Ok(Self { kind })
    }

    fn __add__(&self, right: &Self) -> PyResult<Self> {
        self.add(right)
    }

    fn __sub__(&self, right: &Self) -> PyResult<Self> {
        self.subtract(right)
    }

    fn __neg__(&self) -> PyResult<Self> {
        self.negate()
    }

    fn __mul__(&self, right: &Self) -> PyResult<Self> {
        self.multiply(right)
    }

    fn __repr__(&self) -> String {
        dispatch_ref!(&self.kind, AlgebraElementKind, |element| element
            .print_obj_with(|value| element.algebra().field().scalar_text(value)))
    }

    fn __eq__(&self, other: &Self) -> bool {
        dispatch_pair_ref!(
            &self.kind,
            &other.kind,
            AlgebraElementKind,
            |left, right| left == right,
            false
        )
    }
}
