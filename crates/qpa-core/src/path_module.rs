//! Matrix representations of QPA right modules over quiver algebras.
//!
//! The checked constructor in this module translates the three-argument
//! `RightModuleOverPathAlgebra` method from `lib/module.gi`. In particular,
//! arrow matrices act on row vectors, omitted arrows receive GAP's physical
//! `1`-dimensional stand-ins at zero-dimensional vertices, and every defining
//! relation of a quotient is evaluated before a module is returned.

use std::{
    fmt::{self, Write as _},
    sync::{Arc, OnceLock, RwLock, Weak},
};

use numerica::domains::Field;
use thiserror::Error;

use crate::{
    ArrowIndex, CachedAttribute, DenseMatrix, Generator, MatrixError, Path, PathEndpoint,
    QpaCoefficientField, QpaFieldCharacteristic, QuiverAlgebra, QuiverAlgebraElement,
    QuiverAlgebraElementError, QuiverAlgebraHomomorphism, QuiverAlgebraHomomorphismError,
    QuotientError, StructureConstantsAlgebra, StructureConstantsError, TensorProductAlgebra,
};

// These private source chunks remain in this module through `include!`.
// That preserves every established `qpa_core::path_module::*` path while
// making the QPA/domain boundaries reviewable in isolation.
include!("path_module/objects.rs");
include!("path_module/endomorphisms.rs");
include!("path_module/homomorphisms.rs");
include!("path_module/submodules.rs");
include!("path_module/categorical.rs");
include!("path_module/constructors.rs");
include!("path_module/relation_checking.rs");
include!("path_module/errors.rs");

#[cfg(test)]
#[path = "path_module/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "path_module/action_tests.rs"]
mod action_tests;

#[cfg(test)]
#[path = "path_module/relation_tests.rs"]
mod relation_tests;

#[cfg(test)]
#[path = "path_module/projective_cache_tests.rs"]
mod projective_cache_tests;
