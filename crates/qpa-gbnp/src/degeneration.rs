//! Numerical Auslander--Reiten quivers and finite-type degeneration order.
//!
//! This is a source-faithful translation of QPA 1.37
//! `lib/degorderfinitetype.gi`.  In particular, the Hom-dimension matrix is
//! obtained by constructing QPA's Auslander--Reiten matrix and computing
//! `TransposedMat(Inverse(T))`; no independently chosen degeneration
//! algorithm is substituted.

use std::{fmt, fmt::Write as _};

use numerica::{
    domains::{integer::Integer, rational::Q},
    tensors::matrix::Matrix,
};
use thiserror::Error;

/// A multiplicity or dimension vector over Numerica's arbitrary-precision integers.
pub type NumericalVector = Vec<Integer>;

/// A matrix over Numerica's arbitrary-precision integers.
pub type NumericalMatrix = Vec<NumericalVector>;

/// One mesh in QPA's explicit numerical AR-quiver description.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArMesh {
    /// Direct predecessors of the module at this row, using zero-based indices.
    pub direct_predecessors: Vec<usize>,
    /// The zero-based index of its AR translate, or `None` for a projective.
    pub tau: Option<usize>,
}

impl ArMesh {
    /// Construct one explicit mesh.
    #[must_use]
    pub fn new(direct_predecessors: Vec<usize>, tau: Option<usize>) -> Self {
        Self {
            direct_predecessors,
            tau,
        }
    }
}

/// A position `tau^(-shift)(representative)` in an orbit description.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ArOrbitPosition {
    /// Zero-based orbit number.
    pub orbit: usize,
    /// Nonnegative inverse-translate shift from the representative.
    pub shift: usize,
}

impl ArOrbitPosition {
    /// Construct an orbit position.
    #[must_use]
    pub const fn new(orbit: usize, shift: usize) -> Self {
        Self { orbit, shift }
    }
}

/// One orbit in QPA's compact numerical AR-quiver description.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArOrbit {
    /// Orbit length minus one, exactly QPA's first orbit entry.
    pub length_minus_one: usize,
    /// Direct predecessors of the representative.
    pub direct_predecessors: Vec<ArOrbitPosition>,
}

impl ArOrbit {
    /// Construct one compact orbit description.
    #[must_use]
    pub fn new(length_minus_one: usize, direct_predecessors: Vec<ArOrbitPosition>) -> Self {
        Self {
            length_minus_one,
            direct_predecessors,
        }
    }
}

/// Either accepted encoding of QPA's numerical AR-quiver data.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ArQuiverDescription {
    /// Meshes indexed by indecomposable module.
    Meshes(Vec<ArMesh>),
    /// Compact translate-orbit descriptions.
    Orbits(Vec<ArOrbit>),
}

/// A module argument accepted by QPA's numerical degeneration operations.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum NumericalModule {
    /// A zero-based indecomposable-module index.
    Indecomposable(usize),
    /// A full multiplicity vector in indecomposable order.
    Multiplicities(NumericalVector),
}

impl From<usize> for NumericalModule {
    fn from(value: usize) -> Self {
        Self::Indecomposable(value)
    }
}

impl From<NumericalVector> for NumericalModule {
    fn from(value: NumericalVector) -> Self {
        Self::Multiplicities(value)
    }
}

/// The two input forms accepted specifically by QPA's `ModulesOfDimVect`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DimensionVectorInput {
    /// A zero-based indecomposable-module index whose dimension is used.
    Indecomposable(usize),
    /// A dimension vector of length `NumberOfProjectives`.
    DimensionVector(NumericalVector),
}

impl From<usize> for DimensionVectorInput {
    fn from(value: usize) -> Self {
        Self::Indecomposable(value)
    }
}

impl From<NumericalVector> for DimensionVectorInput {
    fn from(value: NumericalVector) -> Self {
        Self::DimensionVector(value)
    }
}

/// Failures corresponding to QPA's explicit errors or implicit GAP indexing failures.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum DegenerationError {
    /// QPA requires positive numbers of indecomposables and projectives.
    #[error(
        "numbers of indecomposables and projectives must be positive, with projectives no greater than indecomposables"
    )]
    InvalidSize,
    /// The compact or explicit description addressed a nonexistent position.
    #[error("the numerical AR-quiver description contains an out-of-range position")]
    DescriptionIndex,
    /// QPA requires precisely the first declared modules to be projective.
    #[error("indecomposables before the projective boundary must be projective")]
    ProjectivesNotFirst,
    /// QPA's Auslander--Reiten matrix was not invertible.
    #[error("the Auslander-Reiten matrix is singular")]
    SingularArMatrix,
    /// The inverse-transpose did not have integral entries.
    #[error("the inverse-transpose Auslander-Reiten matrix is not integral")]
    NonIntegralHomMatrix,
    /// QPA explicitly requires a nonnegative Hom-dimension matrix.
    #[error("the inverse-transpose Auslander-Reiten matrix is not nonnegative")]
    NegativeHomDimension,
    /// QPA's projective rows did not identify exactly one simple per projective.
    #[error("cannot find the simple modules from the projective Hom-dimension rows")]
    CannotFindSimples,
    /// A module argument was neither an in-range indecomposable nor a full vector.
    #[error("a module must be an in-range indecomposable or a full multiplicity vector")]
    InvalidModule,
    /// `ModulesOfDimVect` received neither an indecomposable nor a dimension vector.
    #[error("the input must be an in-range indecomposable or a projective-length dimension vector")]
    InvalidDimensionVectorInput,
    /// QPA's recursive enumeration requires at least one nonsimple column.
    #[error("the QPA enumeration has no nonsimple column")]
    NoNonsimpleColumn,
    /// QPA's `HowManyTimes` helper received a zero vector.
    #[error("HowManyTimes received a zero vector")]
    ZeroEnumerationColumn,
    /// The requested module was unexpectedly absent from its dimension-vector solutions.
    #[error("the requested module is absent from the generated solution set")]
    ModuleAbsentFromSolutions,
}

/// The result data accumulated by QPA's `TestCodimConjecture` loop.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CodimensionConjectureReport {
    /// Pairs `(predecessor, indecomposable)` violating codimension one.
    pub counterexamples: Vec<(NumericalVector, usize)>,
    /// Total degeneration predecessors visited before termination.
    pub predecessor_count: usize,
    /// Total direct predecessors visited before termination.
    pub direct_predecessor_count: usize,
    /// Whether `break_on_first` stopped the QPA loop at its first counterexample.
    pub stopped_early: bool,
}

/// QPA's numerical finite-type Auslander--Reiten quiver object.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArQuiverNumerical {
    number_of_indecomposables: usize,
    number_of_projectives: usize,
    description: ArQuiverDescription,
    dim_hom_matrix: NumericalMatrix,
    simples: Vec<usize>,
}

impl ArQuiverNumerical {
    /// Construct the object by QPA's AR-matrix inverse-transpose procedure.
    ///
    /// Rust indices are zero-based; QPA's source data are one-based.
    ///
    /// # Errors
    ///
    /// Rejects malformed descriptions, singular AR matrices, nonintegral or
    /// negative Hom dimensions, and data from which QPA cannot locate simples.
    pub fn new(
        number_of_indecomposables: usize,
        number_of_projectives: usize,
        description: ArQuiverDescription,
    ) -> Result<Self, DegenerationError> {
        if number_of_indecomposables == 0
            || number_of_projectives == 0
            || number_of_projectives > number_of_indecomposables
        {
            return Err(DegenerationError::InvalidSize);
        }

        let ar_matrix = match &description {
            ArQuiverDescription::Meshes(meshes) => {
                build_mesh_ar_matrix(number_of_indecomposables, number_of_projectives, meshes)?
            }
            ArQuiverDescription::Orbits(orbits) => {
                build_orbit_ar_matrix(number_of_indecomposables, number_of_projectives, orbits)?
            }
        };
        let dim_hom_matrix = inverse_transpose_integral(&ar_matrix)?;
        let mut simples = Vec::new();
        for projective in 0..number_of_projectives {
            for column in 0..number_of_indecomposables {
                if dim_hom_matrix[projective][column] == 1
                    && (0..number_of_projectives)
                        .all(|other| other == projective || dim_hom_matrix[other][column] == 0)
                {
                    simples.push(column);
                }
            }
        }
        if simples.len() != number_of_projectives {
            return Err(DegenerationError::CannotFindSimples);
        }

        Ok(Self {
            number_of_indecomposables,
            number_of_projectives,
            description,
            dim_hom_matrix,
            simples,
        })
    }

    /// Return QPA's `NumberOfIndecomposables` attribute.
    #[must_use]
    pub const fn number_of_indecomposables(&self) -> usize {
        self.number_of_indecomposables
    }

    /// Return QPA's `NumberOfProjectives` attribute.
    #[must_use]
    pub const fn number_of_projectives(&self) -> usize {
        self.number_of_projectives
    }

    /// Return the raw mesh or orbit description stored by QPA.
    #[must_use]
    pub const fn description(&self) -> &ArQuiverDescription {
        &self.description
    }

    /// Return QPA's matrix `[dim Hom(i,j)]`.
    #[must_use]
    pub fn dim_hom_matrix(&self) -> &[NumericalVector] {
        &self.dim_hom_matrix
    }

    /// Return the zero-based positions of the simple modules found by QPA.
    #[must_use]
    pub fn simples(&self) -> &[usize] {
        &self.simples
    }

    /// Translate QPA's `DimensionVector` dispatch.
    ///
    /// # Errors
    ///
    /// Rejects an out-of-range indecomposable or wrong-length vector.
    pub fn dimension_vector(
        &self,
        module: &NumericalModule,
    ) -> Result<NumericalVector, DegenerationError> {
        self.validate_module(module)?;
        match module {
            NumericalModule::Indecomposable(column) => Ok((0..self.number_of_projectives)
                .map(|row| self.dim_hom_matrix[row][*column].clone())
                .collect()),
            NumericalModule::Multiplicities(vector) => Ok(self.dim_hom_matrix
                [..self.number_of_projectives]
                .iter()
                .map(|row| dot(row, vector))
                .collect()),
        }
    }

    /// Translate QPA's `DimHom` four-way dispatch.
    ///
    /// # Errors
    ///
    /// Rejects out-of-range indecomposables or wrong-length vectors.
    pub fn dim_hom(
        &self,
        source: &NumericalModule,
        target: &NumericalModule,
    ) -> Result<Integer, DegenerationError> {
        self.validate_module(source)?;
        self.validate_module(target)?;
        Ok(match (source, target) {
            (NumericalModule::Indecomposable(left), NumericalModule::Indecomposable(right)) => {
                self.dim_hom_matrix[*left][*right].clone()
            }
            (NumericalModule::Indecomposable(left), NumericalModule::Multiplicities(right)) => {
                dot(&self.dim_hom_matrix[*left], right)
            }
            (NumericalModule::Multiplicities(left), NumericalModule::Indecomposable(right)) => (0
                ..self.number_of_indecomposables)
                .map(|index| &left[index] * &self.dim_hom_matrix[index][*right])
                .fold(Integer::from(0), |sum, value| sum + value),
            (NumericalModule::Multiplicities(left), NumericalModule::Multiplicities(right)) => {
                let mut result = Integer::from(0);
                for (row_index, multiplicity) in left.iter().enumerate() {
                    result += multiplicity * dot(&self.dim_hom_matrix[row_index], right);
                }
                result
            }
        })
    }

    /// Translate QPA's `DimEnd`.
    ///
    /// # Errors
    ///
    /// Propagates module validation failures.
    pub fn dim_end(&self, module: &NumericalModule) -> Result<Integer, DegenerationError> {
        self.dim_hom(module, module)
    }

    /// Translate QPA's orbit-dimension formula.
    ///
    /// # Errors
    ///
    /// Propagates module validation failures.
    pub fn orbit_dimension(&self, module: &NumericalModule) -> Result<Integer, DegenerationError> {
        let dimension = self.dimension_vector(module)?;
        let squares = dimension
            .iter()
            .map(|value| value * value)
            .fold(Integer::from(0), |sum, value| sum + value);
        Ok(squares - self.dim_end(module)?)
    }

    /// Translate QPA's unchecked orbit-codimension formula.
    ///
    /// # Errors
    ///
    /// Propagates module validation failures.
    pub fn orbit_codimension(
        &self,
        source: &NumericalModule,
        target: &NumericalModule,
    ) -> Result<Integer, DegenerationError> {
        Ok(self.dim_end(target)? - self.dim_end(source)?)
    }

    /// Translate QPA's checked `DegOrderLEQ`.
    ///
    /// Unlike GAP's method, a dimension-vector mismatch is returned as `false`
    /// without printing a warning.
    ///
    /// # Errors
    ///
    /// Propagates module validation failures.
    pub fn degeneration_order_leq(
        &self,
        source: &NumericalModule,
        target: &NumericalModule,
    ) -> Result<bool, DegenerationError> {
        if self.dimension_vector(source)? != self.dimension_vector(target)? {
            return Ok(false);
        }
        self.degeneration_order_leq_unchecked(source, target)
    }

    /// Translate QPA's `DegOrderLEQNC`, which omits the dimension-vector check.
    ///
    /// # Errors
    ///
    /// Propagates module validation failures.
    pub fn degeneration_order_leq_unchecked(
        &self,
        source: &NumericalModule,
        target: &NumericalModule,
    ) -> Result<bool, DegenerationError> {
        self.validate_module(source)?;
        self.validate_module(target)?;
        for row in self.number_of_projectives..self.number_of_indecomposables {
            if self.row_hom_into(row, source) > self.row_hom_into(row, target) {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// Translate QPA's recursive `ModulesOfDimVect` generator.
    ///
    /// # Errors
    ///
    /// Propagates module validation and QPA enumeration failures.
    pub fn modules_of_dimension_vector(
        &self,
        which: &DimensionVectorInput,
    ) -> Result<Vec<NumericalVector>, DegenerationError> {
        let dimension = match which {
            DimensionVectorInput::Indecomposable(index)
                if *index < self.number_of_indecomposables =>
            {
                self.dimension_vector(&NumericalModule::Indecomposable(*index))?
            }
            DimensionVectorInput::DimensionVector(vector)
                if vector.len() == self.number_of_projectives =>
            {
                vector.clone()
            }
            _ => return Err(DegenerationError::InvalidDimensionVectorInput),
        };
        self.generate_modules(&dimension, None)
    }

    /// Translate QPA's `DegOrderPredecessors`.
    ///
    /// # Errors
    ///
    /// Propagates module validation and QPA enumeration failures.
    pub fn degeneration_predecessors(
        &self,
        which: &NumericalModule,
    ) -> Result<Vec<NumericalVector>, DegenerationError> {
        self.ordered_neighbours(which, OrderFilter::Predecessor)
    }

    /// Translate QPA's naive direct-predecessor filtering.
    ///
    /// # Errors
    ///
    /// Propagates predecessor-generation failures.
    pub fn direct_degeneration_predecessors(
        &self,
        which: &NumericalModule,
    ) -> Result<Vec<NumericalVector>, DegenerationError> {
        let predecessors = self.degeneration_predecessors(which)?;
        Ok(self.direct_predecessors_from(&predecessors))
    }

    /// Generate predecessors once and return QPA's pair of all and direct ones.
    ///
    /// # Errors
    ///
    /// Propagates predecessor-generation failures.
    pub fn degeneration_predecessors_with_direct(
        &self,
        which: &NumericalModule,
    ) -> Result<(Vec<NumericalVector>, Vec<NumericalVector>), DegenerationError> {
        let predecessors = self.degeneration_predecessors(which)?;
        let direct = self.direct_predecessors_from(&predecessors);
        Ok((predecessors, direct))
    }

    /// Translate QPA's `DegOrderSuccessors`.
    ///
    /// # Errors
    ///
    /// Propagates module validation and QPA enumeration failures.
    pub fn degeneration_successors(
        &self,
        which: &NumericalModule,
    ) -> Result<Vec<NumericalVector>, DegenerationError> {
        self.ordered_neighbours(which, OrderFilter::Successor)
    }

    /// Translate QPA's naive direct-successor filtering.
    ///
    /// # Errors
    ///
    /// Propagates successor-generation failures.
    pub fn direct_degeneration_successors(
        &self,
        which: &NumericalModule,
    ) -> Result<Vec<NumericalVector>, DegenerationError> {
        let successors = self.degeneration_successors(which)?;
        Ok(self.direct_successors_from(&successors))
    }

    /// Generate successors once and return QPA's pair of all and direct ones.
    ///
    /// # Errors
    ///
    /// Propagates successor-generation failures.
    pub fn degeneration_successors_with_direct(
        &self,
        which: &NumericalModule,
    ) -> Result<(Vec<NumericalVector>, Vec<NumericalVector>), DegenerationError> {
        let successors = self.degeneration_successors(which)?;
        let direct = self.direct_successors_from(&successors);
        Ok((successors, direct))
    }

    /// Run QPA's `TestCodimConjecture` case loop and return its data.
    ///
    /// Rust returns structured counters instead of printing progress and a
    /// runtime-dependent summary.  `break_on_first` preserves QPA's early exit.
    ///
    /// # Errors
    ///
    /// Propagates predecessor-generation failures.
    pub fn test_codimension_conjecture(
        &self,
        _print_progress: bool,
        break_on_first: bool,
    ) -> Result<CodimensionConjectureReport, DegenerationError> {
        let mut counterexamples = Vec::new();
        let mut predecessor_count = 0_usize;
        let mut direct_predecessor_count = 0_usize;
        for indecomposable in 0..self.number_of_indecomposables {
            let target = NumericalModule::Indecomposable(indecomposable);
            let (predecessors, direct) = self.degeneration_predecessors_with_direct(&target)?;
            predecessor_count += predecessors.len();
            direct_predecessor_count += direct.len();
            for vector in direct {
                if self
                    .orbit_codimension(&NumericalModule::Multiplicities(vector.clone()), &target)?
                    != 1
                {
                    counterexamples.push((vector, indecomposable));
                    if break_on_first {
                        return Ok(CodimensionConjectureReport {
                            counterexamples,
                            predecessor_count,
                            direct_predecessor_count,
                            stopped_early: true,
                        });
                    }
                }
            }
        }
        Ok(CodimensionConjectureReport {
            counterexamples,
            predecessor_count,
            direct_predecessor_count,
            stopped_early: false,
        })
    }

    fn ordered_neighbours(
        &self,
        which: &NumericalModule,
        filter: OrderFilter,
    ) -> Result<Vec<NumericalVector>, DegenerationError> {
        let dimension = self.dimension_vector(which)?;
        let comparison = (self.number_of_projectives..self.number_of_indecomposables)
            .map(|row| self.row_hom_into(row, which))
            .collect::<Vec<_>>();
        let mut solutions = self.generate_modules(&dimension, Some((filter, &comparison)))?;
        let target = self.multiplicity_vector(which)?;
        let position = solutions
            .iter()
            .position(|solution| solution == &target)
            .ok_or(DegenerationError::ModuleAbsentFromSolutions)?;
        solutions.remove(position);
        Ok(solutions)
    }

    fn generate_modules(
        &self,
        dimension: &[Integer],
        filter: Option<(OrderFilter, &[Integer])>,
    ) -> Result<Vec<NumericalVector>, DegenerationError> {
        let nonsimples = (0..self.number_of_indecomposables)
            .filter(|index| !self.simples.contains(index))
            .collect::<Vec<_>>();
        if nonsimples.is_empty() {
            return Err(DegenerationError::NoNonsimpleColumn);
        }
        let columns = (0..self.number_of_indecomposables)
            .map(|column| {
                (0..self.number_of_projectives)
                    .map(|row| self.dim_hom_matrix[row][column].clone())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut state = GenerationState {
            multiplicities: vec![Integer::from(0); self.number_of_indecomposables],
            collected_dimension: vec![Integer::from(0); self.number_of_projectives],
            solutions: Vec::new(),
        };
        generate_recursive(
            self,
            dimension,
            &columns,
            &nonsimples,
            0,
            filter,
            &mut state,
        )?;
        Ok(state.solutions)
    }

    fn multiplicity_vector(
        &self,
        module: &NumericalModule,
    ) -> Result<NumericalVector, DegenerationError> {
        self.validate_module(module)?;
        Ok(match module {
            NumericalModule::Indecomposable(index) => {
                let mut vector = vec![Integer::from(0); self.number_of_indecomposables];
                vector[*index] = 1.into();
                vector
            }
            NumericalModule::Multiplicities(vector) => vector.clone(),
        })
    }

    fn validate_module(&self, module: &NumericalModule) -> Result<(), DegenerationError> {
        match module {
            NumericalModule::Indecomposable(index) if *index < self.number_of_indecomposables => {
                Ok(())
            }
            NumericalModule::Multiplicities(vector)
                if vector.len() == self.number_of_indecomposables =>
            {
                Ok(())
            }
            _ => Err(DegenerationError::InvalidModule),
        }
    }

    fn row_hom_into(&self, row: usize, module: &NumericalModule) -> Integer {
        match module {
            NumericalModule::Indecomposable(column) => self.dim_hom_matrix[row][*column].clone(),
            NumericalModule::Multiplicities(vector) => dot(&self.dim_hom_matrix[row], vector),
        }
    }

    fn vector_leq(&self, source: &[Integer], target: &[Integer]) -> bool {
        (self.number_of_projectives..self.number_of_indecomposables).all(|row| {
            dot(&self.dim_hom_matrix[row], source) <= dot(&self.dim_hom_matrix[row], target)
        })
    }

    fn direct_predecessors_from(&self, predecessors: &[NumericalVector]) -> Vec<NumericalVector> {
        predecessors
            .iter()
            .enumerate()
            .filter(|(left_index, left)| {
                !predecessors.iter().enumerate().any(|(right_index, right)| {
                    left_index != &right_index && self.vector_leq(left, right)
                })
            })
            .map(|(_, vector)| vector.clone())
            .collect()
    }

    fn direct_successors_from(&self, successors: &[NumericalVector]) -> Vec<NumericalVector> {
        successors
            .iter()
            .enumerate()
            .filter(|(right_index, right)| {
                !successors.iter().enumerate().any(|(left_index, left)| {
                    right_index != &left_index && self.vector_leq(left, right)
                })
            })
            .map(|(_, vector)| vector.clone())
            .collect()
    }
}

impl fmt::Display for ArQuiverNumerical {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "<ARQuiverNumerical with {} indecomposables and {} projectives>",
            self.number_of_indecomposables, self.number_of_projectives
        )
    }
}

/// Return QPA's readable rendering of one multiplicity vector, including its newline.
#[must_use]
pub fn format_multiplicity_vector(vector: &[Integer]) -> String {
    let mut rendered = String::new();
    for (index, multiplicity) in vector.iter().enumerate() {
        if *multiplicity != 0 {
            if !rendered.is_empty() {
                rendered.push_str(" + ");
            }
            write!(&mut rendered, "{multiplicity}*({})", index + 1)
                .expect("writing to a String cannot fail");
        }
    }
    if rendered.is_empty() {
        rendered.push('0');
    }
    rendered.push('\n');
    rendered
}

/// Return QPA's consecutive readable renderings of multiplicity vectors.
#[must_use]
pub fn format_multiplicity_vectors(vectors: &[NumericalVector]) -> String {
    vectors
        .iter()
        .map(|vector| format_multiplicity_vector(vector))
        .collect()
}

#[derive(Clone, Copy)]
enum OrderFilter {
    Predecessor,
    Successor,
}

struct GenerationState {
    multiplicities: NumericalVector,
    collected_dimension: NumericalVector,
    solutions: Vec<NumericalVector>,
}

fn generate_recursive(
    ar: &ArQuiverNumerical,
    dimension: &[Integer],
    columns: &[NumericalVector],
    nonsimples: &[usize],
    nonsimple_position: usize,
    filter: Option<(OrderFilter, &[Integer])>,
    state: &mut GenerationState,
) -> Result<(), DegenerationError> {
    let current = nonsimples[nonsimple_position];
    let column = &columns[current];
    let remaining = dimension
        .iter()
        .zip(&state.collected_dimension)
        .map(|(total, collected)| total - collected)
        .collect::<Vec<_>>();
    let how_many = how_many_times(column, &remaining)?;
    let mut multiplicity = Integer::from(0);
    while multiplicity <= how_many {
        state.multiplicities[current] = multiplicity.clone();
        if nonsimple_position + 1 == nonsimples.len() {
            let mut solution = state.multiplicities.clone();
            for projective in 0..ar.number_of_projectives {
                let used = dot(&ar.dim_hom_matrix[projective], &solution);
                solution[ar.simples[projective]] = &dimension[projective] - used;
            }
            let retained = match filter {
                None => true,
                Some((order, comparison)) => (ar.number_of_projectives
                    ..ar.number_of_indecomposables)
                    .enumerate()
                    .all(|(offset, row)| {
                        let value = dot(&ar.dim_hom_matrix[row], &solution);
                        match order {
                            OrderFilter::Predecessor => value <= comparison[offset],
                            OrderFilter::Successor => value >= comparison[offset],
                        }
                    }),
            };
            if retained {
                state.solutions.push(solution);
            }
        } else {
            let previous = state.collected_dimension.clone();
            for (entry, value) in state.collected_dimension.iter_mut().zip(column) {
                *entry += &multiplicity * value;
            }
            generate_recursive(
                ar,
                dimension,
                columns,
                nonsimples,
                nonsimple_position + 1,
                filter,
                state,
            )?;
            state.collected_dimension = previous;
        }
        multiplicity += 1;
    }
    Ok(())
}

fn how_many_times(column: &[Integer], remaining: &[Integer]) -> Result<Integer, DegenerationError> {
    let first = column
        .iter()
        .position(|entry| *entry != 0)
        .ok_or(DegenerationError::ZeroEnumerationColumn)?;
    let mut count = &remaining[first] / &column[first];
    for index in first + 1..remaining.len() {
        if column[index] != 0 {
            let quotient = &remaining[index] / &column[index];
            if quotient < count {
                count = quotient;
            }
        }
    }
    Ok(count)
}

fn dot(left: &[Integer], right: &[Integer]) -> Integer {
    left.iter()
        .zip(right)
        .map(|(a, b)| a * b)
        .fold(Integer::from(0), |sum, value| sum + value)
}

fn build_mesh_ar_matrix(
    n: usize,
    projectives: usize,
    meshes: &[ArMesh],
) -> Result<Vec<Vec<i64>>, DegenerationError> {
    if meshes.len() < projectives || meshes.len() > n {
        return Err(DegenerationError::DescriptionIndex);
    }
    if meshes[..projectives].iter().any(|mesh| mesh.tau.is_some()) {
        return Err(DegenerationError::ProjectivesNotFirst);
    }
    let mut matrix = identity_i64(n);
    for (row, mesh) in meshes.iter().enumerate() {
        for &predecessor in &mesh.direct_predecessors {
            if predecessor >= n {
                return Err(DegenerationError::DescriptionIndex);
            }
            matrix[row][predecessor] = -1;
        }
        if let Some(tau) = mesh.tau {
            if tau >= n {
                return Err(DegenerationError::DescriptionIndex);
            }
            matrix[row][tau] += 1;
        }
    }
    Ok(matrix)
}

fn build_orbit_ar_matrix(
    n: usize,
    projectives: usize,
    orbits: &[ArOrbit],
) -> Result<Vec<Vec<i64>>, DegenerationError> {
    if orbits.len() < projectives || orbits.len() > n {
        return Err(DegenerationError::DescriptionIndex);
    }
    let representative_count = orbits.len();
    let mut next_extra = representative_count;
    let mut first_extra = Vec::with_capacity(representative_count);
    for orbit in orbits {
        if orbit.length_minus_one == 0 {
            first_extra.push(None);
        } else {
            first_extra.push(Some(next_extra));
            next_extra = next_extra
                .checked_add(orbit.length_minus_one)
                .ok_or(DegenerationError::DescriptionIndex)?;
        }
    }
    if next_extra > n {
        return Err(DegenerationError::DescriptionIndex);
    }
    let mut matrix = identity_i64(n);
    for (orbit_index, orbit) in orbits.iter().enumerate() {
        if orbit_index >= projectives {
            if orbit.length_minus_one == 0 {
                matrix[orbit_index][orbit_index] = 2;
            } else {
                let last = orbit_position(
                    orbits,
                    &first_extra,
                    n,
                    ArOrbitPosition::new(orbit_index, orbit.length_minus_one),
                )?;
                matrix[orbit_index][last] = 1;
            }
        }
        if let Some(first) = first_extra[orbit_index] {
            matrix[first][orbit_index] = 1;
            for (row, entries) in matrix
                .iter_mut()
                .enumerate()
                .skip(first + 1)
                .take(orbit.length_minus_one - 1)
            {
                entries[row - 1] = 1;
            }
        }
    }

    add_orbit_incoming_arrows(&mut matrix, orbits, &first_extra, n, projectives)?;
    add_projective_radical_arrows(&mut matrix, orbits, &first_extra, n, projectives)?;
    Ok(matrix)
}

fn add_orbit_incoming_arrows(
    matrix: &mut [Vec<i64>],
    orbits: &[ArOrbit],
    first_extra: &[Option<usize>],
    n: usize,
    projectives: usize,
) -> Result<(), DegenerationError> {
    let position = |value| orbit_position(orbits, first_extra, n, value);

    for (orbit_index, orbit) in orbits.iter().enumerate() {
        for &predecessor in &orbit.direct_predecessors {
            matrix[orbit_index][position(predecessor)?] = -1;
        }
        for shift in 1..=orbit.length_minus_one {
            let row = position(ArOrbitPosition::new(orbit_index, shift))?;
            for &predecessor in &orbit.direct_predecessors {
                let predecessor_orbit = orbits
                    .get(predecessor.orbit)
                    .ok_or(DegenerationError::DescriptionIndex)?;
                if predecessor.orbit < projectives
                    && shift + predecessor.shift <= predecessor_orbit.length_minus_one
                {
                    let column = position(ArOrbitPosition::new(
                        predecessor.orbit,
                        predecessor.shift + shift,
                    ))?;
                    matrix[row][column] = -1;
                }
                if predecessor.orbit >= projectives {
                    let translated =
                        (predecessor.shift + shift) % (predecessor_orbit.length_minus_one + 1);
                    let column = position(ArOrbitPosition::new(predecessor.orbit, translated))?;
                    matrix[row][column] = -1;
                }
            }
        }
    }
    Ok(())
}

fn add_projective_radical_arrows(
    matrix: &mut [Vec<i64>],
    orbits: &[ArOrbit],
    first_extra: &[Option<usize>],
    n: usize,
    projectives: usize,
) -> Result<(), DegenerationError> {
    let position = |value| orbit_position(orbits, first_extra, n, value);
    for projective in 0..projectives {
        for &radical_summand in &orbits[projective].direct_predecessors {
            let summand_orbit = orbits
                .get(radical_summand.orbit)
                .ok_or(DegenerationError::DescriptionIndex)?;
            for shift in 0..=orbits[projective].length_minus_one {
                let column = position(ArOrbitPosition::new(projective, shift))?;
                if radical_summand.orbit < projectives {
                    let translated = radical_summand.shift + 1 + shift;
                    if translated <= summand_orbit.length_minus_one {
                        let row =
                            position(ArOrbitPosition::new(radical_summand.orbit, translated))?;
                        matrix[row][column] = -1;
                    }
                } else {
                    let translated =
                        (radical_summand.shift + 1 + shift) % (summand_orbit.length_minus_one + 1);
                    let row = position(ArOrbitPosition::new(radical_summand.orbit, translated))?;
                    matrix[row][column] = -1;
                }
            }
        }
    }
    Ok(())
}

fn orbit_position(
    orbits: &[ArOrbit],
    first_extra: &[Option<usize>],
    n: usize,
    position: ArOrbitPosition,
) -> Result<usize, DegenerationError> {
    let orbit = orbits
        .get(position.orbit)
        .ok_or(DegenerationError::DescriptionIndex)?;
    if position.shift == 0 {
        return Ok(position.orbit);
    }
    let first = first_extra[position.orbit].ok_or(DegenerationError::DescriptionIndex)?;
    let result = first
        .checked_add(position.shift - 1)
        .ok_or(DegenerationError::DescriptionIndex)?;
    if result >= n || position.shift > orbit.length_minus_one {
        return Err(DegenerationError::DescriptionIndex);
    }
    Ok(result)
}

fn identity_i64(size: usize) -> Vec<Vec<i64>> {
    (0..size)
        .map(|row| (0..size).map(|column| i64::from(row == column)).collect())
        .collect()
}

fn inverse_transpose_integral(matrix: &[Vec<i64>]) -> Result<NumericalMatrix, DegenerationError> {
    let rational = Matrix::from_nested_vec(
        matrix
            .iter()
            .map(|row| row.iter().copied().map(Into::into).collect())
            .collect(),
        Q,
    )
    .map_err(|_| DegenerationError::DescriptionIndex)?;
    let inverse_transpose = rational
        .inv()
        .map_err(|_| DegenerationError::SingularArMatrix)?
        .into_transposed();
    inverse_transpose
        .row_iter()
        .map(|row| {
            row.iter()
                .map(|entry| {
                    if !entry.is_integer() {
                        Err(DegenerationError::NonIntegralHomMatrix)
                    } else if entry < &0 {
                        Err(DegenerationError::NegativeHomDimension)
                    } else {
                        Ok(entry.numerator())
                    }
                })
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use numerica::domains::integer::Integer;
    use serde_json::{Value, json};

    use super::{
        ArMesh, ArOrbit, ArOrbitPosition, ArQuiverDescription, ArQuiverNumerical,
        DimensionVectorInput, NumericalModule, NumericalVector, format_multiplicity_vector,
        format_multiplicity_vectors,
    };

    fn integers(values: &[i64]) -> Vec<numerica::domains::integer::Integer> {
        values.iter().copied().map(Into::into).collect()
    }

    fn a3_zero() -> ArQuiverNumerical {
        ArQuiverNumerical::new(
            5,
            3,
            ArQuiverDescription::Meshes(vec![
                ArMesh::new(vec![], None),
                ArMesh::new(vec![0], None),
                ArMesh::new(vec![3], None),
                ArMesh::new(vec![1], Some(0)),
                ArMesh::new(vec![2], Some(3)),
            ]),
        )
        .unwrap()
    }

    fn d4_subspace() -> ArQuiverNumerical {
        ArQuiverNumerical::new(
            12,
            4,
            ArQuiverDescription::Meshes(vec![
                ArMesh::new(vec![], None),
                ArMesh::new(vec![0], None),
                ArMesh::new(vec![0], None),
                ArMesh::new(vec![0], None),
                ArMesh::new(vec![1, 2, 3], Some(0)),
                ArMesh::new(vec![4], Some(1)),
                ArMesh::new(vec![4], Some(2)),
                ArMesh::new(vec![4], Some(3)),
                ArMesh::new(vec![5, 6, 7], Some(4)),
                ArMesh::new(vec![8], Some(5)),
                ArMesh::new(vec![8], Some(6)),
                ArMesh::new(vec![8], Some(7)),
            ]),
        )
        .unwrap()
    }

    fn orbit_example() -> ArQuiverNumerical {
        ArQuiverNumerical::new(
            9,
            2,
            ArQuiverDescription::Orbits(vec![
                ArOrbit::new(0, vec![ArOrbitPosition::new(3, 0)]),
                ArOrbit::new(1, vec![ArOrbitPosition::new(2, 2)]),
                ArOrbit::new(
                    3,
                    vec![ArOrbitPosition::new(3, 0), ArOrbitPosition::new(1, 1)],
                ),
                ArOrbit::new(
                    1,
                    vec![ArOrbitPosition::new(2, 1), ArOrbitPosition::new(2, 3)],
                ),
            ]),
        )
        .unwrap()
    }

    fn machine_integer(value: &Integer) -> i64 {
        value.to_i64().unwrap()
    }

    fn machine_vector(vector: &[Integer]) -> Vec<i64> {
        vector.iter().map(machine_integer).collect()
    }

    fn machine_vectors(vectors: &[NumericalVector]) -> Vec<Vec<i64>> {
        vectors
            .iter()
            .map(|vector| machine_vector(vector))
            .collect()
    }

    fn fixture_mesh_case() -> Value {
        let ar = d4_subspace();
        let indecomposable = NumericalModule::Indecomposable(4);
        let decomposable =
            NumericalModule::Multiplicities(integers(&[0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0]));
        let dimension = ar.dimension_vector(&indecomposable).unwrap();
        let modules_number = ar
            .modules_of_dimension_vector(&DimensionVectorInput::Indecomposable(4))
            .unwrap();
        let modules_vector = ar
            .modules_of_dimension_vector(&DimensionVectorInput::DimensionVector(dimension))
            .unwrap();
        let predecessors = ar.degeneration_predecessors(&decomposable).unwrap();
        let direct_predecessors = ar.direct_degeneration_predecessors(&decomposable).unwrap();
        let predecessors_with_direct = ar
            .degeneration_predecessors_with_direct(&decomposable)
            .unwrap();
        let successors = ar.degeneration_successors(&indecomposable).unwrap();
        let direct_successors = ar.direct_degeneration_successors(&indecomposable).unwrap();
        let successors_with_direct = ar
            .degeneration_successors_with_direct(&indecomposable)
            .unwrap();
        json!({
            "number_of_indecomposables": ar.number_of_indecomposables(),
            "number_of_projectives": ar.number_of_projectives(),
            "dim_hom_matrix": machine_vectors(ar.dim_hom_matrix()),
            "simples": ar.simples().iter().map(|index| index + 1).collect::<Vec<_>>(),
            "display": ar.to_string(),
            "dimension_vectors": (0..12).map(|index| {
                machine_vector(&ar.dimension_vector(&NumericalModule::Indecomposable(index)).unwrap())
            }).collect::<Vec<_>>(),
            "dimension_vector_of_vector": machine_vector(&ar.dimension_vector(&decomposable).unwrap()),
            "modules_number": machine_vectors(&modules_number),
            "modules_vector": machine_vectors(&modules_vector),
            "predecessors": machine_vectors(&predecessors),
            "direct_predecessors": machine_vectors(&direct_predecessors),
            "predecessors_with_direct": [
                machine_vectors(&predecessors_with_direct.0),
                machine_vectors(&predecessors_with_direct.1),
            ],
            "successors": machine_vectors(&successors),
            "direct_successors": machine_vectors(&direct_successors),
            "successors_with_direct": [
                machine_vectors(&successors_with_direct.0),
                machine_vectors(&successors_with_direct.1),
            ],
            "dim_hom_dispatch": [
                machine_integer(&ar.dim_hom(&indecomposable, &indecomposable).unwrap()),
                machine_integer(&ar.dim_hom(&indecomposable, &decomposable).unwrap()),
                machine_integer(&ar.dim_hom(&decomposable, &indecomposable).unwrap()),
                machine_integer(&ar.dim_hom(&decomposable, &decomposable).unwrap()),
            ],
            "dim_end": [
                machine_integer(&ar.dim_end(&indecomposable).unwrap()),
                machine_integer(&ar.dim_end(&decomposable).unwrap()),
            ],
            "orbit_dimension": [
                machine_integer(&ar.orbit_dimension(&indecomposable).unwrap()),
                machine_integer(&ar.orbit_dimension(&decomposable).unwrap()),
            ],
            "orbit_codimension": machine_integer(
                &ar.orbit_codimension(&indecomposable, &decomposable).unwrap()
            ),
            "degeneration_leq": [
                ar.degeneration_order_leq(&indecomposable, &decomposable).unwrap(),
                ar.degeneration_order_leq(&decomposable, &indecomposable).unwrap(),
            ],
            "degeneration_leq_nc": [
                ar.degeneration_order_leq_unchecked(&indecomposable, &decomposable).unwrap(),
                ar.degeneration_order_leq_unchecked(&decomposable, &indecomposable).unwrap(),
            ],
            "different_dimension_leq": ar.degeneration_order_leq(
                &NumericalModule::Indecomposable(0),
                &NumericalModule::Indecomposable(1),
            ).unwrap(),
        })
    }

    #[test]
    fn qpa_mesh_constructor_and_numerical_operations() {
        let ar = a3_zero();
        assert_eq!(ar.number_of_indecomposables(), 5);
        assert_eq!(ar.number_of_projectives(), 3);
        assert_eq!(ar.simples(), &[0, 3, 4]);
        assert_eq!(
            ar.dim_hom_matrix(),
            [
                integers(&[1, 1, 0, 0, 0]),
                integers(&[0, 1, 1, 1, 0]),
                integers(&[0, 0, 1, 0, 1]),
                integers(&[0, 0, 1, 1, 0]),
                integers(&[0, 0, 0, 0, 1]),
            ]
        );
        assert_eq!(
            ar.dimension_vector(&NumericalModule::Indecomposable(4))
                .unwrap(),
            integers(&[0, 0, 1])
        );
        assert_eq!(
            ar.dim_hom(
                &NumericalModule::Multiplicities(integers(&[1, 0, 0, 1, 0])),
                &NumericalModule::Indecomposable(4),
            )
            .unwrap(),
            0
        );
    }

    #[test]
    fn qpa_print_multiplicity_vector_format() {
        assert_eq!(
            format_multiplicity_vector(&integers(&[0, 2, 0, 1])),
            "2*(2) + 1*(4)\n"
        );
        assert_eq!(format_multiplicity_vector(&integers(&[0, 0])), "0\n");
        assert_eq!(
            format_multiplicity_vectors(&[integers(&[1, 0, 3]), integers(&[0, 4, 0])]),
            "1*(1) + 3*(3)\n4*(2)\n"
        );
        assert_eq!(
            a3_zero().to_string(),
            "<ARQuiverNumerical with 5 indecomposables and 3 projectives>"
        );
    }

    #[test]
    fn pinned_qpa_degorderfinitetype_fixture_matches() {
        let expected: Value = serde_json::from_str(include_str!(
            "../../../tests/fixtures/degorderfinitetype.schema-v1.json"
        ))
        .unwrap();
        assert_eq!(fixture_mesh_case(), expected["mesh"]);

        let orbit = orbit_example();
        assert_eq!(
            json!({
                "dim_hom_matrix": machine_vectors(orbit.dim_hom_matrix()),
                "simples": orbit.simples().iter().map(|index| index + 1).collect::<Vec<_>>(),
            }),
            expected["orbit"]
        );

        let a3 = a3_zero();
        let report = a3.test_codimension_conjecture(false, false).unwrap();
        let counterexamples = report
            .counterexamples
            .iter()
            .map(|(vector, index)| (machine_vector(vector), index + 1))
            .collect::<Vec<_>>();
        let counts = (0..5)
            .map(|index| {
                let (all, direct) = a3
                    .degeneration_predecessors_with_direct(&NumericalModule::Indecomposable(index))
                    .unwrap();
                [all.len(), direct.len()]
            })
            .collect::<Vec<_>>();
        assert_eq!(
            json!({
                "counterexamples": counterexamples,
                "predecessor_direct_counts": counts,
            }),
            expected["conjecture"]
        );
        assert_eq!(report.predecessor_count, 0);
        assert_eq!(report.direct_predecessor_count, 0);
        assert!(!report.stopped_early);
    }
}
