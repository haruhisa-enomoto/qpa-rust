/// Failure while constructing QPA's regular bimodule over its enveloping
/// algebra.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum EnvelopingModuleError {
    /// A finite quotient must expose the canonical basis prepared by GBNP.
    #[error("a finite-dimensional quotient has no canonical basis vectors")]
    MissingFiniteCanonicalBasis,
    /// A path-algebra image contained a path absent from its canonical basis.
    #[error("an element path is absent from the finite canonical basis")]
    ElementOutsideCanonicalBasis,
    /// The algebra and element variants must remain paired.
    #[error("canonical-basis arithmetic returned another quiver-algebra variant")]
    UnexpectedElementVariant,
    /// The coordinate-swap name must resolve in the original product quiver.
    #[error("cannot find swapped product-quiver generator {name:?}")]
    MissingSwappedProductGenerator {
        /// Generated QPA product name.
        name: Box<str>,
    },
    /// QPA defines algebra Loewy length only for finite-dimensional input.
    #[error("the entered algebra is not finite dimensional")]
    InfiniteLoewyAlgebra,
    /// QPA's method requires an admissible quotient when input is not a path algebra.
    #[error("LoewyLength requires an admissible quotient of a path algebra")]
    NonadmissibleLoewyAlgebra,
    /// `SimpleModules` rejects infinite-dimensional input before constructing
    /// vertex representations.
    #[error("SimpleModules requires a finite-dimensional algebra")]
    InfiniteSimpleModuleAlgebra,
    /// The source requests projective-resolution degree `n + 1`.
    #[error("the global-dimension bound {bound} has no representable successor")]
    GlobalDimensionBoundOverflow {
        /// Entered nonnegative bound.
        bound: usize,
    },
    /// The complexity methods proceed only for admissible quotient algebras
    /// after their cached finite-global-dimension shortcut.
    #[error(
        "entered algebra is not finite dimensional or an admissible quotient of a path algebra"
    )]
    ComplexityRequiresAdmissibleQuotient,
    /// QPA's floating logarithmic expression is undefined for a nonzero
    /// projective-resolution object in these degrees.
    #[error("the complexity estimate is undefined at projective-resolution degree {degree}")]
    UndefinedComplexityEstimate {
        /// Entered arbitrary-precision complex degree.
        degree: Integer,
    },
    /// QPA's `Maximum` requires at least one vertex simple.
    #[error("ComplexityOfAlgebra found no simple modules")]
    ComplexityAlgebraHasNoSimpleModules,
    /// Both inputs to an add(M)-approximation must use the identical acting
    /// algebra.
    #[error("the two modules entered into the add(M)-approximation are not over the same algebra")]
    ApproximationModulesOverDifferentAlgebras,
    /// Every module in QPA's list overload must use the target algebra.
    #[error("not all modules in the approximation list are over the target algebra")]
    ListApproximationModulesOverDifferentAlgebras,
    /// The common-algebra check makes QPA's internal Hom dispatch applicable.
    #[error("HomOverAlgebra unexpectedly rejected add(M)-approximation inputs")]
    ApproximationHomUnavailable,
    /// The number of radical composites must remain addressable.
    #[error("the add(M)-approximation radical-composite count overflows usize")]
    ApproximationMapCountOverflow,
    /// A nonempty quotient basis produces a nonempty direct sum in QPA.
    #[error("an internal nonempty add(M)-approximation direct sum returned fail")]
    EmptyApproximationDirectSum,
    /// Direct-sum structure maps correspond one-for-one with quotient lifts.
    #[error("add(M)-approximation structure-map and generator counts differ")]
    ApproximationMapCountMismatch,
    /// The source takes `Sum` only after a nonempty quotient-basis branch.
    #[error("an internal add(M)-approximation map sum was empty")]
    EmptyApproximationGeneratorList,
    /// QPA's nonempty-list branch indexes direct-sum data even if it selected no map.
    #[error("the list add(M)-approximation selected no source map")]
    EmptyListApproximationDirectSum,
    /// A scanned Hom row must occur in the Hom basis it came from.
    #[error("a list add(M)-approximation row disappeared from its Hom basis")]
    MissingListApproximationHomPosition,
    /// QPA accepts only a one-based index into the entered module list.
    #[error("radical approximation index {index} is outside 1..={length}")]
    RadicalApproximationIndexOutOfRange {
        /// Rejected arbitrary-precision GAP integer.
        index: Integer,
        /// Entered list length.
        length: usize,
    },
    /// QPA's radical branch always direct-sums a nonempty source list.
    #[error("the internal radical add(M)-approximation direct sum returned fail")]
    InternalRadicalApproximationDirectSumFailed,
    /// QPA rejects a negative length before constructing its initial radical
    /// approximation.
    #[error("the entered simple-endo projective-resolution length is negative: {length}")]
    NegativeSimpleEndoProjectiveResolutionLength {
        /// Rejected arbitrary-precision GAP integer.
        length: Integer,
    },
    /// QPA receives `fail` here and then attempts to select its second list
    /// component.
    #[error("CommonDirectSummand returned fail inside the simple-endo projective resolution")]
    SimpleEndoResolutionCommonSummandFail,
    /// The projective cover and its add(M)-approximation share their source.
    #[error("the internal LeftFacMApproximation pushout returned fail")]
    InternalLeftFacPushOutFailed,
    /// The injective envelope and its add(M)-approximation share their range.
    #[error("the internal RightSubMApproximation pullback returned fail")]
    InternalRightSubPullBackFailed,
    /// GAP's negative-bound range is empty and `FiniteComplex` indexes it.
    #[error("HaveFiniteCoresolutionInAddM produced no differentials for bound {bound}")]
    EmptyCoresolutionDifferentialList {
        /// Entered negative GAP integer.
        bound: Integer,
    },
    /// A successful bounded projective-dimension call must bind the attribute.
    #[error("TiltingModule found no cached projective dimension after a successful bounded check")]
    MissingTiltingProjectiveDimensionCache,
    /// GAP's negative-bound range is empty and `FiniteComplex` indexes it.
    #[error("HaveFiniteResolutionInAddM produced no differentials for bound {bound}")]
    EmptyResolutionDifferentialList {
        /// Entered negative GAP integer.
        bound: Integer,
    },
    /// A successful bounded injective-dimension call must bind the attribute.
    #[error("CotiltingModule found no cached injective dimension after a successful bounded check")]
    MissingCotiltingInjectiveDimensionCache,
    /// QPA's first cotilting-approximation check calls an unbound attribute.
    #[error("IsCotiltingModule is not installed on the first approximation argument")]
    CotiltingApproximationPropertyUnbound,
    /// QPA explicitly rejects a cached false cotilting property.
    #[error("the first approximation argument is not a cotilting module")]
    FirstApproximationArgumentNotCotilting,
    /// The approximation methods read QPA's non-computing attribute directly.
    #[error("the cotilting approximation has no cached injective dimension")]
    MissingCotiltingApproximationInjectiveDimensionCache,
    /// At injective dimension zero QPA constructs an empty list and indexes it.
    #[error("the cotilting approximation produced no projective-resolution exact sequences")]
    EmptyCotiltingApproximationExactSequences,
    /// QPA assumes its resolution and pushout objects are isomorphic here.
    #[error("an internal cotilting-approximation module isomorphism returned false")]
    InternalCotiltingApproximationIsomorphismFailed,
    /// QPA indexes the first component of each internal pushout directly.
    #[error("an internal cotilting-approximation pushout returned fail")]
    InternalCotiltingApproximationPushOutFailed,
    /// QPA applies `Minimum` to the complete indecomposable-projective list.
    #[error("FaithfulDimension cannot take Minimum of an empty projective list")]
    EmptyFaithfulDimensionProjectiveList,
    /// A faithful-dimension approximation chain must remain addressable.
    #[error("the FaithfulDimension approximation length overflows usize")]
    FaithfulDimensionOverflow,
    /// A finite faithful dimension is incremented exactly once by QPA.
    #[error("the number of complements overflows usize")]
    NumberOfComplementsOverflow,
    /// A generated GAP resolution degree must fit Rust's addressable lazy
    /// cache in its direction from degree zero.
    #[error("resolution degree is not representable by usize: {degree}")]
    ResolutionDegreeNotRepresentable {
        /// Rejected arbitrary-precision degree.
        degree: Integer,
    },
    /// A panic while generating a lazy differential poisons the shared cache.
    #[error("the lazy projective-resolution cache was poisoned")]
    ProjectiveResolutionCachePoisoned,
    /// Every nonzero projective resolution starts with its degree-zero cover.
    #[error("the nonzero projective-resolution differential cache is empty")]
    EmptyProjectiveResolutionCache,
    /// A panic while generating a lazy injective differential poisons the
    /// shared cache.
    #[error("the lazy injective-resolution cache was poisoned")]
    InjectiveResolutionCachePoisoned,
    /// Every nonzero injective resolution starts with its degree-zero
    /// envelope.
    #[error("the nonzero injective-resolution differential cache is empty")]
    EmptyInjectiveResolutionCache,
    /// QPA's algebra dominant-dimension method requires a boolean
    /// selfinjectivity property after its finite-dimensionality check.
    #[error("selfinjectivity is unavailable for this finite dominant-dimension input")]
    DominantSelfinjectivityUnavailable,
    /// Every indecomposable projective of the checked finite algebra must
    /// remain in the module dominant-dimension method domain.
    #[error("a finite-algebra projective has no dominant-dimension method")]
    FiniteDominantModuleUnavailable,
    /// QPA rejects negative arguments to `NthSyzygy` before constructing the
    /// resolution.
    #[error("entered syzygy degree is negative: {degree}")]
    NegativeSyzygyDegree {
        /// Rejected arbitrary-precision degree.
        degree: Integer,
    },
    /// QPA's two-element `IsNthSyzygy` direct sum cannot be empty.
    #[error("the internal IsNthSyzygy direct sum returned fail")]
    InternalNthSyzygyDirectSumFailed,
    /// QPA's sum of indecomposable-projective socle dimension vectors must
    /// fit Rust's addressable module dimensions.
    #[error("the selfinjective socle multiplicity at vertex {vertex} overflows usize")]
    SelfinjectiveDimensionVectorOverflow {
        /// Zero-based quiver vertex position.
        vertex: usize,
    },
    /// QPA explicitly rejects infinite-dimensional input to its schurian
    /// property method.
    #[error("the entered algebra is not finite dimensional")]
    InfiniteSchurianAlgebra,
    /// A finite selfinjective quotient must construct its regular bimodule.
    #[error("the internal symmetric-algebra regular bimodule construction returned fail")]
    InternalSymmetricRegularModuleFailed,
    /// A finite selfinjective quotient must construct its dual regular
    /// bimodule.
    #[error("the internal symmetric-algebra dual bimodule construction returned fail")]
    InternalSymmetricDualModuleFailed,
    /// The ordered vertex-pair table used by QPA's semicommutativity loop
    /// must fit Rust's address space.
    #[error("the semicommutative vertex-pair count for {vertex_count} vertices overflows usize")]
    SemicommutativePairCountOverflow {
        /// Number of quiver vertices.
        vertex_count: usize,
    },
    /// A nonzero exact determinant must make the rational Cartan matrix
    /// invertible in QPA's following source step.
    #[error("a Cartan matrix with nonzero determinant unexpectedly had no rational inverse")]
    NonzeroCartanDeterminantNotInvertible,
    /// QPA's selfinjectivity call has no applicable method on this quotient.
    #[error("NakayamaPermutation cannot dispatch IsSelfinjectiveAlgebra on this algebra")]
    NakayamaPermutationSelfinjectivityUnavailable,
    /// The dedicated cache can contain only an unbound state, boolean false,
    /// or one list of index images.
    #[error("the NakayamaPermutation attribute cache contains an invalid fail state")]
    InvalidNakayamaPermutationCacheState,
    /// Concurrent attribute evaluation cannot bind both false and a value.
    #[error("the NakayamaPermutation attribute was concurrently cached as false")]
    NakayamaPermutationCacheConflict,
    /// The source closure accepts only modules of total dimension one.
    #[error("not a simple module was entered as an argument for the Nakayama permutation")]
    NakayamaPermutationModuleNotSimple,
    /// A foreign simple module can select beyond the captured projective list.
    #[error("Nakayama permutation module position {position} is outside 1..={vertex_count}")]
    NakayamaPermutationModulePositionOutOfRange {
        /// One-based position of the dimension-vector entry equal to one.
        position: usize,
        /// Number of captured projective socles.
        vertex_count: usize,
    },
    /// QPA's index closure accepts only the algebra's one-based simple indices.
    #[error(
        "incorrect Nakayama permutation index {index}; expected an integer in 1..={vertex_count}"
    )]
    IncorrectNakayamaPermutationIndex {
        /// Rejected arbitrary-precision integer.
        index: Integer,
        /// Number of simple-module indices.
        vertex_count: usize,
    },
    /// Selfinjectivity guarantees that every projective socle is a vertex
    /// simple and therefore has one dimension-vector entry equal to one.
    #[error("projective {projective} has no simple socle position in NakayamaPermutation")]
    SelfinjectiveProjectiveSocleHasNoSimplePosition {
        /// One-based indecomposable-projective position.
        projective: usize,
    },
    /// QPA's selfinjectivity call has no applicable method on this quotient.
    #[error("NakayamaAutomorphism cannot dispatch IsSelfinjectiveAlgebra on this algebra")]
    NakayamaAutomorphismSelfinjectivityUnavailable,
    /// The automorphism cache cannot contain GAP's `fail` sentinel.
    #[error("the NakayamaAutomorphism attribute cache contains an invalid fail state")]
    InvalidNakayamaAutomorphismCacheState,
    /// Concurrent evaluation cannot bind both boolean false and a matrix.
    #[error("the NakayamaAutomorphism attribute was concurrently cached as false")]
    NakayamaAutomorphismCacheConflict,
    /// The source closure requires an element in its captured algebra.
    #[error("the entered argument is not in the algebra")]
    NakayamaAutomorphismElementOutsideAlgebra,
    /// A successful finite selfinjective computation requires a canonical basis.
    #[error("NakayamaAutomorphism requires a finite canonical algebra basis")]
    MissingNakayamaCanonicalBasis,
    /// Every selfinjective indecomposable projective has a nonzero socle generator.
    #[error("projective {projective} has no socle generator in NakayamaAutomorphism")]
    EmptyProjectiveSocleGenerator {
        /// One-based indecomposable-projective position.
        projective: usize,
    },
    /// GAP's Steinitz primitive received inconsistent basis shapes or ranks.
    #[error("BaseSteinitzVectors received invalid Nakayama basis data")]
    InvalidNakayamaSteinitzInput,
    /// The source Steinitz extension must produce a full algebra basis.
    #[error("BaseSteinitzVectors did not complete the Nakayama socle basis")]
    IncompleteNakayamaSteinitzBasis,
    /// Addressable Steinitz column positions cannot overflow.
    #[error("the Nakayama Steinitz column counter overflowed usize")]
    NakayamaSteinitzColumnOverflow,
    /// Addressable Steinitz row positions cannot overflow.
    #[error("the Nakayama Steinitz row counter overflowed usize")]
    NakayamaSteinitzRowOverflow,
    /// A selfinjective algebra's Frobenius pairing must be nondegenerate.
    #[error("the Frobenius pairing constructed by NakayamaAutomorphism is singular")]
    SingularNakayamaPairing,
    /// The cached closure matrix must match the canonical-basis dimension.
    #[error("the cached Nakayama matrix has shape {shape:?}, expected {dimension} x {dimension}")]
    InvalidNakayamaAutomorphismMatrixShape {
        /// Cached matrix shape.
        shape: (usize, usize),
        /// Current canonical-basis dimension.
        dimension: usize,
    },
    /// Internal canonical-basis reconstruction requires equal list lengths.
    #[error(
        "Nakayama linear combination has {coefficients} coefficients for {basis} basis vectors"
    )]
    NakayamaLinearCombinationLength {
        /// Canonical-basis length.
        basis: usize,
        /// Coordinate-list length.
        coefficients: usize,
    },
    /// The Frobenius-form cache cannot contain GAP's `fail` sentinel.
    #[error("the FrobeniusForm attribute cache contains an invalid fail state")]
    InvalidFrobeniusFormCacheState,
    /// Successful Nakayama construction must install the form as a side effect.
    #[error("NakayamaAutomorphism did not install FrobeniusForm")]
    MissingInstalledFrobeniusForm,
    /// Concurrent evaluation cannot bind both false and a form matrix.
    #[error("the FrobeniusForm attribute was concurrently cached as false")]
    FrobeniusFormCacheConflict,
    /// The form closure accepts only elements of its captured quotient.
    #[error("a FrobeniusForm argument is not in the algebra")]
    FrobeniusFormElementOutsideAlgebra,
    /// The cached form matrix must match the canonical-basis dimension.
    #[error("the cached Frobenius form has shape {shape:?}, expected {dimension} x {dimension}")]
    InvalidFrobeniusFormMatrixShape {
        /// Cached matrix shape.
        shape: (usize, usize),
        /// Current canonical-basis dimension.
        dimension: usize,
    },
    /// The Frobenius-functional cache cannot contain GAP's `fail` sentinel.
    #[error("the FrobeniusLinearFunctional cache contains an invalid fail state")]
    InvalidFrobeniusFunctionalCacheState,
    /// Successful Nakayama construction must install the functional.
    #[error("NakayamaAutomorphism did not install FrobeniusLinearFunctional")]
    MissingInstalledFrobeniusFunctional,
    /// Concurrent evaluation cannot bind both false and functional coordinates.
    #[error("the FrobeniusLinearFunctional attribute was concurrently cached as false")]
    FrobeniusFunctionalCacheConflict,
    /// The functional closure accepts only elements of its captured quotient.
    #[error("the FrobeniusLinearFunctional argument is not in the algebra")]
    FrobeniusFunctionalElementOutsideAlgebra,
    /// Cached functional coordinates must match the canonical basis.
    #[error("the cached Frobenius functional has length {actual}, expected {dimension}")]
    InvalidFrobeniusFunctionalLength {
        /// Cached coordinate count.
        actual: usize,
        /// Current canonical-basis dimension.
        dimension: usize,
    },
    /// The order cache cannot contain GAP's `fail` sentinel.
    #[error("the OrderOfNakayamaAutomorphism cache contains an invalid fail state")]
    InvalidNakayamaOrderCacheState,
    /// Path-algebra order dispatch requires a selfinjectivity method.
    #[error("OrderOfNakayamaAutomorphism cannot dispatch IsSelfinjectiveAlgebra")]
    NakayamaOrderSelfinjectivityUnavailable,
    /// Concurrent order evaluation cannot bind both false and a value.
    #[error("OrderOfNakayamaAutomorphism was concurrently cached as false")]
    NakayamaOrderCacheConflict,
    /// QPA's quotient order method attempts to call a false automorphism.
    #[error("OrderOfNakayamaAutomorphism attempted to call boolean false")]
    CalledFalseNakayamaAutomorphism,
    /// GAP's matrix order operation requires a nonempty square matrix here.
    #[error("invalid Nakayama automorphism matrix shape {shape:?} for Order")]
    InvalidNakayamaOrderMatrixShape {
        /// Rejected matrix shape.
        shape: (usize, usize),
    },
    /// The proven cyclotomic-index bound must fit addressable iteration.
    #[error("the rational finite-order bound for dimension {dimension} overflows usize")]
    NakayamaOrderCyclotomicBoundOverflow {
        /// Matrix dimension.
        dimension: usize,
    },
    /// The supported prime characteristic must fit addressable arithmetic.
    #[error("prime characteristic {prime} does not fit usize")]
    NakayamaOrderPrimeTooLarge {
        /// Prime characteristic.
        prime: u64,
    },
    /// The unipotent exponent calculation must fit addressable dimensions.
    #[error("the GL({dimension},{prime}) unipotent exponent bound overflows usize")]
    NakayamaOrderPrimePowerOverflow {
        /// Matrix dimension.
        dimension: usize,
        /// Prime characteristic.
        prime: u64,
    },
    /// An annihilator basis requires QPA's finite canonical algebra basis.
    #[error("AnnihilatorOfModule requires a finite-dimensional algebra")]
    InfiniteAnnihilatorAlgebra,
    /// QPA's finite canonical path basis starts with all vertex idempotents.
    #[error("the finite canonical algebra basis does not start with the quiver vertices")]
    MissingCanonicalVertexBasis,
    /// Positive summand counts make both random projective lists nonempty.
    #[error("an internal RandomModule projective direct sum returned fail")]
    InternalRandomModuleDirectSumFailed,
    /// Both random projective sums act over the identical entered algebra.
    #[error("HomOverAlgebra rejected internal RandomModule projective sums")]
    InternalRandomModuleHomRejected,
    /// `TensorProductOfModules(M, N)` requires `N` over the opposite algebra.
    #[error("the entered modules are not over the appropriate algebras")]
    TensorProductModulesOverWrongAlgebras,
    /// Tensor-presentation dimensions must fit Rust's addressable size.
    #[error("a TensorProductOfModules dimension pairing overflows usize")]
    TensorProductDimensionOverflow,
    /// Every internal indecomposable projective has one distinguished generator.
    #[error("an internal TensorProductOfModules projective has no generator")]
    EmptyTensorProjectiveGenerator,
    /// Distinguished internal projective generators are vertex-uniform.
    #[error("an internal TensorProductOfModules projective generator is not uniform")]
    NonuniformTensorProjectiveGenerator,
    /// The source applies quotient-only nested element extraction to `P1`
    /// generators even when the acting algebra is a path algebra.
    #[error(
        "QPA's nontrivial path-algebra TensorProductOfModules branch cannot extract P1 vertices"
    )]
    TensorProductPathAlgebraP1Extraction,
    /// Both internal projective terms are nonempty in this source branch.
    #[error("an internal TensorProductOfModules projective presentation is empty")]
    EmptyTensorProjectivePresentation,
    /// A nonzero `P1 tensor N` supplies at least one physical relation row.
    #[error("the internal TensorProductOfModules relation matrix is empty")]
    EmptyTensorRelationMatrix,
    /// GAP's Steinitz extension must complete the tensor relation basis.
    #[error(
        "NaturalHomomorphismBySubspace did not produce a complete tensor-product quotient basis"
    )]
    IncompleteTensorProductQuotientBasis,
    /// QPA's `P1 tensor N = 0` closure references a never-assigned local.
    #[error("QPA TensorProductOfModules elementary closure has unassigned local size_matrix")]
    TensorProductUnboundSizeMatrix,
    /// QPA overwrites the closure's captured projection count with the last
    /// `P1` vertex and can consequently index `projectionsP0` out of range.
    #[error(
        "QPA TensorProductOfModules closure requests projection {index}, but only {count} exist"
    )]
    TensorProductClosureProjectionOutOfRange {
        /// One-based projection requested by QPA.
        index: usize,
        /// Number of projective-cover projections.
        count: usize,
    },
    /// The same captured-local bug can make QPA build too few acted summands
    /// for its subsequent `size_matrix[2]` loop.
    #[error(
        "QPA TensorProductOfModules closure requests acted summand {index}, but only {count} exist"
    )]
    TensorProductClosureActedIndexOutOfRange {
        /// One-based acted summand requested by QPA.
        index: usize,
        /// Number of acted summands constructed by the overwritten count.
        count: usize,
    },
    /// The projective cover is surjective, so every entered left element has
    /// a preimage when it belongs to the captured module.
    #[error("a tensor-product left element unexpectedly has no projective-cover preimage")]
    MissingTensorProjectiveCoverPreimage,
    /// `ExtOverAlgebra` requires one common acting algebra.
    #[error("the two modules entered are not modules over the same algebra")]
    ExtModulesOverDifferentAlgebras,
    /// A Hom call cannot reject the already checked Ext modules.
    #[error("HomOverAlgebra rejected modules in the ExtOverAlgebra construction")]
    InternalExtHomRejected,
    /// GAP's Steinitz extension must complete the image basis to the Hom basis.
    #[error("NaturalHomomorphismBySubspace did not produce a complete Ext quotient basis")]
    IncompleteExtQuotientBasis,
    /// The coefficient closure only accepts raw vectors over its original field.
    #[error("the Ext representative map uses a different coefficient field")]
    ExtCoefficientFieldMismatch,
    /// The coefficient closure accepts only maps in `Hom(Omega(M), N)`.
    #[error("the entered map is outside the Hom space used by ExtOverAlgebra")]
    ExtMapOutsideHomSpace,
    /// Internal star Hom calls use modules over the already checked algebra.
    #[error("HomOverAlgebra rejected an internal StarOfModule Hom space")]
    InternalStarHomRejected,
    /// Every indecomposable projective has its distinguished vertex generator.
    #[error("an indecomposable projective has no generator in StarOfModule")]
    EmptyStarProjectiveGenerator,
    /// Arrow multiplication preserves one-vertex support.
    #[error("left multiplication by an arrow did not define a projective map")]
    StarArrowMapNotUniform,
    /// A nonzero star vertex has a raw physical Hom basis.
    #[error("a nonzero StarOfModule vertex has no coordinate basis")]
    MissingStarCoordinateBasis,
    /// Composition with left multiplication remains in the corresponding Hom space.
    #[error("a StarOfModule composite lies outside its Hom basis")]
    StarCompositionOutsideHomBasis,
    /// Both Hom-space families used by the evaluation construction share one algebra.
    #[error("HomOverAlgebra rejected an internal double-star Hom space")]
    InternalDoubleStarHomRejected,
    /// A nonzero double-star vertex has a nonempty raw physical Hom basis.
    #[error("double-star vertex {vertex} has no Hom coordinate basis")]
    MissingDoubleStarHomCoordinateBasis {
        /// One-based vertex position.
        vertex: usize,
    },
    /// Evaluation at a module basis vector is an element of the computed double dual.
    #[error("the evaluation map at vertex {vertex} lies outside the double-star Hom basis")]
    DoubleStarEvaluationOutsideHomBasis {
        /// One-based vertex position.
        vertex: usize,
    },
    /// Finite projective bases require the opposite algebra's canonical basis.
    #[error("the opposite algebra has no canonical basis for double-star evaluation")]
    MissingDoubleStarCanonicalBasis,
    /// Every projective-basis element occurs in the opposite canonical basis.
    #[error(
        "opposite projective {projective}, vertex {vertex}, basis element {position} is absent from the canonical basis"
    )]
    MissingDoubleStarProjectiveBasisElement {
        /// One-based indecomposable-projective position.
        projective: usize,
        /// One-based representation vertex.
        vertex: usize,
        /// One-based position inside that vertex basis.
        position: usize,
    },
    /// Source-generated star actions satisfy the opposite relations.
    #[error("the generated StarOfModule does not satisfy the opposite-algebra relations")]
    GeneratedStarViolatesRelations,
    /// The explicit transpose construction requires QPA's finite canonical basis.
    #[error("TransposeOfModule requires a finite-dimensional algebra")]
    InfiniteTransposeAlgebra,
    /// Every minimal module generator is supported at a quiver vertex.
    #[error("a minimal generator has no vertex support in TransposeOfModule")]
    MissingTransposeGeneratorSupport,
    /// The source's internal projective sums are nonempty in this branch.
    #[error("an internal TransposeOfModule projective direct sum returned fail")]
    InternalTransposeDirectSumFailed,
    /// Projective covers created by the transpose-on-maps construction retain
    /// the direct-sum metadata required by QPA's lifting method.
    #[error("an internal TransposeOfModuleHomomorphism projective lifting returned fail")]
    InternalTransposeHomomorphismLiftFailed,
    /// The lifted projective-cover square commutes by construction.
    #[error("the internal TransposeOfModuleHomomorphism kernel square does not commute")]
    InternalTransposeHomomorphismKernelSquareFailed,
    /// Applying the contravariant star functor preserves the presentation
    /// square used to induce the transpose map.
    #[error("the internal TransposeOfModuleHomomorphism cokernel square does not commute")]
    InternalTransposeHomomorphismCokernelSquareFailed,
    /// Minimal syzygy generators must have exactly one right-support vertex.
    #[error("a minimal syzygy generator is not uniform in TransposeOfModule")]
    NonuniformTransposeSyzygyGenerator,
    /// Every opposite projective summand has its distinguished generator.
    #[error("an opposite projective has no generator in TransposeOfModule")]
    EmptyTransposeProjectiveGenerator,
    /// `DTr(M)` remains over the identical acting algebra of `M`.
    #[error("HomOverAlgebra rejected the internal tau-rigidity Hom space")]
    InternalTauRigidHomRejected,
    /// `SimpleModules` accepts quotients only after QPA's admissibility filter
    /// has been established.
    #[error("SimpleModules requires an admissible quotient of a path algebra")]
    NonadmissibleSimpleModuleAlgebra,
    /// A vertex-simple representation must satisfy every admissible relation.
    #[error("a generated vertex-simple module does not satisfy the algebra relations")]
    GeneratedSimpleViolatesRelations,
    /// A projective action generated from canonical basis multiplication must satisfy relations.
    #[error("a generated indecomposable projective does not satisfy the algebra relations")]
    GeneratedProjectiveViolatesRelations,
    /// Every projective right ideal contains its generating vertex idempotent.
    #[error("the canonical basis of projective {vertex} has no vertex generator")]
    MissingProjectiveVertexGenerator {
        /// One-based projective/vertex position.
        vertex: usize,
    },
    /// QPA's projectivity dimension sum exceeded Rust's addressable size.
    #[error("the projective dimension test overflows usize")]
    ProjectiveDimensionOverflow,
    /// `HomFromProjective` requires an element of the identical target module.
    #[error("the projective generator element is not in the target module")]
    ProjectiveElementInWrongModule,
    /// A projective basis cannot be extracted from infinite-dimensional input.
    #[error("BasisOfProjectives returned fail")]
    ProjectiveBasisUnavailable,
    /// Every nonzero module has at least one top/minimal generator.
    #[error("a nonzero projective cover has no summands")]
    EmptyNonzeroProjectiveCover,
    /// QPA's minimal generators are vertex-uniform.
    #[error("a minimal generator is not supported at exactly one vertex")]
    MinimalGeneratorNotVertexUniform,
    /// Projective lifting validates the source before its other inputs.
    #[error("the source of the second map is not projective")]
    LiftingSourceNotProjective,
    /// Projective lifting requires an epimorphism as its first map.
    #[error("the first map is not surjective")]
    LiftingFirstMapNotSurjective,
    /// Both maps in a projective lifting diagram must end at the same module.
    #[error("the entered maps do not have a common range")]
    LiftingMapsDoNotHaveCommonRange,
    /// Every nonzero direct projective summand has its installed generator.
    #[error("a nonzero projective summand has no minimal generator")]
    EmptyProjectiveSummandGenerators,
    /// Surjectivity guarantees a lift of every projective-generator image.
    #[error("a projective-generator image unexpectedly has no preimage")]
    MissingProjectiveLiftPreimage,
    /// A projective lift preimage must be vertex-uniform.
    #[error("a projective lifting preimage is not vertex-uniform")]
    ProjectiveLiftNotVertexUniform,
    /// Factoring-through computations require a common acting algebra.
    #[error("the factoring modules are not modules over the same algebra")]
    FactoringModulesOverDifferentAlgebras,
    /// QPA explicitly rejects the zero module for its stable endomorphism
    /// quotient operation.
    #[error("zero module given as input for EndModuleProjOverAlgebra")]
    ZeroStableEndomorphismModule,
    /// A Hom call cannot reject the checked projective-cover algebra.
    #[error("HomOverAlgebra rejected the projective cover's common algebra")]
    InternalProjectiveHomRejected,
    /// Flattened physical map data exceeded Rust's addressable size.
    #[error("the physical module-map size overflows usize")]
    PhysicalMapSizeOverflow,
    /// A raw physical map vector must fill every vertex block exactly.
    #[error("the physical module-map vector has length {actual}, expected {expected}")]
    PhysicalMapLength {
        /// Entered flat-vector length.
        actual: usize,
        /// Required total physical block length.
        expected: usize,
    },
    /// Projective positions are one-based quiver-vertex positions.
    #[error("projective vertex {vertex} is outside 1..={vertex_count}")]
    ProjectiveVertexOutOfRange {
        /// Entered vertex position.
        vertex: usize,
        /// Number of vertices/projectives.
        vertex_count: usize,
    },
    /// QPA flattens physical module coordinates against a basis with a zero
    /// placeholder at every zero-dimensional vertex space.
    #[error("projective element has {actual} physical coordinates, expected {expected}")]
    ProjectiveElementPhysicalCoordinateCount {
        /// Entered flattened physical coordinate count.
        actual: usize,
        /// Flattened projective-basis count.
        expected: usize,
    },
    /// The two source methods require an acyclic path algebra or an admissible
    /// quotient, respectively.
    #[error("the entered algebra is not finite dimensional in QPA's projective-element domain")]
    ElementInIndecProjectiveAlgebraNotFinite,
    /// The algebra-element-to-projective conversion requires the identical
    /// path algebra or quotient family.
    #[error("the projective algebra element belongs to another algebra")]
    ForeignProjectiveAlgebraElement,
    /// QPA rejects zero in `ElementIn_vA_AsElementInIndecProj`.
    #[error("the entered projective algebra element is zero")]
    ZeroProjectiveAlgebraElement,
    /// A nonzero element must have exactly one left-support vertex.
    #[error("the projective algebra element is not left uniform")]
    ProjectiveAlgebraElementNotLeftUniform,
    /// QPA's coefficient-bearing successive tips must occur in the selected
    /// projective-basis component.
    #[error(
        "a coefficient-bearing tip is absent from projective basis component ({left_vertex}, {right_vertex})"
    )]
    ProjectiveLeadingTermMissing {
        /// Generating/left-support vertex.
        left_vertex: usize,
        /// Right-support component vertex.
        right_vertex: usize,
    },
    /// A distinguished generator of an indecomposable projective must have
    /// exactly one support vertex.
    #[error("a projective summand generator is not vertex uniform")]
    ProjectiveGeneratorNotVertexUniform,
    /// GAP's `IsMatrix` method domain excludes an empty element matrix here.
    #[error("the projective element matrix must be nonempty")]
    EmptyProjectiveElementMatrix,
    /// GAP matrices are rectangular homogeneous lists.
    #[error("the projective element matrix is ragged")]
    RaggedProjectiveElementMatrix,
    /// Matrix rows/columns and the two entered support lists must agree.
    #[error(
        "projective matrix shape {rows}x{columns} disagrees with {range_vertices} range and {source_vertices} source vertices"
    )]
    ProjectiveMatrixVertexListDimensions {
        /// Matrix row count.
        rows: usize,
        /// Matrix column count.
        columns: usize,
        /// Entered range-vertex count.
        range_vertices: usize,
        /// Entered source-vertex count.
        source_vertices: usize,
    },
    /// Every matrix entry must belong to the identical entered algebra.
    #[error("a projective matrix entry belongs to another algebra")]
    ForeignProjectiveMatrixElement,
    /// Every nonzero matrix entry must be both left and right uniform.
    #[error("a projective matrix entry is not uniform")]
    ProjectiveMatrixElementNotUniform,
    /// All nonzero entries in a row must have one left support.
    #[error("projective matrix row {row} is not homogeneous")]
    ProjectiveMatrixRowNotHomogeneous {
        /// One-based row position.
        row: usize,
    },
    /// All nonzero entries in a column must have one right support.
    #[error("projective matrix column {column} is not homogeneous")]
    ProjectiveMatrixColumnNotHomogeneous {
        /// One-based column position.
        column: usize,
    },
    /// Detected supports must equal the entered lists, including their order.
    #[error("the entered vertex lists disagree with projective matrix support")]
    ProjectiveMatrixSupportDiscrepancy,
    /// Nonempty source/range vertex lists make these internal direct sums
    /// infallible in QPA's method domain.
    #[error("an internal nonempty projective direct sum returned fail")]
    InternalProjectiveDirectSumFailed,
    /// A nonzero uniform matrix entry must produce a projective homomorphism.
    #[error("HomFromProjective rejected an internally uniform matrix entry")]
    InternalProjectiveMatrixHomRejected,
    /// The configured QPA/GBNP finite-basis path failed.
    #[error(transparent)]
    Backend(#[from] PathalgGbnpError),
    /// Enveloping-algebra construction failed.
    #[error(transparent)]
    Tensor(#[from] TensorProductError),
    /// Product-quiver projection failed.
    #[error(transparent)]
    Product(#[from] QuiverProductError),
    /// Opposite-path reconstruction failed.
    #[error(transparent)]
    Path(#[from] PathError),
    /// Arithmetic in `OriginalPathAlgebra(A)` failed.
    #[error(transparent)]
    ParentPathAlgebra(#[from] PathAlgebraError),
    /// Path inclusion in the source algebra failed.
    #[error(transparent)]
    PathTensor(#[from] PathAlgebraTensorError),
    /// Source-algebra multiplication failed.
    #[error(transparent)]
    Element(#[from] QuiverAlgebraElementError),
    /// Quotient canonical-basis coordinates failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Exact matrix construction failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// The checked QPA right-module constructor rejected its generated data.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
    /// Generic QPA module dualization failed.
    #[error(transparent)]
    Functor(#[from] ModuleFunctorError),
}

/// Failure in QPA's generic module-dual functor.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum ModuleFunctorError {
    /// The generated transposed actions should satisfy the opposite relations.
    #[error("the generated dual module does not satisfy the opposite-algebra relations")]
    GeneratedDualViolatesRelations,
    /// Parent path-algebra opposite construction failed.
    #[error(transparent)]
    PathAlgebra(#[from] PathAlgebraError),
    /// Quotient opposite construction through QPA's GBNP backend failed.
    #[error(transparent)]
    Backend(#[from] PathalgGbnpError),
    /// The checked module constructor or reciprocal cache failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
}
