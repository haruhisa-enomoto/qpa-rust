/// Failure in QPA's restriction of modules and maps along an algebra map.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum RestrictionViaAlgebraHomomorphismError {
    /// The entered module must act over the algebra-map range.
    #[error("the entered module is not over the range of the algebra homomorphism")]
    ModuleNotOverRangeAlgebra,
    /// The entered module map must be a map over the algebra-map range.
    #[error("the entered homomorphism is not over the range of the algebra homomorphism")]
    HomomorphismNotOverRangeAlgebra,
    /// A restricted arrow action must land in its target idempotent space.
    #[error("a restricted arrow image lies outside its target vertex space")]
    ArrowImageOutsideRestrictedVertex,
    /// A restricted module-map image must land in its target idempotent space.
    #[error("a restricted map image lies outside its target vertex space")]
    MapImageOutsideRestrictedVertex,
    /// An algebra homomorphism must make the restricted action satisfy the source relations.
    #[error("the restricted module does not satisfy the source-algebra relations")]
    RestrictedModuleRejected,
    /// Algebra generator construction failed.
    #[error(transparent)]
    AlgebraElement(#[from] QuiverAlgebraElementError),
    /// Applying the algebra homomorphism failed.
    #[error(transparent)]
    AlgebraHomomorphism(#[from] QuiverAlgebraHomomorphismError),
    /// Exact row-space and coordinate computations failed.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// Original-module action or checked restricted module/map construction failed.
    #[error(transparent)]
    Module(#[from] PathAlgebraModuleError),
}

/// Failure in QPA's checked right-module constructor.
#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum PathAlgebraModuleError {
    /// QPA declares `LiftIdempotents` but installs no method for it.
    #[error("QPA installs no method for LiftIdempotents")]
    UninstalledLiftIdempotents,
    /// QPA's `Display(M)` passes the integer zero to `PrintArray` for a zero
    /// module over an arrowful quiver and therefore aborts.
    #[error("module Display encountered QPA's scalar-zero arrow matrix")]
    DisplayScalarZeroArrowMatrix,
    /// The dimension vector must have one entry per quiver vertex.
    #[error("dimension vector has length {actual}, expected {expected}")]
    DimensionVectorLength {
        /// Entered length.
        actual: usize,
        /// Required vertex count.
        expected: usize,
    },
    /// QPA requires at least one entered action for every quiver arrow.
    #[error("entered {actual} arrow actions, but the quiver has {expected} arrows")]
    IncompleteInferredArrowActions {
        /// Number of entered list elements.
        actual: usize,
        /// Number of quiver arrows.
        expected: usize,
    },
    /// The two-argument QPA method indexes its first list entry even for an
    /// arrowless quiver.
    #[error("QPA's dimension-inferring constructor cannot dispatch an empty action list")]
    EmptyInferredArrowActions,
    /// An arrow-object input used an index outside the acting quiver.
    #[error("arrow index {index} does not exist in the acting quiver")]
    UnknownArrowIndex {
        /// One-based entered arrow index.
        index: usize,
    },
    /// A labeled or arrow-object input repeated some arrows and omitted these.
    #[error("no inferred action was assigned to arrows {labels:?}")]
    UnassignedInferredArrowActions {
        /// Missing labels in quiver-arrow order.
        labels: Vec<Box<str>>,
    },
    /// A matrix with a zero row or column count cannot be a GAP matrix; QPA
    /// requires its integer-pair notation instead.
    #[error(
        "matrix for arrow {label:?} has zero dimension {dimensions:?}; use QPA's zero-space pair"
    )]
    ZeroDimensionMatrixInInferredConstructor {
        /// Rejected arrow label.
        label: Box<str>,
        /// Entered row and column counts.
        dimensions: (usize, usize),
    },
    /// QPA's integer-pair notation is reserved for an arrow incident with a
    /// zero-dimensional vertex.
    #[error("zero-space dimensions {dimensions:?} for arrow {label:?} contain no zero")]
    InvalidZeroSpaceDimensions {
        /// Rejected arrow label.
        label: Box<str>,
        /// Entered logical source and target dimensions.
        dimensions: (usize, usize),
    },
    /// Two incident arrow inputs inferred different dimensions at one vertex.
    #[error("arrow {arrow:?} infers dimension {entered} at vertex {vertex}, previously {expected}")]
    IncompatibleInferredVertexDimension {
        /// One-based vertex position.
        vertex: usize,
        /// Arrow whose input exposed the conflict.
        arrow: Box<str>,
        /// Newly entered dimension.
        entered: usize,
        /// Previously inferred dimension.
        expected: usize,
    },
    /// The NC source method forgets to replace isolated vertices' `-1`
    /// sentinels before constructing its element family.
    #[error(
        "QPA's unchecked dimension-inferring constructor leaves isolated vertices {vertices:?} at dimension -1"
    )]
    UncheckedInferredIsolatedVertices {
        /// One-based isolated vertex positions.
        vertices: Vec<usize>,
    },
    /// Total module dimension exceeded Rust's addressable size.
    #[error("the total module dimension overflows usize")]
    DimensionOverflow,
    /// A based generic module needs one right-action matrix per quiver vertex.
    #[error("entered {actual} based vertex actions, expected {expected}")]
    BasedModuleVertexActionCount {
        /// Entered matrix count.
        actual: usize,
        /// Quiver vertex count.
        expected: usize,
    },
    /// A based generic module needs one right-action matrix per quiver arrow.
    #[error("entered {actual} based arrow actions, expected {expected}")]
    BasedModuleArrowActionCount {
        /// Entered matrix count.
        actual: usize,
        /// Quiver arrow count.
        expected: usize,
    },
    /// Every based generator action must use the acting algebra's coefficient
    /// field.
    #[error("based action for {generator:?} uses a different coefficient field")]
    BasedModuleActionField {
        /// Vertex or arrow generator description.
        generator: Box<str>,
    },
    /// Every based generator action is an endomorphism of one common entered
    /// vector-space basis.
    #[error("based action for {generator:?} has shape {actual:?}, expected {expected:?}")]
    BasedModuleActionShape {
        /// Vertex or arrow generator description.
        generator: Box<str>,
        /// Entered matrix shape.
        actual: (usize, usize),
        /// Required common square shape.
        expected: (usize, usize),
    },
    /// An entered based action that is not actually a module action can carry
    /// a vertex-basis row outside the target idempotent space.
    #[error("based right action of arrow {arrow:?} lies outside its target vertex space")]
    BasedModuleArrowImageOutsideTarget {
        /// Offending arrow name.
        arrow: Box<str>,
    },
    /// An entered based action that violates quotient relations is rejected by
    /// QPA's final checked matrix-module constructor.
    #[error("the converted based module does not satisfy the entered algebra relations")]
    GeneratedBasedModuleConversionRejectedRelations,
    /// QPA's Hom-space equation matrix exceeded Rust's addressable size.
    #[error("the HomOverAlgebra equation matrix dimensions overflow usize")]
    HomEquationSizeOverflow,
    /// QPA installs its generic three-argument Hom method only for a path
    /// ring, not for a quotient of one.
    #[error("QPA's generic Hom/End method does not dispatch on this algebra")]
    GenericHomPathRingDispatchFailure,
    /// The explicit path-ring argument must be the identical acting algebra
    /// of both matrix-module objects.
    #[error("<ring> must be right acting domain of <M> and <N>")]
    GenericHomActingAlgebraMismatch,
    /// QPA's `HOPF_MatrixModuleHom` indexes the first ordinary basis vector
    /// and therefore fails for a zero source or range module.
    #[error("QPA's generic Hom construction cannot use a zero module")]
    ZeroGenericHomModuleFailure,
    /// GAP loses the intended nullspace width when QPA transposes an empty
    /// equation matrix.
    #[error("QPA's generic Hom construction cannot take this empty matrix nullspace")]
    DegenerateGenericHomNullspaceFailure,
    /// QPA's generic Hom equation or full mapping matrix exceeded Rust's
    /// addressable size.
    #[error("the generic QPA Hom equation matrix dimensions overflow usize")]
    GenericHomEquationSizeOverflow,
    /// A generic basis-matrix mapping accepts elements only from the ambient
    /// family of its source object.
    #[error("the entered element has a foreign generic-Hom source family")]
    ForeignGenericHomomorphismSourceElement,
    /// A parent-module element need not belong to a generic map's submodule
    /// source.
    #[error("the entered element is outside the generic-Hom source object")]
    ElementOutsideGenericHomomorphismSource,
    /// Multiplication in a generic End algebra accepts only endomorphisms of
    /// its identical retained matrix-module object.
    #[error("the entered mapping is outside this generic endomorphism algebra")]
    ForeignGenericEndomorphism,
    /// A supported vertex must have a row block in QPA's Hom unknown vector.
    #[error("supported vertex {vertex} has no HomOverAlgebra unknown block")]
    MissingHomUnknownBlock {
        /// One-based vertex position.
        vertex: usize,
    },
    /// An endomorphism-space query cannot fail its identical algebra check.
    #[error("HomOverAlgebra unexpectedly rejected an endomorphism pair")]
    EndomorphismAlgebraMismatch,
    /// A reverse Hom query for one module map necessarily has a common
    /// acting algebra.
    #[error("HomOverAlgebra unexpectedly rejected a splitting pair")]
    InternalSplittingHomRejected,
    /// Common-summand inputs already passed their common-algebra test.
    #[error("HomOverAlgebra unexpectedly rejected a common-summand pair")]
    InternalCommonSummandHomRejected,
    /// Additive-closure inputs already passed their common-algebra test.
    #[error("HomOverAlgebra unexpectedly rejected an additive-closure pair")]
    InternalAdditiveClosureHomRejected,
    /// A nonzero module has the identity in its endomorphism basis.
    #[error("a minimality test unexpectedly received an empty endomorphism basis")]
    EmptyMinimalityEndomorphismBasis,
    /// Minimality endomorphism inputs necessarily share one acting algebra.
    #[error("HomOverAlgebra unexpectedly rejected a minimality endomorphism pair")]
    InternalMinimalityHomRejected,
    /// Kernel/cokernel-derived modules necessarily share the entered map's
    /// acting algebra during a minimal-version computation.
    #[error("HomOverAlgebra unexpectedly rejected a minimal-version pair")]
    InternalMinimalHomRejected,
    /// QPA's Fitting decomposition requires this canonical kernel inclusion
    /// to have the retraction used by `MoreLeftMinimalVersion`.
    #[error("the internal Fitting kernel inclusion has no right inverse")]
    InternalMinimalKernelNotSplit,
    /// `InverseOfIsomorphism` rejects a map before inspecting its matrices.
    #[error("entered map is not an isomorphism")]
    EnteredMapNotIsomorphism,
    /// Every nonzero component matrix of a module isomorphism is invertible.
    #[error("an isomorphism has a singular nonzero vertex matrix")]
    InternalIsomorphismMatrixSingular,
    /// The two modules in QPA's recursive isomorphism search have already
    /// passed their common-algebra check.
    #[error("HomOverAlgebra unexpectedly rejected an isomorphism-search pair")]
    InternalIsomorphismHomRejected,
    /// QPA's Fitting summands and complements must have canonical retractions.
    #[error("an internal isomorphism-search summand inclusion has no right inverse")]
    InternalIsomorphismSummandNotSplit,
    /// QPA's common-summand route must land in the chosen range summand.
    #[error("a common-summand image has no preimage in the range summand")]
    MissingIsomorphismSummandPreimage,
    /// GAP's splitting solver returns one coefficient per reverse-Hom basis
    /// vector.
    #[error("splitting combination has {actual} coefficients, expected {expected}")]
    SplittingCoefficientCount {
        /// Entered coefficient count.
        actual: usize,
        /// Reverse-Hom basis length.
        expected: usize,
    },
    /// Internal linear-combination construction follows a nonempty-basis
    /// source branch.
    #[error("an internal splitting combination has an empty Hom basis")]
    EmptySplittingHomBasis,
    /// GAP's `DirectSumOfAlgebras` rejects the empty list for the zero module.
    #[error("EndOverAlgebra cannot form a direct sum for a zero module")]
    EmptyEndomorphismDirectSum,
    /// `FromHomMMToEndM` is defined for an endomorphism of one module.
    #[error("FromHomMMToEndM requires structurally equal source and range modules")]
    NonEndomorphismConversion,
    /// `FromEndMToHomMM` requires the module's coefficient field.
    #[error("the endomorphism matrix is over a different coefficient field")]
    EndomorphismMatrixField,
    /// QPA's diagonal-block extraction indexes the full underlying module.
    #[error("endomorphism matrix has shape {actual:?}, but must be at least {minimum:?}")]
    EndomorphismMatrixTooSmall {
        /// Entered matrix shape.
        actual: (usize, usize),
        /// Minimum row and column counts.
        minimum: (usize, usize),
    },
    /// A natural quotient accepts only matrices represented in its source.
    #[error("the matrix is not an element of the represented endomorphism algebra")]
    MatrixOutsideEndomorphismAlgebra,
    /// Every generator of an endomorphism ideal belongs to the parent algebra.
    #[error("an ideal generator is outside the represented endomorphism algebra")]
    EndomorphismIdealGeneratorOutsideAlgebra,
    /// Products of parent and ideal elements must remain in the parent.
    #[error("an endomorphism-ideal product left the represented parent algebra")]
    EndomorphismIdealProductOutsideAlgebra,
    /// Multiplication in a represented endomorphism algebra must remain in its
    /// retained basis before the radical conversion can be applied.
    #[error("an endomorphism-algebra product left the represented basis")]
    EndomorphismAlgebraProductOutsideBasis,
    /// Quotient multiplication requires one identical quotient family.
    #[error("the quotient elements belong to different endomorphism quotient families")]
    ForeignEndomorphismQuotientElement,
    /// Quotient elements carry one coefficient per quotient basis element.
    #[error("endomorphism quotient element has {actual} coordinates, expected {expected}")]
    EndomorphismQuotientElementCoordinateCount {
        /// Entered coordinate count.
        actual: usize,
        /// Quotient dimension.
        expected: usize,
    },
    /// The arrow-kernel construction must define a module satisfying the relations.
    #[error("the generated socle module does not satisfy the algebra relations")]
    GeneratedSocleViolatesRelations,
    /// The all-zero representation must satisfy every supported algebra relation.
    #[error("the generated zero module does not satisfy the algebra relations")]
    GeneratedZeroModuleViolatesRelations,
    /// Every direct summand must act over the same algebra.
    #[error("all direct summands must be modules over the same algebra")]
    DirectSumDifferentAlgebras,
    /// The source block-diagonal action must satisfy the common relations.
    #[error("the generated direct sum does not satisfy the algebra relations")]
    GeneratedDirectSumViolatesRelations,
    /// Direct-sum metadata can be installed only on its fresh module object.
    #[error("direct-sum metadata was already installed")]
    DirectSumMetadataAlreadySet,
    /// Projection and inclusion attributes require a QPA direct sum.
    #[error("the entered module was not built by DirectSumOfQPAModules")]
    NotDirectSumModule,
    /// An explicitly installed minimal generator must belong to the module.
    #[error("a minimal generator does not belong to the entered module")]
    MinimalGeneratorInWrongModule,
    /// Installing a minimal-generator attribute must make it readable.
    #[error("the minimal-generating-set cache could not be initialized")]
    MinimalGeneratingSetCacheFailed,
    /// Every internal use of a direct sum supplies at least one module.
    #[error("an internal nonempty DirectSumOfQPAModules call returned fail")]
    InternalNonemptyDirectSumFailed,
    /// A commutative square must carry the first kernel into the second.
    #[error("the induced kernel image unexpectedly has no preimage")]
    MissingInducedKernelPreimage,
    /// The canonical image projection must lift every image-basis vector.
    #[error("an image basis vector unexpectedly has no preimage")]
    MissingInducedImagePreimage,
    /// The canonical cokernel projection must lift every cokernel-basis vector.
    #[error("a cokernel basis vector unexpectedly has no preimage")]
    MissingInducedCokernelPreimage,
    /// A successful inclusion-containment test must provide every lift.
    #[error("an included basis image unexpectedly has no preimage")]
    MissingInclusionLiftPreimage,
    /// A basis-image map needs at least one image per source basis vector.
    #[error("basis-image list has length {actual}, expected at least {expected}")]
    BasisImageCount {
        /// Entered image count.
        actual: usize,
        /// Required source dimension.
        expected: usize,
    },
    /// QPA requires a nonzero basis image at the source basis vector's vertex.
    #[error("basis image {position} has the wrong vertex support")]
    WrongBasisImageSupport {
        /// One-based position in the entered image list.
        position: usize,
    },
    /// The map-matrix utility cannot infer a map family from no projections.
    #[error("MultiplyListsOfMaps requires a nonempty projection list")]
    EmptyMapProjectionList,
    /// The map-matrix utility cannot return a map from no inclusions.
    #[error("MultiplyListsOfMaps requires a nonempty inclusion list")]
    EmptyMapInclusionList,
    /// QPA indexes one matrix row for every inclusion.
    #[error("map matrix has {actual} rows, expected at least {expected}")]
    MapMatrixRowCount {
        /// Entered row count.
        actual: usize,
        /// Required row count.
        expected: usize,
    },
    /// Each consumed matrix row is dotted with the whole projection list.
    #[error("map matrix row {row} has {actual} columns, expected {expected}")]
    MapMatrixColumnCount {
        /// One-based row position.
        row: usize,
        /// Entered column count.
        actual: usize,
        /// Required projection count.
        expected: usize,
    },
    /// The first binary submodule argument must be a monomorphism.
    #[error("the first submodule argument is not a monomorphism")]
    FirstSubmoduleMapNotInjective,
    /// The second binary submodule argument must be a monomorphism.
    #[error("the second submodule argument is not a monomorphism")]
    SecondSubmoduleMapNotInjective,
    /// Submodule inclusions must have the same structural range.
    #[error("submodule homomorphisms must have the same range")]
    SubmoduleRangesDiffer,
    /// QPA requires a nonempty list of submodule inclusions.
    #[error("the submodule homomorphism list must be nonempty")]
    EmptySubmoduleList,
    /// Every source in a submodule list must act over the same algebra.
    #[error("all submodule-list sources must be modules over the same algebra")]
    SubmoduleListDifferentAlgebras,
    /// Every map in a submodule list must be a monomorphism.
    #[error("not all submodule-list arguments are monomorphisms")]
    SomeSubmoduleMapNotInjective,
    /// Equal ranges make QPA's internal pullback branch infallible.
    #[error("an internal pullback of maps with a common range returned fail")]
    InternalCommonRangePullbackFailed,
    /// `RejectOfModule` requires a common acting algebra.
    #[error("the reject modules are not modules over the same algebra")]
    RejectModulesOverDifferentAlgebras,
    /// A Hom query cannot reject modules after their algebras were checked.
    #[error("an internal HomOverAlgebra call rejected a common algebra")]
    InternalCommonAlgebraHomFailed,
    /// A module element needs one coordinate row per quiver vertex.
    #[error("module element has {actual} vertex rows, expected {expected}")]
    ElementVertexCount {
        /// Entered row count.
        actual: usize,
        /// Required vertex count.
        expected: usize,
    },
    /// A logical coordinate row must match its vertex-space dimension.
    #[error("module element row at vertex {vertex} has length {actual}, expected {expected}")]
    ElementCoordinateCount {
        /// One-based vertex position.
        vertex: usize,
        /// Entered coordinate count.
        actual: usize,
        /// Required logical dimension.
        expected: usize,
    },
    /// A flat logical coordinate list must have total module dimension.
    #[error("module element has {actual} flat coordinates, expected {expected}")]
    ElementFlatCoordinateCount {
        /// Entered coordinate count.
        actual: usize,
        /// Required total dimension.
        expected: usize,
    },
    /// A raw element-space generator must belong to its ambient module family.
    #[error("an element-space generator belongs to a different module family")]
    ForeignElementSpaceGenerator,
    /// Basis vectors and coefficient queries require the ambient element
    /// family.
    #[error("an element basis received an element from a different module family")]
    ForeignElementBasisVector,
    /// `NewBasis` stores only structural basis data in QPA; `BasisNC` or
    /// `Basis` attaches the flattened delegate used for coefficients.
    #[error("the path-module element basis has no delegate coefficient basis")]
    MissingElementBasisDelegate,
    /// QPA's zero-module `CanonicalBasis` reaches an inapplicable `NewBasis`
    /// call before its explicit empty-basis branch.
    #[error("CanonicalBasis of QPA's zero path module fails while constructing NewBasis")]
    ZeroModuleCanonicalBasisFailure,
    /// GAP module-element arithmetic requires an identical element family.
    #[error("module elements belong to different element families")]
    DifferentElementFamilies,
    /// Right action requires an element of the identical acting algebra.
    #[error("the acting element belongs to a different algebra family")]
    ForeignActingElement,
    /// A normalized nonzero algebra element cannot contain QPA's zero path.
    #[error("the acting algebra element contains the zero path")]
    ZeroPathAction,
    /// Module homomorphisms require one common structurally equal algebra.
    #[error("the two modules are not over the same algebra")]
    DifferentHomomorphismAlgebras,
    /// A module homomorphism stores one matrix per vertex.
    #[error("received {actual} module-map matrices, expected {expected}")]
    HomomorphismMatrixCount {
        /// Entered matrix count.
        actual: usize,
        /// Required vertex count.
        expected: usize,
    },
    /// Every vertex map must use the acting algebra's coefficient field.
    #[error("module-map matrix at vertex {vertex} is over a different field")]
    HomomorphismMatrixField {
        /// One-based vertex position.
        vertex: usize,
    },
    /// A vertex map has the wrong physical shape or a nonzero zero-space map.
    #[error(
        "module-map matrix at vertex {vertex} has shape {actual:?}, expected zero matrix of shape {expected:?}"
    )]
    HomomorphismMatrixShape {
        /// One-based vertex position.
        vertex: usize,
        /// Entered shape.
        actual: (usize, usize),
        /// Required physical shape.
        expected: (usize, usize),
    },
    /// The entered vertex maps fail an arrow commutativity square.
    #[error("entered map does not commute with arrow {arrow:?}")]
    NoncommutingHomomorphismSquare {
        /// First failing arrow label.
        arrow: Box<str>,
    },
    /// Addition requires structurally equal source and range modules.
    #[error("the two maps do not belong to the same homomorphism set")]
    DifferentHomomorphismSets,
    /// Composition requires the first range to equal the second source.
    #[error("the first map's range is not the second map's source")]
    NoncomposableHomomorphisms,
    /// Mapping application requires an element of the identical source.
    #[error("the module element does not belong to the mapping source")]
    ForeignMapSourceElement,
    /// Preimage computation requires an element of the identical range.
    #[error("the module element does not belong to the mapping range")]
    ForeignMapRangeElement,
    /// Every generator of a subrepresentation must belong to the input module.
    #[error("a subrepresentation generator is not an element of the input module")]
    ForeignSubrepresentationGenerator,
    /// QPA acts by vertices of the stored original path algebra, which do not
    /// dispatch on outer module elements over a proper quotient.
    #[error("QPA's SubAlgebraModule method fails for modules over proper path-algebra quotients")]
    QuotientSubalgebraModuleFailure,
    /// QPA's identical-family method accepts only ambient module elements.
    #[error("a SubAlgebraModule generator belongs to a different module family")]
    ForeignSubalgebraModuleGenerator,
    /// `SubmoduleAsModule` indexes the first canonical-basis vector of its
    /// zero submodule.
    #[error("QPA cannot convert the zero SubAlgebraModule to a standalone module")]
    ZeroSubmoduleAsModuleFailure,
    /// QPA's generic module-basis delegate attempts to invert a nonsquare
    /// matrix when its retained uniform generators are dependent.
    #[error("QPA cannot construct a canonical basis from dependent SubAlgebraModule generators")]
    DependentSubalgebraModuleGenerators,
    /// `SubmoduleAsModule` feeds zero-space dimension pairs into the
    /// positional two-argument constructor, which QPA misdispatches.
    #[error("SubmoduleAsModule fails at zero-space arrow {arrow:?}")]
    ZeroSpaceSubmoduleAsModuleFailure {
        /// Arrow exposing a zero endpoint dimension.
        arrow: Box<str>,
    },
    /// The empty arrow list does not dispatch QPA's two-argument constructor.
    #[error("SubmoduleAsModule fails for an arrowless path algebra")]
    ArrowlessSubmoduleAsModuleFailure,
    /// A reconstructed submodule arrow image must lie in its target vertex
    /// subspace.
    #[error("an internal SubmoduleAsModule arrow image has no target-basis coordinates")]
    InternalSubmoduleCoordinates,
    /// The named natural mapping requires the entered submodule of `M`.
    #[error("the second argument is not a submodule of the first argument")]
    SubalgebraModuleNotInEnteredParent,
    /// GAP's generic subspace test cannot obtain images for the zero
    /// `SubAlgebraModule` in this method.
    #[error("NaturalHomomorphismBySubAlgebraModule fails for the zero submodule")]
    ZeroNaturalSubmoduleFailure,
    /// QPA feeds quotient zero-space pairs into its positional constructor.
    #[error("NaturalHomomorphismBySubAlgebraModule fails at zero-space arrow {arrow:?}")]
    ZeroSpaceNaturalSubmoduleFailure {
        /// Arrow exposing a zero quotient endpoint.
        arrow: Box<str>,
    },
    /// The empty arrow list does not dispatch QPA's quotient constructor.
    #[error("NaturalHomomorphismBySubAlgebraModule fails for an arrowless path algebra")]
    ArrowlessNaturalSubmoduleFailure,
    /// QPA indexes the first basis vector of its zero quotient.
    #[error("NaturalHomomorphismBySubAlgebraModule cannot construct a zero quotient mapping")]
    ZeroNaturalSubmoduleQuotientFailure,
    /// Basis-image mappings accept only their identical source family.
    #[error("the entered element is outside the basis mapping's source family")]
    ForeignBasisMappingSourceElement,
    /// QPA's canonical subspace basis must remain vertex-uniform.
    #[error("the canonical generated-submodule basis is not vertex-uniform")]
    NonuniformSubrepresentationBasis,
    /// Coordinates in the generated submodule's canonical basis must exist.
    #[error("failed to express a generated-submodule element in its canonical basis")]
    InternalSubrepresentationCoordinates,
    /// A kernel action must remain in the next component kernel.
    #[error("failed to express a kernel action in its component kernel basis")]
    InternalKernelCoordinates,
    /// A range action and every source image must remain in the image basis.
    #[error("failed to express a module-map image in its canonical image basis")]
    InternalImageCoordinates,
    /// Every basis vector of a constructed top must have QPA's chosen lift.
    #[error("a source-top basis vector has no preimage under its top projection")]
    MissingTopBasisPreimage,
    /// A source-generated module unexpectedly failed its defining relations.
    #[error("a module reconstructed by a QPA source method failed the acting-algebra relations")]
    GeneratedModuleRejectedRelations,
    /// QPA's radical method applies only in its finite/admissible source domain.
    #[error(
        "RadicalOfModuleInclusion requires a finite path algebra or finite admissible quotient"
    )]
    UnsupportedRadicalAlgebra,
    /// An entered matrix does not use the acting algebra's field instance.
    #[error("matrix for arrow {label:?} is over a different coefficient field")]
    DifferentCoefficientField {
        /// Arrow label attached to the matrix.
        label: Box<str>,
    },
    /// At least one entered label is absent from the quiver's arrow list.
    #[error("arrow labels {labels:?} do not exist in the acting quiver")]
    UnknownArrowLabels {
        /// Missing labels in entered order, including repetitions.
        labels: Vec<Box<str>>,
    },
    /// QPA forbids a nonzero action incident with a zero-dimensional space.
    #[error("a nonzero matrix was entered for arrow {label:?} at a zero-dimensional space")]
    NonzeroMapAtZeroSpace {
        /// Rejected arrow label.
        label: Box<str>,
    },
    /// GAP cannot retain a zero-row or zero-column entered matrix; such an
    /// arrow must be omitted so the constructor creates its physical stand-in.
    #[error("the zero matrix for arrow {label:?} at a zero-dimensional space must be omitted")]
    ZeroSpaceMapMustBeOmitted {
        /// Rejected arrow label.
        label: Box<str>,
    },
    /// An entered arrow map has dimensions different from its source and target spaces.
    #[error("matrix for arrow {label:?} has shape {actual:?}, expected {expected:?}")]
    WrongArrowMatrixShape {
        /// Rejected arrow label.
        label: Box<str>,
        /// Entered shape.
        actual: (usize, usize),
        /// Source-by-target dimension.
        expected: (usize, usize),
    },
    /// The retained tensor metadata must wrap the identical acting algebra.
    #[error("retained tensor product does not wrap the module's acting algebra")]
    RetainedTensorAlgebraMismatch,
    /// Concurrent construction lost QPA's reciprocal dual-module attribute.
    #[error("concurrent dual-module cache initialization failed")]
    DualCacheRace,
    /// Relation evaluation encountered incompatible explicit matrices.
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    /// Quiver-algebra element construction or arithmetic failed.
    #[error(transparent)]
    AlgebraElement(#[from] QuiverAlgebraElementError),
    /// Path-algebra element construction failed.
    #[error(transparent)]
    PathAlgebra(#[from] crate::PathAlgebraError),
    /// Quotient computation or normal-form reduction failed.
    #[error(transparent)]
    Quotient(#[from] QuotientError),
    /// Exact structure-constant conversion or radical computation failed.
    #[error(transparent)]
    StructureConstants(#[from] StructureConstantsError),
}
