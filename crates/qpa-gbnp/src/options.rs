//! Typed counterparts of the GBNP options used by the Gröbner kernel.

/// Options affecting the source-translated GBNP computation.
///
/// The default matches `GBNP.ClearOptions`: no cancellative-monoid rule and
/// no prefix generators.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct GbnpOptions {
    prefix_generators: usize,
    cancellative_alphabet_size: Option<usize>,
    max_iterations: Option<usize>,
}

impl GbnpOptions {
    /// Number of negative prefix/module generators represented by a tree.
    #[must_use]
    pub const fn prefix_generators(&self) -> usize {
        self.prefix_generators
    }

    /// Set the number of negative prefix/module generators.
    pub const fn set_prefix_generators(&mut self, count: usize) {
        self.prefix_generators = count;
    }

    /// Enable GBNP's cancellative-monoid preprocessing.
    ///
    /// The source requires `Alphabetsize` whenever `CancellativeMonoid` is
    /// bound, although `ReduceCancellation` does not otherwise use the value.
    pub const fn enable_cancellative_monoid(&mut self, alphabet_size: usize) {
        self.cancellative_alphabet_size = Some(alphabet_size);
    }

    /// Disable GBNP's cancellative-monoid preprocessing.
    pub const fn disable_cancellative_monoid(&mut self) {
        self.cancellative_alphabet_size = None;
    }

    /// Return the configured alphabet size when cancellation is enabled.
    #[must_use]
    pub const fn cancellative_alphabet_size(&self) -> Option<usize> {
        self.cancellative_alphabet_size
    }

    /// Set the optional `SGrobner` iteration limit.
    pub const fn set_max_iterations(&mut self, limit: Option<usize>) {
        self.max_iterations = limit;
    }

    /// Optional `SGrobner` iteration limit.
    #[must_use]
    pub const fn max_iterations(&self) -> Option<usize> {
        self.max_iterations
    }
}
