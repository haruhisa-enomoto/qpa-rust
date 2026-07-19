"""Architecture spike for the source-faithful Rust QPA binding.

This package intentionally exposes only a vertical slice. Its job is to prove
the coefficient, ownership, identity, error, and object-lifetime design before
the complete Python binding campaign.
"""

from ._qpa import (
    AlgebraElement,
    Field,
    InvalidInputError,
    MatrixModule,
    Path,
    PathAlgebra,
    QpaError,
    Quiver,
    QuotientAlgebra,
    QuotientElement,
    RelationNotSatisfiedError,
    Scalar,
)

__all__ = [
    "AlgebraElement",
    "Field",
    "InvalidInputError",
    "MatrixModule",
    "Path",
    "PathAlgebra",
    "QpaError",
    "Quiver",
    "QuotientAlgebra",
    "QuotientElement",
    "RelationNotSatisfiedError",
    "Scalar",
]
