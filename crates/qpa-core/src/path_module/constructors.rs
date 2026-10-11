/// Translate QPA's checked two-argument
/// `RightModuleOverPathAlgebra(A, arrow_actions)` method.
///
/// The dimension vector is inferred from every arrow matrix. Isolated
/// vertices receive dimension zero, exactly as in the checked source method.
/// `Ok(None)` is QPA's `fail` when quotient relations are not satisfied.
///
/// # Errors
///
/// Rejects incomplete or malformed action lists, incompatible inferred
/// dimensions, matrices over another field, and invalid zero-space notation;
/// propagates relation-matrix arithmetic failures.
pub fn right_module_over_path_algebra_infer_dimension<F: Field>(
    algebra: QuiverAlgebra<F>,
    actions: &PathAlgebraModuleArrowActions<F>,
) -> Result<Option<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    right_module_over_path_algebra_infer_dimension_impl(algebra, actions, true)
}
/// Translate QPA's unchecked two-argument
/// `RightModuleOverPathAlgebraNC(A, arrow_actions)` method.
///
/// This skips the source method's coefficient-field membership test but keeps
/// its dimension and quotient-relation checks. QPA accidentally leaves
/// isolated-vertex dimensions at `-1`; Rust reports that source failure as a
/// typed error rather than constructing an invalid module family.
///
/// # Errors
///
/// Rejects malformed input, dimension conflicts, and the source method's
/// isolated-vertex bug; propagates relation-matrix arithmetic failures.
pub fn right_module_over_path_algebra_infer_dimension_nc<F: Field>(
    algebra: QuiverAlgebra<F>,
    actions: &PathAlgebraModuleArrowActions<F>,
) -> Result<Option<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    right_module_over_path_algebra_infer_dimension_impl(algebra, actions, false)
}

#[allow(clippy::too_many_lines)]
fn right_module_over_path_algebra_infer_dimension_impl<F: Field>(
    algebra: QuiverAlgebra<F>,
    actions: &PathAlgebraModuleArrowActions<F>,
    checked_field: bool,
) -> Result<Option<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    let quiver = algebra.quiver();
    let arrow_count = quiver.number_of_arrows();
    let entered_count = match actions {
        PathAlgebraModuleArrowActions::Named(entries) => entries.len(),
        PathAlgebraModuleArrowActions::Positional(entries) => entries.len(),
        PathAlgebraModuleArrowActions::Indexed(entries) => entries.len(),
    };
    if entered_count < arrow_count {
        return Err(PathAlgebraModuleError::IncompleteInferredArrowActions {
            actual: entered_count,
            expected: arrow_count,
        });
    }
    if entered_count == 0 {
        return Err(PathAlgebraModuleError::EmptyInferredArrowActions);
    }

    let mut assigned = vec![None; arrow_count];
    match actions {
        PathAlgebraModuleArrowActions::Named(entries) => {
            for entry in entries {
                if checked_field {
                    check_inferred_action_field(entry.label(), entry.action(), algebra.field())?;
                }
                let Some(Generator::Arrow(arrow)) = quiver.generator_by_name(entry.label()) else {
                    return Err(PathAlgebraModuleError::UnknownArrowLabels {
                        labels: vec![entry.label().into()],
                    });
                };
                assigned[arrow.index().zero_based()] = Some(entry.action().clone());
            }
        }
        PathAlgebraModuleArrowActions::Positional(entries) => {
            if checked_field {
                for (position, matrix) in entries.iter().enumerate() {
                    if matrix.field() != algebra.field() {
                        let label = quiver.arrows().nth(position).map_or_else(
                            || format!("position {}", position + 1),
                            |arrow| arrow.name().to_owned(),
                        );
                        return Err(PathAlgebraModuleError::DifferentCoefficientField {
                            label: label.into(),
                        });
                    }
                }
            }
            for (slot, matrix) in assigned.iter_mut().zip(entries) {
                *slot = Some(PathAlgebraModuleArrowAction::Matrix(matrix.clone()));
            }
        }
        PathAlgebraModuleArrowActions::Indexed(entries) => {
            for entry in entries {
                let Some(arrow) = quiver.arrow(entry.arrow()) else {
                    return Err(PathAlgebraModuleError::UnknownArrowIndex {
                        index: entry.arrow().get(),
                    });
                };
                if checked_field {
                    check_inferred_action_field(arrow.name(), entry.action(), algebra.field())?;
                }
                assigned[entry.arrow().zero_based()] = Some(entry.action().clone());
            }
        }
    }

    let missing = quiver
        .arrows()
        .zip(&assigned)
        .filter_map(|(arrow, action)| action.is_none().then(|| Box::<str>::from(arrow.name())))
        .collect::<Vec<_>>();
    if !missing.is_empty() {
        return Err(PathAlgebraModuleError::UnassignedInferredArrowActions { labels: missing });
    }

    let field = algebra.field().clone();
    let mut dimensions = vec![None; quiver.number_of_vertices()];
    let mut matrices = Vec::with_capacity(arrow_count);
    for (arrow, action) in quiver.arrows().zip(assigned) {
        let action = action.expect("missing inferred actions were rejected above");
        let (source_dimension, target_dimension, matrix) =
            normalize_inferred_arrow_action(&action, arrow.name(), field.clone())?;
        assign_inferred_vertex_dimension(
            &mut dimensions,
            arrow.source().index().zero_based(),
            source_dimension,
            arrow.name(),
        )?;
        assign_inferred_vertex_dimension(
            &mut dimensions,
            arrow.target().index().zero_based(),
            target_dimension,
            arrow.name(),
        )?;
        matrices.push(matrix);
    }

    let isolated_vertices = dimensions
        .iter()
        .enumerate()
        .filter_map(|(vertex, dimension)| dimension.is_none().then_some(vertex + 1))
        .collect::<Vec<_>>();
    if !checked_field && !isolated_vertices.is_empty() {
        return Err(PathAlgebraModuleError::UncheckedInferredIsolatedVertices {
            vertices: isolated_vertices,
        });
    }
    let dimension_vector = dimensions
        .into_iter()
        .map(Option::unwrap_or_default)
        .collect::<Vec<_>>();
    let dimension = dimension_vector.iter().try_fold(0_usize, |sum, &entry| {
        sum.checked_add(entry)
            .ok_or(PathAlgebraModuleError::DimensionOverflow)
    })?;

    if dimension != 0
        && algebra.as_quotient().is_some()
        && !relations_are_satisfied(&algebra, &dimension_vector, &matrices, &field)?
    {
        return Ok(None);
    }
    Ok(Some(PathAlgebraModule(Arc::new(
        PathAlgebraModuleData::new(algebra, dimension_vector, dimension, matrices, None),
    ))))
}

fn check_inferred_action_field<F: Field>(
    label: &str,
    action: &PathAlgebraModuleArrowAction<F>,
    field: &F,
) -> Result<(), PathAlgebraModuleError> {
    if let PathAlgebraModuleArrowAction::Matrix(matrix) = action
        && matrix.field() != field
    {
        return Err(PathAlgebraModuleError::DifferentCoefficientField {
            label: label.into(),
        });
    }
    Ok(())
}

fn normalize_inferred_arrow_action<F: Field>(
    action: &PathAlgebraModuleArrowAction<F>,
    label: &str,
    field: F,
) -> Result<(usize, usize, DenseMatrix<F>), PathAlgebraModuleError> {
    match action {
        PathAlgebraModuleArrowAction::Matrix(matrix) => {
            let dimensions = (matrix.nrows(), matrix.ncols());
            if dimensions.0 == 0 || dimensions.1 == 0 {
                return Err(
                    PathAlgebraModuleError::ZeroDimensionMatrixInInferredConstructor {
                        label: label.into(),
                        dimensions,
                    },
                );
            }
            Ok((dimensions.0, dimensions.1, matrix.clone()))
        }
        PathAlgebraModuleArrowAction::ZeroSpaceDimensions { source, target } => {
            if *source != 0 && *target != 0 {
                return Err(PathAlgebraModuleError::InvalidZeroSpaceDimensions {
                    label: label.into(),
                    dimensions: (*source, *target),
                });
            }
            Ok((
                *source,
                *target,
                DenseMatrix::zero(
                    physical_dimension(*source),
                    physical_dimension(*target),
                    field,
                ),
            ))
        }
    }
}

fn assign_inferred_vertex_dimension(
    dimensions: &mut [Option<usize>],
    vertex: usize,
    entered: usize,
    arrow: &str,
) -> Result<(), PathAlgebraModuleError> {
    match dimensions[vertex] {
        None => dimensions[vertex] = Some(entered),
        Some(expected) if expected != entered => {
            return Err(
                PathAlgebraModuleError::IncompatibleInferredVertexDimension {
                    vertex: vertex + 1,
                    arrow: arrow.into(),
                    entered,
                    expected,
                },
            );
        }
        Some(_) => {}
    }
    Ok(())
}

/// Translate QPA's checked three-argument
/// `RightModuleOverPathAlgebra(A, dimension_vector, generators)`.
///
/// `Ok(None)` is QPA's `fail` result when the entered actions do not satisfy
/// the quotient relations. Malformed input is returned as an error, matching
/// the source method's GAP errors.
///
/// # Errors
///
/// Rejects a malformed dimension vector, matrices over another field, unknown
/// labels, nonzero maps incident with a zero space, wrong shapes, and matrix
/// arithmetic failures during relation evaluation.
pub fn right_module_over_path_algebra<F: Field>(
    algebra: QuiverAlgebra<F>,
    dimension_vector: Vec<usize>,
    generators: &[NamedArrowMatrix<F>],
) -> Result<Option<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    right_module_over_path_algebra_impl(algebra, dimension_vector, generators, None, false)
}

/// Construct a checked module while retaining the tensor object underlying
/// its acting algebra. This is the cross-crate hook used by QPA's enveloping
/// algebra attributes.
#[doc(hidden)]
pub fn right_module_over_path_algebra_with_tensor<F: Field>(
    algebra: QuiverAlgebra<F>,
    dimension_vector: Vec<usize>,
    generators: &[NamedArrowMatrix<F>],
    tensor: TensorProductAlgebra<F>,
) -> Result<Option<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    if tensor.algebra() != &algebra {
        return Err(PathAlgebraModuleError::RetainedTensorAlgebraMismatch);
    }
    right_module_over_path_algebra_impl(algebra, dimension_vector, generators, Some(tensor), false)
}

fn right_module_over_path_algebra_impl<F: Field>(
    algebra: QuiverAlgebra<F>,
    dimension_vector: Vec<usize>,
    generators: &[NamedArrowMatrix<F>],
    retained_tensor: Option<TensorProductAlgebra<F>>,
    relations_certified: bool,
) -> Result<Option<PathAlgebraModule<F>>, PathAlgebraModuleError> {
    let quiver = algebra.quiver();
    if dimension_vector.len() != quiver.number_of_vertices() {
        return Err(PathAlgebraModuleError::DimensionVectorLength {
            actual: dimension_vector.len(),
            expected: quiver.number_of_vertices(),
        });
    }
    let dimension = dimension_vector.iter().try_fold(0_usize, |sum, &entry| {
        sum.checked_add(entry)
            .ok_or(PathAlgebraModuleError::DimensionOverflow)
    })?;

    for generator in generators {
        if generator.matrix().field() != algebra.field() {
            return Err(PathAlgebraModuleError::DifferentCoefficientField {
                label: generator.label().into(),
            });
        }
    }
    let missing_labels = generators
        .iter()
        .filter(|generator| {
            !quiver
                .arrows()
                .any(|arrow| arrow.name() == generator.label())
        })
        .map(|generator| Box::<str>::from(generator.label()))
        .collect::<Vec<_>>();
    if !missing_labels.is_empty() {
        return Err(PathAlgebraModuleError::UnknownArrowLabels {
            labels: missing_labels,
        });
    }

    let field = algebra.field().clone();
    let mut arrow_matrices = quiver
        .arrows()
        .map(|arrow| {
            let source = dimension_vector[arrow.source().index().zero_based()];
            let target = dimension_vector[arrow.target().index().zero_based()];
            DenseMatrix::zero(
                physical_dimension(source),
                physical_dimension(target),
                field.clone(),
            )
        })
        .collect::<Vec<_>>();

    for generator in generators {
        let Generator::Arrow(arrow) = quiver
            .generator_by_name(generator.label())
            .expect("membership in the arrow-label list was checked above")
        else {
            unreachable!("an arrow name takes precedence in QPA's quiver record")
        };
        let source = dimension_vector[arrow.source().index().zero_based()];
        let target = dimension_vector[arrow.target().index().zero_based()];
        if source == 0 || target == 0 {
            if !generator.matrix().is_zero() {
                return Err(PathAlgebraModuleError::NonzeroMapAtZeroSpace {
                    label: generator.label().into(),
                });
            }
            return Err(PathAlgebraModuleError::ZeroSpaceMapMustBeOmitted {
                label: generator.label().into(),
            });
        }
        let actual = (generator.matrix().nrows(), generator.matrix().ncols());
        let expected = (source, target);
        if actual != expected {
            return Err(PathAlgebraModuleError::WrongArrowMatrixShape {
                label: generator.label().into(),
                actual,
                expected,
            });
        }
        arrow_matrices[arrow.index().zero_based()] = generator.matrix().clone();
    }

    if !relations_certified && dimension != 0
        && algebra.as_quotient().is_some()
        && !relations_are_satisfied(&algebra, &dimension_vector, &arrow_matrices, &field)?
    {
        return Ok(None);
    }

    Ok(Some(PathAlgebraModule(Arc::new(
        PathAlgebraModuleData::new(
            algebra,
            dimension_vector,
            dimension,
            arrow_matrices,
            retained_tensor,
        ),
    ))))
}

fn physical_dimension(dimension: usize) -> usize {
    dimension.max(1)
}

#[cfg(test)]
fn relations_are_satisfied_uncached<F: Field>(
    algebra: &QuiverAlgebra<F>,
    dimension_vector: &[usize],
    arrow_matrices: &[DenseMatrix<F>],
    field: &F,
) -> Result<bool, MatrixError> {
    for relation in algebra.relations() {
        let mut value: Option<DenseMatrix<F>> = None;
        for term in relation.terms() {
            let term_value =
                path_action_matrix(term.path(), dimension_vector, arrow_matrices, field.clone())?
                    .scale_left(term.coefficient());
            value = Some(match value {
                Some(sum) => sum.add(&term_value)?,
                None => term_value,
            });
        }
        if value.is_some_and(|matrix| !matrix.is_zero()) {
            return Ok(false);
        }
    }
    Ok(true)
}

#[cfg(test)]
fn path_action_matrix<F: Field>(
    path: &Path,
    dimension_vector: &[usize],
    arrow_matrices: &[DenseMatrix<F>],
    field: F,
) -> Result<DenseMatrix<F>, MatrixError> {
    let Some(walk) = path.walk() else {
        unreachable!("normalized path-algebra relations contain no zero-path term")
    };
    let mut walk = walk.iter();
    let Some(first) = walk.next() else {
        let PathEndpoint::Vertex(vertex) = path.source() else {
            unreachable!("a vertex path has a vertex source")
        };
        return Ok(DenseMatrix::identity(
            physical_dimension(dimension_vector[vertex.index().zero_based()]),
            field,
        ));
    };
    let mut result = arrow_matrices[first.zero_based()].clone();
    for arrow in walk {
        result = result.multiply(&arrow_matrices[arrow.zero_based()])?;
    }
    Ok(result)
}
