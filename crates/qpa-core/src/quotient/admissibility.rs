//! Sparse execution of QPA's arrow-ideal/product-space admissibility test.
//!
//! Keep the same descending sequence J, J^2, J^4, ... until dimensions
//! stabilize. Only its representation changes: stream GAP-style echelon rows
//! and cache multiplication by arrows instead of retaining all products and
//! recomputing their normal forms. This is not a finite-dimension shortcut:
//! finite quotients with nonnilpotent arrow ideal must still return false.

use std::collections::{BTreeMap, HashMap};

use numerica::domains::Field;

use super::{
    PathEndpoint, QuotientError, QuotientPathAlgebra, QuotientPathAlgebraElement,
    canonical_basis_path,
};
type Sparse<F> = Vec<(usize, <F as numerica::domains::Set>::Element)>;

/// GAP's incremental semi-echelon span, with absent coordinates stored as zero.
struct Span<F: Field> {
    rows: Vec<Sparse<F>>,
    pivots: Vec<Option<usize>>,
}

impl<F: Field> Span<F> {
    fn new(dimension: usize) -> Self {
        Self {
            rows: Vec::new(),
            pivots: vec![None; dimension],
        }
    }

    fn insert(&mut self, mut row: Sparse<F>, field: &F) {
        while let Some((pivot, coefficient)) = row.first().cloned() {
            if let Some(index) = self.pivots[pivot] {
                row = add_scaled(&row, &self.rows[index], &field.neg(&coefficient), field);
            } else {
                if !field.is_one(&coefficient) {
                    for (_, value) in &mut row {
                        *value = field.div(&*value, &coefficient);
                    }
                }
                self.pivots[pivot] = Some(self.rows.len());
                self.rows.push(row);
                break;
            }
        }
    }
}

fn add_scaled<F: Field>(
    left: &Sparse<F>,
    right: &Sparse<F>,
    scalar: &F::Element,
    field: &F,
) -> Sparse<F> {
    if field.is_zero(scalar) || right.is_empty() {
        return left.clone();
    }
    let mut result = Vec::with_capacity(left.len() + right.len());
    let (mut i, mut j) = (0, 0);
    while i < left.len() || j < right.len() {
        let entry = if j == right.len() || (i < left.len() && left[i].0 < right[j].0) {
            let value = left[i].clone();
            i += 1;
            value
        } else {
            let (column, value) = &right[j];
            let mut value = if field.is_one(scalar) {
                value.clone()
            } else {
                field.mul(value, scalar)
            };
            j += 1;
            if i < left.len() && left[i].0 == *column {
                value = field.add(&left[i].1, &value);
                i += 1;
            }
            (*column, value)
        };
        if !field.is_zero(&entry.1) {
            result.push(entry);
        }
    }
    result
}

struct Products<'a, F: Field> {
    quotient: &'a QuotientPathAlgebra<F>,
    basis: &'a [QuotientPathAlgebraElement<F>],
    indices: HashMap<Vec<usize>, usize>,
    corners: Vec<(usize, usize)>,
    // A nontip's prefix is a nontip. Vertices have no predecessor.
    prefixes: Vec<Option<(usize, usize)>>,
    arrows: Vec<QuotientPathAlgebraElement<F>>,
    arrow_corners: Vec<(usize, usize)>,
    actions: Vec<Option<Sparse<F>>>,
}

fn endpoint(endpoint: PathEndpoint) -> Result<usize, QuotientError> {
    match endpoint {
        PathEndpoint::Vertex(vertex) => Ok(vertex.index().zero_based()),
        PathEndpoint::Zero(_) => Err(QuotientError::InvalidCanonicalBasis),
    }
}

impl<'a, F: Field> Products<'a, F> {
    fn new(
        quotient: &'a QuotientPathAlgebra<F>,
        basis: &'a [QuotientPathAlgebraElement<F>],
    ) -> Result<Self, QuotientError> {
        let paths = basis
            .iter()
            .map(|element| canonical_basis_path(element.representative()))
            .collect::<Result<Vec<_>, _>>()?;
        let indices = paths
            .iter()
            .enumerate()
            .map(|(i, p)| (p.external_representation(), i))
            .collect();
        let corners = paths
            .iter()
            .map(|p| Ok((endpoint(p.source())?, endpoint(p.target())?)))
            .collect::<Result<Vec<_>, QuotientError>>()?;
        let mut vertices = HashMap::new();
        let mut words = HashMap::new();
        for (i, path) in paths.iter().enumerate() {
            if path.is_vertex() {
                vertices.insert(corners[i].0, i);
            } else {
                words.insert(
                    path.walk()
                        .ok_or(QuotientError::InvalidCanonicalBasis)?
                        .to_vec(),
                    i,
                );
            }
        }
        let prefixes = paths
            .iter()
            .enumerate()
            .map(|(i, path)| {
                let word = path.walk().ok_or(QuotientError::InvalidCanonicalBasis)?;
                let Some(last) = word.last() else {
                    return Ok(None);
                };
                let prefix = if word.len() == 1 {
                    vertices.get(&corners[i].0)
                } else {
                    words.get(&word[..word.len() - 1])
                }
                .copied()
                .ok_or(QuotientError::InvalidCanonicalBasis)?;
                Ok(Some((prefix, last.zero_based())))
            })
            .collect::<Result<Vec<_>, QuotientError>>()?;
        let arrow_paths = quotient
            .original_path_algebra()
            .quiver()
            .arrows()
            .map(|a| a.path())
            .collect::<Vec<_>>();
        let arrow_corners = arrow_paths
            .iter()
            .map(|p| Ok((endpoint(p.source())?, endpoint(p.target())?)))
            .collect::<Result<Vec<_>, QuotientError>>()?;
        let arrows = arrow_paths
            .into_iter()
            .map(|p| quotient.path(p))
            .collect::<Result<Vec<_>, _>>()?;
        let actions = vec![None; basis.len() * arrows.len()];
        Ok(Self {
            quotient,
            basis,
            indices,
            corners,
            prefixes,
            arrows,
            arrow_corners,
            actions,
        })
    }

    fn action(&mut self, index: usize, arrow: usize) -> Result<&Sparse<F>, QuotientError> {
        let slot = index * self.arrows.len() + arrow;
        if self.actions[slot].is_none() {
            let mut row = Vec::new();
            if self.corners[index].1 == self.arrow_corners[arrow].0 {
                let product = self.basis[index].multiply(&self.arrows[arrow])?;
                for term in product.representative().terms() {
                    let coordinate = self
                        .indices
                        .get(&term.path().external_representation())
                        .copied()
                        .ok_or(QuotientError::ElementOutsideCanonicalBasis)?;
                    row.push((coordinate, term.coefficient().clone()));
                }
                row.sort_by_key(|(i, _)| *i);
            }
            self.actions[slot] = Some(row);
        }
        Ok(self.actions[slot].as_ref().expect("action was initialized"))
    }

    fn apply_arrow(&mut self, row: &Sparse<F>, arrow: usize) -> Result<Sparse<F>, QuotientError> {
        let field = self.quotient.original_path_algebra().field().clone();
        let mut result = Vec::new();
        for (index, coefficient) in row {
            result = add_scaled(&result, self.action(*index, arrow)?, coefficient, &field);
        }
        Ok(result)
    }

    fn image(
        &mut self,
        left: &Sparse<F>,
        path: usize,
        cache: &mut [Option<Sparse<F>>],
    ) -> Result<Sparse<F>, QuotientError> {
        if let Some(row) = &cache[path] {
            return Ok(row.clone());
        }
        let row = if let Some((prefix, arrow)) = self.prefixes[path] {
            let previous = self.image(left, prefix, cache)?;
            self.apply_arrow(&previous, arrow)?
        } else if self.corners[left[0].0].1 == self.corners[path].0 {
            left.clone()
        } else {
            Vec::new()
        };
        cache[path] = Some(row.clone());
        Ok(row)
    }

    fn square(&mut self, rows: &[Sparse<F>]) -> Result<Vec<Sparse<F>>, QuotientError> {
        let field = self.quotient.original_path_algebra().field().clone();
        let mut span = Span::new(self.basis.len());
        let mut by_source = BTreeMap::<usize, Vec<&Sparse<F>>>::new();
        for row in rows {
            by_source
                .entry(self.corners[row[0].0].0)
                .or_default()
                .push(row);
        }
        for left in rows {
            let Some(rights) = by_source.get(&self.corners[left[0].0].1) else {
                continue;
            };
            let mut cache = vec![None; self.basis.len()];
            for right in rights {
                let mut product = Vec::new();
                for (path, coefficient) in *right {
                    let image = self.image(left, *path, &mut cache)?;
                    product = add_scaled(&product, &image, coefficient, &field);
                }
                span.insert(product, &field);
            }
        }
        Ok(span.rows)
    }
}

pub(super) fn arrow_ideal_is_nilpotent<F: Field>(
    quotient: &QuotientPathAlgebra<F>,
    basis: &[QuotientPathAlgebraElement<F>],
) -> Result<bool, QuotientError> {
    let mut products = Products::new(quotient, basis)?;
    let field = quotient.original_path_algebra().field();
    // The caller checked I subset R^2. Reduction cannot create vertex terms.
    // Thus the nonvertex nontips are exactly a basis of the arrow ideal J:
    // each is an arrow word, and every normal form of an arrow word has no
    // vertex component. This materializes the same ideal as GAP's closure,
    // without enumerating basis * arrows * basis or storing a dense identity.
    let mut current = products
        .prefixes
        .iter()
        .enumerate()
        .filter_map(|(i, prefix)| prefix.map(|_| vec![(i, field.one())]))
        .collect::<Vec<_>>();
    loop {
        let next = products.square(&current)?;
        if next.len() == current.len() {
            return Ok(next.is_empty());
        }
        debug_assert!(next.len() < current.len());
        current = next;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quotient::{PathAlgebraIdeal, quotient_ideal_span, quotient_product_space};
    use crate::{ArrowSpec, PathAlgebra, VertexIndex};
    use numerica::domains::{
        finite_field::{Z2, Zp},
        rational::Q,
    };

    fn compare<F: Field>(
        algebra: PathAlgebra<F>,
        relations: Vec<crate::PathAlgebraElement<F>>,
        expected: bool,
    ) {
        let ideal = PathAlgebraIdeal::new(algebra, relations).unwrap();
        ideal.with_high_level_groebner_basis().unwrap();
        let quotient = QuotientPathAlgebra::by_ideal(&ideal).unwrap();
        let basis = quotient.canonical_basis().unwrap().unwrap();
        let arrows = quotient
            .original_path_algebra()
            .quiver()
            .arrows()
            .map(|a| quotient.path(a.path()).unwrap())
            .collect::<Vec<_>>();
        // Freeze the previous dense implementation as a differential oracle.
        let mut old = quotient_ideal_span(&quotient, &arrows, &basis).unwrap();
        loop {
            let next = quotient_product_space(&quotient, &old, &old, &basis).unwrap();
            if next.len() == old.len() {
                assert_eq!(next.is_empty(), expected);
                break;
            }
            old = next;
        }
        assert_eq!(ideal.is_admissible().unwrap(), expected);
        assert_eq!(ideal.clone().is_admissible().unwrap(), expected);
        assert_eq!(ideal.0.admissible.get(), Some(&expected));
        // Check the cached prefix multiplication against every native product.
        let mut products = Products::new(&quotient, &basis).unwrap();
        let field = quotient.original_path_algebra().field();
        for (i, left) in basis.iter().enumerate() {
            let mut cache = vec![None; basis.len()];
            for (j, right) in basis.iter().enumerate() {
                let row = products
                    .image(&vec![(i, field.one())], j, &mut cache)
                    .unwrap();
                let mut actual = quotient.zero();
                for (k, coefficient) in row {
                    actual = actual
                        .add(&basis[k].scale_left(&coefficient).unwrap())
                        .unwrap();
                }
                assert_eq!(actual, left.multiply(right).unwrap());
            }
        }
    }

    fn loops<F: Field>(field: F) {
        let v = VertexIndex::new(1).unwrap();
        let quiver = crate::Quiver::from_indexed(
            1,
            vec![ArrowSpec::named(v, v, "x"), ArrowSpec::named(v, v, "y")],
        )
        .unwrap();
        let algebra = PathAlgebra::new(field, quiver).unwrap();
        let arrows = algebra
            .quiver()
            .arrows()
            .map(|a| algebra.element(a.path()).unwrap())
            .collect::<Vec<_>>();
        let (x, y) = (&arrows[0], &arrows[1]);
        let commute = x.multiply(y).unwrap().sub(&y.multiply(x).unwrap()).unwrap();
        compare(
            algebra.clone(),
            vec![x.pow(3).unwrap(), y.pow(3).unwrap(), commute],
            true,
        );
        // Nonhomogeneous but nilpotent: cached actions must handle reductions
        // that change ordinary path length, not just homogeneous examples.
        compare(
            algebra.clone(),
            vec![
                x.pow(2).unwrap().sub(&y.pow(3).unwrap()).unwrap(),
                x.multiply(y).unwrap(),
                y.multiply(x).unwrap(),
                y.pow(4).unwrap(),
            ],
            true,
        );
        // Finite-dimensional and all defining terms in R^2, but not admissible.
        compare(
            algebra.clone(),
            vec![
                x.pow(3).unwrap().sub(&x.pow(2).unwrap()).unwrap(),
                y.pow(2).unwrap(),
                x.multiply(y).unwrap(),
                y.multiply(x).unwrap(),
            ],
            false,
        );
    }

    #[test]
    fn sparse_admissibility_matches_dense_in_both_outcomes_and_fields() {
        loops(Q);
        loops(Z2);
        loops(Zp::new(101));
    }

    #[test]
    fn sparse_admissibility_handles_multiple_vertices_and_empty_arrow_ideal() {
        let v = VertexIndex::new(1).unwrap();
        let w = VertexIndex::new(2).unwrap();
        let quiver = crate::Quiver::from_indexed(
            2,
            vec![
                ArrowSpec::named(v, v, "a"),
                ArrowSpec::named(v, w, "b"),
                ArrowSpec::named(w, w, "c"),
            ],
        )
        .unwrap();
        let algebra = PathAlgebra::new(Q, quiver).unwrap();
        let a = algebra
            .quiver()
            .arrows()
            .map(|a| algebra.element(a.path()).unwrap())
            .collect::<Vec<_>>();
        let relations = vec![
            a[0].pow(2).unwrap(),
            a[2].pow(2).unwrap(),
            a[0].multiply(&a[1])
                .unwrap()
                .sub(&a[1].multiply(&a[2]).unwrap())
                .unwrap(),
        ];
        compare(algebra, relations, true);
        let algebra = PathAlgebra::new(Q, crate::Quiver::from_indexed(2, vec![]).unwrap()).unwrap();
        let ideal = PathAlgebraIdeal::new(algebra, vec![]).unwrap();
        let quotient = QuotientPathAlgebra::by_ideal_high_level(&ideal).unwrap();
        let basis = quotient.canonical_basis().unwrap().unwrap();
        assert!(arrow_ideal_is_nilpotent(&quotient, &basis).unwrap());
    }
}
