use super::*;
use crate::{PathAlgebra, Quiver};
use numerica::domains::rational::Q;

#[test]
fn selective_cache_is_weak_and_preserves_live_members_on_full_merge() {
    let algebra =
        QuiverAlgebra::Path(PathAlgebra::new(Q, Quiver::from_indexed(2, vec![]).unwrap()).unwrap());
    let first = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let second = right_module_over_path_algebra(algebra.clone(), vec![0, 1], &[])
        .unwrap()
        .unwrap();
    let weak = first.downgrade();
    assert!(
        algebra
            .merge_indec_projective_module(0, &first)
            .unwrap()
            .is_identical(&first)
    );
    assert!(algebra.cached_indec_projective_modules().is_none());
    assert!(algebra.merge_indec_projective_module(2, &second).is_none());
    let replacement = right_module_over_path_algebra(algebra.clone(), vec![1, 0], &[])
        .unwrap()
        .unwrap();
    let family = algebra.merge_indec_projective_modules(&[replacement, second.clone()]);
    assert!(family[0].is_identical(&first));
    drop(first);
    assert!(weak.upgrade().is_some());
    drop(family);
    assert!(weak.upgrade().is_none());
    assert!(algebra.cached_indec_projective_module(0).is_none());
    assert!(
        algebra
            .cached_indec_projective_module(1)
            .unwrap()
            .is_identical(&second)
    );
    let weak_second = second.downgrade();
    drop(second);
    assert!(weak_second.upgrade().is_none());
    assert!(algebra.cached_indec_projective_modules().is_none());
}
