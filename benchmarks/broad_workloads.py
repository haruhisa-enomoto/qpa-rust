#!/usr/bin/env python3
"""Additional consumer families and scaling probes for the release audit.

Uses the timing/result comparison harness in module_workloads.py. Every case
constructs fresh inputs. Existing consumer code is loaded without editing it.
Requires the consumer's pytest, sympy, numpy and tqdm dependencies. Run with
--consumer-root ../homological-conjectures --repeat 5 --output target/broad.json.
"""
from fractions import Fraction
import runpy
import sys

import qpa
import module_workloads


original_consumer_cases = module_workloads.consumer_cases


def consumer_cases(root):
    cases = original_consumer_cases(root)
    sys.path.insert(0, str(root))
    sys.path.insert(0, str(root / 'lib/python'))
    from lib.ar_knitting_qpa import knit_indecomposables
    from lib.finite_type_g_fan_qpa import finite_type_g_fan
    from lib.natural_endomorphism_qpa import natural_endomorphism_blocks, self_ext_diagnostic
    from lib.ringel_zhang_qpa import build_ringel_zhang_template
    from lib.tensor_words_qpa import verify_word_family
    from lib.cartan_rank_two.qpa import gf2_bipartite_qpa_quotient
    from lib.syzygy_cycles_qpa import bounded_simple_syzygy_screen

    def line(size, field):
        return qpa.PathAlgebra(field, qpa.Quiver.indexed(
            size, [(i, i + 1, f'a{i}') for i in range(1, size)]))

    def knitting(size, field):
        result = knit_indecomposables(line(size, field), expected_count=size * (size + 1) // 2)
        assert result.complete and result.expected_count_matches
        return sorted(list(m.dimension_vector) for m in result.modules)

    def fan():
        result = finite_type_g_fan(line(4, qpa.Field.prime(5)))
        assert len(result.modules) == 10
        return result.fan

    def natural():
        p = qpa.PathAlgebra(qpa.Field.prime(2), qpa.Quiver.indexed(2, [(1, 2, 'a'), (2, 1, 'b')]))
        algebra = p.quotient([p.generator('a') * p.generator('b')])
        projective = algebra.indecomposable_projectives()[1]
        blocks = natural_endomorphism_blocks(qpa.MatrixModule.direct_sum([projective] * 3))
        results = [self_ext_diagnostic(block, 4) for block in blocks]
        assert len(blocks) == 1 and [r['dimension'] for r in results[0]['self_ext']] == [1] * 4
        return [block.dimension for block in blocks], [r['dimension'] for r in results[0]['self_ext']]

    def rational_syzygies():
        parameter = Fraction(3, 2)
        template = build_ringel_zhang_template(parameter)
        module = template.module(parameter)
        regular = template.regular_module()
        dimensions = []
        for _ in range(8):
            dimensions.append(module.ext1(regular).dimension)
            module = module.first_syzygy()
            parameter *= Fraction(3, 2)
            assert module.isomorphic_to(template.module(parameter))
        assert dimensions == [0] * 8
        return dimensions

    presilting = runpy.run_path(str(root / 'tests/test_presilting_window_qpa.py'))
    strings = runpy.run_path(str(root / 'tests/test_string_fan_realization_qpa.py'))
    mesh = runpy.run_path(str(root / 'tests/test_finite_tau_mesh_qpa.py'))
    relative = runpy.run_path(str(root / 'tests/test_relative_generator_ext_qpa.py'))
    cycles = runpy.run_path(str(root / 'tests/test_cartan_rank_two_module_cycles.py'))
    lattice = runpy.run_path(str(root / 'lib/python/subcat/test_qpa_lattice.py'))

    def tensor_words():
        path = line(2, qpa.Field.prime(2))
        algebra = qpa.QuiverAlgebra.from_path_algebra(path)
        projectives = path.indecomposable_projectives()
        report = verify_word_family(algebra, dict(enumerate(projectives)),
                                    [(i, j) for i in range(2) for j in range(2)])
        assert report == {'rank': 4, 'summands': 4, 'checks': [], 'witness': False}
        return report

    def deep_recurrence():
        # The depth-four fixture in test_cartan_rank_two_candidate_verifier.
        relations = (
            ((0, 2, 32), (0, 4, 1)), ((0, 2, 16),), ((0, 2, 8),), ((0, 2, 6),),
            ((1, 2, 32),), ((1, 2, 16),), ((1, 2, 2), (1, 4, 262144)),
            ((0, 3, 68),), ((0, 5, 1),), ((1, 5, 68719476736),),
        )
        algebra = gf2_bipartite_qpa_quotient(relations, (3, 2))
        report = bounded_simple_syzygy_screen(algebra, max_depth=4,
                                             source_dimension_cap=64, target_dimension_cap=256)
        assert report['status'] == 'recurrence_certified'
        assert report['recurrence']['source_dimension_vector'] == [3, 5]
        assert report['recurrence']['target_dimension_vector'] == [23, 30]
        return report

    def window():
        host = presilting['host'].__wrapped__()
        presilting['test_window_enumeration_respects_silting_rank'](host)

    cases.update(
        consumer_ar_a5_binary=lambda: knitting(5, qpa.Field.prime(2)),
        consumer_ar_a7_prime=lambda: knitting(7, qpa.Field.prime(5)),
        consumer_ar_a5_extension=lambda: knitting(5, qpa.Field.extension(2, 2)),
        consumer_fan_a4=fan,
        consumer_presilting_window=window,
        consumer_string_fan=strings['test_representation_infinite_string_algebra_with_and_without_pruning'],
        consumer_natural_endomorphism=natural,
        consumer_rational_syzygies=rational_syzygies,
        consumer_mesh=mesh['test_nonzero_stable_socle_correction'],
        consumer_relative_ext=lambda: relative['test_regular_generator_agrees_with_ordinary_ext'](5),
        consumer_cartan_cycle=cycles['test_coordinate_module_cycle_is_replayed_and_serializable'],
        consumer_cartan_deep_recurrence=deep_recurrence,
        consumer_closure_lattice=lattice['test_peak_a3_collective_generation_fixture'],
        consumer_tensor_words=tensor_words,
    )

    def hom(field):
        module = module_workloads.kronecker(field, 24)
        result = module.hom_dimension_to(module)
        assert result == 24
        return result

    def dense_composition(field):
        size = 32
        algebra = qpa.PathAlgebra(field, qpa.Quiver.indexed(1, []))
        module = algebra.module_from_integer_matrices([size], [])
        values = [field.element(i) for i in range(1, 5)]
        left = [[values[(i * 3 + j) % 4] for j in range(size)] for i in range(size)]
        right = [[values[(i + j * 3) % 4] for j in range(size)] for i in range(size)]
        f = module.homomorphism_to(module, [left])
        g = module.homomorphism_to(module, [right])
        for _ in range(8):
            product = f.compose(g)
        return [[str(entry) for entry in row] for row in product.vertex_matrices[0]]

    cases['hom_prime'] = lambda: hom(qpa.Field.prime(5))
    cases['hom_extension'] = lambda: hom(qpa.Field.extension(2, 2))
    cases['dense_compose_rational'] = lambda: dense_composition(qpa.Field.rational())
    cases['dense_compose_prime'] = lambda: dense_composition(qpa.Field.prime(5))
    cases['dense_compose_binary'] = lambda: dense_composition(qpa.Field.prime(2))
    return cases


if __name__ == '__main__':
    module_workloads.consumer_cases = consumer_cases
    module_workloads.main()
