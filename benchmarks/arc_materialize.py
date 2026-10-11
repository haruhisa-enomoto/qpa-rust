#!/usr/bin/env python3
"""Measure ordinary QPA construction of the ARC input under an external timeout.

Pass every defining relation directly to PathAlgebra.quotient. No manual
partition, preconditioning, or supplied Groebner basis is used. Any structural
optimization is native and applies to arbitrary quivers/presentations.
"""
import argparse
from fractions import Fraction
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import resource
import time

from arc_construction import DEFAULT_INPUT, load_input, validate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--opposite', action='store_true')
    parser.add_argument('--save-basis', type=Path)
    parser.add_argument('--admissibility', action='store_true')
    parser.add_argument('--memory-mib', type=int, default=4096)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (args.memory_mib * 1024**2,) * 2)
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    import qpa
    import qpa._qpa as extension

    def log(event, **data):
        print(json.dumps(dict(event=event, **data)), flush=True)

    def timed(stage, operation):
        log('start', stage=stage)
        start = time.perf_counter()
        result = operation()
        log('timing', stage=stage, seconds=time.perf_counter() - start,
            peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        return result

    source_hash = hashlib.sha256(args.input.read_bytes()).hexdigest()
    log('metadata', construction='ordinary_api', input_sha256=source_hash,
        version=version('qpa-rust'), revision=qpa.PORT_GIT_COMMIT,
        opposite=args.opposite,
        extension_sha256=hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest())
    data = load_input(args.input)
    log('input', **validate(data))
    field = qpa.Field.rational()
    arrows = [(t, s, n) if args.opposite else (s, t, n) for s, t, n in data['arrows']]
    parent = qpa.PathAlgebra(field, qpa.Quiver(data['vertices'], arrows))
    generators = {n: parent.generator(n) for _, _, n in arrows}

    def polynomial(relation):
        result = parent.zero()
        for names, coefficient in relation:
            names = list(reversed(names)) if args.opposite else names
            term = generators[names[0]]
            for name in names[1:]:
                term *= generators[name]
            c = Fraction(coefficient)
            result += term.scale_left(field.element(c.numerator, c.denominator))
        return result

    relations = timed('presentation', lambda: [polynomial(r) for r in
                      data['internal_relations'] + data['connecting_relations']])
    quotient = timed('ordinary_quotient', lambda: parent.quotient(relations))
    basis = quotient.groebner_basis
    assert timed('original_relations_reduce_to_zero',
                 lambda: all(basis.remainder(r).is_zero for r in relations))
    assert int(timed('dimension', quotient.dimension)) == 7708
    assert len(timed('canonical_basis', quotient.canonical_basis)) == 7708
    if args.admissibility:
        assert timed('admissibility', lambda: quotient.is_admissible)
    if args.save_basis:
        payload = dict(input_sha256=source_hash, opposite=args.opposite, field='Q',
                       relations=[[[[a.name for a in p.walk], str(c)] for p, c in r.terms]
                                  for r in quotient.relations])
        args.save_basis.parent.mkdir(parents=True, exist_ok=True)
        args.save_basis.write_text(json.dumps(payload) + '\n')
    log('result', dimension=7708, admissibility_checked=args.admissibility,
        completed_relations=len(quotient.relations))


if __name__ == '__main__':
    main()
