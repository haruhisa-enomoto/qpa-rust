#!/usr/bin/env python3
"""Exact low-degree self-Ext for the recorded 110-dimensional ARC witness.

By default construct the OPPOSITE presentation directly with path.quotient,
from the original defining relations. No supplied basis or preconditioning is
used. --basis is an optional historical checkpoint-reproduction mode.
All construction/module/cover/kernel/Hom/Ext work is native QPA.
Published left-module matrices act on columns; transpose them for QPA right
modules over A^op. Degrees are explicit, never an unbounded Ext assertion.
Run under an external timeout; the process also has an address-space guard.
"""
import argparse
from fractions import Fraction
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import re
import resource
import time

from arc_construction import DEFAULT_INPUT, load_input, validate


WITNESS = Path('../homological-conjectures/results/auslander-reiten-counterexample/arc-expanded-witness.md')


def parse_witness(path):
    text = path.read_text()
    assert 'Field: Q.' in text
    vertices = re.search(r'^Vertices: (.*)\.$', text, re.M).group(1).split(', ')
    dims = list(map(int, re.search(r'^Dimension vector: (.*)\.$', text, re.M).group(1).split(', ')))
    blocks = re.findall(r'```text\n(.*?)\n```', text, re.S)
    assert len(blocks) == 1
    matrices, record = [], None
    for line in blocks[0].splitlines():
        if record is None:
            match = re.fullmatch(r'(\w+): (\w+) -> (\w+); shape (\d+) x (\d+)', line)
            assert match, line
            name, source, target, n, m = match.groups()
            record = dict(name=name, source=source, target=target, shape=[int(n), int(m)], entries=[])
        elif line == 'end':
            matrices.append(record)
            record = None
        else:
            r,c,value = line.split()
            record['entries'].append((int(r),int(c),Fraction(value)))
    assert record is None and len(matrices)==43 and dims == [26,30,40,13,0,0,0,1]
    assert sum(len(r['entries']) for r in matrices)==439
    return vertices, dims, matrices


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--witness', type=Path, default=WITNESS)
    parser.add_argument('--basis', type=Path,
                        help='optional completed-basis checkpoint; omitted for the ordinary API')
    parser.add_argument('--degree', type=int, choices=[0,1,2], default=1,
                        help='0 constructs/checks the module only; 1 or 2 computes through this degree')
    parser.add_argument('--profile-cover', action='store_true',
                        help='profile the reference cover assembly through public APIs; omit Ext queries')
    parser.add_argument('--memory-mib', type=int, default=8192)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (args.memory_mib*1024**2,)*2)
    resource.setrlimit(resource.RLIMIT_CORE, (0,0))
    import qpa
    import qpa._qpa as extension
    report = dict(field='Q', convention='right module over A^op, transposed published left actions',
                  pairs=[['Z','Z']], requested_degrees=list(range(1,args.degree+1)),
                  version=version('qpa-rust'), revision=qpa.PORT_GIT_COMMIT,
                  extension_sha256=hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest(),
                  construction='checkpoint' if args.basis else 'ordinary_api',
                  profile_cover=args.profile_cover,
                  inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [args.input,args.witness,args.basis] if p is not None},
                  events=[], ext_dimensions={}, status='running')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    def event(name, **values):
        value = dict(event=name, **values)
        report['events'].append(value)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(value),flush=True)
    def timed(stage, operation):
        event('start',stage=stage)
        start=time.perf_counter()
        value=operation()
        event('timing',stage=stage,seconds=time.perf_counter()-start,
              peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        return value
    data=load_input(args.input)
    validate(data)
    field=qpa.Field.rational()
    path=qpa.PathAlgebra(field,qpa.Quiver(data['vertices'],[(t,s,n) for s,t,n in data['arrows']]))
    generators={n:path.generator(n) for _,_,n in data['arrows']}
    def polynomial(relation, reverse=False):
        value=path.zero()
        for word,c in relation:
            word=list(reversed(word)) if reverse else word
            term=generators[word[0]]
            for n in word[1:]:term*=generators[n]
            c=Fraction(c)
            value+=term.scale_left(field.element(c.numerator,c.denominator))
        return value
    relations=[polynomial(r,True) for r in data['internal_relations']+data['connecting_relations']]
    if args.basis:
        saved=json.loads(args.basis.read_text())
        assert saved['opposite'] and saved['field']=='Q'
        assert saved['input_sha256']==report['inputs'][str(args.input)]
        complete=timed('load_basis',lambda:[polynomial(r) for r in saved['relations']])
        ideal=path.ideal(relations)
        ideal.with_complete_groebner_basis(complete)
        algebra=timed('quotient',ideal.quotient)
    else:
        algebra=timed('ordinary_quotient',lambda:path.quotient(relations))
    assert timed('verify_input_relations',lambda:all(
        algebra.groebner_basis.remainder(r).is_zero for r in relations))
    assert int(algebra.dimension())==7708
    assert len(algebra.canonical_basis())==7708
    assert timed('admissibility',lambda:algebra.is_admissible)
    def construct():
        vertices,dims,records=parse_witness(args.witness)
        assert vertices==data['vertices']
        actions=[]
        for (s,t,name),r in zip(data['arrows'],records,strict=True):
            assert (s,t,name)==(r['source'],r['target'],r['name'])
            n,m=r['shape']
            assert [n,m]==[dims[vertices.index(s)],dims[vertices.index(t)]]
            matrix=[[field.element(0) for _ in range(max(1,n))] for _ in range(max(1,m))]
            seen=set()
            for row,col,c in r['entries']:
                assert 1<=row<=n and 1<=col<=m and (row,col) not in seen
                seen.add((row,col))
                matrix[col-1][row-1]=field.element(c.numerator,c.denominator)
            if r['entries']:actions.append((name,matrix))
        return algebra.module_from_matrices(dims,actions)
    witness=timed('checked_witness_module',construct)
    assert witness.dimension==110
    event('module',dimension=witness.dimension,dimension_vector=witness.dimension_vector)
    source=witness
    retained=[]
    for degree in range(1,args.degree+1):
        # ProjectiveCover needs this map anyway. Retain its weak cache and
        # time it separately, rather than hiding top construction in the cover.
        top_projection=timed(f'top_projection_{degree}',source.top_projection)
        retained.append(top_projection)
        event('source_top',degree=degree,dimension_vector=top_projection.range.dimension_vector)
        if args.profile_cover:
            # Reference translation, not an alternative optimized workflow.
            generators=timed(f'generators_{degree}',source.minimal_generating_set)
            family=timed(f'projective_family_{degree}',algebra.indecomposable_projectives)
            maps=timed(f'generator_maps_{degree}',lambda:[g.hom_from_projective() for g in generators])
            direct_sum=timed(f'direct_sum_{degree}',lambda:qpa.MatrixModule.direct_sum([m.source for m in maps]))
            projections=timed(f'projections_{degree}',direct_sum.direct_sum_projections)
            composites=timed(f'composites_{degree}',lambda:[p.compose(m) for p,m in zip(projections,maps,strict=True)])
            def add_maps():
                result=composites[0]
                for m in composites[1:]: result=result.add(m)
                return result
            cover=timed(f'add_maps_{degree}',add_maps)
            del family, maps, projections, composites, direct_sum, generators
        else:
            cover=timed(f'cover_{degree}',source.projective_cover)
        inclusion=timed(f'kernel_{degree}',cover.kernel_inclusion)
        retained.extend([cover,inclusion])
        def check_exactness():
            assert cover.is_surjective and inclusion.is_injective
            assert inclusion.compose(cover).is_zero
            assert cover.source.dimension == source.dimension + inclusion.source.dimension
            return True
        assert timed(f'exactness_{degree}',check_exactness)
        event('resolution_term',degree=degree,exact_sequence_checked=True,projective_dimension=cover.source.dimension,
              projective_dimension_vector=cover.source.dimension_vector,
              syzygy_dimension=inclusion.source.dimension,
              syzygy_dimension_vector=inclusion.source.dimension_vector)
        if not args.profile_cover:
            value=timed(f'ext_{degree}',lambda:source.ext1_dimension(witness))
            report['ext_dimensions'][str(degree)]=value
            event('ext',degree=degree,source='Z',target='Z',dimension=value)
        source=inclusion.source
    report['status']='complete'
    event('result',ext_dimensions=report['ext_dimensions'])


if __name__=='__main__':
    main()
