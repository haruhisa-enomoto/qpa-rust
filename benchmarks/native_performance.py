#!/usr/bin/env python3
"""Ordinary API cover, checked-construction and sparse/dense Hom benchmarks.

Independent fresh algebras for each cover sample; rank examples have a known
one-dimensional simultaneous centralizer over every field. Run under timeout.
"""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import statistics
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--repeat', type=int, default=3)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (4*1024**3,)*2)
    resource.setrlimit(resource.RLIMIT_CORE, (0,0))
    import qpa
    import qpa._qpa as extension
    report = dict(revision=qpa.PORT_GIT_COMMIT, repeat=args.repeat,
                  extension_sha256=hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest(),
                  runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  cases=[], status='running')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    def save(): args.output.write_text(json.dumps(report, indent=2)+'\n')
    def timed(f):
        start=time.perf_counter()
        value=f()
        return value,time.perf_counter()-start
    def record(field,case,samples,values):
        row=dict(field=field,case=case,seconds=samples,median_seconds=statistics.median(samples),
                 values=values,process_peak_rss_kib_so_far=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        report['cases'].append(row)
        save()
        print(json.dumps(row),flush=True)
    save()
    for label,field in [('Q',qpa.Field.rational()),('GF2',qpa.Field.prime(2)),
                        ('GF101',qpa.Field.prime(101)),('GF4',qpa.Field.extension(2,2))]:
        for family in ['loop','line']:
            cover_times=[]
            relation_times=[]
            for _ in range(args.repeat):
                if family=='loop':
                    algebra=qpa.QuotientAlgebra.truncated(field,qpa.Quiver.indexed(1,[(1,1,'x'),(1,1,'y')]),4)
                    simple=algebra.simple_modules()[0]
                else:
                    algebra=qpa.PathAlgebra(field,qpa.Quiver.indexed(12,[(i,i+1,f'a{i}') for i in range(1,12)]))
                    simple=algebra.simple_modules()[9]
                source=qpa.MatrixModule.direct_sum([simple]*6)
                top=source.top_projection()
                cover,elapsed=timed(source.projective_cover)
                cover_times.append(elapsed)
                assert cover.is_surjective
                assert source.projective_cover().is_identical(cover)
                dimensions=cover.source.dimension_vector
                vertices=algebra.quiver.vertices
                actions=[(name,matrix) for (s,t,name),matrix in zip(algebra.quiver.arrows,cover.source.arrow_matrices)
                         if dimensions[vertices.index(s)] and dimensions[vertices.index(t)]]
                checked,elapsed=timed(lambda:algebra.module_from_matrices(dimensions,actions))
                relation_times.append(elapsed)
                assert checked.dimension_vector==dimensions
                del top, cover, source, checked, algebra, simple
            record(label,'cover_'+family,cover_times,dimensions)
            record(label,'checked_module_'+family,relation_times,dimensions)
        n=18
        jordan=[[int(j==i+1) for j in range(n)] for i in range(n)]
        corner=[[int(i==j==0) for j in range(n)] for i in range(n)]
        for dense in [False,True]:
            matrices=[[row[:] for row in m] for m in [jordan,corner]]
            if dense:
                # Conjugate by S=I+u*v^t, with u alternating +/-1,
                # v=(1,...,1), v^t*u=0. Thus S^-1=I-u*v^t over
                # EVERY field. This produces dense integer matrices without
                # unbounded random-change-of-basis coefficient growth.
                u=[1 if i%2==0 else -1 for i in range(n)]
                for index,m in enumerate(matrices):
                    right=[[m[i][j]+sum(m[i][k]*u[k] for k in range(n))
                            for j in range(n)] for i in range(n)]
                    matrices[index]=[[right[i][j]-u[i]*sum(right[k][j] for k in range(n))
                                      for j in range(n)] for i in range(n)]
            parent=qpa.PathAlgebra(field,qpa.Quiver.indexed(1,[(1,1,'a'),(1,1,'b')]))
            scalar=field.primitive_element() if label=='GF4' else field.element(1)
            actions=[(name,[[field.element(x)*scalar for x in row] for row in m])
                     for name,m in zip(['a','b'],matrices)]
            module=parent.module_from_matrices([n],actions)
            samples=[]
            for _ in range(args.repeat):
                value,elapsed=timed(lambda:module.hom_dimension_to(module))
                assert value==1
                samples.append(elapsed)
            record(label,'hom_dense' if dense else 'hom_sparse',samples,value)
    report['status']='complete'
    save()


if __name__=='__main__':
    main()
