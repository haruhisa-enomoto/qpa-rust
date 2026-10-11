#!/usr/bin/env python3
"""Guarded stage diagnosis of the recorded ARC presentation, never Ext.

Inputs are read-only. The default input is the tracked expanded bound-quiver
list; --compare-json optionally checks it against the saved structured input.
The dimension 7708 is an inherited expectation, not an assumption of completion.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import resource
import signal
import subprocess
import sys
import time


DEFAULT_INPUT = Path('../homological-conjectures/results/auslander-reiten-counterexample/arc-expanded-bound-quiver.md')


def load_input(path):
    if path.suffix == '.json':
        data = json.loads(path.read_text())
        return {key: data[key] for key in ['vertices', 'arrows', 'internal_relations', 'connecting_relations']}
    text = path.read_text()
    vertices = re.search(r'^Vertices: (.*)\.$', text, re.M).group(1).split(', ')
    blocks = re.findall(r'```text\n(.*?)\n```', text, re.S)
    assert len(blocks) == 2
    arrows = []
    for line in blocks[0].splitlines():
        match = re.fullmatch(r'(\w+): (\w+) -> (\w+)', line)
        assert match, line
        name, source, target = match.groups()
        arrows.append([source, target, name])
    relations = []
    for line in blocks[1].splitlines():
        relation = []
        for term in line.replace(' - ', ' + -').split(' + '):
            sign = -1 if term.startswith('-') else 1
            parts = term.removeprefix('-').split('*')
            coefficient = sign * int(parts.pop(0)) if parts[0].isdigit() else sign
            relation.append([parts, str(coefficient)])
        relations.append(relation)
    return dict(vertices=vertices, arrows=arrows, internal_relations=relations[:62], connecting_relations=relations[62:])


def validate(data):
    vertices, arrows = data['vertices'], data['arrows']
    relations = data['internal_relations'] + data['connecting_relations']
    assert (len(vertices), len(arrows), len(relations)) == (8, 43, 106)
    assert len(set(vertices)) == 8
    endpoints = {name: (s, t) for s, t, name in arrows}
    assert len(endpoints) == 43
    degrees = Counter()
    for relation in relations:
        assert relation
        corners = set()
        for word, coefficient in relation:
            assert Fraction(coefficient) != 0 and len(word) >= 2
            assert all(name in endpoints for name in word)
            assert all(endpoints[a][1] == endpoints[b][0] for a, b in zip(word, word[1:]))
            corners.add((endpoints[word[0]][0], endpoints[word[-1]][1]))
            degrees[len(word)] += 1
        assert len(corners) == 1
    assert all(s in vertices and t in vertices for s, t, _ in arrows)
    return {'vertices': len(vertices), 'arrows': len(arrows), 'relations': len(relations),
            'terms': sum(degrees.values()), 'terms_by_path_length': dict(sorted(degrees.items())),
            'max_terms_per_relation': max(map(len, relations)), 'expected_dimension': 7708}


def emit(event, **data):
    print('BENCH ' + json.dumps(dict(event=event, monotonic=time.monotonic(), **data)), flush=True)


def worker(args):
    cap = args.memory_mib * 1024 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    import qpa
    import qpa._qpa as extension
    from importlib.metadata import version
    emit('metadata', pid=os.getpid(), package=qpa.__file__, version=version('qpa-rust'),
         revision=qpa.PORT_GIT_COMMIT,
         extension_sha256=hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest())
    def timed(stage, operation):
        emit('start', stage=stage)
        start = time.perf_counter()
        value = operation()
        emit('timing', stage=stage, seconds=time.perf_counter()-start,
             cpu_seconds=time.process_time(), peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        return value
    data = timed('parse', lambda: load_input(args.input))
    emit('input', **validate(data))
    if args.compare_json:
        assert data == load_input(args.compare_json), 'Markdown/JSON presentations differ'
        emit('input_match', path=str(args.compare_json))
    def setup():
        field = qpa.Field.rational() if args.field == 'Q' else qpa.Field.prime(101)
        arrows = [(t, s, name) if args.opposite else (s, t, name) for s, t, name in data['arrows']]
        algebra = qpa.PathAlgebra(field, qpa.Quiver(data['vertices'], arrows))
        generators = {name: algebra.generator(name) for _, _, name in arrows}
        relations = []
        for relation in data['internal_relations'] + data['connecting_relations']:
            value = algebra.zero()
            for names, coefficient in relation:
                word = list(reversed(names)) if args.opposite else names
                term = generators[word[0]]
                for name in word[1:]:
                    term = term * generators[name]
                coefficient = Fraction(coefficient)
                value += term.scale_left(field.element(coefficient.numerator, coefficient.denominator))
            relations.append(value)
        return algebra, relations
    algebra, relations = timed('presentation', setup)
    ideal = timed('ideal', lambda: algebra.ideal(relations))
    complete = timed('groebner', ideal.gbnp_groebner_basis)
    emit('completed_basis', relations=len(complete))
    if args.stop_after == 'groebner':
        emit('result', completed_through='groebner', gb_relations=len(complete))
        return
    ideal = timed('completed_ideal', lambda: algebra.ideal(complete))
    timed('attach', lambda: ideal.with_complete_groebner_basis(complete))
    quotient = timed('quotient', ideal.quotient)
    # Intentionally probe these before admissibility so that its product-space
    # work cannot hide successful completion/dimension/basis stages.
    dimension = int(timed('dimension', quotient.dimension))
    emit('dimension', value=dimension, matches_expected=(dimension == 7708))
    basis = timed('basis', quotient.canonical_basis)
    assert len(basis) == dimension
    emit('basis_size', value=len(basis))
    if args.stop_after == 'admissibility':
        assert timed('admissibility', ideal.is_admissible)
    emit('result', dimension=dimension, completed_through=args.stop_after)


def process_usage(pid):
    """Sample the root and descendants; time-series RSS is a sampled lower bound."""
    pids, rss, cpu = [pid], 0, 0.0
    for child in pids:
        try:
            pids.extend(map(int, Path(f'/proc/{child}/task/{child}/children').read_text().split()))
            status = Path(f'/proc/{child}/status').read_text()
            match = re.search(r'^VmRSS:\s+(\d+)', status, re.M)
            rss += int(match.group(1)) if match else 0
            stat = Path(f'/proc/{child}/stat').read_text().rsplit(')', 1)[1].split()
            cpu += (int(stat[11]) + int(stat[12])) / os.sysconf('SC_CLK_TCK')
        except (FileNotFoundError, ProcessLookupError):
            pass
    return {'pids': pids, 'rss_kib': rss, 'cpu_seconds': cpu}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--compare-json', type=Path)
    parser.add_argument('--field', choices=['Q', 'GF101'], default='Q')
    parser.add_argument('--opposite', action='store_true')
    parser.add_argument('--stop-after', choices=['groebner', 'basis', 'admissibility'], default='admissibility')
    parser.add_argument('--stage-timeout', type=float, default=120)
    parser.add_argument('--total-timeout', type=float, default=600)
    parser.add_argument('--memory-mib', type=int, default=4096)
    parser.add_argument('--profile', action='store_true', help='native py-spy sampling at 49 Hz (separate diagnostic run)')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if not args.output or min(args.stage_timeout, args.total_timeout, args.memory_mib) <= 0:
        parser.error('--output and positive limits are required')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stdout_path = args.output.with_suffix('.stdout.log')
    stderr_path = args.output.with_suffix('.stderr.log')
    command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--input', str(args.input.resolve()),
               '--field', args.field, '--stop-after', args.stop_after, '--memory-mib', str(args.memory_mib)]
    if args.compare_json:
        command += ['--compare-json', str(args.compare_json.resolve())]
    if args.opposite:
        command += ['--opposite']
    profile_path = args.output.with_suffix('.speedscope.json')
    if args.profile:
        command = ['py-spy', 'record', '--native', '--rate', '49', '--format', 'speedscope',
                   '--output', str(profile_path), '--', *command]
    report = {'date': datetime.now(timezone.utc).isoformat(), 'platform': platform.platform(),
              'python': sys.version, 'cpu_affinity': sorted(os.sched_getaffinity(0)),
              'input': str(args.input), 'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'field': args.field, 'opposite': args.opposite, 'stop_after': args.stop_after,
              'stage_timeout': args.stage_timeout, 'total_timeout': args.total_timeout,
              'address_space_limit_mib': args.memory_mib, 'command': command,
              'profile': str(profile_path) if args.profile else None, 'status': 'running', 'samples': []}
    start = time.monotonic()
    with stdout_path.open('w') as out, stderr_path.open('w') as err:
        process = subprocess.Popen(command, stdout=out, stderr=err, start_new_session=True)
        try:
            while True:
                now = time.monotonic()
                events = [json.loads(line[6:]) for line in stdout_path.read_text().splitlines()
                          if line.startswith('BENCH ') and line.endswith('}')]
                started = [e for e in events if e['event'] == 'start']
                current = started[-1] if started else {'stage': 'startup', 'monotonic': start}
                usage = process_usage(process.pid)
                report['samples'].append(dict(elapsed_seconds=now-start, stage=current['stage'], **usage))
                report.update(events=events, elapsed_seconds=now-start, last_started_stage=current['stage'],
                              stage_seconds={e['stage']: e['seconds'] for e in events if e['event'] == 'timing'})
                if process.poll() is not None:
                    report.update(returncode=process.returncode,
                                  status='ok' if process.returncode == 0 and any(e['event']=='result' for e in events) else 'error')
                elif now-start > args.total_timeout or now-current['monotonic'] > args.stage_timeout:
                    report.update(status='timeout', timed_out_stage=current['stage'])
                args.output.write_text(json.dumps(report, indent=2) + '\n')
                if report['status'] != 'running':
                    break
                time.sleep(1)
        finally:
            # SIGINT lets py-spy save its profile; always kill any remaining
            # descendants afterwards, including on coordinator interruption.
            if args.profile and process.poll() is None:
                os.killpg(process.pid, signal.SIGINT)
                try:
                    process.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    pass
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
    report['returncode'] = process.returncode
    report['stderr'] = stderr_path.read_text()
    report['sampled_peak_process_tree_rss_kib'] = max(s['rss_kib'] for s in report['samples'])
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['status', 'last_started_stage', 'stage_seconds', 'elapsed_seconds']}, indent=2))


if __name__ == '__main__':
    main()
