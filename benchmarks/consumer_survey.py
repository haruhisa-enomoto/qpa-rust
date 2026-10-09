#!/usr/bin/env python3
"""Time every consumer QPA test, optionally profiling each test call.

Run with this repository's release extension and pytest installed. Fixtures and
assertions are supplied by the consumer; no consumer files are changed. Setup,
call and teardown are recorded separately. Profiling is diagnostic: compare
wall times only between runs without --profile. Run in fresh processes for
independent repetitions; module/class fixtures may be shared within a run.
For the consumer's Cartan verifier, pass --qpa-checkout-root to point its
origin check at the checkout under test (also applied in spawned workers).
"""
import argparse
import cProfile
import hashlib
import json
import os
from pathlib import Path
import platform
import pstats
import subprocess
import sys
import tempfile

import pytest
import qpa
import qpa._qpa as extension


class Survey:
    def __init__(self, profile):
        self.profile = profile
        self.cases = {}

    @pytest.hookimpl(hookwrapper=True)
    def pytest_runtest_call(self, item):
        profiler = cProfile.Profile() if self.profile else None
        if profiler:
            profiler.enable()
        yield
        if profiler:
            profiler.disable()
            stats = pstats.Stats(profiler).stats
            rows = [dict(function=f"{key[0]}:{key[1]}:{key[2]}", calls=value[1],
                         self_seconds=value[2], cumulative_seconds=value[3])
                    for key, value in stats.items()]
            self.cases.setdefault(item.nodeid, {})['profile'] = sorted(
                rows, key=lambda row: row['cumulative_seconds'], reverse=True)[:25]

    def pytest_runtest_logreport(self, report):
        case = self.cases.setdefault(report.nodeid, {})
        case[report.when] = dict(seconds=report.duration, outcome=report.outcome)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--consumer-root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--profile', action='store_true')
    parser.add_argument('--extra-test', action='append', default=[],
                        help='Additional test file relative to the consumer root')
    parser.add_argument('--qpa-checkout-root', type=Path,
                        help='Configure the consumer provenance check for this checkout, including spawned workers')
    args = parser.parse_args()
    root = args.consumer_root.resolve()
    if args.qpa_checkout_root:
        # The consumer normally checks its sibling QPA_rust_python checkout.
        # Point that same origin check at the explicitly selected benchmark
        # checkout, without changing consumer files or falsifying module paths.
        with tempfile.TemporaryDirectory(prefix='qpa-survey-') as directory:
            Path(directory, 'sitecustomize.py').write_text(
                'import os\nfrom pathlib import Path\n'
                'import lib.qpa_provenance as provenance\n'
                'provenance.QPA_CHECKOUT_ROOT = Path(os.environ["QPA_BENCHMARK_CHECKOUT_ROOT"])\n')
            env = dict(os.environ)
            env['QPA_BENCHMARK_CHECKOUT_ROOT'] = str(args.qpa_checkout_root.resolve())
            env['PYTHONPATH'] = os.pathsep.join([directory, str(root), env.get('PYTHONPATH', '')])
            argv = ['--consumer-root', str(root), '--output', str(args.output.resolve())]
            if args.profile:
                argv.append('--profile')
            for path in args.extra_test:
                argv.extend(['--extra-test', path])
            return subprocess.call([sys.executable, str(Path(__file__).resolve()), *argv], env=env)
    sys.path.insert(0, str(root))
    sys.path.insert(0, str(root / 'lib/python'))
    paths = sorted((root / 'tests').glob('*qpa.py'))
    paths += [root / path for path in args.extra_test]
    survey = Survey(args.profile)
    status = pytest.main(['-q', '-p', 'no:cacheprovider', *map(str, paths)], plugins=[survey])
    report = dict(python=sys.version, platform=platform.platform(),
                  qpa_build_revision=qpa.PORT_GIT_COMMIT,
                  extension_sha256=hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest(),
                  cpu_affinity=sorted(os.sched_getaffinity(0)), profiled=args.profile,
                  pytest_exit_code=int(status), cases=survey.cases,
                  provenance_checkout_override=os.environ.get('QPA_BENCHMARK_CHECKOUT_ROOT'),
                  test_source_sha256={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                      for p in paths})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    return int(status)


if __name__ == '__main__':
    raise SystemExit(main())
