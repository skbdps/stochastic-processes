"""Build a read-only Codex handoff; never execute science or change scientific refs."""
from __future__ import annotations
import ast
import gzip
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import stat
import subprocess
from datetime import datetime, timezone
import zipfile

REPO = 'skbdps/stochastic-processes'
FREEZE = '0ff57919f5d8345f52d8689732ff94448ca3ca1d'
MAIN = '4e35ba23d79ec6dc16f10035df27729817e23678'
BRANCHES = {
    'completion': ('validation/week8-completion-20261031', '2a6f5dbf798171f16a23bbd2a47f7ecda378e98b'),
    'runner': ('validation/week8-runner-20261031', 'afcec2b8c82e6f4e1ba7fd58fbde56582a57e1ee'),
    'replay_diagnostic': ('validation/week8-replay-diagnostic-20261031', '1ac7b78fb0217dcf4efd138c43871fa1904b21e9'),
}
ARTIFACT = 13042844395
ARTIFACT_HASH = 'e9a780206925c4085ff4bb47cf988c7ccec66ce40db038f3fdef82703fe52513'
ARTIFACT_BYTES = 35558419
PLAN_BLOB = '2dd700ced7fdb44526e42f477e8f4bbd93b3af41'
PLAN_COMMIT = '088db0a4a05698af34b03d73aa0677b8e1033146'
FONT_SUFFIXES = {'.ttf', '.otf', '.woff', '.woff2', '.eot'}
ROOT = Path.cwd()
OUT = Path(os.environ['RUNNER_TEMP']) / 'OU_Codex_Handoff_20261031'
SOURCES: dict[str, dict] = {}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def command(args: list[str]) -> bytes:
    return subprocess.check_output(args, cwd=ROOT)


def api_bytes(path: str) -> bytes:
    return command(['gh', 'api', f'repos/{REPO}/{path}'])


def api(path: str):
    return json.loads(api_bytes(path))


def put(name: str, data: bytes, source: dict | None = None) -> None:
    relative = PurePosixPath(name)
    if relative.is_absolute() or '..' in relative.parts or '\\' in name:
        raise ValueError('unsafe output path: ' + name)
    if relative.suffix.lower() in FONT_SUFFIXES:
        raise ValueError('font files are not part of the handoff')
    target = OUT.joinpath(*relative.parts)
    if target.exists():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
    SOURCES[name] = source or {'kind': 'handoff-generated'}


def put_json(name: str, value, source: dict | None = None) -> None:
    put(name, (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n').encode(), source)


def tracked(commit: str):
    result = []
    for entry in command(['git', 'ls-tree', '-r', '-z', commit]).split(b'\0'):
        if not entry:
            continue
        header, path = entry.split(b'\t', 1)
        mode, kind, blob = header.decode().split()
        result.append((path.decode(), mode, kind, blob))
    return result


def copy_object(commit: str, path: str, mode: str, kind: str, blob: str, target: str) -> None:
    if kind != 'blob' or mode not in {'100644', '100755'}:
        raise ValueError('only ordinary committed files may be copied: ' + path)
    data = command(['git', 'cat-file', 'blob', blob])
    actual = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
    if actual != blob:
        raise ValueError('Git object integrity failure: ' + path)
    put(target, data, {'kind': 'git', 'commit': commit, 'path': path, 'blob_sha1': blob})


VERIFY_SCRIPT = r'''"""Verify all original handoff bytes without installing packages or editing files."""
from pathlib import Path
import hashlib
import json
import sys

def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / 'MANIFEST.json').read_text(encoding='utf-8'))
    errors = []
    for name, record in manifest['files'].items():
        path = root / name
        if path.is_symlink() or not path.is_file():
            errors.append(name + ': missing/nonregular file')
            continue
        data = path.read_bytes()
        if len(data) != record['bytes'] or hashlib.sha256(data).hexdigest() != record['sha256']:
            errors.append(name + ': size or SHA-256 mismatch')
    if errors:
        print('\n'.join(errors), file=sys.stderr)
        return 1
    print('PASS: ' + str(len(manifest['files'])) + ' handoff files match their recorded bytes.')
    print('This is file-integrity verification, not a new scientific audit or visual review.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
'''

SETUP_SCRIPT = r'''"""Clone all repository history into a NEW directory at the handoff's scientific freeze.

No dependencies are installed, no old checkout is modified, no experiments run,
and no remote reference is changed. Inspect newer remote work before publication.
"""
import argparse
from pathlib import Path
import subprocess
import sys

URL = 'https://github.com/skbdps/stochastic-processes.git'
FREEZE = '0ff57919f5d8345f52d8689732ff94448ca3ca1d'

def run(args):
    return subprocess.check_output(args, text=True).strip()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    destination = args.destination.expanduser().absolute()
    if destination.exists() or destination.is_symlink():
        parser.error('Destination already exists; choose a new directory. Nothing was changed.')
    subprocess.run(['git', 'clone', '-c', 'core.autocrlf=false', '--no-checkout', URL, str(destination)], check=True)
    subprocess.run(['git', '-C', str(destination), 'cat-file', '-e', FREEZE + '^{commit}'], check=True)
    live = run(['git', '-C', str(destination), 'rev-parse', 'origin/weeks5-8-sequential'])
    subprocess.run(['git', '-C', str(destination), 'switch', '-c', 'codex/finish-week8-from-handoff', FREEZE], check=True)
    assert run(['git', '-C', str(destination), 'rev-parse', 'HEAD']) == FREEZE
    print('Created:', destination)
    print('Pinned scientific freeze:', FREEZE)
    print('Research branch at clone time:', live)
    if live != FREEZE:
        print('NOTICE: newer or divergent remote work exists. Inspect it before editing or publishing.')
    print('Read START_HERE.md and the saved numerical environment before running tests or audits.')
    print('The retained full E4 artifact is separate; do not start a replacement simulation.')
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print('Git step failed. The partial new directory was preserved for inspection.', file=sys.stderr)
        raise SystemExit(exc.returncode or 1)
'''


def inspect_e4(data: bytes) -> dict:
    """Validate member bytes, saved manifests and complete raw accounting; no refits."""
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        entries = archive.infolist()
        names = [entry.filename for entry in entries]
        if len(names) != len(set(names)):
            raise ValueError('duplicate ZIP members')
        for entry in entries:
            path = PurePosixPath(entry.filename)
            if path.is_absolute() or '..' in path.parts or '\\' in entry.filename:
                raise ValueError('unsafe ZIP member')
            if stat.S_ISLNK(entry.external_attr >> 16):
                raise ValueError('ZIP symlinks are not permitted')
            if path.suffix.lower() in FONT_SUFFIXES:
                raise ValueError('unexpected font in research artifact')
        bad = archive.testzip()
        if bad:
            raise ValueError('ZIP CRC failure: ' + bad)
        def read(name): return archive.read('week8_e4/' + name)
        def obj(name): return json.loads(read(name))
        metadata = obj('run_metadata.json')
        audit = obj('audit.json')
        design = obj('resolved_design.json')
        figures = obj('figures/figure_manifest.json')
        assert metadata['status'] == 'audited_pending_visual_review'
        assert metadata['provenance']['code_commit'] == FREEZE
        assert metadata['provenance']['uncommitted_dependency_files'] == []
        assert 'week8_e4/visual_review.json' not in names
        assert audit['status'] == 'passed' and audit['full_completion_audit'] is True
        assert audit['planned_paths'] == audit['source_paths_verified'] == 14000
        assert audit['planned_fit_records'] == 42000 and audit['pipeline_errors'] == 0
        assert audit['replayed_fits'] == 42 and audit['replayed_source_paths'] == 14
        assert metadata['audit_sha256'] == sha(read('audit.json'))
        for name, digest in audit['audited_artifact_sha256'].items():
            assert sha(read(name)) == digest, name
        for name, digest in figures['figures'].items():
            assert sha(read('figures/' + name)) == digest, name
        assert len(figures['figures']) == 36
        assert len([n for n in figures['figures'] if n.endswith('.png')]) == 18
        table_manifest = obj('table_manifest.json')
        for name, record in table_manifest.items():
            assert sha(read(name)) == record['sha256'], name
        counts = {'paths': 0, 'fits': 0, 'point_valid': 0, 'interval_valid': 0, 'pipeline_errors': 0, 'partitions': 0}
        methods = design['estimators']
        assert methods == ['exact', 'pfml', 'euler']
        assert design['seed'] == 20261011 and design['replications'] == 2000 and len(design['cells']) == 7
        for cell in design['cells']:
            cid = cell['cell_id']
            manifest = obj('cells/' + cid + '/cell_manifest.json')
            assert manifest['status'] == 'complete' and manifest['completed_replications'] == 2000
            seen = set()
            for part in manifest['partitions']:
                raw = read('cells/' + cid + '/' + part['filename'])
                assert sha(raw) == part['sha256']
                rows = [json.loads(line) for line in gzip.decompress(raw).splitlines()]
                assert len(rows) == part['rows']
                counts['partitions'] += 1
                for row in rows:
                    key = (row['rep'], row['estimator'])
                    assert row['cell_id'] == cid and key not in seen
                    seen.add(key)
                    counts['fits'] += 1
                    counts['point_valid'] += int(row['point_valid'])
                    counts['interval_valid'] += int(row['interval_valid'])
                    counts['pipeline_errors'] += int(bool(row['pipeline_error']))
            assert seen == {(rep, method) for rep in range(2000) for method in methods}
            counts['paths'] += 2000
        assert counts == {'paths': 14000, 'fits': 42000, 'point_valid': 42000, 'interval_valid': 42000, 'pipeline_errors': 0, 'partitions': 56}, counts
        full_summary = json.loads(archive.read('evidence/full_summary.json'))
        exit_record = json.loads(archive.read('evidence/exit.json'))
        assert exit_record['returncode'] == 0 and archive.read('evidence/stderr.txt') == b''
        assert full_summary['status'] == metadata['status']
        assert full_summary['actual_visual_review_pending'] is True
        assert 'week8_e4-figure-renderer/renderer_state.json' in names
        for source, target in [('evidence/full_summary.json', 'evidence/E4_SAVED_FULL_SUMMARY.json'),
                               ('week8_e4/audit.json', 'evidence/E4_SAVED_NUMERICAL_AUDIT.json'),
                               ('evidence/invocation.json', 'evidence/E4_ORIGINAL_INVOCATION.json')]:
            put(target, archive.read(source), {'kind': 'actions-artifact-member', 'artifact_id': ARTIFACT, 'member': source})
        return {'kind': 'handoff_integrity_and_saved_accounting_readback',
                'checked_utc': datetime.now(timezone.utc).isoformat(),
                'artifact_id': ARTIFACT, 'original_zip_sha256': ARTIFACT_HASH,
                'zip_member_count': len(entries), 'zip_crc_check': 'passed',
                'saved_audit_status': audit['status'], 'lifecycle_status': metadata['status'],
                'raw_accounting': counts, 'audited_files_rehashed': len(audit['audited_artifact_sha256']),
                'figure_files_rehashed': len(figures['figures']), 'table_files_rehashed': len(table_manifest),
                'saved_source_array_hashes_verified': audit['source_array_hashes_verified'],
                'saved_schedule_trace_hashes_verified': audit['schedule_trace_hashes_verified'],
                'saved_endpoint_fits_replayed': audit['replayed_fits'],
                'new_sampling_or_fits': False, 'new_independent_numerical_audit': False,
                'new_visual_review': False, 'research_completion_claim': False}


def main() -> None:
    if OUT.exists():
        raise FileExistsError(OUT)
    OUT.mkdir(parents=True)
    live = api('branches?per_page=100')
    heads = {b['name']: b['commit']['sha'] for b in live}
    assert heads['main'] == MAIN and heads['weeks5-8-sequential'] == FREEZE, 'Live scientific/main heads moved; review handoff first.'
    for branch, expected in BRANCHES.values():
        assert heads[branch] == expected, (branch, heads[branch], expected)
    put_json('evidence/github/branches_at_packaging.json', live, {'kind': 'GitHub-GET', 'path': 'branches?per_page=100'})
    for name in ['START_HERE.md', 'CODEX_PROMPT.md', 'DOCUMENT_INDEX.md']:
        put(name, (ROOT/'handoff_codex'/name).read_bytes(), {'kind': 'handoff-source', 'commit': os.environ['GITHUB_SHA'], 'path': 'handoff_codex/'+name})
    for name, code in [('verify_bundle.py', VERIFY_SCRIPT), ('setup_repository.py', SETUP_SCRIPT)]:
        ast.parse(code, filename=name)
        put(name, code.encode())
    root_files = {'README.md', 'LICENSE', 'preregistration.md', 'requirements.txt', 'requirements-lock.txt', 'pytest.ini', 'run_all.py', '.gitignore'}
    prefixes = ('docs/', 'ou_irregular/', 'tests/', 'scripts/', 'configs/',
                'artifacts/week8_evidence/', 'artifacts/week8_validation/',
                'artifacts/week8_validation-figure-renderer/', 'artifacts/week8_source_validation/')
    for path, mode, kind, blob in tracked(FREEZE):
        include = path in root_files or path.startswith(prefixes)
        # Earlier checkpoint roots are orientation only; raw archives arrive with the full clone.
        if path in {f'artifacts/{folder}/{name}' for folder in ['week5_e1', 'week6_e2']
                    for name in ['run_metadata.json', 'audit.json', 'table_manifest.json', 'visual_review.json']}:
            include = True
        if include:
            copy_object(FREEZE, path, mode, kind, blob, 'repository_snapshot/'+path)
    required_docs = ['week8_gate.md', 'week8_exploratory_protocol.md', 'week8_mathematics.md', 'week8_population_targets.md',
                     'week8_implementation.md', 'week8_source.md', 'week8_execution.md', 'week8_reporting.md',
                     'week8_audit.md', 'week8_audit_progress.md', 'week8_plots.md', 'week8_runner.md', 'week8_runner_validation.md']
    assert all((OUT/'repository_snapshot/docs'/name).is_file() for name in required_docs)
    plan = command(['git', 'cat-file', 'blob', PLAN_BLOB])
    assert plan.startswith(b'%PDF')
    assert hashlib.sha1(b'blob '+str(len(plan)).encode()+b'\0'+plan).hexdigest() == PLAN_BLOB
    assert any(blob == PLAN_BLOB for _, _, _, blob in tracked(PLAN_COMMIT))
    put('original_docs/OU_master_plan_v2.pdf', plan, {'kind': 'historical-git-blob', 'blob_sha1': PLAN_BLOB, 'historical_commit': PLAN_COMMIT, 'not_restored_to_repository': True})
    for label, (branch, commit) in BRANCHES.items():
        for path, mode, kind, blob in tracked(commit):
            if path.startswith('.github/workflows/') or ('/' not in path and path.endswith('.json') and any(s in path for s in ['diagnostic', 'result', 'gate'])):
                copy_object(commit, path, mode, kind, blob, 'ci_snapshots/'+label+'/'+path)
    for run in [41198942670, 41200166912, 41200456847]:
        info = api(f'actions/runs/{run}')
        jobs = api(f'actions/runs/{run}/jobs?per_page=100')
        put_json(f'evidence/github/run_{run}.json', info, {'kind': 'GitHub-GET', 'path': f'actions/runs/{run}'})
        put_json(f'evidence/github/jobs_{run}.json', jobs, {'kind': 'GitHub-GET', 'path': f'actions/runs/{run}/jobs'})
        for job in jobs['jobs']:
            put(f'evidence/github/job_{job["id"]}.log', api_bytes(f'actions/jobs/{job["id"]}/logs'), {'kind': 'GitHub-job-log', 'run_id': run, 'job_id': job['id']})
    meta = api(f'actions/artifacts/{ARTIFACT}')
    assert meta['expired'] is False and meta['digest'] == 'sha256:'+ARTIFACT_HASH
    assert meta['size_in_bytes'] == ARTIFACT_BYTES
    data = api_bytes(f'actions/artifacts/{ARTIFACT}/zip')
    assert len(data) == ARTIFACT_BYTES and sha(data) == ARTIFACT_HASH
    put('checkpoints/week8_full_e4_41200456847.zip', data, {'kind': 'original-actions-artifact', 'artifact_id': ARTIFACT, 'run_id': 41200456847})
    put_json('evidence/E4_ARTIFACT_METADATA.json', meta, {'kind': 'GitHub-GET', 'path': f'actions/artifacts/{ARTIFACT}'})
    readback = inspect_e4(data)
    put_json('evidence/E4_READBACK_SUMMARY.json', readback)
    gate = json.loads((OUT/'repository_snapshot/artifacts/week8_evidence/validation_gate.json').read_bytes())
    assert gate['status'] == 'passed' and gate['tests_passed'] == 2055 and gate['all_tests_accounted'] == 2061 and gate['test_failures'] == 0
    after = api('branches?per_page=100')
    after_heads = {b['name']: b['commit']['sha'] for b in after}
    assert after_heads['main'] == MAIN and after_heads['weeks5-8-sequential'] == FREEZE
    report = {'created_utc': datetime.now(timezone.utc).isoformat(), 'packaging_commit': os.environ['GITHUB_SHA'],
              'packaging_run': os.environ['GITHUB_RUN_ID'], 'scientific_freeze': FREEZE, 'main_at_handoff': MAIN,
              'main_and_scientific_refs_unchanged_during_packaging': True,
              'document_files': len(list((OUT/'repository_snapshot/docs').rglob('*'))),
              'original_master_plan_included': True, 'full_E4_original_artifact_included': True,
              'new_scientific_execution': False, 'new_visual_review': False,
              'recorded_tests_passed': 2055, 'recorded_tests_accounted': 2061,
              'E4_readback': readback, 'completion_work_remaining': True}
    put_json('BUILD_REPORT.json', report)
    files = {}
    for path in sorted(OUT.rglob('*')):
        if path.is_file():
            name = path.relative_to(OUT).as_posix()
            payload = path.read_bytes()
            files[name] = {'bytes': len(payload), 'sha256': sha(payload), 'source': SOURCES[name]}
    manifest = {'schema_version': 1, 'scientific_freeze': FREEZE, 'files': files,
                'excluded_self_files': ['MANIFEST.json', 'SHA256SUMS.txt']}
    (OUT/'MANIFEST.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+'\n')
    (OUT/'SHA256SUMS.txt').write_text(''.join(record['sha256']+'  '+name+'\n' for name, record in files.items()))
    subprocess.run([os.sys.executable, str(OUT/'verify_bundle.py')], check=True)
    subprocess.run([os.sys.executable, str(OUT/'setup_repository.py'), '--help'], check=True, stdout=subprocess.DEVNULL)
    print(json.dumps({'bundle_directory': str(OUT), 'files': len(files)+2, 'payload_bytes': sum(r['bytes'] for r in files.values()), 'summary': report}, indent=2))


if __name__ == '__main__':
    main()
