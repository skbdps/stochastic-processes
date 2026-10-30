"""Independent runner lifecycle faults on explicitly literal evidence.

Source records and audit outcomes are labelled test inputs, not simulated OU
observations or numerical completion evidence. Persistence, table generation,
saved-audit binding, rendering, lifecycle transitions and seals are exercised
through the actual components. Every scientific generation/replay entry point
is blocked. Review dictionaries are literal stage inputs, not visual reviews.
"""
from contextlib import contextmanager
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from ou_irregular.experiments import audit as historical_audit
from ou_irregular.experiments import execution as historical_execution
from ou_irregular.experiments.provenance import canonical, file_digest
from ou_irregular.inference.adapters import Week3Estimator
from ou_irregular.week4.simulation import OuSampleSource
from ou_irregular.week7 import execution, floor_audit, path_diagnostics, plots, reporting, runner
from ou_irregular.week7.floor_source import PairedFloorSource

from .test_runner_reference import literal_records, expected_counts


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()


def read_json(path):
    return json.loads(Path(path).read_bytes())


def write_json(path, value):
    Path(path).write_bytes(encoded(value))


def snapshot(output):
    output = Path(output)
    if not output.exists():
        return None
    return {str(path.relative_to(output)): (('link', str(path.readlink())) if path.is_symlink()
            else ('file', hashlib.sha256(path.read_bytes()).hexdigest()) if path.is_file()
            else ('directory',)) for path in sorted(output.rglob('*'))}


def forbidden(*args, **kwargs):
    raise AssertionError('literal runner integrity tests must not simulate or replay science')


@contextmanager
def block_science():
    with pytest.MonkeyPatch.context() as patch:
        for obj, name in ((OuSampleSource, 'sample'), (PairedFloorSource, 'sample'),
                          (Week3Estimator, 'fit'), (historical_execution, 'source_records'),
                          (historical_audit, 'audit_run'), (floor_audit, 'audit_run'),
                          (floor_audit, 'independent_source'), (floor_audit, 'source_records'),
                          (path_diagnostics, 'run_selected'), (path_diagnostics, '_one_selected')):
            patch.setattr(obj, name, forbidden)
        yield patch


@pytest.fixture(autouse=True)
def no_science():
    with block_science():
        yield


def audit_outcome(output, design, kind, rows):
    """An explicitly literal stage result with independent count expectations.

    This is not produced by an auditor and does not establish scientific
    validity. Its fixture field remains in every test copy and completion seal.
    """
    output = Path(output)
    metadata = read_json(output/'run_metadata.json')
    tables = read_json(output/'table_manifest.json')
    endpoints = []
    for cell in design.cells:
        for rep in sorted({0, design.replications-1}):
            row = next(r for r in rows if r['cell_id'] == cell.identity and r['rep'] == rep)
            endpoints.append({'cell_id': cell.identity, 'replication': rep,
                              'times_sha256': row['times_sha256'],
                              'observations_sha256': row['observations_sha256'],
                              'estimators_checked': ['exact', 'pfml', 'euler']})
    report = {
        'fixture': 'literal stage outcome; NOT an independent numerical audit',
        'schema_version': 1, 'status': 'passed', 'full_completion_audit': True,
        'replay_performed': True, 'design_sha256': design.digest,
        'run_identity_sha256': metadata['provenance']['identity_sha256'],
        'code_commit': metadata['provenance']['code_commit'],
        'source_dependency_files_verified': len(metadata['provenance']['source_files']),
        'audited_artifact_sha256': historical_audit.artifact_fingerprint(output),
        'audited_run_metadata_sha256': historical_audit.semantic_metadata_digest(metadata),
        'cells_verified': len(design.cells),
        'partitions_verified': sum(len(read_json(output/'cells'/cell.identity/'cell_manifest.json')['partitions'])
                                  for cell in design.cells),
        **expected_counts(design), 'pipeline_errors': 0,
        'replayed_source_paths': len(endpoints), 'replayed_fits': 3*len(endpoints),
        'replay_records': endpoints,
        'replay_selection': ('replications 0 and R-1 in every physical floor variant, all estimators'
                             if kind == 'floor' else 'replications 0 and R-1 in every cell, all estimators'),
        'tables_verified': {name: entry['rows'] for name, entry in tables.items()
                            if name != 'failed_fits.csv.gz'},
        'tolerances': {'arithmetic_rtol': 2e-10, 'arithmetic_atol': 2e-12,
                       'replay_rtol': 2e-10, 'replay_atol': 2e-12, 'covariance_rtol': 2e-9,
                       'covariance_atol': '2e-12 * max(1, expected matrix Frobenius norm)'},
    }
    if kind == 'floor':
        binding = read_json(output/'epsilon_table_manifest.json')['raw_inputs']
        report.update(audit_contract='week7_independent_paired_floor_audit_v1',
                      scientific_overrides=['paired_floor_source'], source_policy=metadata['source_policy'],
                      raw_input_sha256=hashlib.sha256(canonical(binding).encode()).hexdigest(),
                      planned_source_realizations=12*design.replications,
                      planned_independent_coupling_groups=4*design.replications,
                      source_replay_performed=True, source_paths_verified=12*design.replications,
                      source_array_hashes_verified=24*design.replications)
        report['tables_verified']['epsilon_comparisons.csv'] = 72
        report['tolerances']['source_arrays_and_diagnostics'] = 'exact types, values and SHA-256 array bytes'
    else:
        report['tolerances'].update(
            span_roundoff='8*gamma_(K+8)*max(span,K*nominal_gap)+8 ulps; gamma_m=m*eps/(1-m*eps)',
            minimum_gap_roundoff='4 ulps at the realized-span scale')
    return report


def renderer_sibling(output):
    output = Path(output)
    return output.with_name(output.name+'-figure-renderer')


def copy_archive(source, target):
    shutil.copytree(source, target)
    if renderer_sibling(source).exists():
        shutil.copytree(renderer_sibling(source), renderer_sibling(target))
    return target


def literal_review(output):
    figures = read_json(Path(output)/'figures/figure_manifest.json')
    return {'reviewer': 'literal test review input; no actual visual review claimed',
            'reviewed_pngs': sorted(name for name in figures['figures'] if name.endswith('.png')),
            'all_panels_legible': True,
            'figure_manifest_sha256': file_digest(Path(output)/'figures/figure_manifest.json')}


@pytest.fixture(scope='module', params=('followup', 'floor'))
def literal_case(request, tmp_path_factory):
    kind = request.param
    design = runner.design_for(kind, validation_replications=2)
    rows = literal_records(design, kind)
    by_source = {}
    for row in rows:
        by_source.setdefault((row['cell_id'], row['rep']), []).append(row)
    base = tmp_path_factory.mktemp('runner_literal_integrity_'+kind)
    raw = base/'raw'
    with block_science() as patch:
        patch.setattr(historical_execution, 'source_records',
                      lambda cell, rep, source, services: deepcopy(by_source[cell.identity, rep]))
        execution.run_raw(kind, design, raw, extra_files=runner.WEEK_FILES,
                          workers=1, partition_size=1, progress=False)
    audited = copy_archive(raw, base/'audited')
    with block_science() as patch:
        patch.setattr(historical_audit if kind == 'followup' else floor_audit, 'audit_run',
                      lambda output, run_design, **kwargs: audit_outcome(output, run_design, kind, rows))
        runner.postprocess(kind, audited)
    assert read_json(audited/'run_metadata.json')['status'] == 'audited_pending_visual_review'
    sealed = copy_archive(audited, base/'sealed')
    with block_science():
        runner.seal_run(kind, sealed, literal_review(sealed))
    return SimpleNamespace(kind=kind, design=design, rows=rows, raw=raw, audited=audited, sealed=sealed)


@pytest.fixture
def case(literal_case, tmp_path):
    source = literal_case
    return SimpleNamespace(kind=source.kind, design=source.design, rows=source.rows,
                           root=tmp_path, raw_source=source.raw,
                           output=copy_archive(source.audited, tmp_path/'run'),
                           sealed_source=source.sealed)


def test_audited_and_completed_resumption_is_byte_identical_and_never_reexecutes(case, monkeypatch):
    monkeypatch.setattr(execution, 'run_raw', forbidden)
    monkeypatch.setattr(reporting, 'write_tables', forbidden)
    monkeypatch.setattr(plots, 'write_figures', forbidden)
    for source, name in ((case.output, 'audited'), (case.sealed_source, 'sealed')):
        output = source if name == 'audited' else copy_archive(source, case.root/name)
        before, renderer = snapshot(output), snapshot(renderer_sibling(output))
        result = runner.run_study(case.kind, output, validation_replications=2,
                                  resume=True, progress=False)
        assert result['planned_fit_records'] == 3*len(case.design.cells)*2
        assert snapshot(output) == before
        assert snapshot(renderer_sibling(output)) == renderer


def test_repeated_seal_is_read_only_and_changed_review_is_refused(case):
    output = copy_archive(case.sealed_source, case.root/'sealed')
    before = snapshot(output)
    review = read_json(output/'visual_review.json')
    result = runner.seal_run(case.kind, output, review)
    assert result == read_json(output/'run_metadata.json')
    changed = dict(review, reviewer='different literal reviewer')
    with pytest.raises(ValueError):
        runner.seal_run(case.kind, output, changed)
    assert snapshot(output) == before


@pytest.mark.parametrize('fault', ('duplicate_png', 'missing_png', 'stale_manifest', 'blank_reviewer', 'not_legible'))
def test_visual_review_requires_every_current_png_once(case, fault):
    review = literal_review(case.output)
    if fault == 'duplicate_png':
        review['reviewed_pngs'].append(review['reviewed_pngs'][0])
    elif fault == 'missing_png':
        review['reviewed_pngs'].pop()
    elif fault == 'stale_manifest':
        review['figure_manifest_sha256'] = '0'*64
    elif fault == 'blank_reviewer':
        review['reviewer'] = ' '
    else:
        review['all_panels_legible'] = False
    before = snapshot(case.output)
    with pytest.raises(ValueError):
        runner.seal_run(case.kind, case.output, review)
    assert snapshot(case.output) == before


@pytest.mark.parametrize('mutation', ('missing_runner', 'extra_dependency'))
def test_exact_dependency_membership_is_required_before_resume(case, monkeypatch, mutation):
    metadata = read_json(case.output/'run_metadata.json')
    source_files = metadata['provenance']['source_files']
    if mutation == 'missing_runner':
        del source_files['ou_irregular/week7/runner.py']
    else:
        source_files['README.md'] = file_digest(Path(__file__).parents[2]/'README.md')
    write_json(case.output/'run_metadata.json', metadata)
    monkeypatch.setattr(execution, 'run_raw', forbidden)
    before = snapshot(case.output)
    with pytest.raises(ValueError, match='membership'):
        runner.run_study(case.kind, case.output, validation_replications=2, resume=True)
    assert snapshot(case.output) == before


@pytest.mark.parametrize('field', ('point_valid_records', 'full_completion_audit', 'replay_records'))
def test_saved_full_audit_tampering_is_not_replaced_on_resume(case, monkeypatch, field):
    metadata = read_json(case.output/'run_metadata.json')
    metadata['status'] = 'raw_complete_pending_audit'
    metadata.pop('audit_sha256', None)
    write_json(case.output/'run_metadata.json', metadata)
    audit = read_json(case.output/'audit.json')
    if field == 'point_valid_records':
        audit[field] -= 1
    elif field == 'full_completion_audit':
        audit[field] = False
    else:
        audit[field][0]['observations_sha256'] = '0'*64
    write_json(case.output/'audit.json', audit)
    monkeypatch.setattr(execution, 'run_raw', forbidden)
    monkeypatch.setattr(reporting, 'write_tables', forbidden)
    before = snapshot(case.output)
    with pytest.raises(ValueError):
        runner.run_study(case.kind, case.output, validation_replications=2, resume=True)
    assert snapshot(case.output) == before


def test_partial_figure_resume_reuses_audit_and_preserves_every_existing_byte(case, monkeypatch):
    output = copy_archive(case.raw_source, case.root/'partial')
    module = historical_audit if case.kind == 'followup' else floor_audit
    calls = []
    def audit_boundary(directory, design, **kwargs):
        calls.append('audit')
        return audit_outcome(directory, design, case.kind, case.rows)
    monkeypatch.setattr(module, 'audit_run', audit_boundary)
    original = plots.figure_bytes
    stems = plots.FOLLOWUP_STEMS if case.kind == 'followup' else plots.FLOOR_STEMS
    def interrupted(figure, stem):
        if stem == stems[1]:
            raise RuntimeError('literal interrupted second figure')
        return original(figure, stem)
    monkeypatch.setattr(plots, 'figure_bytes', interrupted)
    with pytest.raises(RuntimeError, match='interrupted second'):
        runner.postprocess(case.kind, output)
    assert calls == ['audit']
    assert read_json(output/'run_metadata.json')['status'] == 'raw_complete_pending_audit'
    assert not (output/'figures/figure_manifest.json').exists()
    existing = {str(path.relative_to(output)): path.read_bytes()
                for path in output.rglob('*') if path.is_file() and path.name != 'run_metadata.json'}
    renderer_before = snapshot(renderer_sibling(output))
    monkeypatch.setattr(module, 'audit_run', forbidden)
    monkeypatch.setattr(execution, 'run_raw', forbidden)
    monkeypatch.setattr(reporting, 'write_tables', forbidden)
    monkeypatch.setattr(plots, 'figure_bytes', original)
    runner.run_study(case.kind, output, validation_replications=2, resume=True)
    assert read_json(output/'run_metadata.json')['status'] == 'audited_pending_visual_review'
    assert all((output/name).read_bytes() == value for name, value in existing.items())
    assert snapshot(renderer_sibling(output)) == renderer_before


def test_root_drift_during_reporting_is_preserved_and_blocks_later_mutation(case, monkeypatch):
    output = copy_archive(case.raw_source, case.root/'drift')
    before = (output/'run_metadata.json').read_bytes()
    original = reporting.write_tables
    observed = []
    def changed_root(*args, **kwargs):
        result = original(*args, **kwargs)
        root = read_json(output/'run_metadata.json')
        root['status'] = 'running_incomplete'
        write_json(output/'run_metadata.json', root)
        observed.append((output/'run_metadata.json').read_bytes())
        return result
    monkeypatch.setattr(reporting, 'write_tables', changed_root)
    with pytest.raises(path_diagnostics.ManifestDriftError):
        runner.postprocess(case.kind, output)
    assert (output/'run_metadata.json').read_bytes() == observed[0]
    incidents = list(output.parent.glob(output.name+'-manifest-incident-*'))
    assert len(incidents) == 1
    assert (incidents[0]/'expected_run_metadata.json').read_bytes() == before
    assert (incidents[0]/'observed_run_metadata.json').read_bytes() == observed[0]
    (output/'run_metadata.json').write_bytes(before)  # injected restoration, not a repair
    frozen = snapshot(output)
    with pytest.raises(ValueError, match='incident'):
        runner.postprocess(case.kind, output)
    assert snapshot(output) == frozen


def test_nonroot_artifact_change_during_visual_review_publication_cannot_be_sealed(case, monkeypatch):
    original = runner.atomic_json
    target = case.output/'parameter_summary.csv'
    mutated = []
    def publication(path, value, **kwargs):
        result = original(path, value, **kwargs)
        if Path(path).name == 'visual_review.json':
            target.write_bytes(target.read_bytes()+b'\n')
            mutated.append(target.read_bytes())
        return result
    monkeypatch.setattr(runner, 'atomic_json', publication)
    with pytest.raises(ValueError):
        runner.seal_run(case.kind, case.output, literal_review(case.output))
    assert mutated and target.read_bytes() == mutated[0]
    assert read_json(case.output/'run_metadata.json')['status'] == 'audited_pending_visual_review'


def test_root_drift_during_visual_review_publication_cannot_be_overwritten(case, monkeypatch):
    original = runner.atomic_json
    observed = []
    def publication(path, value, **kwargs):
        result = original(path, value, **kwargs)
        if Path(path).name == 'visual_review.json':
            root = read_json(case.output/'run_metadata.json')
            root['status'] = 'incomplete'
            write_json(case.output/'run_metadata.json', root)
            observed.append((case.output/'run_metadata.json').read_bytes())
        return result
    monkeypatch.setattr(runner, 'atomic_json', publication)
    with pytest.raises(path_diagnostics.ManifestDriftError):
        runner.seal_run(case.kind, case.output, literal_review(case.output))
    assert (case.output/'run_metadata.json').read_bytes() == observed[0]
    assert list(case.output.parent.glob(case.output.name+'-manifest-incident-*'))


def test_full_run_requires_frozen_dependencies_before_any_output(tmp_path, monkeypatch):
    seen = []
    def unfrozen(design, extra_files=(), **kwargs):
        seen.append((design.scope, kwargs.get('require_frozen'), set(extra_files)))
        raise ValueError('literal uncommitted runner dependency')
    monkeypatch.setattr(execution, 'capture', unfrozen)
    output = tmp_path/'must_not_exist'
    with pytest.raises(ValueError, match='uncommitted'):
        runner.run_study('floor', output, progress=False)
    assert len(seen) == 1 and seen[0][:2] == ('exploratory', True)
    assert 'ou_irregular/week7/runner.py' in seen[0][2]
    assert not output.exists()
