"""Independent saved-data plotting integrity and lifecycle verification.

The floor input is the byte-preserved R=2 validation archive and its already
computed full audit. No test here generates a source or fits an estimator.
Selected diagnostic fixtures, when used, are explicitly literal unit inputs.
"""
from copy import deepcopy
from dataclasses import replace
import csv
import hashlib
import io
import json
from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from ou_irregular.experiments import execution as raw_execution
from ou_irregular.experiments import audit as registered_audit
from ou_irregular.experiments.provenance import ROOT, CORE_FILES, capture, canonical
from ou_irregular.experiments.storage import run_lock
from ou_irregular.inference.adapters import Week3Estimator
from ou_irregular.week4.simulation import OuSampleSource
from ou_irregular.week7 import execution, floor_audit, path_diagnostics, plots
from ou_irregular.week7.floor_source import PairedFloorSource
from ou_irregular.week7.selection import floor_design


FLOOR_NAMES = {stem+'.'+extension for stem in (
    'W7_floor_paired_bias', 'W7_floor_paired_coverage',
    'W7_floor_rmse_ratio', 'W7_floor_spacing') for extension in ('png', 'pdf')}


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_bytes())


def write_json(path, value):
    Path(path).write_bytes(encoded(value))


def snapshot(output):
    output = Path(output)
    if not output.exists():
        return None
    result = {}
    for path in sorted(output.rglob('*')):
        relative = str(path.relative_to(output))
        if path.is_symlink():
            result[relative] = ('symlink', str(path.readlink()))
        elif path.is_file():
            result[relative] = ('file', sha(path))
        elif path.is_dir():
            result[relative] = ('directory',)
    return result


def renderer_state(output):
    return output.with_name(output.name+'-figure-renderer')


def forbidden(*args, **kwargs):
    raise AssertionError('saved-data plotting must not generate, fit, replay, or execute source paths')


@pytest.fixture(autouse=True)
def forbid_scientific_execution(monkeypatch):
    monkeypatch.setattr(OuSampleSource, 'sample', forbidden)
    monkeypatch.setattr(PairedFloorSource, 'sample', forbidden)
    monkeypatch.setattr(Week3Estimator, 'fit', forbidden)
    monkeypatch.setattr(execution, 'run_raw', forbidden)
    monkeypatch.setattr(raw_execution, 'run_raw', forbidden)
    monkeypatch.setattr(raw_execution, 'source_records', forbidden)
    monkeypatch.setattr(floor_audit, 'audit_run', forbidden)
    monkeypatch.setattr(floor_audit, 'source_records', forbidden)
    monkeypatch.setattr(floor_audit, 'independent_source', forbidden)
    monkeypatch.setattr(registered_audit, 'audit_run', forbidden)
    monkeypatch.setattr(path_diagnostics, 'run_selected', forbidden)
    monkeypatch.setattr(path_diagnostics, '_one_selected', forbidden)
    monkeypatch.setattr(path_diagnostics, 'analyze_sample', forbidden)


@pytest.fixture(scope='module')
def saved_floor(tmp_path_factory):
    from .test_plots_reference import saved_floor_archive
    return saved_floor_archive(tmp_path_factory)


@pytest.fixture(scope='module')
def rendered_floor(tmp_path_factory):
    from .test_plots_reference import rendered_floor_archive
    return rendered_floor_archive(tmp_path_factory)


@pytest.fixture
def copied_saved(saved_floor, tmp_path):
    source, design = saved_floor
    output = tmp_path/'run'
    shutil.copytree(source, output)
    return output, design


@pytest.fixture
def copied_rendered(rendered_floor, tmp_path):
    source, design = rendered_floor
    output = tmp_path/'run'
    shutil.copytree(source, output)
    state = renderer_state(source)
    if state.exists():
        shutil.copytree(state, renderer_state(output))
    return output, design


@pytest.fixture
def epsilon_table():
    archive = ROOT/'artifacts/week7_recovery_evidence/floor_audit_current_capture/archive'
    return pd.read_csv(archive/'epsilon_comparisons.csv', float_precision='round_trip')


@pytest.mark.parametrize('fault', ('fractional_common', 'boolean_common', 'negative_common',
                                  'common_exceeds_own', 'coverage_total', 'covered_exceeds_valid'))
def test_complete_epsilon_keys_do_not_allow_false_count_domains(epsilon_table, fault):
    table = epsilon_table.copy(deep=True)
    group = ((table.cv == 1.5) & (table.floor_fraction == 1e-8)
             & (table.estimator == 'exact'))
    first = table.index[group][0]
    if fault in {'fractional_common', 'boolean_common', 'negative_common'}:
        table['n_common_point_valid'] = table['n_common_point_valid'].astype(object)
        value = {'fractional_common': 1.5, 'boolean_common': True, 'negative_common': -1}[fault]
        table.loc[group, 'n_common_point_valid'] = value
    elif fault == 'common_exceeds_own':
        table.loc[group, 'n_point_valid_a'] = 0
    elif fault == 'coverage_total':
        table.loc[first, 'n_neither_cover'] += 1
    else:
        table.loc[first, 'n_covered_a'] = table.loc[first, 'n_interval_valid_a'] + 1
    with pytest.raises(ValueError, match='count|contingency|valid|coverage'):
        plots._epsilon_inputs(table, floor_design().validation(2))


def test_real_png_pdf_rerender_is_deterministic_and_ignores_ambient_rc(epsilon_table):
    archive = ROOT/'artifacts/week7_recovery_evidence/floor_audit_current_capture/archive'
    tables = {name: pd.read_csv(archive/name, float_precision='round_trip')
              for name in ('parameter_summary.csv', 'spacing_summary.csv')}
    design = floor_design().validation(2)
    stem = 'W7_floor_paired_bias'
    results = []
    for settings in (
            {'axes.labelsize': 7, 'savefig.transparent': False, 'text.usetex': False},
            {'axes.labelsize': 31, 'savefig.transparent': True, 'text.usetex': True,
             'font.family': ['serif'], 'mathtext.fontset': 'cm', 'savefig.dpi': 41}):
        with plt.rc_context(settings):
            figures = plots.build_floor_figures(tables, epsilon_table, design)
            try:
                results.append(plots.figure_bytes(figures[stem], stem))
            finally:
                for figure in figures.values():
                    plt.close(figure)
    assert set(results[0]) == {stem+'.png', stem+'.pdf'}
    assert results[0] == results[1]
    assert results[0][stem+'.png'].startswith(b'\x89PNG\r\n\x1a\n')
    assert results[0][stem+'.pdf'].startswith(b'%PDF-')
    assert b'/CreationDate' not in results[0][stem+'.pdf']
    assert b'/ModDate' not in results[0][stem+'.pdf']


def assert_failed_without_writes(output, design, action=None, *, match=None):
    before, state_before = snapshot(output), snapshot(renderer_state(output))
    with pytest.raises((ValueError, FileNotFoundError, FileExistsError), match=match):
        (action or (lambda: plots.write_figures('floor', output, design)))()
    assert snapshot(output) == before
    assert snapshot(renderer_state(output)) == state_before


def refresh_identity(provenance):
    identity = {key: provenance[key] for key in (
        'design', 'design_sha256', 'preregistration_commit', 'source_files', 'environment')}
    provenance['identity_sha256'] = hashlib.sha256(canonical(identity).encode()).hexdigest()


def refresh_audit_artifact_map(output):
    """Fault injection only: do not let an outer hash hide a structural error.

    This deliberately counterfeits a copied saved report. It never claims a
    new numerical audit or changes the shared genuine reference archive.
    """
    report = read_json(output/'audit.json')
    report['audited_artifact_sha256'] = registered_audit.artifact_fingerprint(output)
    write_json(output/'audit.json', report)
    metadata = read_json(output/'run_metadata.json')
    if 'audit_sha256' in metadata:
        metadata['audit_sha256'] = sha(output/'audit.json')
        write_json(output/'run_metadata.json', metadata)


def alter_csv(output, name, mutation):
    rows = list(csv.reader(io.StringIO((output/name).read_text(), newline='')))
    mutation(rows)
    stream = io.StringIO(newline='')
    csv.writer(stream, lineterminator='\n').writerows(rows)
    (output/name).write_text(stream.getvalue())
    manifest = read_json(output/'table_manifest.json')
    manifest[name] = {'rows': len(rows)-1, 'sha256': sha(output/name)}
    write_json(output/'table_manifest.json', manifest)
    refresh_audit_artifact_map(output)


@pytest.mark.parametrize('fault', (
    'duplicate_endpoint', 'missing_endpoint_field', 'boolean_replication',
    'wrong_estimator_scope', 'wrong_source_hash', 'malformed_source_hash',
    'replay_selection', 'replay_tolerance', 'arithmetic_tolerance',
    'not_full', 'wrong_code_commit', 'missing_floor_override', 'wrong_source_policy',
    'false_source_count', 'false_fit_count', 'false_table_count'))
def test_saved_full_audit_claim_requires_exact_bound_scope(copied_saved, fault):
    output, design = copied_saved
    report = read_json(output/'audit.json')
    record = report['replay_records'][0]
    if fault == 'duplicate_endpoint':
        report['replay_records'][1] = deepcopy(record)
    elif fault == 'missing_endpoint_field':
        del record['times_sha256']
    elif fault == 'boolean_replication':
        record['replication'] = False
    elif fault == 'wrong_estimator_scope':
        record['estimators_checked'] = ['exact', 'pfml', 'pfml']
    elif fault == 'wrong_source_hash':
        record['observations_sha256'] = '0'*64
    elif fault == 'malformed_source_hash':
        record['times_sha256'] = 'invalid hash'
    elif fault == 'replay_selection':
        report['replay_selection'] = 'first source only'
    elif fault == 'replay_tolerance':
        report['tolerances']['replay_rtol'] = .1
    elif fault == 'arithmetic_tolerance':
        report['tolerances']['arithmetic_atol'] = .1
    elif fault == 'not_full':
        report['full_completion_audit'] = False
    elif fault == 'wrong_code_commit':
        report['code_commit'] = '0'*40
    elif fault == 'missing_floor_override':
        report['scientific_overrides'] = []
    elif fault == 'wrong_source_policy':
        report['source_policy']['root_seed'] += 1
    elif fault == 'false_source_count':
        report['source_array_hashes_verified'] -= 1
    elif fault == 'false_fit_count':
        report['replayed_fits'] -= 1
    else:
        report['tables_verified']['parameter_summary.csv'] -= 1
    write_json(output/'audit.json', report)
    assert_failed_without_writes(output, design, match='audit|endpoint|tolerance')


@pytest.mark.parametrize('name', tuple(floor_audit.TABLE_KEYS))
def test_every_standard_table_requires_complete_unique_keys_before_state_creation(copied_saved, name):
    output, design = copied_saved
    keys = floor_audit.TABLE_KEYS[name]
    def duplicate_key(rows):
        indices = [rows[0].index(key) for key in keys]
        for index in indices:
            rows[2][index] = rows[1][index]
    alter_csv(output, name, duplicate_key)
    assert_failed_without_writes(output, design, match='key|duplicate|group|incomplete')


@pytest.mark.parametrize('fault', ('spacing_header', 'spacing_denominator', 'summary_coordinate', 'summary_header',
                                  'outcome_total', 'numerical_coordinate'))
def test_rehashed_tables_still_require_complete_preflight(copied_saved, fault):
    output, design = copied_saved
    name = ('spacing_summary.csv' if fault.startswith('spacing') else
            'parameter_summary.csv' if fault.startswith('summary') else
            'outcome_counts.csv' if fault.startswith('outcome') else 'numerical_diagnostics.csv')
    def change(rows):
        if fault == 'spacing_header':
            rows[0][rows[0].index('min_gap_mcse')] = 'unrecognized_mcse'
        elif fault == 'summary_header':
            rows[0][rows[0].index('operational_wilson_lower')] = 'unrecognized_lower'
        elif fault == 'spacing_denominator':
            rows[1][rows[0].index('min_gap_n_available')] = '1'
        elif fault == 'outcome_total':
            rows[1][rows[0].index('count')] = '1'
        else:
            rows[1][rows[0].index('n')] = '250'
    alter_csv(output, name, change)
    assert_failed_without_writes(output, design, match='header|source|R |coordinate|count|fraction|outcome')


def test_real_manifest_binds_epsilon_separately_and_does_not_upgrade_raw_provenance(rendered_floor):
    output, design = rendered_floor
    metadata = read_json(output/'run_metadata.json')
    manifest = plots.verify_figures('floor', output, design)
    state = read_json(renderer_state(output)/'renderer_state.json')
    assert metadata['status'] == 'raw_complete_pending_audit'
    assert 'audit_sha256' not in metadata
    assert set(manifest['figures']) == FLOOR_NAMES
    assert manifest['renderer'] == state
    assert manifest['renderer_state_sha256'] == hashlib.sha256(canonical(state).encode()).hexdigest()
    assert set(manifest['input_table_sha256']) == set(plots.reporting.STANDARD_FILES)
    assert manifest['epsilon_inputs'] == {
        name: sha(output/name) for name in ('epsilon_comparisons.csv', 'epsilon_table_manifest.json')}
    assert manifest['input_binding']['saved_audit_full_completion'] is True
    assert manifest['new_sources_or_fits'] is False
    assert manifest['headline_results_replaced'] is False
    assert manifest['visual_review_performed'] is False
    assert 'ou_irregular/week7/plots.py' not in metadata['provenance']['source_files']
    assert set(state['provenance']['source_files']) == set(CORE_FILES) | set(plots.FILES)
    assert state['descriptor']['renderer_environment']['pdf_creation_date'] is None
    assert state['descriptor']['renderer_environment']['pdf_modification_date'] is None
    assert state['descriptor']['require_frozen'] is False


def test_complete_readback_is_pure_self_contained_and_does_not_recapture(copied_rendered, monkeypatch):
    output, design = copied_rendered
    shutil.rmtree(renderer_state(output))
    # A complete archive is deliberately independent of optional sidecar state.
    renderer_state(output).mkdir()
    (renderer_state(output)/'unrelated.txt').write_text('optional external material')
    before, state_before = snapshot(output), snapshot(renderer_state(output))
    mtimes = {str(p.relative_to(output)): p.stat().st_mtime_ns for p in output.rglob('*') if p.is_file()}
    for name in ('capture', 'run_lock', 'atomic_bytes', 'atomic_json', 'figure_bytes'):
        monkeypatch.setattr(plots, name, forbidden)
    first = plots.write_figures('floor', output, design)
    second = plots.verify_figures('floor', output, design)
    assert first == second
    assert snapshot(output) == before
    assert snapshot(renderer_state(output)) == state_before
    assert {str(p.relative_to(output)): p.stat().st_mtime_ns for p in output.rglob('*') if p.is_file()} == mtimes


@pytest.mark.parametrize('fault', ('missing_png', 'changed_pdf', 'extra_name', 'nested_directory',
                                  'image_symlink', 'manifest_missing_file', 'wrong_epsilon_hash'))
def test_complete_declared_products_are_never_recreated_or_overwritten(copied_rendered, fault, monkeypatch):
    output, design = copied_rendered
    images = output/'figures'
    png, pdf = images/'W7_floor_paired_bias.png', images/'W7_floor_paired_bias.pdf'
    if fault == 'missing_png':
        png.unlink()
    elif fault == 'changed_pdf':
        pdf.write_bytes(pdf.read_bytes()+b'changed')
    elif fault == 'extra_name':
        (images/'extra.png').write_bytes(b'extra')
    elif fault == 'nested_directory':
        (images/'nested').mkdir()
    elif fault == 'image_symlink':
        png.unlink()
        png.symlink_to(pdf)
    else:
        manifest = read_json(images/'figure_manifest.json')
        if fault == 'manifest_missing_file':
            del manifest['figures'][png.name]
        else:
            manifest['epsilon_inputs']['epsilon_comparisons.csv'] = '0'*64
        write_json(images/'figure_manifest.json', manifest)
    monkeypatch.setattr(plots, 'figure_bytes', forbidden)
    assert_failed_without_writes(output, design)


@pytest.mark.parametrize('fault', ('missing_dependency', 'extra_dependency', 'changed_dependency',
                                  'changed_environment', 'renderer_environment', 'descriptor_digest'))
def test_embedded_renderer_provenance_is_exact_and_verified(copied_rendered, fault):
    output, design = copied_rendered
    path = output/'figures/figure_manifest.json'
    manifest = read_json(path)
    state = manifest['renderer']
    provenance = state['provenance']
    if fault == 'missing_dependency':
        del provenance['source_files']['ou_irregular/week7/plots.py']
    elif fault == 'extra_dependency':
        provenance['source_files']['README.md'] = sha(ROOT/'README.md')
    elif fault == 'changed_dependency':
        provenance['source_files']['ou_irregular/week7/plots.py'] = '0'*64
    elif fault == 'changed_environment':
        provenance['environment']['machine'] = 'invented-platform'
    elif fault == 'renderer_environment':
        state['descriptor']['renderer_environment']['pillow'] = 'invented-version'
    else:
        state['descriptor_sha256'] = '0'*64
    refresh_identity(provenance)
    manifest['renderer_state_sha256'] = hashlib.sha256(canonical(state).encode()).hexdigest()
    write_json(path, manifest)
    assert_failed_without_writes(output, design, match='dependenc|environment|state|identity')


def cached_pairs(rendered_floor):
    output, _ = rendered_floor
    return {stem: {stem+'.'+ext: (output/'figures'/(stem+'.'+ext)).read_bytes()
                   for ext in ('png', 'pdf')} for stem in plots.FLOOR_STEMS}


def make_partial(output, keep=('W7_floor_paired_bias.png',)):
    for path in (output/'figures').iterdir():
        if path.name not in keep:
            path.unlink()


def test_partial_byte_equal_files_are_adopted_with_original_state(copied_rendered, rendered_floor, monkeypatch):
    output, design = copied_rendered
    pairs = cached_pairs(rendered_floor)
    make_partial(output)
    retained = output/'figures/W7_floor_paired_bias.png'
    before = retained.read_bytes(), retained.stat().st_mtime_ns
    root_bytes = (output/'run_metadata.json').read_bytes()
    state_bytes = (renderer_state(output)/'renderer_state.json').read_bytes()
    visited = []
    def saved_serialization(figure, stem):
        visited.append(stem)
        return pairs[stem]
    monkeypatch.setattr(plots, 'figure_bytes', saved_serialization)
    monkeypatch.setattr(plots, 'capture', forbidden)
    manifest = plots.write_figures('floor', output, design)
    assert visited == list(plots.FLOOR_STEMS)
    assert (retained.read_bytes(), retained.stat().st_mtime_ns) == before
    assert (output/'run_metadata.json').read_bytes() == root_bytes
    assert (renderer_state(output)/'renderer_state.json').read_bytes() == state_bytes
    assert set(manifest['figures']) == FLOOR_NAMES


@pytest.mark.parametrize('fault', ('missing_state', 'changed_partial', 'extra_state_file',
                                  'state_symlink', 'state_descriptor', 'empty_state_directory'))
def test_partial_outputs_fail_closed_without_original_matching_state(copied_rendered, rendered_floor, monkeypatch, fault):
    output, design = copied_rendered
    pairs = cached_pairs(rendered_floor)
    make_partial(output)
    state_dir = renderer_state(output)
    state_path = state_dir/'renderer_state.json'
    if fault == 'missing_state':
        shutil.rmtree(state_dir)
    elif fault == 'changed_partial':
        (output/'figures/W7_floor_paired_bias.png').write_bytes(b'foreign partial bytes')
    elif fault == 'extra_state_file':
        (state_dir/'extra.txt').write_text('unexpected')
    elif fault == 'state_symlink':
        stored = state_path.read_bytes()
        state_path.unlink()
        target = output.parent/'state-target.json'
        target.write_bytes(stored)
        state_path.symlink_to(target)
    elif fault == 'state_descriptor':
        state = read_json(state_path)
        state['descriptor']['input_binding']['audit_sha256'] = '0'*64
        write_json(state_path, state)
    else:
        state_path.unlink()
    monkeypatch.setattr(plots, 'figure_bytes', lambda figure, stem: pairs[stem])
    assert_failed_without_writes(output, design)


def test_incomplete_render_owns_the_run_lock_before_images_or_state(copied_saved, monkeypatch):
    output, design = copied_saved
    monkeypatch.setattr(plots, 'figure_bytes', forbidden)
    with run_lock(output):
        before, state_before = snapshot(output), snapshot(renderer_state(output))
        with pytest.raises(RuntimeError, match='lock'):
            plots.write_figures('floor', output, design)
        assert snapshot(output) == before
        assert snapshot(renderer_state(output)) == state_before


def test_renderer_state_cannot_be_nested_inside_an_existing_archive(copied_saved):
    output, design = copied_saved
    write_json(output.parent/'run_metadata.json', {'status': 'unrelated_existing_archive'})
    assert_failed_without_writes(output, design, match='archiv|nested|inside')


@pytest.mark.parametrize('raise_after_change', (False, True))
def test_external_root_drift_is_preserved_and_blocks_automatic_resume(copied_rendered, rendered_floor, monkeypatch, raise_after_change):
    output, design = copied_rendered
    make_partial(output, keep=())
    original = (output/'run_metadata.json').read_bytes()
    foreign = read_json(output/'run_metadata.json')
    foreign['external_drift_literal_fixture'] = True
    foreign_bytes = encoded(foreign)
    pairs = cached_pairs(rendered_floor)
    def drift(figure, stem):
        (output/'run_metadata.json').write_bytes(foreign_bytes)
        if raise_after_change:
            raise KeyboardInterrupt('literal interruption after external write')
        return pairs[stem]
    monkeypatch.setattr(plots, 'figure_bytes', drift)
    with pytest.raises(path_diagnostics.ManifestDriftError, match='manifest drift'):
        plots.write_figures('floor', output, design)
    assert (output/'run_metadata.json').read_bytes() == foreign_bytes
    assert not (output/'figures/figure_manifest.json').exists()
    incidents = [path for path in output.parent.iterdir()
                 if path.name.startswith(output.name+'-manifest-incident-')]
    assert len(incidents) == 1
    incident = incidents[0]
    assert (incident/'expected_run_metadata.json').read_bytes() == original
    assert (incident/'observed_run_metadata.json').read_bytes() == foreign_bytes
    record = read_json(incident/'incident.json')
    assert record['observed_manifest_overwritten_by_handler'] is False
    # Even a human restoring the expected bytes does not silently erase the
    # unresolved incident or permit an automatic continuation.
    (output/'run_metadata.json').write_bytes(original)
    before_incident = snapshot(incident)
    assert_failed_without_writes(output, design, match='incident')
    assert snapshot(incident) == before_incident


@pytest.mark.parametrize('target', ('epsilon_comparisons.csv', 'audit.json', 'renderer_state'))
def test_nonroot_input_or_renderer_state_drift_blocks_image_commit(copied_rendered, rendered_floor, monkeypatch, target):
    output, design = copied_rendered
    make_partial(output, keep=())
    root = (output/'run_metadata.json').read_bytes()
    pairs = cached_pairs(rendered_floor)
    path = renderer_state(output)/'renderer_state.json' if target == 'renderer_state' else output/target
    foreign = path.read_bytes()+b' '
    def drift(figure, stem):
        path.write_bytes(foreign)
        return pairs[stem]
    monkeypatch.setattr(plots, 'figure_bytes', drift)
    with pytest.raises(ValueError, match='changed|state|input'):
        plots.write_figures('floor', output, design)
    assert path.read_bytes() == foreign
    assert (output/'run_metadata.json').read_bytes() == root
    assert not list((output/'figures').iterdir())


def test_actual_recorded_git_bytes_are_checked_without_requiring_current_head(copied_saved):
    _, design = copied_saved
    provenance = capture(design, require_frozen=True)
    # All core dependencies predate the renderer and are unchanged at this
    # genuine older commit; unrelated later commits are not a new code identity.
    provenance['code_commit'] = '517f5664d14ad06650b0a05e6a77c5de2699032f'
    plots._frozen_provenance(provenance, design, required=set(CORE_FILES),
                             exact=set(CORE_FILES), frozen=True)
    provenance['code_commit'] = '0'*40
    with pytest.raises(ValueError, match='recorded frozen commit'):
        plots._frozen_provenance(provenance, design, required=set(CORE_FILES),
                                 exact=set(CORE_FILES), frozen=True)


def test_full_renderer_freeze_cannot_be_counterfeited_by_clean_dirty_list(copied_saved):
    output, design = copied_saved
    prepared = replace(plots._load_run_inputs('floor', output, design), frozen=True)
    # This is an isolated freeze-gate fixture, not a relabelled full study.
    descriptor = plots._descriptor(prepared)
    renderer_design = plots._RendererDesign(canonical(descriptor))
    provenance = capture(renderer_design, plots.FILES, require_frozen=False)
    provenance['uncommitted_dependency_files'] = []
    provenance['code_commit'] = '517f5664d14ad06650b0a05e6a77c5de2699032f'
    state = {'schema_version': 1, 'rendering_contract': plots.CONTRACT,
             'descriptor': descriptor, 'descriptor_sha256': renderer_design.digest,
             'provenance': provenance}
    before = snapshot(output)
    with pytest.raises(ValueError, match='recorded frozen commit'):
        plots._validate_renderer_state(state, descriptor, ())
    assert snapshot(output) == before
    assert not renderer_state(output).exists()


def seal_literal_completed(output):
    """A lifecycle fixture only; this is never a visual-review deliverable."""
    manifest = read_json(output/'figures/figure_manifest.json')
    review = {'all_panels_legible': True,
              'reviewed_pngs': sorted(name for name in manifest['figures'] if name.endswith('.png')),
              'figure_manifest_sha256': sha(output/'figures/figure_manifest.json'),
              'reviewer': 'literal lifecycle fixture; not an actual visual certification'}
    write_json(output/'visual_review.json', review)
    metadata = read_json(output/'run_metadata.json')
    metadata.update(status='completed', audit_sha256=sha(output/'audit.json'),
                    visual_review_sha256=sha(output/'visual_review.json'))
    metadata['artifact_sha256'] = {str(path.relative_to(output)): sha(path)
                                  for path in output.rglob('*') if path.is_file()
                                  and path.name not in {'run_metadata.json', '.execution.lock'}}
    write_json(output/'run_metadata.json', metadata)


def test_legitimate_later_lifecycle_does_not_invalidate_semantic_figure_binding(copied_rendered, monkeypatch):
    output, design = copied_rendered
    original = read_json(output/'figures/figure_manifest.json')
    seal_literal_completed(output)
    before = snapshot(output), snapshot(renderer_state(output))
    for name in ('figure_bytes', 'capture', 'run_lock', 'atomic_bytes', 'atomic_json'):
        monkeypatch.setattr(plots, name, forbidden)
    assert plots.write_figures('floor', output, design) == original
    assert plots.verify_figures('floor', output, design) == original
    assert (snapshot(output), snapshot(renderer_state(output))) == before


@pytest.mark.parametrize('fault', ('missing_image', 'missing_manifest', 'missing_all_figures'))
def test_completed_archive_never_recreates_missing_figures(copied_rendered, monkeypatch, fault):
    output, design = copied_rendered
    seal_literal_completed(output)
    if fault == 'missing_image':
        (output/'figures/W7_floor_spacing.png').unlink()
    elif fault == 'missing_manifest':
        (output/'figures/figure_manifest.json').unlink()
    else:
        shutil.rmtree(output/'figures')
    monkeypatch.setattr(plots, 'figure_bytes', forbidden)
    assert_failed_without_writes(output, design)


@pytest.mark.parametrize('fault', ('audit_seal', 'root_source_policy', 'wrong_design', 'raw_dependency_omission'))
def test_saved_run_preflight_cannot_relabel_original_provenance(copied_saved, fault):
    output, design = copied_saved
    metadata = read_json(output/'run_metadata.json')
    if fault == 'audit_seal':
        metadata.update(status='audited_pending_visual_review', audit_sha256='0'*64)
    elif fault == 'root_source_policy':
        metadata['source_policy']['stationary_initialization'] = False
    elif fault == 'raw_dependency_omission':
        del metadata['provenance']['source_files']['ou_irregular/week7/floor_source.py']
        refresh_identity(metadata['provenance'])
    else:
        design = replace(design, replications=3)
    write_json(output/'run_metadata.json', metadata)
    assert_failed_without_writes(output, design)


@pytest.fixture(scope='module')
def literal_selected_archive(tmp_path_factory):
    """All 53 fixed keys, but explicitly literal unavailable diagnostics.

    Neither original E1/E2 records nor any original replay are fabricated as
    successful evidence. Every payload retains a literal replay failure and
    the closed lifecycle is diagnostic_review_required. This fixture covers
    archive placement/binding only; reference tests own evaluated artist data.
    """
    output = tmp_path_factory.mktemp('plots_literal_selected')/'closed_selected'
    (output/'sources').mkdir(parents=True)
    selection = read_json(ROOT/'artifacts/week7_evidence/selection.json')
    selection['literal_validation_fixture'] = True
    selection_bytes = encoded(selection)
    descriptor = path_diagnostics._diagnostic_design(selection, selection_bytes)
    provenance = capture(descriptor, path_diagnostics.FILES, require_frozen=False)
    metadata = {'schema_version': 1, 'project': 'ou-week7-path-diagnostics',
                'status': 'diagnostic_review_required', 'require_frozen': False,
                'provenance': provenance, 'source_files': {}}
    (output/'selection.json').write_bytes(selection_bytes)
    write_json(output/'diagnostic_design.json', descriptor.to_dict())
    originals = {}
    for item in selection['sources']:
        original = {method: {'cell_id': item['cell_id'], 'rep': item['rep'], 'estimator': method,
                             'point_valid': False, 'interval_valid': False,
                             'point_status': 'literal_unavailable', 'interval_status': 'literal_unavailable',
                             'candidate_theta': 123.0, 'theta': None, 'mu': None, 'sigma': None,
                             'literal_validation_fixture': True}
                    for method in ('exact', 'pfml', 'euler')}
        originals[path_diagnostics._source_key(item)] = original
        payload = {'schema_version': 1, 'run_identity': provenance['identity_sha256'],
                   'selection': item, 'status': 'original_replay_failed',
                   'headline_results_replaced': False, 'literal_validation_fixture': True,
                   'original_records': original,
                   'original_replay': {'status': 'failed',
                                       'error': 'Literal unavailable replay; no original E2 source was generated'}}
        filename = path_diagnostics._filename(item)
        write_json(output/'sources'/filename, payload)
        metadata['source_files'][filename] = sha(output/'sources'/filename)
    write_json(output/'run_metadata.json', metadata)
    report = path_diagnostics._audit_selected_files(output, selection, originals, metadata)
    assert report['status'] == 'review_required'
    assert report['planned_sources'] == report['retained_sources'] == 53
    assert report['verified_fit_replays'] == report['evaluated_sources'] == 0
    assert len(report['failures']) == 53
    write_json(output/'audit.json', report)
    metadata['audit_sha256'] = sha(output/'audit.json')
    write_json(output/'run_metadata.json', metadata)
    return output


@pytest.fixture
def copied_selected(literal_selected_archive, tmp_path):
    output = tmp_path/'closed_selected'
    shutil.copytree(literal_selected_archive, output)
    return output, tmp_path/'external_figures'


def test_selected_preflight_binds_all_53_payloads_and_both_fixed_cases(copied_selected):
    source, output = copied_selected
    before = snapshot(source)
    prepared = plots._load_selected_inputs(source, output)
    assert prepared.metadata['status'] == 'diagnostic_review_required'
    assert prepared.frozen is False
    assert len(prepared.binding['source_payload_sha256']) == 53
    assert len(prepared.binding['saved_audit_failures']) == 53
    assert set(prepared.payloads) == {'W7_selected_E2_cv0p5_rep000465', 'W7_selected_E2_cv0_rep001424'}
    assert prepared.binding['saved_audit_status'] == 'review_required'
    assert prepared.binding['source_payload_sha256'] == {
        path.name: sha(path) for path in (source/'sources').iterdir()}
    assert snapshot(source) == before
    assert not output.exists()
    assert not renderer_state(output).exists()


@pytest.mark.parametrize('relation', ('inside', 'same', 'ancestor', 'symlink'))
def test_selected_output_must_be_disjoint_from_closed_input(copied_selected, relation):
    source, output = copied_selected
    if relation == 'inside':
        output = source/'figures'
    elif relation == 'same':
        output = source
    elif relation == 'ancestor':
        output = source.parent
    else:
        output.symlink_to(source, target_is_directory=True)
    before = snapshot(source.parent)
    with pytest.raises(ValueError, match='outside|disjoint|symbolic'):
        plots.write_selected_figures(source, output)
    assert snapshot(source.parent) == before


@pytest.mark.parametrize('fault', ('missing_unplotted_source', 'changed_unplotted_source', 'missing_fixed_case',
                                  'unknown_source_file', 'payload_symlink', 'false_passed_state', 'nested_nonfinite'))
def test_selected_loader_checks_every_retained_payload_not_only_two_plots(copied_selected, fault):
    source, output = copied_selected
    selection = read_json(source/'selection.json')
    fixed = {(plots.Cell(250, 2., .5).identity, 465), (plots.Cell(250, 2., 0.).identity, 1424)}
    unplotted = next(item for item in selection['sources'] if (item['cell_id'], item['rep']) not in fixed)
    path = source/'sources'/path_diagnostics._filename(unplotted)
    if fault == 'missing_unplotted_source':
        path.unlink()
    elif fault == 'changed_unplotted_source':
        path.write_bytes(path.read_bytes()+b' ')
    elif fault == 'missing_fixed_case':
        item = next(item for item in selection['sources'] if (item['cell_id'], item['rep']) in fixed)
        (source/'sources'/path_diagnostics._filename(item)).unlink()
    elif fault == 'unknown_source_file':
        (source/'sources/extra.json').write_text('{}\n')
    elif fault == 'payload_symlink':
        data = path.read_bytes()
        path.unlink()
        other = source.parent/'foreign_payload.json'
        other.write_bytes(data)
        path.symlink_to(other)
    elif fault == 'nested_nonfinite':
        value = read_json(path)
        value['original_records']['exact']['point_diagnostics_json'] = '{"invalid":NaN}'
        write_json(path, value)
    else:
        metadata = read_json(source/'run_metadata.json')
        metadata['status'] = 'completed_diagnostic_evidence'
        write_json(source/'run_metadata.json', metadata)
    before = snapshot(source.parent)
    with pytest.raises(ValueError):
        plots.write_selected_figures(source, output)
    assert snapshot(source.parent) == before


def test_selected_unavailable_figures_are_external_and_retain_review_required_scope(copied_selected, monkeypatch):
    source, output = copied_selected
    before = snapshot(source)
    seen = []
    actual_save = plots.figure_bytes
    def save_literal(figure, stem):
        text = '\n'.join(item.get_text() for item in figure.texts)
        assert 'LITERAL VALIDATION FIXTURE' in text
        assert 'not an original E2 replay' in text
        assert all(not axis.lines and not axis.collections for axis in figure.axes)
        seen.append(stem)
        return actual_save(figure, stem)
    monkeypatch.setattr(plots, 'figure_bytes', save_literal)
    manifest = plots.write_selected_figures(source, output)
    assert seen == ['W7_selected_E2_cv0p5_rep000465', 'W7_selected_E2_cv0_rep001424']
    assert manifest['input_binding']['saved_audit_status'] == 'review_required'
    assert len(manifest['input_binding']['saved_audit_failures']) == 53
    assert len(manifest['input_binding']['source_payload_sha256']) == 53
    assert manifest['visual_review_performed'] is False
    assert not (source/'figures').exists()
    assert snapshot(source) == before
    before_output = snapshot(output), snapshot(renderer_state(output))
    for name in ('figure_bytes', 'capture', 'run_lock', 'atomic_bytes', 'atomic_json'):
        monkeypatch.setattr(plots, name, forbidden)
    assert plots.write_selected_figures(source, output) == manifest
    assert plots.verify_selected_figures(source, output) == manifest
    assert (snapshot(output), snapshot(renderer_state(output))) == before_output
    assert snapshot(source) == before

