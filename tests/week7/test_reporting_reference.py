"""Saved-table scientific references using literal, non-simulated raw fixtures.

These fixtures bind reporting dependencies before storage, but their numerical
records and descriptive source hashes are fabricated unit-test inputs. They
are never represented as simulated experiments or source-replay evidence.
"""
import gc
import hashlib
import json
import shutil
import weakref

import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments import execution as historical_execution
from ou_irregular.experiments.analysis import csv_bytes, write_tables as historical_tables
from ou_irregular.experiments.design import ESTIMATORS
from ou_irregular.experiments.execution import load_records
from ou_irregular.inference.adapters import Week3Estimator
from ou_irregular.week4.simulation import OuSampleSource
from ou_irregular.week7 import execution, reporting
from ou_irregular.week7.floor_comparisons import COLUMNS, epsilon_comparisons
from ou_irregular.week7.floor_source import PairedFloorSource
from ou_irregular.week7.selection import floor_design, followup_design

from .test_floor_comparisons_reference import fabricated_records


def close(value):
    return pytest.approx(value, rel=1e-12, abs=2e-14)


def _augment(record):
    """Literal spacing/diagnostic fields required by the unchanged table writer."""
    row = dict(record)
    mean, cv, rep = row['mean_gap'], row['cv'], row['rep']
    floored = cv >= row['floor_cv_threshold']
    row.update(family='equidistant' if cv == 0 else 'gamma', target_mean_gap=mean,
               requested_cv=cv, floor_applied=floored,
               floor_value=row['floor_fraction']*mean if floored else 0.,
               clip_count=0, clip_fraction=0., raw_mean_gap=mean, realized_mean_gap=mean,
               raw_cv=cv, realized_cv=cv, realized_span=(row['n']-1)*mean,
               min_gap=mean if cv == 0 else mean/2,
               expected_mean_gap_after_floor=mean, expected_clip_fraction=0.,
               point_reason='literal arithmetic fixture', interval_reason='literal interval status',
               candidate_note='literal candidate retained without replacing it by zero')
    if row['point_valid']:
        row.update(hessian_condition_number=float(3+rep) if row['interval_valid'] else float(200+rep),
                   hessian_relative_difference=.001 if row['interval_valid'] else .01,
                   hessian_whitened_difference=.001 if row['interval_valid'] else .02,
                   standardized_score=.001 if row['interval_valid'] else .02)
    else:
        row.update(hessian_condition_number=np.nan, hessian_relative_difference=np.nan,
                   hessian_whitened_difference=np.nan, standardized_score=np.nan)
    return row


def _floor_literal_records():
    _, records = fabricated_records(
        [1., 2., 4., 100., 999.], [2., 1., 3., 999., 20.],
        point_a=[True, True, True, True, False], point_b=[True, True, True, False, True],
        interval_a=[True, True, True, False, False], interval_b=[True, False, True, False, True],
        cover_a=[True, True, False, False, False], cover_b=[True, False, True, False, True],
        cvs=(0., 1.5, 2., 3.))
    # One entire estimator/variant has no valid points. It must survive every
    # table with missing point metrics, while all-R coverage remains defined.
    unavailable = ((records.cv == 3.) & (records.floor_fraction == 1e-4)
                   & (records.estimator == 'euler'))
    records.loc[unavailable, ['point_valid', 'optimizer_success', 'interval_valid']] = False
    records.loc[unavailable, 'point_status'] = 'boundary_supremum'
    records.loc[unavailable, 'interval_status'] = 'invalid_point_fit'
    for p in ('theta', 'mu', 'sigma'):
        records.loc[unavailable, p+'_se'] = np.nan
    return [_augment(row) for row in records.to_dict(orient='records')]


def _followup_literal_records(design):
    rows = []
    for cell in design.cells:
        words = tuple(int(cell.identity[i:i+8], 16) for i in range(0, 64, 8))
        for rep in range(2):
            key = (4, *words, rep)
            valid = cell.cv == 0.
            for name in ESTIMATORS:
                row = {**cell.to_dict(), 'rep': rep, 'estimator': name,
                       'source_seed': design.seed,
                       'source_spawn_key': json.dumps(key, separators=(',', ':')),
                       'gap_spawn_key': json.dumps((*key, 0), separators=(',', ':')),
                       'path_spawn_key': json.dumps((*key, 1), separators=(',', ':')),
                       'times_sha256': hashlib.sha256(f'literal-times-{cell.identity}-{rep}'.encode()).hexdigest(),
                       'observations_sha256': hashlib.sha256(f'literal-states-{cell.identity}-{rep}'.encode()).hexdigest(),
                       'source_valid': True, 'source_error': '', 'pipeline_error': '',
                       'point_valid': valid, 'optimizer_success': valid,
                       'point_status': 'valid_interior' if valid else 'boundary_supremum',
                       'interval_valid': valid, 'interval_status': 'valid' if valid else 'invalid_point_fit',
                       'nll': 0.}
                for p, truth in zip(('theta', 'mu', 'sigma'), (1., 0., .5)):
                    error = float(rep+1) if valid else float(rep+9)
                    row.update({p: truth+error, p+'_se': 2. if valid else np.nan,
                                p+'_lower': truth-.125, p+'_upper': truth+2*error+.125})
                rows.append(_augment(row))
    return rows


def _raw_snapshot(output):
    paths = [output/'run_metadata.json', output/'resolved_design.json']
    paths.extend(p for p in (output/'cells').rglob('*') if p.is_file())
    return {str(path.relative_to(output)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(paths)}


def _forbidden(*args, **kwargs):
    raise AssertionError('reporter must not generate observations or refit any estimator')


@pytest.fixture(scope='module', params=('followup', 'floor'))
def saved_report(request, tmp_path_factory):
    kind = request.param
    design = followup_design().validation(2) if kind == 'followup' else floor_design().validation(5)
    literal = _followup_literal_records(design) if kind == 'followup' else _floor_literal_records()
    by_source = {}
    for row in literal:
        by_source.setdefault((row['cell_id'], row['rep']), []).append(row)
    directory = tmp_path_factory.mktemp('reporting_literal_reference_'+kind)
    output = directory/'guarded_report'
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(historical_execution, 'source_records',
                      lambda cell, rep, source, services: [dict(row) for row in by_source[cell.identity, rep]])
        patch.setattr(OuSampleSource, 'sample', _forbidden)
        patch.setattr(PairedFloorSource, 'sample', _forbidden)
        patch.setattr(Week3Estimator, 'fit', _forbidden)
        execution.run_raw(kind, design, output, extra_files=reporting.FILES,
                          workers=1, partition_size=2, progress=False)
    records = list(load_records(output, design))
    legacy = directory/'unchanged_standard_reference'
    shutil.copytree(output, legacy)
    before = _raw_snapshot(output)
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(OuSampleSource, 'sample', _forbidden)
        patch.setattr(PairedFloorSource, 'sample', _forbidden)
        patch.setattr(Week3Estimator, 'fit', _forbidden)
        patch.setattr(historical_execution, 'source_records', _forbidden)
        patch.setattr(execution, 'run_raw', _forbidden)
        result = reporting.write_tables(kind, design, output)
        reference = historical_tables(legacy, design)
    assert _raw_snapshot(output) == before
    return {'kind': kind, 'design': design, 'output': output, 'records': records,
            'result': result, 'legacy': legacy, 'standard_reference': reference,
            'raw_before': before}


def test_standard_composition_and_every_raw_byte_are_preserved(saved_report):
    case = saved_report
    assert case['result']['standard'] == case['standard_reference']
    assert set(case['result']['standard']) == set(reporting.STANDARD_FILES)
    for name in (*reporting.STANDARD_FILES, 'table_manifest.json'):
        assert (case['output']/name).read_bytes() == (case['legacy']/name).read_bytes(), name
    assert _raw_snapshot(case['output']) == case['raw_before']
    metadata = json.loads((case['output']/'run_metadata.json').read_text())
    assert metadata['status'] == 'raw_complete_pending_audit'
    assert not (case['output']/'audit.json').exists()
    assert not (case['output']/'figures').exists()
    assert not (case['output']/'visual_review.json').exists()
    if case['kind'] == 'followup':
        assert case['result']['epsilon'] is None
        assert not (case['output']/reporting.EPSILON_FILE).exists()
        assert not (case['output']/reporting.EPSILON_MANIFEST).exists()


def test_standard_counts_missing_values_and_failure_candidates_survive(saved_report):
    case = saved_report
    raw = pd.DataFrame(case['records'])
    summary = pd.read_csv(case['output']/'parameter_summary.csv', float_precision='round_trip')
    assert len(summary) == 9*len(case['design'].cells)
    assert summary[summary.parameter == 'mu'].relative_bias.isna().all()
    assert summary[summary.parameter == 'mu'].relative_rmse.isna().all()
    assert summary.n_planned.eq(case['design'].replications).all()
    missing_cell = (next(c for c in case['design'].cells if c.cv == .5).identity
                    if case['kind'] == 'followup' else
                    next(c for c in case['design'].cells if c.cv == 3. and c.floor_fraction == 1e-4).identity)
    missing = summary[(summary.cell_id == missing_cell)
                      & (True if case['kind'] == 'followup' else summary.estimator.eq('euler'))]
    assert len(missing) == (9 if case['kind'] == 'followup' else 3)
    assert missing.n_point_valid.eq(0).all() and missing.n_interval_valid.eq(0).all()
    for column in ('bias', 'rmse', 'mcse_bias', 'mcse_rmse', 'coverage_conditional',
                   'empirical_sd_point_valid', 'mean_reported_se'):
        assert missing[column].isna().all(), column
    assert missing.coverage_operational.eq(0.).all()
    assert missing.operational_mcse.eq(0.).all()
    failed = pd.read_csv(case['output']/'failed_fits.csv.gz', float_precision='round_trip')
    mask = ~raw.point_valid | ~raw.interval_valid | raw.pipeline_error.ne('')
    keys = ['cell_id', 'rep', 'estimator']
    assert set(failed[keys].itertuples(index=False, name=None)) == set(raw.loc[mask, keys].itertuples(index=False, name=None))
    assert len(failed) == (6 if case['kind'] == 'followup' else 75)
    expected = raw.loc[mask].sort_values(keys).reset_index(drop=True)
    actual = failed.sort_values(keys).reset_index(drop=True)
    for name in ('theta', 'mu', 'sigma', 'nll', 'theta_lower', 'theta_upper',
                 'sigma_lower', 'sigma_upper', 'theta_se', 'point_valid', 'interval_valid'):
        np.testing.assert_array_equal(actual[name].to_numpy(), expected[name].to_numpy())
    assert actual.candidate_note.eq('literal candidate retained without replacing it by zero').all()
    outcomes = pd.read_csv(case['output']/'outcome_counts.csv')
    assert outcomes.groupby(['cell_id', 'estimator'])['count'].sum().eq(case['design'].replications).all()


def test_saved_standard_metrics_use_own_valid_fits_and_keep_failed_hessian_diagnostics(saved_report):
    if saved_report['kind'] != 'floor':
        pytest.skip('paired-floor reference applies only to the floor fixture')
    output, design = saved_report['output'], saved_report['design']
    cell = next(c for c in design.cells if c.cv == 1.5 and c.floor_fraction == 1e-8)
    summary = pd.read_csv(output/'parameter_summary.csv', float_precision='round_trip')
    row = summary[(summary.cell_id == cell.identity) & (summary.estimator == 'exact')
                  & (summary.parameter == 'theta')].iloc[0]
    assert row.n_point_valid == 4 and row.n_interval_valid == 3
    assert row.bias == close(107/4)
    assert row.rmse == close(np.sqrt((1+4+16+10000)/4))
    assert row.coverage_operational == .4
    assert row.coverage_conditional == close(2/3)
    diagnostics = pd.read_csv(output/'numerical_diagnostics.csv', float_precision='round_trip')
    value = diagnostics[(diagnostics.cell_id == cell.identity) & (diagnostics.estimator == 'exact')].iloc[0]
    assert value.n_interval_valid == 3
    assert value.hessian_condition_number_n_finite == 4
    assert value.hessian_condition_number_n_nonfinite == 1
    assert value.hessian_condition_number_max == 203.
    assert value.hessian_condition_number_median == 4.5


def test_bounded_epsilon_iterator_releases_raw_groups_and_matches_full_arithmetic(saved_report, monkeypatch):
    if saved_report['kind'] != 'floor':
        pytest.skip('paired-floor reference applies only to the floor fixture')
    output, design = saved_report['output'], saved_report['design']
    metadata = json.loads((output/'run_metadata.json').read_text())
    identity = metadata['provenance']['identity_sha256']
    reader, comparator = reporting.read_cell, reporting.epsilon_comparisons
    loaded, raw_refs, comparison_inputs = [], [], []
    def tracked_read(output, design, cell, identity):
        frame = reader(output, design, cell, identity)
        loaded.append((cell.cv, cell.floor_fraction, len(frame)))
        raw_refs.append(weakref.ref(frame))
        return frame
    def tracked_compare(records, partial_design):
        assert len(partial_design.cells) == 3
        assert len({c.cv for c in partial_design.cells}) == 1
        assert len(records) == 3*3*design.replications
        comparison_inputs.append(weakref.ref(records))
        return comparator(records, partial_design)
    monkeypatch.setattr(reporting, 'read_cell', tracked_read)
    monkeypatch.setattr(reporting, 'epsilon_comparisons', tracked_compare)
    pieces = []
    for index, piece in enumerate(reporting.iter_epsilon_tables(output, design, identity)):
        assert len(piece) == 18 and tuple(piece.columns) == COLUMNS
        assert len(loaded) == (index+1)*3
        gc.collect()
        assert all(reference() is None for reference in raw_refs)
        assert all(reference() is None for reference in comparison_inputs)
        pieces.append(piece)
    assert len(pieces) == 4
    assert loaded == [(cv, epsilon, 3*design.replications)
                      for cv in (0., 1.5, 2., 3.) for epsilon in (1e-8, 1e-6, 1e-4)]
    combined = pd.concat(pieces, ignore_index=True)
    full = epsilon_comparisons(saved_report['records'], design)
    assert csv_bytes(combined) == csv_bytes(full)
    assert (output/reporting.EPSILON_FILE).read_bytes() == csv_bytes(full)


def test_saved_epsilon_hand_references_missing_domains_and_complete_schema(saved_report):
    if saved_report['kind'] != 'floor':
        pytest.skip('paired-floor reference applies only to the floor fixture')
    case = saved_report
    table = pd.read_csv(case['output']/reporting.EPSILON_FILE, float_precision='round_trip')
    assert len(table) == 72 and tuple(table.columns) == COLUMNS and len(COLUMNS) == 71
    row = table[(table.cv == 1.5) & (table.floor_fraction == 1e-8)
                & (table.estimator == 'exact') & (table.parameter == 'theta')].iloc[0]
    assert row.n_planned == 5 and row.n_common_point_valid == 3
    assert row.bias_a_common == close(7/3) and row.bias_b_common == 2.
    assert row.bias_difference == close(1/3) and row.bias_difference_mcse == close(2/3)
    assert row.mse_difference == close(7/3) and row.mse_difference_mcse == close(np.sqrt(76)/3)
    assert row.rmse_ratio == close(np.sqrt(1.5)) and row.rmse_ratio_log_mcse == close(5/28)
    assert row.operational_coverage_difference == close(-.2)
    assert row.operational_coverage_difference_mcse == close(np.sqrt(.14))
    assert (row.n_both_cover, row.n_a_only_cover, row.n_b_only_cover, row.n_neither_cover) == (1, 1, 2, 1)
    missing = table[(table.cv == 3.) & (table.floor_fraction == 1e-4) & (table.estimator == 'euler')]
    assert len(missing) == 3 and missing.n_common_point_valid.eq(0).all()
    for column in ('bias_a_common', 'bias_b_common', 'bias_difference', 'mse_difference',
                   'rmse_difference', 'rmse_ratio', 'rmse_ratio_mcse'):
        assert missing[column].isna().all(), column
    assert missing.operational_coverage_difference.eq(-.6).all()
    np.testing.assert_allclose(missing.operational_coverage_difference_mcse, np.sqrt(.06), rtol=1e-12)
    relative = [name for name in table if name.startswith('relative_')]
    assert table.loc[table.parameter == 'mu', relative].isna().all().all()
    assert table.loc[table.cv == 0., 'rmse_ratio'].eq(1.).all()
    assert table.loc[table.cv == 0., 'rmse_ratio_mcse'].eq(0.).all()
    assert not table.floor_fraction.eq(1e-6).any()
    manifest = case['result']['epsilon']
    assert manifest['table']['columns'] == list(COLUMNS) and manifest['table']['rows'] == 72
    assert set(manifest['raw_inputs']['cells']) == {c.identity for c in case['design'].cells}
    assert manifest['raw_inputs']['source_policy']['planned_source_realizations'] == 60
    assert manifest['raw_inputs']['source_policy']['planned_independent_coupling_groups'] == 20
    assert manifest['raw_inputs']['source_policy']['planned_fit_records'] == 180
    assert 'run_metadata_sha256' not in manifest['raw_inputs']
    assert set(case['result']['standard']) == set(reporting.STANDARD_FILES)
    assert reporting.EPSILON_FILE not in case['result']['standard']
