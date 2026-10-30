"""Independent accounting boundaries for the final Week 7 orchestrator.

The shared lifecycle fixtures are explicitly literal: durable partitions and
table composition are real, but the source records, numerical audit outcomes,
and figure stage outcomes are test inputs. No OU observations, scientific
random draws, fit replays, or selected-source diagnostics are generated here.
Such fixtures are never independent evidence of scientific completion.
"""
from copy import deepcopy
from dataclasses import replace
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from ou_irregular.experiments import execution as historical_execution
from ou_irregular.experiments.design import Cell, ESTIMATORS
from ou_irregular.experiments.provenance import canonical, file_digest
from ou_irregular.inference.adapters import Week3Estimator
from ou_irregular.week4.simulation import OuSampleSource
from ou_irregular.week7 import execution, floor_audit, path_diagnostics, reporting
from ou_irregular.week7.floor_source import COUPLING_POLICY, PairedFloorSource
from ou_irregular.week7.selection import floor_design, followup_design


PARAMETERS = ('theta', 'mu', 'sigma')


def snapshot(directory):
    """Literal byte map, deliberately independent of production fingerprints."""
    root = Path(directory)
    return {str(path.relative_to(root)): path.read_bytes()
            for path in sorted(root.rglob('*')) if path.is_file()}


def state(directory):
    return json.loads((Path(directory)/'run_metadata.json').read_text())


def write_json(path, value):
    Path(path).write_text(canonical(value)+'\n')


def forbidden(*args, **kwargs):
    raise AssertionError('literal runner references must not simulate or replay scientific inputs')


def literal_records(design, kind, *, all_valid=False):
    """Complete records with a known six-record validity pattern per cell.

    For R=2 each cell has six records: five valid points, three valid intervals,
    and three retained failures. The unavailable point has a finite candidate;
    unavailable intervals have finite covering endpoints to ensure that their
    recorded validity remains part of all-R coverage. Numerical diagnostics
    are descriptive fixture values, not valid inference replay evidence.
    """
    rows = []
    for cell in design.cells:
        anchor = replace(cell, floor_fraction=1e-6) if kind == 'floor' else cell
        words = [int(anchor.identity[i:i+8], 16) for i in range(0, 64, 8)]
        physical = anchor.identity if kind == 'floor' and cell.cv == 0 else cell.identity
        floor_applied = cell.cv >= cell.floor_cv_threshold
        mean = cell.q/cell.theta
        for rep in range(design.replications):
            key = [4, *words, rep]
            source = {
                **cell.to_dict(), 'rep': rep, 'source_seed': design.seed,
                'source_spawn_key': canonical(key), 'gap_spawn_key': canonical([*key, 0]),
                'path_spawn_key': canonical([*key, 1]),
                'times_sha256': hashlib.sha256(f'literal:times:{physical}:{rep}'.encode()).hexdigest(),
                'observations_sha256': hashlib.sha256(f'literal:states:{physical}:{rep}'.encode()).hexdigest(),
                'source_valid': True, 'source_error': '', 'pipeline_error': '',
                'family': 'equidistant' if cell.cv == 0 else 'gamma',
                'target_mean_gap': mean, 'requested_cv': cell.cv,
                'floor_applied': floor_applied,
                'floor_value': cell.floor_fraction*mean if floor_applied else 0.,
                'clip_count': 0, 'clip_fraction': 0., 'raw_mean_gap': mean,
                'realized_mean_gap': mean, 'raw_cv': cell.cv, 'realized_cv': cell.cv,
                'realized_span': (cell.n-1)*mean,
                'min_gap': mean if cell.cv == 0 else mean/2,
                'expected_mean_gap_after_floor': mean, 'expected_clip_fraction': 0.,
            }
            if kind == 'floor':
                source.update(actual_cell_id=cell.identity, seed_anchor_cell_id=anchor.identity,
                              seed_anchor_floor_fraction=1e-6, source_coupling_policy=COUPLING_POLICY)
            for method in ESTIMATORS:
                point = bool(all_valid or not (rep % 2 == 0 and method == 'euler'))
                interval = bool(all_valid or (point and
                    ((rep % 2 == 0 and method == 'exact') or
                     (rep % 2 == 1 and method != 'exact'))))
                row = {**source, 'estimator': method,
                       'point_valid': point, 'optimizer_success': point,
                       'point_status': 'valid_interior' if point else 'boundary',
                       'point_reason': 'literal finite candidate fixture',
                       'interval_valid': interval,
                       'interval_status': ('valid' if interval else
                                           'step_disagreement' if point else 'invalid_point_fit'),
                       'interval_reason': 'literal interval availability fixture',
                       'nll': 1.,
                       'point_diagnostics_json': '{"fixture":"literal"}',
                       'interval_diagnostics_json': '{"fixture":"literal"}',
                       'covariance_type': 'inverse_observed_hessian_model_based',
                       'hessian_condition_number': 3. if point else np.nan,
                       'hessian_relative_difference': .001 if interval else .02 if point else np.nan,
                       'hessian_whitened_difference': .001 if interval else .02 if point else np.nan,
                       'standardized_score': .001 if point else np.nan}
                for parameter, truth in zip(PARAMETERS, (1., 0., .5)):
                    row.update({parameter: truth+.1*(rep+1),
                                parameter+'_se': .5 if interval else np.nan,
                                parameter+'_lower': truth-1., parameter+'_upper': truth+1.})
                # These uncomputed arrays are deliberately unavailable. The
                # audit callback below is a stage-boundary double, never a
                # claim that literal records pass a numerical Hessian audit.
                for i, parameter in enumerate(PARAMETERS):
                    row[f'gradient_{i}'] = np.nan
                    for j in range(i, 3):
                        row[f'cov_{parameter}_{PARAMETERS[j]}'] = np.nan
                        row[f'working_cov_{i}_{j}'] = np.nan
                        row[f'working_hessian_{i}_{j}'] = np.nan
                rows.append(row)
    return rows


def expected_counts(design):
    """Count the literal pattern without applying any production summarizer."""
    even = (design.replications+1)//2
    odd = design.replications//2
    cells = len(design.cells)
    return {'planned_paths': cells*design.replications,
            'planned_fit_records': 3*cells*design.replications,
            'point_valid_records': cells*(2*even+3*odd),
            'interval_valid_records': cells*(even+2*odd),
            'retained_failure_records': cells*(2*even+odd)}


def write_failure_csv(output, rows, columns):
    """Build a literal CSV using csv.writer so duplicate headers survive."""
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, lineterminator='\n')
    writer.writerow(columns)
    for row in rows:
        writer.writerow([row.get(name, '') for name in columns])
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output/'failed_fits.csv.gz').write_bytes(gzip.compress(stream.getvalue().encode(), mtime=0))


@pytest.mark.parametrize('kind', ('followup', 'floor'))
@pytest.mark.parametrize('validation_replications', (None, 1, 4))
def test_runner_designs_retain_fixed_protocol_and_separate_validation_namespace(
        kind, validation_replications):
    from ou_irregular.week7 import runner
    design = runner.design_for(kind, validation_replications=validation_replications)
    expected_cells = ([(250, 2., cv, 1e-6) for cv in (0., .5)] if kind == 'followup'
                      else [(1000, .5, cv, floor) for cv in (0., 1.5, 2., 3.)
                            for floor in (1e-8, 1e-6, 1e-4)])
    assert [(cell.n, cell.q, cell.cv, cell.floor_fraction) for cell in design.cells] == expected_cells
    assert all((cell.theta, cell.mu, cell.sigma, cell.floor_cv_threshold) == (1., 0., .5, 1.5)
               for cell in design.cells)
    full_name = 'week7_e3_followup' if kind == 'followup' else 'week7_floor_sensitivity'
    full_seed = 20261007 if kind == 'followup' else 20261009
    if validation_replications is None:
        assert (design.name, design.scope, design.seed, design.replications) == (
            full_name, 'exploratory', full_seed, 10000 if kind == 'followup' else 2000)
    else:
        assert (design.name, design.scope, design.seed, design.replications) == (
            full_name+'_validation', 'validation', full_seed+1000000, validation_replications)
    assert tuple(design.to_dict()['estimators']) == ('exact', 'pfml', 'euler')


@pytest.mark.parametrize('replications', (0, -1, True, 1.5))
def test_validation_replication_count_is_a_positive_integer(replications):
    from ou_irregular.week7 import runner
    with pytest.raises(ValueError):
        runner.design_for('followup', validation_replications=replications)


@pytest.mark.parametrize('kind', ('followup', 'floor'))
@pytest.mark.parametrize('empty', (False, True))
def test_exact_failure_header_survives_even_when_no_fits_failed(tmp_path, kind, empty):
    from ou_irregular.week7 import runner
    design = (followup_design() if kind == 'followup' else floor_design()).validation(2)
    rows = literal_records(design, kind, all_valid=empty)
    columns = sorted(rows[0])
    assert len(columns) == (84 if kind == 'followup' else 88)
    failures = [row for row in rows if not row['point_valid'] or not row['interval_valid']]
    assert bool(failures) is not empty
    write_failure_csv(tmp_path, failures, columns)
    before = snapshot(tmp_path)
    runner.strict_failure_header(tmp_path, kind)
    assert snapshot(tmp_path) == before


@pytest.mark.parametrize('corruption', ('duplicate', 'missing', 'extra', 'reordered'))
@pytest.mark.parametrize('empty', (False, True))
def test_followup_failure_header_rejects_set_equivalent_or_incomplete_csv(
        tmp_path, corruption, empty):
    from ou_irregular.week7 import runner
    rows = literal_records(followup_design().validation(2), 'followup', all_valid=empty)
    columns = sorted(rows[0])
    if corruption == 'duplicate':
        columns.append(columns[0])
    elif corruption == 'missing':
        columns.remove('covariance_type')
    elif corruption == 'extra':
        columns.append('literal_extra_column')
    else:
        columns[:2] = columns[1::-1]
    failures = [row for row in rows if not row['point_valid'] or not row['interval_valid']]
    write_failure_csv(tmp_path, failures, columns)
    before = snapshot(tmp_path)
    with pytest.raises(ValueError):
        runner.strict_failure_header(tmp_path, 'followup')
    assert snapshot(tmp_path) == before
