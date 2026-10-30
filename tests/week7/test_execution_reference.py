"""Scientific composition and deterministic scheduling of Week 7 validation.

Only the exact fixed designs' separate validation namespace with R=2 runs.
The original observation counts/grids are retained. Direct references compose
the historical source_records and registered inference services, without any
substitute estimator, source randomness or numerical gate.
"""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess

import pytest

from ou_irregular.experiments.design import ESTIMATORS
from ou_irregular.experiments.execution import (
    PROCEDURE, load_records, registered_wald, source_records,
)
from ou_irregular.experiments.provenance import canonical
from ou_irregular.week4.simulation import OuSampleSource
from ou_irregular.week7.execution import require_design, run_raw
from ou_irregular.week7.floor_source import PairedFloorSource
from ou_irregular.week7.selection import floor_design, followup_design


BASELINE = '11d4a13187dfca35df1d8d9d449fb12f865a479f'
ROOT = Path(__file__).resolve().parents[2]
HISTORICAL_FILES = (
    'requirements.txt', 'requirements-lock.txt', 'preregistration.md', 'docs/plan_amendment_week3_validation.md',
    'docs/week7_exploratory_protocol.md', 'ou_irregular/ou.py', 'ou_irregular/spacing.py',
    'ou_irregular/estimators.py', 'ou_irregular/inference/contracts.py',
    'ou_irregular/inference/coordinates.py', 'ou_irregular/inference/hessian.py',
    'ou_irregular/inference/wald.py', 'ou_irregular/inference/adapters.py',
    'ou_irregular/week4/config.py', 'ou_irregular/week4/simulation.py',
    'ou_irregular/week4/summaries.py', 'ou_irregular/experiments/design.py',
    'ou_irregular/experiments/execution.py', 'ou_irregular/experiments/storage.py',
    'ou_irregular/experiments/provenance.py', 'ou_irregular/experiments/summaries.py',
    'ou_irregular/experiments/comparisons.py', 'ou_irregular/week7/selection.py',
    'ou_irregular/week7/path_diagnostics.py',
    'ou_irregular/week7/floor_source.py', 'ou_irregular/week7/floor_comparisons.py',
)
FACTORIES = {'followup': followup_design, 'floor': floor_design}
LAYOUTS = (('workers1_partition1', 1, 1),
           ('workers2_partition1', 2, 1),
           ('workers2_partition2', 2, 2))


def ordered(records):
    return sorted(records, key=lambda row: (row['cell_id'], row['rep'], ESTIMATORS.index(row['estimator'])))


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def partition_files(output):
    return {str(path.relative_to(output)): file_hash(path)
            for path in sorted(Path(output).glob('cells/*/part-*.jsonl.gz'))}


@pytest.fixture(scope='module', params=('followup', 'floor'))
def executed(request, tmp_path_factory):
    kind = request.param
    design = FACTORIES[kind]().validation(2)
    assert design.scope == 'validation'
    assert design.seed == (21261007 if kind == 'followup' else 21261009)
    directory = tmp_path_factory.mktemp('execution_reference_'+kind)
    runs = []
    for label, workers, partition_size in LAYOUTS:
        output = directory/label
        metadata = run_raw(kind, design, output, workers=workers,
                           partition_size=partition_size, progress=False)
        rows = ordered(load_records(output, design))
        runs.append({'label': label, 'output': output, 'metadata': metadata,
                     'rows': rows, 'partitions': partition_files(output),
                     'workers': workers, 'partition_size': partition_size})
    # Independent integration route: direct old source_records, with the
    # original registered services constructed exactly once per scientific cell.
    source = OuSampleSource(design.seed) if kind == 'followup' else PairedFloorSource(design.seed)
    direct = []
    for cell in design.cells:
        services = {name: registered_wald(name) for name in ESTIMATORS}
        for rep in range(2):
            direct.extend(source_records(cell, rep, source, services))
    return {'kind': kind, 'design': design, 'runs': runs, 'direct': ordered(direct)}


def test_frozen_scientific_definitions_and_validation_namespaces():
    followup, floor = followup_design(), floor_design()
    assert followup.seed == 20261007 and followup.replications == 10000
    assert followup.scope == 'exploratory'
    assert {(c.n, c.q, c.cv, c.floor_fraction) for c in followup.cells} == {
        (250, 2., 0., 1e-6), (250, 2., .5, 1e-6)}
    assert floor.seed == 20261009 and floor.replications == 2000
    assert floor.scope == 'exploratory'
    assert {(c.n, c.q, c.cv, c.floor_fraction) for c in floor.cells} == {
        (1000, .5, cv, epsilon)
        for cv, epsilon in itertools.product((0., 1.5, 2., 3.), (1e-8, 1e-6, 1e-4))}
    for kind, full in (('followup', followup), ('floor', floor)):
        # Requiring a full design here only inspects metadata; it generates no
        # observations. All calls to run_raw above use validation(2).
        assert require_design(kind, full) == full
        validation = full.validation(2)
        assert require_design(kind, validation) == validation
        assert validation.seed == full.seed+1000000
        assert validation.cells == full.cells and validation.replications == 2
        assert validation.scope == 'validation'
        assert all((c.theta, c.mu, c.sigma, c.floor_cv_threshold) == (1., 0., .5, 1.5)
                   for c in validation.cells)


def test_historical_scientific_services_and_runners_remain_byte_identical():
    for path in HISTORICAL_FILES:
        previous = subprocess.run(['git', 'show', BASELINE+':'+path], cwd=ROOT,
                                  capture_output=True, check=True).stdout
        assert (ROOT/path).read_bytes() == previous, path


def test_every_record_is_identical_across_workers_partitions_and_direct_composition(executed):
    expected = canonical(executed['direct'])
    for run in executed['runs']:
        # Compare every serialized field, including failures/candidates, raw
        # source diagnostics, optimizer statuses, full Hessians, scores,
        # covariances, SEs, interval endpoints and nested diagnostic JSON.
        assert canonical(run['rows']) == expected, run['label']
    serial, parallel, alternate = executed['runs']
    assert serial['partitions'] == parallel['partitions']
    assert set(serial['partitions']) != set(alternate['partitions'])
    for run in executed['runs']:
        for cell in executed['design'].cells:
            manifest = json.loads((run['output']/'cells'/cell.identity/'cell_manifest.json').read_text())
            assert manifest['completed_replications'] == 2 and manifest['status'] == 'complete'
            expected_ranges = [(0, 1), (1, 2)] if run['partition_size'] == 1 else [(0, 2)]
            assert [(p['start'], p['stop']) for p in manifest['partitions']] == expected_ranges
            assert sum(p['rows'] for p in manifest['partitions']) == 6


def test_complete_keys_source_addresses_and_independent_group_accounting(executed):
    kind, design = executed['kind'], executed['design']
    rows = executed['runs'][0]['rows']
    expected = set(itertools.product((cell.identity for cell in design.cells), range(2), ESTIMATORS))
    assert {(row['cell_id'], row['rep'], row['estimator']) for row in rows} == expected
    assert len(rows) == (12 if kind == 'followup' else 72)
    assert len({(row['cell_id'], row['rep']) for row in rows}) == (4 if kind == 'followup' else 24)
    cells = {cell.identity: cell for cell in design.cells}
    addresses = set()
    for row in rows:
        cell = cells[row['cell_id']]
        if kind == 'floor':
            anchor = next(c for c in design.cells if c.cv == cell.cv and c.floor_fraction == 1e-6)
            assert row['actual_cell_id'] == cell.identity
            assert row['seed_anchor_cell_id'] == anchor.identity
            assert row['seed_anchor_floor_fraction'] == 1e-6
            assert row['source_coupling_policy'] == 'week7_paired_floor_v1_shared_gamma_stationary_x0_and_standard_normal_innovations'
        else:
            anchor = cell
            assert 'actual_cell_id' not in row and 'seed_anchor_cell_id' not in row
            assert 'source_coupling_policy' not in row
        words = tuple(int(anchor.identity[i:i+8], 16) for i in range(0, 64, 8))
        expected_key = (4, *words, row['rep'])
        assert row['source_seed'] == design.seed
        assert tuple(json.loads(row['source_spawn_key'])) == expected_key
        assert tuple(json.loads(row['gap_spawn_key'])) == (*expected_key, 0)
        assert tuple(json.loads(row['path_spawn_key'])) == (*expected_key, 1)
        addresses.add((row['source_seed'], expected_key))
    assert len(addresses) == (4 if kind == 'followup' else 8)


def test_provenance_inference_policy_and_honest_raw_completion(executed):
    kind, design = executed['kind'], executed['design']
    identities = set()
    for run in executed['runs']:
        metadata, rows = run['metadata'], run['rows']
        identities.add(metadata['provenance']['identity_sha256'])
        assert metadata['provenance']['design'] == design.to_dict()
        assert metadata['provenance']['design_sha256'] == design.digest
        expected_procedure = dict(PROCEDURE)
        if kind == 'floor':
            expected_procedure['source_address'] = ('SeedSequence(root,spawn_key='
                '(4,*8_registered_floor_anchor_cell_words,rep)).spawn(2)')
        assert metadata['procedure'] == expected_procedure
        assert metadata['scientific_overrides'] == ([] if kind == 'followup' else ['paired_floor_source'])
        assert metadata['execution_kind'] == kind
        assert metadata['execution_contract'] == 'week7_guarded_cell_execution_v1'
        assert metadata['require_frozen'] is False
        assert metadata['resume_events'] == []
        assert len(metadata['execution_attempts']) == 1
        assert metadata['execution_attempts'][0]['workers'] == run['workers']
        assert metadata['execution_attempts'][0]['partition_size'] == run['partition_size']
        anchor_map = {}
        for cell in design.cells:
            anchor = (next(c for c in design.cells if c.cv == cell.cv and c.floor_fraction == 1e-6)
                      if kind == 'floor' else cell)
            anchor_map[cell.identity] = anchor.identity
        expected_policy = {
            'name': ('week7_paired_floor_v1_shared_gamma_stationary_x0_and_standard_normal_innovations'
                     if kind == 'floor' else 'registered_ou_source_v1'),
            'root_seed': design.seed, 'actual_to_seed_anchor': anchor_map,
            'seed_anchor_floor_fraction': 1e-6 if kind == 'floor' else None,
            'planned_source_realizations': 24 if kind == 'floor' else 4,
            'planned_independent_coupling_groups': 8 if kind == 'floor' else 4,
            'planned_fit_records': 72 if kind == 'floor' else 12,
            'stationary_initialization': True, 'renormalize_floored_gaps': False,
        }
        assert metadata['source_policy'] == expected_policy
        assert set(metadata['completed_cells']) == {c.identity for c in design.cells}
        for cell in design.cells:
            completion = metadata['completed_cells'][cell.identity]
            assert completion['cell_id'] == cell.identity
            assert completion['paths'] == 2 and completion['rows'] == 6
        errors = sum(bool(row['pipeline_error']) for row in rows)
        assert metadata['pipeline_errors'] == errors
        assert metadata['status'] == ('implementation_review_required' if errors else 'raw_complete_pending_audit')
        if kind == 'followup':
            assert metadata['project'] == 'ou-irregular-experiments'
        else:
            assert metadata['project'] == 'ou-week7-paired-floor'
        assert not (run['output']/'completion.json').exists()
        assert not (run['output']/'audit.json').exists()
    assert len(identities) == 1
