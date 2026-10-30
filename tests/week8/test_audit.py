"""Independent E4 audit checks: literal arithmetic, archived schedules and R=4.

The native fixture is exactly the fixed development design (28 paths/84 fits),
never the full seed. One retained fixture is shared within a test invocation.
Mutation tests change copies and refresh local digests deliberately: rejection
must come from the substantive check, not just a stale checksum. The fixture is
module-integration evidence, not the later complete-runner validation or E4 data.
"""
from copy import deepcopy
from dataclasses import replace
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

from ou_irregular.experiments import audit as core
from ou_irregular.experiments import analysis as producer_analysis
from ou_irregular.experiments import execution as engine
from ou_irregular.experiments.provenance import CORE_FILES, ROOT, canonical
from ou_irregular.week8 import audit, execution, reporting, source
from ou_irregular.week8.design import full_design, validation_design
from .test_reporting import literal_records


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def snapshot(root):
    return {str(p.relative_to(root)): ('dir' if p.is_dir() else ('link:'+str(p.readlink()))
                                      if p.is_symlink() else sha(p)) for p in root.rglob('*')}


def write_json(path, value):
    Path(path).write_text(canonical(value)+'\n')


def read_raw(output, design):
    rows = []
    for cell in design.cells:
        for path in sorted((output/'cells'/cell.identity).glob('part-*.jsonl.gz')):
            with gzip.open(path, 'rt') as stream:
                rows.extend(json.loads(line) for line in stream)
    return rows


def forbidden(*args, **kwargs):
    raise AssertionError('producer expected-data helper must not be called')


@pytest.fixture(scope='module')
def native(tmp_path_factory):
    design = validation_design()
    assert design.seed == 21261011 and design.replications == 4 and len(design.cells) == 7
    assert all(c.n == 1000 for c in design.cells)
    parent = Path(os.environ['OU_E4_AUDIT_FIXTURE_ROOT']) if 'OU_E4_AUDIT_FIXTURE_ROOT' in os.environ else tmp_path_factory.mktemp('e4_audit_native')
    output = parent/'native_archive'
    assert not output.exists(), 'native fixture cannot overwrite a previous invocation'
    execution.run_raw(design, output, extra_files=audit.FILES, workers=1, partition_size=2, progress=False)
    reporting.write_tables(design, output, extra_files=audit.FILES)
    before = snapshot(output)
    report = audit.audit_run(output, design, replay=True)
    assert snapshot(output) == before
    write_json(parent/'native_audit.json', report)
    return {'output': output, 'design': design, 'report': report, 'rows': read_raw(output, design)}


@pytest.fixture
def copy_archive(native, tmp_path):
    output = tmp_path/'mutation'
    shutil.copytree(native['output'], output)
    return output, native['design']


def refresh_table(output, name):
    path = output/'table_manifest.json'
    value = json.loads(path.read_bytes())
    value[name]['sha256'] = sha(output/name)
    write_json(path, value)


def mutate_table(output, name, callback):
    path = output/name
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        columns, rows = reader.fieldnames, list(reader)
    callback(columns, rows)
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    refresh_table(output, name)


def mutate_partition(output, cell, callback):
    root = output/'cells'/cell.identity
    path = root/'cell_manifest.json'
    manifest = json.loads(path.read_bytes())
    part = manifest['partitions'][0]
    raw = root/part['filename']
    with gzip.open(raw, 'rt') as stream:
        rows = [json.loads(line) for line in stream]
    callback(rows)
    raw.write_bytes(gzip.compress((''.join(canonical(r)+'\n' for r in rows)).encode(), mtime=0))
    part['sha256'] = sha(raw)
    write_json(path, manifest)
    # Also refresh the raw-to-table binding, so numerical/schema checks must
    # reject a coherent checksum rewrite rather than merely a stale digest.
    path = output/'table_binding.json'
    binding = json.loads(path.read_bytes())
    binding['raw_inputs']['cells'][cell.identity] = {'manifest_sha256':sha(root/'cell_manifest.json'),
                                                   'partitions':manifest['partitions']}
    binding['raw_input_sha256'] = hashlib.sha256(canonical(binding['raw_inputs']).encode()).hexdigest()
    write_json(path, binding)


@pytest.mark.parametrize('payload', ['{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}',
                                     '{"a":{"x":1,"x":2}}', '{"x":1e999}'])
def test_strict_json_rejects_ambiguous_or_nonfinite_evidence(payload):
    with pytest.raises(ValueError): audit.strict_json(payload)


@pytest.mark.parametrize('mode', [1, None, 'true'])
def test_mode_type_rejected_before_source_or_output_creation(tmp_path, mode, monkeypatch):
    monkeypatch.setattr(audit, 'independent_source', forbidden)
    with pytest.raises(ValueError, match='Boolean'):
        audit.audit_run(tmp_path/'missing', validation_design(), replay=mode)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('change', ['seed','order','scope','replications','subset'])
def test_only_complete_fixed_designs_are_accepted(tmp_path, change):
    d = validation_design()
    updates = {'seed': {'seed': d.seed+1}, 'order': {'cells':tuple(reversed(d.cells))},
               'scope':{'scope':'exploratory'}, 'replications':{'replications':3}, 'subset':{'cells':d.cells[:1]}}
    with pytest.raises(ValueError): audit.audit_run(tmp_path/'missing', replace(d, **updates[change]))
    assert list(tmp_path.iterdir()) == []


def test_source_reconstruction_matches_all_archived_sequences_without_ou_or_producer_calls(monkeypatch):
    for name in ('generate_gaps', 'gap_statistics'):
        monkeypatch.setattr(source, name, forbidden)
    monkeypatch.setattr(source.E4Source, 'sample', forbidden)
    monkeypatch.setattr(source.E4Source, 'gaps', forbidden)
    monkeypatch.setattr(audit, 'simulate_ou', forbidden)
    cells = {c.identity:c for c in full_design().cells}
    folder = ROOT/'artifacts/week8_source_validation'
    files = sorted(folder.glob('sequences-*.jsonl.gz'))
    assert len(files) == 6
    count = 0
    for path in files:
        with gzip.open(path, 'rt') as stream:
            rows = [json.loads(line) for line in stream]
        assert len(rows) == 256
        for row in rows:
            times, observed, child = audit.reconstruct_schedule(cells[row['cell_id']], row['rep'], 22261011)
            assert times.shape == (1000,) and np.all(np.diff(times)>0)
            assert set(observed) == source.SOURCE_FIELDS-{'observations_sha256'}
            for key, value in observed.items(): audit._compare(row[key], value, key)
            assert tuple(child.spawn_key)[-1] == 1
            count += 1
    assert count == 1536


def test_regular_schedule_is_exact_and_uses_no_schedule_randomness(monkeypatch):
    class NoDraws:
        random = exponential = forbidden
    monkeypatch.setattr(np.random, 'default_rng', lambda child: NoDraws())
    times, diagnostic, _ = audit.reconstruct_schedule(full_design().cells[-1], 0, 21261011)
    np.testing.assert_array_equal(times, np.arange(1000)*.5)
    assert diagnostic['source_initial_phase'] == -1
    assert diagnostic['source_phase_reset_count'] == diagnostic['source_switch_count'] == diagnostic['clip_count'] == 0


def test_literal_independent_contrasts_and_source_uncertainty_do_not_call_producer_calculators(monkeypatch):
    design = validation_design()
    cells = {c.identity:c for c in design.cells}
    rows = [r for c in design.cells for rep in range(4) for r in literal_records(c,rep)]
    monkeypatch.setattr(reporting, 'independent_comparisons', forbidden)
    monkeypatch.setattr(reporting, 'source_summary', forbidden)
    monkeypatch.setattr(producer_analysis, 'read_cell', forbidden)
    summary = []
    for cell in design.cells:
        group = [r for r in rows if r['cell_id']==cell.identity]
        summary.extend(audit._enrich(core._parameter_summaries(group,cell,design), cell))
    matched = audit.independent_contrasts(summary,design,'matched_iid')
    regular = audit.independent_contrasts(summary,design,'regular')
    assert len(matched)==27 and len(regular)==54
    for row in matched:
        assert cells[row['reference_cell_id']].mechanism=='iid'
        assert row['ratio']==row['reference_ratio']
        if row['estimator']=='exact':
            assert (row['cell_n_point_valid'],row['reference_n_point_valid'])==(3,2)
            assert row['excess_bias']==-1
            assert row['excess_bias_mcse']==pytest.approx(2/math.sqrt(3))
            assert row['rmse_ratio']==pytest.approx(math.sqrt(7/15))
            assert row['rmse_ratio_log_mcse']==pytest.approx(math.sqrt(61)/20)
            assert row['operational_coverage_difference']==.25
            assert row['operational_coverage_difference_mcse']==pytest.approx(math.sqrt(7)/8)
            assert row['conditional_coverage_difference']==.5
        if row['parameter']=='mu': assert math.isnan(row['relative_excess_bias'])
    cell = design.cells[2]
    group = [r for r in rows if r['cell_id']==cell.identity]
    result, = audit.independent_source_summary(group,cell,design)
    assert result['n_unique_sources']==4 and result['n_fit_records']==12
    assert result['raw_mean_gap_mean']==4
    assert result['raw_mean_gap_mcse']==pytest.approx(math.sqrt(5/3))
    assert result['raw_second_moment_mean']==22
    assert result['raw_second_moment_mcse']==pytest.approx(4*math.sqrt(7))
    assert result['clip_count_total']==10
    assert result['initial_phase_burst_mcse']==pytest.approx(1/math.sqrt(12))


def test_declared_source_closure_covers_imported_repository_modules():
    required = set(CORE_FILES)|set(audit.FILES)
    modules = []
    for name,module in list(sys.modules.items()):
        path = getattr(module,'__file__',None)
        if name.startswith('ou_irregular.') and path and Path(path).resolve().is_relative_to(ROOT):
            relative = str(Path(path).resolve().relative_to(ROOT))
            modules.append(relative)
    assert set(modules) <= required, set(modules)-required


def test_native_complete_source_geometry_tables_and_fixed_endpoint_replay(native):
    report = native['report']
    assert report['status']=='passed' and report['full_completion_audit'] is True
    assert report['planned_paths']==report['source_paths_verified']==28
    assert report['planned_fit_records']==84
    assert report['source_array_hashes_verified']==56
    assert report['schedule_trace_hashes_verified']==84
    assert report['replayed_source_paths']==14 and report['replayed_fits']==42
    assert {(r['cell_id'],r['replication']) for r in report['replay_records']} == {
        (c.identity,rep) for c in native['design'].cells for rep in (0,3)}
    assert report['tables_verified']['parameter_summary.csv']==63
    assert report['tables_verified']['matched_iid_comparisons.csv']==27
    assert report['tables_verified']['regular_comparisons.csv']==54
    assert report['pipeline_errors']==0


def test_readback_reconstructs_sources_without_endpoint_fitting_or_producer_expected_data(native, monkeypatch):
    for name in ('write_tables','verify_tables','source_summary','independent_comparisons','raw_input_binding'):
        monkeypatch.setattr(reporting,name,forbidden)
    monkeypatch.setattr(producer_analysis,'read_cell',forbidden)
    monkeypatch.setattr(source,'generate_gaps',forbidden)
    monkeypatch.setattr(source,'gap_statistics',forbidden)
    monkeypatch.setattr(source.E4Source,'sample',forbidden)
    monkeypatch.setattr(audit,'registered_wald',forbidden)
    before = snapshot(native['output'])
    result = audit.audit_run(native['output'],native['design'],replay=False)
    assert result['full_completion_audit'] is False and result['replay_performed'] is False
    assert result['replayed_fits']==0 and result['replay_records']==[] and result['source_paths_verified']==28
    assert snapshot(native['output'])==before


@pytest.mark.parametrize('table,column', [
    ('parameter_summary.csv','bias'), ('parameter_summary.csv','mcse_bias'),
    ('parameter_summary.csv','operational_wilson_lower'),
    ('paired_comparisons.csv','bias_difference'), ('paired_comparisons.csv','rmse_ratio'),
    ('regular_comparisons.csv','excess_bias_mcse'), ('regular_comparisons.csv','rmse_ratio'),
    ('matched_iid_comparisons.csv','excess_bias'), ('matched_iid_comparisons.csv','rmse_ratio_log_mcse'),
    ('spacing_summary.csv','realized_mean_gap_mean'), ('source_summary.csv','raw_second_moment_mcse'),
    ('source_summary.csv','source_switch_count_total'), ('outcome_counts.csv','wilson_lower'),
    ('numerical_diagnostics.csv','hessian_condition_number_p95'),
])
def test_rehashed_wrong_table_arithmetic_is_rejected_readonly(copy_archive,table,column):
    output,design = copy_archive
    def change(columns,rows):
        row = next(r for r in rows if r[column]!='')
        row[column]=str(float(row[column])+.125)
    mutate_table(output,table,change)
    before = snapshot(output)
    with pytest.raises(ValueError,match='table (arithmetic|integer) changed'):
        audit.audit_run(output,design,replay=False)
    assert snapshot(output)==before


@pytest.mark.parametrize('defect', ['duplicate_key','missing_row','extra_column','column_order','wrong_reference','missing_mu'])
def test_rehashed_structural_table_corruption_is_rejected(copy_archive,defect):
    output,design = copy_archive
    name = 'matched_iid_comparisons.csv' if defect=='wrong_reference' else 'parameter_summary.csv'
    def change(columns,rows):
        if defect=='duplicate_key': rows[-1]=dict(rows[0])
        elif defect=='missing_row': rows.pop()
        elif defect=='extra_column':
            columns.append('unknown')
            for row in rows: row['unknown']='0'
        elif defect=='column_order': columns.reverse()
        elif defect=='wrong_reference': rows[0]['reference_cell_id']=design.cells[-1].identity
        else: next(r for r in rows if r['parameter']=='mu')['relative_bias']='0'
    mutate_table(output,name,change)
    before = snapshot(output)
    with pytest.raises(ValueError): audit.audit_run(output,design,replay=False)
    assert snapshot(output)==before


@pytest.mark.parametrize('field,value', [
    ('source_seed',21261012),('raw_gaps_sha256','0'*64),('phase_trace_sha256','1'*64),
    ('source_initial_phase',1),('raw_second_moment',900.),('source_switch_count',999999),
    ('times_sha256','2'*64),('observations_sha256','3'*64),('source_spawn_key','[8,0]'),
])
def test_rehashed_raw_source_corruption_cannot_change_scientific_evidence(copy_archive,field,value):
    output,design = copy_archive
    cell=design.cells[-1] if field=='source_initial_phase' else design.cells[0]
    def change(rows):
        for row in rows:
            if row['rep']==0: row[field]=value
    mutate_partition(output,cell,change)
    before=snapshot(output)
    with pytest.raises(ValueError): audit.audit_run(output,design,replay=False)
    assert snapshot(output)==before


@pytest.mark.parametrize('defect', ['bool_flag','int_as_float','extra_raw','duplicate_row','bad_nested','covariance','endpoint','score'])
def test_rehashed_raw_inference_and_schema_corruption_rejected(copy_archive,defect):
    output,design=copy_archive
    def change(rows):
        row=rows[0]
        if defect=='bool_flag': row['point_valid']=1
        elif defect=='int_as_float': row['rep']=0.
        elif defect=='extra_raw': row['new_field']=0
        elif defect=='duplicate_row': rows[1]=deepcopy(row)
        elif defect=='bad_nested': row['point_diagnostics_json']='{"x":1,"x":2}'
        elif defect=='covariance': row['cov_theta_theta']*=2
        elif defect=='endpoint': row['theta_lower']+=.1
        else: row['standardized_score']=5.
    mutate_partition(output,design.cells[-1],change)
    before=snapshot(output)
    with pytest.raises(ValueError): audit.audit_run(output,design,replay=False)
    assert snapshot(output)==before


@pytest.mark.parametrize('defect', ['missing_cell','unknown_root','unknown_cell_file','linked_file','missing_table','bad_binding','wrong_scope'])
def test_layout_and_provenance_fail_closed_without_mutation(copy_archive,defect):
    output,design=copy_archive
    if defect=='missing_cell': shutil.rmtree(output/'cells'/design.cells[0].identity)
    elif defect=='unknown_root': (output/'extra').write_text('unexpected')
    elif defect=='unknown_cell_file': (output/'cells'/design.cells[0].identity/'extra').write_text('unexpected')
    elif defect=='linked_file':
        p=output/'parameter_summary.csv'; saved=p.read_bytes(); p.unlink()
        target=output.parent/'linked.csv'; target.write_bytes(saved); p.symlink_to(target)
    elif defect=='missing_table': (output/'spacing_summary.csv').unlink()
    elif defect=='bad_binding':
        p=output/'table_binding.json'; v=json.loads(p.read_bytes()); v['raw_input_sha256']='f'*64; write_json(p,v)
    else:
        p=output/'run_metadata.json'; v=json.loads(p.read_bytes()); v['require_frozen']=True; write_json(p,v)
    before=snapshot(output)
    with pytest.raises(ValueError): audit.audit_run(output,design,replay=False)
    assert snapshot(output)==before


def test_empty_failure_export_keeps_complete_header_and_rejects_rehashed_omission(copy_archive):
    output,design=copy_archive
    p=output/'failed_fits.csv.gz'
    with gzip.open(p,'rt') as stream: rows=list(csv.reader(stream))
    assert rows[0]==sorted(audit.RAW_FIELDS)
    # Independent membership/count checks apply whether the fixture has failures
    # or not; do not assume favorable numerical outcomes in a test fixture.
    index=rows[0].index('phase_trace_sha256')
    for row in rows: row.pop(index)
    import io
    text=io.StringIO(newline=''); csv.writer(text,lineterminator='\n').writerows(rows)
    p.write_bytes(gzip.compress(text.getvalue().encode(),mtime=0)); refresh_table(output,p.name)
    before=snapshot(output)
    with pytest.raises(ValueError,match='failure export'): audit.audit_run(output,design,replay=False)
    assert snapshot(output)==before


def test_closed_audit_requires_full_bound_report_and_readback_cannot_claim_replay(copy_archive,native):
    output,design=copy_archive
    write_json(output/'audit.json',native['report'])
    metadata=json.loads((output/'run_metadata.json').read_bytes())
    metadata.update(status='audited_pending_visual_review',audit_sha256=sha(output/'audit.json'))
    write_json(output/'run_metadata.json',metadata)
    before=snapshot(output)
    result=audit.audit_run(output,design,replay=False)
    assert not result['full_completion_audit'] and result['replayed_fits']==0
    assert snapshot(output)==before
    report=json.loads((output/'audit.json').read_bytes()); report['source_paths_verified']-=1
    write_json(output/'audit.json',report)
    metadata['audit_sha256']=sha(output/'audit.json'); write_json(output/'run_metadata.json',metadata)
    before=snapshot(output)
    with pytest.raises(ValueError,match='invariant'): audit.audit_run(output,design,replay=False)
    assert snapshot(output)==before


def test_closed_archive_cannot_use_a_partial_replay_report(copy_archive,native):
    output,design=copy_archive
    report=audit.audit_run(output,design,replay=False)
    write_json(output/'audit.json',report)
    metadata=json.loads((output/'run_metadata.json').read_bytes())
    metadata.update(status='audited_pending_visual_review',audit_sha256=sha(output/'audit.json'))
    write_json(output/'run_metadata.json',metadata)
    before=snapshot(output)
    with pytest.raises(ValueError,match='full endpoint audit'): audit.audit_run(output,design,replay=False)
    assert snapshot(output)==before


def test_endpoint_replay_compares_complete_nested_diagnostics(copy_archive,monkeypatch):
    output,design=copy_archive
    original=audit.source_records
    def changed(*args,**kwargs):
        rows=original(*args,**kwargs)
        diag=json.loads(rows[0]['point_diagnostics_json']); diag['injected_extra']=1
        rows[0]['point_diagnostics_json']=canonical(diag)
        return rows
    monkeypatch.setattr(audit,'source_records',changed)
    before=snapshot(output)
    with pytest.raises(ValueError,match='endpoint'): audit.audit_run(output,design,replay=True)
    assert snapshot(output)==before


def test_invalid_point_candidates_and_nonempty_failure_exports_are_not_discarded(native,tmp_path):
    """Construct an explicit invalid-point unit case from a retained real fit.

    This does not alter the retained archive or an estimator's decision. Only
    the unit input's eligibility is changed; all candidate numbers remain.
    """
    from ou_irregular.inference.contracts import PointFit
    cell=native['design'].cells[-1]
    sample,diagnostic=audit.independent_source(cell,0,native['design'].seed)
    failures=[]
    for method in audit.METHODS:
        saved=next(r for r in native['rows'] if r['cell_id']==cell.identity and r['rep']==0 and r['estimator']==method)
        pointdiag=json.loads(saved['point_diagnostics_json'])
        pointdiag.update(valid_for_point_summary=False,fit_status='boundary_suspected')
        point=PointFit(estimator=method,theta=saved['theta'],mu=saved['mu'],sigma=saved['sigma'],nll=saved['nll'],
                       valid=False,optimizer_success=saved['optimizer_success'],status='boundary_suspected',
                       reason='constructed invalid-point unit input',diagnostics=pointdiag)
        assert point.valid is False and point.optimizer_success is saved['optimizer_success']
        result=engine.registered_wald(method).infer(sample,point_fit=point)
        row=json.loads(canonical({**cell.to_dict(),**diagnostic,'rep':0,'source_valid':True,
                                  'source_error':'','pipeline_error':'',**engine.inference_record(result)}))
        audit._validate_row(row,cell,0,method)
        assert not row['point_valid'] and not row['interval_valid']
        assert row['theta']==saved['theta'] and row['nll']==saved['nll']
        assert row['theta_lower'] is None and row['theta_se'] is None
        failures.append(row)
    path=tmp_path/'failed.csv.gz'
    def save(rows):
        with gzip.open(path,'wt',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=sorted(audit.RAW_FIELDS),lineterminator='\n')
            writer.writeheader();writer.writerows(rows)
    failures.sort(key=lambda r:r['estimator'])
    save(failures)
    with audit.FailureExport(path) as checker:
        checker.check_rows(failures); assert checker.finish(3)==3
    save(failures[:-1])
    with pytest.raises(ValueError,match='omitted'):
        with audit.FailureExport(path) as checker: checker.check_rows(failures)
    save(list(reversed(failures)))
    with pytest.raises(ValueError):
        with audit.FailureExport(path) as checker: checker.check_rows(failures)


def test_constructed_unavailable_interval_is_retained_but_inconsistent_gate_rejected(native):
    """A fabricated rejected-interval case, not a modified scientific outcome."""
    cell=native['design'].cells[-1]
    row=deepcopy(next(r for r in native['rows'] if r['cell_id']==cell.identity and r['rep']==0 and r['estimator']=='exact'))
    row.update(interval_valid=False,interval_status='unstable_hessian',interval_reason='constructed unit rejection',
               hessian_relative_difference=.006,hessian_condition_number=None,standardized_score=None)
    for k in row:
        if k.startswith(('cov_','working_cov_')) or k.endswith(('_lower','_upper','_se')): row[k]=None
    diag=json.loads(row['interval_diagnostics_json'])
    diag.update(hessian_status='unstable_hessian',hessian_relative_difference=.006)
    for key in ('standard_errors','critical_value','hessian_eigenvalues','hessian_condition_number','standardized_score'): diag.pop(key,None)
    row['interval_diagnostics_json']=canonical(diag)
    audit._validate_row(row,cell,0,'exact')
    assert row['point_valid'] is True and row['interval_valid'] is False
    row['hessian_relative_difference']=.001;diag['hessian_relative_difference']=.001
    row['interval_diagnostics_json']=canonical(diag)
    with pytest.raises(ValueError,match='without a recorded failed'): audit._validate_row(row,cell,0,'exact')


def test_saved_full_endpoint_descriptors_are_checked_even_when_fitting_is_skipped(copy_archive,native):
    output,design=copy_archive
    record=deepcopy(native['report']);record['replay_records'][0]['observations_sha256']='0'*64
    write_json(output/'audit.json',record)
    metadata=json.loads((output/'run_metadata.json').read_bytes())
    metadata.update(status='audited_pending_visual_review',audit_sha256=sha(output/'audit.json'))
    write_json(output/'run_metadata.json',metadata)
    before=snapshot(output)
    with pytest.raises(ValueError,match='endpoint descriptors'): audit.audit_run(output,design,replay=False)
    assert snapshot(output)==before
