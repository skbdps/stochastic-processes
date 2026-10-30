"""Deterministic E3 selection from sealed headline evidence, with no reruns.

Selection is outcome-informed and explicitly exploratory. The rule is fixed
before affected data are generated; it must not be changed to favor a result.
Every failed source is retained. Extremum ties use scientific ID then replication.
"""
from collections import defaultdict
import gzip
import json
import math
from pathlib import Path

import pandas as pd

from ..experiments.audit import artifact_fingerprint, semantic_metadata_digest
from ..experiments.design import Cell, RunDesign
from ..experiments.provenance import ROOT, file_digest, verify
from ..experiments.workflow import load_saved_design, verify_figures
from ..week5.config import load_config as load_e1
from ..week6.config import load_config as load_e2


CONTROL_METRICS = ('hessian_whitened_difference', 'standardized_score')


def select_records(records):
    """Select failed sources and accepted extremes, retaining tie determinism."""
    reasons = defaultdict(set)
    maxima = {}
    failures = []
    seen = set()
    for row in records:
        key = (row['cell_id'], int(row['rep']))
        fit_key = (*key, row['estimator'])
        if fit_key in seen:
            raise ValueError('duplicate source/estimator key during selection')
        seen.add(fit_key)
        if not row['point_valid'] or not row['interval_valid'] or row.get('pipeline_error'):
            reasons[key].add('unavailable:' + row['estimator'] + ':' + row['interval_status'])
            failures.append(dict(cell_id=key[0], rep=key[1], estimator=row['estimator'],
                                 point_status=row['point_status'], interval_status=row['interval_status']))
        if row['interval_valid']:
            for metric in CONTROL_METRICS:
                value = row[metric]
                if value is not None and math.isfinite(float(value)):
                    group = (row['estimator'], metric)
                    candidate = (-float(value), *key)
                    if group not in maxima or candidate < maxima[group]:
                        maxima[group] = candidate
    for (method, metric), (_, cid, rep) in maxima.items():
        reasons[cid, rep].add('accepted_maximum:' + method + ':' + metric)
    return ([dict(cell_id=cid, rep=rep, reasons=sorted(why))
             for (cid, rep), why in sorted(reasons.items())],
            sorted(failures, key=lambda x: (x['cell_id'], x['rep'], x['estimator'])))


def sealed_input(run_dir):
    """Verify the seal and all retained bytes without regenerating observations."""
    path = Path(run_dir)
    metadata = json.loads((path / 'run_metadata.json').read_text())
    design = load_saved_design(path)
    audit = json.loads((path / 'audit.json').read_text())
    if metadata['status'] != 'completed' or audit['status'] != 'passed' or not audit['full_completion_audit']:
        raise ValueError('selection requires completed, fully audited headlines')
    verify(metadata['provenance'], expected_design=design)
    if (file_digest(path / 'audit.json') != metadata['audit_sha256']
            or semantic_metadata_digest(metadata) != audit['audited_run_metadata_sha256']
            or artifact_fingerprint(path) != audit['audited_artifact_sha256']):
        raise ValueError('headline seal/audit binding changed')
    figures = verify_figures(path, design)
    verify_completion_seal(path, metadata, figures)
    return design, dict(code_commit=metadata['provenance']['code_commit'],
                        design_sha256=design.digest, audit_sha256=file_digest(path / 'audit.json'),
                        metadata_sha256=file_digest(path / 'run_metadata.json'))


def verify_completion_seal(path, metadata, figures):
    """Require a complete byte map and the actual bound visual-review record."""
    path = Path(path)
    actual = {str(p.relative_to(path)): file_digest(p) for p in sorted(path.rglob('*'))
              if p.is_file() and p != path / 'run_metadata.json' and p.name != '.execution.lock'}
    if metadata.get('artifact_sha256') != actual:
        raise ValueError('complete sealed artifact map changed')
    review_path = path / 'visual_review.json'
    if not review_path.is_file() or file_digest(review_path) != metadata.get('visual_review_sha256'):
        raise ValueError('visual review seal changed')
    review = json.loads(review_path.read_text())
    expected_pngs = {name for name in figures['figures'] if name.endswith('.png')}
    if (not expected_pngs or review.get('all_panels_legible') is not True
            or set(review.get('reviewed_pngs', [])) != expected_pngs
            or review.get('figure_manifest_sha256') != file_digest(path / 'figures/figure_manifest.json')
            or not isinstance(review.get('reviewer'), str) or not review['reviewer'].strip()):
        raise ValueError('visual review content does not bind the completed figures')


def saved_rows(path):
    for partition in sorted((Path(path) / 'cells').glob('*/*.jsonl.gz')):
        with gzip.open(partition, 'rt') as handle:
            for line in handle:
                yield json.loads(line)


def build_selection(e1, e2):
    result = {'scope': 'exploratory_selection_after_complete_headlines', 'inputs': {},
              'sources': [], 'failed_fit_records': [], 'unresolved_parameter_rows': []}
    failure_cells = {}
    for label, directory in [('E1', e1), ('E2', e2)]:
        design, identity = sealed_input(directory)
        registered = (load_e1(ROOT / 'configs/week5_e1.yaml') if label == 'E1'
                      else load_e2(ROOT / 'configs/week6_e2.yaml'))
        if design.to_dict() != registered.to_dict():
            raise ValueError('E3 selection requires the exact registered E1/E2 headline designs')
        result['inputs'][label] = identity
        sources, failures = select_records(saved_rows(directory))
        result['sources'].extend(dict(experiment=label, **s) for s in sources)
        result['failed_fit_records'].extend(dict(experiment=label, **f) for f in failures)
        cells = {c.identity: c for c in design.cells}
        for f in failures:
            failure_cells[f['cell_id']] = cells[f['cell_id']]
        summary = pd.read_csv(Path(directory) / 'parameter_summary.csv')
        for row in summary[summary.practical_classification == 'unresolved'].to_dict('records'):
            result['unresolved_parameter_rows'].append(dict(experiment=label, **{
                k: row[k] for k in ['cell_id', 'n', 'theta_mean_gap', 'cv', 'estimator', 'parameter']}))
    expected = {Cell(n=250, q=2, cv=c).identity for c in (0, .5)}
    if set(failure_cells) != expected:
        raise ValueError('published inputs differ from the dated two-cell E3 follow-up protocol')
    result['followup_cells'] = [c.to_dict() for c in sorted(failure_cells.values(), key=lambda c: c.cv)]
    return result


def followup_design():
    return RunDesign('week7_e3_followup', 20261007, 10000,
                     tuple(Cell(n=250, q=2, cv=c) for c in (0, .5)), 'exploratory')


def floor_design():
    # CV=0 at all three floors is an exact shared-source negative control.
    return RunDesign('week7_floor_sensitivity', 20261009, 2000,
                     tuple(Cell(n=1000, q=.5, cv=c, floor_fraction=epsilon)
                           for c in (0, 1.5, 2, 3) for epsilon in (1e-8, 1e-6, 1e-4)), 'exploratory')
