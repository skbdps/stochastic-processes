"""Read back the additive review from the repository root; print its evidence."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

import pandas as pd

from ou_irregular.week10.analytical_review import verify


ROOT = Path.cwd()
REVIEW = Path('artifacts/week10_analytical_review')
PACKAGE = Path('artifacts/week10_supplement')
DOCUMENTS = [
    'README.md', 'artifacts/README.md', 'manuscript/README.md',
    'manuscript/week10_source_map.md', 'manuscript/week10_draft3.md',
    'docs/week10_mechanism.md', 'docs/week10_pfml_boundary_moments.md',
    'docs/week10_analytical_review.md', 'docs/week10_reproduction.md',
    'docs/week10_supplement_completion.md', str(REVIEW/'README.md'),
]
SOURCES = [
    'ou_irregular/week10/analytical_review.py',
    'tests/week10/test_analytical_review.py',
    'tests/week10/test_boundary_argument.py',
    'ou_irregular/week10/theory.py',
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


before = json.loads((REVIEW/'before_review.json').read_text())
for name, expected in before['files'].items():
    assert digest(REVIEW/'before_review'/name) == expected, name
assert digest('manuscript/week10_draft2.md') == before['files']['manuscript/week10_draft2.md']
assert digest(PACKAGE/'completion_seal.json') == before['original_seal_sha256']
initial = json.loads((REVIEW/'initial_source/manifest.json').read_text())
for name, expected in initial['files'].items():
    assert digest(REVIEW/'initial_source'/name) == expected, name
assert digest(REVIEW/'calibration_initial/manifest.json') == initial['calibration_manifest_sha256']
initial_manifest = json.loads((REVIEW/'calibration_initial/manifest.json').read_text())
for name, expected in initial_manifest['files'].items():
    assert digest(REVIEW/'calibration_initial'/name) == expected, name

calibration = verify(PACKAGE, REVIEW/'calibration')
assert calibration == json.loads((REVIEW/'calibration_verify.json').read_text())
original = json.loads((REVIEW/'original_verify.stdout').read_text())
assert all(record['status'] == 'passed' for record in original.values())
assert original['seal']['seal_sha256'] == before['original_seal_sha256']
for name in ('original_verify.stderr', 'calibration_verify.stderr', 'week10_tests.stderr'):
    assert not (REVIEW/name).read_bytes(), name
suites = ET.parse(REVIEW/'week10_tests.xml').getroot().findall('testsuite')
counts = {key: sum(int(s.attrib[key]) for s in suites)
          for key in ('tests', 'failures', 'errors', 'skipped')}
assert counts == {'tests': 132, 'failures': 0, 'errors': 0, 'skipped': 0}, counts
numerical = json.loads((REVIEW/'numerical_validation/critical_numerical_audit.json').read_text())
assert numerical['max_relative_error'] < 5e-16
assert numerical['current_theory_sha256'] == digest('ou_irregular/week10/theory.py')

text = Path('docs/week10_analytical_review.md').read_text()
draft = Path('manuscript/week10_draft3.md').read_text()
numeric_checks = 0
for filename, design in [('original_sd_comparisons.csv', 'Original fixed N'),
                         ('matched_prediction_comparisons.csv', 'New fixed T')]:
    table = pd.read_csv(REVIEW/'calibration'/filename, float_precision='round_trip')
    for row in table[table.ratio.isin([5, 20])].itertuples(index=False):
        values = [row.predicted_sd_ratio, row.observed_pfml_sd_ratio, row.observed_pfml_rmse_ratio]
        assert f'| {design} | {int(row.ratio)} | ' + ' | '.join(f'{x:.6f}' for x in values) + ' |' in text
        assert f'| {int(row.ratio)} | {design} | ' + ' | '.join(f'{x:.4f}' for x in values) + ' |' in draft
        numeric_checks += 6
paired = pd.read_csv(REVIEW/'calibration/paired_coverage.csv', float_precision='round_trip')
cells = pd.read_csv(REVIEW/'calibration/cell_calibration.csv', float_precision='round_trip')
for mechanism in ('iid', 'mmpp'):
    row = paired[(paired.mechanism == mechanism) & (paired.ratio == 20) & (paired.contrast == 'C_minus_B')].iloc[0]
    cell = cells[(cells.mechanism == mechanism) & (cells.ratio == 20)].iloc[0]
    values = [f'{100*cell[label+"_coverage"]:.2f}%' for label in ('A', 'B', 'C')]
    interval = f'{100*row.coverage_difference:.3f} [{100*row.pointwise_mc_lower:.3f}, {100*row.pointwise_mc_upper:.3f}]'
    assert f'| {mechanism.upper()}, r=20 | ' + ' | '.join(values) + ' | ' + interval + ' |' in text
    assert f'{100*cell.exact_mixture_known_width_coverage_estimate:.4f}%' in text
    numeric_checks += 7
assert (paired[paired.contrast == 'B_minus_A'].coverage_difference.abs() <= .0005+1e-15).all()
zero = paired[(paired.first_only_covered == 0) & (paired.second_only_covered == 0)]
assert len(zero) == 2 and (zero.pointwise_mc_lower < 0).all() and (zero.pointwise_mc_upper > 0).all()

links = 0
for name in DOCUMENTS:
    content = Path(name).read_text()
    assert content.count('```') % 2 == 0, name
    assert content.count('$$') % 2 == 0, name
    for target in re.findall(r'\]\(([^)]+)\)', content):
        if target.startswith(('http:', 'https:', '#', 'mailto:')):
            continue
        target = target.split('#', 1)[0]
        resolved = (Path(name).parent/target).resolve()
        if resolved == (ROOT/REVIEW/'final_validation.json').resolve():
            continue  # This command prints the record written to that destination.
        assert resolved.exists(), (name, target)
        links += 1
subprocess.run(['git', 'diff', '--check'], check=True, capture_output=True, text=True)
files = {name: digest(name) for name in DOCUMENTS+SOURCES}
excluded = {'final_validation.json', 'final_validation.stderr'}
evidence = {str(path.relative_to(REVIEW)): digest(path)
            for path in sorted(REVIEW.rglob('*')) if path.is_file()
            and path.name not in excluded and '__pycache__' not in path.parts}
result = {
    'status': 'passed', 'verified_utc': datetime.now(timezone.utc).isoformat(),
    'scope': 'analytical correction and saved-data diagnostics; no new research replications',
    'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'committed_or_pushed': False, 'week10_tests': counts,
    'original_readback': original, 'calibration_readback': calibration,
    'draft2_unchanged': True, 'pre_review_copies_verified': len(before['files']),
    'initial_calibration_and_source_preserved': True,
    'numerical_document_entries_checked': numeric_checks,
    'local_document_links_checked': links, 'git_diff_whitespace_check': 'passed',
    'source_and_document_sha256': files, 'review_evidence_sha256': evidence,
    'limits': ['No academic review or publication decision is claimed.',
               'Fitted-mean empirical RMSE summaries do not establish ideal population moments.',
               'Original generated narrative is historical; revised interpretation is in Draft 3.'],
}
print(json.dumps(result, indent=2, allow_nan=False))
