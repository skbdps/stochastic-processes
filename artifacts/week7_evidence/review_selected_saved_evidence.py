"""Read-only independent arithmetic/disposition of closed selected diagnostics.

No OU project modules are imported. No observations, model fits, profile
evaluations, Hessian stencils, or random values are generated. This checks
saved source bindings, diagnostic decisions and algebra on saved arrays.
"""
from collections import Counter
from decimal import Decimal, getcontext
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh


getcontext().prec = 50
ROOT = Path('/tmp/ou_stochastic_recovery_20261030/artifacts/week7_selected')
OUTPUT = Path('/tmp/ou_resume_20261030/selected_saved_evidence_independent_review.json')
CHECKS = Counter()
ERRORS = []
MAX_ERROR = {}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot():
    return {str(path.relative_to(ROOT)): sha(path) for path in sorted(ROOT.rglob('*')) if path.is_file()}


def require(value, label):
    CHECKS[label] += 1
    if not value:
        ERRORS.append(label)


def near(actual, expected, label, *, rtol=2e-8, atol=5e-10):
    a, b = np.asarray(actual, float), np.asarray(expected, float)
    difference = float(np.max(np.abs(a-b)))
    MAX_ERROR[label] = max(MAX_ERROR.get(label, 0.), difference)
    require(bool(np.allclose(a, b, rtol=rtol, atol=atol)), label)


def number(record, field):
    return Decimal(record[field+'_extended'])


def check_ladder(ladder):
    entries = ladder['entries']
    require([e['relative_step'] for e in entries] == [.001, .0005, .00025, .000125], 'fixed_hessian_steps')
    for entry in entries:
        hessian = np.asarray(entry['hessian'], float)
        gradient = np.asarray(entry['gradient'], float)
        require(hessian.shape == (3, 3) and gradient.shape == (3,), 'saved_derivative_shapes')
        require(bool(np.array_equal(hessian, hessian.T)), 'matrix_symmetry')
        eigen = eigvalsh(hessian, driver='evr')
        near(eigen, entry['eigenvalues'], 'saved_eigenvalues')
        require(entry['positive_definite'] == bool(np.min(eigen) > 0), 'positive_definite_decision')
        near(np.linalg.svd(hessian, compute_uv=False)[0]/np.linalg.svd(hessian, compute_uv=False)[-1],
             entry['condition_number'], 'condition_number', atol=2e-6)
        near(hessian[0, 0], entry['theta_hessian_diagonal'], 'theta_diagonal')
        schur = hessian[0, 0]-hessian[0, 1:]@np.linalg.solve(hessian[1:, 1:], hessian[1:, 0])
        near(schur, entry['theta_schur_complement_at_stored_point'], 'schur_complement')
        if entry['positive_definite']:
            score = math.sqrt(gradient@np.linalg.solve(hessian, gradient))
            near(score, entry['standardized_score'], 'standardized_score')
        if ladder['precision'] == 'longdouble':
            near([[float(Decimal(v)) for v in row] for row in entry['hessian_extended']],
                 hessian, 'extended_matrix_display_conversion', rtol=0, atol=0)
    for coarse, fine, declared in zip(entries[:-1], entries[1:], ladder['adjacent']):
        a, b = np.asarray(coarse['hessian']), np.asarray(fine['hessian'])
        near(np.linalg.norm(a-b, 'fro')/np.linalg.norm(b, 'fro'),
             declared['relative_frobenius'], 'adjacent_frobenius')
        if fine['positive_definite']:
            # Generalized eigenvalues independently check the producer's
            # Cholesky-whitened spectral-radius calculation.
            spectral = float(np.max(np.abs(eigvalsh(a-b, b, driver='gvd'))))
            near(spectral, declared['whitened_spectral'], 'adjacent_generalized_eigenvalue')


before = snapshot()
metadata = json.loads((ROOT/'run_metadata.json').read_text())
audit = json.loads((ROOT/'audit.json').read_text())
selection = json.loads((ROOT/'selection.json').read_text())
require(metadata['status'] == 'completed_diagnostic_evidence', 'closed_diagnostics')
require(audit['status'] == 'passed' and audit['failures'] == [], 'saved_replay_audit_passed')
require(audit['independent_recalculation_of_all_diagnostic_math'] is False, 'original_audit_scope_retained')
require(metadata['audit_sha256'] == sha(ROOT/'audit.json'), 'root_audit_binding')
require(len(metadata['source_files']) == len(selection['sources']) == 53, 'full_53_source_selection')

profiles, boundary, key_cases, profile_rows = [], [], {}, []
counts = Counter()
for filename, digest in sorted(metadata['source_files'].items()):
    path = ROOT/'sources'/filename
    require(sha(path) == digest == audit['source_files'][filename], 'source_payload_hash')
    payload = json.loads(path.read_text())
    require(payload['status'] == 'evaluated' and payload['original_replay']['status'] == 'passed', 'evaluated_and_replayed')
    require(payload['headline_results_replaced'] is False, 'headline_records_preserved')
    diagnostic = payload['diagnostic']
    ar = diagnostic['ar1']
    original = payload['original_records']
    row = original['exact']
    counts['sources'] += 1
    phi = number(ar, 'phi')
    require(number(ar, 'innovation_variance') > 0, 'positive_regression_residual_variance')
    expected = ('PFML_supremum_at_phi_zero_no_finite_interior_estimate' if phi <= 0 else
                'PFML_supremum_at_phi_one_lower_theta_boundary' if phi >= 1 else
                'PFML_interior_positive_correlation')
    require(ar['conclusion'] == expected, 'ols_boundary_conclusion')
    require(ar['exact_boundary_proof_applicable'] == ar['exactly_equal_gaps'], 'boundary_proof_scope')
    old_phi = json.loads(original['pfml']['point_diagnostics_json'])['unrestricted_phi']
    near(float(phi), old_phi, 'independent_ols_vs_original_pfml_slope', rtol=2e-10, atol=2e-12)
    if phi <= 0:
        counts['nonpositive_ols_sources'] += 1
        require(original['pfml']['point_valid'] is False and original['pfml']['point_status'] == 'boundary',
                'pfml_boundary_retained')
        counts['pfml_boundary_records'] += 1
        if ar['exactly_equal_gaps']:
            require(original['exact']['point_valid'] is False and original['exact']['point_status'] == 'boundary',
                    'equal_gap_exact_boundary_retained')
            counts['exact_equal_gap_boundary_records'] += 1
        boundary.append({'source': filename, 'cv': row['cv'], 'rep': row['rep'],
                         'phi': float(phi), 'exactly_equal_gaps': ar['exactly_equal_gaps'],
                         'exact_point_valid': original['exact']['point_valid']})
    for method, fit in diagnostic['fits'].items():
        if original[method]['point_valid']:
            counts['valid_point_diagnostic_targets'] += 1
            check_ladder(fit['float64_ladder'])
            counts['saved_hessian_matrices'] += 4
            if method == 'exact':
                check_ladder(fit['independent_longdouble_ladder'])
                counts['saved_hessian_matrices'] += 4
        else:
            require(fit['status'] == 'skipped_original_point_invalid', 'invalid_candidate_not_reinterpreted')
            counts['invalid_point_diagnostic_skips'] += 1
    profile = diagnostic['exact_profile']
    if profile['status'] == 'evaluated':
        counts['exact_profiles'] += 1
        candidates = [c for c in profile['candidates'] if c['status'] == 'finite']
        best = min(number(c, 'nll') for c in candidates)
        require(number(profile['best_finite'], 'nll') == best, 'profile_best_saved_candidate')
        resolution = Decimal('1e-7')+Decimal('1e-9')*abs(best)
        require(abs(Decimal(profile['resolution_extended'])-resolution) < Decimal('1e-21'), 'profile_resolution_formula')
        expected_differences = {
            'stored_full_minus_best': number(profile['stored_full'], 'nll')-best,
            'stored_profile_minus_best': number(profile['stored_profile'], 'nll')-best,
            'stored_full_minus_stored_profile': number(profile['stored_full'], 'nll')-number(profile['stored_profile'], 'nll'),
            'brownian_limit_minus_best': number(profile['limits']['brownian'], 'nll')-best,
            'iid_limit_minus_best': number(profile['limits']['iid'], 'nll')-best,
        }
        for name, expected_value in expected_differences.items():
            # Serialized longdouble operands round at roughly 19 significant
            # digits. Allow their decimal conversion roundoff, not a change
            # in the much larger fixed scientific diagnostic resolution.
            bound = Decimal('2e-18')*max(Decimal(1), abs(best), abs(expected_value))
            require(abs(number(profile['differences'], name)-expected_value) <= bound,
                    'profile_extended_subtraction')
        bounds = profile['bounds']
        grid = np.asarray(profile['grid_log_theta'])
        require(grid[0] == bounds['lower'] and grid[-1] == bounds['upper'], 'profile_grid_endpoints')
        require(bool(np.all(np.diff(grid) > 0) and np.max(np.diff(grid)) <= .05), 'profile_grid_max_spacing')
        require([s['offset'] for s in profile['local_searches']] == [-2, -1, 0, 1, 2], 'five_fixed_local_starts')
        require(len(profile['grid_nll_extended']) == len(grid), 'complete_profile_grid_objectives')
        for search in profile['local_searches']:
            require(search['success'] is True, 'local_solver_reported_convergence')
            start = max(bounds['lower'], min(bounds['upper'], profile['stored_profile']['log_theta']+search['offset']))
            near(search['x0'][0], start, 'local_start_address', rtol=0, atol=1e-14)
        tail_unresolved = not bounds['lower_tail_reached'] or not bounds['upper_tail_reached']
        limiting = min(number(x, 'nll') for x in profile['limits'].values()) <= best+resolution
        improvement = expected_differences['stored_full_minus_best'] > resolution
        expected_disposition = ('cap_limited_tail_unresolved' if tail_unresolved else
                                'limiting_likelihood_competitive_or_numerically_unresolved' if limiting else
                                'finite_search_improves_stored_full_nll_above_diagnostic_resolution' if improvement else
                                'stored_full_fit_matches_finite_search_within_diagnostic_resolution')
        require(profile['disposition'] == expected_disposition, 'profile_disposition_rule')
        require(profile['global_optimum_proved'] is False and profile['replacement_fit_or_interval'] is False,
                'finite_search_scope_limit')
        profiles.append(profile['disposition'])
        local = [{'offset': s['offset'], 'success': s['success'],
                  'log_theta': s['candidate']['log_theta'],
                  'nll_minus_best': float(number(s['candidate'], 'nll')-best),
                  'above_resolution': number(s['candidate'], 'nll')-best > resolution}
                 for s in profile['local_searches']]
        counts['local_starts'] += len(local)
        counts['local_starts_above_best_resolution'] += sum(r['above_resolution'] for r in local)
        profile_rows.append({'source': filename, 'n': row['n'], 'q': row['theta_mean_gap'],
                             'cv': row['cv'], 'rep': row['rep'],
                             'stored_full_minus_best': profile['differences']['stored_full_minus_best'],
                             'resolution': float(resolution),
                             'iid_limit_minus_best': profile['differences']['iid_limit_minus_best'],
                             'brownian_limit_minus_best': profile['differences']['brownian_limit_minus_best'],
                             'both_tail_targets_reached': not tail_unresolved, 'local_starts': local})
        if payload['selection']['experiment'] == 'E2' and row['n'] == 250 and row['theta_mean_gap'] == 2 and (row['cv'], row['rep']) in [(0., 1424), (.5, 465)]:
            name = 'irregular_rejected_rep465' if row['cv'] else 'regular_accepted_rep1424'
            exact = diagnostic['fits']['exact']
            ladders = {}
            for precision in ['float64_ladder', 'independent_longdouble_ladder']:
                ladder = exact[precision]
                ladders[precision] = {
                    'adjacent_whitened': [a['whitened_spectral'] for a in ladder['adjacent']],
                    'adjacent_frobenius': [a['relative_frobenius'] for a in ladder['adjacent']],
                    'entries': [{k: e[k] for k in ['relative_step', 'eigenvalues', 'condition_number',
                                                 'standardized_score', 'theta_hessian_diagonal',
                                                 'theta_schur_complement_at_stored_point']}
                                for e in ladder['entries']]}
            ld = ladders['independent_longdouble_ladder']['adjacent_whitened']
            key_cases[name] = {
                'source': filename, 'original': {k: row[k] for k in ['theta', 'mu', 'sigma', 'point_valid',
                    'interval_valid', 'interval_status', 'hessian_relative_difference', 'hessian_whitened_difference']},
                'ar1': ar, 'precision': diagnostic['precision'], 'profile': profile_rows[-1],
                'representation_differences': diagnostic['representation_differences'], 'ladders': ladders,
                'longdouble_adjacent_decay_ratios': [ld[i]/ld[i+1] for i in range(2)],
                'interpretation': 'Weak curvature amplifies finite-difference truncation: independent extended arithmetic retains the initial discrepancy and halving h gives approximately fourfold decay; cancellation alone does not explain the original gate outcome. No replacement interval is constructed.'}
    else:
        require(profile['status'] == 'skipped_original_point_invalid', 'invalid_exact_profile_skipped')

require(snapshot() == before, 'complete_archive_byte_preservation')
result = {
    'status': 'passed' if not ERRORS else 'review_required',
    'scope': 'independent algebra/decision review of saved selected diagnostic payloads; no samples, fits, objectives or stencils recalculated',
    'input_archive': str(ROOT), 'input_sha256': before,
    'recorded_code_commit': metadata['provenance']['code_commit'],
    'original_replay_audit': {k: v for k, v in audit.items() if k != 'source_files'},
    'checks': dict(CHECKS), 'total_checks': sum(CHECKS.values()), 'failures': ERRORS,
    'matrix_review_comparison_tolerance': {'relative': 2e-8, 'absolute': 5e-10,
        'condition_absolute': 2e-6, 'purpose': 'algebra-recalculation roundoff only; does not change any registered fit/inference gate'},
    'matrix_review_arithmetic': 'Binary64 NumPy/SciPy linear algebra on saved numerical matrix entries, including the binary64 displays of longdouble matrices. Extended decimal entries are verified to round to those displays; spectral algebra is not independently repeated at longdouble precision.',
    'maximum_absolute_recalculation_differences': MAX_ERROR,
    'counts': dict(counts), 'profile_dispositions': dict(Counter(profiles)),
    'nonpositive_ols_cases': boundary, 'profile_cases': profile_rows, 'key_cases': key_cases,
    'limitations': [
        'Observed arrays are not stored in these selected diagnostic payloads, so this review does not independently recompute OLS slopes or likelihoods from raw observations.',
        'OLS decisions are checked against the saved independently calculated slope and archived unrestricted PFML coefficient; the equal-gap boundary implication is analytical.',
        'Saved full-Hessian matrices are algebraically checked; finite-difference objective evaluations are not repeated.',
        'Matrix spectral checks use binary64 arithmetic on saved numerical entries, not an independent longdouble or arbitrary-precision reconstruction of the Hessians.',
        'All 14 finite searches match the stored exact fit within their fixed resolution; finite searches do not prove global optima.',
        'Boundary and interval-unavailable headline records remain unchanged; the diagnostics are selected exploratory evidence and do not estimate population rates.',
        'No newly verified registered-estimator/inference implementation defect is established by these selected outcomes. Historical unattributed runtime incidents remain a separate operational limitation.'
    ],
}
OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
print(json.dumps({'status': result['status'], 'counts': result['counts'], 'checks': result['total_checks'],
                  'failures': ERRORS, 'profile_dispositions': result['profile_dispositions'],
                  'review_path': str(OUTPUT)}, sort_keys=True))
