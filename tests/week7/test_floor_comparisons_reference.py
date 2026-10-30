"""Literal paired-floor references: no observations or random draws generated.

The records are deliberately fabricated arithmetic fixtures with canonical
scientific/address metadata. Independent rational moments and covariance
expansions test the paired formulas instead of reusing production helpers.
"""
import hashlib
import itertools
import json
import math
from statistics import NormalDist

import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments.design import Cell, ESTIMATORS, RunDesign
from ou_irregular.week7.floor_comparisons import epsilon_comparisons


EPSILONS = (1e-8, 1e-6, 1e-4)
PARAMETERS = ('theta', 'mu', 'sigma')
Z = NormalDist().inv_cdf(.975)
POLICY = 'week7_paired_floor_v1_shared_gamma_stationary_x0_and_standard_normal_innovations'


def close(value):
    return pytest.approx(value, rel=1e-12, abs=2e-14)


def fabricated_records(errors_a, errors_b, *, point_a=None, point_b=None,
                       interval_a=None, interval_b=None, cover_a=None, cover_b=None,
                       cvs=(1.5,)):
    """Complete three-floor fixtures; CV0 always copies the reference evidence."""
    count = len(errors_a)
    assert len(errors_b) == count
    point_a = list(point_a if point_a is not None else [True]*count)
    point_b = list(point_b if point_b is not None else [True]*count)
    interval_a = list(interval_a if interval_a is not None else point_a)
    interval_b = list(interval_b if interval_b is not None else point_b)
    cover_a = list(cover_a if cover_a is not None else [True]*count)
    cover_b = list(cover_b if cover_b is not None else [True]*count)
    design = RunDesign('paired_floor_literal_reference', 21261009, count,
                       tuple(Cell(1000, .5, cv, floor_fraction=epsilon)
                             for cv in cvs for epsilon in EPSILONS), 'validation')
    rows = []
    for cell in design.cells:
        is_a = cell.floor_fraction == 1e-8 and cell.cv != 0
        errors = errors_a if is_a else errors_b
        point = point_a if is_a else point_b
        interval = interval_a if is_a else interval_b
        hits = cover_a if is_a else cover_b
        anchor = Cell(1000, .5, cell.cv, floor_fraction=1e-6)
        words = tuple(int(anchor.identity[i:i+8], 16) for i in range(0, 64, 8))
        for rep in range(count):
            key = (4, *words, rep)
            encode = lambda x: json.dumps(x, separators=(',', ':'))
            # These strings identify fabricated evidence; they are not hashes
            # of generated observations. Irregular variants may differ, while
            # the regular negative control must have identical evidence.
            hash_id = (cell.cv, 1e-6 if cell.cv == 0 else cell.floor_fraction, rep)
            source = {'actual_cell_id': cell.identity, 'seed_anchor_cell_id': anchor.identity,
                      'seed_anchor_floor_fraction': 1e-6, 'source_coupling_policy': POLICY,
                      'source_seed': design.seed, 'source_spawn_key': encode(key),
                      'gap_spawn_key': encode((*key, 0)), 'path_spawn_key': encode((*key, 1)),
                      'times_sha256': hashlib.sha256(repr(('literal-times', hash_id)).encode()).hexdigest(),
                      'observations_sha256': hashlib.sha256(repr(('literal-values', hash_id)).encode()).hexdigest(),
                      'source_valid': True, 'source_error': '', 'pipeline_error': ''}
            for estimator in ESTIMATORS:
                row = {**cell.to_dict(), **source, 'rep': rep, 'estimator': estimator,
                       'point_valid': point[rep], 'optimizer_success': point[rep],
                       'point_status': 'valid_interior' if point[rep] else 'boundary_supremum',
                       'nll': 0., 'interval_valid': interval[rep],
                       'interval_status': ('valid' if interval[rep] else
                                           'step_disagreement' if point[rep] else 'invalid_point_fit')}
                for parameter, truth in zip(PARAMETERS, (1., 0., .5)):
                    estimate = truth+errors[rep]
                    if interval[rep]:
                        width = abs(errors[rep])+.125 if hits[rep] else abs(errors[rep])/2
                        assert width > 0
                        lower, upper, se = estimate-width, estimate+width, width/Z
                    else:
                        # Finite endpoints deliberately cover the truth even
                        # though interval_valid=False: availability must enter
                        # the operational coverage indicator explicitly.
                        lower, upper, se = truth-1, truth+1, np.nan
                    row.update({parameter: estimate, f'{parameter}_lower': lower,
                                f'{parameter}_upper': upper, f'{parameter}_se': se})
                rows.append(row)
    return design, pd.DataFrame(rows)


def selected(result, *, cv=1.5, epsilon=1e-8, estimator='exact', parameter='theta'):
    rows = result[(result.cv == cv) & (result.floor_fraction == epsilon)
                  & (result.estimator == estimator) & (result.parameter == parameter)]
    assert len(rows) == 1
    return rows.iloc[0]


@pytest.fixture(scope='module')
def rational_result():
    design, records = fabricated_records(
        [1., 2., 4., 100., 999.], [2., 1., 3., 999., 20.],
        point_a=[True, True, True, True, False], point_b=[True, True, True, False, True],
        interval_a=[True, True, True, False, False], interval_b=[True, False, True, False, True],
        cover_a=[True, True, False, False, False], cover_b=[True, False, True, False, True])
    # Reversal gives deterministic out-of-order input without constructing RNG.
    return epsilon_comparisons(records.iloc[::-1].reset_index(drop=True), design)


@pytest.mark.parametrize('estimator,parameter', tuple(itertools.product(ESTIMATORS, PARAMETERS)))
def test_rational_common_valid_point_contrasts_and_full_delta_covariance(
        rational_result, estimator, parameter):
    row = selected(rational_result, estimator=estimator, parameter=parameter)
    assert row.n_planned == 5
    assert row.n_point_valid_a == row.n_point_valid_b == 4
    assert row.n_common_point_valid == 3
    assert row.n_interval_valid_a == row.n_interval_valid_b == 3
    assert row.n_common_interval_valid == 2
    assert row.bias_a_common == close(7/3)
    assert row.bias_b_common == close(2.)
    assert row.bias_a_common != close(107/4)  # A's own-valid mean.
    assert row.bias_b_common != close(26/4)   # B's own-valid mean.
    assert row.mse_a_common == close(7.)
    assert row.mse_b_common == close(14/3)
    assert row.rmse_a_common == close(math.sqrt(7))
    assert row.rmse_b_common == close(math.sqrt(14/3))
    assert row.bias_difference == close(1/3)
    assert row.bias_difference_mcse == close(2/3)
    assert row.mse_difference == close(7/3)
    assert row.mse_difference_mcse == close(math.sqrt(76)/3)
    assert row.squared_error_covariance_common == close(27.)
    # Sample Var(A²)=63, Var(B²)=49/3, Cov(A²,B²)=27. Expanding the
    # gradient (1/(2*RMSE_A),-1/(2*RMSE_B)) retains the paired covariance.
    expected_variance = (63/(4*7)+(49/3)/(4*(14/3))
                         -27/(2*math.sqrt(7)*math.sqrt(14/3)))/3
    assert row.rmse_difference == close(math.sqrt(7)-math.sqrt(14/3))
    assert row.rmse_difference_mcse**2 == close(expected_variance)
    assert row.rmse_ratio == close(math.sqrt(1.5))
    assert row.rmse_ratio_log_mcse == close(5/28)
    assert row.rmse_ratio_mcse == close(math.sqrt(1.5)*5/28)
    assert row.bias_difference_mc_lower == close(1/3-Z*2/3)
    assert row.bias_difference_mc_upper == close(1/3+Z*2/3)
    assert row.rmse_ratio_mc_lower == close(math.sqrt(1.5)*math.exp(-Z*5/28))
    assert row.rmse_ratio_mc_upper == close(math.sqrt(1.5)*math.exp(Z*5/28))
    assert row.point_subset == 'intersection_of_valid_point_fits'
    assert row.difference_direction == 'epsilon_minus_registered_floor'
    assert row.ratio_direction == 'epsilon_over_registered_floor'


@pytest.mark.parametrize('parameter', PARAMETERS)
def test_rational_operational_coverage_retains_all_planned_pairs(rational_result, parameter):
    row = selected(rational_result, parameter=parameter)
    assert (row.n_both_cover, row.n_a_only_cover, row.n_b_only_cover, row.n_neither_cover) == (1, 1, 2, 1)
    assert row.n_covered_a == 2 and row.n_covered_b == 3
    assert row.operational_coverage_a == .4 and row.operational_coverage_b == .6
    assert row.operational_coverage_difference == close(-.2)
    assert row.operational_coverage_difference_mcse == close(math.sqrt(.14))
    assert row.operational_coverage_difference_mc_lower == close(-.2-Z*math.sqrt(.14))
    assert row.operational_coverage_difference_mc_upper == close(-.2+Z*math.sqrt(.14))
    assert row.coverage_subset == 'all_planned_paths_unavailable_counts_as_not_covering'


def test_zero_truth_relative_metrics_remain_undefined_but_absolute_mu_ratio_exists(rational_result):
    mu = selected(rational_result, parameter='mu')
    relative_columns = [name for name in rational_result if name.startswith('relative_')]
    assert relative_columns
    assert mu[relative_columns].isna().all()
    assert mu.bias_difference == close(1/3)
    assert mu.rmse_ratio == close(math.sqrt(1.5))
    sigma = selected(rational_result, parameter='sigma')
    assert sigma.relative_bias_difference == close(2/3)
    assert sigma.relative_bias_difference_mcse == close(4/3)


def test_proportional_errors_cancel_log_ratio_variance_through_covariance():
    design, records = fabricated_records([2., 4., 6.], [1., 2., 3.])
    row = selected(epsilon_comparisons(records, design))
    assert row.squared_error_covariance_common > 0
    assert row.rmse_ratio == 2.
    assert row.rmse_ratio_log_mcse == row.rmse_ratio_mcse == 0.
    assert row.rmse_ratio_mc_lower == row.rmse_ratio_mc_upper == 2.
    assert row.bias_difference_mcse == close(1/math.sqrt(3))
    # Observed proportionality cancels this sample's ratio influence; it does
    # not assert a deterministic population ratio for irregular floor variants.


def test_regular_control_is_identical_with_zero_contrast_uncertainty_when_defined():
    design, records = fabricated_records([1., 2., 3.], [1., 2., 3.],
                                         cover_a=[True, False, True], cover_b=[True, False, True],
                                         cvs=(0.,))
    result = epsilon_comparisons(records, design)
    assert len(result) == 18
    for metric in ('bias_difference', 'mse_difference', 'rmse_difference', 'operational_coverage_difference'):
        assert result[metric].eq(0.).all()
        assert result[metric+'_mcse'].eq(0.).all()
        assert result[metric+'_mc_lower'].eq(0.).all()
        assert result[metric+'_mc_upper'].eq(0.).all()
    assert result.rmse_ratio.eq(1.).all()
    assert result.rmse_ratio_log_mcse.eq(0.).all()
    assert result.rmse_ratio_mcse.eq(0.).all()
    assert result.rmse_ratio_mc_lower.eq(1.).all() and result.rmse_ratio_mc_upper.eq(1.).all()


@pytest.mark.parametrize('pattern', ['all_invalid', 'disjoint', 'one_common', 'one_planned',
                                     'zero_numerator', 'zero_denominator', 'both_zero'])
def test_missing_and_singular_metric_cases_keep_their_mathematical_domains(pattern):
    kwargs = {}
    errors_a, errors_b = [1., 2.], [2., 1.]
    if pattern == 'all_invalid':
        kwargs = {'point_a': [False, False], 'point_b': [False, False]}
    elif pattern == 'disjoint':
        kwargs = {'point_a': [True, False], 'point_b': [False, True]}
    elif pattern == 'one_common':
        errors_a, errors_b = [2., 99.], [1., 99.]
        kwargs = {'point_a': [True, False], 'point_b': [True, False]}
    elif pattern == 'one_planned':
        errors_a, errors_b = [2.], [1.]
    elif pattern == 'zero_numerator':
        errors_a, errors_b = [0., 0.], [1., 2.]
    elif pattern == 'zero_denominator':
        errors_a, errors_b = [1., 2.], [0., 0.]
    elif pattern == 'both_zero':
        errors_a = errors_b = [0., 0.]
    design, records = fabricated_records(errors_a, errors_b, **kwargs)
    row = selected(epsilon_comparisons(records, design))
    assert np.isnan(row.rmse_ratio_log_mcse)
    assert np.isnan(row.rmse_ratio_mcse)
    assert np.isnan(row.rmse_ratio_mc_lower) and np.isnan(row.rmse_ratio_mc_upper)
    if pattern in ('all_invalid', 'disjoint'):
        assert row.n_common_point_valid == 0
        for name in ('bias_a_common', 'bias_b_common', 'mse_a_common', 'mse_b_common',
                     'bias_difference', 'mse_difference', 'rmse_difference', 'rmse_ratio'):
            assert np.isnan(row[name])
    if pattern in ('one_common', 'one_planned'):
        assert row.n_common_point_valid == 1 and row.rmse_ratio == 2.
        for name in ('bias_difference_mcse', 'mse_difference_mcse', 'rmse_difference_mcse',
                     'squared_error_covariance_common'):
            assert np.isnan(row[name])
    if pattern == 'one_planned':
        assert np.isnan(row.operational_coverage_difference_mcse)
    if pattern == 'one_common':
        assert row.operational_coverage_difference_mcse == 0.
    if pattern == 'zero_numerator':
        assert row.rmse_ratio == 0.
        assert row.rmse_difference == close(-math.sqrt(2.5))
        assert np.isnan(row.rmse_difference_mcse)
    if pattern in ('zero_denominator', 'both_zero'):
        assert np.isnan(row.rmse_ratio) and np.isnan(row.rmse_difference_mcse)
    if pattern == 'both_zero':
        assert row.bias_difference == row.mse_difference == row.rmse_difference == 0.
        assert row.bias_difference_mcse == row.mse_difference_mcse == 0.
    if pattern == 'all_invalid':
        assert row.operational_coverage_difference == row.operational_coverage_difference_mcse == 0.
    if pattern == 'disjoint':
        assert row.operational_coverage_difference == 0.
        assert row.operational_coverage_difference_mcse == 1.
        assert row.operational_coverage_difference_mc_lower == close(-Z)
        assert row.operational_coverage_difference_mc_upper == close(Z)


def test_regular_singleton_keeps_sample_uncertainty_undefined():
    design, records = fabricated_records([1.], [1.], cvs=(0.,))
    result = epsilon_comparisons(records, design)
    assert result.rmse_ratio.eq(1.).all()
    for name in ('bias_difference_mcse', 'mse_difference_mcse', 'rmse_difference_mcse',
                 'rmse_ratio_log_mcse', 'rmse_ratio_mcse', 'operational_coverage_difference_mcse'):
        assert result[name].isna().all()


def test_complete_frozen_grid_produces_exactly_72_directional_rows():
    design, records = fabricated_records([2., 4., 6.], [1., 2., 3.], cvs=(0., 1.5, 2., 3.))
    result = epsilon_comparisons(records.iloc[::-1].reset_index(drop=True), design)
    assert len(result) == 72
    expected = {(cell.identity, method, parameter)
                for cell in design.cells if cell.floor_fraction != 1e-6
                for method in ESTIMATORS for parameter in PARAMETERS}
    assert set(result[['cell_id', 'estimator', 'parameter']].itertuples(index=False, name=None)) == expected
    assert result.groupby('cv').size().to_dict() == {0.: 18, 1.5: 18, 2.: 18, 3.: 18}
    assert not result.floor_fraction.eq(1e-6).any()
    assert result.reference_floor_fraction.eq(1e-6).all()
    for cv in (1.5, 2., 3.):
        assert selected(result, cv=cv).rmse_ratio == 2.
    assert selected(result, cv=0.).rmse_ratio == 1.
