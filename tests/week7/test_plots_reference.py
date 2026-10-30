"""Scientific display references; no sources, estimators or audits are rerun.

One real floor input is copied from the already published R=2 validation
archive and bound to its byte-preserved, previously computed full audit.  The
other figures use explicitly literal table/payload inputs, never purported
original E2 replay evidence.  Literal fixtures expose missing domains, paired
denominators, extended-decimal subtraction and saved uncertainty endpoints.
"""
from contextlib import contextmanager
from decimal import Decimal
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
from statistics import NormalDist

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments import audit as historical_audit
from ou_irregular.experiments import execution as historical_execution
from ou_irregular.experiments.analysis import numerical_diagnostics
from ou_irregular.experiments.comparisons import (
    paired_estimator_comparisons, regular_reference_comparisons,
)
from ou_irregular.experiments.design import Cell, ESTIMATORS
from ou_irregular.experiments.provenance import ROOT
from ou_irregular.experiments.summaries import (
    ExperimentSummarizer, outcome_counts, spacing_summary,
)
from ou_irregular.inference.adapters import Week3Estimator
from ou_irregular.week4.simulation import OuSampleSource
from ou_irregular.week7 import execution, floor_audit, path_diagnostics
from ou_irregular.week7.floor_comparisons import epsilon_comparisons
from ou_irregular.week7.floor_source import PairedFloorSource
from ou_irregular.week7.selection import floor_design, followup_design

from .test_floor_comparisons_reference import fabricated_records
from .test_reporting_reference import _augment


METHODS = ('exact', 'pfml', 'euler')
PARAMETERS = ('theta', 'mu', 'sigma')
FLOORS = (1e-8, 1e-6, 1e-4)
STEPS = (.001, .0005, .00025, .000125)
Z = NormalDist().inv_cdf(.975)
_SAVED_FLOOR_ARCHIVE = None
_RENDERED_FLOOR_ARCHIVE = None


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot(directory):
    return {str(path.relative_to(directory)): sha(path)
            for path in sorted(directory.rglob('*')) if path.is_file()}


def forbidden(*args, **kwargs):
    raise AssertionError('saved-data plotting must not generate sources, refit, or rerun an audit')


@contextmanager
def forbid_science():
    """Instrument every relevant public generation/replay entry point."""
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(OuSampleSource, 'sample', forbidden)
        patch.setattr(PairedFloorSource, 'sample', forbidden)
        patch.setattr(Week3Estimator, 'fit', forbidden)
        patch.setattr(historical_execution, 'source_records', forbidden)
        patch.setattr(execution, 'run_raw', forbidden)
        patch.setattr(historical_audit, 'audit_run', forbidden)
        patch.setattr(floor_audit, 'audit_run', forbidden)
        patch.setattr(floor_audit, 'independent_source', forbidden)
        patch.setattr(path_diagnostics, 'run_selected', forbidden)
        patch.setattr(path_diagnostics, 'analyze_sample', forbidden)
        yield


def saved_floor_archive(tmp_path_factory):
    """One copied, genuinely audited input shared through pytest's week7 module.

    The original raw-complete lifecycle remains unchanged: the actual runner
    saves its full numerical audit before creating figures, then advances the
    root lifecycle. Raw, table, dependency-provenance and audit bytes retain
    their prior identities.
    This function neither executes nor fabricates a scientific audit.
    """
    global _SAVED_FLOOR_ARCHIVE
    if _SAVED_FLOOR_ARCHIVE is None:
        source = ROOT/'artifacts/week7_recovery_evidence/floor_audit_current_capture/archive'
        report = ROOT/'artifacts/week7_recovery_evidence/plots_reference_inputs/floor_audit.json'
        before = snapshot(source)
        assert sha(report) == '6acdc696a073886d73c6650831c588f20916b7f950f0c1262c0935a57d926e2c'
        output = tmp_path_factory.mktemp('plots_saved_floor_reference')/'archive'
        shutil.copytree(source, output)
        shutil.copyfile(report, output/'audit.json')
        metadata = json.loads((output/'run_metadata.json').read_text())
        assert metadata['status'] == 'raw_complete_pending_audit'
        assert {key: value for key, value in snapshot(output).items()
                if key != 'audit.json'} == before
        assert snapshot(source) == before
        _SAVED_FLOOR_ARCHIVE = output, floor_design().validation(2)
    return _SAVED_FLOOR_ARCHIVE


def rendered_floor_archive(tmp_path_factory):
    """One real four-figure PNG/PDF set; all science/replay calls forbidden."""
    global _RENDERED_FLOOR_ARCHIVE
    if _RENDERED_FLOOR_ARCHIVE is None:
        from ou_irregular.week7 import plots
        source, design = saved_floor_archive(tmp_path_factory)
        before = snapshot(source)
        output = tmp_path_factory.mktemp('plots_rendered_floor_reference')/'archive'
        shutil.copytree(source, output)
        tables = {name: pd.read_csv(output/name, float_precision='round_trip')
                  for name in ('parameter_summary.csv', 'spacing_summary.csv')}
        epsilon = pd.read_csv(output/'epsilon_comparisons.csv', float_precision='round_trip')
        visited = []
        save = plots.figure_bytes
        def checked_save(figure, stem):
            # Observe the actual figure about to be serialized by the public
            # writer, so the integrated loader cannot silently route a
            # different table while only pure-builder tests remain correct.
            assert_floor_figure_matches_saved(figure, stem, design, tables, epsilon)
            visited.append(stem)
            return save(figure, stem)
        with forbid_science(), pytest.MonkeyPatch.context() as patch:
            patch.setattr(plots, 'figure_bytes', checked_save)
            plots.write_figures('floor', output, design)
        assert visited == ['W7_floor_paired_bias', 'W7_floor_paired_coverage',
                           'W7_floor_rmse_ratio', 'W7_floor_spacing']
        for key, digest in before.items():
            assert sha(output/key) == digest, key
        assert snapshot(source) == before
        _RENDERED_FLOOR_ARCHIVE = output, design
    return _RENDERED_FLOOR_ARCHIVE


def _tables(records, design):
    data = pd.DataFrame(records)
    return {
        'parameter_summary.csv': ExperimentSummarizer().summarize(data, design),
        'regular_comparisons.csv': regular_reference_comparisons(data, design),
        'paired_comparisons.csv': paired_estimator_comparisons(data, design),
        'spacing_summary.csv': spacing_summary(data, design),
        'outcome_counts.csv': outcome_counts(data, design),
        'numerical_diagnostics.csv': numerical_diagnostics(data, design),
        'failed_fits.csv.gz': data.loc[~data.point_valid | ~data.interval_valid].copy(),
    }


def literal_followup_tables():
    """Complete fixed cell grid with deliberate, parameter-shared failures."""
    design = followup_design().validation(5)
    counts = {0.: ((5, 4), (4, 3), (3, 2)), .5: ((4, 2), (3, 1), (0, 0))}
    records = []
    for cell in design.cells:
        for method, (points, intervals) in zip(METHODS, counts[cell.cv]):
            for rep in range(design.replications):
                point, interval = rep < points, rep < intervals
                row = {**cell.to_dict(), 'rep': rep, 'estimator': method,
                       'point_valid': point, 'optimizer_success': point,
                       'point_status': 'valid_interior' if point else 'boundary_supremum',
                       'interval_valid': interval,
                       'interval_status': 'valid' if interval else
                                          'step_disagreement' if point else 'invalid_point_fit',
                       'nll': float(rep), 'pipeline_error': ''}
                for parameter, truth in zip(PARAMETERS, (1., 0., .5)):
                    error = float(rep+1)
                    # The first available interval misses theta, the last
                    # misses sigma; all available mu intervals cover.
                    hit = (parameter == 'mu' or
                           parameter == 'theta' and rep > 0 or
                           parameter == 'sigma' and rep < intervals-1)
                    width = error+.25 if hit else error/2
                    row.update({parameter: truth+error,
                                parameter+'_se': width/Z if interval else np.nan,
                                parameter+'_lower': truth+error-width,
                                parameter+'_upper': truth+error+width})
                records.append(_augment(row))
    return design, _tables(records, design)


def literal_floor_tables():
    """Complete 12-cell, 72-comparison literal grid, including S=0/1 and RMSE=0."""
    rational = dict(
        errors_a=[1., 2., 4., 100., 999.], errors_b=[2., 1., 3., 999., 20.],
        point_a=[True, True, True, True, False], point_b=[True, True, True, False, True],
        interval_a=[True, True, True, False, False], interval_b=[True, False, True, False, True],
        cover_a=[True, True, False, False, False], cover_b=[True, False, True, False, True])
    groups = []
    for cv in (0., 1.5):
        _, rows = fabricated_records(**rational, cvs=(cv,))
        groups.append(rows)
    _, zero = fabricated_records([0.]*5, [1., 2., 3., 4., 5.],
                                 interval_a=[True, True, True, True, False],
                                 cover_a=[True, True, True, True, False],
                                 cover_b=[False, False, False, False, True], cvs=(2.,))
    groups.append(zero)
    _, unavailable = fabricated_records([1.]*5, [1., 2., 3., 4., 5.],
                                        point_a=[False]*5, cvs=(3.,))
    mask = ((unavailable.floor_fraction == 1e-4)
            & ((unavailable.estimator == 'pfml') & (unavailable.rep > 0)
               | (unavailable.estimator == 'euler')))
    unavailable.loc[mask, ['point_valid', 'optimizer_success', 'interval_valid']] = False
    unavailable.loc[mask, 'point_status'] = 'boundary_supremum'
    unavailable.loc[mask, 'interval_status'] = 'invalid_point_fit'
    for parameter in PARAMETERS:
        unavailable.loc[mask, parameter+'_se'] = np.nan
    groups.append(unavailable)
    records = []
    for source_row in pd.concat(groups, ignore_index=True).to_dict('records'):
        row = _augment(source_row)
        rep, cv, epsilon = row['rep'], row['cv'], row['floor_fraction']
        if cv:
            row['realized_mean_gap'] = .5*(1+epsilon*(rep+1))
            row['realized_span'] = 999*row['realized_mean_gap']
            row['realized_cv'] = cv+.01*rep
            row['clip_count'] = rep+1
            row['clip_fraction'] = (rep+1)/999
            row['min_gap'] = .5*epsilon*(rep+1)
        records.append(row)
    design = floor_design().validation(5)
    epsilon = epsilon_comparisons(records, design)
    return design, _tables(records, design), epsilon


def _literal_ladder(*, extended=False, indefinite=False):
    entries = []
    for index, step in enumerate(STEPS):
        first = (.01, .1, 0., -.2)[index] if indefinite else .01
        bad = first <= 0
        eigenvalues = [first, 2., 10.]
        score = None if bad else (0. if index % 2 == 0 else .002)
        entries.append({
            'relative_step': step, 'steps': [step]*3,
            'status': 'evaluated', 'hessian': np.diag(eigenvalues).tolist(),
            'gradient': [0. if score is None else score*math.sqrt(eigenvalues[0]), 0., 0.],
            'eigenvalues': eigenvalues,
            'positive_definite': not bad,
            'condition_number': None if bad else 10./first,
            'standardized_score': score,
            'theta_hessian_diagonal': eigenvalues[0],
            'theta_schur_complement_at_stored_point': eigenvalues[0],
            'schur_status': 'evaluated',
            'schur_profile_interpretation': 'only at nuisance stationarity in the same representation',
        })
    adjacent = []
    for index, (coarse, fine) in enumerate(zip(STEPS[:-1], STEPS[1:])):
        a, b = entries[index]['eigenvalues'][0], entries[index+1]['eigenvalues'][0]
        adjacent.append({
            'coarse_relative_step': coarse, 'fine_relative_step': fine,
            'relative_frobenius': abs(a-b)/math.sqrt(b*b+104.),
            'whitened_spectral': abs(a-b)/b if b > 0 else None,
            'status': 'fine_hessian_not_positive_definite' if b <= 0 else 'evaluated',
        })
    return {'status': 'evaluated', 'precision': 'longdouble' if extended else 'registered_float64',
            'entries': entries, 'adjacent': adjacent, 'replacement_interval_constructed': False}


def literal_selected_payload(cv):
    """Direct builder input only: deliberately no original-record or replay claim."""
    cell = Cell(250, 2., cv)
    rep = 465 if cv == .5 else 1424
    best = '1000000.000000000000'
    extended_grid = ['1000000.000000000020', '1000000.000000000040',
                     '1000000.000000000010', best,
                     '999999.999999999990', '1000000.000000000030']
    scalar = lambda value: {'nll': float(value), 'nll_extended': value}
    profile = {
        'status': 'evaluated', 'grid_log_theta': [-1., 0., .25, .5, 1., 2.],
        'grid_nll': list(map(float, extended_grid)), 'grid_nll_extended': extended_grid,
        'best_finite': {'log_theta': 1., **scalar(best)},
        'stored_profile': {'log_theta': .25, **scalar('1000000.000000000010')},
        'stored_full': scalar('1000000.000000000050'),
        'limits': {'brownian': scalar('1000000.000000000060'),
                   'iid': scalar('999999.999999999980')},
        'resolution': .0010001, 'resolution_extended': '.0010001',
        'bounds': {'lower': -1., 'upper': 2., 'lower_tail_reached': False,
                   'upper_tail_reached': False, 'lower_cap_limited': True,
                   'upper_cap_limited': True},
        'disposition': 'cap_limited_tail_unresolved', 'global_optimum_proved': False,
        'verified_software_defect': False, 'replacement_fit_or_interval': False,
    }
    deltas = {'stored_profile_minus_best': '0.000000000010',
              'stored_full_minus_best': '0.000000000050',
              'stored_full_minus_stored_profile': '0.000000000040',
              'brownian_limit_minus_best': '0.000000000060',
              'iid_limit_minus_best': '-0.000000000020'}
    profile['differences'] = {key: float(value) for key, value in deltas.items()}
    profile['differences'].update({key+'_extended': value for key, value in deltas.items()})
    fits = {}
    for method in METHODS:
        rejected = cv == .5 and method == 'exact'
        record = {
            'status': 'evaluated', 'original_point_valid': True,
            'original_point_status': 'valid_interior',
            'original_interval_valid': not rejected,
            'original_interval_status': 'step_disagreement' if rejected else 'valid',
            'working_point': [.25, 0., 0.],
            'float64_ladder': _literal_ladder(indefinite=rejected),
        }
        if method == 'exact':
            record['independent_longdouble_ladder'] = _literal_ladder(extended=True)
        fits[method] = record
    return {
        'selection': {'experiment': 'E2', 'cell_id': cell.identity, 'rep': rep,
                      'reasons': ['literal validation fixture, not selected observations']},
        'status': 'evaluated', 'headline_results_replaced': False,
        'literal_validation_fixture': True,
        'diagnostic': {
            'status': 'evaluated', 'fits': fits, 'exact_profile': profile,
            'coordinates': {'state_offset': .2, 'state_scale': 2., 'time_scale': 4.,
                            'source': 'observed_data_only'},
            'precision': {'dtype': 'float128', 'storage_bytes': 16, 'fraction_bits': 63,
                          'significand_bits': 64, 'decimal_precision': 18,
                          'epsilon': '1.084202172485504434e-19',
                          'extended_beyond_float64': True, 'arbitrary_precision': False},
            'representation_differences': {
                'registered_ladder': 'lift_exact_float64_normalized_values_and_diff_times',
                'profile': 'normalize_original_float64_values_and_times_in_longdouble',
                'maximum_absolute_value_difference': 2e-17,
                'maximum_absolute_gap_difference': 3e-17,
            },
            'headline_results_replaced': False,
        },
    }


def axes_by_gid(figure):
    result = {axis.get_gid(): axis for axis in figure.axes if axis.get_gid()}
    assert len(result) == len(figure.axes)
    return result


def artist_by_gid(artists, gid):
    matches = [artist for artist in artists if artist.get_gid() == gid]
    assert len(matches) == 1, (gid, [artist.get_gid() for artist in artists])
    return matches[0]


def figure_text(figure):
    return '\n'.join([text.get_text() for text in figure.texts]
                     + [text for axis in figure.axes
                        for text in (axis.get_title(), axis.get_xlabel(), axis.get_ylabel())]
                     + [text.get_text() for axis in figure.axes for text in axis.texts])


def assert_series(axis, gid, x, value, lower, upper, *, ratio=False):
    """Inspect actual line coordinates and CI segments, not an auxiliary spec."""
    x = np.asarray(x, float)
    value, lower, upper = [np.asarray(array, float) for array in (value, lower, upper)]
    line = artist_by_gid(axis.lines, gid)
    displayed = value.copy()
    if ratio:
        displayed[value == 0] = np.nan
    np.testing.assert_allclose(line.get_xdata(), x, rtol=0, atol=1e-14)
    np.testing.assert_allclose(line.get_ydata(), displayed, rtol=1e-13, atol=2e-14, equal_nan=True)
    eligible = np.isfinite(value) & np.isfinite(lower) & np.isfinite(upper)
    if ratio:
        eligible &= (value > 0) & (lower > 0) & (upper > 0)
    expected = np.asarray([[[xx, lo], [xx, hi]]
                           for xx, lo, hi in zip(x[eligible], lower[eligible], upper[eligible])])
    collection = artist_by_gid(axis.collections, gid+':ci')
    actual = np.asarray(collection.get_segments())
    assert actual.shape == expected.shape
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=2e-14)


def wilson(hits, total):
    """Independent closed-form Wilson reference, including its empty domain."""
    if not total:
        return np.nan, np.nan
    fraction = hits/total
    denominator = 1+Z*Z/total
    center = (fraction+Z*Z/(2*total))/denominator
    half = Z*math.sqrt(fraction*(1-fraction)/total+Z*Z/(4*total*total))/denominator
    return max(0., center-half), min(1., center+half)


def assert_floor_figure_matches_saved(figure, stem, design, tables, epsilon):
    """Independently map the actual serialization input to its saved CSV fields."""
    axes = axes_by_gid(figure)
    if stem == 'W7_floor_spacing':
        spacing = tables['spacing_summary.csv'].copy()
        mapping = {cell.identity: cell for cell in design.cells}
        spacing['epsilon'] = spacing.cell_id.map(lambda cid: mapping[cid].floor_fraction)
        fields = ('realized_mean_gap', 'clip_fraction', 'min_gap', 'realized_span')
        assert set(axes) == {'spacing:'+field for field in fields}
        for field in fields:
            scale = 100. if field == 'clip_fraction' else 1.
            for index, floor in enumerate(FLOORS):
                rows = spacing[spacing.epsilon == floor].set_index('cv').loc[[0., 1.5, 2., 3.]]
                value, se = scale*rows[field+'_mean'], scale*rows[field+'_mcse']
                assert_series(axes['spacing:'+field], f'epsilon:{floor:g}',
                              np.arange(4)+(index-1)*.055, value, value-Z*se, value+Z*se)
        return
    kind = {'W7_floor_paired_bias': 'bias', 'W7_floor_paired_coverage': 'coverage',
            'W7_floor_rmse_ratio': 'ratio'}[stem]
    assert set(axes) == {f'{kind}:{floor:g}:{p}' for floor in (1e-8, 1e-4) for p in PARAMETERS}
    for floor in (1e-8, 1e-4):
        for parameter in PARAMETERS:
            metric = ('bias_difference' if parameter == 'mu' else 'relative_bias_difference')
            if kind == 'coverage':
                metric = 'operational_coverage_difference'
            elif kind == 'ratio':
                metric = 'rmse_ratio'
            scale = 1. if kind == 'ratio' or kind == 'bias' and parameter == 'mu' else 100.
            for index, method in enumerate(METHODS):
                rows = epsilon[(epsilon.floor_fraction == floor) & (epsilon.parameter == parameter)
                               & (epsilon.estimator == method)].set_index('cv').loc[[0., 1.5, 2., 3.]]
                assert_series(axes[f'{kind}:{floor:g}:{parameter}'], method,
                              np.arange(4)+(index-1)*.055, scale*rows[metric],
                              scale*rows[metric+'_mc_lower'], scale*rows[metric+'_mc_upper'],
                              ratio=kind == 'ratio')


@pytest.fixture(scope='module')
def followup_figures():
    from ou_irregular.week7 import plots
    design, tables = literal_followup_tables()
    # Reversing every table makes ordering independent of serialization order.
    inputs = {name: value.iloc[::-1].reset_index(drop=True) for name, value in tables.items()}
    before = {name: table.copy(deep=True) for name, table in inputs.items()}
    with forbid_science():
        figures = plots.build_followup_figures(inputs, design)
    for name in inputs:
        pd.testing.assert_frame_equal(inputs[name], before[name])
    yield design, tables, figures
    for figure in figures.values():
        plt.close(figure)


@pytest.fixture(scope='module')
def floor_figures():
    from ou_irregular.week7 import plots
    design, tables, epsilon = literal_floor_tables()
    inputs = {name: value.iloc[::-1].reset_index(drop=True) for name, value in tables.items()}
    paired = epsilon.iloc[::-1].reset_index(drop=True)
    before = {name: table.copy(deep=True) for name, table in inputs.items()}
    paired_before = paired.copy(deep=True)
    with forbid_science():
        figures = plots.build_floor_figures(inputs, paired, design)
    for name in inputs:
        pd.testing.assert_frame_equal(inputs[name], before[name])
    pd.testing.assert_frame_equal(paired, paired_before)
    yield design, tables, epsilon, figures
    for figure in figures.values():
        plt.close(figure)


def test_followup_all_coverage_points_and_wilson_endpoints_match_independent_counts(followup_figures):
    design, tables, figures = followup_figures
    assert set(figures) == {'W7_followup_coverage', 'W7_followup_availability'}
    axes = axes_by_gid(figures['W7_followup_coverage'])
    assert set(axes) == {f'coverage:{scope}:{p}' for scope in ('operational', 'conditional')
                         for p in PARAMETERS}
    intervals_by_method = {'exact': (4, 2), 'pfml': (3, 1), 'euler': (2, 0)}
    for parameter in PARAMETERS:
        for scope in ('operational', 'conditional'):
            for index, method in enumerate(METHODS):
                intervals = intervals_by_method[method]
                hits = tuple(count if parameter == 'mu' else max(0, count-1) for count in intervals)
                totals = (5, 5) if scope == 'operational' else intervals
                estimates = [100*hit/total if total else np.nan for hit, total in zip(hits, totals)]
                bounds = [wilson(hit, total) for hit, total in zip(hits, totals)]
                assert_series(axes[f'coverage:{scope}:{parameter}'], method,
                              np.arange(2)+(index-1)*.055, estimates,
                              [100*value[0] for value in bounds], [100*value[1] for value in bounds])
    text = figure_text(figures['W7_followup_coverage'])
    assert 'pointwise' in text and 'model-based' in text
    assert 'n=250 observations' in text and 'R=5' in text and 'VALIDATION ONLY' in text
    assert 'valid intervals only' in text and 'all R planned paths' in text
    assert tables['parameter_summary.csv'].query("parameter == 'mu'").relative_bias.isna().all()


def test_followup_availability_and_disjoint_failure_bars_count_sources_once(followup_figures):
    _, _, figures = followup_figures
    axes = axes_by_gid(figures['W7_followup_availability'])
    assert set(axes) == {'availability:point', 'availability:interval', 'failures:point_invalid',
                         'failures:interval_unavailable_given_point'}
    point = {'exact': (5, 4), 'pfml': (4, 3), 'euler': (3, 0)}
    interval = {'exact': (4, 2), 'pfml': (3, 1), 'euler': (2, 0)}
    for kind, counts in (('point', point), ('interval', interval)):
        for index, method in enumerate(METHODS):
            bounds = [wilson(count, 5) for count in counts[method]]
            assert_series(axes['availability:'+kind], method, np.arange(2)+(index-1)*.055,
                          100*np.array(counts[method])/5,
                          [100*value[0] for value in bounds], [100*value[1] for value in bounds])
    for category in ('point_invalid', 'interval_unavailable_given_point'):
        for method in METHODS:
            bars = [bar for bar in axes['failures:'+category].patches if bar.get_gid() == method]
            assert len(bars) == 2
            expected = ([5-count for count in point[method]] if category == 'point_invalid' else
                        [p-i for p, i in zip(point[method], interval[method])])
            assert [bar.get_height() for bar in bars] == expected
    assert 'not tripled' in figure_text(figures['W7_followup_availability'])


def test_floor_every_paired_point_and_saved_interval_matches_all_72_keys(floor_figures):
    _, _, epsilon, figures = floor_figures
    assert len(epsilon) == 72 and len(epsilon.columns) == 71
    assert set(figures) == {'W7_floor_paired_bias', 'W7_floor_paired_coverage',
                            'W7_floor_rmse_ratio', 'W7_floor_spacing'}
    for kind, stem in (('bias', 'W7_floor_paired_bias'), ('coverage', 'W7_floor_paired_coverage'),
                        ('ratio', 'W7_floor_rmse_ratio')):
        axes = axes_by_gid(figures[stem])
        assert set(axes) == {f'{kind}:{floor:g}:{parameter}' for floor in (1e-8, 1e-4)
                            for parameter in PARAMETERS}
        for floor in (1e-8, 1e-4):
            for parameter in PARAMETERS:
                metric = ('bias_difference' if parameter == 'mu' else 'relative_bias_difference')
                if kind == 'coverage':
                    metric = 'operational_coverage_difference'
                elif kind == 'ratio':
                    metric = 'rmse_ratio'
                scale = 1. if kind == 'ratio' or kind == 'bias' and parameter == 'mu' else 100.
                for index, method in enumerate(METHODS):
                    rows = epsilon[(epsilon.floor_fraction == floor) & (epsilon.parameter == parameter)
                                   & (epsilon.estimator == method)].set_index('cv').loc[[0., 1.5, 2., 3.]]
                    assert_series(axes[f'{kind}:{floor:g}:{parameter}'], method,
                                  np.arange(4)+(index-1)*.055, scale*rows[metric],
                                  scale*rows[metric+'_mc_lower'], scale*rows[metric+'_mc_upper'],
                                  ratio=kind == 'ratio')


def test_floor_rational_pairing_mu_units_and_log_ratio_limits_are_not_substituted(floor_figures):
    _, tables, epsilon, figures = floor_figures
    for parameter, scale in (('theta', 100.), ('mu', 1.), ('sigma', 200.)):
        axis = axes_by_gid(figures['W7_floor_paired_bias'])[f'bias:1e-08:{parameter}']
        value = artist_by_gid(axis.lines, 'exact').get_ydata()[1]
        assert value == pytest.approx(scale/3, rel=1e-13)
        bars = artist_by_gid(axis.collections, 'exact:ci').get_segments()
        np.testing.assert_allclose(bars[1][:, 1], scale*np.array([1/3-Z*2/3, 1/3+Z*2/3]), rtol=1e-13)
    # The own-valid mean difference is 107/4-26/4=81/4, rather than 1/3.
    assert artist_by_gid(axes_by_gid(figures['W7_floor_paired_bias'])['bias:1e-08:theta'].lines,
                         'exact').get_ydata()[1] != pytest.approx(100*81/4)
    assert epsilon.loc[epsilon.parameter == 'mu', 'relative_bias_difference'].isna().all()
    ratio = math.sqrt(1.5)
    lo, hi = ratio*math.exp(-Z*5/28), ratio*math.exp(Z*5/28)
    axis = axes_by_gid(figures['W7_floor_rmse_ratio'])['ratio:1e-08:mu']
    assert artist_by_gid(axis.lines, 'exact').get_ydata()[1] == pytest.approx(ratio)
    endpoints = artist_by_gid(axis.collections, 'exact:ci').get_segments()[1][:, 1]
    np.testing.assert_allclose(endpoints, [lo, hi], rtol=1e-13)
    assert not np.allclose(endpoints, [ratio-Z*ratio*5/28, ratio+Z*ratio*5/28])
    note = figure_text(figures['W7_floor_paired_bias'])
    assert 'No practical-breach threshold applies to a contrast' in note
    assert 'CV=1.5: 3/3/3' in note and 'CV=3: 5/1/0' in note


def test_floor_all_r_coverage_limits_are_paired_and_not_clipped(floor_figures):
    _, _, _, figures = floor_figures
    axis = axes_by_gid(figures['W7_floor_paired_coverage'])['coverage:1e-08:theta']
    line = artist_by_gid(axis.lines, 'exact')
    # First positive-CV literal: d=(0,1,-1,0,-1), mean=-1/5, sample SE=sqrt(.14).
    assert line.get_ydata()[1] == pytest.approx(-20.)
    bounds = artist_by_gid(axis.collections, 'exact:ci').get_segments()[1][:, 1]
    np.testing.assert_allclose(bounds, 100*np.array([-.2-Z*math.sqrt(.14), -.2+Z*math.sqrt(.14)]))
    # CV2 literal d=(1,1,1,1,-1) has mean=.6, sample SE=.4. Its normal upper
    # endpoint is >1, which remains a visible approximation rather than clipping.
    assert line.get_ydata()[2] == pytest.approx(60.)
    bounds = artist_by_gid(axis.collections, 'exact:ci').get_segments()[2][:, 1]
    np.testing.assert_allclose(bounds, 100*np.array([.6-Z*.4, .6+Z*.4]))
    assert bounds[1] > 100 and axis.get_ylim()[1] >= bounds[1]
    text = figure_text(figures['W7_floor_paired_coverage'])
    assert 'All R=5 planned paired paths' in text and 'pointwise' in text


def test_floor_zero_missing_and_singleton_ratio_domains_are_visible(floor_figures):
    _, _, epsilon, figures = floor_figures
    axes = axes_by_gid(figures['W7_floor_rmse_ratio'])
    low = axes['ratio:1e-08:theta']
    high = axes['ratio:0.0001:theta']
    # CV0 is a shared control: unit ratio and collapsed interval on the same S.
    assert artist_by_gid(low.lines, 'exact').get_ydata()[0] == 1.
    np.testing.assert_array_equal(artist_by_gid(low.collections, 'exact:ci').get_segments()[0][:, 1], [1., 1.])
    # CV2 numerator RMSE is zero; the point is explicitly labelled off the log
    # scale, and no epsilon, unity, or artificial delta interval is substituted.
    assert np.isnan(artist_by_gid(low.lines, 'exact').get_ydata()[2])
    assert any(text.get_text() == '0' and text.get_position()[0] == pytest.approx(2-.055)
               for text in low.texts)
    assert np.isnan(artist_by_gid(low.lines, 'exact').get_ydata()[3])
    assert any(text.get_text() == 'NA' and text.get_position()[0] == pytest.approx(3-.055)
               for text in low.texts)
    # CV3 high-floor PFML has S=1: finite unit ratio, unavailable MC uncertainty.
    assert artist_by_gid(high.lines, 'pfml').get_ydata()[3] == 1.
    assert all(segment[0, 0] != 3. for segment in artist_by_gid(high.collections, 'pfml:ci').get_segments())
    assert 'one shared control' in figure_text(figures['W7_floor_rmse_ratio'])


def test_spacing_uses_each_source_once_and_saved_mean_mcse_without_minimum_substitution(floor_figures):
    design, tables, _, figures = floor_figures
    axes = axes_by_gid(figures['W7_floor_spacing'])
    spacing = tables['spacing_summary.csv'].copy()
    mapping = {cell.identity: cell for cell in design.cells}
    spacing['epsilon'] = spacing.cell_id.map(lambda cid: mapping[cid].floor_fraction)
    expected_fields = ('realized_mean_gap', 'clip_fraction', 'min_gap', 'realized_span')
    assert set(axes) == {'spacing:'+field for field in expected_fields}
    assert spacing.n_planned.eq(5).all()
    for field in expected_fields:
        scale = 100. if field == 'clip_fraction' else 1.
        for index, floor in enumerate(FLOORS):
            rows = spacing[spacing.epsilon == floor].set_index('cv').loc[[0., 1.5, 2., 3.]]
            value, se = scale*rows[field+'_mean'], scale*rows[field+'_mcse']
            assert_series(axes['spacing:'+field], f'epsilon:{floor:g}', np.arange(4)+(index-1)*.055,
                          value, value-Z*se, value+Z*se)
    # Literal mean path minimum is three times its global minimum over sources.
    axis = axes['spacing:min_gap']
    assert artist_by_gid(axis.lines, 'epsilon:1e-08').get_ydata()[1] == pytest.approx(1.5e-8)
    assert artist_by_gid(axis.lines, 'epsilon:1e-08').get_ydata()[1] != pytest.approx(.5e-8)
    assert axis.get_yscale() == 'symlog'
    # Span is built from K=999 transitions, not n=1000 observations.
    means = artist_by_gid(axes['spacing:realized_mean_gap'].lines, 'epsilon:1e-08').get_ydata()
    spans = artist_by_gid(axes['spacing:realized_span'].lines, 'epsilon:1e-08').get_ydata()
    np.testing.assert_allclose(spans, 999*means, rtol=1e-14)
    text = figure_text(figures['W7_floor_spacing'])
    assert 'not three estimator copies' in text and 'mean of path minima' in text
    assert 'K=n-1=999' in text


def test_cv_axes_explicitly_describe_equally_spaced_design_cells(followup_figures, floor_figures):
    # CV=0→1.5 and CV=1.5→2 occupy equal categorical spacings. The axis must
    # disclose that representation rather than imply numerical CV distances.
    for figures in (followup_figures[-1], floor_figures[-1]):
        for figure in figures.values():
            assert all('equally spaced design cells' in axis.get_xlabel() for axis in figure.axes)


def literal_constant_spacing_table():
    """R=2000 scalar summary arithmetic only, with no OU/Gamma observations.

    The constant path-minimum values are the three actual nominal floors.
    NumPy's finite summation can put their computed mean a few ulps above the
    identical observed minimum/maximum; even its tiny computed SD is retained.
    """
    design = floor_design().validation(2000)
    rows = []
    for cell in design.cells:
        row = {'cell_id': cell.identity, 'n': 1000, 'theta_mean_gap': .5,
               'cv': cell.cv, 'n_planned': 2000}
        values = {'realized_mean_gap': .5, 'realized_cv': cell.cv,
                  'realized_span': 499.5, 'clip_fraction': 0. if not cell.cv else 1/999,
                  'min_gap': .5 if not cell.cv else .5*cell.floor_fraction}
        for name, value in values.items():
            constants = np.full(2000, value, dtype=np.float64)
            row.update({name+'_n_available': 2000,
                        name+'_mean': float(np.mean(constants)),
                        name+'_mcse': float(np.std(constants, ddof=1)/math.sqrt(2000)),
                        name+'_min': float(np.min(constants)),
                        name+'_max': float(np.max(constants))})
        rows.append(row)
    return design, {'spacing_summary.csv': pd.DataFrame(rows)}


def test_r2000_constant_minima_roundoff_is_accepted_without_altering_any_saved_scalar():
    from ou_irregular.week7 import plots
    design, tables = literal_constant_spacing_table()
    original = tables['spacing_summary.csv'].copy(deep=True)
    positive = original[original.cv > 0]
    assert len(positive) == 9
    # This actual arithmetic counterexample rejects an exact min<=mean<=max
    # preflight; it is not an invented adjustment to the source or its values.
    assert positive.min_gap_mean.gt(positive.min_gap_max).all()
    assert positive.min_gap_mcse.gt(0.).all()
    assert set(positive.min_gap_min) == {5e-9, 5e-7, 5e-5}
    with forbid_science():
        accepted = plots._spacing_inputs(tables, design)
    pd.testing.assert_frame_equal(tables['spacing_summary.csv'], original, check_exact=True)
    pd.testing.assert_frame_equal(accepted.drop(columns='floor_fraction'), original, check_exact=True)
    assert accepted.min_gap_mean.to_numpy().tobytes() == original.min_gap_mean.to_numpy().tobytes()
    assert accepted.min_gap_mcse.to_numpy().tobytes() == original.min_gap_mcse.to_numpy().tobytes()


@pytest.mark.parametrize('direction', (-1, 1))
def test_r2000_roundoff_allowance_does_not_accept_a_material_extrema_violation(direction):
    from ou_irregular.week7 import plots
    design, tables = literal_constant_spacing_table()
    frame = tables['spacing_summary.csv']
    # A relative 1e-8 deviation is many orders larger than binary64 summation
    # roundoff at R2000, even for the very small 5e-9 physical floor.
    target = frame.index[(frame.cv == 1.5) & frame.min_gap_min.eq(5e-9)][0]
    frame.loc[target, 'min_gap_mean'] = 5e-9*(1+direction*1e-8)
    with forbid_science(), pytest.raises(ValueError, match='extrema/mean'):
        plots._spacing_inputs(tables, design)


def test_zero_extrema_do_not_gain_an_absolute_roundoff_tolerance():
    from ou_irregular.week7 import plots
    design, tables = literal_constant_spacing_table()
    frame = tables['spacing_summary.csv']
    target = frame.index[frame.cv == 0.][0]
    assert frame.loc[target, 'clip_fraction_min'] == frame.loc[target, 'clip_fraction_max'] == 0.
    frame.loc[target, 'clip_fraction_mean'] = 1e-300
    with forbid_science(), pytest.raises(ValueError, match='extrema/mean'):
        plots._spacing_inputs(tables, design)


@pytest.fixture(scope='module')
def selected_figures():
    from ou_irregular.week7 import plots
    result = {}
    with forbid_science():
        for cv in (.5, 0.):
            payload = literal_selected_payload(cv)
            before = copy.deepcopy(payload)
            figure = plots.build_selected_figure(payload, validation=True)
            assert payload == before
            result[cv] = (payload, figure)
    yield result
    for _, figure in result.values():
        plt.close(figure)


def test_selected_profile_subtracts_extended_absolutes_before_float_conversion(selected_figures):
    from ou_irregular.week7 import plots
    for payload, figure in selected_figures.values():
        profile = payload['diagnostic']['exact_profile']
        # Every binary64 absolute is identical: subtracting those would erase
        # every signed difference used by this independent decimal reference.
        assert len(set(profile['grid_nll'])) == 1
        exact = np.array([float(Decimal(value)-Decimal(profile['best_finite']['nll_extended']))
                          for value in profile['grid_nll_extended']])
        with forbid_science():
            series = plots.profile_series(payload)
        # The platform's 64-significand-bit arithmetic rounds absolute values
        # by at most about 3e-14 near 1e6; it is not arbitrary precision.
        np.testing.assert_allclose(series['delta_nll'], exact, rtol=0, atol=6e-14)
        assert series['delta_nll'][4] < 0 and series['delta_nll'][0] > 0
        axes = axes_by_gid(figure)
        for gid, mask in (('profile:global', np.ones(6, dtype=bool)),
                          ('profile:zoom', series['zoom_mask'])):
            curve = artist_by_gid(axes[gid].lines, 'profile-grid')
            np.testing.assert_array_equal(curve.get_xdata(), np.asarray(profile['grid_log_theta'])[mask])
            np.testing.assert_array_equal(curve.get_ydata(), series['delta_nll'][mask])
        assert 'theta_u' in axes['profile:global'].get_xlabel()


def test_selected_stored_and_limit_markers_keep_saved_representation_and_signed_values(selected_figures):
    for payload, figure in selected_figures.values():
        axes = axes_by_gid(figure)
        for key in ('profile:global', 'profile:zoom'):
            axis = axes[key]
            np.testing.assert_allclose(artist_by_gid(axis.collections, 'profile-best').get_offsets(), [[1., 0.]])
            for name, expected in (('stored_profile_minus_best', 1e-11), ('stored_full_minus_best', 5e-11)):
                actual = artist_by_gid(axis.collections, name).get_offsets()
                np.testing.assert_allclose(actual, [[.25, expected]], rtol=1e-14, atol=0)
            np.testing.assert_array_equal(artist_by_gid(axis.lines, 'profile-resolution').get_ydata(), [.0010001]*2)
        for name, expected in (('brownian_limit_minus_best', 6e-11), ('iid_limit_minus_best', -2e-11)):
            np.testing.assert_array_equal(artist_by_gid(axes['profile:global'].lines, name).get_ydata(), [expected]*2)
        text = figure_text(figure)
        assert 'No replacement fit, interval or covariance' in text
        assert 'Finite searches do not prove a global optimum' in text
        assert 'False/False' in text and 'cap limited tail unresolved' in text


def test_selected_all_four_steps_and_three_pairs_retain_zero_negative_and_missing_diagnostics(selected_figures):
    for payload, figure in selected_figures.values():
        axes = axes_by_gid(figure)
        assert set(axes) == {'profile:global', 'profile:zoom', 'hessian:eigenvalues',
                             'hessian:condition', 'hessian:score', 'hessian:adjacent'}
        ladders = {method: payload['diagnostic']['fits'][method]['float64_ladder'] for method in METHODS}
        ladders['exact_longdouble'] = payload['diagnostic']['fits']['exact']['independent_longdouble_ladder']
        for method, ladder in ladders.items():
            for key, field in (('hessian:condition', 'condition_number'), ('hessian:score', 'standardized_score')):
                expected = [np.nan if entry[field] is None else entry[field] for entry in ladder['entries']]
                line = artist_by_gid(axes[key].lines, method)
                np.testing.assert_array_equal(line.get_xdata(), np.arange(4))
                np.testing.assert_allclose(line.get_ydata(), expected, rtol=0, atol=0, equal_nan=True)
                assert axes[key].get_yscale() == 'symlog'
            for field in ('relative_frobenius', 'whitened_spectral'):
                expected = [np.nan if entry[field] is None else entry[field] for entry in ladder['adjacent']]
                line = artist_by_gid(axes['hessian:adjacent'].lines, method+':'+field)
                np.testing.assert_array_equal(line.get_xdata(), np.arange(3))
                np.testing.assert_allclose(line.get_ydata(), expected, rtol=0, atol=0, equal_nan=True)
        for method in ('exact', 'exact_longdouble'):
            for index in range(3):
                expected = [entry['eigenvalues'][index] for entry in ladders[method]['entries']]
                line = artist_by_gid(axes['hessian:eigenvalues'].lines, f'{method}:lambda{index+1}')
                np.testing.assert_array_equal(line.get_xdata(), np.arange(4))
                np.testing.assert_array_equal(line.get_ydata(), expected)
        assert axes['hessian:eigenvalues'].get_yscale() == 'symlog'
        assert axes['hessian:adjacent'].get_yscale() == 'symlog'
        assert [label.get_text() for label in axes['hessian:score'].get_xticklabels()] == [
            '.001', '.0005', '.00025', '.000125']
        assert [label.get_text() for label in axes['hessian:adjacent'].get_xticklabels()] == [
            'h → h/2', 'h/2 → h/4', 'h/4 → h/8']
    failed = axes_by_gid(selected_figures[.5][1])
    np.testing.assert_array_equal(artist_by_gid(failed['hessian:eigenvalues'].lines, 'exact:lambda1').get_ydata(),
                                  [.01, .1, 0., -.2])
    assert artist_by_gid(failed['hessian:score'].lines, 'pfml').get_ydata()[0] == 0.
    assert np.isnan(artist_by_gid(failed['hessian:score'].lines, 'exact').get_ydata()[2])


def test_selected_fixed_cases_are_literal_and_do_not_promote_finer_steps_or_storage_width(selected_figures):
    for cv, (payload, figure) in selected_figures.items():
        text = figure_text(figure)
        assert 'LITERAL VALIDATION FIXTURE' in text and 'not an original E2 replay' in text
        assert f'CV={cv:g}, replication {465 if cv else 1424}' in text
        assert '63 fraction bits' in text and '64 significant binary digits' in text and '16 storage bytes' in text
        assert 'quadruple' not in text and '128 significant' not in text
        assert 'LD profile normalizes original values/gaps separately' in text
        assert 'h/2 is the registered matrix candidate subject to all original gates' in text
        assert 'h/4 and h/8 are diagnostic only' in text
        original = 'step_disagreement' if cv else 'valid'
        assert 'valid_interior / '+original in text
        assert payload['headline_results_replaced'] is False


def test_selected_missing_grid_values_remain_gaps_and_literal_scope_cannot_be_promoted():
    from ou_irregular.week7 import plots
    payload = literal_selected_payload(.5)
    payload['diagnostic']['exact_profile']['grid_nll_extended'][2] = None
    with forbid_science():
        series = plots.profile_series(payload)
        assert np.isnan(series['delta_nll'][2])
        assert series['delta_nll'][4] < 0
        with pytest.raises(ValueError, match='literal'):
            plots.build_selected_figure(payload, validation=False)


def test_selected_retained_failure_is_an_unavailable_figure_not_a_zero_diagnostic():
    from ou_irregular.week7 import plots
    payload = {'selection': literal_selected_payload(.5)['selection'],
               'status': 'original_replay_failed',
               'diagnostic': {'status': 'not_run', 'error': 'literal validation replay-failure fixture'},
               'headline_results_replaced': False, 'literal_validation_fixture': True}
    with forbid_science():
        figure = plots.build_selected_figure(payload, validation=True)
    try:
        assert len(figure.axes) == 6
        assert all(not axis.lines and not axis.collections for axis in figure.axes)
        text = figure_text(figure)
        assert text.count('Saved diagnostics unavailable') == 6
        assert 'literal validation replay-failure fixture' in text
        assert 'saved failure remains unresolved' in text
        assert 'LITERAL VALIDATION FIXTURE' in text
    finally:
        plt.close(figure)


@pytest.mark.parametrize('fault', ('missing_step', 'wrong_pair', 'replacement_interval'))
def test_selected_full_ladder_display_contract_rejects_incomplete_or_relabelled_steps(fault):
    from ou_irregular.week7 import plots
    diagnostic = literal_selected_payload(.5)['diagnostic']
    ladder = diagnostic['fits']['exact']['float64_ladder']
    if fault == 'missing_step':
        ladder['entries'].pop()
    elif fault == 'wrong_pair':
        ladder['adjacent'][1]['fine_relative_step'] = .000125
    else:
        ladder['replacement_interval_constructed'] = True
    with forbid_science(), pytest.raises(ValueError):
        plots._ladder_traces(diagnostic)


def test_actual_saved_floor_render_is_bound_to_existing_audit_and_does_not_change_scientific_bytes(tmp_path_factory):
    from ou_irregular.week7 import plots
    source, design = saved_floor_archive(tmp_path_factory)
    before = snapshot(source)
    output, actual_design = rendered_floor_archive(tmp_path_factory)
    assert actual_design == design
    after_render = snapshot(output)
    with forbid_science():
        manifest = plots.verify_figures('floor', output, design)
    assert snapshot(output) == after_render and snapshot(source) == before
    metadata = json.loads((output/'run_metadata.json').read_text())
    audit = json.loads((output/'audit.json').read_text())
    assert metadata['status'] == 'raw_complete_pending_audit'
    assert metadata['source_policy']['planned_source_realizations'] == 24
    assert metadata['source_policy']['planned_independent_coupling_groups'] == 8
    assert manifest['input_binding']['audit_sha256'] == sha(output/'audit.json')
    assert manifest['input_binding']['raw_input_sha256'] == audit['raw_input_sha256']
    assert manifest['design_sha256'] == design.digest
    assert manifest['new_sources_or_fits'] is False
    assert manifest['headline_results_replaced'] is False
    assert manifest['visual_review_performed'] is False
    assert set(manifest['figures']) == {stem+'.'+extension for stem in (
        'W7_floor_paired_bias', 'W7_floor_paired_coverage', 'W7_floor_rmse_ratio', 'W7_floor_spacing')
        for extension in ('png', 'pdf')}
    for name, digest in manifest['figures'].items():
        assert sha(output/'figures'/name) == digest
    assert not (output/'visual_review.json').exists()
    assert plots.renderer_state_path(output).is_file()

