"""Independent direct Gamma/OU references for the fixed paired-floor source.

Only the separately specified validation seed 21261009 and replications 0/1
are used. These fixtures are not exploratory results. The references construct
the seed address, Gamma draws, timestamps, stationary X0, and exact OU recursion
without calling any production source, spacing, or simulation function.
"""
import hashlib
import itertools
import json

import numpy as np
import pytest
from scipy.stats import gamma

from ou_irregular.experiments.design import Cell, ESTIMATORS
from ou_irregular.experiments.execution import registered_wald, source_records
from ou_irregular.experiments.provenance import canonical
from ou_irregular.week7.floor_source import PairedFloorSource, anchor_cell


VALIDATION_SEED = 21261009
CVS = (0., 1.5, 2., 3.)
EPSILONS = (1e-8, 1e-6, 1e-4)
REPS = (0, 1)
COUPLING = 'week7_paired_floor_v1_shared_gamma_stationary_x0_and_standard_normal_innovations'


def manual_cell(cv, epsilon):
    """Independent canonical scientific descriptor, including actual epsilon."""
    fields = {'n': 1000, 'theta_mean_gap': .5, 'mean_gap': .5, 'cv': float(cv),
              'true_theta': 1., 'true_mu': 0., 'true_sigma': .5,
              'floor_fraction': float(epsilon), 'floor_cv_threshold': 1.5}
    encoded = json.dumps(fields, sort_keys=True, separators=(',', ':'), allow_nan=False)
    return {'cell_id': hashlib.sha256(encoded.encode()).hexdigest(), **fields}


def independent_hash(values):
    return hashlib.sha256(np.asarray(values, dtype='<f8').tobytes()).hexdigest()


def independent_moments(values):
    """Population moments in extended arithmetic, independently of scaled _moments."""
    values = np.asarray(values, dtype=np.longdouble)
    mean = np.sum(values)/len(values)
    cv = np.sqrt(np.sum((values-mean)**2)/len(values))/mean
    return float(mean), float(cv)


def manual_reference(cv, epsilon, rep):
    actual, anchor = manual_cell(cv, epsilon), manual_cell(cv, 1e-6)
    words = tuple(int(anchor['cell_id'][start:start+8], 16) for start in range(0, 64, 8))
    sequence = np.random.SeedSequence(VALIDATION_SEED, spawn_key=(4, *words, rep))
    gap_seed, path_seed = sequence.spawn(2)
    gap_rng, path_rng = np.random.default_rng(gap_seed), np.random.default_rng(path_seed)
    if cv == 0:
        raw = np.full(999, .5)
        floored = raw.copy()
        times = np.arange(1000, dtype=float)*.5
        floor_value, clip_count, expected_clip, expected_mean = 0., 0, 0., .5
        raw_mean = realized_mean = .5
        raw_cv = realized_cv = 0.
    else:
        # Fixed rational shape/scale values follow k=1/CV², s=.5*CV².
        shape, scale = {1.5: (4/9, 9/8), 2.: (1/4, 2.), 3.: (1/9, 9/2)}[cv]
        raw = gap_rng.gamma(shape, scale, size=999)
        floor_value = epsilon*.5
        clip_count = int(np.count_nonzero(raw < floor_value))
        floored = np.maximum(raw, floor_value)
        times = np.concatenate(([0.], np.cumsum(floored)))
        raw_mean, raw_cv = independent_moments(raw)
        realized_mean, realized_cv = independent_moments(np.diff(times))
        expected_clip = float(gamma.cdf(floor_value, a=shape, scale=scale))
        # Integrating d*f_k(d)=k*s*f_{k+1}(d) above the floor gives this
        # independent distribution-level mean; the sample is never rescaled.
        expected_mean = float(floor_value*expected_clip
                              +.5*gamma.sf(floor_value, a=shape+1, scale=scale))
    # The model consumes the stationary draw first, then K standard normals.
    # Its transitions use differences of final float64 timestamps, preserving
    # the historical cumulative-rounding semantics instead of substituting raw
    # floored gaps into the simulator.
    x0 = path_rng.normal(0., np.sqrt(.125))
    innovations = path_rng.standard_normal(999)
    realized_gaps = np.diff(times)
    correlation = np.exp(-realized_gaps)
    innovation_sd = np.sqrt(.125*(-np.expm1(-2*realized_gaps)))
    values = np.empty(1000)
    values[0] = x0
    for i in range(999):
        values[i+1] = values[i]*correlation[i]+innovation_sd[i]*innovations[i]
    encode_key = lambda key: json.dumps(key, separators=(',', ':'))
    diagnostics = {
        'family': 'equidistant' if cv == 0 else 'gamma',
        'target_mean_gap': .5, 'requested_cv': cv,
        'floor_applied': cv >= 1.5, 'floor_value': floor_value,
        'clip_count': clip_count, 'clip_fraction': clip_count/999,
        'raw_mean_gap': raw_mean, 'realized_mean_gap': realized_mean,
        'raw_cv': raw_cv, 'realized_cv': realized_cv,
        'realized_span': float(times[-1]), 'min_gap': float(np.min(realized_gaps)),
        'expected_mean_gap_after_floor': expected_mean, 'expected_clip_fraction': expected_clip,
        'source_seed': VALIDATION_SEED, 'source_spawn_key': encode_key(sequence.spawn_key),
        'gap_spawn_key': encode_key(gap_seed.spawn_key),
        'path_spawn_key': encode_key(path_seed.spawn_key),
        'times_sha256': independent_hash(times), 'observations_sha256': independent_hash(values),
        'actual_cell_id': actual['cell_id'], 'seed_anchor_cell_id': anchor['cell_id'],
        'seed_anchor_floor_fraction': 1e-6, 'source_coupling_policy': COUPLING,
    }
    return {'cell': actual, 'anchor': anchor, 'raw': raw, 'floored': floored,
            'times': times, 'values': values, 'x0': x0, 'innovations': innovations,
            'diagnostics': diagnostics}


@pytest.fixture(scope='module')
def validation_sources():
    source = PairedFloorSource(VALIDATION_SEED)
    result = {}
    for cv, epsilon, rep in itertools.product(CVS, EPSILONS, REPS):
        reference = manual_reference(cv, epsilon, rep)
        sample, diagnostics = source.sample(reference['cell'], rep)
        result[(cv, epsilon, rep)] = (sample, diagnostics, reference)
    return result


@pytest.mark.parametrize('cv,epsilon,rep', tuple(itertools.product(CVS, EPSILONS, REPS)))
def test_every_frozen_variant_matches_direct_seed_gamma_ou_arrays_and_diagnostics(
        validation_sources, cv, epsilon, rep):
    sample, diagnostics, reference = validation_sources[(cv, epsilon, rep)]
    np.testing.assert_array_equal(sample.times, reference['times'])
    np.testing.assert_array_equal(sample.values, reference['values'])
    assert independent_hash(sample.times) == reference['diagnostics']['times_sha256']
    assert independent_hash(sample.values) == reference['diagnostics']['observations_sha256']
    assert set(diagnostics) == set(reference['diagnostics'])
    for name, expected in reference['diagnostics'].items():
        actual = diagnostics[name]
        if isinstance(expected, float):
            assert actual == pytest.approx(expected, rel=5e-13, abs=5e-15), name
        else:
            assert type(actual) is type(expected) and actual == expected, name
    assert 'cell_id' not in diagnostics
    # Equality to this direct, unrescaled timestamp construction excludes any
    # hidden normalization of the observed mean or realized span.
    if cv > 0:
        assert not np.isclose(sample.times[-1], 999*.5, rtol=0., atol=1e-8)


def test_canonical_anchor_ids_and_namespace_cardinality():
    actual_ids, anchor_ids = set(), set()
    for cv, epsilon in itertools.product(CVS, EPSILONS):
        actual = manual_cell(cv, epsilon)
        canonical_cell = Cell(n=1000, q=.5, cv=cv, floor_fraction=epsilon)
        assert actual == canonical_cell.to_dict()
        expected_anchor = Cell(n=1000, q=.5, cv=cv, floor_fraction=1e-6)
        assert expected_anchor.to_dict() == manual_cell(cv, 1e-6)
        assert anchor_cell(actual) == expected_anchor
        actual_ids.add(actual['cell_id'])
        anchor_ids.add(expected_anchor.identity)
    assert len(actual_ids) == 12
    assert len(anchor_ids) == 4


@pytest.mark.parametrize('cv,rep', tuple(itertools.product(CVS, REPS)))
def test_common_random_numbers_are_preserved_while_physical_variants_follow_floor(
        validation_sources, cv, rep):
    baseline = validation_sources[(cv, 1e-6, rep)]
    sample0, diagnostics0, reference0 = baseline
    for epsilon in EPSILONS:
        sample, diagnostics, reference = validation_sources[(cv, epsilon, rep)]
        np.testing.assert_array_equal(reference['raw'], reference0['raw'])
        np.testing.assert_array_equal(reference['innovations'], reference0['innovations'])
        assert sample.values[0] == reference0['x0'] == sample0.values[0]
        for key in ('source_spawn_key', 'gap_spawn_key', 'path_spawn_key', 'seed_anchor_cell_id'):
            assert diagnostics[key] == diagnostics0[key]
        gaps = np.diff(sample.times)
        recovered = ((sample.values[1:]-sample.values[:-1]*np.exp(-gaps))
                     /np.sqrt(.125*(-np.expm1(-2*gaps))))
        np.testing.assert_allclose(recovered, reference0['innovations'], rtol=2e-10, atol=3e-11)
    for epsilon_a, epsilon_b in itertools.combinations(EPSILONS, 2):
        sample_a, _, reference_a = validation_sources[(cv, epsilon_a, rep)]
        sample_b, _, reference_b = validation_sources[(cv, epsilon_b, rep)]
        if np.array_equal(reference_a['floored'], reference_b['floored']):
            # Equal finite realizations, including high-CV no-change pairs,
            # are legitimate. No difference is forced to create a result.
            np.testing.assert_array_equal(sample_a.times, sample_b.times)
            np.testing.assert_array_equal(sample_a.values, sample_b.values)
        else:
            assert not np.array_equal(sample_a.times, sample_b.times)
            assert not np.array_equal(sample_a.values, sample_b.values)


@pytest.mark.parametrize('rep', REPS)
def test_regular_control_has_identical_complete_engine_fit_results(rep):
    """CV0 is one shared sample, not three independent control datasets."""
    source = PairedFloorSource(VALIDATION_SEED)
    reference = None
    for epsilon in EPSILONS:
        cell = Cell(n=1000, q=.5, cv=0., floor_fraction=epsilon)
        _, diagnostic = source.sample(cell.to_dict(), rep)
        rows = source_records(cell, rep, source,
                              {name: registered_wald(name) for name in ESTIMATORS})
        assert len(rows) == 3
        assert {row['estimator'] for row in rows} == set(ESTIMATORS)
        assert all(row['source_valid'] and not row['source_error'] and not row['pipeline_error']
                   for row in rows)
        metadata = set(cell.to_dict()) | set(diagnostic)
        # Compare every remaining field, including all Hessian/covariance/SE/
        # interval numbers, missing values, validity/status and embedded JSON.
        # Only physical-design/source labels are excluded from this comparison.
        numerical_records = canonical([{k: v for k, v in row.items() if k not in metadata}
                                       for row in rows])
        if reference is None:
            reference = numerical_records
        else:
            assert numerical_records == reference
