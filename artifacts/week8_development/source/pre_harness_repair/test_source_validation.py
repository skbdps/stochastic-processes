"""Independent source-validation arithmetic and artifact-contract checks.

All inputs are literal arrays, analytic moment references, or explicitly mocked
sequence rows. No random schedule, OU observation, fit, or actual 256-sequence
scientific validation is generated. The artifact tests retain the fixed plan
size while substituting a mock source and mock provenance at their boundaries.
"""
from copy import deepcopy
from fractions import Fraction
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from ou_irregular.week8 import source as e4_source
from ou_irregular.week8.design import E4Cell


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("week8_validation_under_independent_test", ROOT / "scripts/week8_validate_source.py")
validation = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = validation
SPEC.loader.exec_module(validation)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


@pytest.fixture(autouse=True)
def no_scientific_generation(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("this test must not generate random schedules or OU observations")
    monkeypatch.setattr(e4_source, "generate_gaps", forbidden)
    monkeypatch.setattr(e4_source, "simulate_ou", forbidden)
    monkeypatch.setattr(np.random, "default_rng", forbidden)


def address(cell, rep):
    words = [int(cell.identity[i:i + 8], 16) for i in range(0, 64, 8)]
    parent = [8, *words, rep]
    return {"source_seed": 22261011, "source_spawn_key": canonical(parent),
            "gap_spawn_key": canonical(parent + [0]), "path_spawn_key": canonical(parent + [1])}


def literal_gap_result(cell, raw):
    """Construct arrays directly; these are deterministic test data, not draws."""
    raw = np.array(raw, dtype=float)
    floor = float(cell.floor_fraction * cell.q / cell.theta if cell.cv >= cell.floor_cv_threshold else 0)
    transformed = np.maximum(raw, floor)
    phases = np.zeros(len(raw), dtype=np.int64)
    return e4_source.GapResult(cell.identity, raw, transformed,
                              np.concatenate(([0.0], np.cumsum(transformed))),
                              phases, phases, phases, floor)


def test_fixed_source_only_plan_and_complete_scientific_cell_hashes():
    plan = validation.validation_plan()
    assert plan["scope"] == "source_only_validation"
    assert plan["seed"] == 22261011 and plan["source_namespace"] == 8
    assert plan["schedule_child"] == 0
    assert plan["replications_per_cell"] == 256
    assert plan["gaps_per_sequence"] == 999
    assert plan["adjacent_pairs_per_sequence"] == 998
    assert plan["planned_sequences"] == 1536
    assert plan["planned_raw_gaps"] == 1534464
    assert plan["planned_adjacent_pairs"] == 1532928
    assert plan["cdf_t_over_m"] == [.1, .5, 1., 2., 5.]
    assert plan["cdf_indicator"] == "gap<=t"
    assert plan["moment_cdf_gate"] == {
        "mcse_multiplier": 6., "absolute_tolerance": 1e-10,
        "mcse": "across_sequence_sample_sd_ddof1_div_sqrt256"}
    assert plan["phase_gate"] == {
        "observations_per_cell": 256, "phase_observation": "first_phase_only",
        "standard_error": "theoretical_binomial", "multiplier": 6.}
    assert plan["clipping_gate"] == "none_report_only"
    assert plan["new_ou_observations"] == plan["new_fits"] == 0
    assert not {"estimators", "planned_fit_records"} & set(plan)
    assert [(c["mechanism"], c["ratio"]) for c in plan["cells"]] == [
        ("mmpp", 1.), ("mmpp", 5.), ("mmpp", 20.),
        ("iid", 1.), ("iid", 5.), ("iid", 20.)]
    for cell in plan["cells"]:
        assert cell["cell_id"] == hashlib.sha256(canonical({k: v for k, v in cell.items() if k != "cell_id"}).encode()).hexdigest()
    plan["cells"].clear()
    assert len(validation.validation_plan()["cells"]) == 6


@pytest.mark.parametrize("ratio,second,mmpp_adjacent", [
    (1., Fraction(1, 2), Fraction(1, 4)),
    (5., Fraction(2, 3), Fraction(5, 16)),
    (20., Fraction(290, 219), Fraction(11105, 21316)),
])
def test_raw_targets_match_hand_derived_rational_moments(ratio, second, mmpp_adjacent):
    for mechanism in ("mmpp", "iid"):
        cell = E4Cell(n=1000, q=.5, mechanism=mechanism, ratio=ratio)
        targets = validation.expected_statistics(cell)
        assert len(targets) == 16
        assert targets["raw_mean_gap"] == pytest.approx((.5, 999), rel=2e-14)
        assert targets["raw_second_moment"] == pytest.approx((float(second), 999), rel=2e-14)
        adjacent = float(mmpp_adjacent) if mechanism == "mmpp" else .25
        assert targets["raw_adjacent_product"] == pytest.approx((adjacent, 998), rel=2e-14)
        for name, (_, terms) in targets.items():
            assert terms == (998 if name.endswith("adjacent_product") else 999)
        if ratio == 1:
            for t_over_m, label in ((.1, "0p1"), (.5, "0p5"), (1., "1"), (2., "2"), (5., "5")):
                reference = -math.expm1(-t_over_m)
                for kind in ("raw", "transformed"):
                    assert targets[f"{kind}_cdf_t_over_m_{label}"] == pytest.approx((reference, 999), rel=2e-14)
            for suffix in ("mean_gap", "second_moment", "adjacent_product"):
                assert targets["transformed_" + suffix] == pytest.approx(targets["raw_" + suffix], rel=2e-14)


def test_999_gap_moments_and_998_adjacent_products_have_independent_denominators():
    cell = E4Cell(n=1000, q=.5, mechanism="mmpp", ratio=1.)
    result = literal_gap_result(cell, [1.] * 998 + [3.])
    observed = validation.gap_statistics(cell, result)
    # Sums are 998+3=1001, 998+9=1007, and 997*1+1*3=1000.
    for kind in ("raw", "transformed"):
        assert observed[kind + "_mean_gap"] == pytest.approx(1001 / 999, rel=2e-15)
        assert observed[kind + "_second_moment"] == pytest.approx(1007 / 999, rel=2e-15)
        assert observed[kind + "_adjacent_product"] == pytest.approx(1000 / 998, rel=2e-15)
        assert observed[kind + "_cdf_t_over_m_2"] == pytest.approx(998 / 999, rel=2e-15)
        assert observed[kind + "_cdf_t_over_m_5"] == pytest.approx(998 / 999, rel=2e-15)


def test_literal_floor_statistics_cdf_equality_and_clip_accounting():
    # Off-grid deterministic fixture makes the floor effect visible without
    # changing any fixed validation cell or drawing rare-event observations.
    cell = E4Cell(n=5, q=2., mechanism="mmpp", ratio=20., floor_fraction=.5)
    result = literal_gap_result(cell, [.5, 1., 2., 3.])
    observed = validation.gap_statistics(cell, result)
    assert observed["raw_mean_gap"] == pytest.approx(13 / 8)
    assert observed["raw_second_moment"] == pytest.approx(57 / 16)
    assert observed["raw_adjacent_product"] == pytest.approx(17 / 6)
    assert observed["transformed_mean_gap"] == pytest.approx(7 / 4)
    assert observed["transformed_second_moment"] == pytest.approx(15 / 4)
    assert observed["transformed_adjacent_product"] == pytest.approx(3.)
    assert observed["clip_count"] == observed["raw_below_floor_count"] == 1
    assert observed["clip_fraction"] == .25
    assert observed["transformed_at_floor_count"] == 2
    assert observed["transformed_below_floor_count"] == 0
    for kind in ("raw", "transformed"):
        assert [observed[f"{kind}_cdf_t_over_m_{label}"] for label in ("0p1", "0p5", "1", "2", "5")] == [0., .5, .75, 1., 1.]


def test_cluster_mcse_uses_256_sequence_statistics_and_ddof_one():
    values = [-1.] * 128 + [1.] * 128
    result = validation.comparison(values, 0.)
    # Sum of centered squares=256; sample variance=256/255; SE^2=1/255.
    assert result["observed"] == 0.
    assert result["standard_error"] == pytest.approx(1 / math.sqrt(255), rel=2e-15)
    assert result["tolerance"] == pytest.approx(6 / math.sqrt(255) + 1e-10, rel=2e-15)
    assert result["standard_error_method"] == "across_sequence_sample_sd_ddof1_over_sqrt256"
    assert result["passed"]
    threshold = 6 / math.sqrt(255) + 1e-10
    assert validation.comparison(values, threshold - 1e-8)["passed"]
    assert not validation.comparison(values, threshold + 1e-8)["passed"]


def test_zero_mcse_retains_absolute_tolerance_and_inclusive_boundary():
    values = [0.] * 256
    assert validation.comparison(values, 1e-10)["passed"]
    assert not validation.comparison(values, math.nextafter(1e-10, math.inf))["passed"]
    result = validation.comparison(values, .001)
    assert result["standard_error"] == 0. and result["tolerance"] == 1e-10
    assert not result["passed"]


@pytest.mark.parametrize("values", [[0.] * 255, [0.] * 257, [[0.]] * 256, [float("nan")] * 256, [float("inf")] * 256])
def test_comparison_rejects_wrong_sequence_count_or_nonfinite_values(values):
    with pytest.raises(ValueError):
        validation.comparison(values, 0.)


def test_initial_phase_gate_uses_theoretical_binomial_se_with_n_256():
    inside = validation.comparison([1.] * 243 + [0.] * 13, .8, theoretical_probability=.8)
    outside = validation.comparison([1.] * 244 + [0.] * 12, .8, theoretical_probability=.8)
    assert inside["standard_error"] == pytest.approx(.025, rel=2e-15)
    assert inside["tolerance"] == pytest.approx(.15, rel=2e-15)
    assert inside["standard_error_method"] == "theoretical_binomial_first_phase_N256"
    assert inside["passed"] and not outside["passed"]
    constant = validation.comparison([1.] * 256, .8, theoretical_probability=.8)
    assert constant["standard_error"] == pytest.approx(.025, rel=2e-15)
    with pytest.raises(ValueError):
        validation.comparison([.5] * 256, .8, theoretical_probability=.8)


@pytest.fixture(scope="module")
def literal_rows():
    """Synthetic algebra fixtures for schema/gates; not generated data.

    Their moment/CDF values are assigned targets to exercise the evaluator and
    storage boundaries. Independent target arithmetic is checked separately
    above. They must never be described as source-validation evidence.
    """
    rows = []
    for cell in validation.validation_cells():
        targets = validation.expected_statistics(cell)
        m, floor = cell.q / cell.theta, e4_source.floor_value(cell)
        base = {name: 0 for name in e4_source.GAP_FIELDS}
        base.update({name: expected for name, (expected, _) in targets.items()})
        for kind in ("raw", "transformed"):
            first, second = base[kind + "_mean_gap"], base[kind + "_second_moment"]
            cv = math.sqrt(second - first * first) / first
            base.update({kind + "_cv": cv, kind + "_span": 999 * first,
                         kind + "_min_gap": .01,
                         "expected_" + kind + "_mean_gap": first,
                         "expected_" + kind + "_second_moment": second,
                         "expected_" + kind + "_cv": cv,
                         "expected_" + kind + "_adjacent_product": base[kind + "_adjacent_product"]})
        base.update(family="event_stationary_mmpp" if cell.mechanism == "mmpp" else "iid_full_palm_marginal",
                    target_mean_gap=m, requested_cv=cell.cv, floor_applied=cell.ratio == 20., floor_value=floor,
                    expected_mean_gap_after_floor=base["transformed_mean_gap"],
                    expected_clip_fraction=e4_source.cdf(cell, floor), clip_fraction=0.,
                    realized_mean_gap=base["transformed_mean_gap"], realized_cv=base["transformed_cv"],
                    realized_span=base["transformed_span"], min_gap=.01,
                    raw_gaps_sha256="0" * 64, transformed_gaps_sha256="0" * 64, phase_trace_sha256="0" * 64,
                    source_phase_reset_count=1 if cell.mechanism == "mmpp" else 999)
        calm_count = round(256 * 5 / (cell.ratio + 5))
        for rep in range(256):
            phase = int(rep >= calm_count)
            row = {**cell.to_dict(), **base, **address(cell, rep), "rep": rep,
                   "status": "evaluated", "error": "", "n_gaps": 999, "adjacent_pairs": 998,
                   "times_sha256": "0" * 64, "source_initial_phase": phase,
                   "source_calm_start_count": 999 if phase == 0 else 0,
                   "source_burst_start_count": 999 if phase == 1 else 0,
                   "source_calm_end_count": 999 if phase == 0 else 0,
                   "source_burst_end_count": 999 if phase == 1 else 0}
            rows.append(row)
    return tuple(rows)


def test_complete_evaluation_preserves_all_108_check_rows_and_has_no_clipping_gate(literal_rows):
    checks = validation.evaluate_rows(list(literal_rows))
    assert len(checks) == 108 and all(row["passed"] for row in checks)
    assert len({(row["cell_id"], row["metric"]) for row in checks}) == 108
    assert all(row["n_sequences"] == 256 for row in checks)
    assert not any("clip" in row["metric"] for row in checks)
    assert all(set(row) == set(validation.CHECK_COLUMNS) for row in checks)
    for row in checks:
        expected_terms = 1 if row["metric"].startswith("initial_phase_") else 998 if row["metric"].endswith("adjacent_product") else 999
        assert row["terms_per_sequence"] == expected_terms


def test_phase_comparison_does_not_pool_dependent_event_phases(literal_rows):
    rows = deepcopy(literal_rows)
    for row in rows:
        # Deliberately extreme counts after the first event cannot change a
        # first-phase-only comparison. They remain conformable counts.
        row["source_calm_start_count"] = 999 if row["source_initial_phase"] == 0 else 998
        row["source_burst_start_count"] = 0 if row["source_initial_phase"] == 0 else 1
    checks = validation.evaluate_rows(rows)
    phase_checks = [row for row in checks if row["metric"].startswith("initial_phase_")]
    assert len(phase_checks) == 12 and all(row["passed"] for row in phase_checks)
    for row in phase_checks:
        p = row["expected"]
        assert row["standard_error"] == pytest.approx(math.sqrt(p * (1 - p) / 256), rel=2e-15)


@pytest.mark.parametrize("change", [
    {"rep": False}, {"rep": 0.0}, {"n_gaps": 998}, {"adjacent_pairs": 999},
    {"mechanism": "iid"}, {"ratio": 5.}, {"source_seed": 20261011},
    {"source_spawn_key": "[8,0]"}, {"source_initial_phase": 2},
    {"source_initial_phase": False}, {"source_initial_phase": 0.0},
    {"status": "unknown"}, {"undeclared": "extra"},
])
def test_evaluator_rejects_invalid_identity_schema_counts_or_phase(literal_rows, change):
    rows = list(literal_rows)
    rows[0] = {**rows[0], **change}
    with pytest.raises(ValueError):
        validation.evaluate_rows(rows)


def test_evaluator_rejects_missing_duplicate_and_absent_required_fields(literal_rows):
    with pytest.raises(ValueError):
        validation.evaluate_rows(list(literal_rows[:-1]))
    with pytest.raises(ValueError):
        validation.evaluate_rows([*literal_rows[:-1], literal_rows[0]])
    rows = list(literal_rows)
    rows[0] = {k: v for k, v in rows[0].items() if k != "n_gaps"}
    with pytest.raises(ValueError):
        validation.evaluate_rows(rows)


def test_source_failure_retains_planned_denominator_and_unavailable_checks(literal_rows):
    rows = list(literal_rows)
    rows[0] = {**rows[0], "status": "source_failure", "error": "LiteralFailure: retained fixture"}
    checks = validation.evaluate_rows(rows)
    affected = [row for row in checks if row["cell_id"] == rows[0]["cell_id"]]
    assert len(checks) == 108 and len(affected) == 18
    assert all(not row["passed"] and row["n_sequences"] == 256
               and row["status"] == "source_failure_or_unfinished_sequence"
               and row["observed"] is None for row in affected)
    assert all(row["passed"] for row in checks if row["cell_id"] != rows[0]["cell_id"])


@pytest.fixture
def mocked_run(monkeypatch, literal_rows):
    """Mock source/provenance only; production row writing and checks remain real."""
    lookup = {(row["cell_id"], row["rep"]): row for row in literal_rows}
    calls = {"addresses": [], "gaps": [], "capture": [], "verify": []}
    settings = {"source_failure": None}
    cells = {cell.identity: cell for cell in validation.validation_cells()}

    def capture(identity, extra_files=(), require_frozen=False):
        assert require_frozen is True
        calls["capture"].append(identity.to_dict())
        return {"source_files": {name: "0" * 64 for name in set(validation.CORE_FILES) | set(extra_files)},
                "identity_sha256": "f" * 64, "code_commit": "MOCK_ONLY_NOT_A_PUBLISHED_FREEZE"}

    def verify(provenance, expected_design):
        assert provenance["identity_sha256"] == "f" * 64
        assert expected_design.to_dict() == validation.validation_plan()
        calls["verify"].append(True)

    class MockSource:
        def __init__(self, seed):
            assert seed == 22261011

        def address(self, mapping, rep):
            calls["addresses"].append((mapping["cell_id"], rep))
            return address(cells[mapping["cell_id"]], rep)

        def gaps(self, mapping, rep):
            key = (mapping["cell_id"], rep)
            calls["gaps"].append(key)
            if key == settings["source_failure"]:
                raise RuntimeError("literal mocked generation failure; do not retry")
            return SimpleNamespace(key=key, raw_gaps=np.ones(999), times=np.arange(1000, dtype=float))

        def sample(self, *args, **kwargs):
            raise AssertionError("source-only validation must not request OU samples")

    def statistics(cell, gaps):
        row = lookup[gaps.key]
        assert row["cell_id"] == cell.identity
        return {name: row[name] for name in e4_source.GAP_FIELDS}

    monkeypatch.setattr(validation, "capture", capture)
    monkeypatch.setattr(validation, "verify", verify)
    monkeypatch.setattr(validation, "E4Source", MockSource)
    monkeypatch.setattr(validation, "gap_statistics", statistics)
    return settings, calls


def test_freeze_failure_precedes_output_or_source_allocation(tmp_path, monkeypatch):
    output = tmp_path / "blocked"
    def not_frozen(*args, **kwargs):
        assert kwargs["require_frozen"] is True
        raise ValueError("literal uncommitted-dependency fixture")
    monkeypatch.setattr(validation, "capture", not_frozen)
    monkeypatch.setattr(validation, "E4Source", lambda *args: pytest.fail("source allocated before freeze"))
    with pytest.raises(ValueError, match="uncommitted-dependency"):
        validation.run_validation(output)
    assert not output.exists()


def test_mock_run_writes_complete_immutable_rows_and_no_ou_or_fits(tmp_path, mocked_run):
    settings, calls = mocked_run
    output = tmp_path / "mock-only"
    report = validation.run_validation(output)
    assert report["status"] == "passed"
    assert report["planned_sequences"] == report["retained_sequence_rows"] == report["evaluated_sequences"] == 1536
    assert report["check_rows"] == 108
    assert report["new_ou_observations"] == report["new_fits"] == 0
    assert report["clipping_test"] == "reported_only_no_empirical_normal_gate"
    expected_keys = [(cell.identity, rep) for cell in validation.validation_cells() for rep in range(256)]
    assert calls["addresses"] == calls["gaps"] == expected_keys
    assert len(calls["capture"]) == len(calls["verify"]) == 1
    partitions = sorted(output.glob("sequences-*.jsonl.gz"))
    assert len(partitions) == 6
    all_keys = set()
    for path in partitions:
        with gzip.open(path, "rt") as stream:
            rows = [json.loads(line) for line in stream]
        assert len(rows) == 256
        all_keys.update((row["cell_id"], row["rep"]) for row in rows)
    assert all_keys == set(expected_keys)
    assert json.loads((output / "validation.json").read_text()) == report
    for name, expected in report["artifacts"].items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == expected
    before = {path.name: path.read_bytes() for path in output.iterdir()}
    with pytest.raises(FileExistsError):
        validation.run_validation(output)
    assert {path.name: path.read_bytes() for path in output.iterdir()} == before


def test_mock_source_failure_is_retained_once_with_fixed_plan(tmp_path, mocked_run):
    settings, calls = mocked_run
    key = (validation.validation_cells()[0].identity, 0)
    settings["source_failure"] = key
    output = tmp_path / "mock-failure"
    report = validation.run_validation(output)
    assert report["status"] == "failed" and report["evaluated_sequences"] == 1535
    assert report["retained_sequence_rows"] == 1536 and report["check_rows"] == 108
    assert len(report["source_failures"]) == 1
    failure = report["source_failures"][0]
    assert (failure["cell_id"], failure["rep"]) == key
    assert failure["status"] == "source_failure" and "do not retry" in failure["error"]
    assert calls["gaps"].count(key) == 1
    assert len(report["failed_checks"]) == 18
    with gzip.open(output / f"sequences-{key[0]}.jsonl.gz", "rt") as stream:
        retained = json.loads(next(stream))
    assert retained["source_seed"] == 22261011 and retained["n_gaps"] == 999
    assert retained["status"] == "source_failure" and retained["error"] == failure["error"]
    assert json.loads((output / "validation.json").read_text())["status"] == "failed"


def test_mock_check_construction_exception_still_preserves_failed_report(tmp_path, mocked_run, monkeypatch):
    def broken_checks(rows):
        raise ValueError("literal comparison-construction failure")
    monkeypatch.setattr(validation, "evaluate_rows", broken_checks)
    output = tmp_path / "mock-check-failure"
    report = validation.run_validation(output)
    assert report["status"] == "failed" and report["retained_sequence_rows"] == 1536
    assert "comparison-construction failure" in report["fatal_error"]
    assert (output / "sequence_summary.csv").is_file()
    assert json.loads((output / "validation.json").read_text())["status"] == "failed"


def test_output_rejects_symlinks_and_scientific_archive_descendants(tmp_path):
    actual = tmp_path / "actual"
    actual.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(actual, target_is_directory=True)
    with pytest.raises(ValueError, match="symbolic"):
        validation.run_validation(alias / "new")
    (actual / "validation_plan.json").write_text("{}\n")
    with pytest.raises(ValueError, match="scientific evidence"):
        validation.run_validation(actual / "new")


def test_cli_exposes_no_seed_replication_tolerance_resume_or_overwrite_override(tmp_path, monkeypatch):
    monkeypatch.setattr(validation, "run_validation", lambda *args: pytest.fail("invalid CLI must not execute"))
    for option in ("--seed", "--replications", "--tolerance", "--resume", "--overwrite"):
        with pytest.raises(SystemExit) as error:
            validation.main(["--output", str(tmp_path / "never"), option, "1"])
        assert error.value.code == 2
