import json
import shutil

import pandas as pd
import pytest

from ou_irregular.experiments.analysis import write_tables
from ou_irregular.experiments.execution import run_raw
from ou_irregular.experiments.provenance import ROOT, file_digest
from ou_irregular.week6.config import load_config
from ou_irregular.week6.reporting import load_tables, write_grid_reports


@pytest.fixture(scope="module")
def full_grid_tables(tmp_path_factory):
    root = tmp_path_factory.mktemp("e2_reporting_validation")
    design = load_config(ROOT / "configs/week6_e2.yaml").validation(2)
    run_raw(design, root, progress=False)
    write_tables(root, design)
    return root, design


def replace_table(root, name, table):
    table.to_csv(root/name, index=False, float_format="%.17g")
    path = root/"table_manifest.json"
    manifest = json.loads(path.read_text())
    manifest[name] = {"rows": len(table), "sha256": file_digest(root/name)}
    path.write_text(json.dumps(manifest))


def test_all_grid_groups_and_zero_failure_cells_are_retained(tmp_path, full_grid_tables):
    source, design = full_grid_tables
    report = write_grid_reports(source, design, output_dir=tmp_path)
    assert report["reports"]["classification_table.csv"]["rows"] == 540
    assert report["reports"]["classification_counts.csv"]["rows"] == 108
    assert report["reports"]["failed_cells.csv"]["rows"] == 180
    counts = pd.read_csv(tmp_path/"classification_counts.csv")
    assert counts.groupby(["n", "parameter", "estimator"]).n_cells.sum().eq(20).all()
    assert counts.query("parameter == 'mu' and practical_classification == 'not_applicable_zero_truth'").n_cells.eq(20).all()
    failures = pd.read_csv(tmp_path/"failed_cells.csv")
    assert (~failures.has_failure).any()
    assert failures.n_planned.eq(2).all()
    first = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    assert write_grid_reports(source, design, output_dir=tmp_path) == report
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == first


def test_scientific_noncoverage_is_not_a_numerical_failure(tmp_path, full_grid_tables):
    source, design = full_grid_tables
    tables, _ = load_tables(source, design)
    summary = tables["parameter_summary.csv"]
    selected = summary[(summary.n_interval_valid == 2) & (summary.coverage_operational == 0)]
    assert len(selected) > 0
    write_grid_reports(source, design, output_dir=tmp_path)
    failures = pd.read_csv(tmp_path/"failed_cells.csv").set_index(["cell_id", "estimator"])
    for row in selected.itertuples():
        assert not failures.loc[(row.cell_id, row.estimator), "has_failure"]


@pytest.mark.parametrize("corruption", ["missing_whole_group", "duplicate_with_same_length"])
def test_report_keys_checked_against_the_whole_plan(tmp_path, full_grid_tables, corruption):
    source, design = full_grid_tables
    shutil.copytree(source, tmp_path, dirs_exist_ok=True)
    table = pd.read_csv(tmp_path/"parameter_summary.csv", float_precision="round_trip")
    if corruption == "missing_whole_group":
        table = table.iloc[3:]
    else:
        table = pd.concat([table.iloc[1:], table.iloc[[1]]], ignore_index=True)
    replace_table(tmp_path, "parameter_summary.csv", table)
    with pytest.raises(ValueError, match="full planned key set"):
        write_grid_reports(tmp_path, design)


def test_failed_cell_log_avoids_double_counting(tmp_path, full_grid_tables):
    source, design = full_grid_tables
    shutil.copytree(source, tmp_path, dirs_exist_ok=True)
    cid, method = design.cells[0].identity, "exact"
    # Synthetic reporting fixture: one point-invalid path and one valid point
    # without an interval. This does not claim to be an audited scientific run.
    summary = pd.read_csv(tmp_path/"parameter_summary.csv", float_precision="round_trip")
    mask = (summary.cell_id == cid) & (summary.estimator == method)
    summary.loc[mask, "n_point_valid"] = 1
    summary.loc[mask, "n_interval_valid"] = 0
    replace_table(tmp_path, "parameter_summary.csv", summary)
    numerical = pd.read_csv(tmp_path/"numerical_diagnostics.csv", float_precision="round_trip")
    numerical.loc[(numerical.cell_id == cid) & (numerical.estimator == method), "n_interval_valid"] = 0
    replace_table(tmp_path, "numerical_diagnostics.csv", numerical)
    outcomes = pd.read_csv(tmp_path/"outcome_counts.csv", float_precision="round_trip")
    mask = (outcomes.cell_id == cid) & (outcomes.estimator == method)
    a = outcomes[mask].iloc[0].copy();b = a.copy()
    a["point_status"], a["interval_status"], a["count"] = "boundary_candidate", "point_fit_invalid", 1
    b["point_status"], b["interval_status"], b["count"] = "valid_interior", "hessian_unstable", 1
    replace_table(tmp_path, "outcome_counts.csv", pd.concat([outcomes[~mask], pd.DataFrame([a,b])], ignore_index=True))
    write_grid_reports(tmp_path, design)
    log = pd.read_csv(tmp_path/"figures/failed_cells.csv").set_index(["cell_id", "estimator"])
    row = log.loc[(cid, method)]
    assert row.has_failure and row.n_point_invalid == 1 and row.n_interval_unavailable == 2
    assert row.n_point_valid_without_interval == 1 and row.failure_fraction == 1
    assert sum(x["count"] for x in json.loads(row.status_counts_json)) == 2


def test_missing_entire_outcome_group_is_rejected(tmp_path, full_grid_tables):
    source, design = full_grid_tables
    shutil.copytree(source, tmp_path, dirs_exist_ok=True)
    table = pd.read_csv(tmp_path/"outcome_counts.csv", float_precision="round_trip")
    first = table.iloc[0]
    table = table[~((table.cell_id == first.cell_id) & (table.estimator == first.estimator))]
    replace_table(tmp_path, "outcome_counts.csv", table)
    with pytest.raises(ValueError, match="all planned cell/methods"):
        write_grid_reports(tmp_path, design)


def test_unverified_input_bytes_are_rejected(tmp_path, full_grid_tables):
    source, design = full_grid_tables
    shutil.copytree(source, tmp_path, dirs_exist_ok=True)
    with (tmp_path/"parameter_summary.csv").open("a") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="input digest mismatch"):
        write_grid_reports(tmp_path, design)
