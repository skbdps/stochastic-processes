import json

import pytest

from ou_irregular.experiments.analysis import write_tables
from ou_irregular.experiments.execution import run_raw
from ou_irregular.week5.plots import render_figures
from test_execution import small_design


def test_all_requested_and_companion_figures_regenerate_exactly(tmp_path):
    design = small_design()
    run_raw(design, tmp_path, progress=False)
    write_tables(tmp_path, design)
    first = render_figures(tmp_path, design)
    second = render_figures(tmp_path, design)
    assert first == second and len(first["figures"]) == 10
    for name in first["figures"]:
        data = (tmp_path / "figures" / name).read_bytes()
        assert data.startswith(b"\x89PNG") if name.endswith("png") else data.startswith(b"%PDF")


def test_plotting_rejects_changed_input_tables(tmp_path):
    design = small_design()
    run_raw(design, tmp_path, progress=False)
    write_tables(tmp_path, design)
    path = tmp_path / "parameter_summary.csv"
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="digest"):
        render_figures(tmp_path, design)
