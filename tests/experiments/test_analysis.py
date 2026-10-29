from pathlib import Path

import pandas as pd
import pytest

from ou_irregular.experiments.analysis import write_tables
from ou_irregular.experiments.execution import run_raw
from test_execution import small_design


def test_every_cell_survives_exports_and_exports_regenerate_exactly(tmp_path):
    design = small_design()
    run_raw(design, tmp_path, progress=False)
    one = write_tables(tmp_path, design)
    two = write_tables(tmp_path, design)
    assert one == two
    assert one["parameter_summary.csv"]["rows"] == 18
    assert one["regular_comparisons.csv"]["rows"] == 18
    assert one["paired_comparisons.csv"]["rows"] == 18
    assert one["spacing_summary.csv"]["rows"] == 2
    summary = pd.read_csv(tmp_path / "parameter_summary.csv")
    assert summary.n_planned.eq(3).all()
    assert summary.query("parameter == 'mu'").relative_bias.isna().all()
    assert summary.query("parameter == 'mu'").relative_rmse.isna().all()


def test_missing_entire_cell_cannot_be_silently_omitted(tmp_path):
    design = small_design()
    run_raw(design, tmp_path, progress=False)
    directory = tmp_path / "cells" / design.cells[0].identity
    directory.rename(tmp_path / "removed_cell")
    with pytest.raises(ValueError, match="planned cell set"):
        write_tables(tmp_path, design)
