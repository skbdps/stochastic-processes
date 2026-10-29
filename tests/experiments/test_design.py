from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest
import yaml

from ou_irregular.experiments.design import Cell, RunDesign
from ou_irregular.week4.config import load_config as week4_config
from ou_irregular.week5.config import load_config

ROOT = Path(__file__).resolve().parents[2]


def test_registered_design_and_historical_identity():
    design = load_config(ROOT / "configs/week5_e1.yaml")
    assert len(design.cells) == 7
    assert len(design.cells) * design.replications * 3 == 42000
    historical = week4_config(ROOT / "configs/week4_validation.yaml")
    for cell in design.cells:
        assert cell.to_dict() == historical.scientific_cell({"q": .5, "cv": cell.cv, "n": 1000})
    assert design.validation().scope == "validation"
    assert design.validation().seed != design.seed
    assert design.validation().cells == design.cells
    with pytest.raises(FrozenInstanceError):
        design.cells[0].q = 1


@pytest.mark.parametrize("changes", [{"seed": 7}, {"replications": 1000}, {"q": [.1]},
                                    {"n": [250]}, {"floor_fraction": 1e-8}, {"seed": True}])
def test_registration_changes_rejected(tmp_path, changes):
    data = yaml.safe_load((ROOT / "configs/week5_e1.yaml").read_text())
    data.update(changes)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError):
        load_config(path)


@pytest.mark.parametrize("changes", [{"n": True}, {"q": 0}, {"cv": -1}, {"sigma": 0},
                                    {"floor_fraction": 1}, {"mu": float("nan")}])
def test_cell_domains(changes):
    arguments = {"n": 100, "q": .5, "cv": 1}
    arguments.update(changes)
    with pytest.raises(ValueError):
        Cell(**arguments)


def test_missing_duplicate_cells_and_invalid_counts():
    cell = Cell(1000, .5, 0)
    for cells, reps in [((), 2), ((cell, cell), 2), ((cell,), True), ((cell,), 0)]:
        with pytest.raises(ValueError):
            RunDesign("test", 1, reps, cells)
