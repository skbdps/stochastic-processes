import json

import pytest

from ou_irregular.experiments.provenance import ROOT
from ou_irregular.week5.config import load_config as load_e1
from ou_irregular.week6.runner import _saved_design, main


def test_headline_replication_override_is_rejected_before_output(tmp_path):
    target=tmp_path/"unregistered"
    with pytest.raises(SystemExit):
        main(["run", "--output", str(target), "--replications", "2"])
    assert not target.exists()


def test_invalid_validation_replications_are_rejected_before_output(tmp_path):
    target=tmp_path/"invalid"
    with pytest.raises(ValueError):
        main(["run", "--validation", "--replications", "0", "--output", str(target)])
    assert not target.exists()


def test_e1_archive_cannot_be_claimed_as_e2(tmp_path):
    design=load_e1(ROOT/"configs/week5_e1.yaml")
    (tmp_path/"resolved_design.json").write_text(json.dumps(design.to_dict()))
    with pytest.raises(ValueError,match="frozen E2"):
        _saved_design(tmp_path)
