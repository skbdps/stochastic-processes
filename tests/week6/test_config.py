from copy import deepcopy
from itertools import product

import pytest
import yaml

from ou_irregular.experiments.provenance import ROOT
from ou_irregular.week6.config import REGISTERED, load_config


def test_full_registered_cartesian_design_and_separate_validation_stream():
    design = load_config(ROOT / "configs/week6_e2.yaml")
    assert design.seed == 20261002 and design.replications == 2000
    expected = set(product([250, 1000, 4000], [.1, .5, 1., 2.], [0., .5, 1., 2., 3.]))
    assert {(c.n, c.q, c.cv) for c in design.cells} == expected
    assert len(design.cells) == len({c.identity for c in design.cells}) == 60
    assert len(design.cells)*design.replications == 120000
    shared = next(c for c in design.cells if (c.n, c.q, c.cv) == (1000, .5, 3.))
    assert shared.identity == "b75093677ab1ae4fcca15bdf6b3992b2249ec307bba67b83042a01710a530e9e"
    validation = design.validation(2)
    assert validation.seed == 21261002 and validation.scope == "validation"
    assert validation.cells == design.cells and validation.replications == 2


@pytest.mark.parametrize("field,value", [
    ("seed",20261001), ("replications",10000), ("q",[.1,.5,1.]),
    ("cv",[0.,.5,1.,2.]), ("n",[250,1000]), ("floor_fraction",1e-8),
    ("floor_cv_threshold",2.), ("experiment","E1"),
    ("truth",{"theta":1.,"mu":0.,"sigma":1.}), ("extra",True),
])
def test_config_cannot_shrink_or_amend_headline(tmp_path,field,value):
    data=deepcopy(REGISTERED);data[field]=value
    path=tmp_path/"altered.yaml";path.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError,match="frozen registration"):
        load_config(path)


def test_duplicate_configuration_keys_rejected(tmp_path):
    path=tmp_path/"duplicate.yaml"
    path.write_text((ROOT/"configs/week6_e2.yaml").read_text()+"seed: 20261002\n")
    with pytest.raises(ValueError):
        load_config(path)
