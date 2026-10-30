"""Author checks of independent-audit boundaries; no source generation here."""
from copy import deepcopy
from dataclasses import replace

import pytest

from ou_irregular.week7 import floor_audit
from ou_irregular.week7.selection import floor_design


@pytest.mark.parametrize("payload", [
    '{"x":1,"x":2}', '{"diagnostic":{"x":1,"x":2}}',
    '{"x":NaN}', '{"x":Infinity}', '{"x":-Infinity}',
    '{"x":1e999}', '{"x":[0,{"y":-1e999}]}',
])
def test_strict_json_rejects_ambiguous_or_nonfinite_evidence(payload):
    with pytest.raises(ValueError):
        floor_audit.strict_json(payload)


def test_nested_diagnostic_strings_are_strict_before_legacy_parsing():
    record = {"point_diagnostics_json": '{"inner_json":"{\\"x\\":1,\\"x\\":2}"}'}
    with pytest.raises(ValueError, match="duplicate"):
        floor_audit._nested_json(record)
    with pytest.raises(ValueError, match="nonfinite"):
        floor_audit._nested_json({"interval_diagnostics_json": '{"hessian_diagnostics":{"relative_difference":1e999}}'})


@pytest.mark.parametrize("change", ["seed", "cell_order", "scope", "name", "one_floor"])
def test_fixed_design_requires_exact_namespace_order_and_complete_floor_groups(change):
    design = floor_design().validation(2)
    if change == "seed":
        changed = replace(design, seed=design.seed+1)
    elif change == "cell_order":
        changed = replace(design, cells=tuple(reversed(design.cells)))
    elif change == "scope":
        changed = replace(design, scope="exploratory")
    elif change == "name":
        changed = replace(design, name="unregistered_floor_name")
    else:
        changed = replace(design, cells=(design.cells[0],))
    with pytest.raises(ValueError):
        floor_audit._require_design(changed)
    if change == "one_floor":
        with pytest.raises(ValueError):
            floor_audit._require_design(changed, subset=True)


def test_fixed_complete_cv_subdesign_is_allowed_only_for_pure_epsilon_arithmetic():
    design = floor_design().validation(2)
    group = replace(design, cells=tuple(c for c in design.cells if c.cv == 2.))
    assert floor_audit._require_design(group, subset=True) is group
    with pytest.raises(ValueError):
        floor_audit._require_design(group)


def test_complete_replay_comparison_preserves_missingness_types_and_nested_keys():
    expected = {"rep": 1, "valid": True, "value": 2., "missing": None,
                "point_diagnostics_json": '{"gradient":[1.0,2.0],"status":"valid"}'}
    floor_audit._compare(deepcopy(expected), expected, "fixture")
    tolerated = deepcopy(expected)
    tolerated["value"] += 1e-11
    floor_audit._compare(tolerated, expected, "fixture")
    with pytest.raises(ValueError):
        floor_audit._compare(tolerated, expected, "source fixture", exact=True)
    mutations = [{**expected, "rep": True}, {**expected, "value": 2},
                 {**expected, "missing": 0.}, {**expected, "unexpected": 0},
                 {**expected, "point_diagnostics_json": '{"gradient":[1.0,2.0],"status":"valid","extra":0}'}]
    for changed in mutations:
        with pytest.raises(ValueError):
            floor_audit._compare(changed, expected, "fixture")


@pytest.mark.parametrize("replay", [1, "true", None])
def test_invalid_replay_option_is_rejected_before_any_output_creation(tmp_path, replay):
    output = tmp_path/"never_created"
    with pytest.raises(ValueError, match="Boolean"):
        floor_audit.audit_run(output, floor_design().validation(2), replay=replay)
    assert not output.exists()


def test_missing_archive_is_not_created_or_regenerated(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("an absent archive must not generate a source")
    monkeypatch.setattr(floor_audit, "independent_source", forbidden)
    output = tmp_path/"absent_archive"
    with pytest.raises(ValueError, match="existing directory"):
        floor_audit.audit_run(output, floor_design().validation(2))
    assert list(tmp_path.iterdir()) == []
