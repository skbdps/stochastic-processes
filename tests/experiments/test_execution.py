from dataclasses import replace
import json

import numpy as np
import pandas as pd
import pytest

from ou_irregular.experiments.design import Cell, RunDesign
from ou_irregular.experiments.execution import (
    registered_wald, source_records, run_raw, load_records,
)
from ou_irregular.experiments.storage import CellStore
from ou_irregular.experiments.summaries import ExperimentSummarizer
from ou_irregular.inference.adapters import OuLikelihood, Week3Estimator
from ou_irregular.inference.wald import WaldInference
from ou_irregular.week4.simulation import OuSampleSource


def small_design():
    return RunDesign("engine_validation", 22222222, 3, (Cell(20, .5, 0), Cell(20, .5, 2)), "validation")


def test_tracing_preserves_registered_inference_exactly():
    design = small_design()
    for cell in design.cells:
        sample, _ = OuSampleSource(design.seed).sample(cell.to_dict(), 0)
        for name in ("exact", "pfml", "euler"):
            traced = registered_wald(name)
            a = traced.infer(sample)
            b = WaldInference(Week3Estimator(name), OuLikelihood(name)).infer(sample)
            assert a.valid == b.valid and a.status == b.status
            np.testing.assert_array_equal(a.point_fit.parameters, b.point_fit.parameters)
            np.testing.assert_array_equal(a.covariance, b.covariance)
            assert traced.hessian.last is not None


def test_raw_run_resume_does_not_resimulate_saved_paths(tmp_path, monkeypatch):
    design = small_design()
    manifest = run_raw(design, tmp_path, partition_size=2, progress=False)
    assert manifest["status"] == "raw_complete_pending_audit"
    rows = pd.DataFrame(load_records(tmp_path, design))
    assert len(rows) == 18
    summary = ExperimentSummarizer().summarize(rows, design)
    assert len(summary) == 18
    def forbidden(*args, **kwargs):
        raise AssertionError("completed source must not be resimulated on resume")
    monkeypatch.setattr(OuSampleSource, "sample", forbidden)
    run_raw(design, tmp_path, resume=True, partition_size=1, progress=False)
    pd.testing.assert_frame_equal(rows, pd.DataFrame(load_records(tmp_path, design)))


def test_interrupt_keeps_prefix_and_resumes_same_inputs(tmp_path, monkeypatch):
    design = replace(small_design(), cells=small_design().cells[:1])
    original = CellStore.append
    def interrupt(self, records, stop):
        if self.next_replication > 0:
            raise KeyboardInterrupt("simulated interruption")
        return original(self, records, stop)
    monkeypatch.setattr(CellStore, "append", interrupt)
    with pytest.raises(KeyboardInterrupt):
        run_raw(design, tmp_path, partition_size=1, progress=False)
    state = json.loads((tmp_path / "run_metadata.json").read_text())
    assert state["status"] == "incomplete"
    monkeypatch.setattr(CellStore, "append", original)
    run_raw(design, tmp_path, resume=True, partition_size=1, progress=False)
    clean = tmp_path.parent / (tmp_path.name + "_clean")
    run_raw(design, clean, partition_size=1, progress=False)
    pd.testing.assert_frame_equal(pd.DataFrame(load_records(tmp_path, design)),
                                  pd.DataFrame(load_records(clean, design)))


def test_source_and_inference_exceptions_retain_planned_candidates(tmp_path):
    class FailingSource:
        def sample(self, cell, rep):
            raise ArithmeticError("injected source failure")
    design = small_design()
    state = run_raw(design, tmp_path, source=FailingSource(), progress=False)
    assert state["status"] == "implementation_review_required" and state["pipeline_errors"] == 18
    rows = pd.DataFrame(load_records(tmp_path, design))
    summary = ExperimentSummarizer().summarize(rows, design)
    assert summary.n_planned.eq(3).all() and summary.coverage_operational.eq(0).all()
    class BadInterval:
        estimator = Week3Estimator("exact")
        def infer(self, *args, **kwargs):
            raise RuntimeError("injected inference error")
    services = {name: registered_wald(name) for name in ("exact", "pfml", "euler")}
    services["exact"] = BadInterval()
    records = source_records(design.cells[0], 0, OuSampleSource(design.seed), services)
    assert records[0]["point_valid"] and not records[0]["interval_valid"]
    assert np.isfinite(records[0]["theta"])


def test_headline_requires_frozen_design_and_no_injection(tmp_path):
    design = replace(small_design(), scope="headline")
    with pytest.raises(ValueError, match="frozen entry-point"):
        run_raw(design, tmp_path, progress=False)
