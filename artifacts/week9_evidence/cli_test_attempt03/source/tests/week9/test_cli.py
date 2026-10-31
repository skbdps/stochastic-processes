"""A real saved-evidence build and negative integrity/lifecycle checks.

No positive visual-review record is fabricated here. Review-shaped payloads
are explicitly blocked, uninspected test data and can never complete a seal.
Only the temporary derivative output is modified by corruption tests.
"""
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
import builtins
import importlib
import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit

import pytest

from ou_irregular.week9 import cli
from ou_irregular.week9.inputs import load_sources


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT/"configs/week9_publication.json"


def forbidden(*args, **kwargs):
    raise AssertionError("Week 9 publication must not execute simulation, model fitting, or numerical source replay")


def forbid_scientific_execution(patch):
    targets = {
        "ou_irregular.ou": ("simulate_ou", "fit"),
        "ou_irregular.estimators": ("fit_estimator",),
        "ou_irregular.week8.source": ("generate_gaps",),
        "ou_irregular.week8.audit": ("reconstruct_schedule", "audit_run"),
        "ou_irregular.experiments.audit": ("audit_run",),
    }
    for module_name, names in targets.items():
        module = sys.modules.get(module_name)
        if module is None:
            continue
        for name in names:
            patch.setattr(module, name, forbidden)
    for module_name, class_name in (("ou_irregular.week4.simulation", "OuSampleSource"),
                                   ("ou_irregular.week7.floor_source", "PairedFloorSource"),
                                   ("ou_irregular.week8.source", "E4Source")):
        module = sys.modules.get(module_name)
        if module is not None:
            patch.setattr(getattr(module, class_name), "sample", forbidden)
    original_import = builtins.__import__
    scientific = (*targets, "ou_irregular.spacing", "ou_irregular.inference",
                  "ou_irregular.week4.simulation", "ou_irregular.week7.floor_source",
                  "ou_irregular.week7.execution", "ou_irregular.week7.path_diagnostics",
                  "ou_irregular.week8.execution", "ou_irregular.week5.runner",
                  "ou_irregular.week6.runner", "ou_irregular.week7.runner", "ou_irregular.week8.runner")

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        absolute = importlib.util.resolve_name("."*level+name, globals["__package__"]) if level else name
        if any(absolute == prefix or absolute.startswith(prefix+".") for prefix in scientific):
            forbidden()
        return original_import(name, globals, locals, fromlist, level)

    # Prevent the publication layer from loading scientific execution modules
    # at all. Already-loaded entry points above are also guarded. This keeps
    # the saved-data test independent of historical-only dependencies.
    patch.setattr(builtins, "__import__", guarded_import)


@pytest.fixture(scope="module")
def publication(tmp_path_factory):
    """One real build; original scientific inputs remain byte-identical."""
    with pytest.MonkeyPatch.context() as patch:
        forbid_scientific_execution(patch)
        bundle = load_sources(CONFIG, root=ROOT)
        before = bundle.snapshot()
        output = tmp_path_factory.mktemp("week9_cli")/"pending_package"
        result = cli.build(CONFIG, output, root=ROOT)
        assert result["verified"] and result["status"] == "pending_visual_review"
        assert result["new_simulated_sources"] == result["new_model_fits"] == 0
        assert bundle.snapshot() == before
        yield output
        assert bundle.snapshot() == before
        assert cli.verify(output, root=ROOT)["status"] == "pending_visual_review"
        assert not (output/"visual_review.json").exists()


@contextmanager
def mutated_member(output, relative, payload, *, rebind=False):
    """Temporarily corrupt one derivative file, optionally rebinding its hashes.

    Rebinding demonstrates that deterministic source reconstruction of the
    derivative view checks more than its user-editable outer manifest.
    """
    member = output/relative
    names = [member, output/"manifest.json", output/"package_state.json"]
    originals = {path: path.read_bytes() for path in names}
    try:
        member.write_bytes(payload)
        if rebind:
            manifest = json.loads(originals[output/"manifest.json"])
            manifest["files"][relative] = cli.sha(payload)
            (output/"manifest.json").write_bytes(cli.json_bytes(manifest))
            state = json.loads(originals[output/"package_state.json"])
            state["manifest_sha256"] = cli.digest(output/"manifest.json")
            if relative == "figure_manifest.json":
                state["figure_manifest_sha256"] = cli.digest(output/"figure_manifest.json")
            (output/"package_state.json").write_bytes(cli.json_bytes(state))
        yield
    finally:
        for path, original in originals.items():
            path.write_bytes(original)


def blocked_uninspected_review(output):
    """Deliberately invalid payload; zero hashes are test placeholders only."""
    figures = cli.read_json(output/"figure_manifest.json")
    pdfs = [f for f in figures["figures"] if f["kind"] == "new_week9_layout" and f["path"].endswith(".pdf")]
    return {
        "test_payload_only": True,
        "reviewer": "TEST ONLY: no visual inspection performed",
        "review_type": "actual_ai_raster_and_pdf_review",
        "reviewed_utc": "2000-01-01T00:00:00+00:00",
        "figure_manifest_sha256": cli.digest(output/"figure_manifest.json"),
        "reviewed_png_sha256": {f["path"]: f["sha256"] for f in figures["figures"]
                                 if f["kind"] == "new_week9_layout" and f["path"].endswith(".png")},
        "reviewed_pdf_pages": {f["path"]: {"page": 1, "rendered_page_sha256": "0"*64,
            "source_pdf_sha256": f["sha256"], "render_command": "TEST PLACEHOLDER: no render run",
            "renderer": "TEST PLACEHOLDER: no renderer executed"} for f in pdfs},
        "all_panels_legible": False,
        "blocking_issues": ["Deliberately invalid uninspected test payload; completion forbidden"],
        "notes": ["No image was inspected by this test. This payload must never seal a package."],
    }


def test_real_build_generates_bound_tables_and_visible_qualifications(publication):
    record = cli.read_json(publication/"build_record.json")
    assert record["source_preservation_passed"]
    assert record["new_simulated_sources"] == record["new_model_fits"] == 0
    assert record["tables"]["T1_e2.csv"]["rows"] == 360
    assert record["tables"]["T1_e1.csv"]["rows"] == 42
    assert record["tables"]["unresolved.csv"]["rows"] == 6
    assert record["tables"]["f1_ecdf.csv"]["rows"] == 126000
    assert (publication/"tables/f1_ecdf.csv.gz").is_file()
    assert (publication/"source_catalog.json").read_bytes() == CONFIG.read_bytes()
    claims = cli.read_json(publication/"claims.json")
    readable = (publication/"claims.md").read_text()
    for claim in claims:
        assert claim["text"] in readable
        assert claim["qualifications"] in readable
        assert claim["status"] in readable
    assert (publication/"claims.csv").read_bytes() == cli.claims_csv_bytes(claims)
    captions = (publication/"captions.md").read_text()
    assert "pointwise" in captions and "realized-span" in captions
    assert "no continuous boundary" in captions.lower() or "continuous boundary" in captions.lower()


def test_pending_package_cannot_be_reported_completed(publication):
    with pytest.raises(ValueError, match="awaiting actual visual review"):
        cli.verify(publication, root=ROOT, require_completed=True)


def test_readable_t1_preserves_scope_and_all_unresolved_cases(publication):
    text = (publication/"tables/T1.md").read_text()
    assert "W does not certify nominal 95% coverage" in text
    assert "No continuous or maximum-safe-CV boundary" in text
    assert "Mu's zero-truth relative-bias criterion is inapplicable" in text
    assert "All six unresolved E1/E2 classifications" in text
    unresolved = cli.pd.read_csv(publication/"tables/unresolved.csv", float_precision="round_trip")
    assert len(unresolved) == 6
    for row in unresolved.itertuples():
        assert f"| {row.study} | {int(row.n)} | {row.theta_mean_gap:g} | {row.cv:g} | {row.estimator} | {row.parameter} |" in text


def test_copied_notes_have_portable_resolving_links(publication):
    links = 0
    for note in sorted((publication/"notes").glob("*.md")):
        text = note.read_text()
        for match in re.finditer(r"\[[^\]]*\]\((?:<([^>]+)>|([^\s)]+))\)", text):
            target = match.group(1) or match.group(2)
            links += 1
            parsed = urlsplit(target)
            assert parsed.scheme in ("", "http", "https")
            assert "/Users/" not in target and "file://" not in target
            if parsed.scheme:
                assert parsed.netloc
                continue
            if not parsed.path:
                continue
            path = (note.parent/unquote(parsed.path)).resolve()
            assert path.is_relative_to(publication.resolve()), (note.name, target)
            assert path.is_file(), (note.name, target)
    assert links > 0, "copied publication notes must retain their evidence links"


def test_changed_member_is_rejected_without_rebinding(publication):
    path = "claims.json"
    with mutated_member(publication, path, (publication/path).read_bytes()+b" "):
        with pytest.raises(ValueError, match="publication file changed"):
            cli.verify(publication, root=ROOT)


@pytest.mark.parametrize("relative,expected", [
    ("tables/T1.md", "readable T1 differs"),
    ("notes/week9_protocol.md", "publication note/link bindings changed"),
    ("claims.csv", "claim CSV differs"),
    ("claims.md", "readable claims differ"),
    ("captions.md", "captions differ"),
])
def test_rebound_human_views_are_regenerated_from_sources(publication, relative, expected):
    with mutated_member(publication, relative, (publication/relative).read_bytes()+b"\nFORGED DERIVATIVE TEXT\n", rebind=True):
        with pytest.raises(ValueError, match=expected):
            cli.verify(publication, root=ROOT)


def test_rebound_source_catalog_cannot_redefine_accepted_inputs(publication):
    catalog = cli.read_json(publication/"source_catalog.json")
    catalog["baseline_commit"] = "0"*40
    with mutated_member(publication, "source_catalog.json", cli.json_bytes(catalog), rebind=True):
        with pytest.raises(ValueError, match="catalog digest|baseline"):
            cli.verify(publication, root=ROOT)


@pytest.mark.parametrize("kind", ["extra_file", "symlink"])
def test_compare_checks_actual_membership_even_when_manifests_match(publication, tmp_path, kind):
    unexpected = publication/"unexpected_member.txt"
    if kind == "extra_file":
        unexpected.write_text("not manifested")
    else:
        target = tmp_path/"harmless_target.txt"
        target.write_text("symlinks must be rejected")
        unexpected.symlink_to(target)
    try:
        with pytest.raises(ValueError, match="membership changed|symlink"):
            cli.compare(publication, publication, root=ROOT)
    finally:
        unexpected.unlink()


@pytest.mark.parametrize("mutation,expected", [
    ("missing_png", "omits or changes a new PNG"),
    ("missing_pdf", "omits a new PDF"),
    ("forged_figure_manifest", "different figures"),
    ("nonhex_render_hash", "PDF page review binding"),
    ("forged_pdf_source", "PDF page review binding"),
    ("future_timestamp", "timestamp"),
    ("blocking_issues", "blocking issues"),
])
def test_invalid_review_payloads_cannot_certify_uninspected_figures(publication, mutation, expected):
    review = blocked_uninspected_review(publication)
    if mutation == "missing_png":
        review["reviewed_png_sha256"].pop(next(iter(review["reviewed_png_sha256"])))
    elif mutation == "missing_pdf":
        review["reviewed_pdf_pages"].pop(next(iter(review["reviewed_pdf_pages"])))
    elif mutation == "forged_figure_manifest":
        review["figure_manifest_sha256"] = "f"*64
    elif mutation == "nonhex_render_hash":
        next(iter(review["reviewed_pdf_pages"].values()))["rendered_page_sha256"] = "not a hexadecimal SHA-256"
    elif mutation == "forged_pdf_source":
        next(iter(review["reviewed_pdf_pages"].values()))["source_pdf_sha256"] = "f"*64
    elif mutation == "future_timestamp":
        review["reviewed_utc"] = (datetime.now(timezone.utc)+timedelta(days=1)).isoformat()
    with pytest.raises(ValueError, match=expected):
        cli.validate_review(publication, review)


def test_finalize_rejects_invalid_review_without_lifecycle_mutation(publication, tmp_path):
    review = blocked_uninspected_review(publication)
    review["reviewed_png_sha256"] = {}
    path = tmp_path/"deliberately_invalid_uninspected_review.json"
    path.write_bytes(cli.json_bytes(review))
    before = (publication/"package_state.json").read_bytes()
    with pytest.raises(ValueError, match="omits or changes a new PNG"):
        cli.finalize(publication, path, root=ROOT)
    assert (publication/"package_state.json").read_bytes() == before
    assert not (publication/"visual_review.json").exists()
    assert not (publication/".package_state.tmp").exists()
