"""Record test identities/outcomes without changing pytest behavior."""
import json, os
from pathlib import Path
_reports = []
def pytest_collection_finish(session):
    Path(os.environ["OU_RELEASE_COLLECTION_OUTPUT"]).write_text(
        json.dumps([item.nodeid for item in session.items], indent=2) + "\n")
def pytest_runtest_logreport(report):
    _reports.append({"nodeid": report.nodeid, "when": report.when,
                     "outcome": report.outcome})
def pytest_sessionfinish(session, exitstatus):
    Path(os.environ["OU_RELEASE_REPORT_OUTPUT"]).write_text(
        json.dumps({"exitstatus": int(exitstatus), "reports": _reports},
                   indent=2) + "\n")
