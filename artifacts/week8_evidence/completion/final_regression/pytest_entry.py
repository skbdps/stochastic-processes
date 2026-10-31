"""Spawn-safe release test entry; no scientific/test monkeypatches."""
import json, os
from pathlib import Path
import sys
sys.path.insert(0, str(Path.cwd()))
def main():
    import ou_irregular.runner
    import pytest
    os.environ.pop("MKL_NUM_THREADS", None)
    os.environ.pop("VECLIB_MAXIMUM_THREADS", None)
    from ou_irregular.experiments.provenance import environment
    expected = json.loads(Path(os.environ["OU_RELEASE_RUNTIME"]).read_text())
    if environment() != expected:
        raise RuntimeError("test-entry runtime differs from recorded current runtime")
    return pytest.main(sys.argv[1:])
if __name__ == "__main__":
    raise SystemExit(main())
