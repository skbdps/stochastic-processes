"""Run tests and regenerate the notebook, CSVs and figure with this interpreter.

The notebook uses synchronous Python cells. Execute them in an in-process
IPython shell, capturing text and rich outputs without a separate kernel server.
"""
import os
from pathlib import Path
import subprocess
import sys
import time

# Small likelihood evaluations do not benefit from large BLAS thread pools.
# Set before NumPy is imported, and pass the same environment to pytest.
for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"
os.environ["MPLBACKEND"] = "Agg"

import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output
from traitlets.config import Config

ROOT = Path(__file__).resolve().parent


def execute_notebook(path):
    notebook = nbformat.read(path, as_version=4)
    for cell in notebook.cells:
        if cell.cell_type == "code":
            cell.outputs = []
            cell.execution_count = None
    config = Config()
    config.HistoryManager.hist_file = ":memory:"
    shell = InteractiveShell.instance(config=config)
    original_directory = Path.cwd()
    os.chdir(path.parent)
    try:
        for number, cell in enumerate(notebook.cells, 1):
            if cell.cell_type != "code":
                continue
            print(f"Executing notebook cell {number}", flush=True)
            with capture_output() as captured:
                result = shell.run_cell(cell.source, store_history=True)
            cell.execution_count = result.execution_count
            for name, value in (("stdout", captured.stdout), ("stderr", captured.stderr)):
                if value:
                    cell.outputs.append(nbformat.v4.new_output("stream", name=name, text=value))
                    print(value, end="", file=sys.stdout if name == "stdout" else sys.stderr)
            for output in captured.outputs:
                cell.outputs.append(nbformat.v4.new_output(
                    "display_data", data=output.data, metadata=output.metadata))
            if not result.success:
                error = result.error_before_exec or result.error_in_exec
                cell.outputs.append(nbformat.v4.new_output(
                    "error", ename=type(error).__name__, evalue=str(error),
                    traceback=[str(error)]))
                raise RuntimeError(f"Notebook cell {number} failed") from error
        notebook.metadata["language_info"] = {
            "name": "python", "version": sys.version.split()[0],
            "mimetype": "text/x-python", "file_extension": ".py"}
        nbformat.validate(notebook)
    finally:
        # On a failed run, retain the failure and clear later stale outputs.
        nbformat.write(notebook, path)
        os.chdir(original_directory)
        InteractiveShell.clear_instance()


def main():
    started = time.perf_counter()
    print("1/2: unit tests", flush=True)
    subprocess.run([sys.executable, "-m", "pytest", "-q", str(ROOT / "tests")],
                   cwd=ROOT, check=True)
    print("2/2: notebook, replication CSVs and teaser figure", flush=True)
    execute_notebook(ROOT / "notebooks" / "01_exact_mle_demo.ipynb")
    print(f"ALL PASSED - elapsed {time.perf_counter() - started:.1f} s", flush=True)


if __name__ == "__main__":
    main()
