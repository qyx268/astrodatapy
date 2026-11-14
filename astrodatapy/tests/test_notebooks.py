"""Execute the example notebooks in astrodatapy/utils against the installed package.

These need the [notebooks] extra (matplotlib, scipy, astropy, nbformat, nbclient,
ipykernel) and are skipped when it is not installed:

    pip install "astrodatapy[test,notebooks]"
    pytest --pyargs astrodatapy
"""
import json
import os
import sys

import pytest

nbformat = pytest.importorskip("nbformat")
nbclient = pytest.importorskip("nbclient")
pytest.importorskip("ipykernel")
matplotlib = pytest.importorskip("matplotlib")

UTILS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "utils")
TIMEOUT = 1800  # seconds per cell


def _version_tuple(version):
    return tuple(int(part) for part in version.split(".")[:2] if part.isdigit())


def _execute(name, workdir):
    """Run notebook `name` with a kernel that uses this interpreter, in `workdir`."""
    from jupyter_client.kernelspec import KernelSpecManager
    from jupyter_client.manager import KernelManager

    # A private kernel spec: a "python3" kernel registered elsewhere (e.g. ~/.local)
    # might point to a different interpreter and test the wrong astrodatapy.
    kernel_dir = os.path.join(str(workdir), "kernels", "astrodatapy-test")
    os.makedirs(kernel_dir)
    with open(os.path.join(kernel_dir, "kernel.json"), "w") as f:
        json.dump({"argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
                   "display_name": "astrodatapy test", "language": "python"}, f)
    spec_manager = KernelSpecManager(kernel_dirs=[os.path.dirname(kernel_dir)])
    kernel_manager = KernelManager(kernel_name="astrodatapy-test", kernel_spec_manager=spec_manager)

    nb = nbformat.read(os.path.join(UTILS, name), as_version=4)
    client = nbclient.NotebookClient(nb, km=kernel_manager, timeout=TIMEOUT,
                                     resources={"metadata": {"path": str(workdir)}})
    client.execute()  # raises CellExecutionError if any cell fails
    return nb


def _n_figures(nb):
    return sum("image/png" in out.get("data", {})
               for cell in nb.cells if cell.cell_type == "code" for out in cell.outputs)


def test_plots_notebook(tmp_path):
    nb = _execute("plots.ipynb", tmp_path)
    assert _n_figures(nb) >= 10


@pytest.mark.skipif(sys.version_info < (3, 7), reason="combine.ipynb uses dataclasses (Python >= 3.7)")
@pytest.mark.skipif(_version_tuple(matplotlib.__version__) < (3, 4),
                    reason="combine.ipynb uses Figure.supxlabel (matplotlib >= 3.4)")
def test_combine_notebook(tmp_path):
    pytest.importorskip("scipy")
    pytest.importorskip("astropy")
    nb = _execute("combine.ipynb", tmp_path)
    assert _n_figures(nb) >= 2

    expected = {  # output file: redshifts that must have a combined constraint
        "combined_oiii_lf.json": {"6", "7", "8"},
        "combined_glf_uv.json": {"5", "6", "7", "8", "9", "10"},
        "combined_qlf_uv.json": {"5", "6", "7"},
    }
    for filename, redshifts in expected.items():
        path = tmp_path / "combined_constraints" / filename
        assert path.exists(), filename
        combined = json.loads(path.read_text())
        assert redshifts <= set(combined), (filename, sorted(combined))
        for z, entry in combined.items():
            n = len(entry["x"])
            assert n > 0 and all(len(entry[key]) == n for key in ("y", "yerr_lo", "yerr_hi")), (filename, z)
            values = entry["y"] + entry["yerr_lo"] + entry["yerr_hi"]
            assert all(v == v and abs(v) < 1e3 for v in values), (filename, z)  # finite, no NaN
            assert all(v >= 0 for v in entry["yerr_lo"] + entry["yerr_hi"]), (filename, z)
