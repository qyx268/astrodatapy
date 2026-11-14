# Sphinx configuration for astrodatapy.
#
# Build locally with:
#     pip install ".[docs]"
#     sphinx-build -b html docs docs/_build/html
import configparser
import datetime
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# -- Project information (from setup.cfg) -------------------------------------
setup_cfg = configparser.ConfigParser()
setup_cfg.read(os.path.join(ROOT, "setup.cfg"))
project = setup_cfg["metadata"]["name"]
author = setup_cfg["metadata"]["author"]
copyright = "{0}, {1}".format(datetime.datetime.now().year, author)

try:
    import astrodatapy
    release = astrodatapy.__version__ or "unknown"
except ImportError:
    release = "unknown"
version = release.split("+")[0]

# The example notebooks in docs/astrodatapy are symlinks to astrodatapy/utils/*.ipynb
# (shipped with the package and executed by the test suite).

# -- General configuration -----------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx_automodapi.automodapi",
    "myst_nb",
]
nb_execution_mode = "off"  # show the stored outputs; the notebooks are executed by the tests
automodapi_toctreedirnm = "api"
exclude_patterns = ["_build", "**.ipynb_checkpoints", "_templates"]
templates_path = ["_templates"]
highlight_language = "python3"

# -- HTML output ----------------------------------------------------------------
html_theme = "alabaster"
html_title = "{0} v{1}".format(project, release)
htmlhelp_basename = project + "doc"
