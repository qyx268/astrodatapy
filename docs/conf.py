# Sphinx configuration for astrodatapy.
#
# Build locally with:
#     pip install ".[docs]"
#     sphinx-build -b html docs docs/_build/html
import configparser
import datetime
import os
import re
from urllib.parse import quote

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
# (shipped with the package and executed by the test suite). They are executed again
# at every docs build, so the docs always show the current code and data, whether or
# not the outputs stored in the notebooks are up to date.

# -- Dataset catalogue -----------------------------------------------------------
# docs/astrodatapy/datasets.rst is written here, at every build, from the info.txt
# files, so that it always lists the data in the repository.
DATA = os.path.join(ROOT, "astrodatapy", "data")
FEATURE_GROUPS = [
    ("Number densities", "astrodatapy.number_density", [
        ("GSMF", "Galaxy stellar mass function -- all"),
        ("GSMF_Blue", "Galaxy stellar mass function -- blue"),
        ("GSMF_Red", "Galaxy stellar mass function -- red"),
        ("GSMF_Quiescent", "Galaxy stellar mass function -- quiescent"),
        ("GSMF_Bulge", "Galaxy stellar mass function -- bulge"),
        ("GSMF_Disk", "Galaxy stellar mass function -- disk"),
        ("SFRF", "Star formation rate function"),
        ("GLF_UV", "Galaxy luminosity function -- UV"),
        ("GLF_OIII", "Galaxy luminosity function -- [OIII] 5008"),
        ("QLF_UV", "Quasar luminosity function -- UV"),
        ("QLF_optical", "Quasar luminosity function -- optical"),
        ("QLF_bolometric", "Quasar luminosity function -- bolometric"),
        ("QLF_Xray", "Quasar luminosity function -- X-ray"),
        ("BHMF", "Black hole mass function"),
        ("BHM", "Black hole mass"),
    ]),
    ("Correlations", "astrodatapy.correlation", [
        ("Magorrian", "Black hole - galaxy bulge mass relation"),
        ("Tully_Fisher", "Mass - velocity of spiral galaxies"),
        ("DiskSize_StellarMass", "Disk size - stellar mass"),
        ("GasFraction_StellarMass", "Gas fraction - stellar mass"),
        ("sSFR_StellarMass_Blue", "sSFR - stellar mass -- blue"),
        ("HaloMass_StellarMass", "Halo mass - stellar mass"),
        ("HaloMass_StellarMass_Blue", "Halo mass - stellar mass -- blue"),
        ("HaloMass_StellarMass_Red", "Halo mass - stellar mass -- red"),
    ]),
    ("Clustering", "astrodatapy.clustering", [
        ("QC_2PTCF", "Quasar two-point correlation function"),
    ]),
    ("Redshift evolution of the IGM", "astrodatapy.evolution", [
        ("IGM_Emissivity", "Ionizing emissivity, returned in s^-1 cMpc^-3"),
        ("IGM_Gamma", "H I photoionization rate, returned in units of 1e-12 s^-1"),
        ("IGM_MFP", "Mean free path of ionizing photons at 912 A, returned in proper Mpc"),
        ("IGM_T0", "IGM temperature at the mean density, returned in K"),
        ("IGM_xHI", "Volume-averaged neutral hydrogen fraction, returned as a fraction (0-1)"),
    ]),
]


def _read_info(feature):
    """Rows of <feature>/info.txt, as dicts keyed by the column names of its #Name line."""
    header, rows = None, []
    with open(os.path.join(DATA, feature, "info.txt")) as f:
        for line in f:
            if line.startswith("#Name"):
                header = line[1:].split()
            elif header and line.strip() and not line.startswith("#"):
                rows.append(dict(zip(header, line.split())))
    return rows


def _redshifts(feature, name):
    """Redshifts of the z*.dat files of a dataset, e.g. z10pt5.dat -> 10.5; for a dataset of
    evolution (one file <name>.dat of rows starting with z), its redshift range and size."""
    directory = os.path.join(DATA, feature, name)
    single = os.path.join(DATA, feature, name + ".dat")
    if os.path.isfile(single):
        zs = [float(l.split()[0]) for l in open(single) if l.strip() and not l.startswith("#")]
        return "%g-%g (%d)" % (min(zs), max(zs), len(zs)) if len(zs) > 1 else "%g" % zs[0]
    if not os.path.isdir(directory):
        return ""
    zs = []
    for fname in os.listdir(directory):
        m = re.match(r"^z(\d+)pt(\d+)\.dat$", fname)
        if m:
            zs.append(float(m.group(1) + "." + m.group(2)))
    return ", ".join("%g" % z for z in sorted(zs))


def _reference(bib):
    """An ADS link for a bibcode, else the text as it is."""
    if re.match(r"^\d{4}\S{15}$", bib):
        return "`%s <https://ui.adsabs.harvard.edu/abs/%s/abstract>`__" % (bib, quote(bib))
    return bib


def write_dataset_catalogue(path):
    listed = set(f for _, _, features in FEATURE_GROUPS for f, _ in features)
    unlisted = sorted(d for d in os.listdir(DATA)
                      if os.path.isfile(os.path.join(DATA, d, "info.txt")) and d not in listed)
    groups = list(FEATURE_GROUPS)
    if unlisted:
        groups.append(("Other", "astrodatapy", [(f, "") for f in unlisted]))
    body, n_total = [], 0
    for title, module, features in groups:
        body += [title, "-" * len(title), "", "Read with ``%s``." % module, ""]
        for feature, description in features:
            rows = _read_info(feature)
            n_total += len(rows)
            heading = feature + (": " + description if description else "")
            columns = ["Dataset", "Redshifts"]
            columns += [c for c in ("Type", "Method", "h", "IMF", "BinWidth") if any(c in r for r in rows)]
            columns += ["Reference"]
            body += [heading, "~" * len(heading), "",
                     ".. list-table::", "   :header-rows: 1", ""]
            body += ["   * - " + columns[0]] + ["     - " + c for c in columns[1:]]
            for r in rows:
                cells = ["``%s``" % r["Name"], _redshifts(feature, r["Name"])]
                cells += [r.get(c, "") for c in columns[2:-1]]
                cells += [_reference(r.get("Bibliographic", ""))]
                body += ["   * - " + cells[0]] + ["     - " + c for c in cells[1:]]
            body.append("")
    head = ["Datasets", "========", "",
            "The %d datasets in the package, with the redshifts of their data files, the"
            % n_total,
            "``h`` of the source paper and the ADS entry of the paper. For the redshift",
            "evolution of the IGM, the redshift range and the number of points are given;",
            "``data_sim`` marks simulation (model) curves, not observations. This page is",
            "generated from the ``info.txt`` files when the documentation is built.", ""]
    text = "\n".join(head + body) + "\n"
    old = open(path).read() if os.path.exists(path) else None
    if text != old:  # do not touch the file (and trigger a rebuild) when nothing changed
        with open(path, "w") as f:
            f.write(text)


write_dataset_catalogue(os.path.join(HERE, "astrodatapy", "datasets.rst"))

# -- General configuration -----------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx_automodapi.automodapi",
    "myst_nb",
]
nb_execution_mode = "force"  # run the notebooks at every build (needs the notebooks' packages)
nb_execution_timeout = 900  # seconds per cell
nb_execution_raise_on_error = True  # a failing notebook fails the build
nb_execution_in_temp = True  # do not write the notebooks' output files into docs/
napoleon_use_ivar = True  # Attributes sections as fields, not as separate entries
autodoc_preserve_defaults = True  # show folder=fileDir+"/data/", not the path of the build machine
automodapi_toctreedirnm = "api"
exclude_patterns = ["_build", "**.ipynb_checkpoints", "_templates"]
templates_path = ["_templates"]
highlight_language = "python3"

# -- HTML output ----------------------------------------------------------------
html_theme = "alabaster"
html_title = "{0} v{1}".format(project, release)
htmlhelp_basename = project + "doc"
