# Contributing to astrodatapy

Thank you for helping. New datasets, corrections to existing ones, bug reports and code
improvements are all welcome.

## Reporting a problem or asking for help

Open an issue at <https://github.com/qyx268/astrodatapy/issues>. For a data problem, name
the feature and dataset (e.g. `GSMF` / `Baldry2012`) and, if you can, the table or figure of
the paper the value should match. For a bug, include the code that shows it, the error, and
your astrodatapy, Python, numpy and pandas versions. Questions about using the package are
welcome as issues too.

## Adding or correcting data

Each feature is a directory in `astrodatapy/data/` with an `info.txt` that describes its
datasets (data format, log or linear, normalisation, Hubble parameter, IMF, type, method,
bibcode) and one file per redshift (or one file, for the redshift-evolution features) per
dataset. The [data format](astrodatapy/data/README.rst) explains every column, the units,
the h and IMF conversions, and how upper and lower limits are written.

To add a dataset:

1. Check that it is not already there, perhaps under its arXiv name: datasets are named
   first author + publication year (e.g. `Harikane2024`), with a suffix for several models
   of one paper (`Qin2026_Fiducial`).
2. Keep the values in the units and conventions of the paper, and describe them in
   `info.txt`; astrodatapy converts them when they are loaded.
3. Start each data file with a header that cites the paper, the table or figure the values
   come from, and anything unusual (error type, h convention, values read off a figure).
4. Add the paper's bibcode to `info.txt` and a BibTeX entry, keyed by the bibcode, to
   `astrodatapy/data/bibliography.bib`.
5. Run the tests (below); they load every dataset of every feature.

Simulation (model) curves are welcome as well: mark them `data_sim` in `info.txt` and say in
the header how the run was calibrated.

## Code changes

Work on a branch and open a pull request against `master`. Please

* keep `python_requires >= 3.6` working and the sources ASCII;
* add or update tests in `astrodatapy/tests/` for the behaviour you change;
* update the docstrings, the README or `astrodatapy/data/README.rst` when the behaviour or
  the data format changes.

## Running the tests and building the docs

```console
$ pip install -e ".[test,notebooks]"
$ pytest --pyargs astrodatapy         # loads every dataset and executes the example notebooks
$ pip install -e ".[docs]"
$ sphinx-build -W -b html docs docs/_build/html
```

The continuous integration runs the tests on Python 3.7 to 3.14 and builds the docs for
every pull request.
