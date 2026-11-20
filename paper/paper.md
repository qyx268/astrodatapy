---
title: 'astrodatapy: a Python compilation of observational constraints for galaxy formation and reionization models'
tags:
  - Python
  - astronomy
  - galaxy formation
  - luminosity functions
  - reionization
  - intergalactic medium
authors:
  - name: Yuxiang Qin
    orcid: 0000-0000-0000-0000   # CONFIRM: your ORCID
    affiliation: 1
affiliations:
  - name: The Australian National University, Australia   # CONFIRM: school/centre, and any co-authors
    index: 1
date: 20 November 2026
bibliography: paper.bib
---

# Summary

Models of galaxy formation and cosmic reionization are calibrated and tested against
statistical measurements of the galaxy and quasar populations and of the intergalactic
medium (IGM): how many galaxies there are per unit volume at a given stellar mass or
luminosity, how their properties correlate, and how the ionization state and temperature of
the IGM evolve with redshift. `astrodatapy` is a Python package that collects these
measurements in one place. Version 2.3.0 holds 307 datasets from 187 papers in 29 quantities:
stellar mass, star formation rate, ultraviolet (UV) and [O III] luminosity functions of
galaxies; quasar luminosity functions in the UV, optical, X-ray and bolometric luminosity and
black hole mass functions; scaling relations such as the black hole–bulge, Tully–Fisher,
size–mass, gas fraction–mass and stellar-to-halo mass relations; quasar clustering; and the
evolution of the IGM ionizing emissivity, H I photoionization rate, mean free path of ionizing
photons, temperature and neutral fraction. A single call loads every dataset of a quantity
near a target redshift and returns it converted to the Hubble parameter, initial mass function
(IMF) and units the user asks for (\autoref{fig:example}).

# Statement of need

The measurements a model is compared with are published in heterogeneous forms. Number
densities are given per Mpc$^3$ or per $h^{-3}$ Mpc$^3$ for different values of the Hubble
parameter $h$; stellar masses and star formation rates assume a @Salpeter1955, @Kroupa2001 or
@Chabrier2003 IMF; luminosity functions use different magnitude systems and bands; quantities
are tabulated in linear or logarithmic form with symmetric, asymmetric or logarithmic errors;
some points are upper or lower limits; lengths are comoving or proper; some results are only
given as fitted Schechter or power-law parameters, and some only as figures. Every modelling
group re-collects and re-converts these data, conversions are a common source of silent
errors, and the exact data behind a published comparison are rarely recorded.

`astrodatapy` stores each dataset once, in the units and conventions of its paper, with those
conventions recorded as metadata, and performs the conversions when the data are loaded. It
is aimed at modellers who calibrate semi-analytic, semi-numerical or hydrodynamical models of
galaxy formation and reionization and at observers who place a new measurement among
previous ones. Because every dataset carries a bibcode and a header that cites the table or
figure its values come from, a comparison made with a given `astrodatapy` version can be
reproduced and traced back to the literature.

# State of the field

Several simulation codes ship the observational data they are compared with: the
semi-analytic models Shark [@Lagos2018], Galacticus [@Benson2012] and L-Galaxies
[@Henriques2015], the empirical model UniverseMachine [@Behroozi2019], and the analysis
pipeline of the SWIFT code [@Schaller2024] (its `velociraptor-comparison-data` repository).
These collections serve one code: the data are converted to that code's assumptions in
advance, or the corrections for $h$ and the IMF are written into the script that reads each
dataset, and they are not distributed as libraries for other codes. Inference frameworks such
as 21CMMC [@Greig2015; @Park2019] bundle the few datasets their likelihoods use.
Compilations of a single observable collect one kind of measurement, such as the quasar
luminosity functions homogenised by @Shen2020 or the published upper limits on the 21-cm power
spectrum from reionization and cosmic dawn in `eor_limits` [@eorlimits], which provides them
as human-readable files with plotting tools. The closest general-purpose tool is CoReCon [@Garaldi2023], a pip-installable
collection of reionization constraints (ionized fraction, photoionization rate, mean free
path, IGM temperature, UV luminosity functions and others) that records references and
limits. CoReCon stores each quantity in one fixed unit, does not convert between Hubble
parameters or IMFs, and does not cover stellar or black hole mass functions or scaling
relations. Archives such as VizieR [@Ochsenbein2000], accessed with astroquery
[@Ginsburg2019], serve published catalogues and tables rather than curated, homogenised
constraints.

`astrodatapy` has been public since 2018, five years before CoReCon. Its scope (galaxy,
black hole and quasar statistics, scaling relations and IGM evolution in one package) and its
central design choice (store published values with their conventions and convert on request)
do not fit the fixed-unit structure of an existing collection, so it is developed as a
separate package.

# Software design

The package separates data from code. Each quantity ("feature") is a directory of plain-text
files: an `info.txt` with one row per dataset (data format, logarithmic or linear values,
normalisation, Hubble parameter and its power, IMF, comoving or proper units, dataset type,
measurement method and bibcode), and one data file per dataset and redshift with a header
that documents the source. Plain text was chosen over a binary or database format so that
every value can be read, diffed and reviewed in a pull request, and so that contributors can
add data without writing code. Datasets keep their published values rather than values
pre-converted to one convention: the provenance stays checkable against the paper, and any
target convention can be requested.

Four loader classes cover the shapes of the data: `number_density` for distribution functions
at a given redshift, `correlation` for scaling relations, `clustering` for correlation lengths
and `evolution` for quantities tabulated against redshift. Each returns a `pandas` [@McKinney2010]
table of datasets whose data are `numpy` [@Harris2020] arrays of values with upper and lower
bounds, together with flags for upper and lower limits, the dataset type and the method. The
conversions are applied in a fixed order (in `number_density`: evaluation of parametric fits,
IMF, magnitude band, error format, limits, logarithm, normalisation and Hubble parameter), so
a new dataset needs only a row of metadata. Upper limits are ordinary rows whose errors are written as
not-a-number; they survive every conversion and are flagged at the end. Simulation curves are
stored in the same format and marked `data_sim`, so they can be drawn next to the
measurements without being mistaken for them.

The test suite loads every dataset of every feature, checks unit conversions against
hand-computed values, and executes the example notebooks; it runs on Python 3.7 to 3.14 in
continuous integration, and the documentation, built on every change, executes the
notebooks and generates a catalogue of all datasets from the `info.txt` files. One notebook
combines the high-redshift galaxy and quasar UV and [O III] luminosity functions into binned
constraints, archived on Zenodo [@Qin2026zenodo].

![Four of the 29 quantities in `astrodatapy` 2.3.0, as loaded with $h = 0.6751$: (a) the
galaxy stellar mass function at $z = 0$, converted to a @Chabrier2003 IMF; (b) the galaxy UV
luminosity function at $z = 7$; (c) the quasar UV luminosity function at $z = 6$; (d) the
volume-averaged neutral hydrogen fraction of the IGM against redshift, with simulation curves
in grey. Arrows mark upper and lower limits. Produced by `paper/make_figure.py` with Matplotlib
[@Hunter2007].\label{fig:example}](figure.png)

# Research impact statement

`astrodatapy` has been used since 2018 to calibrate and validate the Meraxes semi-analytic
model of the DRAGONS programme [@Mutch2016]. The Meraxes analysis package `dragons` requires
it for its comparisons with the observed stellar mass, star formation rate, UV luminosity,
quasar luminosity and black hole mass functions, and the Meraxes calibration driver `mhysa`
lists it as a dependency. Published studies that used it include the Meraxes work on the
sizes and morphologies of high-redshift galaxies [@Marshall2019] and on the evolution of
black hole scaling relations [@Marshall2020], whose public analysis scripts load the
comparison data with `astrodatapy`, and the study of $z \gtrsim 12$ galaxies observed by
JWST of @Qin2023, which used it to convert the observations to a common IMF. Outside this
collaboration, the stellar mass function plotting tool of the MP-Gadget simulation code has
used it since 2019, and the Markov chain Monte Carlo calibration framework of the SAGE
semi-analytic model (`sage-mcmc`) takes its observational constraints from it. The
high-redshift UV and [O III] luminosity-function constraints combined with the package are
archived on Zenodo [@Qin2026zenodo]. Five contributors besides the author have added data
or code through pull requests, and issues from users have shaped the package; for example,
the per-point bin widths it has returned since version 2.1.0 were requested for likelihood
analyses.

# AI usage disclosure

<!-- CONFIRM (author): that you reviewed every change, and that nothing before November 2026 used generative AI
     (no other AI tools, e.g. Copilot or chat assistants). JOSS treats an inaccurate disclosure as an ethics issue. -->

The work of November 2026 that led to versions 1.0.0 to 2.3.0 was done with the assistance
of Claude Code (Anthropic), using the Claude Opus 5.5 model; the 74 commits it co-authored
carry a `Co-Authored-By` trailer. It assisted with porting the package to Python 3 and modern
packaging, writing tests, continuous integration and documentation, auditing existing
datasets against their source papers, transcribing tables and reading curves from published
figures for new datasets, and drafting this paper. The author reviewed every change; from
version 1.0.0 onwards each was proposed as a pull request and merged by the author. New and
corrected values were checked against the source papers, and the scripts and checks used to
read values from figures are kept with the data provenance. The design of the package and of
its data format, and the choice of datasets, are the author's. The code and data written
between 2018 and July 2026 were produced without generative AI.

# Acknowledgements

We thank Balu Sreedhar, Fatma Shaban, Madeline Marshall, Simon Mutch and Varun Nikam for
their contributions to the package, and Manodeep Sinha and Emanuele Maria Ventura for issues
and discussions that improved it.
<!-- CONFIRM (author): funding (grants, centres such as ASTRO 3D) and computing facilities to acknowledge, and
     whether any sponsor had a role; JOSS requires financial support to be acknowledged. If any of the people
     above become co-authors, remove them here. -->

# References
