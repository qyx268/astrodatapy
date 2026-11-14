Data format
===========

This directory holds the data that astrodatapy loads. Each feature (``GLF_UV``,
``GSMF``, ``Magorrian``, ...) has its own directory::

    <feature>/info.txt                  one line per dataset
    <feature>/<dataset>/z5pt0.dat       one file per redshift
    <feature>/<dataset>_ulimit/z5pt0.dat  upper limits from the same paper, if any

Only the datasets listed in ``info.txt`` are loaded. The redshift is read from the
file name (``z5pt0.dat`` is z = 5.0, ``z10pt5.dat`` is z = 10.5), so name each file
after the centre of the paper's redshift bin. A dataset is named after its first
author and year of publication, e.g. ``Harikane2024``; when an arXiv paper is
published in a later year, rename it.

The data files are read with ``numpy.loadtxt``, so lines starting with ``#`` are
comments. Start every new file with a header giving the paper, the table or figure
the values come from, the redshift range, the units and anything left out.

The BibTeX entries of all papers are collected in ``bibliography.bib``, keyed by
their ADS bibcode.

Number densities
----------------

These are the features read by ``astrodatapy.number_density``: ``GSMF*``, SFRF,
BHMF, BHM, GLF_UV, GLF_OIII and ``QLF_*``. The columns of their ``info.txt`` are separated by
whitespace:

``#Name``
    Name of the dataset, i.e. of the directory holding its files.
``MAG``
    (GLF and QLF) What x is. ``M1450``, ``M1500``, ``M1600``, ``MABB`` (AB magnitude
    in the B band), ``OIII`` (log10 L[OIII] in erg/s), ``LX_2_10keV`` (log10 LX in
    erg/s) and ``Bolometric`` are used as they are; ``MB2MABB``, ``Mg2MB``,
    ``M14502MB`` and ``LBsol2M1450`` convert between bands.
``IMF``
    (``GSMF*`` and SFRF) IMF of the stellar masses or star formation rates:
    ``Salpeter``, ``DietSalpeter``, ``Kroupa`` or ``Chabrier``. The data are
    converted to the ``IMF_out`` requested by the user.
``h``
    H0/(100 km/s/Mpc) of the cosmology the paper used.
``NormalizationOfPhi``
    Factor the tabulated phi is multiplied by, e.g. ``1e-5`` for a table in units of
    10^-5 Mpc^-3 mag^-1. It must be 1 when ``Log`` is 1.
``Errors``
    Layout of the error columns, see below.
``Log``
    1 if phi and its errors are log10 values (errors in dex), 0 if they are linear.
``HinPhi``
    1 if phi is in h^3 Mpc^-3, 0 if it is in Mpc^-3.
``Type``
    ``data``, ``dataULimit`` (upper limits), ``data_sim`` (a model, drawn as a line
    with a band), or a parametric fit: ``Schechter``, ``DoubleSchechter`` or
    ``Hopkins2007``.
``hPowerM``
    (BHM and BHMF) Power of h that the black hole mass scales with, which depends on
    how the masses were estimated: 0 for pitch angles, Sersic indices or velocity
    dispersions, -1 for virial or dynamical masses, -2 for bulge luminosities.
    Without this column, BHMF masses are scaled with -1.
``BinWidth``
    (QLF_Xray) Width of the luminosity bins in dex.
``Bibliographic``
    ADS bibcode of the paper, e.g. ``2024ApJ...966...74W``, or
    ``2025arXiv251011373M`` while it is only on arXiv.

Each line of a data file holds x, phi and then the errors, as given by ``Errors``:

============  ====================================================
``Errors``    Columns after x and phi
------------  ----------------------------------------------------
``Delta``     sigma (symmetric)
``ULDeltas``  sigma_upper, sigma_lower (both positive)
``LUDeltas``  sigma_lower, sigma_upper (both positive)
``ULLimits``  phi_upper, phi_lower (values, not offsets)
``LULimits``  phi_lower, phi_upper (values, not offsets)
``None``      none; used for upper limits, whose file holds x and the limit
============  ====================================================

x is a log10 mass, luminosity or star formation rate, or an AB magnitude (see the
feature table in the main README). For BHM, x is the redshift and the second column
is the black hole mass. phi is per dex, or per magnitude for magnitudes.

A parametric fit (``Type`` other than ``data``) has a single line with its
parameters, in the order of the matching method in ``number_density.py``. For example
a ``Schechter`` fit in magnitudes is ``xmin xmax M_star phi_star alpha``, and one in
stellar mass is ``xmin xmax log10M_star alpha phi_star``. The function is evaluated
between xmin and xmax.

Units and h
~~~~~~~~~~~

The loader converts every dataset to linear phi in Mpc^-3 dex^-1 (Mpc^-3 mag^-1 for
magnitudes) and to the ``h`` requested by the user, using the ``h`` of the dataset:

+---------------------------------------------+------------------------------------------------+
| Quantity                                    | Conversion from h_obs to h                     |
+=============================================+================================================+
| phi                                         | ``* (h/h_obs)^3``, or ``* h^3`` if ``HinPhi``  |
|                                             | is 1                                           |
+---------------------------------------------+------------------------------------------------+
| stellar mass, SFR or luminosity (``GSMF*``, | ``log10 x - 2 log10(h/h_obs)``, as L ~ h^-2    |
| SFRF, GLF_OIII, QLF_bolometric, QLF_Xray)   |                                                |
+---------------------------------------------+------------------------------------------------+
| magnitude (GLF_UV, QLF_UV, QLF_optical)     | ``M + 5 log10(h/h_obs)``                       |
+---------------------------------------------+------------------------------------------------+
| black hole mass (BHMF, BHM)                 | ``* (h/h_obs)^hPowerM``                        |
+---------------------------------------------+------------------------------------------------+

Only h is rescaled; differences in the other cosmological parameters are ignored.
With the default ``h=1`` the results are in h = 1 units, e.g. phi in h^3 Mpc^-3.

Correlations
------------

The features read by ``astrodatapy.correlation`` (Magorrian, Tully_Fisher, ...) have
these ``info.txt`` columns, separated by two tabs:

``DataFormat``
    The x and y columns and their errors, e.g. ``xy`` (no errors), ``xyULDeltas``
    (errors on y only) or ``xLUDeltasyLUDeltas``. The error formats are those of the
    number densities above.
``Log``
    Which of the two is in log10: ``x``, ``y``, ``xy`` or ``None``.
``NormalizationOfx``, ``NormalizationOfy``
    Factors x and y are multiplied by, e.g. ``1e9``.
``Bibliographic``
    ADS bibcode, as above.

These data are returned in the units of the papers; ``h`` is not applied to them yet
(issue #7).

Clustering
----------

``QC_2PTCF`` holds power-law fits of the quasar correlation function. Each file has
one line with the parameters of the ``Type`` (``PowerLaw`` or
``PowerLaw_2COMPONENTS``), in the order of the matching method in ``clustering.py``.
The correlation length is converted to the requested ``h``.

Adding a dataset
----------------

#. Put the files in ``<feature>/<dataset>/``, with the comment header described
   above, and any upper limits in ``<feature>/<dataset>_ulimit/``.
#. Add a line for each to ``<feature>/info.txt``.
#. Add the paper's BibTeX entry to ``bibliography.bib``.
#. Run the tests, which load every dataset of every feature::

       pytest --pyargs astrodatapy

#. Plot it next to the other data, e.g. with ``astrodatapy/utils/plots.ipynb``.

Notes on individual features
----------------------------

QLF_Xray
~~~~~~~~

The X-ray luminosity function is tabulated against ``log10(LX/erg s^-1)`` in
the rest-frame 2-10 keV band, with ``Phi = dPhi/dlogLX`` in
``Mpc^-3 dex^-1``. Because the abscissa is a luminosity rather than a
magnitude, it is rescaled as ``h^-2``, as for QLF_bolometric.

=====================  ============  =======================================================
**Dataset**            **Redshift**  **Reference**
---------------------  ------------  -------------------------------------------------------
BarlowHall2022         6.05          Barlow-Hall et al. 2023, MNRAS 519, 6055 (Swift/ExSeSS)
BarlowHall2022_ulimit  6.05          as above, 1 sigma upper limit
Meraxes_obscured       5.007, 6.088  Meraxes (DRAGONS) model line, after obscuration
Pouliasis2024          5.5           Pouliasis et al. 2024, A&A 685, A97 (Chandra + XMM)
Pouliasis2024_ulimit   5.5           as above, 1 sigma upper limits
Roster2026_DR2         0.1 - 4.9     Roster et al. 2026, arXiv:2607.27887 (eROSITA DR2)
Roster2026_eFEDS       0.1 - 4.9     Roster et al. 2026, arXiv:2607.27887 (eROSITA eFEDS)
Wolf2021               6.05          Wolf et al. 2021, A&A 647, A5 (eROSITA/eFEDS)
=====================  ============  =======================================================

The redshift column gives the value to pass as ``z_target``. Each observational
file is tabulated at the centre of its redshift bin: 5.0-6.0 for Pouliasis2024 and
5.7-6.4 for Wolf2021 and BarlowHall2022. Roster2026 provides nine redshifts,
0.1, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0 and 4.9, one file each, so use a small
``z_tol`` to select a single one.

``info.txt`` also gives the width of each dataset's luminosity bins in dex
(column ``BinWidth``: 0.5 for Pouliasis2024 and Wolf2021, 1.0 for Barlow-Hall and Roster, 0.2 for
the model line), so that a point can be drawn with an x error bar of half that width.

``Meraxes_obscured`` is a model line, not data (type ``data_sim``, drawn as a line with a band, like
``Qin2017_Tiamat`` in the other quasar luminosity functions). It is the X-ray luminosity function of the
Meraxes (DRAGONS) semi-analytic model after obscuration, i.e. the observed 2-10 keV luminosity, not absorption
corrected, in 0.2 dex bins at z = 5.007 and 6.088 (the Genesis trees end at z = 5). Its band is the
Bernoulli duty-cycle uncertainty that Meraxes saves, which does not include Poisson noise or cosmic
variance, and only bins expected to hold at least 10 AGN are listed. Its ``h`` is the model's, 0.6751.

The ``h`` of each source paper is used when rescaling: 0.68 for Wolf2021, 0.6751 for Meraxes_obscured and
0.70 for the others. Every data file carries a header giving the figure or
table its values came from (or, for the model line, how it was made), the energy band,
the column-density range and the cosmology.
