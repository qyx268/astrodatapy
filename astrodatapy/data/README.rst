Data format
===========

This directory holds the data that astrodatapy loads. Each feature (``GLF_UV``,
``GSMF``, ``Magorrian``, ...) has its own directory::

    <feature>/info.txt                  one line per dataset
    <feature>/<dataset>/z5pt0.dat       one file per redshift

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
    10^-5 Mpc^-3 mag^-1 (for BHM, the factor that gives the black hole masses in
    Msun, e.g. ``1e6`` for a table in 10^6 Msun). It must be 1 when ``Log`` is 1.
``Errors``
    Layout of the error columns, see below.
``Log``
    1 if phi and its errors are log10 values (errors in dex), 0 if they are linear.
``HinPhi``
    1 if phi is in h^3 Mpc^-3, 0 if it is in Mpc^-3.
``Type``
    ``data``, ``data_sim`` (a model, drawn as a line with a band), or a parametric
    fit: ``Schechter``, ``DoubleSchechter`` or ``Hopkins2007``. (Older versions kept
    upper limits as separate datasets of Type ``dataULimit``; these are still read,
    with every row an upper limit.)
``hPowerM``
    (BHM and BHMF) Power of h that the black hole mass scales with, which depends on
    how the masses were estimated: 0 for pitch angles, Sersic indices or velocity
    dispersions, -1 for virial or dynamical masses, -2 for bulge luminosities.
    Without this column, BHMF masses are scaled with -1.
``BinWidth``
    (Optional) Width of the bins of x, in the units of x as returned (dex or mag): a
    number if every row has the same width, ``column`` if the data file has an extra last
    column with each row's width (for bins of different widths), or ``None`` if unknown.
    The loader returns the widths in the ``BinWidth`` column of ``target_observation``
    (``nan`` when unknown) without converting them: the ``MAG`` band conversions and the h
    conversion only shift x, which leaves widths unchanged, except ``LBsol2M1450``, which
    turns luminosities into magnitudes, so give its widths in magnitudes.
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
``None``      none (parametric fits)
============  ====================================================

x is a log10 mass, luminosity or star formation rate, or an AB magnitude (see the
feature table in the main README). For BHM, x is the redshift and the second column
is the black hole mass, which ``NormalizationOfPhi`` converts to Msun. phi is per
dex, or per magnitude for magnitudes.

Upper limits go in the same file as the other points of the paper at that redshift,
as rows whose error columns are ``nan``: ``x  limit  nan`` for ``Delta``, ``x  limit
nan  nan`` for the other formats. The loader marks them in the ``ULimit`` column of
``target_observation`` and sets their upper and lower bounds to the limit. A redshift
with only upper limits gets a file of such rows.

Models (``data_sim``) are tabulated like ``data``.

Bin widths are optional (``BinWidth`` above). Please give them for new datasets when the
paper does. Most existing datasets don't record their bin widths yet (QLF_Xray and
Baldry+2012 in GSMF do),
so their ``BinWidth`` is ``nan``; widths can be added to them later, checked against the
papers (a ``BinWidth`` column in ``info.txt`` for constant widths, or a last column of
widths in the data files), without changing anything else.
A parametric fit (``Type`` ``Schechter``, ``DoubleSchechter`` or ``Hopkins2007``) has
instead a single line with its parameters, in the order of the matching method in
``number_density.py``. For example
a ``Schechter`` fit in magnitudes is ``xmin xmax M_star phi_star alpha``, and one in
stellar mass is ``xmin xmax log10M_star alpha phi_star``. The function is evaluated
between xmin and xmax.

Units and h
~~~~~~~~~~~

The loader converts phi (or, for BHM, the black hole mass) to linear values in
Mpc^-3 dex^-1 (Mpc^-3 mag^-1 for magnitudes), and then converts every dataset from its
own ``h``, h_obs, to the ``h`` requested by the user. x stays as it is stored, a log10
value or a magnitude, and is shifted:

+-----------------------------------------------+----------------------------------------------+
| Column                                        | Conversion from h_obs to h                   |
+===============================================+==============================================+
| phi (linear)                                  | ``phi * (h/h_obs)^3``, or ``phi * h^3`` if   |
|                                               | ``HinPhi`` is 1                              |
+-----------------------------------------------+----------------------------------------------+
| x = log10 of a stellar mass, SFR or           | ``x - 2 log10(h/h_obs)``, as L ~ h^-2        |
| luminosity (``GSMF*``, SFRF, GLF_OIII,        |                                              |
| QLF_bolometric, QLF_Xray)                     |                                              |
+-----------------------------------------------+----------------------------------------------+
| x = magnitude (GLF_UV, QLF_UV, QLF_optical)   | ``x + 5 log10(h/h_obs)``                     |
+-----------------------------------------------+----------------------------------------------+
| x = log10 M_BH (BHMF)                         | ``x + hPowerM log10(h/h_obs)``               |
+-----------------------------------------------+----------------------------------------------+
| black hole mass, linear (BHM, second to       | ``M_BH * (h/h_obs)^hPowerM``                 |
| fourth columns; x is the redshift)            |                                              |
+-----------------------------------------------+----------------------------------------------+

Only h is rescaled; differences in the other cosmological parameters are ignored.
With the default ``h=1`` the results are in h = 1 units, e.g. phi in h^3 Mpc^-3.

Correlations
------------

The features read by ``astrodatapy.correlation`` (Magorrian, Tully_Fisher, ...) have
these ``info.txt`` columns, separated by whitespace:

``DataFormat``
    The x and y columns and their errors, e.g. ``xy`` (no errors), ``xyULDeltas``
    (errors on y only) or ``xLUDeltasyLUDeltas``. The error formats are those of the
    number densities above.
``Log``
    Which of the two is in log10: ``x``, ``y``, ``xy`` or ``None``.
``NormalizationOfx``, ``NormalizationOfy``
    Factors x and y are multiplied by, e.g. ``1e9``.
``h``
    H0/(100 km/s/Mpc) that the tabulated values assume; 1 if they are tabulated in
    h-scaled units such as h^-1 Msun.
``hPowerx``, ``hPowery``
    Powers of h that x and y scale with: -2 for stellar, gas or bulge masses derived
    from luminosities or fluxes, -1 for sizes in kpc and for halo masses, 0 for
    velocities and for ratios such as gas fractions and sSFR; for black hole masses,
    as ``hPowerM`` above.
``IMF``
    IMF of the stellar masses (``Salpeter``, ``DietSalpeter``, ``Kroupa`` or
    ``Chabrier``), or ``None`` to leave them as they are. The stellar-mass axis
    (x, or y for Tully_Fisher) is converted to the ``IMF_out`` requested by the user,
    and the gas fraction of GasFraction_StellarMass, M_gas/M*, the other way.
``Bibliographic``
    ADS bibcode, as above.

x and y are converted to linear values (``Log``) and multiplied by their
normalisations, then by (h/h_obs)^hPowerx and (h/h_obs)^hPowery, where h_obs is the
``h`` of the dataset. A dataset without ``h`` or ``IMF`` is returned as tabulated.

Clustering
----------

``QC_2PTCF`` holds power-law fits of the quasar correlation function. Each file has
one line with the parameters of the ``Type`` (``PowerLaw`` or
``PowerLaw_2COMPONENTS``), in the order of the matching method in ``clustering.py``.
The correlation length is converted to the requested ``h``.

Redshift evolution (IGM)
------------------------

The features read by ``astrodatapy.evolution`` (``IGM_Emissivity``, ``IGM_Gamma``,
``IGM_MFP``, ``IGM_T0``) have **one file per dataset**, ``<feature>/<dataset>.dat``, whose
rows start with the redshift; the file name carries no redshift. Their ``info.txt`` columns
(whitespace separated):

``DataFormat``
    As for the correlations, with x the redshift: e.g. ``xLULimitsyLUDeltas`` (z, lower and
    upper edge of the redshift bin, value, lower and upper error), ``xyLUDeltas`` (no bin),
    ``xLULimitsyDelta`` (symmetric error), ``xyULLimits`` (upper and lower bound as values),
    ``xy`` (no errors, e.g. a simulation curve).
``Log``
    ``y`` if the values (and bounds) are log10, else ``None``.
``NormalizationOfy``
    Factor the linear values are multiplied by: e.g. ``1e51`` for an emissivity in
    1e51 s^-1 cMpc^-3, ``1e12`` for a log10 photoionization rate in s^-1 (returned in
    1e-12 s^-1), ``1e4`` for a temperature in 1e4 K.
``h``, ``hPowery``
    The value scales as h^hPowery: ``-1`` for lengths tabulated in h^-1 Mpc (with ``h = 1``)
    or h70^-1 Mpc (with ``h = 0.7``); ``0`` when the paper gives physical units for its
    cosmology (``h`` is then for information).
``Comoving``
    ``1`` for a length tabulated in comoving units, which the loader divides by 1+z (the
    mean free path is returned in proper Mpc).
``Type``
    ``data`` for an observation, ``data_sim`` for a simulation or model curve.
``Method``
    How the values were obtained, e.g. ``LyC_stack``, ``Lya_forest_tau_eff``,
    ``QSO_proximity``, ``Lya_forest_power_spectrum``, or the name of a simulation.
``Bibliographic``
    ADS bibcode.

Upper limits are rows whose errors are ``nan``, lower limits rows whose errors are
``inf``; the loader flags them (``ULimit``, ``LLimit``) and gives them the value as both
bounds. Store the values as published, with a header saying where they come from.

The bounds usually bracket the value, but not always: a posterior given as its MAP
(best-fitting) model with the 16th/84th percentiles (``Qin2025``) can have the MAP outside
them, so ``value - lower`` or ``upper - value`` can be negative (``plt.errorbar`` refuses
that; draw such a dataset as a band, e.g. ``fill_between``, with the value as a line).

Adding a dataset
----------------

#. Put the files in ``<feature>/<dataset>/``, with the comment header described
   above; upper limits are rows with ``nan`` errors in the same files. Give the bin
   widths (``BinWidth``) if the paper does.
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
BarlowHall2022         6.05          Barlow-Hall et al. 2023, MNRAS 519, 6055 (Swift/ExSeSS);
                                     with a 1 sigma upper limit
Meraxes_obscured       5.007, 6.088  Meraxes (DRAGONS) model line, after obscuration
Pouliasis2024          5.5           Pouliasis et al. 2024, A&A 685, A97 (Chandra + XMM);
                                     with 1 sigma upper limits
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
