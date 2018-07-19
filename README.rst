AstroDataPy
===========

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge

Python tools to collect astronomical data.

This is a python package dedicated to collect up-to-date astronomical 
data from both observational and modelling campaigns, with a focus on 
galaxy properties. Current package includes 1) statistical properties 
such as the number densities of galaxies (as functions of stellar mass, 
UV magnitude, star formation rate) and AGN (as functions of black hole 
mass, quasar UV/optical/bolometric luminosities); 2) correlations between 
galaxies properties such as the Magorrian relation, Tully Fisher relation, 
Disk size - stellar mass relation, and halo - stellar mass relation, etc; 
and 3) clustering of quasars such as the two point correlation function.

Current modelling results include DRAGONS (Meraxes).

Installation
============

pip install git+https://github.com/qyx268/astrodatapy

Usage
=====

Example 1
---------
read quasar two point correlation function at z = 4 with a redshift range of [3.5, 4.5]

.. code:: python

    from astrodatapy.clustering import clustering
    obs = clustering(feature = 'QC_2PTCF', z_target = 4.0, z_tol = 0.5)
    # show all available data of QC_2PTCF
    print(obs.available_observation)
    # show redshifts of all available data of QC_2PTCF
    print(obs.z_available_observation)
    # show the target data of QC_2PTCF at z = 4
    print(obs.target_observation)
    

Example 2
---------
read Magorrian Relation at redshift 0, output with h=0.678, and do not show information

.. code:: python

  from astrodatapy.correlation import correlation
  obs = correlation(feature = 'Magorrian', z_target = 0, quiet = 1)


Example 3
---------
plot galaxy stellar mass function at redshift 5 and show labels

.. code:: python

  import matplotlib.pyplot as plt
  from astrodatapy.number_density import number_density
  fig, ax = plt.subplots(1, 1)
  obs = number_density(feature = 'GSMF', z_target = 5.0)
  for ii in range(obs.n_target_observation):
    data       = obs.target_observation['Data'][ii]
    label      = obs.target_observation.index[ii]
    # error bar plot
    ax.errorbar(data[:,0],  data[:,1], yerr = [data[:,1]-data[:,3], data[:,2]- data[:,1]], label = label)

    # uncertainty range plot
    #ax.plot(data[:,0], data[:,1], label = label)
    #ax.fill_between(data[:,0], data[:,2], data[:,3], alpha=0.5)

More examples can be found in astrodatapy/utils/plots.ipynb and Documentation.

Documentation
=============

http://astrodatapy.readthedocs.io

Features
============

Number density
--------------

==============             ==========================================
**Features**               **Descriptions**
--------------             ------------------------------------------
BHM                        Black Hole Mass
BHMF                       Black Hole Mass Function
GLF_UV                     Galaxy Luminosity Function -- UV
GSMF                       Galaxy Stellar Mass Function -- all
GSMF_Blue                  Galaxy Stellar Mass Function -- blue
GSMF_Bulge                 Galaxy Stellar Mass Function -- bulge
GSMF_Disk                  Galaxy Stellar Mass Function -- disk
GSMF_Quiescent             Galaxy Stellar Mass Function -- quiescent
GSMF_Red                   Galaxy Stellar Mass Function -- red
QLF_bolometric             Quasar Luminosity Function -- bolometric
QLF_optical                Quasar Luminosity Function -- optical
QLF_UV                     Quasar Luminosity Function -- UV
SFRF                       Star Formation Rate Function
==============             ==========================================

Correlation
-----------


=========================  ================================================
**Features**               **Descriptions**
-------------------------  ------------------------------------------------
BHM                        Black Hole Mass
Magorrian                  Black Hole - Galaxy Bulge Mass Scaling Relation
Tully_Fisher               Mass - Velocity of Spiral Galaxies
DiskSize_StellarMass       DiskSize - StellarMass
GasFraction_StellarMass    GasFraction - StellarMass
sSFR_StellarMass_Blue      sSFR - StellarMass -- blue
HaloMass_StellarMass       HaloMass - StellarMass
HaloMass_StellarMass_Blue  HaloMass - StellarMass -- blue
HaloMass_StellarMass_Red   HaloMass - StellarMass -- red
=========================  ================================================

Clustering
----------

==============             =================================================
**Features**               **Descriptions**
--------------             -------------------------------------------------
QC_2PTCF                   Quasar Clustering -- 2 point correlation function
==============             =================================================

License
=======

* Free software: BSD license

* This project is Copyright (c) Yuxiang Qin and licensed under the terms of the BSD 3-Clause license. See the licenses folder for more information.

Contributors
============

* Yuxiang Qin (The University of Melbourne)
