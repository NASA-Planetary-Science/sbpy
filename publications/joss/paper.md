---
title: 'sbpy: A Python module for small-body planetary astronomy'
tags:
  - python
  - astronomy
  - solar system
  - planetary science
  - small bodies
  - asteroids
  - comets
  - meteoroids
  - trojans
  - centaurs
  - kuiper belt objects
  - trans-neptunian objects
authors:
  - name: Michael Mommert
    orcid: 0000-0002-8132-778X
    affiliation: 1
  - name: Michael S. P. Kelley
    orcid: 0000-0002-6702-7676
    affiliation: 2
  - name: Miguel de Val-Borro
    orcid: 0000-0002-0455-9384	
    affiliation: 3
  - name: Jian-Yang Li
    orcid: 0000-0003-3841-9977
    affiliation: 3
  - name: Giannina Guzman
    orcid: 0000-0001-6340-8220
    affiliation: 4
  - name: Brigitta Sipőcz
    orcid:  0000-0002-3713-6337
    affiliation: 5
  - name: Josef Ďurech
    affiliation: 6
  - name: Mikael Granvik
    orcid: 0000-0002-5624-1888
    affiliation: 7
  - name: Will Grundy
    orcid: 0000-0002-8296-6540
    affiliation: 1
  - name: Nick Moskovitz
    orcid: 0000-0001-6765-6336
    affiliation: 1
  - name: Antti Penttilä
    orcid: 0000-0001-7403-1721
    affiliation: 7
  - name: Nalin Samarasinha
    orcid: 0000-0001-8925-7010
    affiliation: 3
affiliations:
 - name: Lowell Observatory, US
   index: 1
 - name: University of Maryland, US
   index: 2
 - name: Planetary Science Institute, US
   index: 3
 - name: Villanova University, US
   index: 4
 - name: DIRAC Institute, Department of Astronomy, University of Washington
   index: 5
 - name: Charles University, Prague, Czech Republic
   index: 6
 - name: University of Helsinki, Finland
   index: 7

date: 31 March 2019
bibliography: paper.bib
---

# Summary

Planetary astronomy - the study of Solar System objects with
telescopic observations from the ground or from space - utilizes a
wide range of methods that are common in observational
astronomy. However, some aspects, including the planning of
observations, as well as the analysis and interpretation of the
results, require tailored techniques and models that are unique and
disparate from those used in most other fields of
astronomy. Currently, there is no single open source software package
available to support small-body planetary astronomers in their study
of asteroids and comets in the same way in which ``Astropy``
[@astropy; @astropy2] supports the general astronomy community.

``sbpy`` is a community effort to build a Python package for
small-body planetary astronomy in the form of an ``Astropy`` affiliated
package. The goal is to collect and implement well-tested and
well-documented code for the scientific study of asteroids and comets,
including (but not limited to):

* observation planning tools tailored to moving objects;
* photometric models for resolved and unresolved observations;
* wrappers and tools for astrometry and orbit fitting;
* spectroscopy analysis tools and models for reflected solar light and
  emission from gas;
* cometary gas and dust coma simulation and analysis tools;
* asteroid thermal models for flux estimation and size/albedo estimation;
* image enhancement tools for comet comae and PSF subtraction tools;
* lightcurve and shape analysis tools;
* access tools for various databases containing orbital and physical data,
  as well as ephemerides services.

``sbpy`` is available and being maintained as a github repository at
[github.com/NASA-Planetary-Science/sbpy](https://github.com/NASA-Planetary-Science/sbpy);
documentation is available on
[sbpy.readthedocs.io](https://sbpy.readthedocs.io/en/latest/)

All functionality provided as part of ``sbpy`` has been tested against
published results in order to ensure its correctness. An internal
reference tracking system enables users to query a list of appropriate
references depending on the functionality that has been utilized. In
order to improve the performance of computationally intensive
functions, like thermal modeling of atmosphereless bodies,
C-extensions are utilized.

Designed as an ``Astropy`` affiliated package, ``sbpy`` utilizes a
wide range of functionality of ``Astropy``, including its table,
constants, units, and modeling submodules, and other affiliated
packages like ``astroquery`` [@astroquery]. ``sbpy``'s API is deliberately emulating
that of ``Astropy`` wherever possible in order to enable user-friendly
and consistent access to its functionality.

``sbpy`` is set up in a highly modular class-based fashion, roughly
separated by the individual tasks listed above. Figure 1 provides
an overview of these modules and their
interrelations as indicated through colored symbols. The main modules,
shown as black boxes, interact through data containers for
ephemerides and observations (red box symbol), orbits (blue triangle),
physical properties (green diamond), and target names (orange
cross). All modules make use of the reference tracking system (purple
circle).

![``sbpy`` module structure.](structure.png)

``sbpy`` is designed to support both professional researchers and
interested citizen scientists in their activities.

# Acknowledgments

The development of ``sbpy`` is supported by NASA PDART Grant
No. 80NSSC18K0987.

# References
