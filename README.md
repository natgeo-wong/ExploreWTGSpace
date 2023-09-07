# **<div align="center">ExploreWTGSpace</div>**

<p align="center">
  <a href="https://www.repostatus.org/#active">
    <img alt="Repo Status" src="https://www.repostatus.org/badges/latest/active.svg?style=flat-square" />
  </a>
  <a href="https://mit-license.org">
    <img alt="MIT License" src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square">
  </a>
  <a href="https://natgeo-wong.github.io/ExploreWTGSpace/dev/">
    <img alt="Latest Documentation" src="https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square">
  </a>
</p>

**Authored By:** 
* Nathanael Wong (nathanaelwong@fas.harvard.edu)
* Professor Kuang Zhiming

The Weak Temperature Gradient (WTG) approximation is
a simplified framework for atmospheric dynamics in the deep tropics where the
Coriolis force is weak.

This project aims to explore the different implementations of the WTG Approximation
in small-domain cloud-resolving models and how these schemes interact with the vast
parameter space with variables such as:
* Horizontal domain size and resolution
* Sea-surface temperature
* Interactive and non-interactive radiation, presence of a diurnal cycle

Over the past few decades, there are two popular frameworks that have emerged that
implement the Weak Temperature Gradient approximation:
1. Temperature Gradient Relaxation (TGR)
2. Damped Gravity Wave (DGW)

This repository contains files that generate the prm files required to run our
experiments in the System of Atmospheric Modelling v6.11.8 as modified by the Kuang Lab
at Harvard University, along with notebooks that will aid in the analysis of our model
results.

The resulting sounding (snd) files from RCE spinups found in `exp/snd` can also be used
in future experiments.

## Related Repositories:

* [2023GL104350](https://github.com/natgeo-wong/2023GL104350)

  This project is a subset of ExploreWTGSpace, containing the relevant experiments,
  notebooks and scripts required to generate the data for the paper submission to GRL
  with the ID `2023GL104350`.

* [SelfAggregation](https://github.com/natgeo-wong/SelfAggregation)

  A complement to the ExploreWTGSpace project.  Investigates the self-aggregation of
  convection in the System of Atmospheric Modelling for similar horizontal resolution
  allowing for direct comparison to WTG results.

## Installation

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project.  To (locally) reproduce this project, do the
following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> ] activate .
    Activating environment at `path/to/this/project`

   (ExploreWTGSpace) pkg> instantiate
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.

## Other Acknowledgements
> Project Repository Template generated using [DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl) created by George Datseris.