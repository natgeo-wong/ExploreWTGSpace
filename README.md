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

In this project, we explore the instability caused by implementing the Weak Temperature Gradient approximation onto a small-domain Radiative-Convective Equilibrium simulation under a series of different parameter spaces including:
* Horizontal domain size and resolution
* Sea-surface temperature
* Vertical domain resolution
* Interactive and non-interactive radiation

## Installation

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project.  To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> ] activate .
    Activating environment at `path/to/this/project`

   (PiPWV) pkg> instantiate
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.

## Other Acknowledgements
> Project Repository Template generated using [DrWatson.jl](https://github.com/JuliaDynamics/DrWatson.jl) created by George Datseris.