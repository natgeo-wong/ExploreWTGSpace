# [ExploreWTGSpace](https://github.com/natgeo-wong/ExploreWTGSpace)
*Exploring the different implementations and parameters of the WTG Approximation*

The Weak Temperature Gradient (WTG) approximation is
a simplified framework for large-scale atmospheric dynamics in the deep tropics
where the Coriolis force is weak.

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

!!! note
    The Temperature Gradient Relaxation framework is also known as the Weak
    Temperature Gradient framework and has been called as such in several different
    studies ([Daleu2015](@citet), [Romps2012](@citet)). I use Temperature Gradient Relaxation in order to distinguish it from
    the overall Weak Temperature Gradient approximation.