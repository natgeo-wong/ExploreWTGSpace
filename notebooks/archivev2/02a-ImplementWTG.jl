### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 681658b0-5914-11eb-0d65-bbace277d145
begin
	using Pkg; Pkg.activate()
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 6dce35fc-5914-11eb-0ce2-0d4e164e1898
begin
	@quickactivate "ExploreWTGSpace"
	using PlutoUI
	using Printf
	using SpecialFunctions
	using StatsBase
	
	using PyCall, LaTeXStrings
	using PNGFiles, ImageShow
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the ExploreWTGSpace project..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 02a. Implementing the WTG Approximation in SAM

In this notebook, we discuss the different implementations of the WTG approximation, namely the:
1. DGW Implementation of Blossey et al. [2009]
2. WTG Implementation of Raymond & Zeng [2005] and their derivatives

Furthermore, we also develop a way to implement the WTG approximation gradually in the System of Atmospheric Modelling, that gradually increases the strength of the WTG momentum-damping parameter from a near/pseudo-RCE state.
"

# ╔═╡ 63f1f675-61db-4996-ac2a-157701c3da8b
md"
### A. Implementing the WTG Approximation in Models

#### i. The Damped Gravity Wave (DGW) Implementation

Blossey et al. [2009] were the first to offically implement the WTG approximation in SAM, using the DGW method. The strength of the DGW implementation is affected by the momentum-damping coefficient $a_m$ in the below formula:

$$\frac{\partial}{\partial p}
\left( \frac{f^2+a_m^2}{a_m}\frac{\partial\omega'}{\partial p} \right)
\approx \frac{k^2R_d}{p} T_v'$$

In our CRM simulations, we assume that $f$, the coriolis parameter, is 0, within the tropics (and therefore negligible), and that  $a_m$ is constant in height (which was noted in Blossey et al. [2009], though they favoured that $a_m \propto p/p_\text{ref}$).  Thus, the formula simplifies itself to

$$\frac{\partial^2\omega'}{\partial p^2} \approx \frac{k^2}{a_m} \frac{R_d}{p}T_v'$$

Therefore, we see that for a bigger $a_m$, the induced $\omega'$ (or the WTG-induced vertical adjustment) is smaller.  Thus, bigger $a_m$ would result in less deviation from an RCE solution.  In other words, the magnitude of the WTG forcing is proportional to $a_m^{-1}$.

#### ii. The Weak Temperature Gradient (WTG) Implementation

The  assumes that the vertical advection of virtual potential temperature $w\partial_z\theta_v$ removes the differences in the buoyancy between the CRM and large-scale environment over a time-scale $\tau$, such that at a height in the free troposphere $z_i$:

$$\left. w_\text{wtg}(z_i) 
\frac{\partial{\overline{\theta}}}{\partial{z}} \right|_{z=z_i}
= \frac{\overline{\theta}(z_i)-\theta_0(z_i)}{\tau} \cdot f(z_i)$$

Where $\theta$ is the potential temperature, and $\overline{\theta}$ is its horizontal average.  This method of implementation in a CRM was first done by Raymond and Zeng [2005], and has since become the most popular method of implementing the WTG approximation due to the straightforward conceptual picture it provides, though there are other variations.  We would expect that as $\tau$ decreases, the degree to which the WTG approximation must be obeyed in the troposphere must increase,

There are different forms of $f(z_i)$ that have been used:
* From Raymond and Zeng [2005]: $f(z_i) = \sin\frac{\pi z_i}{h}$
* From Daleu et al. [2015]: $f(z_i) = 1$

We also implement an alternative: $f(z_i) = \tau_a * a * \sin\frac{\pi z_i}{h} + \tau_b * b * \sin\frac{2\pi z_i}{h}$, where $a$ and $b$ are terms that vary with every timestep depending on the magnitudes of the sine and half-sine curve decompositions of the $\overline{\theta}(z_i)-\theta_0(z_i)$ curve.
"


# ╔═╡ de98b3f4-d550-4f00-966d-9e55a034ec51
md"
### B. Implementing a smooth transition

We implement a transitionary period of the first 100 days (out of 500 model run days) in our WTG runs, such that the model scales up from a pseudo-RCE state to the user-defined strength of the WTG approximation.  Let us take $a_m$ (or the DGW implementation) as an example:

* The final momentum-damping strength is $a_{m.0}$
* We take $t_{wtg}$ to be the time $t$ where $a_m = a_{m.0}$
* At time $t=0$, we have that $a_m \rightarrow \infty$

All in all, if we want $w_{wtg}$ to increase linearly over time, we must have that

$$w_{wtg}(t) \propto \frac{1}{a_m(t)} = \frac{1}{a_{m.0}} * \frac{t}{t_{wtg}}$$
$$\therefore a_m(t) = a_{m.0} * \frac{t_{wtg}}{t}, t < t_{wtg};\quad a_m(t) = a_{m.0}, t \geq t_{wtg}$$
"

# ╔═╡ f66daa4b-fe38-47a8-82e5-adfdc91a6a98
function t2am(t::Real, twtg::Real, am0::Real)

	if t < twtg
		  return am0 * twtg / t
	else; return am0
	end
	
end

# ╔═╡ c5e8f839-6ace-4a86-bbb0-e0c725f25c70
function am2wtg(am)
	return 1 / am
end

# ╔═╡ d1eda967-dea3-4b84-9a46-98ad6cdde41e
time = 0 : 0.001 : 1; time = time[2:end]

# ╔═╡ 215c1a45-22f3-473f-93ac-ae98c5aa2f5d
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=2,aspect=1.5,sharey=0)
	clr = pplt.Colors("Blues",9)
	amvec = [0.5,0.5*sqrt(2),1,sqrt(2),2]
	for ii in eachindex(amvec)
		amt = t2am.(time,0.2,amvec[ii])
		axs[1].plot(time,amt,c=clr[ii+2])
		axs[2].plot(
			time,am2wtg.(amt),c=clr[ii+2],
			label=L"$a_{m.0}$ = " * "$(@sprintf("%.2f",amvec[ii]))",
			legend="r",legend_kw=Dict("ncol"=>1,"frame"=>false)
		)
	end

	axs[1].format(yscale="log",ylabel=L"$a_m$",ylim=(0.2,50))
	axs[2].format(ylabel=L"$w_{wtg}\propto {a_m}^{-1}$",ylim=(0,2.5),ytickloc="right")

	for ax in axs
		ax.plot([0.2,0.2],[0,50],c="k",linestyle="--")
		ax.format(
			xlim=(0,1),xlabel="Nondimensionalized Time",
			suptitle="Implementing the WTG Approximation"
		)
	end
	
	fig.savefig(plotsdir("02a-ImplementWTG.png"),dpi=400,transparent=false)
	load(plotsdir("02a-ImplementWTG.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─63f1f675-61db-4996-ac2a-157701c3da8b
# ╟─de98b3f4-d550-4f00-966d-9e55a034ec51
# ╠═f66daa4b-fe38-47a8-82e5-adfdc91a6a98
# ╠═c5e8f839-6ace-4a86-bbb0-e0c725f25c70
# ╠═d1eda967-dea3-4b84-9a46-98ad6cdde41e
# ╟─215c1a45-22f3-473f-93ac-ae98c5aa2f5d
