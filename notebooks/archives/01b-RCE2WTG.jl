### A Pluto.jl notebook ###
# v0.16.1

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
	using Printf
	using SpecialFunctions
	using Statistics
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 6b. Transitioning from RCE to WTG

In this notebook, we investigate and develop a way to implement the WTG forcing gradually in the System of Atmospheric Modelling.  The sudden introduction of WTG large-scale forcing often causes a model to enter a \"shocked\" state that unnaturally forces the model into a different state.  Here, we develop a method that gradually increases the strength of the WTG momentum-damping parameter from a near/pseudo-RCE state.
"

# ╔═╡ 63f1f675-61db-4996-ac2a-157701c3da8b
md"
### A. Implementing the WTG Approximation in Models

Blossey et al. [2009] were the first to implement to WTG approximation in SAM.  The strength of the WTG is affected by the momentum-damping coefficient $a_m$ in the below formula:

$$\frac{\partial}{\partial p}
\left( \frac{f^2+a_m^2}{a_m}\frac{\partial\omega'}{\partial p} \right)
\approx \frac{k^2R_d}{p} T_v'$$

In our CRM simulations, we assume that $f$, the coriolis parameter, is 0, within the tropics (and therefore negligible).  Thus, the formula simplifies itself to

$$\frac{\partial}{\partial p}
\left( a_m\frac{\partial\omega'}{\partial p} \right)
\approx \frac{k^2R_d}{p} T_v'$$

If we do a further simplification by assuming that $a_m$ is constant in height (which was noted in Blossey et al. [2009], though they favoured that $a_m \propto p/p_\text{ref}$), this further simplifies the equation down to

$$\frac{\partial^2\omega'}{\partial p^2} \approx \frac{k^2}{a_m} \frac{R_d}{p}T_v'$$

Therefore, we see that for a bigger $a_m$, the induced $\omega'$ (or the WTG-induced vertical adjustment) is smaller.  Thus, bigger $a_m$ would result in less deviation from an RCE solution.  In other words, the magnitude of the WTG forcing is proportional to $a_m^{-1}$.
"


# ╔═╡ 5c4cb1a2-c566-4208-b799-7e9644f197dd
wtgstrength(am::Real) = 1/am

# ╔═╡ de98b3f4-d550-4f00-966d-9e55a034ec51
md"
### B. Implementing a smooth transition

Implementing the Weak-Temperature Gradient approximation suddenly can cause a \"shock\" to the model, which will then transfer into a different regime state.  Therefore, we implement a smooth transition from a pseudo-RCE state ($a_m\gg1$) to a WTG state ($a_m = a_{m.0}$).  This was by varying $a_m$ in the form of an error function.  This will allow the increase in the strength of the WTG approximation to taper off gently, instead of increasing more and more rapidly, as it approaches $a_{m.0}$.

The error function is given by

$$\text{erf}\> t = \frac{2}{\sqrt{t}} \int_0^t e^{-x^2} \>\text{d}{x}$$

We use the error function to vary the momentum damping parameter $a_m$, or more specifically, the parameter $x$ where $x = \log_2a_m$.

The formula for $x$ as a function of time $t$, is given by

$$x(t) = \left( \frac{10-\log_2a_{m.0}}{2} \right)
\text{erf}\> \left( \frac{t_\text{max}/2-t}{t_\text{max}/5} \right)
+ 5 + \frac{1}{2}\log_2a_{m.0}$$
$$a_m(t) = 2^{x(t)}$$

Where $t_\text{max}$ is a scaling parameter, that defines the time when $a_m \approx a_{m.0}$ as a fraction of the total time that the WTG approximation is applied.  So if the WTG approximation is applied over 100 days, then if $t_\text{max} = 0.25$, $a_m \approx a_{m.0}$ occurs at day 25.
"

# ╔═╡ 2171952b-9551-40a8-8838-7d79d8977303
twtg_scale = 0.8

# ╔═╡ f66daa4b-fe38-47a8-82e5-adfdc91a6a98
begin
	pplt.close(); f,axs = pplt.subplots(ncols=3,aspect=1,axwidth=2,sharey=0);

	for am in 2 .^(0:10)
		tmax = 1; twtg_max = twtg_scale * tmax
		t = 0:0.001:tmax;
		test = erf.((twtg_max/2 .-t)/(twtg_max/5)) * (10-log10(am)/log10(2))/2 .+ 5 .+ log(2,am)/2
		amc = 2 .^(test)
		wtg = wtgstrength.(amc)

		axs[1].plot(t,test,c="k",lw=1)
		axs[2].plot(t,amc,c="k",lw=1)
		axs[3].plot(t,wtg,c="k",lw=1)
	end

	axs[1].plot([1,1]*twtg_scale,[0,10],c="r")
	axs[2].plot([1,1]*twtg_scale,[0,1024],c="r")
	axs[3].plot([1,1]*twtg_scale,[0,1],c="r")

	axs[1].format(ylim=(0,10),ltitle=L"(a) $\log_2 (a_m)$")
	axs[2].format(ylim=(0,1024),ltitle=L"(b) $a_m$")
	axs[3].format(ylim=(0,1),ltitle=L"(c) WTG $\propto$ $a_m^{-1}$",
		xlabel="Nondimensionalized time")
	
	for ax in axs
		ax.format(xlim=(0,1))
	end

	f.savefig(plotsdir("wtgcoeff.png"),transparent=false,dpi=200)
	load(plotsdir("wtgcoeff.png"))
end


# ╔═╡ 08173a11-bb0d-4566-9602-017815306de4
md"
### C. Implementation in SAM

We encoded this transition in WTG strength via an error function into the `forcing.f90` file of SAM.  Alongside this, we made the following adjustments to the SAM model:
* Ensemble perturbations were added to the initialized **thermodynamic** profile rather than the **vertical** profile, which allows consistency in the vertical pressure grid constructed between ensemble models.
* When both vertical and pressure coordinates are given, the temperature constructed from the potential temperature in the sounding file is based off **both given vertical and pressure coordinates**, rather than the reconstructed pressure profile, which differs from the original.

The number of ensemble members ran depends on the number of members needed for both the dry and wet states to appear.  This can range from 5, to 15.  Regardless, we adjust our analysis code to allow for varying number of ensemble members for each final $a_m$.
"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
expname = "P064km301d7"

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	configvec = [
		# "damping001",
		"damping002","damping004","damping006","damping008",
		"damping010","damping011","damping012","damping016","damping020",
		"damping024","damping023","damping028","damping032","damping040",
		"damping045","damping048","damping056","damping064","damping080",
		"damping090","damping096","damping112","damping128","damping160",
		"damping192","damping224","damping256","damping512",
	]
	ncon = length(configvec)
	blues = pplt.Colors("Blues",(ncon+2))
	reds  = pplt.Colors("Reds",(ncon+2))
	teals = pplt.Colors("Greens",(ncon+2))
	lgd = Dict("frame"=>false,"ncols"=>1)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 223b4286-8811-11eb-0e67-4da65e1999a5
function daymean(data)
	
	return dropdims(mean(reshape(data,1,:),dims=1),dims=1)
	
end

# ╔═╡ 55230f4a-7661-11eb-1c37-8b022b95e08e
begin
	pplt.close()
	fts,ats = pplt.subplots(nrows=3,aspect=3,axwidth=4,hspace=0.2,sharey=0)
	
	for ic in 1 : ncon
		config = configvec[ic]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		imem = 0
		
		while imem < 100; imem += 1
			fnc = outstatname(expname,configvec[ic],false,true,imem)
			if isfile(fnc)
				_,_,t = retrievedims(fnc); t = t .- floor(t[1])
				td = daymean(t)
				pr = retrievevar("PREC",fnc)
				pa = retrievevar("AREAPREC",fnc)
				pa = daymean(pr ./ pa)
				pr = daymean(pr)
				sw = daymean(retrievevar("SWNS",fnc))
				lw = daymean(retrievevar("LWNS",fnc))
				sh = daymean(retrievevar("SHF",fnc))
				lh = daymean(retrievevar("LHF",fnc))
				pw = daymean(retrievevar("PW",fnc))
				seb = sw .- lw .- sh .- lh
				ats[1].plot(td,pr,lw=1,color=blues[ic+1])
				ats[3].plot(td,seb,lw=1,color=reds[ic+1])
				if imem == 1
					constr = @sprintf("%d",config)
					ats[2].plot(
						td,pw,lw=1,color=teals[ic+1],
						label=(L"$a_m =$" * " $(constr)"),
						legend="r",legend_kw=lgd
					)
				else
					ats[2].plot(td,pw,lw=1,color=teals[ic+1])
				end
			end
		end
		
	end
	
	ats[1].format(
		ylabel=L"Rainfall / mm day$^{-1}$",yscale="log",
		ylocator=10. .^(-3:3),
		# yscale_kw=Dict("linthresh"=>0.1),
		ylim=(0.005,200),
		suptitle=expname,ultitle="(a)"
	)
	
	ats[2].format(
		ylim=(0,100),ylabel="PW / mm",
		suptitle=expname,ultitle="(b)"
	)
	
	ats[3].format(
		xlim=(00,500),xlabel="Time / Days",
		ylim=(-250,250),ylabel=L"SEB / W m$^{-2}$",ylocator=(-3:3)*100,
		suptitle=expname,ultitle="(c)"
	)
	
	fts.savefig(plotsdir(
		"rce2wtg-$(expname).png"),
		transparent=false,dpi=150
	)
	load(plotsdir("rce2wtg-$(expname).png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─63f1f675-61db-4996-ac2a-157701c3da8b
# ╠═5c4cb1a2-c566-4208-b799-7e9644f197dd
# ╟─de98b3f4-d550-4f00-966d-9e55a034ec51
# ╠═2171952b-9551-40a8-8838-7d79d8977303
# ╠═f66daa4b-fe38-47a8-82e5-adfdc91a6a98
# ╟─08173a11-bb0d-4566-9602-017815306de4
# ╠═d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─223b4286-8811-11eb-0e67-4da65e1999a5
# ╟─55230f4a-7661-11eb-1c37-8b022b95e08e
