### A Pluto.jl notebook ###
# v0.15.1

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
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 1c. Momentum Damping Strength

Previous studies (e.g. Emanuel et al. [2014]) have shown that imposing the WTG approximation onto a model that has reached RCE causes it to transition into one of two regimes: a wet regime and a dry regime.  These regimes are analogues to the wet and dry regimes found in a large-area domain, where self-aggregation of convection naturally occurs in RCE.

In this notebook, we explore the characteristics of these wet and dry regimes under different $a_m$ strength.
"

# ╔═╡ b6892634-9199-11eb-38d5-8dda8da16ed7
md"
### A. The Wet and Dry Regimes in WTG Ensemble Runs

We recall the momentum damping equation (assuming $a_m$ is constant with height):

$$\frac{\partial^2\omega'}{\partial p^2} \approx \frac{k^2}{a_m} \frac{R_d}{p}T_v'$$

Recall that $k$ is the wavenumber.  Therefore, increasing $a_m$ increases the wavenumber required for the WTG response to be the same.  Or in other words:

$$k' = \frac{k}{\sqrt{a_m}}$$

The pseudo-wavenumber $k'$ is smaller than the actual wavenumber $k$, which implies that the wavelength is much, much larger.  This is an analogue to the domain being farther away from the baseline RCE domain, which is taken to be a large-scale domain average.  So as $a_m$ increases, we should see the dry and wet states converge back into the initial RCE state.
"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
expname = "P064km305d0"

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	configvec = [
		"damping001","damping002","damping004","damping006","damping008",
		"damping010","damping011","damping012","damping016","damping020",
		"damping024","damping023","damping028","damping032","damping040",
		"damping045","damping048","damping056","damping064","damping080",
		"damping090","damping096","damping112","damping128","damping160",
		"damping192","damping224","damping256","damping512",
	]
	ncon = length(configvec)
	blues = pplt.Colors("Blues",(ncon+2))
	grns  = pplt.Colors("Teal",(ncon+2))
	brwns = pplt.Colors("Brown",(ncon+2))
	lgd = Dict("frame"=>false,"ncols"=>4)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 7160ccfe-2189-4341-9e83-e700773f3d3e
ndy = 250

# ╔═╡ 55230f4a-7661-11eb-1c37-8b022b95e08e
begin
	pplt.close()
	fts,ats = pplt.subplots(ncols=5,aspect=0.5,axwidth=1.2,sharex=0)
	cr = []
	rh = zeros(64)
	pl = zeros(64)
	inc = 0
	
	for ic in 1 : ncon
		config = configvec[ic]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		imem = 0
		
		while imem < 100; imem += 1
			fnc = outstatname(expname,configvec[ic],false,true,imem)
			if isfile(fnc); global inc += 1
				_,p,_ = retrievedims(fnc); global pl += p; p = p * 100
				pw = retrievevar("PW",fnc)[(end-(ndy-1)):end]
				ta = retrievevar("TABS",fnc)[:,(end-(ndy-1)):end]
				qv = retrievevar("QV",fnc)[:,(end-(ndy-1)):end] / 1000
				ri = calcrh(qv,ta,p)
				sw = calcswp(ri,qv,p)
				global cr = vcat(cr,pw ./ sw / 10)
				global rh = cat(rh,ri,dims=2)
			end
		end
		
	end
	
	rh = rh[:,2:end]
	pl = pl / inc
	
	md"Loading data for the last $ndy days ..."
	
end

# ╔═╡ f706782e-6b0e-40fa-ba21-ff87f2594f5a
begin
	bin = 0:2:100; nbin = length(bin); binp = (bin[1:(end-1)].+bin[2:end])/2
	rmean = zeros(64,nbin-1) * NaN
	for ibin = 1 : (nbin-1)
		ind = (cr .>= bin[ibin]) .& (cr .<= bin[ibin+1])
		for ipre = 1 : 64
			if !iszero(sum(ind)); rmean[ipre,ibin] = mean(rh[ipre,ind]) end
		end
	end
	# rnan = reshape(.!isnan.(sum(rmean,dims=1)),:)
	# binp  = collect(binp)[rnan]
	# rmean = rmean[:,rnan]
end

# ╔═╡ d3c6ff52-e86c-4eed-a892-545f92a0d305
begin
	pplt.close(); f2,a2 = pplt.subplots(aspect=2,axwidth=4)
	
	c2 = a2[1].pcolormesh(binp,pl,rmean,levels=5:10:95,extend="both",cmap="Blues")
	a2[1].format(
		xlim=(20,100),
		ylim=(1010,20),yscale="log"
	)
	a2[1].colorbar(c2,loc="r")
	a2[1].format(
		grid=true,
		ylabel="Pressure / hPa",
		xlabel="Column Relative Humidity / %",
		suptitle="$(expname) | WTG Simulations"
	)
	
	f2.savefig(plotsdir("$(expname)-WTGsimCRHvsP.png"),transparent=false,dpi=200)
	load(plotsdir("$(expname)-WTGsimCRHvsP.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─b6892634-9199-11eb-38d5-8dda8da16ed7
# ╠═d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╠═7160ccfe-2189-4341-9e83-e700773f3d3e
# ╟─55230f4a-7661-11eb-1c37-8b022b95e08e
# ╠═f706782e-6b0e-40fa-ba21-ff87f2594f5a
# ╟─d3c6ff52-e86c-4eed-a892-545f92a0d305
