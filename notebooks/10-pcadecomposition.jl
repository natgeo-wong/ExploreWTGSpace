### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 681658b0-5914-11eb-0d65-bbace277d145
begin
	using Pkg; Pkg.activate()
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 6dce35fc-5914-11eb-0ce2-0d4e164e1898
begin
	@quickactivate "ExploreWTGSpace"
	using MultivariateStats
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the ExploreWTGSpace paper ..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 3. Time-Series and Power Spectrum

In this notebook, we investigate and develop a way to implement the WTG forcing gradually in the System of Atmospheric Modelling.  The sudden introduction of WTG large-scale forcing often causes a model to enter a \"shocked\" state that unnaturally forces the model into a different state.  Here, we develop a method that gradually increases the strength of the WTG momentum-damping parameter from a near/pseudo-RCE state.
"

# ╔═╡ f188190f-81bf-4b29-b479-38e9a85a997c
@bind wtgscheme Select([
	"DGW" => "(DGW) Damped Gravity Wave [Blossey et al., 2009]",
	"KGW" => "(KGW) Damped Gravity Waves [Kuang et al., 2008]",
	"TDG" => "(TDG) Time-Dependent Damped Gravity Waves [Kuang et al., 2008]",
	"TGR" => "(TGR) Temperature Gradient Relaxation [Raymond and Zeng, 2005]",
	"SPC" => "(SPC) Spectral TGR [Herman and Raymond, 2014]",
])

# ╔═╡ b7a79d4e-4007-4c55-99cd-33abe6ee9f32
@bind prefix Select([
	"P" => "Perpetual Insolation (P)",
	"D" => "Diurnal Insolation (D)",
	"T" => "Temperature Tendency (T)",
])

# ╔═╡ 292ff637-7f96-4d9b-beeb-8b3d7b28a218
if prefix == "D"
	md"Toggle Domain Size: $(@bind hres PlutoUI.Slider(1:2,default=1))"
else
	md"Toggle Domain Size: $(@bind hres PlutoUI.Slider(0.5:0.5:1,default=1))"
end

# ╔═╡ 026110d9-55be-484a-b962-1aed19528933
md"Domain Size = $(128*hres) x $(128*hres) km"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
begin
	expname = "$(prefix)$(@sprintf("%03d",128*hres))2km300V64"

md"**Experiment Set:** $expname"
end

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	if checkschemeDGW(wtgscheme)
		configWTG = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
	else
		configWTG = [
			0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
			1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),
			10,10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2)
		]
	end
	nconWTG = length(configWTG)
	blues_WTG = pplt.get_colors("Blues",(nconWTG))
	lgd_WTG = Dict("frame"=>false,"ncols"=>3)
	md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 228c6070-0682-48ec-9870-d83ba44cf128
md"WTG Strength: $(@bind iconWTG PlutoUI.Slider(1:nconWTG,default=1))"

# ╔═╡ dbdffc08-71d1-46e8-b583-7a18fa351d0b
configWTG[iconWTG]

# ╔═╡ dc06d1f6-b3ab-4c7e-9f30-cfe1bb5e9cbd
begin
	pplt.close()
	fts,ats = pplt.subplots(aspect=3,axwidth=5)

	if checkschemeDGW(wtgscheme)
		fnc = joinpath(wtgscheme,expname,"$(dampingstrprnt(configWTG[iconWTG])).nc")
	else
		fnc = joinpath(wtgscheme,expname,"$(relaxscalestrprnt(configWTG[iconWTG])).nc")
	end
	
	ds_dgwprcp = NCDataset(datadir("wwtg",fnc))
	dt   = ds_dgwprcp["time"][(end-4799):end]
	p    = ds_dgwprcp["p"][:,1]
	wwtg = ds_dgwprcp["wwtg"][:,(end-4799):end,11] * 1000
	ptrp = ds_dgwprcp["ptrop"][(end-4799):end,11]
	close(ds_dgwprcp)

	ats[1].pcolormesh(dt,p,wwtg)
	ats[1].plot(dt,ptrp)

	for ax in ats
		ax.format(
			xlim=(225,250),yscale="log",ylim=(1000,25),
			ylabel="Pressure / hPa",xlabel="Days"
		)
	end

	fts.savefig(
		plotsdir("10-wwtg-$(expname)-$wtgscheme.png"),
		transparent=false,dpi=400
	)
		
	load(plotsdir("10-wwtg-$(expname)-$wtgscheme.png"))
end

# ╔═╡ 2dfdd378-a118-4335-9fc9-98675b027757
M = fit(PCA, wwtg; pratio=0.99)

# ╔═╡ 98c7a489-2f8b-4c73-b96e-f396a0bc573a
evec = eigvecs(M);

# ╔═╡ 5426c985-b6f9-4143-9ea5-87e31ea9f784
begin
	pplt.close(); f2,a2 = pplt.subplots(aspect=0.5,axwidth=1)
	
	a2[1].plot(evec[:,1],p)
	# a2[1].plot(evec[:,2],p)
	# a2[1].plot(evec[:,3],p)
	# a2[1].plot(evec[:,4],p)
	a2[1].plot([-1,1]*0.1,ones(2)*mean(ptrp))
	# a2[1].plot(evec[:,3],p)
	a2[1].format(ylim=(1000,25),yscale="log",xlim=(-0.5,0.5))
	
	f2.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ 082f10dc-e9e6-4d1c-8186-97b771a2d8ad
size(ptrp)

# ╔═╡ 09260a1d-4c76-43ce-b8c7-267fa534f0b4
Y = predict(M,wwtg)

# ╔═╡ 84b4c52c-1820-4be9-9bf2-037b91d8c897
begin
	pplt.close()
	f3,a3 = pplt.subplots(aspect=3,axwidth=5)

	a3[1].pcolormesh(dt,p,reconstruct(M,Y))

	for ax in a3
		ax.format(
			xlim=(225,250),yscale="log",ylim=(1000,25),
			ylabel="Pressure / hPa",xlabel="Days"
		)
	end

	f3.savefig(
		plotsdir("10-wwtg-$(expname)-$wtgscheme.png"),
		transparent=false,dpi=400
	)
		
	load(plotsdir("10-wwtg-$(expname)-$wtgscheme.png"))
end

# ╔═╡ 9b866d0e-3ab0-4673-b482-56e9125a42ae
begin
	pplt.close()
	f4,a4 = pplt.subplots(aspect=3,axwidth=5)

	a4[1].plot(dt,Y')

	for ax in a4
		ax.format(
			xlim=(225,250),xlabel="Days"
		)
	end

	f4.savefig(
		plotsdir("10-wwtg-$(expname)-$wtgscheme.png"),
		transparent=false,dpi=400
	)
		
	load(plotsdir("10-wwtg-$(expname)-$wtgscheme.png"))
end

# ╔═╡ 834fb00a-2e75-486b-a265-4d7b94e1a816
principalvars(M) ./ (tprincipalvar(M) + tresidualvar(M))

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─f188190f-81bf-4b29-b479-38e9a85a997c
# ╟─b7a79d4e-4007-4c55-99cd-33abe6ee9f32
# ╟─292ff637-7f96-4d9b-beeb-8b3d7b28a218
# ╟─026110d9-55be-484a-b962-1aed19528933
# ╟─d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─228c6070-0682-48ec-9870-d83ba44cf128
# ╟─dbdffc08-71d1-46e8-b583-7a18fa351d0b
# ╠═dc06d1f6-b3ab-4c7e-9f30-cfe1bb5e9cbd
# ╠═2dfdd378-a118-4335-9fc9-98675b027757
# ╟─98c7a489-2f8b-4c73-b96e-f396a0bc573a
# ╠═5426c985-b6f9-4143-9ea5-87e31ea9f784
# ╠═082f10dc-e9e6-4d1c-8186-97b771a2d8ad
# ╠═09260a1d-4c76-43ce-b8c7-267fa534f0b4
# ╠═84b4c52c-1820-4be9-9bf2-037b91d8c897
# ╠═9b866d0e-3ab0-4673-b482-56e9125a42ae
# ╠═834fb00a-2e75-486b-a265-4d7b94e1a816
