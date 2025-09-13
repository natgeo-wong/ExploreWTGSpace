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
	"KGW" => "(KGW) Damped Gravity Wave [Kuang et al., 2008]",
	"TDG" => "(TDG) Time-Dependent Damped Gravity Wave [Kuang et al., 2008]",
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
	nmode   = 5
	blues_WTG = pplt.get_colors("Blues_r",(nmode+2))
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
	fts,ats = pplt.subplots(nrows=2,aspect=3,axwidth=5)

	if checkschemeDGW(wtgscheme)
		fnc = joinpath(wtgscheme,expname,"$(dampingstrprnt(configWTG[iconWTG])).nc")
	else
		fnc = joinpath(wtgscheme,expname,"$(relaxscalestrprnt(configWTG[iconWTG])).nc")
	end

	imem = 14
	ds_dgwprcp = NCDataset(datadir("wwtg",fnc))
	dt   = ds_dgwprcp["time"][:]
	p    = ds_dgwprcp["p"][:,imem]
	wwtg = ds_dgwprcp["wwtg"][:,:,imem]
	ptrp = ds_dgwprcp["ptrop"][:,imem]
	ztrp = ds_dgwprcp["ztrop"][:,imem]
	wₙ   = ds_dgwprcp["wₙ"][:,:,imem] ./ sum(abs.(ds_dgwprcp["wₙ"][:,:,imem]),dims=2)
	close(ds_dgwprcp)

	ats[1].pcolormesh(dt,p,wwtg,levels=-0.05:0.005:0.05)
	ats[1].plot(dt,ptrp)
	ats[2].plot(dt,wₙ)
	pnl = ats[2].panel("r")
	pnl.scatter(zeros(10),mean(abs.(wₙ),dims=1)')

	ats[1].format(ylim=(1000,25),yscale="log")
	ats[2].format(ylim=(-1,1))
	for ax in ats
		ax.format(
			xlim=(240,250),xlabel="Days"
		)
	end

	fts.savefig(
		plotsdir("11-wn-$(expname)-$wtgscheme.png"),
		transparent=false,dpi=400
	)
		
	load(plotsdir("11-wn-$(expname)-$wtgscheme.png"))
end

# ╔═╡ a08e9b0e-fb3b-4d84-b9ae-7c618941dd3f
md"
### A. Compilation
"

# ╔═╡ 59cdcd12-394b-4f66-bb83-be11959da38c
begin
	pplt.close()
	f2,a2 = pplt.subplots(ncols=1,aspect=0.4,axwidth=0.9)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE","P1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:]
	prcpRCEμ_P = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)
	ds_rceprcp = NCDataset(datadir("precipitation","RCE","T1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:]
	prcpRCEμ_T = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	function ncname(wtgvec,ii)

		iconWTG = wtgvec[ii]

		if checkschemeDGW(wtgscheme)
			return "$(dampingstrprnt(iconWTG)).nc"
		else
			return "$(relaxscalestrprnt(iconWTG)).nc"
		end 
		
	end

	for iconfig in 1 : nconWTG

		fnc = ncname(configWTG,iconfig)
			
		ds = NCDataset(datadir("wwtg","SPC","P1282km300V64",fnc))
		wₙ = ds["wₙ"][:,:,:]
		close(ds)
		
		ds = NCDataset(datadir("precipitation","SPC","P1282km300V64",fnc))
		prcp = ds["precipitation"][1001:2000,:]
		close(ds)
		
		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_P / 0.9)
		jj = prcp .> (prcpRCEμ_P * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				a2[1].scatter(-mean(wₙ[imode,ii])./mean(wₙ[1,ii]),configWTG[iconfig],c=blues_WTG[imode+1],s=5,zorder=5)
			end
			if !iszero(length(jj))
				a2[1].scatter(mean(wₙ[imode,jj])./mean(wₙ[1,jj]),configWTG[iconfig],c=blues_WTG[imode+1],s=5,zorder=5)
			end
		end

	end

	
	for ax in a2
		ax.format(
			ylim=(0.005,2000),yscale="log",xlim=(-1.1,1.1)
		)
	end

	f2.savefig(
		plotsdir("11-compiledwn-$(expname)-$wtgscheme.png"),
		transparent=false,dpi=400
	)
		
	load(plotsdir("11-compiledwn-$(expname)-$wtgscheme.png"))
end

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
# ╟─a08e9b0e-fb3b-4d84-b9ae-7c618941dd3f
# ╠═59cdcd12-394b-4f66-bb83-be11959da38c
