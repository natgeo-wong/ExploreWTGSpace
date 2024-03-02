### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
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
	using DSP
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
# 2a. Transitioning from RCE to WTG

In this notebook, we investigate and develop a way to implement the WTG forcing gradually in the System of Atmospheric Modelling.  The sudden introduction of WTG large-scale forcing often causes a model to enter a \"shocked\" state that unnaturally forces the model into a different state.  Here, we develop a method that gradually increases the strength of the WTG momentum-damping parameter from a near/pseudo-RCE state.
"

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	configWTG = [0.02,0.05,0.1,0.2,0.5,1,2,5]
	nconWTG = length(configWTG)
	blues_WTG = pplt.get_colors("Blues",(nconWTG))
	lgd_WTG = Dict("frame"=>false,"ncols"=>3)
	md"Loading time dimension and defining the damping experiments ..."
	prcp = zeros(Float32,32,32,480,8)
	md"Some preallocation and setup"
end

# ╔═╡ be7a062d-7424-48d9-bea2-a4c58b132030
md"Do Animation? $(@bind isanim PlutoUI.Slider(0:1,default=0))"

# ╔═╡ dc06d1f6-b3ab-4c7e-9f30-cfe1bb5e9cbd
begin
	if isone(isanim)
		for ic in 1 : nconWTG
	
			prcpii = @view prcp[:,:,:,ic]
		
			fnc = datadir(
				"DGW","T0642km300V64",dampingstrprnt(configWTG[ic]),"OUT_2D",
				"DGW_ExploreWTGSpace-T0642km300V64-member01-OUT2D_32.2Dcom_1.nc"
			)
			ds_dgwprcp = NCDataset(datadir(fnc))
			NCDatasets.load!(ds_dgwprcp["TB"].var,prcpii,:,:,:)
			close(ds_dgwprcp)
	
		end
	
		for ii = 1 : 480
			pplt.close()
			fts,ats = pplt.subplots(nrows=2,ncols=4,axwidth=1)
	
			lvls = 200 : 5 : 300
	
			for ic = 1 : nconWTG
				ats[ic].pcolormesh(1:2:63,1:2:63,prcp[:,:,ii,ic],levels=lvls,cmap="Blues",extend="both")
			end
		
			for ax in ats
				ax.format(
					xlim=(0,64),ylim=(0,64),
				)
			end
		
			fts.savefig(plotsdir("02Doutput-$ii.png"),transparent=false,dpi=150)
		end
	end
end

# ╔═╡ 42024543-700e-42bf-acd2-b9edabe73cf6
load(plotsdir("02Doutput-1.png"))

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─be7a062d-7424-48d9-bea2-a4c58b132030
# ╟─dc06d1f6-b3ab-4c7e-9f30-cfe1bb5e9cbd
# ╟─42024543-700e-42bf-acd2-b9edabe73cf6
