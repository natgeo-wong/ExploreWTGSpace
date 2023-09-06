### A Pluto.jl notebook ###
# v0.19.16

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
	using NCDatasets
	using PlutoUI
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
# 2a. Transitioning from RCE to WTG

In this notebook, we investigate and develop a way to implement the WTG forcing gradually in the System of Atmospheric Modelling.  The sudden introduction of WTG large-scale forcing often causes a model to enter a \"shocked\" state that unnaturally forces the model into a different state.  Here, we develop a method that gradually increases the strength of the WTG momentum-damping parameter from a near/pseudo-RCE state.
"

# ╔═╡ f188190f-81bf-4b29-b479-38e9a85a997c
@bind wtgscheme Select([
	"WTG" => "(WTG) Weak Temperature Gradient [Raymond and Zeng, 2005]",
	"CON" => "(CON) Constant WTG [Daleu et al. 2015]",
	"SIN" => "(SIN) Full Sine-Curve",
	"DCM" => "(DCM) Decomposed WTG Sine-Curve",
])

# ╔═╡ b7a79d4e-4007-4c55-99cd-33abe6ee9f32
@bind prefix Select([
	"P" => "Perpetual Insolation (P)",
	"D" => "Diurnal Insolation (D)",
	"T" => "Non-interactive Radiation (T)",
	"S" => "Bulk-surface Fluxes (S)",
])

# ╔═╡ fc813056-4cb6-4db8-ba14-2fd2243b7c0d
md"Toggle Domain Size $(@bind islarge PlutoUI.Slider(64:64:128,default=128,show_value=true)) km"

# ╔═╡ 60a1f908-7b2d-4ed2-bfc1-3f3a416c7c88
md"Toggle Horizontal Resolution: $(@bind hres PlutoUI.Slider(-1:1,default=0))"

# ╔═╡ bdf1a99a-24f9-4826-af25-3603686f23ad
md"Sea Surface Temperature: $(@bind sst PlutoUI.Slider(295:5:305,default=300, show_value=true))"

# ╔═╡ 292ff637-7f96-4d9b-beeb-8b3d7b28a218
md"Coarse Vertical Grid? $(@bind iscvg PlutoUI.Slider(0:1))"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
begin
	domsize = @sprintf("%03d",islarge)

	if wtgscheme == "WTG"

		if islarge == 64
			  res = Int(2. ^(hres-1)*2)
		else; res = Int(2. ^hres*2)
		end

		if iszero(iscvg)
			  vgrd = 64
		else; vgrd = 28
		end
		expname = "$(prefix)$(domsize)$(res)km$(sst)V$(vgrd)"

	else

		if prefix == "T"
			expname = "T1282km300V64"
		elseif prefix == "S"
			expname = "S1284km300V64"
		else
			expname = "P1282km300V64"
		end

	end

md"**WTG Scheme:** $wtgscheme | **Experiment Set:** $expname"
end

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	configDGW = [0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000]
	configWTG = [
		0.01,0.02,0.02*sqrt(2.5),0.05,0.05*sqrt(2),
		0.1,0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
		1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),
		10,10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),100
	]
	nconDGW = length(configDGW)
	nconWTG = length(configWTG)
    blues_DGW = pplt.get_colors("Blues",(nconDGW+2))
	blues_WTG = pplt.get_colors("Blues",(nconWTG+2))
	lgd_DGW = Dict("frame"=>false,"ncols"=>2)
	lgd_WTG = Dict("frame"=>false,"ncols"=>3)
	md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 4581e38f-a680-4692-96ce-45d5e8799953
md"Create Image? $(@bind createimage PlutoUI.Slider(0:1))"

# ╔═╡ dc06d1f6-b3ab-4c7e-9f30-cfe1bb5e9cbd
begin
	if isone(createimage)
		pplt.close()
		fts,ats = pplt.subplots(nrows=2,aspect=2,axwidth=3.5)
	
		for ic in 1 : nconDGW
			config = "damping$(dampingstrprnt(configDGW[ic]))"
			imem = 0
	
			while imem < 15; imem += 1
				fnc = outstatname("DGW",expname,config,false,true,imem)
				if isfile(fnc)
					_,p,t = retrievedims_fnc(fnc); t = t .- floor(t[1])
					pr = retrievevar_fnc("PREC",fnc) / 24
					nt = length(t); nt = nt - mod(nt,8)
					# t  = dropdims(mean(reshape(t[1:nt],8,:),dims=1),dims=1)
					# pr = dropdims(mean(reshape(pr[1:nt],8,:),dims=1),dims=1)
					if imem == 1
						constr = @sprintf("%.e",configDGW[ic])
						ats[1].plot(
							t,pr,color=blues_DGW[ic+1],
							label=(L"$a_m =$" * " $(constr)" * L" day$^{-1}$"),
							legend="r",legend_kw=lgd_DGW
						)
					else
						ats[1].plot(t,pr,color=blues_DGW[ic+1])
					end
				end
			end
	
		end
		
		for ic in 1 : nconWTG
			config = relaxscalestrprnt(configWTG[ic])
			imem = 0
	
			while imem < 15; imem += 1
				fnc = outstatname(wtgscheme,expname,config,false,true,imem)
				if isfile(fnc)
					_,p,t = retrievedims_fnc(fnc); t = t .- floor(t[1])
					pr = retrievevar_fnc("PREC",fnc) / 24
					nt = length(t); nt = nt - mod(nt,8)
					# t  = dropdims(mean(reshape(t[1:nt],8,:),dims=1),dims=1)
					# pr = dropdims(mean(reshape(pr[1:nt],8,:),dims=1),dims=1)
					if imem == 1
						constr = @sprintf("%.1e",configWTG[ic])
						ats[2].plot(
							t,pr,color=blues_WTG[ic+1],
							label=(L"$\tau =$" * " $(constr) hr"),
							legend="r",legend_kw=lgd_WTG
						)
					else
						ats[2].plot(t,pr,color=blues_WTG[ic+1])
					end
				end
			end
	
		end
	
		# for imem = 1 : 10
		# 	fnc = outstatname("RCE",expname,"",false,true,imem)
		# 	if isfile(fnc)
		# 		_,p,t = retrievedims_fnc(fnc); t = t .- floor(t[1])
		# 		pr = retrievevar_fnc("PREC",fnc) / 24
		# 		if imem == 1
		# 			ats[1].plot(
		# 				t,pr,color="k",label=("RCE"),
		# 				legend="r",legend_kw=lgd_DGW
		# 			)
		# 			ats[2].plot(
		# 				t,pr,color="k",label=("RCE"),
		# 				legend="r",legend_kw=lgd_WTG
		# 			)
		# 		else
		# 			ats[1].plot(t,pr,color="k")
		# 			ats[2].plot(t,pr,color="k")
		# 		end
		# 	end
		# end

		for ax in ats
			ax.format(
				xlim=(000,250),#yscale="symlog",yscale_kw=Dict("linthresh"=>0.01),
				ylim=(0,2),
				ylabel=L"Rainfall Rate / mm hr$^{-1}$",xlabel="Days"
			)
		end

		ats[1].format(ultitle="(a) Precipitation Time-Series (DGW)")
		ats[2].format(ultitle="(b) Precipitation Time-Series ($wtgscheme)")

		fts.savefig(
			plotsdir("02b-rce2wtg-$(expname)-$wtgscheme.png"),
			transparent=false,dpi=400
		)
		
	end
	load(plotsdir("02b-rce2wtg-$(expname)-$wtgscheme.png"))
end

# ╔═╡ d8db43e6-e537-4348-97f5-da3f7dc76e0f
begin
	if isone(createimage)
		pplt.close()
		fstd,astd = pplt.subplots(ncols=2,aspect=3,axwidth=2)
	
		for ic in 1 : nconDGW
			config = "damping$(dampingstrprnt(configDGW[ic]))"
			imem = 0
	
			while imem < 15; imem += 1
				fnc = outstatname("DGW",expname,config,false,true,imem)
				if isfile(fnc)
					pr = retrievevar_fnc("PREC",fnc) / 24
					if length(pr) == 2000
						astd[1].errorbar(configDGW[ic],0,std(pr[(end-799):end]),c="b")
					end
				end
			end
	
		end
		
		for ic in 1 : nconWTG
			config = relaxscalestrprnt(configWTG[ic])
			imem = 0
	
			while imem < 15; imem += 1
				fnc = outstatname(wtgscheme,expname,config,false,true,imem)
				if isfile(fnc)
					pr = retrievevar_fnc("PREC",fnc) / 24
					if length(pr) == 2000
						astd[2].errorbar(configWTG[ic],0,std(pr[(end-799):end]),c="b")
					end
				end
			end
	
		end

		for ax in astd
			ax.format(xscale="log")
		end

		astd[1].format(urtitle="(a) Standard Deviation (DGW)")
		astd[2].format(urtitle="(b) Standard Deviation ($wtgscheme)")

		fstd.savefig(
			plotsdir("02b-stddeviation-$(expname)-$wtgscheme.png"),
			transparent=false,dpi=400
		)
		
	end
	load(plotsdir("02b-stddeviation-$(expname)-$wtgscheme.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─f188190f-81bf-4b29-b479-38e9a85a997c
# ╟─b7a79d4e-4007-4c55-99cd-33abe6ee9f32
# ╟─fc813056-4cb6-4db8-ba14-2fd2243b7c0d
# ╟─60a1f908-7b2d-4ed2-bfc1-3f3a416c7c88
# ╟─bdf1a99a-24f9-4826-af25-3603686f23ad
# ╟─292ff637-7f96-4d9b-beeb-8b3d7b28a218
# ╟─d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─4581e38f-a680-4692-96ce-45d5e8799953
# ╟─dc06d1f6-b3ab-4c7e-9f30-cfe1bb5e9cbd
# ╠═d8db43e6-e537-4348-97f5-da3f7dc76e0f
